classdef PulseStreamerHttpNoJava < handle
    %PULSESTREAMERHTTPNOJAVA JSON-RPC over HTTP/1.1 using tcpclient, for MATLAB without Java.
    %
    %   h = PulseStreamerHttpNoJava('132.64.56.113', 8050);
    %   reply = h.post('/json-rpc', jsonText);
    %
    % Why it exists: from R2026b MATLAB no longer installs a Java runtime, and urlread2
    % (PulseStreamer's transport since it was written) is built on java.net.
    % PulseStreamer.httpRequest keeps using urlread2 whenever Java is available -- the
    % proven path on the hardware -- and uses this class only when it is not.
    %
    % It sends the request line and headers urlread2's Java HttpURLConnection sends and,
    % like Java, keeps the connection open between calls: opening a tcpclient costs
    % ~60-100 ms, a request on an open connection ~2 ms, and the Pulse Streamer is polled
    % in a loop. Keeping a connection open is made safe by:
    %   - closing it after IDLE_CLOSE_SECONDS without use (Java does the same after ~5 s),
    %     so a connection the server may have timed out is never reused, and an idle
    %     connection never ties up the server;
    %   - on R2026b+, where tcpclient can report a closed peer (Connected), checking it
    %     before reuse;
    %   - if a REUSED connection fails before a single reply byte arrived, the server
    %     never processed the request, so it is sent once more on a fresh connection. It
    %     is never re-sent once reply bytes have arrived, so a command cannot run twice.
    % Replies with Content-Length, chunked encoding, or a body that ends when the server
    % closes are all handled. A non-2xx status is returned, not thrown, like urlread2.
    %
    % Base MATLAB only (tcpclient, timer); no Java and no toolbox.

    properties (Constant)
        IDLE_CLOSE_SECONDS = 2      % close a kept-alive connection after this long unused
    end

    properties (SetAccess = private)
        host                        % char
        port                        % double
        timeout                     % seconds, for connecting and for each reply
    end

    properties (Access = private)
        t = []                      % tcpclient while a connection is open
        idleTimer = []              % closes the connection when unused
        busy = false                % true while a request is in flight
        replyBytes = 0              % bytes of the current reply received so far
    end

    methods
        function obj = PulseStreamerHttpNoJava(host, port, timeout)
            if nargin < 3 || isempty(timeout)
                timeout = 10;
            end
            obj.host = char(host);
            obj.port = port;
            obj.timeout = timeout;
            obj.idleTimer = timer('ExecutionMode', 'singleShot', ...
                                  'StartDelay', obj.IDLE_CLOSE_SECONDS, ...
                                  'TimerFcn', @(~, ~) obj.closeIfIdle(), ...
                                  'Name', 'PulseStreamerHttpNoJava idle close');
        end

        function [body, status] = post(obj, path, payload)
            stop(obj.idleTimer);
            obj.busy = true;
            cleanupBusy = onCleanup(@() obj.setIdle());
            bytes = unicode2native(char(payload), 'UTF-8');
            request = [uint8(sprintf(['POST %s HTTP/1.1\r\n' ...
                                      'Content-Type: application/x-www-form-urlencoded\r\n' ...
                                      'Host: %s:%d\r\n' ...
                                      'Accept: text/html, image/gif, image/jpeg, *; q=.2, */*; q=.2\r\n' ...
                                      'Connection: keep-alive\r\n' ...
                                      'Content-Length: %d\r\n\r\n'], ...
                                     path, obj.host, obj.port, numel(bytes))), ...
                       reshape(bytes, 1, [])];
            for attempt = 1:2
                reused = obj.isOpen();
                if ~reused
                    obj.open();
                end
                obj.replyBytes = 0;
                try
                    write(obj.t, request);
                    [body, status, keepAlive] = obj.readResponse();
                catch err
                    obj.close();
                    if reused && obj.replyBytes == 0 && attempt == 1
                        continue        % stale connection, request never processed: retry once
                    end
                    rethrow(err);
                end
                if keepAlive
                    start(obj.idleTimer);
                else
                    obj.close();
                end
                return
            end
        end

        function close(obj)
            if ~isempty(obj.t)
                try
                    delete(obj.t);      % deleting a tcpclient closes the connection
                catch
                end
                obj.t = [];
            end
        end

        function delete(obj)
            obj.close();
            if ~isempty(obj.idleTimer) && isvalid(obj.idleTimer)
                stop(obj.idleTimer);
                delete(obj.idleTimer);
            end
        end
    end

    methods (Access = private)
        function closeIfIdle(obj)
            if ~obj.busy
                obj.close();
            end
        end

        function setIdle(obj)
            obj.busy = false;
        end

        function tf = isOpen(obj)
            tf = ~isempty(obj.t) && isvalid(obj.t);
            if tf && isprop(obj.t, 'Connected') && ~obj.t.Connected   % R2026b+
                obj.close();
                tf = false;
            end
        end

        function open(obj)
            args = {'Timeout', obj.timeout, 'ConnectTimeout', obj.timeout, 'EnableTransferDelay', false};
            try
                % R2026b+: let tcpclient notice a server-closed connection (see isOpen).
                obj.t = tcpclient(obj.host, obj.port, args{:}, 'DetectServerDisconnect', true);
            catch
                obj.t = tcpclient(obj.host, obj.port, args{:});
            end
        end

        function [body, status, keepAlive] = readResponse(obj)
            started = tic;
            raw = obj.readMore(started);
            crlf2 = strfind(char(raw), sprintf('\r\n\r\n'));
            while isempty(crlf2)
                raw = [raw, obj.readMore(started)]; %#ok<AGROW>
                crlf2 = strfind(char(raw), sprintf('\r\n\r\n'));
            end
            head = char(raw(1:crlf2(1)-1));
            rest = raw(crlf2(1)+4:end);
            tok = regexp(head, '^HTTP/(\d\.\d)\s+(\d{3})', 'tokens', 'once');
            if isempty(tok)
                error('PulseStreamerHttpNoJava:badReply', 'Not an HTTP reply from %s:%d: %s', ...
                      obj.host, obj.port, head(1:min(end, 80)));
            end
            status = str2double(tok{2});
            connHdr = regexp(head, '(?im)^connection:\s*([^\r\n]*)', 'tokens', 'once');
            if strcmp(tok{1}, '1.0')
                keepAlive = ~isempty(connHdr) && contains(lower(connHdr{1}), 'keep-alive');
            else
                keepAlive = isempty(connHdr) || ~contains(lower(connHdr{1}), 'close');
            end
            lenTok = regexp(head, '(?im)^content-length:\s*(\d+)', 'tokens', 'once');
            isChunked = ~isempty(regexp(head, '(?im)^transfer-encoding:[^\r\n]*chunked', 'once'));

            if ~isempty(lenTok)
                n = str2double(lenTok{1});
                while numel(rest) < n
                    rest = [rest, obj.readMore(started)]; %#ok<AGROW>
                end
                data = rest(1:n);
            elseif isChunked
                data = uint8([]);
                while true
                    lineEnd = strfind(char(rest), sprintf('\r\n'));
                    while isempty(lineEnd)
                        rest = [rest, obj.readMore(started)]; %#ok<AGROW>
                        lineEnd = strfind(char(rest), sprintf('\r\n'));
                    end
                    n = hex2dec(strtrim(strtok(char(rest(1:lineEnd(1)-1)), ';')));
                    rest = rest(lineEnd(1)+2:end);
                    if n == 0
                        break                               % trailers, if any, are ignored
                    end
                    while numel(rest) < n + 2
                        rest = [rest, obj.readMore(started)]; %#ok<AGROW>
                    end
                    data = [data, rest(1:n)]; %#ok<AGROW>
                    rest = rest(n+3:end);                   % skip the chunk and its CRLF
                end
            else
                % No length given: the body ends when the server closes the connection.
                keepAlive = false;
                data = rest;
                idle = tic;
                while toc(idle) < 0.05
                    if obj.t.NumBytesAvailable > 0
                        data = [data, obj.readMore(started)]; %#ok<AGROW>
                        idle = tic;
                    else
                        pause(0.0002);
                    end
                end
            end
            body = native2unicode(data, 'UTF-8');
        end

        function b = readMore(obj, started)
            % Wait for at least one byte (polling: ~1 ms faster per request than a
            % blocking read), then take everything already waiting.
            while obj.t.NumBytesAvailable == 0
                if toc(started) > obj.timeout
                    error('PulseStreamerHttpNoJava:timeout', ...
                          'No complete reply from the Pulse Streamer at %s:%d within %g s.', ...
                          obj.host, obj.port, obj.timeout);
                end
                pause(0.0002);
            end
            b = read(obj.t, obj.t.NumBytesAvailable, 'uint8');
            obj.replyBytes = obj.replyBytes + numel(b);
        end
    end
end
