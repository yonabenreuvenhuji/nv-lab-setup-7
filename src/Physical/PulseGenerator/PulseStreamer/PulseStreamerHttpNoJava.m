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
    % One connection per request. The Pulse Streamer closes the connection after every
    % reply and sends no Content-Length, so the reply ends when the server closes
    % (measured on setup 7, R2026b). Reusing connections is therefore impossible, and
    % this class sends "Connection: close" and never tries.
    %
    % Finding the end of a close-delimited reply:
    %   - R2026b+: tcpclient's DetectServerDisconnect reports the close (Connected
    %     becomes false), so the reply is complete as soon as it happens (~59 ms per
    %     request on setup 7, nearly all of it opening the connection). The close is
    %     also reported through ErrorOccurredFcn, whose default prints "An unexpected
    %     error was encountered. Additional Information: End of file"; a silent
    %     handler is installed because here the close is the normal end of a reply.
    %   - Older releases: the reply is taken as complete after 50 ms without new bytes.
    % Replies with Content-Length or chunked encoding are handled too. A non-2xx
    % status is returned, not thrown, like urlread2.
    %
    % Base MATLAB only (tcpclient); no Java and no toolbox.

    properties (Constant)
        QUIET_SECONDS = 0.05        % no server-close detection: reply done after this silence
    end

    properties (SetAccess = private)
        host                        % char
        port                        % double
        timeout                     % seconds, for connecting and for each reply
    end

    methods
        function obj = PulseStreamerHttpNoJava(host, port, timeout)
            if nargin < 3 || isempty(timeout)
                timeout = 10;
            end
            obj.host = char(host);
            obj.port = port;
            obj.timeout = timeout;
        end

        function [body, status] = post(obj, path, payload)
            bytes = unicode2native(char(payload), 'UTF-8');
            request = [uint8(sprintf(['POST %s HTTP/1.1\r\n' ...
                                      'Content-Type: application/x-www-form-urlencoded\r\n' ...
                                      'Host: %s:%d\r\n' ...
                                      'Accept: text/html, image/gif, image/jpeg, *; q=.2, */*; q=.2\r\n' ...
                                      'Connection: close\r\n' ...
                                      'Content-Length: %d\r\n\r\n'], ...
                                     path, obj.host, obj.port, numel(bytes))), ...
                       reshape(bytes, 1, [])];
            [t, detectsClose] = obj.open();
            closeConnection = onCleanup(@() delete(t));    % deleting a tcpclient closes it
            write(t, request);
            [body, status] = obj.readResponse(t, detectsClose);
        end
    end

    methods (Access = private)
        function [t, detectsClose] = open(obj)
            args = {'Timeout', obj.timeout, 'ConnectTimeout', obj.timeout, 'EnableTransferDelay', false};
            try
                t = tcpclient(obj.host, obj.port, args{:}, 'DetectServerDisconnect', true);   % R2026b+
                detectsClose = true;
            catch err
                if ~strcmp(err.identifier, 'MATLAB:InputParser:UnmatchedParameter') ...
                        && ~contains(err.message, 'DetectServerDisconnect')
                    rethrow(err)        % a real connection failure, not an older release
                end
                t = tcpclient(obj.host, obj.port, args{:});
                detectsClose = false;
            end
            if detectsClose
                % The server closing is how every reply ends; don't print it as an error.
                % Set before the request is written, so before the server can close.
                t.ErrorOccurredFcn = @(~, ~) [];
            end
        end

        function [body, status] = readResponse(obj, t, detectsClose)
            started = tic;
            raw = uint8([]);
            crlf2 = [];
            while isempty(crlf2)
                [b, closed] = obj.readMore(t, detectsClose, started);
                if closed
                    error('PulseStreamerHttpNoJava:closed', ...
                          'The Pulse Streamer at %s:%d closed the connection without a reply.', ...
                          obj.host, obj.port);
                end
                raw = [raw, b]; %#ok<AGROW>
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
            lenTok = regexp(head, '(?im)^content-length:\s*(\d+)', 'tokens', 'once');
            isChunked = ~isempty(regexp(head, '(?im)^transfer-encoding:[^\r\n]*chunked', 'once'));

            if ~isempty(lenTok)
                n = str2double(lenTok{1});
                while numel(rest) < n
                    rest = [rest, obj.readOrFail(t, detectsClose, started)]; %#ok<AGROW>
                end
                data = rest(1:n);
            elseif isChunked
                data = uint8([]);
                while true
                    lineEnd = strfind(char(rest), sprintf('\r\n'));
                    while isempty(lineEnd)
                        rest = [rest, obj.readOrFail(t, detectsClose, started)]; %#ok<AGROW>
                        lineEnd = strfind(char(rest), sprintf('\r\n'));
                    end
                    n = hex2dec(strtrim(strtok(char(rest(1:lineEnd(1)-1)), ';')));
                    rest = rest(lineEnd(1)+2:end);
                    if n == 0
                        break                               % trailers, if any, are ignored
                    end
                    while numel(rest) < n + 2
                        rest = [rest, obj.readOrFail(t, detectsClose, started)]; %#ok<AGROW>
                    end
                    data = [data, rest(1:n)]; %#ok<AGROW>
                    rest = rest(n+3:end);                   % skip the chunk and its CRLF
                end
            elseif detectsClose
                % The Pulse Streamer's style: the body ends when the server closes.
                data = rest;
                closed = false;
                while ~closed
                    [b, closed] = obj.readMore(t, detectsClose, started);
                    data = [data, b]; %#ok<AGROW>
                end
            else
                % Same, without close detection: done after QUIET_SECONDS of silence.
                data = rest;
                quiet = tic;
                while toc(quiet) < obj.QUIET_SECONDS
                    if t.NumBytesAvailable > 0
                        data = [data, read(t, t.NumBytesAvailable, 'uint8')]; %#ok<AGROW>
                        quiet = tic;
                    else
                        pause(0.0002);
                    end
                end
            end
            body = native2unicode(data, 'UTF-8');
        end

        function [b, closed] = readMore(obj, t, detectsClose, started)
            % Wait for at least one byte (polling: ~1 ms faster per request than a
            % blocking read), then take everything already waiting. closed is true, with
            % no bytes, once the server has closed and everything it sent was read.
            closed = false;
            while t.NumBytesAvailable == 0
                if detectsClose && ~t.Connected && t.NumBytesAvailable == 0
                    b = uint8([]);
                    closed = true;
                    return
                end
                if toc(started) > obj.timeout
                    error('PulseStreamerHttpNoJava:timeout', ...
                          'No complete reply from the Pulse Streamer at %s:%d within %g s.', ...
                          obj.host, obj.port, obj.timeout);
                end
                pause(0.0002);
            end
            b = read(t, t.NumBytesAvailable, 'uint8');
        end

        function b = readOrFail(obj, t, detectsClose, started)
            [b, closed] = obj.readMore(t, detectsClose, started);
            if closed
                error('PulseStreamerHttpNoJava:truncated', ...
                      'The Pulse Streamer at %s:%d closed the connection in the middle of a reply.', ...
                      obj.host, obj.port);
            end
        end
    end
end
