function results = awgTransportBench(awg, resources, nQueries)
%AWGTRANSPORTBENCH Measure per-transaction latency and arm time per VISA transport.
%
%   results = awgTransportBench(awg)
%   results = awgTransportBench(awg, {'USB0::...', 'TCPIP0::...::inst0::INSTR'})
%   results = awgTransportBench(awg, resources, 40)
%
% WHY: the sweep spends most of its time inside Run(), and Run() polls
% AWGControl:RSTate? in a loop until the AWG reports armed. So the measured
% ~450-600 ms "arm latency" is (true firmware arm time) + (number of polls x
% per-query latency). Those two are indistinguishable from the outside until you
% change the per-query latency and see whether the arm time moves.
%
% VXI-11 (…inst0::INSTR) is RPC-based and typically the slowest of the three
% transports an AWG70002A offers; USBTMC and raw SOCKET are usually much lower
% latency. If the arm time tracks the query latency, the arm is mostly transport
% and switching resource is a real win. If it does not move, the 450 ms really is
% firmware and the only fix is to stop re-arming per point.
%
% Reports per resource:
%   openTime      time to construct the session and read *IDN?
%   queryMedian   median AWGControl:RSTate? round trip
%   queryP90      90th percentile, because the tail is what a polling loop pays
%   armMedian     median stop-then-Run() cycle, i.e. what each sweep point costs
%
% SAFETY: this stops and re-arms the AWG a handful of times, exactly as a normal
% sweep does. The pulse generator is not triggering, so no MW is emitted. It does
% NOT touch loaded waveforms. It restores the original resource on the way out,
% including if it errors.
%
% Requires awg.useVisadev = true (the legacy interface cannot open a SOCKET
% resource cleanly and is what we are trying to get away from).

    if ~exist('nQueries', 'var') || isempty(nQueries)
        nQueries = 30;
    end
    if ~exist('resources', 'var') || isempty(resources)
        resources = { 'USB0::0x0699::0x0503::B040645::0::INSTR', ...
                      'TCPIP0::132.64.56.37::inst0::INSTR', ...
                      'TCPIP0::132.64.56.37::4000::SOCKET' };
    end
    if ~iscell(resources)
        resources = {resources};
    end
    if ~awg.useVisadev
        error('Set awg.useVisadev = true first -- this benchmark is about the modern interface.');
    end

    original = awg.getResource();
    fprintf('\n=== awgTransportBench (original resource: %s) ===\n\n', original);
    fprintf('%-46s %-10s %-12s %-10s %-10s\n', ...
            'resource', 'open (s)', 'query med', 'query p90', 'arm med');

    results = struct('resource',{},'openTime',{},'queryMedian',{}, ...
                     'queryP90',{},'armMedian',{},'ok',{});

    for r = 1:length(resources)
        res = resources{r};
        entry = struct('resource', res, 'openTime', NaN, 'queryMedian', NaN, ...
                       'queryP90', NaN, 'armMedian', NaN, 'ok', false);
        try
            tOpen = tic;
            awg.setResource(res);
            idn = awg.queryRaw('*IDN?', 10);
            entry.openTime = toc(tOpen);
            if ~contains(idn, 'AWG70002A')
                fprintf('%-46s  SKIPPED - *IDN? returned "%s"\n', shorten(res), strtrim(idn));
                results(end+1) = entry;  %#ok<AGROW>
                continue
            end

            % --- per-query latency, the thing a polling loop multiplies --------
            t = zeros(1, nQueries);
            for k = 1:nQueries
                tq = tic;
                awg.queryRaw('AWGControl:RSTate?', 10);
                t(k) = toc(tq);
            end
            entry.queryMedian = median(t);
            entry.queryP90    = prctileSimple(t, 90);

            % --- full stop-then-arm cycle, what each sweep point pays ----------
            nArms = 5;
            ta = zeros(1, nArms);
            for k = 1:nArms
                awg.sendCommand('awgcontrol:stop:immediate');
                tArm = tic;
                awg.Run;
                ta(k) = toc(tArm);
            end
            entry.armMedian = median(ta);
            entry.ok = true;

            fprintf('%-46s %-10.3f %-12.1f %-10.1f %-10.0f\n', shorten(res), ...
                    entry.openTime, 1e3*entry.queryMedian, 1e3*entry.queryP90, ...
                    1e3*entry.armMedian);
        catch err
            fprintf('%-46s  FAILED - %s\n', shorten(res), err.message);
        end
        results(end+1) = entry;  %#ok<AGROW>
    end

    % --- restore, whatever happened ---------------------------------------
    try
        awg.setResource(original);
    catch err
        warning('Could not restore the original resource %s: %s', original, err.message);
    end

    fprintf('\n  (query times in ms, arm times in ms)\n');
    good = results([results.ok]);
    if length(good) > 1
        [~, iFast] = min([good.armMedian]);
        [~, iSlow] = max([good.armMedian]);
        fprintf('  fastest arm: %s at %.0f ms\n', shorten(good(iFast).resource), 1e3*good(iFast).armMedian);
        fprintf('  slowest arm: %s at %.0f ms\n', shorten(good(iSlow).resource), 1e3*good(iSlow).armMedian);
        spreadArm   = good(iSlow).armMedian   / max(good(iFast).armMedian, eps);
        spreadQuery = good(iSlow).queryMedian / max(good(iFast).queryMedian, eps);
        fprintf('  arm spread %.2fx, query spread %.2fx\n', spreadArm, spreadQuery);
        if spreadArm > 1.3
            fprintf('  => the arm time TRACKS transport latency. Switching resource is worth it:\n');
            fprintf('     put awg.setResource(''%s'') in the speed-knobs cell.\n', good(iFast).resource);
        else
            fprintf('  => the arm time barely moves with transport, so it really is firmware.\n');
            fprintf('     Only eliminating the per-point re-arm will help.\n');
        end
    end
    fprintf('\n');
end

function s = shorten(res)
    s = res;
    if length(s) > 44
        s = ['...' s(end-40:end)];
    end
end

function p = prctileSimple(x, pct)
    % Avoids depending on the Statistics Toolbox
    x = sort(x(:));
    if isempty(x); p = NaN; return; end
    idx = max(1, min(numel(x), ceil(pct/100 * numel(x))));
    p = x(idx);
end
