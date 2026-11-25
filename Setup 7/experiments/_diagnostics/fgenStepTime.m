function out = fgenStepTime(opts)
%FGENSTEPTIME How long the FGEN tone takes to settle after a frequency change.
%
%   out = fgenStepTime
%   out = fgenStepTime(struct('steps', [2e6 20e6], 'repeats', 5))
%
% Settles the one question that decides whether FGEN mode is usable for ESR. No NVs,
% no field, no laser -- just the AWG, the switch, the amplifier and the scope.
%
% ===================================================================
% THE QUESTION
% ===================================================================
%
% FGEN mode was abandoned for ESR after four rounds of debugging produced 0.039%
% contrast against AWG mode's 3.9%. The last standing explanation was that the
% AWG70000's built-in function generator accepts 2.87 GHz, stores it, reports it back
% and emits nothing useful. That explanation is now DEAD: every awgFgenTone and
% calibResonator run since has read real power at the scope across 2.5-3.3 GHz, with
% the spectral peak exactly at the drive frequency. FGEN emits.
%
% But calibResonator turned up a replacement suspect. It needs THREE SECONDS of settling
% after each frequency change -- at 0.3 s roughly half the acquisitions came back as
% flat traces with no tone at all. If that three seconds is the cost of a frequency
% write, then a 96-point ESR sweep spends five minutes settling against thirty seconds
% for the entire working AWG-mode run, FGEN is an order of magnitude SLOWER rather than
% faster, and every FGEN ESR point was measured mid-transient -- which would explain the
% missing contrast exactly.
%
% There is a strong reason to think the three seconds is not the frequency write at all.
% calibResonator's setTone does an unconditional full bring-up per point: mode assert,
% RUN:IMMEDIATE, poll the run state, outputs on. And separately, randomize is on by
% default, so consecutive frequencies there are often 800 MHz apart, where an ESR sweep
% steps 2 MHz sequentially. Either could dominate. So this measures:
%
%   'freq'  ONE command -- FGEN:CHANn:FREQuency -- and nothing else. This is what an
%           FGEN ESR sweep would actually do per point.
%   'full'  the whole setTone bring-up, which is what calibResonator does.
%
% at several step sizes spanning the ESR case (2 MHz) and the randomised resonator case
% (up to 800 MHz). If 'freq' settles in milliseconds and 'full' takes seconds, then FGEN
% ESR is viable and calibResonator was simply paying for a restart it did not need.
%
% ===================================================================
% HOW IT MEASURES
% ===================================================================
%
% Write the new frequency, then acquire as fast as the scope allows, timestamping each
% acquisition, and track the amplitude at the new frequency and where the spectral peak
% actually is. The settling curve falls straight out.
%
% RESOLUTION IS THE PROBE CADENCE, and it is measured rather than assumed: a control
% run with no frequency change at all gives both the cadence and the amplitude noise
% floor before any conclusion is drawn. Expect tens of milliseconds. That cannot resolve
% a 1 ms settle, and it does not need to -- the question is milliseconds against three
% seconds, a factor of a thousand.
%
% RECORD LENGTH IS SET BY THE SMALLEST STEP, not by convenience. A 100 ns record has
% 10 MHz frequency bins, so a 2 MHz step is invisible in it -- the DFT at the new
% frequency would read full amplitude while the synthesiser still sat at the old one,
% and the measurement would report "settled instantly" for the wrong reason. Default
% 2 us gives 0.5 MHz bins, which separates a 2 MHz step cleanly.
%
% ===================================================================
% OPTIONS (struct fields, all optional)
%   channel       AWG channel to drive (default 1)
%   scopeChannel  scope input carrying it (default = channel)
%   scopeResource VISA resource (default: the EXR/MXR/DSO from visadevlist)
%   f0            starting frequency, Hz (default 2.83e9)
%   steps         frequency steps to test, Hz (default [2e6 20e6 200e6 800e6])
%   modes         {'freq'} / {'full'} / both (default both)
%   vpp           amplitude, volts peak-to-peak (default 0.01 -- see the power budget:
%                 0.5 Vpp puts ~2.5 Vrms on the scope and ~12 W in the 20 dB pad)
%   recordTime    seconds per acquisition (default 2e-6; see RECORD LENGTH above)
%   maxWait       seconds to keep probing after the write (default 5)
%   nProbe        maximum acquisitions per trial (default 120)
%   repeats       trials per (mode, step) (default 3)
%   band          fractional amplitude band that counts as settled (default 0.05)
%   preSettle     seconds to wait before a trial's baseline (default 3)
%   range         scope full-scale volts to start from (default 0.5)
%   mwChannel     pulse-generator channel holding the switch open (default MW / MW2)
%
% Writes nothing. Restores the bench on exit, including on error.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct( ...
        'channel', 1, 'scopeChannel', [], 'scopeResource', '', ...
        'f0', 2.83e9, 'steps', [2e6 20e6 200e6 800e6], 'modes', {{'freq', 'full'}}, ...
        'vpp', 0.01, 'recordTime', 2e-6, 'maxWait', 5, 'nProbe', 120, ...
        'repeats', 3, 'band', 0.05, 'preSettle', 3, 'range', 0.5, 'mwChannel', ''));
    if isempty(opts.scopeChannel); opts.scopeChannel = opts.channel; end
    if isempty(opts.mwChannel)
        if opts.channel == 1; opts.mwChannel = 'MW'; else; opts.mwChannel = 'MW2'; end
    end
    ch  = opts.channel;
    sch = opts.scopeChannel;
    modes = opts.modes;
    if ischar(modes) || isstring(modes); modes = cellstr(modes); end

    % Record length must resolve the smallest step, or a "settled instantly" verdict
    % would be an artefact of coarse bins rather than a measurement.
    binHz = 1/opts.recordTime;
    if min(opts.steps) < 3*binHz
        need = 3/min(opts.steps);
        warning(['fgenStepTime:resolution', newline, ...
                 'recordTime %.3g s gives %.3g MHz bins, which cannot separate the %.3g MHz ' ...
                 'step. Raising recordTime to %.3g s.'], ...
                opts.recordTime, binHz/1e6, min(opts.steps)/1e6, need);
        opts.recordTime = need;
        binHz = 1/opts.recordTime;
    end

    fprintf('\n=== fgenStepTime ===\n');
    fprintf('AWG ch%d -> scope ch%d, %g Vpp, from %.4f GHz\n', ch, sch, opts.vpp, opts.f0/1e9);
    fprintf('steps: %s MHz\n', strjoin(arrayfun(@(s) sprintf('%g', s/1e6), ...
            opts.steps, 'UniformOutput', false), ', '));
    fprintf('modes: %s\n', strjoin(modes, ', '));
    fprintf('record %.3g us -> %.3g MHz bins\n', opts.recordTime*1e6, binHz/1e6);

    awg = findAwg();
    sc  = openScope(opts.scopeResource);
    pg  = getObjByName(PulseGenerator.NAME);

    out = struct('opts', opts, 'trials', [], 'control', [], 'summary', []);

    try
        pg.on(opts.mwChannel);
        setupScope(sc, sch, opts);

        % ---- bring up at f0 and fix the vertical scale --------------------------
        bringUp(awg, ch, opts.f0, opts.vpp);
        fprintf('\nsettling %g s before baseline\n', opts.preSettle);
        pause(opts.preSettle);

        [xinc, yinc, yorg, yref] = rangeAndScale(sc, sch, opts, opts.f0);

        % ---- control: probe cadence and amplitude noise, with NOTHING changing ----
        % This has to come first. Every number below is quoted against this cadence,
        % and a settle time shorter than the cadence can only be reported as an upper
        % bound -- which is a real answer, but only if the bound is measured.
        fprintf('\n--- control (no frequency change) ---\n');
        cn = min(30, opts.nProbe);
        cA = nan(1, cn); cT = nan(1, cn);
        tc = tic;
        for k = 1:cn
            cT(k) = toc(tc);
            q = probeOnce(sc, sch, opts.f0, xinc, yinc, yorg, yref);
            cA(k) = q.A;
        end
        cadence = median(diff(cT));
        A0   = median(cA, 'omitnan');
        noise = 1.4826 * median(abs(cA - A0), 'omitnan') / max(A0, eps);
        fprintf('  probe cadence      %.1f ms  <-- this is the time resolution\n', cadence*1e3);
        fprintf('  baseline amplitude %.4g mV pk\n', 1e3*A0);
        fprintf('  amplitude scatter  %.2f%% (MAD)\n', 100*noise);
        out.control = struct('t', cT, 'A', cA, 'cadence', cadence, ...
                             'A0', A0, 'relNoise', noise);
        if A0 < 1e-4
            error(['Baseline amplitude is %.4g mV -- there is no tone to watch settle. ' ...
                   'Run scopeProbe(%g) first and fix the signal path before this.'], ...
                  1e3*A0, opts.f0);
        end
        % A settled band narrower than the noise floor would never be reached.
        band = max(opts.band, 3*noise);
        if band > opts.band
            fprintf('  settled band widened %.0f%% -> %.0f%% (3x the noise floor)\n', ...
                    100*opts.band, 100*band);
        end

        % ---- the trials ---------------------------------------------------------
        T = struct('mode', {}, 'step', {}, 'rep', {}, 't', {}, 'A', {}, 'peak', {}, ...
                   'Afinal', {}, 'tFirst', {}, 'tSettle', {}, 'firstProbeSettled', {}, ...
                   'fTarget', {});
        for im = 1:numel(modes)
            for is = 1:numel(opts.steps)
                for ir = 1:opts.repeats
                    f1 = opts.f0 + opts.steps(is);
                    % Always start each trial from f0 with a full bring-up, so every
                    % trial begins from an identical, known-settled state regardless of
                    % which mode is being timed.
                    bringUp(awg, ch, opts.f0, opts.vpp);
                    pause(opts.preSettle);

                    tr = runTrial(awg, sc, ch, sch, modes{im}, f1, opts, ...
                                  xinc, yinc, yorg, yref, band);
                    tr.mode = modes{im}; tr.step = opts.steps(is); tr.rep = ir;
                    T(end+1) = tr;                                     %#ok<AGROW>
                    fprintf('  %-5s step %6.0f MHz  rep %d:  first tone %s, settled %s%s\n', ...
                            modes{im}, opts.steps(is)/1e6, ir, ...
                            tStr(tr.tFirst, cadence), tStr(tr.tSettle, cadence), ...
                            ternary(tr.firstProbeSettled, '  (already settled at the first probe)', ''));
                end
            end
        end
        out.trials = T;

        % ---- summary and verdict ------------------------------------------------
        out.summary = summarise(T, modes, opts, cadence);
        report(out, opts, cadence);
        plotIt(out, opts, cadence);

    catch err
        shutDown(awg, pg, sc);
        rethrow(err);
    end
    shutDown(awg, pg, sc);
end

% ======================================================================

function tr = runTrial(awg, sc, ch, sch, mode, f1, opts, xinc, yinc, yorg, yref, band)
% Issue the frequency change, then probe as fast as possible while timestamping.
%
% The timestamp is taken IMMEDIATELY BEFORE :DIGitize, because that is when the 2 us
% snapshot is taken. The transfer that follows adds latency to the next sample, not to
% this one, so the time axis stays honest even though the cadence is coarse.
    n = opts.nProbe;
    t = nan(1, n); A = nan(1, n); pk = nan(1, n);

    t0 = tic;
    switch lower(mode)
        case 'freq'
            % ONE command. No mode assert, no RUN, no output toggle, no readback.
            % This is exactly what an FGEN ESR sweep would issue per point.
            awg.sendCommand(sprintf('FGEN:CHAN%d:FREQuency %.10g', ch, f1));
        case 'full'
            % What calibResonator's setTone does: the entire bring-up, every point.
            bringUp(awg, ch, f1, opts.vpp);
        otherwise
            error('mode must be ''freq'' or ''full''.');
    end

    for k = 1:n
        if toc(t0) > opts.maxWait; break; end
        t(k) = toc(t0);
        q = probeOnce(sc, sch, f1, xinc, yinc, yorg, yref);
        A(k) = q.A; pk(k) = q.peakFreq;
    end
    keep = ~isnan(t);
    t = t(keep); A = A(keep); pk = pk(keep);

    % Final amplitude from the tail, so "settled" means settled relative to where it
    % ended up rather than to the pre-step baseline -- the two differ if the chain's
    % response is genuinely different at the new frequency, which for a resonator sweep
    % it certainly is.
    nTail  = min(numel(A), max(3, floor(numel(A)/4)));
    Afinal = median(A(end-nTail+1:end), 'omitnan');

    tFirst = NaN;
    i1 = find(A >= 0.5*Afinal, 1);
    if ~isempty(i1); tFirst = t(i1); end

    % Settled = this probe and EVERY later probe inside the band. Taken from the end
    % backwards, so a brief excursion after an apparent settle correctly disqualifies it.
    inBand  = abs(A - Afinal) <= band*Afinal;
    tSettle = NaN; iS = [];
    for k = 1:numel(A)
        if all(inBand(k:end)); iS = k; break; end
    end
    if ~isempty(iS); tSettle = t(iS); end

    tr = struct('mode', '', 'step', NaN, 'rep', NaN, 't', t, 'A', A, 'peak', pk, ...
                'Afinal', Afinal, 'tFirst', tFirst, 'tSettle', tSettle, ...
                'firstProbeSettled', ~isempty(iS) && iS == 1, 'fTarget', f1);
end

function q = probeOnce(sc, sch, fEval, xinc, yinc, yorg, yref)
% One acquisition, one amplitude, one spectral peak. The preamble is NOT re-read: the
% scaling was captured once when the range was fixed and nothing here changes it, and
% skipping that query is most of the reason the cadence is tens rather than hundreds of
% milliseconds. The *OPC? stays -- :DIGitize does not reliably block on this instrument,
% and the losses look like a flat trace, which is precisely the artefact under
% investigation.
    writeline(sc, sprintf(':DIGitize CHANnel%d', sch));
    try
        query1(sc, '*OPC?');
    catch
        pause(0.05);
    end
    writeline(sc, ':WAVeform:DATA?');
    raw = double(readbinblock(sc, 'int16'));
    try; readline(sc); catch; end

    v = (raw - yref)*yinc + yorg;
    n = numel(v);
    t = (0:n-1)'*xinc;
    w = 0.5*(1 - cos(2*pi*(0:n-1)'/(n-1)));
    q.A = 2*abs(sum(v(:).*w.*exp(-2i*pi*fEval*t)))/sum(w);
    V = abs(fft(v(:).*w)); V = V(1:floor(n/2));
    [~, ip] = max(V(2:end)); ip = ip + 1;
    q.peakFreq = (ip-1)/(n*xinc);
    q.vrms = sqrt(mean(v.^2));
    q.n = n;
end

function [xinc, yinc, yorg, yref] = rangeAndScale(sc, sch, opts, fEval)
% Range once, then freeze it and cache the preamble scaling. Freezing matters twice
% over: a range change mid-trial would put a step in the amplitude trace that looks
% exactly like a settling transient, and re-reading the preamble every probe would
% roughly double the cadence.
    rng = opts.range;
    for attempt = 1:5
        rng = setRange(sc, sch, rng);
        [xinc, yinc, yorg, yref] = readPreamble(sc, sch);
        q = probeOnce(sc, sch, fEval, xinc, yinc, yorg, yref);
        % Fill from VOLTS against the channel range, not from ADC counts: the EXR is a
        % 10-bit part and the count-to-full-scale mapping is not the 8-bit one.
        fill = 2*q.A/rng;
        if fill > 0.85
            if rng >= 8; break; end
            rng = min(rng*2.5, 8);
        elseif fill < 0.06 && rng > 0.008
            rng = max(rng/2.5, 0.008);
        else
            break
        end
    end
    [xinc, yinc, yorg, yref] = readPreamble(sc, sch);
    fprintf('  scope range fixed at %.4g V full scale\n', rng);
end

function [xinc, yinc, yorg, yref] = readPreamble(sc, sch)
    writeline(sc, sprintf(':WAVeform:SOURce CHANnel%d', sch));
    pre = str2double(strsplit(strtrim(query1(sc, ':WAVeform:PREamble?')), ','));
    xinc = pre(5); yinc = pre(8); yorg = pre(9); yref = pre(10);
end

function setupScope(sc, sch, opts)
    writeline(sc, ':SYSTem:HEADer OFF');
    drainErrors(sc);
    cmds = { ':ACQuire:MODE RTIMe', ':ACQuire:AVERage OFF', ...
             ':ACQuire:SRATe:ANALog MAX', ':ACQuire:POINts:ANALog AUTO', ...
             ':ACQuire:INTerpolate OFF', ':TRIGger:SWEep AUTO', ...
             sprintf(':TIMebase:SCALe %g', opts.recordTime/10), ...
             ':TIMebase:POSition 0', sprintf(':CHANnel%d:DISPlay ON', sch), ...
             ':WAVeform:FORMat WORD', ':WAVeform:BYTeorder LSBFirst' };
    refused = {};
    for i = 1:numel(cmds)
        if ~sendChecked(sc, cmds{i}); refused{end+1} = cmds{i}; end     %#ok<AGROW>
    end
    if ~isempty(refused)
        fprintf('  (instrument refused: %s)\n', strjoin(refused, ', '));
    end
    setInputFifty(sc, sch);
    writeline(sc, sprintf(':CHANnel%d:OFFSet 0', sch));
end

function bringUp(awg, ch, freqHz, vpp)
% The full FGEN bring-up, in the ONLY order that works on this instrument: configure,
% RUN, poll until the run state is non-zero, THEN outputs on. Enabling outputs before
% the run fails silently -- no SCPI error, clean error queue, nothing on the connector.
% And the run state reads 0 for about a second after a working RUN, so it is polled
% rather than asked once.
    if ~strcmpi(strtrim(awg.queryRaw('INST:MODE?')), 'FGEN')
        awg.setValue('mode', 'FGEN');
    end
    awg.sendCommandBatch({ sprintf('FGEN:CHAN%d:TYPE SINE', ch), ...
                           sprintf('FGEN:CHAN%d:FREQuency %.10g', ch, freqHz), ...
                           sprintf('FGEN:CHAN%d:AMPLITUDE %.6g', ch, vpp), ...
                           sprintf('FGEN:CHAN%d:PHASE 0', ch), ...
                           sprintf('FGEN:CHAN%d:OFFSET 0', ch) });
    awg.sendCommand('awgcontrol:run:immediate');
    t0 = tic;
    while toc(t0) < 10
        if ~strcmp(strtrim(awg.queryRaw('AWGCONTROL:RSTATE?')), '0'); break; end
        pause(0.05);
    end
    awg.workChannels = ch;
end

function s = summarise(T, modes, opts, cadence)
    s = struct('mode', {}, 'step', {}, 'tSettleMed', {}, 'tFirstMed', {}, ...
               'allFirstProbe', {});
    for im = 1:numel(modes)
        for is = 1:numel(opts.steps)
            m = strcmp({T.mode}, modes{im}) & [T.step] == opts.steps(is);
            if ~any(m); continue; end
            ts = [T(m).tSettle]; tf = [T(m).tFirst];
            s(end+1) = struct('mode', modes{im}, 'step', opts.steps(is), ...   %#ok<AGROW>
                'tSettleMed', median(ts, 'omitnan'), 'tFirstMed', median(tf, 'omitnan'), ...
                'allFirstProbe', all([T(m).firstProbeSettled]));
        end
    end
end

function report(out, opts, cadence)
    S = out.summary;
    fprintf('\n--- settling time ---\n');
    fprintf('%-6s %10s %14s %14s\n', 'mode', 'step/MHz', 'first tone', 'settled');
    for i = 1:numel(S)
        fprintf('%-6s %10.0f %14s %14s\n', S(i).mode, S(i).step/1e6, ...
                tStr(S(i).tFirstMed, cadence), tStr(S(i).tSettleMed, cadence));
    end

    fprintf('\n--- what this means for ESR ---\n');
    iF = strcmp({S.mode}, 'freq');
    if ~any(iF)
        fprintf('  (mode ''freq'' was not run, so nothing can be said about an ESR sweep)\n');
        return
    end
    SF = S(iF);
    % An ESR sweep steps a couple of MHz sequentially, so judge it on the SMALLEST step
    % measured -- and say so, rather than quietly using the best number available.
    [~, iSmall] = min([SF.step]);
    tS = SF(iSmall).tSettleMed;
    bounded = SF(iSmall).allFirstProbe;
    nPts = 96;
    if bounded
        fprintf(['  A %g MHz step was already settled at the FIRST probe, so settling is\n' ...
                 '  BELOW the %.0f ms resolution. For %d ESR points that is under %.1f s of\n' ...
                 '  settling in total, against ~30 s for the whole working AWG-mode run.\n' ...
                 '  Frequency switching is NOT what broke FGEN ESR.\n'], ...
                SF(iSmall).step/1e6, cadence*1e3, nPts, nPts*cadence);
    elseif isfinite(tS)
        fprintf(['  A %g MHz step settles in %.0f ms. For %d ESR points that is %.1f s of\n' ...
                 '  settling alone, against ~30 s for the whole working AWG-mode run.\n'], ...
                SF(iSmall).step/1e6, tS*1e3, nPts, nPts*tS);
        if nPts*tS > 30
            fprintf('  So FGEN would be SLOWER than the path it was meant to replace.\n');
        else
            fprintf('  So the speed case survives, and the missing contrast needs another cause.\n');
        end
    else
        fprintf('  The %g MHz step never settled inside %g s. That is the answer.\n', ...
                SF(iSmall).step/1e6, opts.maxWait);
    end

    % Does the step size matter? If big steps are slow and small ones fast, that alone
    % explains calibResonator's three seconds, because randomize puts 800 MHz between
    % consecutive points there while ESR steps 2 MHz.
    if numel(SF) > 1
        ts = [SF.tSettleMed]; st = [SF.step];
        [~, iB] = max(st);
        if isfinite(ts(iB)) && isfinite(tS) && ts(iB) > max(4*tS, 4*cadence)
            fprintf(['\n  Step size matters: %g MHz takes %.0f ms against %.0f ms for %g MHz.\n' ...
                     '  calibResonator randomises its frequency order, so consecutive points\n' ...
                     '  there are often ~800 MHz apart -- which would account for its 3 s on its\n' ...
                     '  own, with no implication for a sequential 2 MHz ESR sweep.\n'], ...
                    st(iB)/1e6, ts(iB)*1e3, tS*1e3, SF(iSmall).step/1e6);
        end
    end

    % And is the bring-up the real cost?
    iU = strcmp({S.mode}, 'full');
    if any(iU)
        SU = S(iU);
        [~, jSmall] = min([SU.step]);
        tU = SU(jSmall).tSettleMed;
        if isfinite(tU) && isfinite(tS) && tU > max(4*tS, 4*cadence)
            fprintf(['\n  The BRING-UP is the expensive part: %.0f ms for the full setTone\n' ...
                     '  sequence against %.0f ms for a bare frequency write at the same step.\n' ...
                     '  calibResonator can drop its settle time a long way by not restarting\n' ...
                     '  the generator every point.\n'], tU*1e3, tS*1e3);
        elseif isfinite(tU) && isfinite(tS)
            fprintf(['\n  The bring-up is not the expensive part (%.0f ms full against %.0f ms\n' ...
                     '  bare), so calibResonator''s 3 s is the frequency change itself.\n'], ...
                    tU*1e3, tS*1e3);
        end
    end
end

function plotIt(out, opts, cadence)
    T = out.trials;
    if isempty(T); return; end
    figure('Name', 'fgenStepTime', 'Color', 'w');
    modes = unique({T.mode}, 'stable');
    for im = 1:numel(modes)
        subplot(numel(modes), 1, im); hold on; box on;
        m = strcmp({T.mode}, modes{im});
        steps = unique([T(m).step]);
        co = lines(numel(steps));
        for is = 1:numel(steps)
            k = find(m & [T.step] == steps(is));
            for j = k
                y = T(j).A / max(T(j).Afinal, eps);
                plot(max(T(j).t, cadence/2)*1e3, 100*y, '.-', 'Color', co(is,:), ...
                     'MarkerSize', 8, 'HandleVisibility', ternary(j == k(1), 'on', 'off'), ...
                     'DisplayName', sprintf('%g MHz', steps(is)/1e6));
            end
        end
        yline(100, 'k:', 'HandleVisibility', 'off');
        yline(100*(1 - opts.band), 'k--', 'HandleVisibility', 'off');
        xline(cadence*1e3, 'r-', 'HandleVisibility', 'off');
        set(gca, 'XScale', 'log');
        xlabel('ms since the frequency was written  (red line = probe cadence = resolution)');
        ylabel('% of final amplitude');
        title(sprintf('mode ''%s''', modes{im}));
        legend('Location', 'southeast');
    end
end

% ======================================================================

function s = tStr(t, cadence)
    if ~isfinite(t); s = 'never'; return; end
    if t <= cadence*1.5; s = sprintf('<%.0f ms', cadence*1e3); return; end
    if t < 1; s = sprintf('%.0f ms', t*1e3); else; s = sprintf('%.2f s', t); end
end

function v = ternary(c, a, b)
    if c; v = a; else; v = b; end
end

function shutDown(awg, pg, sc)
% Microwave off, switch closed, AWG back in the mode the experiments expect. Runs even
% on the error path, because an abort would otherwise leave the amplifier driven.
    try; pg.off; catch; end
    try; awg.workChannels = []; catch; end
    try; awg.setValue('mode', 'AWG'); awg.invalidateLoaded; catch; end
    try; drainErrors(sc); catch; end
    fprintf('\n(bench restored: switch off, AWG outputs off, AWG mode)\n');
end

function awg = findAwg()
% getObjByName('AWG') does not work -- the object is named after
% FrequencyGeneratorTektronixAWG.TYPE plus its serial, e.g. 'TektronixAWG-SN-NA'. This
% mirrors ExperimentAWG's own lookup.
    awg = getObjByName('TektronixAWG');
    if ~isempty(awg); return; end
    fgs = FrequencyGenerator.getFG();
    for i = 1:numel(fgs)
        if contains(fgs{i}.name, 'TektronixAWG') || contains(fgs{i}.name, 'AWGdummy')
            awg = fgs{i}; return
        end
    end
    error('No Tektronix AWG found.');
end

function sc = openScope(resource)
    if isempty(resource)
        L = visadevlist;
        M = upper(string(L.Model));
        hit = startsWith(M, "EXR") | startsWith(M, "MXR") | contains(M, "DSO");
        if ~any(hit)
            error(['No Keysight scope found in visadevlist (looked for EXR*, MXR*, ' ...
                   '*DSO*). Pass opts.scopeResource explicitly.']);
        end
        idx = find(hit);
        isUsb = startsWith(upper(string(L.ResourceName(idx))), "USB");
        if any(isUsb); idx = idx(find(isUsb, 1)); else; idx = idx(1); end
        resource = L.ResourceName(idx);
        fprintf('scope resource: %s (%s)\n', resource, L.Model(idx));
    end
    sc = visadev(resource);
    sc.Timeout = 20;
    fprintf('scope: %s\n', strtrim(query1(sc, '*IDN?')));
end

function ok = setInputFifty(sc, c)
% 50 ohm, verified by readback. The EXR inputs are switchable and the 1 Mohm path rolls
% off a few hundred MHz up, so a GHz tone into it reads as almost nothing -- which here
% would look exactly like a tone that never settles.
    ok = false;
    attempts = { sprintf(':CHANnel%d:INPut DC50', c),      sprintf(':CHANnel%d:INPut?', c); ...
                 sprintf(':CHANnel%d:IMPedance FIFTy', c), sprintf(':CHANnel%d:IMPedance?', c) };
    for k = 1:size(attempts, 1)
        try; writeline(sc, attempts{k,1}); catch; end
        r = upper(queryOpt(sc, attempts{k,2}));
        if contains(r, "DC50") || contains(r, "DCF") || contains(r, "FIFT") || contains(r, "50")
            fprintf('  scope ch%d input: %s (50 ohm confirmed)\n', c, r);
            ok = true; return
        end
    end
    warning(['fgenStepTime:inputImpedance', newline, ...
             'Could not confirm 50 ohm on scope channel %d. On 1 Mohm the tone is rolled ' ...
             'off by the input path and every amplitude here is meaningless.'], c);
end

function actual = setRange(sc, c, want)
    try
        writeline(sc, sprintf(':CHANnel%d:RANGe %g', c, want));
        actual = str2double(query1(sc, sprintf(':CHANnel%d:RANGe?', c)));
    catch
        actual = NaN;
    end
    if ~isfinite(actual); actual = want; end
end

function ok = sendChecked(sc, cmd)
    ok = true;
    try
        writeline(sc, cmd);
    catch
        ok = false; return
    end
    e = queryOpt(sc, ':SYSTem:ERRor? STRing');
    if ~isempty(e) && ~startsWith(e, '0') && ~startsWith(e, '+0') && ...
       ~strcmp(e, '(not supported)') && ~strcmp(e, '(empty response)')
        ok = false;
    end
end

function o = withDefaults(o, d)
% Same local helper the other calib scripts carry. Note it treats an EMPTY field as
% absent, which is why 'scopeChannel', [] and 'mwChannel', '' work as "work it out".
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(o, fn{i}) || isempty(o.(fn{i})); o.(fn{i}) = d.(fn{i}); end
    end
end

function s = query1(sc, cmd)
    writeline(sc, cmd);
    s = char(readline(sc));
end

function s = queryOpt(sc, cmd)
% A query the instrument is allowed to refuse. An unsupported mnemonic returns nothing
% and the read blocks, so the timeout drops to 2 s and the connection is flushed after,
% otherwise the next legitimate read collects this query's late answer.
    old = sc.Timeout;
    s = '(not supported)';
    try
        sc.Timeout = 2;
        writeline(sc, cmd);
        s = strtrim(char(readline(sc)));
        if isempty(s); s = '(empty response)'; end
    catch
        try; flush(sc); catch; end
        drainErrors(sc);
    end
    try; sc.Timeout = old; catch; end
end

function drainErrors(sc)
    old = sc.Timeout;
    try
        sc.Timeout = 2;
        for k = 1:20
            writeline(sc, ':SYSTem:ERRor? STRing');
            r = strtrim(char(readline(sc)));
            if isempty(r) || startsWith(r, '0') || startsWith(r, '+0'); break; end
        end
    catch
    end
    try; sc.Timeout = old; catch; end
end
