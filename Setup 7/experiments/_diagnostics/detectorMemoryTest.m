function out = detectorMemoryTest(opts)
%DETECTORMEMORYTEST Does the photodiode + digitizer remember the previous window?
%
%   out = detectorMemoryTest                         % ESR defaults: 125 us on, 10 us gap
%   out = detectorMemoryTest(struct('onTime', 125, 'gap', 10))
%   out = detectorMemoryTest(struct('replay', out))  % re-analyse a saved result, no hardware
%
% NO MICROWAVE. Only the green laser is switched; the AWG is not used.
%
% ===================================================================
% THE QUESTION
% ===================================================================
%
% In CW ESR the "without MW" (reference) window shows a small copy of the ESR dip,
% about 0.11-0.14 of the signal-window dip. One candidate is the DETECTOR: if the
% photodiode / amplifier / digitizer chain has a slow tail, the reference window,
% which starts only delayBetweenDetections (10 us) after the signal window, still
% carries part of the signal window's lower level.
%
% ===================================================================
% WHY A LASER STEP AND NOT MICROWAVE
% ===================================================================
%
% After the laser switches off, the real light is gone within ~30 ns (level 1
% measured tau = 27 ns on both laser edges). NV fluorescence lives ~12 ns and there
% is no NV physics in the dark that emits in the detection band. So whatever the
% digitizer still reads 1 us, 10 us, 100 us after the laser went off IS the
% detector's memory -- nothing else can produce it.
%
% A microwave test cannot separate the two: a MW-induced change also carries SPIN
% memory (NVs that have not repolarised yet), which looks exactly like detector
% memory in the reference window. The laser step removes the spins from the question.
% It is also a 100% step instead of a ~5% one, so it is ~20x more sensitive.
%
% The level-1 laser-edge calibration had the laser on for only 3 us, which is too
% short to charge up a slow component (tens of us). Here the laser stays on for the
% ESR window length, exactly as long as the signal window's lower level lasts.
%
% ===================================================================
% THE NUMBER THAT ANSWERS IT
% ===================================================================
%
% r(t) = (reading - dark) / (bright - dark), t = time after the laser went off.
% For a linear detector, the reference-window deficit caused by detector memory is
%
%     M = mean of r(t) over [gap, gap + window]
%
% times the signal-window deficit. So M compares directly with the observed
% (reference dip) / (signal dip), 0.11-0.14:
%
%     M << 0.1   the detector is not the cause
%     M ~ 0.1    the detector explains the reference dip
%
% If the detector is clean, the leading suspect is spin memory: calibDetectionWindow
% found repolarisation is NOT a single exponential (tau_pol = 12.2 us, chi2 = 47,
% still 1% left at 56.5 us), so slowly repolarising NVs (e.g. in the weak edges of
% the laser spot) can carry the depolarisation into a reference window 10 us later.
%
% ===================================================================
% HOW IT RUNS
% ===================================================================
%
% ExpAWGDigitizerCalibrateDelays (the level-1 experiment) with only 'greenLaser':
%   [darkTime dark, detector starts] [onTime laser on] [darkTime dark] [5 us dark]
% and the digitizer records the whole thing as a trace at 100 MS/s.
% Run it with the SAME laser power and detector input configuration as the ESR.
%
% A NOTE ON THE TRACE. Experiment.processData normalises photodiode data as
% (odd column)/(even column), which is meant for signal/reference pairs but is also
% applied to full traces, so it divides NEIGHBOURING SAMPLES. The even samples come
% out as plain averages; the odd ones are ratio-scaled and biased in the dark where
% the level is small. This analysis therefore uses the even samples only (50 MS/s
% effective, 20 ns per point -- far finer than anything that matters here).
%
% OPTIONS (struct fields, all optional)
%   onTime     us laser on. Default 125 (the ESR window). Use the ESR window you ran.
%   darkTime   us dark before and after. Default 300. Must be >= gap + window.
%   gap        us, the ESR delayBetweenDetections to evaluate. Default 10.
%   window     us, the ESR detection window to evaluate. Default = onTime.
%   repeats    default 100
%   averages   default 20. One average moves ~7 M samples per channel (725 us at
%              100 MS/s x 100 repeats); if the digitizer refuses, lower repeats or darkTime.
%   esrRatio   observed (ref dip)/(signal dip) to compare with. Default [0.11 0.14].
%   replay     a previous 'out' (or struct with fields t, Y): analyse without hardware
%   plot       default true
%
% OUTPUT out.M, out.Merr, out.checkpoints (table of r at fixed times), out.verdict,
%   plus the raw trace (out.t, out.Y) so it can be replayed or saved.

    if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
    opts = withDefaults(opts, struct('onTime', 125, 'darkTime', 300, 'gap', 10, ...
        'window', [], 'repeats', 100, 'averages', 20, 'esrRatio', [0.11 0.14], ...
        'replay', [], 'plot', true), 'detectorMemoryTest');
    windowGiven = ~isempty(opts.window);
    if ~isempty(opts.replay) && isfield(opts.replay, 'opts')
        % Analyse with the timing the data was TAKEN with.
        opts.onTime   = opts.replay.opts.onTime;
        opts.darkTime = opts.replay.opts.darkTime;
    end
    if ~windowGiven; opts.window = opts.onTime; end
    % processData pairs neighbouring samples, so the trace needs an EVEN number of
    % samples at 100 MS/s: keep every duration on a 0.02 us grid.
    opts.onTime   = round(opts.onTime * 50) / 50;
    opts.darkTime = round(opts.darkTime * 50) / 50;
    if opts.gap + opts.window > opts.darkTime - 5
        error(['detectorMemoryTest: darkTime (%g us) must exceed gap + window (%g us) plus a ' ...
               'margin, or the reference-window equivalent is not inside the recorded trace.'], ...
              opts.darkTime, opts.gap + opts.window);
    end

    % ---- acquire (or replay) ------------------------------------------------
    if isempty(opts.replay)
        [t, Y, info] = acquire(opts);
    else
        t = opts.replay.t(:);
        Y = opts.replay.Y;
        info = struct('source', 'replay');
    end

    out = analyse(t, Y, opts);
    out.t = t; out.Y = Y; out.opts = opts; out.info = info;

    report(out, opts);
    if opts.plot; plotIt(out, opts); end
end

% =========================================================================

function [t, Y, info] = acquire(opts)
    c = ExpAWGDigitizerCalibrateDelays;
    c.listOfChannels           = {'greenLaser'};     % laser only: no MW, no AWG trigger
    c.channelDelay             = opts.onTime;
    c.detectionStartStopOffset = opts.darkTime;
    c.repeats                  = opts.repeats;
    c.averages                 = opts.averages;
    c.isTracking               = 0;

    spcm = getObjByName(Spcm.NAME);
    info = struct('source', 'hardware');
    try; info.inputConfiguration = spcm.inputConfiguration; catch; info.inputConfiguration = '?'; end
    fprintf(['\n=== detectorMemoryTest ===\nlaser on %g us, dark %g us each side, %d repeats x %d ' ...
             'averages, detector input: %s\n'], opts.onTime, opts.darkTime, opts.repeats, ...
            opts.averages, char(info.inputConfiguration));
    fprintf('Use the same laser power and input configuration as the ESR you are explaining.\n');

    c.run;

    nAvg = c.currIter;
    if nAvg < 1
        error('detectorMemoryTest: no completed averages to analyse.');
    end
    t = c.time(:);
    Y = reshape(c.signal(1, :, 1:nAvg), [], nAvg);      % samples x averages
    if size(Y, 1) ~= numel(t)
        % The digitizer can add samples for its channel delay. The edges are found
        % from the data anyway, so rebuild the axis from the sample count.
        fprintf('NOTE trace has %d samples, time axis %d: rebuilding the axis at %g MS/s.\n', ...
                size(Y,1), numel(t), c.sampleRate);
        t = (0:size(Y,1)-1)' / c.sampleRate;
    end
end

% =========================================================================

function out = analyse(t, Y, opts)
    % Even samples only -- see the header.
    keep = 2:2:numel(t);
    t = t(keep);
    Y = Y(keep, :);
    y = mean(Y, 2);

    onCmd  = opts.darkTime;
    offCmd = opts.darkTime + opts.onTime;

    % Rough levels from the commanded timing, then the real edges from the data.
    dark0  = median(y(t > 2 & t < onCmd - 2));
    bright = median(y(t > offCmd - min(30, opts.onTime/3) & t < offCmd - 0.5));
    if bright - dark0 <= 0
        error('detectorMemoryTest: no laser step in the trace (bright %.3g <= dark %.3g).', bright, dark0);
    end
    half = (dark0 + bright) / 2;
    iOff = find(t > offCmd - 1 & y < half, 1);          % first sample below half after the off command
    iOn  = find(t > onCmd - 1 & y > half, 1);
    if isempty(iOff) || isempty(iOn)
        error('detectorMemoryTest: could not find the laser edges near the commanded times.');
    end
    tOff = t(iOff); tOn = t(iOn);

    % Final levels. dark: the pre-pulse dark, more than (darkTime + 5 + darkTime) after
    % the previous pulse. darkLate: the end of the post-pulse dark, as a consistency check.
    preWin  = t > 2 & t < tOn - 2;
    lateWin = t > t(end) - min(50, opts.darkTime/4);
    brWin   = t > tOff - min(30, opts.onTime/3) & t < tOff - 0.5;
    dark     = mean(y(preWin));
    darkLate = mean(y(lateWin));
    bright   = mean(y(brWin));
    step     = bright - dark;

    after = t >= tOff;
    tau = t(after) - tOff;
    r   = (y(after) - dark) / step;

    % The ESR number, per average so it has an error bar.
    wESR = tau >= opts.gap & tau <= opts.gap + opts.window;
    nA = size(Y, 2);
    Mk = zeros(1, nA);
    for k = 1:nA
        dk = mean(Y(preWin, k));
        sk = mean(Y(brWin, k)) - dk;
        rk = (Y(after, k) - dk) / sk;
        Mk(k) = mean(rk(wESR));
    end
    M = mean(r(wESR));
    if nA > 1; Merr = std(Mk) / sqrt(nA); else; Merr = NaN; end

    % Noise of a single (averaged) point in the dark, as a scale for "zero". From
    % point-to-point differences, so a slow slope in the dark does not inflate it.
    noisePt = std(diff(y(preWin))) / sqrt(2) / step;

    % Is the dark settled? If a slow tail from the PREVIOUS pulse is still decaying
    % through the pre-pulse dark, the baseline is too high and M comes out too low.
    pre = find(preWin);
    h = floor(numel(pre)/2);
    preSlope = (mean(y(pre(h+1:end))) - mean(y(pre(1:h)))) / step;

    % r at fixed times after the edge, averaged over a +/-20% log bin.
    cpT = [0.1 0.3 1 3 10 30 100 300];
    cpT = cpT(cpT <= tau(end) * 0.8);
    cpR = nan(size(cpT)); cpE = nan(size(cpT));
    for j = 1:numel(cpT)
        b = tau >= cpT(j)/1.2 & tau <= cpT(j)*1.2;
        if any(b)
            cpR(j) = mean(r(b));
            rkb = zeros(1, nA);
            for k = 1:nA
                dk = mean(Y(preWin, k)); sk = mean(Y(brWin, k)) - dk;
                yk = Y(after, k);
                rkb(k) = mean((yk(b) - dk) / sk);
            end
            if nA > 1; cpE(j) = std(rkb) / sqrt(nA); end
        end
    end

    out = struct('M', M, 'Merr', Merr, 'Mk', Mk, 'dark', dark, 'darkLate', darkLate, ...
                 'bright', bright, 'step', step, 'tOff', tOff, 'tOn', tOn, ...
                 'tauAfter', tau, 'rAfter', r, 'noisePerPoint', noisePt, ...
                 'checkpoints', struct('t_us', cpT, 'r', cpR, 'err', cpE), ...
                 'lateResidual', (darkLate - dark) / step, 'preDarkDrift', preSlope);

    % Verdict against the observed reference/signal dip ratio.
    need = min(opts.esrRatio);
    share = M ./ opts.esrRatio;
    sig = abs(M) / max(Merr, eps);
    if M < 0 && sig > 3
        v = sprintf(['NEGATIVE memory (M = %.2g): the reading UNDERSHOOTS after the light goes off ' ...
                     '(AC coupling / baseline restore). That would make the reference look ' ...
                     'BRIGHTER, the opposite of the dip you see. Not the cause.'], M);
    elseif abs(M) < 0.1 * need
        v = sprintf(['Detector memory is negligible: M = %.2g +/- %.1g, i.e. %.1f%%-%.1f%% of the ' ...
                     'observed reference dip. The photodiode/digitizer is NOT the cause.'], ...
                    M, Merr, 100*min(share), 100*max(share));
    elseif M > 0.5 * need
        v = sprintf(['Detector memory is LARGE: M = %.2g +/- %.1g, i.e. %.0f%%-%.0f%% of the observed ' ...
                     'reference dip. The detector explains most of it.'], ...
                    M, Merr, 100*min(share), 100*max(share));
    else
        v = sprintf(['Detector memory is PARTIAL: M = %.2g +/- %.1g, i.e. %.0f%%-%.0f%% of the ' ...
                     'observed reference dip. The rest is something else (spin memory?).'], ...
                    M, Merr, 100*min(share), 100*max(share));
    end
    out.verdict = v;
end

% =========================================================================

function report(out, opts)
    fprintf('\nlaser on edge at %.3f us, off edge at %.3f us (commanded %g and %g)\n', ...
            out.tOn, out.tOff, opts.darkTime, opts.darkTime + opts.onTime);
    fprintf('levels: dark %.4g, bright %.4g, step %.4g\n', out.dark, out.bright, out.step);
    fprintf('\nfraction of the step still read after the laser went off:\n');
    fprintf('   %8s   %12s\n', 'time', 'r');
    for j = 1:numel(out.checkpoints.t_us)
        fprintf('   %6g us   %+10.2e  +/- %.1e\n', out.checkpoints.t_us(j), ...
                out.checkpoints.r(j), out.checkpoints.err(j));
    end
    fprintf('   (one averaged point in the dark scatters by %.1e of the step)\n', out.noisePerPoint);
    fprintf(['\nESR-equivalent memory M = mean r over [%g, %g] us after the edge\n' ...
             '   M = %.3g +/- %.2g   (needed to explain the reference dip: %s)\n'], ...
            opts.gap, opts.gap + opts.window, out.M, out.Merr, mat2str(opts.esrRatio));
    % 1e-3 of the step is 1% of what it would take to explain the reference dip.
    if abs(out.lateResidual) > 1e-3 || abs(out.preDarkDrift) > 1e-3
        fprintf(['NOTE the dark has not settled: the end of the post-pulse dark sits %.2g of the\n' ...
                 '     step above the pre-pulse dark, which itself drifts by %.2g. There is a\n' ...
                 '     component slower than ~%g us, so the baseline and M are biased. Rerun with\n' ...
                 '     a longer darkTime (e.g. %g).\n'], out.lateResidual, out.preDarkDrift, opts.darkTime, ...
                4*opts.darkTime);
    end
    fprintf('\nVERDICT: %s\n\n', out.verdict);
end

% =========================================================================

function plotIt(out, opts)
    figure('Name', 'detectorMemoryTest', 'Color', 'w');
    subplot(2,1,1);
    keep = 2:2:numel(out.t);
    plot(out.t(keep), mean(out.Y(keep, :), 2), 'b-');
    hold on;
    xline(out.tOn, ':k'); xline(out.tOff, ':k');
    xlabel('time in trace (us)'); ylabel('reading');
    title(sprintf('laser on %g us, no MW', opts.onTime));
    grid on;

    subplot(2,1,2);
    tau = out.tauAfter; r = out.rAfter;
    m = tau > 0;
    pos = m & r > 0; neg = m & r <= 0;
    loglog(tau(pos), r(pos), '.', 'Color', [0 0.45 0.74], 'MarkerSize', 3); hold on;
    loglog(tau(neg), -r(neg), '.', 'Color', [0.85 0.33 0.1], 'MarkerSize', 3);
    yl = [max(min(abs(r(m & r ~= 0))), 1e-6), 2];
    patch([opts.gap, opts.gap + opts.window, opts.gap + opts.window, opts.gap], ...
          [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    loglog(xlim, min(opts.esrRatio) * [1 1], 'k--');
    loglog([opts.gap, opts.gap + opts.window], max(abs(out.M), 1e-6) * [1 1], 'r-', 'LineWidth', 2);
    ylim(yl);
    xlabel('time after laser off (us)'); ylabel('|r| = |reading - dark| / step');
    title(sprintf('M = %.2g +/- %.1g  (grey: ESR reference window; dashed: needed to explain the dip)', ...
                  out.M, out.Merr));
    legend({'r > 0', 'r < 0 (undershoot)'}, 'Location', 'southwest');
    grid on;
end

% =========================================================================

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d); given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('%s: unknown option "%s". Valid: %s.', who, bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || (isempty(o.(known{i})) && ~strcmp(known{i}, 'replay'))
            o.(known{i}) = d.(known{i});
        end
    end
end
