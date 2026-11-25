function out = laserRiseFall(opts)
%LASERRISEFALL Rise and fall time of the green laser, read on the REFERENCE photodiode only.
%
%   out = laserRiseFall                                   % defaults below
%   out = laserRiseFall(struct('averages', 50))
%   out = laserRiseFall(struct('replay', out))            % refit a saved result, no hardware
%
% DIAGNOSTIC. Writes nothing. Leaves the laser OFF (as every experiment does) and
% leaves fixDelays as it found it.
%
% ===================================================================
% WHAT RUNS -- AND WHAT DOES NOT
% ===================================================================
%
% The PulseStreamer switches exactly two channels:
%   'greenLaser'  the gate to the AOM driver: laser on for onTime, then off
%   'detector'    the digitizer trigger, one rising edge per record
% No MW switch, no AWG trigger, no microwave. The AOM level is the static value on the
% laser panel (NI-DAQ ao0); this script does not touch it, so run it at the power you
% care about.
%
% The digitizer acquires and reads ONE channel: 'input_2' of setupInfo.json, the
% reference photodiode that samples the beam. The fluorescence photodiode (input_1)
% is not read, and nothing goes through Experiment.processData or
% spcm.readFromExperiment (whose 'referenced' mode would return input_1 divided by
% input_2, and whose photodiode normalisation divides neighbouring samples of a trace).
% The digitizer is driven directly, with zero-suppression OFF: readExperimentData's
% default turns every exact 0 V sample into NaN, which on a photodiode sitting near 0 V
% in the dark would bias the dark level. Records that were missed entirely (all zeros)
% are dropped instead.
%
% fixDelays = 0: the stored greenLaser on/off delays are NOT applied, so the edges sit
% at their true (absolute) delay after the command. The delays are printed, but the
% numbers this script is for are the rise and fall times, which do not depend on that.
%
% ===================================================================
% 1 ns SAMPLING FROM A 10 ns DIGITIZER  (equivalent-time sampling)
% ===================================================================
%
% The digitizer samples every 10 ns; an AOM edge is a few tens of ns, i.e. a handful of
% samples. So every repeat holds nShifts blocks, and in block k the laser command is
% delayed by k*shiftStep (default 1 ns, the PulseStreamer's resolution) relative to the
% digitizer trigger:
%
%   block k:  [detector, darkBefore + k ns] [greenLaser, onTime] [dark, darkAfter - k ns + deadTime]
%
% Each block is one digitizer record. Plotting record k against t - k ns puts the ten
% sample grids between each other: one point per ns. The block length is the same for
% every k, so the laser duty cycle is the same as without the shifts.
%
% ===================================================================
% TRIGGER ALIGNMENT (triggerJitter)
% ===================================================================
%
% The digitizer runs on its own clock (externalRefClock = 0), not the PulseStreamer's,
% so each record starts on the first digitizer sample after the trigger: a random 0-10 ns
% offset per record. Averaged, the recorded edge is the laser edge smeared by a 10 ns
% wide box. Every experiment on this setup sees that same smear. The fit includes the
% box, and reports both:
%   laser        the edge with the box taken out  (the headline numbers)
%   as recorded  with the box, as the digitizer sees it
% For a 50 ns edge the difference is < 1 ns; for a 10 ns edge it is a few ns.
% triggerJitter = 0 turns the box off (e.g. if the clocks are ever locked).
%
% The times always include the reference photodiode, its amplifier and the digitizer.
% A photodiode with bandwidth B has 10-90% ~ 0.35/B on its own; if that is not much
% shorter than the result, the photodiode is part of what you are measuring.
%
% ===================================================================
% THE FIT -- one edge at a time
% ===================================================================
%
% Each edge is fitted in its own window (edgeWindow, default 1 us before to 2 us after
% its 50% point):
%
%   rise:  y = lo + (hi - lo) *      G(t)
%   fall:  y = lo + (hi - lo) * (1 - G(t))
%
%   G = Gaussian-smoothed Heaviside x exponential (as level 1, calibChannelDelays):
%       sigma -> 0 is a pure exponential (RC-like), tau -> 0 a pure erf (the acoustic
%       wave crossing a Gaussian beam in the AOM). 5 parameters.
%   'double' adds a second, slower exponential with its own weight w (w < 0 is an
%       overshoot that settles back): G = (1-w) S(tau) + w S(rho*tau). Kept by 'auto'
%       only when it is really needed (F > 20, as in fitPinholeFocus).
%
% Weighted by the scatter across averages. Levenberg-Marquardt written here (base MATLAB,
% no toolbox), several starts. From the fitted G:
%   t10, t50, t90   when the edge is 10/50/90% COMPLETE, in ns after the command. For the
%                   fall, 90% complete = 10% of the light left (the Setup 7 delay
%                   convention: 90% on a turn-on, 10% left on a turn-off).
%   rise time = t90 - t10 of the rise (10-90%), fall time = t90 - t10 of the fall (90-10%).
% Errors from the covariance, inflated by sqrt(chi2/dof) when that is above 1.
%
% MODEL-FREE CHECK. The same crossings read straight off the averaged data (5 ns running
% mean, levels from the window ends), with the error from the spread over averages. If
% the fit and the data disagree, the fit does not describe the edge, and the script says
% so -- trust the data crossings then, and look at the residuals.
%
% ===================================================================
% OPTIONS (struct fields, all optional)
% ===================================================================
%   onTime         us laser on (default 10). Long enough for the light to settle.
%   darkBefore     us dark before the on command, inside the record (default 3)
%   darkAfter      us dark after the off command, inside the record (default 5)
%   deadTime       us between records, for the digitizer to re-arm (default 5)
%   repeats        repeats per average (default 200); each repeat has nShifts records
%   averages       default 20
%   shiftStep      us between the laser shifts (default 0.001 = 1 ns)
%   nShifts        number of shifts (default []: fill one digitizer sample, 10).
%                  1 = no equivalent-time sampling.
%   edgeWindow     [before after] the 50% point, us, for the fit (default [1 2])
%   triggerJitter  'auto' (default: one digitizer sample if its clock is internal), or a
%                  number in us, 0 = none
%   model          'auto' (default), 'single' or 'double'
%   plot           default true
%   replay         a previous 'out': refit without hardware. Its timing is used;
%                  edgeWindow, triggerJitter, model and plot can be changed.
%
% OUTPUT
%   out.rise, out.fall   t10/t50/t90 (ns after the command), time (10-90 / 90-10, ns),
%                        the same 'asRecorded', errors, fit parameters, chi2/dof, model,
%                        and the model-free values (out.rise.data)
%   out.levels           dark, on, step (V), plateau change during the pulse, dark
%                        recovery, clipping
%   out.tau, out.Y       equivalent-time axis (us after the trigger, unshifted) and one
%                        averaged trace per average (V) -- enough to replay
%   out.opts, out.info
%
% Save it with  save('laserRiseFall_<date>.mat', 'out')  if you want to keep it.

    if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
    opts = withDefaults(opts, struct( ...
        'onTime', 10, 'darkBefore', 3, 'darkAfter', 5, 'deadTime', 5, ...
        'repeats', 200, 'averages', 20, 'shiftStep', 0.001, 'nShifts', [], ...
        'edgeWindow', [1 2], 'triggerJitter', 'auto', 'model', 'auto', ...
        'plot', true, 'replay', []), 'laserRiseFall');

    replay = ~isempty(opts.replay);
    if replay
        R = opts.replay;
        if ~isstruct(R) || ~all(isfield(R, {'tau', 'Y', 'opts', 'info'}))
            error('laserRiseFall: replay must be an earlier output (fields tau, Y, opts, info).');
        end
        % Analyse with the timing the data was TAKEN with.
        keep = {'onTime', 'darkBefore', 'darkAfter', 'deadTime', 'repeats', 'averages', ...
                'shiftStep', 'nShifts'};
        for i = 1:numel(keep)
            opts.(keep{i}) = R.opts.(keep{i});
        end
    end

    opts = checkTiming(opts);
    if ~any(strcmpi(opts.model, {'auto', 'single', 'double'}))
        error('laserRiseFall: model must be ''auto'', ''single'' or ''double''.');
    end
    opts.model = lower(opts.model);

    % ---- acquire (or replay) ------------------------------------------------
    if replay
        tau  = R.tau(:);
        Y    = R.Y;
        info = R.info;
        info.source = 'replay';
        fprintf('\n=== laserRiseFall: REPLAY of %d averages, no hardware touched ===\n', size(Y, 2));
    else
        [tau, Y, info] = acquire(opts);
    end
    if ischar(opts.triggerJitter)
        if ~strcmpi(opts.triggerJitter, 'auto')
            error('laserRiseFall: triggerJitter must be ''auto'' or a number of us.');
        end
        J = info.triggerJitterAuto;
    else
        J = opts.triggerJitter;
    end
    info.triggerJitter = J;

    % ---- analyse ------------------------------------------------------------
    out = analyse(tau, Y, opts, info, J);
    out.tau = tau;  out.Y = Y;  out.opts = opts;  out.info = info;

    report(out);
    if opts.plot
        plotIt(out);
    end
end

% =========================================================================
% TIMING
% =========================================================================

function opts = checkTiming(opts)
    % Whole digitizer samples (10 ns) for everything the record contains, whole
    % PulseStreamer ns for the shifts.
    r10 = @(x) round(x * 100) / 100;
    opts.onTime     = r10(opts.onTime);
    opts.darkBefore = r10(opts.darkBefore);
    opts.darkAfter  = r10(opts.darkAfter);
    opts.deadTime   = r10(opts.deadTime);
    opts.shiftStep  = max(1, round(opts.shiftStep * 1e3)) / 1e3;
    if isempty(opts.nShifts)
        opts.nShifts = max(1, round(0.01 / opts.shiftStep));     % fill one 10 ns sample
    end
    opts.nShifts = round(opts.nShifts);
    w = opts.edgeWindow;
    if ~(isnumeric(w) && numel(w) == 2 && all(w > 0))
        error('laserRiseFall: edgeWindow must be [before after] in us, both > 0.');
    end
    maxShift = (opts.nShifts - 1) * opts.shiftStep;
    % The laser edges arrive ~0.5 us after the command on Setup 7; leave room for that
    % plus the fit window on each side.
    if opts.darkBefore < max(1, w(1))
        error('laserRiseFall: darkBefore (%g us) must be at least max(1, edgeWindow(1)) = %g us.', ...
              opts.darkBefore, max(1, w(1)));
    end
    if opts.darkAfter < w(2) + 1 + maxShift
        error('laserRiseFall: darkAfter (%g us) must be at least edgeWindow(2) + 1 us = %g us.', ...
              opts.darkAfter, w(2) + 1 + maxShift);
    end
    if opts.onTime < sum(w) + 1
        error('laserRiseFall: onTime (%g us) must be at least sum(edgeWindow) + 1 us = %g us.', ...
              opts.onTime, sum(w) + 1);
    end
    if opts.deadTime < 2
        error('laserRiseFall: deadTime must be at least 2 us (digitizer re-arm).');
    end
    if maxShift >= opts.darkAfter
        error('laserRiseFall: the largest shift (%g us) must be shorter than darkAfter.', maxShift);
    end
end

% =========================================================================
% ACQUISITION
% =========================================================================

function [tau, Y, info] = acquire(opts)
    pg   = getObjByName(PulseGenerator.NAME);
    spcm = getObjByName(Spcm.NAME);
    digi = getObjByName(Digitizer.NAME);
    if isempty(pg) || isempty(spcm) || isempty(digi)
        error(['laserRiseFall: pulse generator, photodiode or digitizer object not found. ' ...
               'Start the setup (the usual init) first.']);
    end
    try
        names = spcm.digitizerChannelNames;
    catch
        error('laserRiseFall: the detector (%s) is not the digitizer photodiode.', class(spcm));
    end
    refChan = find(strcmp(names, 'input_2'), 1);
    if isempty(refChan)
        error(['laserRiseFall: no ''input_2'' (reference photodiode) in the digitizer channel ' ...
               'names of setupInfo.json: %s'], strjoin(names, ', '));
    end

    fs   = digi.MAX_SAMPLE_RATE;            % samples per us
    K    = opts.nShifts;
    s    = (0:K-1) * opts.shiftStep;        % us, laser shift per block
    tPre = opts.darkBefore;  tOn = opts.onTime;  tPost = opts.darkAfter;  tDead = opts.deadTime;
    segDur = tPre + tOn + tPost;
    nSamp  = round(segDur * fs);
    nSeg   = opts.repeats * K;
    ns = @(x) round(x * 1e3) / 1e3;         % whole PulseStreamer ns

    S = Sequence;
    for k = 1:K
        S.addEvent(ns(tPre + s(k)),         'detector');      % rising edge = record start
        S.addEvent(ns(tOn),                 'greenLaser');
        S.addEvent(ns(tPost - s(k) + tDead), '');
    end

    jitterAuto = 0;
    try
        if ~digi.externalRefClock; jitterAuto = 1 / fs; end
    catch
        jitterAuto = 1 / fs;
    end
    maxV = NaN; offV = 0;
    try maxV = digi.maxVoltage(refChan); catch; end
    try offV = digi.voltageOffset(refChan); catch; end

    fprintf('\n=== laserRiseFall: green laser on the REFERENCE photodiode only ===\n');
    fprintf(['PulseStreamer switches greenLaser and detector (digitizer trigger) only; digitizer ' ...
             'reads channel %d (input_2) only.\n'], refChan);
    fprintf(['laser on %g us, %g us dark before and %g after inside each record, %g us between ' ...
             'records\n'], tOn, tPre, tPost, tDead);
    if K > 1
        fprintf('%d laser shifts of %g ns: %g ns effective sampling (digitizer %g ns)\n', ...
                K, 1e3*opts.shiftStep, 1e3/fs/K, 1e3/fs);
    end
    fprintf('%d repeats x %d records per average, %d averages; fixDelays = 0 (absolute delays)\n', ...
            opts.repeats, K, opts.averages);
    fprintf('Uses the AOM level currently set on the laser panel.\n');

    oldFix = pg.fixDelays;
    cleanup = onCleanup(@() restoreState(pg, spcm, digi, oldFix)); %#ok<NASGU>

    pg.setSequence(S);
    pg.repeats   = opts.repeats;
    pg.fixDelays = 0;
    spcm.setSPCMEnable(true);               % photodiode supply gate, as every experiment does
    digi.prepareAcquisition(refChan, nSeg, segDur, 1, fs);

    t = (0:nSamp-1)' / fs;                  % us after the trigger
    tauGrid = t - s;                        % nSamp x K, equivalent time
    [tau, order] = sort(tauGrid(:));

    Y = nan(numel(tau), opts.averages);
    nMissed = 0;  nClip = 0;  nTotal = 0;  vMax = -Inf;  vMin = Inf;
    nAvg = 0;
    for a = 1:opts.averages
        V = [];
        for trial = 1:3
            try
                digi.startExperiment;
                pg.run;
                V = digi.readExperimentData(refChan, true, false);
                break
            catch err
                warning('laserRiseFall:trial', 'average %d, trial %d failed: %s', a, trial, err.message);
                try; digi.stopAcquisition; catch; end
                V = [];
            end
        end
        if isempty(V)
            if nAvg == 0
                error('laserRiseFall: the first average failed three times (see warnings above).');
            end
            fprintf('Stopping after %d averages: average %d failed three times.\n', nAvg, a);
            break
        end
        V = reshape(V, size(V, 2), size(V, 3));            % records x samples
        if size(V, 1) ~= nSeg
            error('laserRiseFall: got %d records, expected %d.', size(V, 1), nSeg);
        end
        if size(V, 2) < nSamp
            error('laserRiseFall: records have %d samples, expected %d.', size(V, 2), nSamp);
        end
        V = V(:, 1:nSamp);

        missed = all(V == 0, 2) | any(isnan(V), 2);
        nMissed = nMissed + nnz(missed);
        V(missed, :) = NaN;

        good = V(~missed, :);
        vMax = max(vMax, max(good(:)));  vMin = min(vMin, min(good(:)));
        if isfinite(maxV)
            nClip = nClip + nnz(good >= offV + 0.995*maxV | good <= offV - 0.995*maxV);
        end
        nTotal = nTotal + numel(good);

        % Records come in play order: (repeat 1: shift 1..K), (repeat 2: ...), ...
        V3 = reshape(V', nSamp, K, opts.repeats);
        V3 = V3(:, :, 2:end);       % first repeat: the laser had been off for a long time
        Yk = nmean(V3, 3);          % nSamp x K
        Y(:, a) = Yk(order);
        nAvg = a;
        fprintf('.');
        if mod(a, 50) == 0; fprintf('\n'); end
    end
    fprintf('\n');
    Y = Y(:, 1:nAvg);

    info = struct('source', 'hardware', 'refChannel', refChan, 'sampleRate', fs, ...
                  'triggerJitterAuto', jitterAuto, 'tCmdOn', tPre, 'tCmdOff', tPre + tOn, ...
                  'maxVoltage', maxV, 'voltageOffset', offV, 'vMax', vMax, 'vMin', vMin, ...
                  'clippedFraction', nClip / max(nTotal, 1), 'missedRecords', nMissed, ...
                  'nAverages', nAvg, 'stored', storedLaserDelays(), 'when', datestr(now, 31));
    if nMissed > 0
        fprintf('NOTE %d of %d records came back empty (missed trigger) and were dropped.\n', ...
                nMissed, nSeg * nAvg);
    end
end

function restoreState(pg, spcm, digi, oldFix)
    try; digi.stopAcquisition; catch; end
    try; pg.fixDelays = oldFix; catch; end
    try; spcm.setSPCMEnable(false); catch; end
end

function s = storedLaserDelays()
    % greenLaser on/off delays in setupInfo.json, for orientation only.
    s = struct('on', NaN, 'off', NaN);
    try
        j = JsonInfoReader.getJson();
        L = j.lasers;
        if isstruct(L); L = num2cell(L); end
        for i = 1:numel(L)
            if isfield(L{i}, 'switch') && strcmp(L{i}.switch.switchChannelName, 'greenLaser')
                s.on  = L{i}.switch.onDelay;
                s.off = L{i}.switch.offDelay;
            end
        end
    catch
    end
end

% =========================================================================
% ANALYSIS
% =========================================================================

function out = analyse(tau, Y, opts, info, J)
    good = any(isfinite(Y), 1);
    Y = Y(:, good);
    A = size(Y, 2);
    if A == 0
        error('laserRiseFall: no data.');
    end
    y = nmean(Y, 2);
    if A >= 3
        % Scatter across averages, pooled over +/-10 ns: a variance from a handful of
        % averages is itself noisy, and weights 1/s^2 from it would favour the points that
        % happen to scatter little (chi2 then comes out ~(A-1)/(A-3) instead of 1).
        v = nstd(Y, 2).^2;
        v = smoothRun(v, max(1, round(0.021 / median(diff(tau)))));
        e = sqrt(v / A);
    else
        e = [];
    end
    % Floor the error bars: a point with ~0 scatter (or A < 3) must not dominate the fit.
    noise = robustNoise(y);
    if isempty(e)
        e = noise * ones(size(y));
    else
        e = max(e, 0.3 * noise);
        e(~isfinite(e)) = noise;
    end

    tOn = info.tCmdOn;  tOff = info.tCmdOff;
    darkPre = nmedian(y(tau > 0.1 & tau < tOn - 0.05));
    onLevel = nmedian(y(tau > tOn + 0.5*opts.onTime & tau < tOff));
    step = onLevel - darkPre;
    if ~(abs(step) > 5 * noise)
        error(['laserRiseFall: no laser step on the reference photodiode (dark %.4g V, on %.4g V, ' ...
               'noise %.2g V per point).\nCheck that the reference photodiode is on digitizer ' ...
               'input 2, that the laser is on and the AOM level is not at its minimum.'], ...
              darkPre, onLevel, noise);
    end
    sgn = sign(step);
    ys  = smoothRun(y, max(1, round(0.005 / median(diff(tau)))));
    half = (darkPre + onLevel) / 2;
    t50on  = crossTime(tau, sgn * ys, sgn * half, tOn - 0.05, +1);
    t50off = crossTime(tau, sgn * ys, sgn * half, tOff - 0.05, -1);
    if ~isfinite(t50on) || ~isfinite(t50off)
        error('laserRiseFall: could not find both laser edges after their commands.');
    end

    w = opts.edgeWindow;
    rise = fitEdge(tau, y, e, Y, 'rise', tOn,  t50on,  [max(tau(1), t50on - w(1)),  min(t50on + w(2), tOff - 0.05)], J, opts.model);
    fall = fitEdge(tau, y, e, Y, 'fall', tOff, t50off, [max(t50off - w(1), t50on + 0.2), min(t50off + w(2), tau(end))], J, opts.model);

    % Levels and checks.
    lv.dark   = darkPre;
    lv.on     = onLevel;
    lv.step   = step;
    lv.noisePerPoint = noise;
    early = tau > t50on + w(2) & tau < t50on + w(2) + 1;
    late  = tau > tOff - 1 & tau < tOff - 0.02;
    lv.plateauChange = (nmedian(y(late)) - nmedian(y(early))) / step;   % fraction of the step
    lv.darkEnd = nmedian(y(tau > tau(end) - min(1, opts.darkAfter/3)));
    lv.darkRecovery = (lv.darkEnd - darkPre) / step;                     % fraction of the step
    lv.clippedFraction = getf(info, 'clippedFraction', 0);
    lv.rangeUsed = NaN;
    if isfinite(getf(info, 'maxVoltage', NaN))
        lv.rangeUsed = max(abs([getf(info, 'vMax', onLevel), getf(info, 'vMin', darkPre)] ...
                           - getf(info, 'voltageOffset', 0))) / info.maxVoltage;
    end

    out = struct('rise', rise, 'fall', fall, 'levels', lv, 'nAverages', A, 'triggerJitter', J);
end

% -------------------------------------------------------------------------

function E = fitEdge(tau, y, e, Y, kind, tCmd, t50, win, J, modelChoice)
    isRise = strcmp(kind, 'rise');
    m = tau >= win(1) & tau <= win(2) & isfinite(y);
    t = tau(m);  yy = y(m);  ee = e(m);
    Yw = Y(m, :);

    % Levels from the two ends of the window.
    span = win(2) - win(1);
    pre  = t < t50 - max(0.15, 0.25 * (t50 - win(1)));
    post = t > t50 + 0.5 * (win(2) - t50);
    yPre = nmedian(yy(pre));  yPost = nmedian(yy(post));
    if isRise; lo = yPre; hi = yPost; else; hi = yPre; lo = yPost; end
    E = struct('kind', kind, 'commanded', tCmd, 'window', win);

    % ---- model-free crossings (data) ----------------------------------------
    sw = max(1, round(0.005 / median(diff(t))));
    E.data = dataCrossings(t, Yw, yy, sw, lo, hi, isRise, tCmd, pre, post);

    % ---- fit ----------------------------------------------------------------
    wr = E.data.time / 1e3;                         % us, rough 10-90 width
    if ~(isfinite(wr) && wr > 0); wr = 0.05; end
    wr = min(max(wr, 0.002), 0.5 * span);
    d = hi - lo;
    hStep = [1e-4*abs(d), 1e-4*abs(d), 1e-5, 1e-5, 1e-5, 1e-3, 1e-2];

    % Single model: p = [lo hi t0 tau sg], several starts on (sg, tau).
    starts = [wr/2.56, wr/50; wr/50, wr/2.2; wr/4, wr/4.4; wr/2, wr/8];
    best = [];
    for i = 1:size(starts, 1)
        sg0 = starts(i, 1);  ta0 = starts(i, 2);
        t0 = t50 - offsetAt(0.5, [0 1 0 ta0 sg0], 0);
        p0 = [lo hi t0 ta0 sg0];
        clampFn = @(p) clampP(p, win);
        resid = @(p) (edgeModel(t, p, isRise, J) - yy) ./ ee;
        [p, Jac, chi2] = lmFit(resid, p0, hStep(1:5), clampFn);
        if isempty(best) || chi2 < best.chi2
            best = struct('p', p, 'J', Jac, 'chi2', chi2);
        end
    end
    sgl = best;
    sgl.n = numel(t);  sgl.k = 5;

    dbl = [];
    if ~strcmp(modelChoice, 'single')
        ps = sgl.p;
        st2 = [0.05 10; -0.05 10; 0.1 50; 0.03 300];      % [w rho] starts
        for i = 1:size(st2, 1)
            p0 = [ps st2(i, :)];
            clampFn = @(p) clampP(p, win);
            resid = @(p) (edgeModel(t, p, isRise, J) - yy) ./ ee;
            [p, Jac, chi2] = lmFit(resid, p0, hStep, clampFn);
            if isempty(dbl) || chi2 < dbl.chi2
                dbl = struct('p', p, 'J', Jac, 'chi2', chi2);
            end
        end
        dbl.n = numel(t);  dbl.k = 7;
    end

    F = NaN;
    useDouble = strcmp(modelChoice, 'double');
    if ~isempty(dbl)
        F = ((sgl.chi2 - dbl.chi2) / 2) / (dbl.chi2 / max(dbl.n - dbl.k, 1));
        if strcmp(modelChoice, 'auto') && F > 20 && abs(dbl.p(6)) > 0.005
            useDouble = true;
        end
    end
    if useDouble; fit = dbl; E.model = 'double'; else; fit = sgl; E.model = 'single'; end
    E.F = F;

    p = fit.p;
    dof = max(fit.n - fit.k, 1);
    E.chi2nu = fit.chi2 / dof;
    birge = max(1, E.chi2nu);
    [C, free] = covariance(fit.J, p, win);
    C = C * birge;

    E.params = p;
    E.paramNames = {'lo', 'hi', 't0', 'tau', 'sigma', 'w', 'rho'};
    E.paramNames = E.paramNames(1:numel(p));
    E.lo = p(1);  E.hi = p(2);
    E.tau_ns = 1e3 * p(4);  E.sigma_ns = 1e3 * p(5);
    if numel(p) == 7
        E.w = p(6);  E.tauSlow_ns = 1e3 * p(4) * p(7);
    else
        E.w = 0;  E.tauSlow_ns = NaN;
    end

    % Derived times, laser (no box) and as recorded (box), with errors by propagation.
    E.laser      = derivedTimes(p, C, free, tCmd, 0, win);
    E.asRecorded = derivedTimes(p, C, free, tCmd, J, win);
    E.time = E.laser.time;   E.timeErr = E.laser.timeErr;
    E.t10 = E.laser.t(1);    E.t50 = E.laser.t(2);   E.t90 = E.laser.t(3);
    % The position of the edge also wanders from average to average (each record starts
    % at a random point of the digitizer clock), which the point-by-point fit errors do
    % not see. The spread of the per-average data crossings does: take the larger.
    E.tErr = max(E.laser.tErr, nz2(E.data.tErr));

    E.fitT = t;  E.fitY = yy;  E.fitE = ee;
    E.fitCurve = edgeModel(t, p, isRise, J);
    E.laserCurve = edgeModel(t, p, isRise, 0);
    E.resid = (yy - E.fitCurve) ./ ee;

    % Does the fit agree with the data crossings?
    E.agree = true;
    if isfinite(E.data.time)
        tol = max(3 * hypot(E.asRecorded.timeErr, nz(E.data.timeErr)), 0.1 * E.asRecorded.time);
        tol = max(tol, 1);                                         % ns
        E.agree = abs(E.asRecorded.time - E.data.time) <= tol;
    end
end

function D = dataCrossings(t, Yw, y, sw, lo, hi, isRise, tCmd, pre, post)
    % 10/50/90% completion read off the data, on the mean and per average.
    D = struct('t', nan(1, 3), 'tErr', nan(1, 3), 'time', NaN, 'timeErr', NaN);
    c = crossings(t, smoothRun(y, sw), lo, hi, isRise);
    D.t = 1e3 * (c - tCmd);
    D.time = D.t(3) - D.t(1);
    A = size(Yw, 2);
    if A >= 3
        per = nan(A, 3);
        for a = 1:A
            ya = Yw(:, a);
            if ~any(isfinite(ya)); continue; end
            % Each average against its own levels (laser power drifts between averages),
            % taken from the same stretches as the mean's.
            if isRise
                la = nmedian(ya(pre)); ha = nmedian(ya(post));
            else
                ha = nmedian(ya(pre)); la = nmedian(ya(post));
            end
            per(a, :) = crossings(t, smoothRun(ya, sw), la, ha, isRise);
        end
        ok = all(isfinite(per), 2);
        if nnz(ok) >= 3
            D.tErr = 1e3 * std(per(ok, :), 0, 1) / sqrt(nnz(ok));
            D.timeErr = 1e3 * std(per(ok, 3) - per(ok, 1)) / sqrt(nnz(ok));
        end
    end
end

function c = crossings(t, y, lo, hi, isRise)
    % Completion r rises 0 -> 1 on both edges. The 50% point is the first crossing; the
    % 10% point is the LAST time before it that r was below 10%, the 90% point the FIRST
    % time after it that r reached 90% -- robust against noise far from the edge.
    if isRise; r = (y - lo) / (hi - lo); else; r = (hi - y) / (hi - lo); end
    c = nan(1, 3);
    i50 = find(r >= 0.5, 1);
    if isempty(i50) || i50 < 2; return; end
    c(2) = interpX(t, r, i50 - 1, i50, 0.5);
    k = find(r(1:i50) < 0.1, 1, 'last');
    if ~isempty(k) && k < numel(t); c(1) = interpX(t, r, k, k + 1, 0.1); end
    k = find(r(i50:end) >= 0.9, 1);
    if ~isempty(k)
        k = k + i50 - 1;
        c(3) = interpX(t, r, k - 1, k, 0.9);
    end
end

function x = interpX(t, r, i, j, level)
    if r(j) == r(i)
        x = t(j);
    else
        x = t(i) + (level - r(i)) * (t(j) - t(i)) / (r(j) - r(i));
    end
end

% -------------------------------------------------------------------------
% The model

function y = edgeModel(t, p, isRise, J)
    G = shapeG(t, p, J);
    if isRise
        y = p(1) + (p(2) - p(1)) * G;
    else
        y = p(1) + (p(2) - p(1)) * (1 - G);
    end
end

function G = shapeG(t, p, J)
    % Normalised edge (0 -> 1), smeared by a centred box of width J (trigger alignment).
    % The box average is exact: (I(t + J/2) - I(t - J/2)) / J with I the integral of the
    % step (stepInt), so it costs two evaluations instead of a sum over many shifts.
    if J > 0
        G = (bareI(t + J/2, p) - bareI(t - J/2, p)) / J;
    else
        G = bareG(t, p);
    end
end

function G = bareG(t, p)
    G = stepResp(t, p(3), p(4), p(5));
    if numel(p) >= 7 && p(6) ~= 0
        G = (1 - p(6)) * G + p(6) * stepResp(t, p(3), p(4) * p(7), p(5));
    end
end

function I = bareI(t, p)
    I = stepInt(t, p(3), p(4), p(5));
    if numel(p) >= 7 && p(6) ~= 0
        I = (1 - p(6)) * I + p(6) * stepInt(t, p(3), p(4) * p(7), p(5));
    end
end

function s = stepResp(t, t0, tau, sg)
% H(t-t0)*(1 - exp(-(t-t0)/tau)) convolved with a Gaussian of width sg (calibChannelDelays):
%   S = Phi(x/sg) - Q,   Q = exp(-x/tau + sg^2/(2 tau^2)) * Phi(x/sg - sg/tau),   x = t - t0
    x = t - t0;
    if sg < 1e-6 * tau
        s = (x > 0) .* (1 - exp(-max(x, 0) / tau));
        return
    end
    s = 0.5 * erfc(-x / (sg*sqrt(2))) - qTerm(x, tau, sg);
end

function I = stepInt(t, t0, tau, sg)
% Integral of stepResp from -inf to t:
%   I = (x - tau) Phi(x/sg) + sg phi(x/sg) + tau Q
% (d/dx gives Phi - Q = S, since dQ/dx = -Q/tau + phi(x/sg)/sg).
    x = t - t0;
    if sg < 1e-6 * tau
        xp = max(x, 0);
        I = (x > 0) .* (xp - tau + tau * exp(-xp / tau));
        return
    end
    I = (x - tau) .* (0.5 * erfc(-x / (sg*sqrt(2)))) + sg * exp(-x.^2 / (2*sg^2)) / sqrt(2*pi) ...
        + tau * qTerm(x, tau, sg);
end

function Q = qTerm(x, tau, sg)
% exp(-x/tau + sg^2/(2 tau^2)) * Phi(x/sg - sg/tau), in the form that cannot overflow:
% directly when z = x/sg - sg/tau >= 0, otherwise through erfcx.
    z = x / sg - sg / tau;
    Q = zeros(size(x));
    pos = z >= 0;
    Q(pos)  = exp(-x(pos)/tau + sg^2/(2*tau^2)) .* 0.5 .* erfc(-z(pos)/sqrt(2));
    Q(~pos) = 0.5 * exp(-x(~pos).^2 / (2*sg^2)) .* erfcx(-z(~pos)/sqrt(2));
end

function q = offsetAt(frac, p, J)
    % Time after t0 at which G first reaches frac (p(3) is ignored and taken as 0).
    p(3) = 0;
    span = 10 * p(5) + 60 * p(4);
    if numel(p) >= 7; span = span + 60 * p(4) * p(7); end
    tg = linspace(-10 * p(5) - J - 1e-3, span + J, 4001)';
    G = shapeG(tg, p, J);
    i = find(G >= frac, 1);
    if isempty(i) || i == 1; q = NaN; return; end
    q = interpX(tg, G, i - 1, i, frac);
end

function p = clampP(p, win)
    p(3) = min(max(p(3), win(1) - 0.5), win(2));
    p(4) = min(max(p(4), 1e-4), 5);          % tau: 0.1 ns .. 5 us
    p(5) = min(max(p(5), 1e-4), 2);          % sigma
    if numel(p) >= 7
        p(6) = min(max(p(6), -0.5), 0.95);   % weight of the slow part
        p(7) = min(max(p(7), 2), 1e4);       % slow tau / fast tau
    end
end

function [C, free] = covariance(Jac, p, win)
    % Parameters sitting on a bound are held fixed; the rest from a column-scaled J'J
    % (a plain inverse can lose a direction when the scales differ by 1e6).
    lo = [-Inf -Inf win(1)-0.5 1e-4 1e-4 -0.5 2];
    hi = [ Inf  Inf win(2)     5    2    0.95 1e4];
    n = numel(p);
    atB = abs(p - lo(1:n)) <= 1e-9 * max(1, abs(lo(1:n))) | abs(p - hi(1:n)) <= 1e-9 * max(1, abs(hi(1:n)));
    free = find(~atB);
    Jf = Jac(:, free);
    sc = sqrt(sum(Jf.^2, 1));  sc(sc == 0) = 1;
    Js = Jf ./ sc;
    Cs = pinv(Js' * Js);
    C = zeros(n);
    C(free, free) = Cs ./ (sc' * sc);
end

function T = derivedTimes(p, C, free, tCmd, J, win)
    % t10/t50/t90 (ns after the command) and the 10-90 time, with errors.
    f = @(pp) completion(pp, J, win);
    v = f(p);
    T.t = 1e3 * (v - tCmd);
    T.time = 1e3 * (v(3) - v(1));
    D = zeros(4, numel(p));
    for k = free(:)'
        h = max(1e-6, 1e-4 * abs(p(k)));
        if k <= 2; h = 1e-4 * abs(p(2) - p(1)); end
        pp = p; pp(k) = pp(k) + h;  pm = p; pm(k) = pm(k) - h;
        pp = clampP(pp, win);  pm = clampP(pm, win);
        hh = pp(k) - pm(k);
        if hh == 0; continue; end
        vp = f(pp);  vm = f(pm);
        D(1:3, k) = (vp - vm)' / hh;
        D(4, k)   = ((vp(3) - vp(1)) - (vm(3) - vm(1))) / hh;
    end
    Cd = D * C * D';
    sd = 1e3 * sqrt(max(diag(Cd), 0))';
    T.tErr = sd(1:3);
    T.timeErr = sd(4);
end

function v = completion(p, J, win)
    % Absolute times (us) where G first reaches 10, 50, 90%: a 1 ns grid to bracket the
    % crossing, then a 5 ps grid inside the bracket (vectorised, no scalar loop).
    a = min(win(1), p(3) - 10*p(5) - J) - 1e-3;
    b = win(2) + 5 * max(p(4), 0) + J;
    if numel(p) >= 7; b = max(b, p(3) + 10*p(5) + 12 * p(4) * p(7)); end
    tg = (a : 1e-3 : b)';
    if numel(tg) > 2e5; tg = linspace(a, b, 2e5)'; end
    G = shapeG(tg, p, J);
    lev = [0.1 0.5 0.9];
    v = nan(1, 3);
    for i = 1:3
        k = find(G >= lev(i), 1);
        if isempty(k) || k == 1; continue; end
        tf = linspace(tg(k-1), tg(k), 201)';
        Gf = shapeG(tf, p, J);
        j = find(Gf >= lev(i), 1);
        if isempty(j) || j == 1; v(i) = tf(end); continue; end
        v(i) = interpX(tf, Gf, j - 1, j, lev(i));
    end
end

% -------------------------------------------------------------------------
% Levenberg-Marquardt (as calibChannelDelays)

function [p, J, chi2] = lmFit(resid, p, h, clampFn)
    ids = {'MATLAB:singularMatrix', 'MATLAB:nearlySingularMatrix'};
    if exist('OCTAVE_VERSION', 'builtin')
        ids = {'Octave:singular-matrix', 'Octave:nearly-singular-matrix'};
    end
    old = {};
    for k = 1:numel(ids)
        try
            old{end+1} = warning('query', ids{k}); %#ok<AGROW>
            warning('off', ids{k});
        catch
        end
    end
    cleanup = onCleanup(@() restoreWarnings(old)); %#ok<NASGU>

    p = clampFn(p);
    r = resid(p);
    chi2 = r' * r;
    lambda = 1e-3;
    small = 0;
    J = jacobianOf(resid, p, h);
    for it = 1:200
        A = J' * J;
        g = J' * r;
        dA = diag(A);
        dA(dA <= 0) = max(dA) * 1e-12 + eps;
        stepped = false;
        while lambda <= 1e12
            dp = -(A + lambda * diag(dA)) \ g;
            pn = clampFn(p + dp');
            rn = resid(pn);
            chin = rn' * rn;
            if isfinite(chin) && chin < chi2
                stepped = true;
                break
            end
            lambda = lambda * 10;
        end
        if ~stepped; break; end
        rel = (chi2 - chin) / chi2;
        p = pn;  r = rn;  chi2 = chin;
        lambda = max(lambda / 10, 1e-12);
        J = jacobianOf(resid, p, h);
        if rel < 1e-8
            small = small + 1;
            if small >= 3; break; end
        else
            small = 0;
        end
    end
end

function restoreWarnings(old)
    for k = 1:numel(old)
        warning(old{k}.state, old{k}.identifier);
    end
end

function J = jacobianOf(resid, p, h)
    r0 = resid(p);
    J = zeros(numel(r0), numel(p));
    for k = 1:numel(p)
        dp = zeros(size(p));  dp(k) = h(k);
        J(:, k) = (resid(p + dp) - resid(p - dp)) / (2*h(k));
    end
end

% =========================================================================
% REPORT AND PLOT
% =========================================================================

function report(out)
    o = out.opts;  lv = out.levels;  info = out.info;
    fprintf('\n%d averages', out.nAverages);
    if strcmp(info.source, 'hardware'); fprintf(', taken %s', info.when); end
    fprintf('.  Reference photodiode: dark %.4g V, laser on %.4g V (step %.4g V)\n', ...
            lv.dark, lv.on, lv.step);
    if out.triggerJitter > 0
        fprintf(['Trigger alignment: a %g ns box (digitizer clock not locked to the PulseStreamer) ' ...
                 'is removed in the "laser" numbers.\n'], 1e3 * out.triggerJitter);
    end

    fprintf('\n                 10-90%% time         completion after the command (ns)\n');
    fprintf('                                     10%%                50%%               90%%\n');
    printEdge('RISE', out.rise);
    printEdge('FALL', out.fall);
    fprintf('  (fall: 90%% complete = 10%% of the light left)\n');

    for E = [out.rise, out.fall]
        if strcmp(E.model, 'double')
            m = sprintf('erf x exp + slow part (w %+.3f, tau %.0f ns)', E.w, E.tauSlow_ns);
        else
            m = 'erf x exp';
        end
        fprintf('\n%s fit: %s, sigma %.1f ns, tau %.1f ns, chi2/dof %.2f (%d points)', ...
                upper(E.kind), m, E.sigma_ns, E.tau_ns, E.chi2nu, numel(E.fitT));
        if out.triggerJitter > 0
            fprintf('\n   as recorded (with the box): %.1f +/- %.1f ns', E.asRecorded.time, E.asRecorded.timeErr);
        end
        fprintf('\n   straight from the data:     %.1f +/- %.1f ns  (10%% at %.1f, 90%% at %.1f ns)', ...
                E.data.time, E.data.timeErr, E.data.t(1), E.data.t(3));
        if ~E.agree
            fprintf(['\n   NOTE the fit and the data crossings disagree: the model does not describe ' ...
                     'this edge.\n        Trust the data value, and look at the residual panel.']);
        end
        if E.chi2nu > 3
            fprintf('\n   NOTE chi2/dof %.1f: structure the model does not follow (see residuals).', E.chi2nu);
        end
    end
    fprintf('\n');

    % Checks.
    if lv.clippedFraction > 0
        fprintf(['\nWARNING %.2g%% of the samples sit at the edge of the +/-%g V input range: the ' ...
                 'reference signal is CLIPPED.\n        Clipping flattens the top of the edge and ' ...
                 'makes the times wrong. Lower the light on the\n        reference photodiode ' ...
                 '(or its gain) and rerun.\n'], 100 * lv.clippedFraction, info.maxVoltage);
    elseif isfinite(lv.rangeUsed) && lv.rangeUsed > 0.9
        fprintf('\nNOTE the signal uses %.0f%% of the input range: close to clipping.\n', 100 * lv.rangeUsed);
    elseif isfinite(lv.rangeUsed) && lv.rangeUsed < 0.02
        fprintf('\nNOTE the step is only %.1f%% of the input range: small signal, larger errors.\n', ...
                100 * abs(lv.step) / info.maxVoltage);
    end
    if abs(lv.plateauChange) > 0.01
        fprintf(['NOTE the light changes by %+.1f%% of the step during the %g us pulse (AOM heating ' ...
                 'or laser): not part of the edge.\n'], 100 * lv.plateauChange, o.onTime);
    end
    if abs(lv.darkRecovery) > 3e-3 && abs(lv.darkRecovery * lv.step) > 3 * lv.noisePerPoint
        fprintf(['NOTE %g us after the off command the reading is still %+.2g of the step away from ' ...
                 'the dark\n     before the pulse: a slow tail (photodiode or light). A longer ' ...
                 'darkAfter shows how slow.\n'], o.darkAfter, lv.darkRecovery);
    end

    st = getf(info, 'stored', struct('on', NaN, 'off', NaN));
    if isfinite(st.on)
        fprintf(['\nFor orientation: setupInfo.json greenLaser onDelay %.1f ns, offDelay %.1f ns ' ...
                 '(90%% convention, measured\non the NV fluorescence at level 1). Reference ' ...
                 'photodiode 90%% points: on %.1f ns, off %.1f ns.\nThe two need not agree: ' ...
                 'different detector, cable and the NV response. Do not edit the json from this.\n'], ...
                1e3 * st.on, 1e3 * st.off, out.rise.t90, out.fall.t90);
    end
    fprintf('\n');
end

function printEdge(name, E)
    fprintf('  %s   %7.1f +/- %-5.2g ns  %7.1f +/- %-5.2g %7.1f +/- %-5.2g %7.1f +/- %-5.2g\n', ...
            name, E.time, E.timeErr, E.t10, E.tErr(1), E.t50, E.tErr(2), E.t90, E.tErr(3));
end

function plotIt(out)
    figure('Name', 'laserRiseFall', 'Color', 'w');
    tau = out.tau;  y = nmean(out.Y, 2);
    info = out.info;

    subplot(3, 2, [1 2]);
    plot(tau, y, 'b-'); hold on;
    yl = [min(y) max(y)] + [-0.05 0.05] * (max(y) - min(y));
    plot(info.tCmdOn  * [1 1], yl, 'k:');
    plot(info.tCmdOff * [1 1], yl, 'k:');
    ylim(yl); grid on;
    xlabel('time after the digitizer trigger (us)'); ylabel('reference PD (V)');
    title(sprintf('green laser on the reference photodiode: rise %.1f ns, fall %.1f ns (10-90%%)', ...
                  out.rise.time, out.fall.time));

    edges = {out.rise, out.fall};
    for i = 1:2
        E = edges{i};
        x = 1e3 * (E.fitT - E.commanded);
        subplot(3, 2, 2 + i);
        plot(x, E.fitY, '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 4); hold on;
        plot(x, E.fitCurve, 'r-', 'LineWidth', 1.2);
        if out.triggerJitter > 0
            plot(x, E.laserCurve, 'b--', 'LineWidth', 1);
        end
        lims = [min(E.fitY) max(E.fitY)];
        for k = [1 3]
            plot(E.laser.t(k) * [1 1], lims, 'k:');
        end
        c50 = E.t50;
        xlim(c50 + [-1 1] * max(4 * max(E.time, 5), 30));
        grid on;
        xlabel(sprintf('ns after the %s command', ternary(i == 1, 'on', 'off')));
        ylabel('V');
        if out.triggerJitter > 0
            legend({'data', 'fit (as recorded)', 'laser (box removed)'}, 'Location', 'best');
        else
            legend({'data', 'fit'}, 'Location', 'best');
        end
        title(sprintf('%s %s: %.1f +/- %.1f ns (data %.1f)', upper(E.kind), ...
                      ternary(i == 1, '10-90%', '90-10%'), E.time, E.timeErr, E.data.time));

        subplot(3, 2, 4 + i);
        plot(x, E.resid, '.', 'MarkerSize', 4); hold on;
        plot([min(x) max(x)], [0 0], 'k-');
        grid on;
        xlabel(sprintf('ns after the %s command', ternary(i == 1, 'on', 'off')));
        ylabel('residual / error');
        title(sprintf('chi2/dof %.2f, model %s', E.chi2nu, E.model));
    end
end

% =========================================================================
% SMALL HELPERS
% =========================================================================

function tc = crossTime(t, y, level, after, direction)
% First time at or after 'after' where y crosses 'level' going in 'direction'.
    tc = NaN;
    idx = find(t >= after);
    if numel(idx) < 2; return; end
    v = direction * (y - level);
    for k = idx(1)+1 : numel(t)
        if v(k-1) < 0 && v(k) >= 0
            tc = interpX(t, v, k - 1, k, 0);
            return
        end
    end
end

function ys = smoothRun(y, n)
    % Running mean over n points (n = 1: unchanged). NaN-tolerant.
    if n <= 1; ys = y; return; end
    k = ones(n, 1);
    ok = isfinite(y);
    y0 = y;  y0(~ok) = 0;
    ys = conv(y0, k, 'same') ./ max(conv(double(ok), k, 'same'), 1);
    ys(~ok) = NaN;
end

function s = robustNoise(y)
    % Noise of one point, from point-to-point differences (so the step does not count).
    d = diff(y(isfinite(y)));
    s = 1.4826 * median(abs(d - median(d))) / sqrt(2);
    if ~(s > 0); s = std(d) / sqrt(2); end      % e.g. a clipped, perfectly flat plateau
    if ~(s > 0); s = eps; end
end

function m = nmean(X, dim)
    ok = isfinite(X);
    X(~ok) = 0;
    m = sum(X, dim) ./ sum(ok, dim);
end

function s = nstd(X, dim)
    m = nmean(X, dim);
    ok = isfinite(X);
    D = X - m;  D(~ok) = 0;
    s = sqrt(sum(D.^2, dim) ./ max(sum(ok, dim) - 1, 1));
end

function m = nmedian(x)
    x = x(isfinite(x));
    if isempty(x); m = NaN; else; m = median(x); end
end

function v = getf(s, f, d)
    if isstruct(s) && isfield(s, f); v = s.(f); else; v = d; end
end

function x = nz(x)
    if ~isfinite(x); x = 0; end
end

function x = nz2(x)
    x(~isfinite(x)) = 0;
end

function out = ternary(c, a, b)
    if c; out = a; else; out = b; end
end

function o = withDefaults(o, d, who)
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
