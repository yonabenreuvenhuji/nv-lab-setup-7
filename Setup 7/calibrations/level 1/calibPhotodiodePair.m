function out = calibPhotodiodePair(opts)
%CALIBPHOTODIODEPAIR Timing and noise of the two photodiodes: fluorescence (input_1) and laser monitor (input_2).
%
%   out = calibPhotodiodePair                               % both parts, defaults below
%   out = calibPhotodiodePair(struct('doNoise', false))     % edges only
%   out = calibPhotodiodePair(struct('replay', out))        % re-analyse a saved result, no hardware
%
% LEVEL 1, cell 1.6. No microwave, no field. Writes nothing: the one number it can suggest,
% spcm.channel_delay(2) in setupInfo.json, is a hand edit like every other delay in that
% file. Leaves the laser off, the photodiode gate off and fixDelays as it found them.
%
% ===================================================================
% WHAT THE QUESTION IS
% ===================================================================
%
% The digitizer samples both inputs on the same trigger, in the same segment, so there is
% no "read delay" between them today. There is a knob for one: spcm.channel_delay in
% setupInfo.json, per digitizer input, in us, in whole 10 ns samples (negative allowed;
% Digitizer.readExperimentData crops each input with its own offset). It is [0 0 0 0] and
% nothing has ever measured what it should be.
%
% The two signals do arrive at different times: the monitor sees the beam after the AOM,
% the NVs' light has the optical path, the ~12 ns NV lifetime and its own photodiode,
% amplifier and cable. A FIXED offset does not matter -- 'referenced' mode divides each
% window by the monitor's window mean normalised to its mean over the read, so a constant
% ratio divides out. What can matter is the edge in a PULSED readout window: the digitizer
% clock is not locked to the PulseStreamer, so every record starts 0-10 ns at random and
% the laser edge moves inside the window. If both inputs see the edge at the same place
% the movement cancels in S/R; if they are offset it does not. CW windows have no edge.
%
% And there is a second question the same data answers, which is probably the bigger one:
% does dividing by the monitor actually remove laser noise? 'referenced' divides RAW
% input_1 by RAW input_2, with no dark level taken off either. The residual laser-noise
% coupling is then (light fraction of input_1) - (light fraction of input_2), not zero.
% Level 0 found the fluorescence linear in power, so subtracting the darks first would make
% the cancellation exact -- if the monitor's own noise does not eat the gain.
%
% ===================================================================
% PART 1 -- THE EDGES  (fixDelays = 0, absolute times)
% ===================================================================
%
% laserRiseFall's measurement, with BOTH rows kept: every record holds input_1 and input_2
% from the same trigger, so the 0-10 ns record-start ambiguity is common to both and drops
% out of their difference. Equivalent-time sampling (nShifts laser shifts of 1 ns) gives
% 1 ns resolution from the 10 ns digitizer. Each edge of each row is fitted with the same
% model as laserRiseFall and level 1 (Gaussian-smoothed Heaviside x exponential, box for the
% trigger alignment); the fluorescence's NV lifetime simply comes out in its tau.
%
%   offset = t50(monitor) - t50(fluorescence), on the rise and on the fall
%
% from the fits and, model-free, from the per-average 50% crossings PAIRED average by
% average (the error comes from that paired spread). Positive = the monitor arrives later.
% Suggested channel_delay(2) = current + the rise offset rounded to 10 ns. If the rise and
% fall offsets differ by more than a sample, the two detectors have different bandwidths
% and no single shift lines both edges up; part 2 then decides.
%
% ===================================================================
% PART 2 -- THE PAYOFF  (fixDelays = 1, the readout window experiments use)
% ===================================================================
%
% Each repeat is two identical halves, as a balanced experiment's signal and reference:
%   [dark gap] [laser on laserInitializationDuration, window = detectionDuration at its start]
% The digitizer segment starts `margin` before the window and ends `margin` after it, so the
% monitor's window can be slid against the fluorescence's in software, +/- maxShift samples,
% on the same data. A dark run (same sequence, laser never opened) gives both dark levels.
%
% For every window k: S_k (input_1 mean over the window), R_k(s) (input_2, window shifted s
% samples). Three ways to normalise, each at every shift s:
%   none      S_k                                       ('nomral' input configuration)
%   raw       S_k / (R_k / <R>)                         what 'referenced' does today
%   dark      (S_k - dS) / ((R_k - dR) / <R - dR>)      darks off both first
% <..> over the read (one average), as readFromExperiment does. Then the ratio of the two
% halves of each repeat -- the quantity a balanced experiment's contrast is built from --
% and its spread over all repeats. Reported relative to `none`, with bootstrap errors
% (paired: same repeats, so the comparison is much tighter than either number alone).
% The `dark` spread is multiplied by the fluorescence's light fraction fS: taking dS off
% the numerator makes the spin contrast 1/fS larger and the noise with it, so only the
% scaled number is noise per unit of real contrast, comparable with the other two.
%
% It also prints kappa = cov(ln S, ln R)/var(ln R) per window: 1 would mean the monitor
% tracks the fluorescence one for one; dividing by R^kappa is the best a pure power law can do.
%
% ===================================================================
% OPTIONS (struct fields, all optional)
% ===================================================================
%   doEdges, doNoise   run each part (default true, true)
%   PART 1   onTime 10, darkBefore 3, darkAfter 5, deadTime 5 (us), edgeRepeats 200,
%            edgeAverages 20, shiftStep 0.001 us, nShifts [] (fill one sample = 10),
%            edgeWindow [0.5 1] us before/after the 50% point, triggerJitter 'auto' or us,
%            model 'auto' | 'single' | 'double'
%   PART 2   window [] (json detectionDuration), laserOn [] (json laserInitializationDuration),
%            gap 1 us of dark before each half, margin 0.1 us, maxShift 5 samples,
%            noiseRepeats 1000, noiseAverages 10, darkAverages 2, bootstrap 300
%   jsonLocation  expProperties.json (default setup7Paths().json, only for the two timings)
%   plot          default true
%   replay        a previous `out`: re-analyse, no hardware (its timing is kept)
%
% OUTPUT
%   out.edges   .fl / .ref (rise and fall fits, levels), .offset (ns, fit and paired data),
%               .suggestedChannelDelay (us, absolute), .currentChannelDelay
%   out.noise   .table (config x shift: sigma of the half ratio, relative to `none`, errors),
%               .bestShift, .kappa, .darks, .levels
%   out.verdict text lines, as printed
%   out.raw     what replay needs

    if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
    opts = withDefaults(opts, struct( ...
        'doEdges', true, 'doNoise', true, ...
        'onTime', 10, 'darkBefore', 3, 'darkAfter', 5, 'deadTime', 5, ...
        'edgeRepeats', 200, 'edgeAverages', 20, 'shiftStep', 0.001, 'nShifts', [], ...
        'edgeWindow', [0.5 1], 'triggerJitter', 'auto', 'model', 'auto', ...
        'window', [], 'laserOn', [], 'gap', 1, 'margin', 0.1, 'maxShift', 5, ...
        'noiseRepeats', 1000, 'noiseAverages', 10, 'darkAverages', 2, 'bootstrap', 300, ...
        'jsonLocation', '', 'plot', true, 'replay', []), 'calibPhotodiodePair');

    replay = ~isempty(opts.replay);
    out = struct('opts', opts, 'edges', [], 'noise', [], 'verdict', {{}}, 'raw', struct());
    if replay
        R = opts.replay;
        if ~isstruct(R) || ~isfield(R, 'raw')
            error('calibPhotodiodePair: replay must be an earlier output (with a raw field).');
        end
        keep = {'onTime', 'darkBefore', 'darkAfter', 'deadTime', 'edgeRepeats', 'edgeAverages', ...
                'shiftStep', 'nShifts', 'window', 'laserOn', 'gap', 'margin', 'maxShift', ...
                'noiseRepeats', 'noiseAverages', 'darkAverages'};
        for i = 1:numel(keep)
            if isfield(R.opts, keep{i}); opts.(keep{i}) = R.opts.(keep{i}); end
        end
        opts.doEdges = isfield(R.raw, 'edges') && ~isempty(R.raw.edges);
        opts.doNoise = isfield(R.raw, 'noise') && ~isempty(R.raw.noise);
        fprintf('\n=== calibPhotodiodePair: REPLAY, no hardware touched ===\n');
    else
        fprintf('\n=== calibPhotodiodePair: the two photodiodes, laser only ===\n');
    end
    opts = checkTiming(opts);
    if ~any(strcmpi(opts.model, {'auto', 'single', 'double'}))
        error('calibPhotodiodePair: model must be ''auto'', ''single'' or ''double''.');
    end
    opts.model = lower(opts.model);
    out.opts = opts;

    % ---- hardware session (one cleanup for both parts) ------------------------
    dev = [];
    if ~replay
        dev = getDevices();
        oldFix = dev.pg.fixDelays;
        cleanup = onCleanup(@() restoreState(dev, oldFix)); %#ok<NASGU>
        if opts.doNoise
            opts = timingsFromJson(opts);
            out.opts = opts;
        end
    end

    % ---- part 1 -----------------------------------------------------------------
    if opts.doEdges
        if replay
            E = R.raw.edges;
        else
            E = acquireEdges(dev, opts);
        end
        out.raw.edges = E;
        if ischar(opts.triggerJitter)
            if ~strcmpi(opts.triggerJitter, 'auto')
                error('calibPhotodiodePair: triggerJitter must be ''auto'' or a number of us.');
            end
            J = E.info.triggerJitterAuto;
        else
            J = opts.triggerJitter;
        end
        out.edges = analyseEdges(E, opts, J);
    end

    % ---- part 2 -----------------------------------------------------------------
    if opts.doNoise
        if replay
            N = R.raw.noise;
        else
            N = acquireNoise(dev, opts);
        end
        out.raw.noise = N;
        out.noise = analyseNoise(N, opts);
    end

    out.verdict = verdict(out);
    report(out);
    if opts.plot
        try
            plotIt(out);
        catch err
            fprintf('(plot failed: %s)\n', err.message);
        end
    end
end

% =========================================================================
% TIMING AND SET-UP
% =========================================================================

function opts = checkTiming(opts)
    r10 = @(x) round(x * 100) / 100;            % whole 10 ns digitizer samples
    opts.onTime     = r10(opts.onTime);
    opts.darkBefore = r10(opts.darkBefore);
    opts.darkAfter  = r10(opts.darkAfter);
    opts.deadTime   = r10(opts.deadTime);
    opts.shiftStep  = max(1, round(opts.shiftStep * 1e3)) / 1e3;
    if isempty(opts.nShifts)
        opts.nShifts = max(1, round(0.01 / opts.shiftStep));
    end
    opts.nShifts = round(opts.nShifts);
    w = opts.edgeWindow;
    if ~(isnumeric(w) && numel(w) == 2 && all(w > 0))
        error('calibPhotodiodePair: edgeWindow must be [before after] in us, both > 0.');
    end
    maxShift = (opts.nShifts - 1) * opts.shiftStep;
    if opts.darkBefore < max(1, w(1))
        error('calibPhotodiodePair: darkBefore must be at least max(1, edgeWindow(1)) us.');
    end
    if opts.darkAfter < w(2) + 1 + maxShift
        error('calibPhotodiodePair: darkAfter must be at least edgeWindow(2) + 1 us.');
    end
    if opts.onTime < sum(w) + 1
        error('calibPhotodiodePair: onTime must be at least sum(edgeWindow) + 1 us.');
    end
    if opts.deadTime < 2
        error('calibPhotodiodePair: deadTime must be at least 2 us (digitizer re-arm).');
    end
    opts.margin   = r10(opts.margin);
    opts.maxShift = round(opts.maxShift);
    if opts.maxShift < 0
        error('calibPhotodiodePair: maxShift must be >= 0 samples.');
    end
    if opts.margin * 100 < opts.maxShift
        error('calibPhotodiodePair: margin (%g us) must cover maxShift (%d samples of 10 ns).', ...
              opts.margin, opts.maxShift);
    end
    if opts.gap < opts.margin + 0.05
        error('calibPhotodiodePair: gap must be at least margin + 0.05 us.');
    end
end

function opts = timingsFromJson(opts)
    % The readout window and the laser time of the real experiments, unless given.
    if ~isempty(opts.window) && ~isempty(opts.laserOn); return; end
    loc = opts.jsonLocation;
    if isempty(loc)
        loc = setup7Paths().json;
    end
    j = jsondecode(fileread(loc));
    if isempty(opts.window)
        opts.window = j.detectionDuration;
        fprintf('readout window %g us (detectionDuration, %s)\n', opts.window, loc);
    end
    if isempty(opts.laserOn)
        opts.laserOn = j.laserInitializationDuration;
        fprintf('laser on %g us per half (laserInitializationDuration)\n', opts.laserOn);
    end
    opts.window  = round(opts.window * 100) / 100;
    opts.laserOn = round(opts.laserOn * 100) / 100;
    if opts.laserOn < opts.window + opts.margin
        error('calibPhotodiodePair: laserOn must be at least window + margin.');
    end
end

function dev = getDevices()
    dev.pg   = getObjByName(PulseGenerator.NAME);
    dev.spcm = getObjByName(Spcm.NAME);
    dev.digi = getObjByName(Digitizer.NAME);
    if isempty(dev.pg) || isempty(dev.spcm) || isempty(dev.digi)
        error(['calibPhotodiodePair: pulse generator, photodiode or digitizer object not found. ' ...
               'Start the setup (the usual init) first.']);
    end
    try
        names = dev.spcm.digitizerChannelNames;
        vals  = dev.spcm.digitizerChannelValues;
    catch
        error('calibPhotodiodePair: the detector (%s) is not the digitizer photodiode.', class(dev.spcm));
    end
    i1 = find(strcmp(names, 'input_1'), 1);
    i2 = find(strcmp(names, 'input_2'), 1);
    if isempty(i1) || isempty(i2)
        error('calibPhotodiodePair: setupInfo.json must name both input_1 and input_2 (has: %s).', ...
              strjoin(names, ', '));
    end
    % BOTH inputs, always: on Setup 7 acquiring input_2 on its own returns zeros (the same
    % fault as the panel's 'reference' option), while [1 2] -- what every experiment
    % acquires -- works.
    dev.chans = [vals(i1), vals(i2)];
    dev.fs = dev.digi.MAX_SAMPLE_RATE;          % samples per us
    dev.currentDelay = [0 0];
    try
        cd = dev.digi.channelDelay;
        dev.currentDelay = [cd(dev.chans(1)), cd(dev.chans(2))];
    catch
    end
    dev.maxV = [NaN NaN];  dev.offV = [0 0];
    try; dev.maxV = dev.digi.maxVoltage(dev.chans); catch; end
    try; dev.offV = dev.digi.voltageOffset(dev.chans); catch; end
end

function restoreState(dev, oldFix)
    try; dev.digi.stopAcquisition; catch; end
    try; dev.pg.fixDelays = oldFix; catch; end
    try; dev.spcm.setSPCMEnable(false); catch; end
end

function V = readTwo(dev, nSeg, nSamp)
    % One acquisition, both rows, no zero suppression; records that never triggered (all
    % zeros in both rows) come back as NaN.
    V = [];
    for trial = 1:3
        try
            dev.digi.startExperiment;
            dev.pg.run;
            V = dev.digi.readExperimentData(dev.chans, true, false);
            break
        catch err
            warning('calibPhotodiodePair:trial', 'trial %d failed: %s', trial, err.message);
            try; dev.digi.stopAcquisition; catch; end
            V = [];
        end
    end
    if isempty(V); return; end
    if size(V, 1) ~= 2 || size(V, 2) ~= nSeg || size(V, 3) < nSamp
        error('calibPhotodiodePair: digitizer returned %s, expected [2 x %d x >=%d].', ...
              mat2str(size(V)), nSeg, nSamp);
    end
    V = V(:, :, 1:nSamp);
    z1 = reshape(all(V(1, :, :) == 0, 3), 1, nSeg);
    z2 = reshape(all(V(2, :, :) == 0, 3), 1, nSeg);
    hasNaN = reshape(any(any(isnan(V), 3), 1), 1, nSeg);
    missed = (z1 & z2) | hasNaN;
    V(:, missed, :) = NaN;
end

% =========================================================================
% PART 1: ACQUISITION
% =========================================================================

function E = acquireEdges(dev, opts)
    fs = dev.fs;
    K = opts.nShifts;
    s = (0:K-1) * opts.shiftStep;
    tPre = opts.darkBefore;  tOn = opts.onTime;  tPost = opts.darkAfter;  tDead = opts.deadTime;
    segDur = tPre + tOn + tPost;
    nSamp = round(segDur * fs);
    nSeg  = opts.edgeRepeats * K;
    ns = @(x) round(x * 1e3) / 1e3;

    S = Sequence;
    for k = 1:K
        S.addEvent(ns(tPre + s(k)),          'detector');
        S.addEvent(ns(tOn),                  'greenLaser');
        S.addEvent(ns(tPost - s(k) + tDead), '');
    end

    fprintf('\n--- part 1: the edges, both photodiodes in the same records ---\n');
    fprintf(['laser %g us on, %g dark before, %g after; %d shifts of %g ns (%g ns effective ' ...
             'sampling); %d repeats x %d averages; fixDelays = 0\n'], tOn, tPre, tPost, K, ...
            1e3 * opts.shiftStep, 1e3 / fs / K, opts.edgeRepeats, opts.edgeAverages);

    dev.pg.setSequence(S);
    dev.pg.repeats = opts.edgeRepeats;
    dev.pg.fixDelays = 0;
    dev.spcm.setSPCMEnable(true);
    dev.digi.prepareAcquisition(dev.chans, nSeg, segDur, 1, fs);

    t = (0:nSamp-1)' / fs;
    tauGrid = t - s;
    [tau, order] = sort(tauGrid(:));
    Yfl  = nan(numel(tau), opts.edgeAverages);
    Yref = nan(numel(tau), opts.edgeAverages);
    nMissed = 0;  nAvg = 0;  vMax = -Inf(1, 2);  vMin = Inf(1, 2);
    for a = 1:opts.edgeAverages
        V = readTwo(dev, nSeg, nSamp);
        if isempty(V)
            if nAvg == 0; error('calibPhotodiodePair: the first average failed three times.'); end
            fprintf('Stopping after %d averages.\n', nAvg);
            break
        end
        nMissed = nMissed + nnz(all(isnan(squeeze(V(1, :, :))), 2));
        for r = 1:2
            Vr = reshape(V(r, :, :), nSeg, nSamp);
            vMax(r) = max(vMax(r), max(Vr(:)));  vMin(r) = min(vMin(r), min(Vr(:)));
            V3 = reshape(Vr', nSamp, K, opts.edgeRepeats);
            V3 = V3(:, :, 2:end);           % first repeat follows a long laser-off
            Yk = nmean(V3, 3);
            if r == 1; Yfl(:, a) = Yk(order); else; Yref(:, a) = Yk(order); end
        end
        nAvg = a;
        fprintf('.');
        if mod(a, 50) == 0; fprintf('\n'); end
    end
    fprintf('\n');
    jitterAuto = 0;
    try
        if ~dev.digi.externalRefClock; jitterAuto = 1 / fs; end
    catch
        jitterAuto = 1 / fs;
    end
    E = struct('tau', tau, 'Yfl', Yfl(:, 1:nAvg), 'Yref', Yref(:, 1:nAvg));
    E.info = struct('tCmdOn', tPre, 'tCmdOff', tPre + tOn, 'sampleRate', fs, ...
                    'triggerJitterAuto', jitterAuto, 'missedRecords', nMissed, ...
                    'currentDelay', dev.currentDelay, 'maxVoltage', dev.maxV, ...
                    'voltageOffset', dev.offV, 'vMax', vMax, 'vMin', vMin, ...
                    'when', datestr(now, 31)); %#ok<TNOW1,DATST>
    if nMissed > 0
        fprintf('NOTE %d records came back empty (missed trigger) and were dropped.\n', nMissed);
    end
end

% =========================================================================
% PART 1: ANALYSIS
% =========================================================================

function out = analyseEdges(E, opts, J)
    info = E.info;
    fl  = analyseRow(E.tau, E.Yfl,  info, opts, J, 'fluorescence (input_1)', 1);
    ref = analyseRow(E.tau, E.Yref, info, opts, J, 'laser monitor (input_2)', 2);
    out = struct('fl', fl, 'ref', ref, 'triggerJitter', J, 'info', info);

    % Offsets, monitor minus fluorescence, at 10/50/90% of each edge.
    off.riseFit  = ref.rise.laser.t - fl.rise.laser.t;          % ns, [10 50 90]
    off.fallFit  = ref.fall.laser.t - fl.fall.laser.t;
    off.riseFitErr = hypot(ref.rise.laser.tErr, fl.rise.laser.tErr);
    off.fallFitErr = hypot(ref.fall.laser.tErr, fl.fall.laser.tErr);
    % Model-free, paired average by average: the record-start ambiguity is common to both
    % rows of an average, so it cancels in the difference rather than adding to the error.
    [off.rise50Data, off.rise50DataErr] = pairedOffset(ref.rise.data.per, fl.rise.data.per, 2);
    [off.fall50Data, off.fall50DataErr] = pairedOffset(ref.fall.data.per, fl.fall.data.per, 2);
    out.offset = off;

    % Suggested channel_delay(2): the rise 50% offset, in whole samples, on top of what is
    % already set (the measurement sees the inputs AFTER the current delays are applied).
    d50 = off.rise50Data;
    if ~isfinite(d50); d50 = off.riseFit(2); end
    step = 1e3 / info.sampleRate;                               % ns per sample
    nSampShift = round(d50 / step);
    cur = getf(info, 'currentDelay', [0 0]);
    out.residualSamples = nSampShift;
    out.currentChannelDelay = cur(2);
    out.suggestedChannelDelay = round((cur(2) + nSampShift * step / 1e3) * 100) / 100;
    out.riseFallMismatch = off.riseFit(2) - off.fallFit(2);    % ns
end

function [m, e] = pairedOffset(perRef, perFl, col)
    m = NaN;  e = NaN;
    if isempty(perRef) || isempty(perFl); return; end
    d = 1e3 * (perRef(:, col) - perFl(:, col));
    d = d(isfinite(d));
    if numel(d) >= 3
        m = mean(d);
        e = std(d) / sqrt(numel(d));
    elseif ~isempty(d)
        m = mean(d);
    end
end

function R = analyseRow(tau, Y, info, opts, J, label, rowIdx)
    good = any(isfinite(Y), 1);
    Y = Y(:, good);
    A = size(Y, 2);
    if A == 0
        error('calibPhotodiodePair: no data for the %s.', label);
    end
    y = nmean(Y, 2);
    noise = robustNoise(y);
    if A >= 3
        v = nstd(Y, 2).^2;
        v = smoothRun(v, max(1, round(0.021 / median(diff(tau)))));
        e = sqrt(v / A);
        e = max(e, 0.3 * noise);
        e(~isfinite(e)) = noise;
    else
        e = noise * ones(size(y));
    end
    tOn = info.tCmdOn;  tOff = info.tCmdOff;
    dark = nmedian(y(tau > 0.1 & tau < tOn - 0.05));
    on   = nmedian(y(tau > tOn + 0.5*opts.onTime & tau < tOff));
    step = on - dark;
    if ~(abs(step) > 5 * noise)
        error(['calibPhotodiodePair: no laser step on the %s (dark %.4g V, on %.4g V, noise ' ...
               '%.2g V per point).\nCheck the photodiode, its gate, and that the laser is on.'], ...
              label, dark, on, noise);
    end
    sgn = sign(step);
    ys = smoothRun(y, max(1, round(0.005 / median(diff(tau)))));
    half = (dark + on) / 2;
    t50on  = crossTime(tau, sgn * ys, sgn * half, tOn - 0.05, +1);
    t50off = crossTime(tau, sgn * ys, sgn * half, tOff - 0.05, -1);
    if ~isfinite(t50on) || ~isfinite(t50off)
        error('calibPhotodiodePair: could not find both laser edges on the %s.', label);
    end
    w = opts.edgeWindow;
    rise = fitEdge(tau, y, e, Y, 'rise', tOn,  t50on,  [max(tau(1), t50on - w(1)), min(t50on + w(2), tOff - 0.05)], J, opts.model);
    fall = fitEdge(tau, y, e, Y, 'fall', tOff, t50off, [max(t50off - w(1), t50on + 0.2), min(t50off + w(2), tau(end))], J, opts.model);

    lv.dark = dark;  lv.on = on;  lv.step = step;  lv.noisePerPoint = noise;
    lv.rangeUsed = NaN;
    mv = getf(info, 'maxVoltage', [NaN NaN]);
    if numel(mv) >= rowIdx && isfinite(mv(rowIdx))
        vmx = getf(info, 'vMax', [on on]);  vmn = getf(info, 'vMin', [dark dark]);
        of = getf(info, 'voltageOffset', [0 0]);
        lv.rangeUsed = max(abs([vmx(rowIdx), vmn(rowIdx)] - of(rowIdx))) / mv(rowIdx);
    end
    R = struct('label', label, 'rise', rise, 'fall', fall, 'levels', lv, 'nAverages', A, ...
               'tau', tau, 'y', y);
end

% =========================================================================
% PART 2: ACQUISITION
% =========================================================================

function N = acquireNoise(dev, opts)
    fs = dev.fs;
    W = opts.window;  L = opts.laserOn;  m = opts.margin;  g = opts.gap;
    segDur = W + 2 * m;
    nSamp = round(segDur * fs);
    nSeg = 2 * opts.noiseRepeats;
    ns = @(x) round(x * 1e3) / 1e3;

    fprintf('\n--- part 2: the payoff, at the experiments'' readout window ---\n');
    fprintf(['two halves per repeat: %g us dark, laser %g us with a %g us window at its start ' ...
             '(segment %g us);\n%d repeats x %d averages, +%d dark averages; fixDelays = 1\n'], ...
            g, L, W, segDur, opts.noiseRepeats, opts.noiseAverages, opts.darkAverages);

    dev.pg.fixDelays = 1;
    dev.spcm.setSPCMEnable(true);
    % The detector line stays high for margin + window and the segment is window + 2 margins
    % long, so its falling edge lands INSIDE the segment and cannot start a second one.
    [Sd, Rd] = runHalves(dev, opts, buildHalves(g, m, W, L, false, ns), nSeg, nSamp, ...
                         opts.darkAverages, 'dark');
    [Sl, Rl] = runHalves(dev, opts, buildHalves(g, m, W, L, true, ns), nSeg, nSamp, ...
                         opts.noiseAverages, 'laser');
    N = struct('S', {Sl}, 'R', {Rl}, 'Sdark', {Sd}, 'Rdark', {Rd}, ...
               'shifts', -opts.maxShift:opts.maxShift, 'sampleRate', fs, ...
               'currentDelay', dev.currentDelay, 'when', datestr(now, 31)); %#ok<TNOW1,DATST>
end

function S = buildHalves(g, m, W, L, laserOn, ns)
    las = 'greenLaser';
    if ~laserOn; las = ''; end
    S = Sequence;
    for h = 1:2
        S.addEvent(ns(g - m), '');
        S.addEvent(ns(m), 'detector');
        if laserOn
            S.addEvent(ns(W), {las, 'detector'});
        else
            S.addEvent(ns(W), 'detector');
        end
        S.addEvent(ns(L - W), las);
    end
end

function [Sc, Rc] = runHalves(dev, opts, Sq, nSeg, nSamp, nAv, what)
    % Per average: S (windows x 1) and R (windows x shifts), first repeat dropped.
    dev.pg.setSequence(Sq);
    dev.pg.repeats = opts.noiseRepeats;
    dev.digi.prepareAcquisition(dev.chans, nSeg, nSamp / dev.fs, 1, dev.fs);
    mS = round(opts.margin * dev.fs);
    wS = round(opts.window * dev.fs);
    shifts = -opts.maxShift:opts.maxShift;
    Sc = {};  Rc = {};
    fprintf('%s: ', what);
    for a = 1:nAv
        V = readTwo(dev, nSeg, nSamp);
        if isempty(V)
            if isempty(Sc); error('calibPhotodiodePair: the first %s average failed three times.', what); end
            break
        end
        V1 = reshape(V(1, :, :), nSeg, nSamp);
        V2 = reshape(V(2, :, :), nSeg, nSamp);
        Sa = mean(V1(:, mS + (1:wS)), 2);
        Ra = nan(nSeg, numel(shifts));
        for k = 1:numel(shifts)
            Ra(:, k) = mean(V2(:, mS + shifts(k) + (1:wS)), 2);
        end
        Sc{end+1} = Sa(3:end);                  %#ok<AGROW>   drop the first repeat
        Rc{end+1} = Ra(3:end, :);               %#ok<AGROW>
        fprintf('.');
    end
    fprintf('\n');
end

% =========================================================================
% PART 2: ANALYSIS
% =========================================================================

function out = analyseNoise(N, opts)
    shifts = N.shifts(:)';
    k0 = find(shifts == 0, 1);
    dS = nmean(vertcat(N.Sdark{:}), 1);
    Rd = vertcat(N.Rdark{:});
    dR = nmean(Rd(:, k0), 1);
    nA = numel(N.S);
    if nA == 0
        error('calibPhotodiodePair: no part-2 data.');
    end

    % Per read normalisation, as readFromExperiment does, then the half ratios.
    qNone = [];  qRaw = [];  qDark = [];
    lS = [];  lR = [];  lRd = [];
    for a = 1:nA
        S = N.S{a};  R = N.R{a};
        n = 2 * floor(numel(S) / 2);
        S = S(1:n);  R = R(1:n, :);
        ok = isfinite(S) & all(isfinite(R), 2);
        ok = ok(1:2:end) & ok(2:2:end);              % whole repeats only
        ok = reshape([ok'; ok'], [], 1);
        S = S(ok);  R = R(ok, :);
        if numel(S) < 4; continue; end
        Rn  = R ./ mean(R, 1);
        Rdn = (R - dR) ./ mean(R - dR, 1);
        xNone = S;
        xRaw  = S ./ Rn;
        xDark = (S - dS) ./ Rdn;
        qNone = [qNone; xNone(1:2:end) ./ xNone(2:2:end)];                %#ok<AGROW>
        qRaw  = [qRaw;  xRaw(1:2:end, :) ./ xRaw(2:2:end, :)];            %#ok<AGROW>
        qDark = [qDark; xDark(1:2:end, :) ./ xDark(2:2:end, :)];          %#ok<AGROW>
        lS  = [lS;  log(S / mean(S))];                                     %#ok<AGROW>
        lR  = [lR;  log(Rn)];                                              %#ok<AGROW>
        lRd = [lRd; log(max(Rdn, eps))];                                   %#ok<AGROW>
    end
    if numel(qNone) < 20
        error('calibPhotodiodePair: only %d usable repeats in part 2.', numel(qNone));
    end

    % One outlier mask for every configuration (a glitch is a glitch whichever way the
    % window is normalised), from the un-normalised ratio: > 8 robust sigma.
    med = median(qNone);  mad = 1.4826 * median(abs(qNone - med));
    keep = abs(qNone - med) <= 8 * max(mad, eps);
    keep = keep & all(isfinite(qRaw), 2) & all(isfinite(qDark), 2);
    qNone = qNone(keep);  qRaw = qRaw(keep, :);  qDark = qDark(keep, :);
    nPair = numel(qNone);

    % SAME SIGNAL, NOT SAME NUMBERS. Taking the dark off the fluorescence makes every
    % relative change in the NV light -- the spin contrast included -- 1/fS times larger in
    % the ratio, and its noise grows by the same factor, so the raw spread of a dark-
    % subtracted ratio is not comparable with the others. Scaled by fS it is: all three are
    % then noise per unit of true contrast.
    Sall0 = vertcat(N.S{:});
    fS = (nmean(Sall0, 1) - dS) / nmean(Sall0, 1);
    sNone = std(qNone);
    sRaw  = std(qRaw, 0, 1);
    sDark = std(qDark, 0, 1) * fS;

    % Bootstrap over repeats, paired: every configuration resampled with the same indices.
    B = max(0, round(opts.bootstrap));
    rRaw = nan(B, numel(shifts));  rDark = nan(B, numel(shifts));
    if B > 0
        % A fixed seed, so a replay gives the same error bars, and the caller's generator
        % state put back afterwards (experiments draw their randperm orders from it).
        % 'twister' named, because a bare rng(seed) errors if the session is in legacy mode.
        rngOld = rng;
        rng(12345, 'twister');
        for b = 1:B
            ix = ceil(rand(nPair, 1) * nPair);
            s0 = std(qNone(ix));
            rRaw(b, :)  = std(qRaw(ix, :), 0, 1) / s0;
            rDark(b, :) = std(qDark(ix, :), 0, 1) * fS / s0;
        end
        rng(rngOld);
    end
    T = struct();
    T.shifts = shifts;
    T.sigmaNone = sNone;
    T.sigmaRaw  = sRaw;
    T.sigmaDark = sDark;
    T.relRaw  = sRaw / sNone;
    T.relDark = sDark / sNone;
    T.relRawErr  = nstd(rRaw, 1);
    T.relDarkErr = nstd(rDark, 1);
    % Paired differences against the current configuration (raw, shift 0).
    T.darkVsRawErr = nstd(rDark(:, k0) - rRaw(:, k0), 1);

    % Best shift for each normalisation, preferring 0 when within one error of the minimum.
    [~, iR] = min(T.relRaw);  [~, iD] = min(T.relDark);
    if T.relRaw(k0) - T.relRaw(iR) <= max(T.relRawErr(k0), eps) * 1; iR = k0; end
    if T.relDark(k0) - T.relDark(iD) <= max(T.relDarkErr(k0), eps) * 1; iD = k0; end
    out = struct('table', T, 'k0', k0, 'bestShiftRaw', shifts(iR), 'bestShiftDark', shifts(iD), ...
                 'iBestRaw', iR, 'iBestDark', iD, 'nPairs', nPair, 'nAverages', nA, ...
                 'darks', [dS dR], 'outliers', nnz(~keep));
    % Shift-gain error from the bootstrap of the paired difference.
    if B > 0
        out.shiftGainDarkErr = nstd(rDark(:, iD) - rDark(:, k0), 1);
        out.shiftGainRawErr  = nstd(rRaw(:, iR) - rRaw(:, k0), 1);
    else
        out.shiftGainDarkErr = NaN;  out.shiftGainRawErr = NaN;
    end

    % Levels and kappa at shift 0.
    Sall = vertcat(N.S{:});  Rall = vertcat(N.R{:});
    out.levels = struct('S', nmean(Sall, 1), 'R', nmean(Rall(:, k0), 1), 'darkS', dS, 'darkR', dR);
    out.levels.lightFractionS = (out.levels.S - dS) / out.levels.S;
    out.levels.lightFractionR = (out.levels.R - dR) / out.levels.R;
    out.kappaRaw  = regressSlope(lS, lR(:, k0));
    out.kappaDark = regressSlope(lS, lRd(:, k0));
    out.relNoiseS = relNoiseOf(Sall);
    out.relNoiseR = relNoiseOf(Rall(:, k0));
end

function k = regressSlope(y, x)
    ok = isfinite(x) & isfinite(y);
    x = x(ok) - mean(x(ok));  y = y(ok) - mean(y(ok));
    k = sum(x .* y) / max(sum(x .^ 2), eps);
end

function r = relNoiseOf(x)
    % Window-to-window relative noise from consecutive differences (slow drift excluded).
    x = x(isfinite(x));
    r = std(diff(x)) / sqrt(2) / abs(mean(x));
end

% =========================================================================
% VERDICT AND REPORT
% =========================================================================

function V = verdict(out)
    V = {};
    if ~isempty(out.edges)
        E = out.edges;
        o = E.offset;
        d = o.rise50Data;  de = o.rise50DataErr;
        if ~isfinite(d); d = o.riseFit(2); de = o.riseFitErr(2); end
        if E.residualSamples == 0
            V{end+1} = sprintf(['EDGES: the monitor and the fluorescence rise within half a sample ' ...
                '(%+.1f +/- %.1f ns). channel_delay(2) = %g us stays.'], d, de, E.currentChannelDelay);
        else
            V{end+1} = sprintf(['EDGES: the monitor rises %+.1f +/- %.1f ns after the fluorescence ' ...
                '(%+d samples). To line them up: spcm.channel_delay(2) = %g us in setupInfo.json ' ...
                '(now %g) -- only if part 2 shows it pays.'], d, de, E.residualSamples, ...
                E.suggestedChannelDelay, E.currentChannelDelay);
        end
        if abs(E.riseFallMismatch) > 10
            V{end+1} = sprintf(['EDGES: the fall offset differs from the rise offset by %.1f ns: the ' ...
                'two detectors have different bandwidths (or the NV lifetime shows), so no single ' ...
                'shift lines up both edges.'], E.riseFallMismatch);
        end
    end
    if ~isempty(out.noise)
        N = out.noise;  T = N.table;  k0 = N.k0;
        raw = T.relRaw(k0);  rawE = T.relRawErr(k0);
        if raw < 1 - 2 * rawE
            V{end+1} = sprintf(['NOISE: referencing as done today (raw) cuts the half-ratio noise to ' ...
                '%.3f +/- %.3f of no referencing. It helps.'], raw, rawE);
        elseif raw > 1 + 2 * rawE
            V{end+1} = sprintf(['NOISE: referencing as done today (raw) RAISES the half-ratio noise to ' ...
                '%.3f +/- %.3f of no referencing. The monitor adds more noise than it removes; ' ...
                '''nomral'' would be quieter.'], raw, rawE);
        else
            V{end+1} = sprintf(['NOISE: referencing as done today (raw) makes no measurable difference ' ...
                '(%.3f +/- %.3f of no referencing).'], raw, rawE);
        end
        dk = T.relDark(k0);  gain = raw - dk;  gE = T.darkVsRawErr;
        if gain > 2 * gE && gain > 0.01
            V{end+1} = sprintf(['NOISE: subtracting the darks first would lower it further, to %.3f ' ...
                '(%.1f%% better than today, %.1f sigma). That needs a small change in the ' ...
                'photodiode class''s ''referenced'' path.'], dk, 100 * gain / raw, gain / max(gE, eps));
        elseif gain < -2 * gE && gain < -0.01
            V{end+1} = sprintf(['NOISE: subtracting the darks first would make it WORSE (%.3f against ' ...
                '%.3f today, %.1f sigma): the monitor''s own noise, no longer diluted by its dark ' ...
                'offset, outweighs the better cancellation. Keep the code as it is.'], dk, raw, ...
                -gain / max(gE, eps));
        else
            V{end+1} = sprintf(['NOISE: subtracting the darks first makes no measurable difference ' ...
                '(%.3f vs %.3f, %.1f sigma). Keep the code as it is.'], dk, raw, gain / max(gE, eps));
        end
        best = N.bestShiftDark;  sg = T.relDark(k0) - T.relDark(N.iBestDark);
        if best ~= 0 && sg > 2 * N.shiftGainDarkErr
            V{end+1} = sprintf(['NOISE: sliding the monitor window by %+d sample(s) lowers the noise by ' ...
                '%.3f (%.1f sigma), dark-subtracted. That is the delay to set.'], best, sg, ...
                sg / max(N.shiftGainDarkErr, eps));
        else
            V{end+1} = sprintf(['NOISE: no shift of the monitor window within +/-%d samples beats 0 ' ...
                'by more than its error. A channel_delay would not pay at this window.'], ...
                max(abs(T.shifts)));
        end
    end
end

function report(out)
    if ~isempty(out.edges)
        E = out.edges;
        fprintf('\n%d averages. Trigger box %g ns removed in the fitted times.\n', ...
                E.fl.nAverages, 1e3 * E.triggerJitter);
        fprintf('\n                          rise 10-90%%   t10      t50      t90    | fall 90-10%%  t10      t50      t90   (ns after the command)\n');
        printRow(E.fl);  printRow(E.ref);
        o = E.offset;
        fprintf('\nmonitor minus fluorescence (ns):  rise  10%% %+6.1f  50%% %+6.1f +/- %.1f  90%% %+6.1f\n', ...
                o.riseFit(1), o.riseFit(2), o.riseFitErr(2), o.riseFit(3));
        fprintf('                                  fall  10%% %+6.1f  50%% %+6.1f +/- %.1f  90%% %+6.1f\n', ...
                o.fallFit(1), o.fallFit(2), o.fallFitErr(2), o.fallFit(3));
        fprintf('  paired per average, 50%% from the data:  rise %+.1f +/- %.1f   fall %+.1f +/- %.1f\n', ...
                o.rise50Data, o.rise50DataErr, o.fall50Data, o.fall50DataErr);
        fprintf('  (measured with the current channel_delay = [%g %g] us applied)\n', ...
                E.info.currentDelay(1), E.info.currentDelay(2));
        for R = [E.fl, E.ref]
            if isfinite(R.levels.rangeUsed) && R.levels.rangeUsed > 0.9
                fprintf('NOTE %s uses %.0f%% of its input range: close to clipping.\n', R.label, 100 * R.levels.rangeUsed);
            end
        end
    end
    if ~isempty(out.noise)
        N = out.noise;  T = N.table;
        fprintf('\n%d repeat pairs over %d averages (%d outliers dropped).\n', N.nPairs, N.nAverages, N.outliers);
        fprintf(['levels: fluorescence %.5g V (dark %.5g, light fraction %.3f), monitor %.5g V ' ...
                 '(dark %.5g, light fraction %.3f)\n'], N.levels.S, N.levels.darkS, ...
                N.levels.lightFractionS, N.levels.R, N.levels.darkR, N.levels.lightFractionR);
        fprintf('window-to-window noise: fluorescence %.3g, monitor %.3g (relative)\n', ...
                N.relNoiseS, N.relNoiseR);
        fprintf('kappa (ln S on ln R per window): raw %.3f, dark-subtracted %.3f\n', ...
                N.kappaRaw, N.kappaDark);
        fprintf('\nspread of the half ratio, relative to no referencing (%.4g):\n', T.sigmaNone);
        fprintf('  shift (samples)   raw (today)          dark-subtracted\n');
        for k = 1:numel(T.shifts)
            tag = '';
            if T.shifts(k) == 0; tag = '   <- current'; end
            fprintf('  %+5d             %.4f +/- %.4f    %.4f +/- %.4f%s\n', T.shifts(k), ...
                    T.relRaw(k), T.relRawErr(k), T.relDark(k), T.relDarkErr(k), tag);
        end
    end
    fprintf('\n==================== WHAT IT SAYS ====================\n');
    for i = 1:numel(out.verdict)
        fprintf('%s\n', wrapText(out.verdict{i}, 92));
    end
    fprintf('Nothing was written. channel_delay lives in setupInfo.json (spcm section), edited by hand.\n');
    fprintf('======================================================\n');
end

function printRow(R)
    fprintf('  %-24s %6.1f +/- %-4.1f %7.1f %8.1f %8.1f  | %6.1f +/- %-4.1f %7.1f %8.1f %8.1f\n', R.label, ...
            R.rise.time, R.rise.timeErr, R.rise.t10, R.rise.t50, R.rise.t90, ...
            R.fall.time, R.fall.timeErr, R.fall.t10, R.fall.t50, R.fall.t90);
end

function s = wrapText(s, n)
    words = strsplit(s, ' ');
    lines = {};  cur = '';
    for i = 1:numel(words)
        if isempty(cur)
            cur = words{i};
        elseif numel(cur) + 1 + numel(words{i}) > n
            lines{end+1} = cur; %#ok<AGROW>
            cur = ['    ' words{i}];
        else
            cur = [cur ' ' words{i}]; %#ok<AGROW>
        end
    end
    lines{end+1} = cur;
    s = strjoin(lines, newline);
end

function plotIt(out)
    nP = double(~isempty(out.edges)) * 2 + double(~isempty(out.noise));
    if nP == 0; return; end
    figure('Name', 'calibPhotodiodePair', 'Color', 'w');
    p = 0;
    if ~isempty(out.edges)
        E = out.edges;
        kinds = {'rise', 'fall'};
        for i = 1:2
            p = p + 1;
            subplot(nP, 1, p); hold on;
            for R = [E.fl, E.ref]
                ed = R.(kinds{i});
                x = 1e3 * (ed.fitT - ed.commanded);
                yn = (ed.fitY - ed.lo) / (ed.hi - ed.lo);
                plot(x, yn, '.', 'MarkerSize', 4);
            end
            c = [E.fl.(kinds{i}).t50, E.ref.(kinds{i}).t50];
            xlim(mean(c) + [-1 1] * max(80, 3 * abs(diff(c))));
            ylim([-0.1 1.1]); grid on;
            legend({'fluorescence', 'monitor'}, 'Location', 'best');
            xlabel(sprintf('ns after the %s command', ternary(i == 1, 'on', 'off')));
            ylabel('normalised');
            title(sprintf('%s: monitor - fluorescence at 50%% = %+.1f ns', kinds{i}, diff(c)));
        end
    end
    if ~isempty(out.noise)
        T = out.noise.table;
        p = p + 1;
        subplot(nP, 1, p); hold on;
        errorbar(T.shifts, T.relRaw, T.relRawErr, 'o-');
        errorbar(T.shifts, T.relDark, T.relDarkErr, 's-');
        plot([min(T.shifts) max(T.shifts)], [1 1], 'k:');
        grid on;
        legend({'raw (today)', 'dark-subtracted', 'no referencing'}, 'Location', 'best');
        xlabel('monitor window shift (10 ns samples)');
        ylabel('half-ratio spread / unreferenced');
        title('does dividing by the monitor pay, and at which shift');
    end
end

% =========================================================================
% EDGE FIT (as laserRiseFall: one edge at a time)
% =========================================================================

function E = fitEdge(tau, y, e, Y, kind, tCmd, t50, win, J, modelChoice)
    isRise = strcmp(kind, 'rise');
    m = tau >= win(1) & tau <= win(2) & isfinite(y);
    t = tau(m);  yy = y(m);  ee = e(m);
    Yw = Y(m, :);
    span = win(2) - win(1);
    pre  = t < t50 - max(0.1, 0.25 * (t50 - win(1)));
    post = t > t50 + 0.5 * (win(2) - t50);
    yPre = nmedian(yy(pre));  yPost = nmedian(yy(post));
    if isRise; lo = yPre; hi = yPost; else; hi = yPre; lo = yPost; end
    E = struct('kind', kind, 'commanded', tCmd, 'window', win);

    sw = max(1, round(0.005 / median(diff(t))));
    E.data = dataCrossings(t, Yw, yy, sw, lo, hi, isRise, tCmd, pre, post);

    wr = E.data.time / 1e3;
    if ~(isfinite(wr) && wr > 0); wr = 0.05; end
    wr = min(max(wr, 0.002), 0.5 * span);
    d = hi - lo;
    hStep = [1e-4*abs(d), 1e-4*abs(d), 1e-5, 1e-5, 1e-5, 1e-3, 1e-2];

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
    sgl = best;  sgl.n = numel(t);  sgl.k = 5;

    dbl = [];
    if ~strcmp(modelChoice, 'single')
        ps = sgl.p;
        st2 = [0.05 10; -0.05 10; 0.1 50; 0.03 300];
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
    [C, free] = covariance(fit.J, p, win);
    C = C * max(1, E.chi2nu);
    E.params = p;
    E.lo = p(1);  E.hi = p(2);
    E.tau_ns = 1e3 * p(4);  E.sigma_ns = 1e3 * p(5);
    E.laser      = derivedTimes(p, C, free, tCmd, 0, win);
    E.asRecorded = derivedTimes(p, C, free, tCmd, J, win);
    E.time = E.laser.time;   E.timeErr = E.laser.timeErr;
    E.t10 = E.laser.t(1);    E.t50 = E.laser.t(2);   E.t90 = E.laser.t(3);
    E.fitT = t;  E.fitY = yy;  E.fitE = ee;
    E.fitCurve = edgeModel(t, p, isRise, J);
end

function D = dataCrossings(t, Yw, y, sw, lo, hi, isRise, tCmd, pre, post)
    % 10/50/90% completion off the data: on the mean, and per average (kept, so the two
    % photodiodes can be differenced average by average).
    D = struct('t', nan(1, 3), 'tErr', nan(1, 3), 'time', NaN, 'timeErr', NaN, 'per', []);
    c = crossings(t, smoothRun(y, sw), lo, hi, isRise);
    D.t = 1e3 * (c - tCmd);
    D.time = D.t(3) - D.t(1);
    A = size(Yw, 2);
    per = nan(A, 3);
    for a = 1:A
        ya = Yw(:, a);
        if ~any(isfinite(ya)); continue; end
        if isRise
            la = nmedian(ya(pre)); ha = nmedian(ya(post));
        else
            ha = nmedian(ya(pre)); la = nmedian(ya(post));
        end
        per(a, :) = crossings(t, smoothRun(ya, sw), la, ha, isRise);
    end
    D.per = per;
    ok = all(isfinite(per), 2);
    if nnz(ok) >= 3
        D.tErr = 1e3 * std(per(ok, :), 0, 1) / sqrt(nnz(ok));
        D.timeErr = 1e3 * std(per(ok, 3) - per(ok, 1)) / sqrt(nnz(ok));
    end
end

function c = crossings(t, y, lo, hi, isRise)
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

function y = edgeModel(t, p, isRise, J)
    G = shapeG(t, p, J);
    if isRise
        y = p(1) + (p(2) - p(1)) * G;
    else
        y = p(1) + (p(2) - p(1)) * (1 - G);
    end
end

function G = shapeG(t, p, J)
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
    x = t - t0;
    if sg < 1e-6 * tau
        s = (x > 0) .* (1 - exp(-max(x, 0) / tau));
        return
    end
    s = 0.5 * erfc(-x / (sg*sqrt(2))) - qTerm(x, tau, sg);
end

function I = stepInt(t, t0, tau, sg)
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
    z = x / sg - sg / tau;
    Q = zeros(size(x));
    pos = z >= 0;
    Q(pos)  = exp(-x(pos)/tau + sg^2/(2*tau^2)) .* 0.5 .* erfc(-z(pos)/sqrt(2));
    Q(~pos) = 0.5 * exp(-x(~pos).^2 / (2*sg^2)) .* erfcx(-z(~pos)/sqrt(2));
end

function q = offsetAt(frac, p, J)
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
    p(4) = min(max(p(4), 1e-4), 5);
    p(5) = min(max(p(5), 1e-4), 2);
    if numel(p) >= 7
        p(6) = min(max(p(6), -0.5), 0.95);
        p(7) = min(max(p(7), 2), 1e4);
    end
end

function [C, free] = covariance(Jac, p, win)
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
% SMALL HELPERS
% =========================================================================

function tc = crossTime(t, y, level, after, direction)
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
    if n <= 1; ys = y; return; end
    k = ones(n, 1);
    ok = isfinite(y);
    y0 = y;  y0(~ok) = 0;
    ys = conv(y0, k, 'same') ./ max(conv(double(ok), k, 'same'), 1);
    ys(~ok) = NaN;
end

function s = robustNoise(y)
    d = diff(y(isfinite(y)));
    s = 1.4826 * median(abs(d - median(d))) / sqrt(2);
    if ~(s > 0); s = std(d) / sqrt(2); end
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
        if ~isfield(o, known{i}) || (isempty(o.(known{i})) && ~any(strcmp(known{i}, {'replay', 'window', 'laserOn', 'nShifts'})))
            o.(known{i}) = d.(known{i});
        end
    end
end
