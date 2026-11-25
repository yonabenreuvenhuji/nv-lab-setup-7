function results = calibLastDelay(opts)
%CALIBLASTDELAY The gap between the microwave ending and the detection window opening, FITTED.
%
%   results = calibLastDelay()
%   results = calibLastDelay(opts)
%   results = calibLastDelay(struct('data', q))     % refit a saved result, no hardware
%
% LEVEL 2 CALIBRATION, run before calibDetectionWindow (see level2.m for why that order).
%
% ===================================================================
% WHAT IS MEASURED
% ===================================================================
%
% Sequence (AWGExpLastDelayDuration): trigger, saturating MW burst, GAP, then laser + detector.
% The MW switch stays open through the gap, so what ends at the start of the gap is the AWG
% waveform. While the gap is shorter than the time the microwave takes to actually stop at
% the NVs, MW overlaps the start of the detection window and costs contrast. Past that point
% the contrast no longer depends on the gap.
%
% ===================================================================
% WHY A FIT, AND WHAT MODEL
% ===================================================================
%
% The old analysis looked for a "knee": the first gap after which every point sat within two
% sigma of a plateau taken from the last third of the scan. On the 2026-11-03 run that gave
% 0.705 us, and it was wrong for two reasons that compound:
%
%   1. The curve has TWO parts. A sharp step near 0.1 us (the microwave stopping), then a slow,
%      straight creep of about +3.5% of the contrast per microsecond that does not flatten
%      anywhere in the scan. The "plateau" was taken on the creep, so it was never flat, and
%      the knee was simply where the creep crossed an arbitrary line.
%   2. Its noise estimate was the scatter of that sloping tail about a constant, which
%      includes the slope. It came out twice the true point-to-point noise.
%
% So this fits the two parts separately:
%
%       C(g) = background(g) + h * Phi((g - g0) / w)
%
%   Phi        the normal CDF: a smoothed step of height h, centre g0, width w
%   background a polynomial in g whose degree (0..3) is chosen by BIC. Degree 0 is a plain
%              flat plateau; the data decide whether a slope is needed.
%
% lastDelay is then where the STEP is finished, not where the background stops rising: the
% 99% point of the step, tc = g0 + 2.326 w, taken at its one-sided 95% upper bound so a noisy
% run errs long rather than short. The bound comes from a profile chi-square over (tc, w), so
% it stays honest when the edge is sharper than the scan step and w is not resolved.
%
% ONLY THE SHORT-GAP REGION IS FITTED (opts.stepWindow, default 0.6 us), and the sequence is
% BALANCED (opts.balancedSequence, default 1). Both learned the hard way on 2026-11-04:
%
%   The long scans (0-150 us) showed contrast rising ~0.6% with the gap, identically ON and OFF
%   resonance, with the MW off, with the AWG never triggered, with the switches mirrored. The
%   cause was the sequence: Experiment's balancedSequence defaults to 0, so the reference window
%   came straight after continuous laser while the signal window followed TRIG + burst + gap of
%   dark. Fluorescence right after a dark period is lower than under steady light -- by 0.14%
%   after ~1 us and ~0.7% after 100 us here (NV charge state or laser turn-on) -- and that read as
%   "contrast". Balanced, both windows follow the same dark and it cancels. Rabi, Ramsey and Echo
%   in the run script and levels 3-4 already set balancedSequence = 1; T1's contrast is |0> - |1>,
%   both after the same dark, so it cancels there too.
%
%   A step + cubic fitted over the whole 150 us once returned a "step" 551 us BEFORE the scan
%   with a height of 10^11 %, so the fit is confined to the step and checked for sanity.
%
% Class bugs fixed the same day (AWGExpLastDelayDuration): the AWG sequence alternated zeros and
% burst, so every other repeat had no microwave (half the contrast -- confirmed: the contrast
% doubled once fixed); the balanced reference used tau(end) instead of piTime for its dark.
%
% ===================================================================
% ERROR BARS
% ===================================================================
%
% The contrast is formed PER AVERAGE, 1 - S1./S2, and the error is the standard error of that
% across averages. Laser and focus drift between averages are common to signal and reference,
% so they cancel in the ratio; propagating S1 and S2 errors separately (the old way) charged
% that drift as noise and inflated the error ~20x. With fewer than 3 averages, or on a replay
% of a result that has no per-average data, the noise is estimated from the point-to-point
% scatter instead (robust to the slope).
%
% ===================================================================
% OPTIONS (struct fields, all optional)
%   data             a previous result struct (gaps + contrast, optionally contrastPerAvg):
%                    refit it, touch no hardware
%   jsonLocation     path to expProperties.json for loadJsonParams
%   gaps             gaps to scan, us (default 0.005:0.01:0.5 -- the step is at ~0.1 us)
%   detectionDuration  us (default 0.2; level2 passes a few us)
%   laserInitializationDuration  us (default detectionDuration + 5*tauPol)
%   mwBurst          saturating pulse length, us (default 1)
%   repeats, averages, frequency (MHz, 2870), amplitude (dBm, -10), fixDelays (1)
%   completion       fraction of the step that counts as finished (default 0.99)
%   confidence       one-sided confidence on the completion point (default 0.95)
%   maxBackgroundDegree  default 2 (inside the step window)
%   stepWindow       only gaps up to this are fitted, us (default 0.6)
%   analyse          false: measure and return the curve, no fit, no plot (calibDarkTime)
%   mwOff, gapSwitch, refSwitch, noTrigger   passed to the class: see AWGExpLastDelayDuration
%   minLastDelay     what to recommend when there is NO step (default 0.05 us): the MW is gone
%                    before the shortest gap, and 50 ns is a margin for delay drift that costs
%                    nothing (0.1% of a repeat)
%   balancedSequence  default 1: the reference window gets the same dark as the signal window.
%                     0 reproduces the old (wrong for this purpose) unbalanced scan.
%
% OUTPUT (adds to the old fields gaps/contrast/contrastErr/nIter/opts)
%   lastDelay        recommended value, us (NaN if no step could be located)
%   fit              struct: g0, w, h, tc, tcUpper, degree, coef, chi2nu, bic, stepSignificance
%   toWrite          struct('lastDelay', value) for expPropsWrite, or empty
%
% Writes nothing. level2.m cell 2.3 writes, via expPropsWrite.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('data', [], 'jsonLocation', '', ...
                        'gaps', 0.005:0.01:0.5, 'detectionDuration', 0.2, ...
                        'mwBurst', 1, 'repeats', 1000, 'averages', 20, ...
                        'frequency', 2870, 'amplitude', -10, 'fixDelays', 1, ...
                        'nChannels', 2, 'AWGorSRSswitch', 0, ...
                        'GIMeas', 0, 'timeDelay', 0, 'amplitudeDiff', [0 0], ...
                        'tauPol', 6.3, 'laserInitializationDuration', [], ...
                        'completion', 0.99, 'confidence', 0.95, ...
                        'maxBackgroundDegree', 2, 'stepWindow', 0.6, 'analyse', true, ...
                        'mwOff', [], 'gapSwitch', [], 'refSwitch', [], 'noTrigger', [], ...
        'balancedSequence', 1, 'minLastDelay', 0.05, 'plot', true));

    % ---- where the calibration lives ------------------------------------
    % Tolerant on purpose: a replay needs no json, and the offline tests run with no Setup 7
    % tree on the path.
    if isempty(opts.jsonLocation)
        try
            opts.jsonLocation = setup7Paths().json;
        catch
            opts.jsonLocation = '';
        end
    end

    if ~isempty(opts.data)
        % ================= REPLAY: refit saved data ======================
        d = opts.data;
        g  = d.gaps(:)';
        Ck = [];
        if isfield(d, 'contrastPerAvg') && ~isempty(d.contrastPerAvg)
            Ck = double(d.contrastPerAvg);
        end
        C0 = d.contrast(:);
        nIter = NaN; if isfield(d, 'nIter'); nIter = d.nIter; end
        if isfield(d, 'opts')
            % Keep the run's own settings for the printout, but not its data/plot flags.
            keep = {'detectionDuration', 'laserInitializationDuration', 'mwBurst', ...
                    'repeats', 'averages', 'frequency', 'amplitude'};
            for i = 1:numel(keep)
                if isfield(d.opts, keep{i}); opts.(keep{i}) = d.opts.(keep{i}); end
            end
        end
        fprintf('\n=== calibLastDelay (REPLAY of saved data, no hardware) ===\n');
    else
        % ================= MEASURE ======================================
        if isempty(opts.laserInitializationDuration)
            opts.laserInitializationDuration = opts.detectionDuration + 5*opts.tauPol;
        end
        g = opts.gaps(:)';
        fprintf('\n=== calibLastDelay ===\n');
        fprintf('%d gaps from %g to %g us, detection window %g us, MW burst %g us, %d x %d\n', ...
                numel(g), g(1), g(end), opts.detectionDuration, opts.mwBurst, ...
                opts.repeats, opts.averages);
        fprintf('laser initialization %g us.  ZERO FIELD required (%g MHz); the Helmholtz is not touched.\n', ...
                opts.laserInitializationDuration, opts.frequency);
        if ~opts.balancedSequence
            fprintf(['UNBALANCED sequence: the reference window follows laser light, the signal window\n' ...
                     '   follows dark. The result mixes in the fluorescence-after-dark effect.\n']);
        end

        c = AWGExpLastDelayDuration;
        if ~isempty(opts.jsonLocation)
            c.loadJsonParams('all', opts.jsonLocation);
        end
        % Setup 7 hardware configuration, set explicitly: AWGExpLastDelayDuration never
        % overrides ExperimentAWG's AWGorSRSswitch = 1 (a channel Setup 7 does not have), and
        % nChannels defaults to 1 where Setup 7 runs two.
        c.AWGorSRSswitch             = opts.AWGorSRSswitch;
        c.nChannels                  = opts.nChannels;
        c.GIMeas                     = opts.GIMeas;
        c.timeDelay                  = opts.timeDelay;
        c.amplitudeDiff              = opts.amplitudeDiff;

        c.tau                        = g;
        c.repeats                    = opts.repeats;
        c.averages                   = opts.averages;
        c.detectionDuration          = opts.detectionDuration;
        % Must equal the detection duration: the digitizer acquires both windows with one set
        % of segment settings (Experiment.prepareInternal rejects anything else).
        c.referenceDetectionDuration = opts.detectionDuration;
        c.laserInitializationDuration = opts.laserInitializationDuration;
        c.frequency                  = opts.frequency;
        c.amplitude                  = opts.amplitude + zeros(1, opts.nChannels);
        c.piTime                     = opts.mwBurst;
        c.isTracking                 = 0;
        c.fixDelays                  = opts.fixDelays;
        % BALANCED: the reference window follows the same dark (TRIG + burst + gap) as the signal
        % window. Experiment's default is 0, and with 0 the reference follows continuous laser --
        % then the lower fluorescence after a dark period reads as a slow "contrast" rise with the
        % gap (0.14% -> 0.7% over 0-150 us on 2026-11-04), which is what all the long scans saw.
        c.balancedSequence           = opts.balancedSequence;
        % Isolation switches, only when asked for (and only if the class has them).
        iso = {'mwOff', 'gapSwitch', 'refSwitch', 'noTrigger'};
        for i = 1:numel(iso)
            if ~isempty(opts.(iso{i}))
                if ~isprop(c, iso{i})
                    error(['calibLastDelay: this AWGExpLastDelayDuration has no "%s" property. ' ...
                           'Update the class from 2026-11-04 or later.'], iso{i});
                end
                c.(iso{i}) = opts.(iso{i});
            end
        end

        % ---- PRE-FLIGHT: will the delay compensation survive these gaps? ---
        % fixSequenceDelays turns each pulse of length D into D + onDelay - offDelay; the gap is
        % a pulse on the MW channels, so a gap below max(offDelay - onDelay) goes negative and
        % Pulse/set.duration throws. Drop those gaps with one sentence instead of 200 lines.
        if opts.fixDelays
            [gSafe, gBad, limit, limCh] = safeGaps(g);
            if ~isempty(gBad)
                fprintf(['\nPRE-FLIGHT: %d of %d gaps are below %.4f us = offDelay - onDelay on "%s".\n' ...
                         'Dropped: %s\n'], numel(gBad), numel(g), limit, limCh, numList(gBad));
                if isempty(gSafe)
                    error(['calibLastDelay: every gap is below the %.4f us floor set by "%s".\n' ...
                           'Run pgDelayAudit(%.4f) to see the arithmetic per channel. The MW OFF edge\n' ...
                           'is NV re-polarisation, not a switch, so an off delay stored for it is\n' ...
                           'almost certainly wrong; setting it to 0 fixes this.'], limit, limCh, min(g));
                end
                g = gSafe;
                c.tau = g;
            end
        end

        c.run;

        if isempty(c.signal); error('The experiment produced no data.'); end
        if ~any(isfinite(c.signal(:)) & c.signal(:) ~= 0)
            error(['calibLastDelay: every point is zero or NaN, so the run aborted rather than ' ...
                   'measured. Scroll up for the first error; if it was a negative Pulse duration, ' ...
                   'run pgDelayAudit.']);
        end
        nIter = c.currIter;
        S1 = reshape(c.signal(1, :, 1:nIter), [], nIter);
        S2 = reshape(c.signal(2, :, 1:nIter), [], nIter);
        n  = min([numel(g), size(S1, 1)]);
        g  = g(1:n); S1 = S1(1:n, :); S2 = S2(1:n, :);
        Ck = 1 - S1 ./ S2;                       % contrast PER AVERAGE -- drift cancels here
        Ck(S1 == 0 & S2 == 0) = NaN;             % averages that never ran
        C0 = nanMean2(Ck);
    end

    % ---- point errors -----------------------------------------------------
    g = g(:);
    [C, dC, noiseSrc] = pointErrors(C0, Ck);
    results = struct('gaps', g, 'contrast', C, 'contrastErr', dC, 'nIter', nIter, ...
                     'noiseSource', noiseSrc, 'opts', rmfieldSafe(opts, 'data'));
    if ~isempty(Ck); results.contrastPerAvg = single(Ck); end
    if ~opts.analyse; return; end            % calibDarkTime does its own analysis

    % ---- the fit ------------------------------------------------------------
    F = fitStep(g, C, dC, opts);
    results.fit = F;
    results.plateau = F.plateau;
    results.lastDelay = F.recommended;
    results.toWrite = struct();
    if isfinite(F.recommended); results.toWrite.lastDelay = F.recommended; end

    [sv, haveStored] = storedVal('lastDelay', opts.jsonLocation);
    printReport(g, C, dC, F, opts, noiseSrc, sv, haveStored);
    if opts.plot
        plotFit(g, C, dC, F, sv, haveStored);
    end
end

% ======================================================================
%  NOISE
% ======================================================================

function [C, dC, src] = pointErrors(C0, Ck)
    C = C0(:);
    if ~isempty(Ck) && size(Ck, 2) >= 3
        K  = sum(isfinite(Ck), 2);
        dC = nanStd2(Ck) ./ sqrt(max(K, 1));
        src = sprintf('standard error over %d averages of the per-average contrast', ...
                      round(median(K)));
        % One wild average must not set a point's weight to near-infinite either way.
        dC = max(dC, 0.25*median(dC));
    else
        % Point-to-point scatter, from first differences so a slope does not inflate it.
        dd  = diff(C);
        sig = 1.4826 * median(abs(dd - median(dd))) / sqrt(2);
        dC  = max(sig, eps) * ones(size(C));
        src = 'point-to-point scatter (no per-average data)';
    end
end

% ======================================================================
%  THE FIT
% ======================================================================

function F = fitStep(g, C, dC, opts)
% Erf step + polynomial background, by profile chi-square on a (tc, w) grid. For fixed
% (tc, w) every other parameter is linear, so each grid point is one weighted least-squares
% solve -- no start values, nothing to converge, no local minima in the linear part.
    z  = sqrt(2) * erfinv(2*opts.completion - 1);     % 2.326 for 99%
    dq = chi2inv1(2*opts.confidence - 1);             % one-sided bound: 2.706 for 95%
    % Only the step region. Past it the slow hardware rise (see header) is not a background a
    % low-order polynomial can follow, and fitting it is how the absurd fits happened.
    k = g <= opts.stepWindow;
    if sum(k) < 8
        F = emptyFit(); F.status = 'coarse'; F.nInWindow = sum(k); return
    end
    gAll = g;
    g = g(k); C = C(k); dC = dC(k);
    span = g(end) - g(1);
    dgMin = min(diff(g));
    gc = (g - mean(g)) / max(span, eps);              % centred, scaled: well-conditioned
    wv = 1 ./ dC.^2;
    N  = numel(g);

    wGrid  = exp(linspace(log(1e-3), log(max(span/3, 2e-3)), 60));
    tcStep = max(min(dgMin/5, 0.002), 5e-4);
    tcGrid = (g(1) : tcStep : g(end))';

    best = struct('bic', Inf);
    for deg = 0:opts.maxBackgroundDegree
        if N < deg + 6; break; end
        Pb = gc .^ (0:deg);                           % background basis
        chi = profileGrid(g, C, wv, Pb, tcGrid, wGrid, z);
        [cmin, k] = min(chi(:));
        if ~isfinite(cmin); continue; end
        kPar = (deg + 1) + 3;                         % background, h, tc, w
        bic  = cmin + kPar * log(N);
        % The same background with no step, for the significance test.
        [~, c2none] = wls(Pb, C, wv);
        if bic < best.bic
            [it, iw] = ind2sub(size(chi), k);
            best = struct('bic', bic, 'deg', deg, 'chi', chi, 'cmin', cmin, ...
                          'it', it, 'iw', iw, 'kPar', kPar, 'c2none', c2none);
        end
    end
    if ~isfinite(best.bic)
        F = emptyFit(); return
    end

    % ---- best point, refined off the grid -----------------------------------
    deg = best.deg;  Pb = gc .^ (0:deg);
    tc  = tcGrid(best.it);  w = wGrid(best.iw);
    obj = @(p) stepChi(p, g, C, wv, Pb, z);
    try
        p = fminsearch(obj, [tc, log(w)], optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6));
        if obj(p) <= best.cmin; tc = p(1); w = exp(p(2)); end
    catch
    end
    g0 = tc - z*w;
    st = 0.5 * (1 + erf((g - g0) / (sqrt(2)*w)));
    [x, cmin] = wls([Pb, st], C, wv);
    dof    = max(N - best.kPar, 1);
    chi2nu = cmin / dof;
    s2     = max(chi2nu, 1);                          % scale errors up, never down

    % ---- profile over tc: the upper bound on the completion point -------------
    prof = min(best.chi, [], 2);
    thr  = min(prof) + dq * s2;
    in   = prof <= thr;
    % Connected region around the best point only.
    iu = best.it; while iu < numel(tcGrid) && in(iu+1); iu = iu + 1; end
    il = best.it; while il > 1 && in(il-1); il = il - 1; end
    tcUpper = tcGrid(iu);  tcLower = tcGrid(il);
    % Width: profile over w for the error printout.
    profW = min(best.chi, [], 1);
    inW = profW <= min(profW) + 1*s2;
    wHi = wGrid(find(inW, 1, 'last'));
    wResolved = wGrid(find(inW, 1, 'first')) > wGrid(1) * 1.01;

    % Step significance: chi2 drop from adding the step to the same background. And whether
    % the contrast is flat at all -- "no step" is only an answer when nothing rises, because a
    % step that lies past the end of the scan shows up as a rising background, not as a step.
    dchiStep = (best.c2none - cmin) / s2;
    [~, c2flat] = wls(ones(N, 1), C, wv);
    risesAnyway = (c2flat - best.c2none) / s2 > 25;

    bgFun   = @(gg) (((gg(:) - mean(g)) / max(span, eps)) .^ (0:deg)) * x(1:end-1);
    stepFun = @(gg) x(end) * 0.5 * (1 + erf((gg(:) - g0) / (sqrt(2)*w)));
    model   = @(gg) bgFun(gg) + stepFun(gg);
    % Background slope at the completion point, in contrast per us.
    hgg = 1e-4;
    slope = (bgFun(tc + hgg) - bgFun(tc - hgg)) / (2*hgg);
    plateau = model(tc);

    F = struct('g0', g0, 'w', w, 'wUpper', wHi, 'wResolved', wResolved, 'h', x(end), ...
               'tc', tc, 'tcLower', tcLower, 'tcUpper', tcUpper, 'degree', deg, ...
               'coef', x, 'chi2nu', chi2nu, 'bic', best.bic, 'stepDchi2', dchiStep, ...
               'slope', slope, 'plateau', plateau, 'model', model, 'background', bgFun, ...
               'step', stepFun, 'z', z, 'completion', opts.completion, ...
               'recommended', NaN, 'status', '');

    F.gMax = g(end);  F.nInWindow = numel(g);  F.gAllEnd = gAll(end);

    % ---- sanity: a step that is not a step is no answer -----------------------
    range = max(C) - min(C) + 3*median(dC);
    if g0 < g(1) - 0.5*span || g0 > g(end) || w > span/3 || x(end) > 3*range
        F.status = 'nofit'; F.recommended = NaN; return
    end

    % ---- the recommendation, with the guards --------------------------------
    if dchiStep < 25 && (risesAnyway || g(end) < 0.4)
        % Rises without a separable step, or flat over a window too short to tell "MW already
        % gone" from "MW not gone yet" -- both read as flat. 0.4 us settles it: no switch or
        % amplifier on this setup takes that long (the step is at ~0.1 us).
        F.status = 'unclear';
    elseif dchiStep < 25
        F.status = 'nostep';
        % Flat everywhere: the MW is gone before the shortest gap. Recommend a small margin, not
        % the shortest gap itself, and show the no-step fit rather than the rejected step.
        F.recommended = max(ceil(g(1)/0.005)*0.005, opts.minLastDelay);
        [xb, ~] = wls(Pb, C, wv);
        F.background = @(gg) (((gg(:) - mean(g)) / max(span, eps)) .^ (0:deg)) * xb;
        F.model = F.background;  F.h = 0;
        F.tcLower = NaN; F.tcUpper = NaN; F.tc = NaN;
    elseif tcUpper >= g(end) - 2*dgMin
        F.status = 'beyond';            % still stepping at the end: cannot place it
    elseif tcLower <= g(1)
        F.status = 'atstart';           % step at the very start: MW gone before g(1)
        F.recommended = ceil(max(tcUpper, g(1))/0.005)*0.005;
    else
        F.status = 'ok';
        F.recommended = ceil(tcUpper/0.005)*0.005;
    end
end

function chi = profileGrid(g, C, wv, Pb, tcGrid, wGrid, z)
% chi-square at every (tc, w), with the background and the step height solved exactly.
% For fixed w the background is projected out once and every tc becomes one column, so the
% whole tc axis is a handful of matrix products instead of thousands of separate solves.
    sw  = sqrt(wv(:));
    Pw  = Pb .* sw;                     % weighted background basis
    yw  = C(:) .* sw;
    [Q, ~] = qr(Pw, 0);                 % orthonormal basis of the background space
    ry  = yw - Q*(Q'*yw);               % data with the background removed
    rr  = ry' * ry;
    chi = Inf(numel(tcGrid), numel(wGrid));
    for iw = 1:numel(wGrid)
        w  = wGrid(iw);
        G0 = tcGrid(:)' - z*w;                                   % 1 x nTc step centres
        S  = 0.5 * (1 + erf((g(:) - G0) / (sqrt(2)*w)));        % N x nTc step columns
        Sw = S .* sw;
        Sr = Sw - Q*(Q'*Sw);                                     % step with background removed
        num = Sr' * ry;  den = sum(Sr.^2, 1)';
        h   = num ./ max(den, eps);
        c2  = rr - num.^2 ./ max(den, eps);
        c2(h <= 0 | den < 1e-12) = Inf;                          % a step DOWN is not the MW clearing
        chi(:, iw) = c2;
    end
end

function c2 = stepChi(p, g, C, wv, Pb, z)
    w = exp(p(2));
    if w < 1e-4 || w > (g(end)-g(1)); c2 = Inf; return; end
    st = 0.5 * (1 + erf((g - (p(1) - z*w)) / (sqrt(2)*w)));
    [x, c2] = wls([Pb, st], C, wv);
    if x(end) <= 0; c2 = Inf; end
end

function [x, c2] = wls(A, y, wv)
    sw = sqrt(wv(:));
    Aw = A .* sw;  yw = y(:) .* sw;
    x  = Aw \ yw;
    r  = yw - Aw*x;
    c2 = r' * r;
end

function q = chi2inv1(p)
% Inverse chi-square CDF for ONE degree of freedom, without the Statistics Toolbox.
    q = 2 * erfinv(p)^2;
end

function F = emptyFit()
    F = struct('g0', NaN, 'w', NaN, 'wUpper', NaN, 'wResolved', false, 'h', NaN, ...
               'tc', NaN, 'tcLower', NaN, 'tcUpper', NaN, 'degree', NaN, 'coef', [], ...
               'chi2nu', NaN, 'bic', NaN, 'stepDchi2', 0, 'slope', NaN, 'plateau', NaN, ...
               'model', @(gg) NaN(size(gg)), 'background', @(gg) NaN(size(gg)), ...
               'step', @(gg) NaN(size(gg)), 'z', NaN, 'completion', NaN, ...
               'recommended', NaN, 'status', 'nofit');
end

% ======================================================================
%  REPORT
% ======================================================================

function printReport(g, C, dC, F, opts, noiseSrc, sv, haveStored)
    fprintf('\n--- fit:  C(g) = background + h * step((g - g0)/w) ---\n');
    fprintf('   %d gaps, %.3f to %.3f us.  Errors: %s (median %.4f%%).\n', ...
            numel(g), g(1), g(end), noiseSrc, 100*median(dC));
    switch F.status
        case 'coarse'
            fprintf(['   NO FIT: only %d gaps at or below %.2f us, where the step is. This scan does not\n' ...
                     '   resolve it. For lastDelay use fine gaps, e.g. 0.005:0.01:0.5 (the default). For the\n' ...
                     '   long-gap behaviour use calibDarkTime. Nothing to write.\n\n'], F.nInWindow, opts.stepWindow);
            return
        case 'nofit'
            fprintf(['   NO FIT: the step model does not describe the data below %.2f us (the best fit put\n' ...
                     '   the step outside the window, or made it wider or taller than the data allow).\n' ...
                     '   Look at the curve. Nothing to write.\n\n'], opts.stepWindow);
            return
    end
    fprintf('   fitted: the %d gaps up to %.2f us (opts.stepWindow); the rest is shown, not fitted.\n', ...
            F.nInWindow, F.gMax);
    bgName = {'flat', 'linear', 'quadratic', 'cubic'};
    fprintf('   background %-9s (degree chosen by BIC)   chi2_nu = %.2f\n', ...
            bgName{F.degree+1}, F.chi2nu);
    if F.chi2nu > 2
        fprintf('   chi2_nu above 2: the model misses structure or the errors are optimistic; the\n');
        fprintf('   bound below has been widened by sqrt(chi2_nu) to compensate.\n');
    end
    if strcmp(F.status, 'nostep')
        fprintf('   step      none: adding one improves the fit by only %.1f sigma\n', sqrt(max(F.stepDchi2, 0)));
    end
    if ~strcmp(F.status, 'nostep')
    fprintf('   step      height %.4f%% (%.1f%% of the contrast), %.0f sigma\n', ...
            100*F.h, 100*F.h/max(F.plateau, eps), sqrt(max(F.stepDchi2, 0)));
    end
    if strcmp(F.status, 'nostep')
        % nothing to say about a centre or width
    elseif F.wResolved
        fprintf('             centre %.4f us, width %.4f us\n', F.g0, F.w);
    else
        fprintf(['             centre %.4f us, width < %.4f us -- sharper than the scan step, so\n' ...
                 '             not resolved; the bound below allows for any width up to that.\n'], ...
                F.g0, F.wUpper);
    end

    fprintf('\n--- lastDelay ---\n');
    switch F.status
        case 'ok'
            fprintf('   MW gone (%.0f%% of the step) at %.4f us; %.0f%% upper bound %.4f us\n', ...
                    100*opts.completion, F.tc, 100*opts.confidence, F.tcUpper);
            fprintf('   ->  lastDelay = %.3f us\n', F.recommended);
        case 'atstart'
            fprintf('   The step sits at the very start of the scan, so the MW is gone by %.4f us.\n', F.tcUpper);
            fprintf('   ->  lastDelay = %.3f us\n', F.recommended);
        case 'nostep'
            fprintf(['   No step: the contrast does not depend on the gap anywhere from %.3f us, so the\n' ...
                     '   MW is already gone when the window opens, even at the shortest gap.\n' ...
                     '   ->  lastDelay = %.3f us  (opts.minLastDelay: a margin for delay drift; there\n' ...
                     '       is nothing to gain below it)\n'], g(1), F.recommended);
        case 'beyond'
            fprintf(['   The step is still going at the end of the scan (%.3f us), so it cannot be\n' ...
                     '   placed. Extend opts.gaps past %.2f us and re-run. Nothing to write.\n'], ...
                    g(end), 2*g(end));
        case 'unclear'
            if F.gMax < 0.4
                fprintf(['   No step inside a window that stops at %.3f us. Flat data that short cannot tell\n' ...
                         '   "MW already gone" from "MW not gone yet". Scan to 0.5 us. Nothing to write.\n'], F.gMax);
            else
                fprintf(['   The contrast rises through the window but no step can be separated from the\n' ...
                         '   rise. Look at the curve; nothing to write.\n']);
            end
    end
    if haveStored && isfinite(F.recommended)
        fprintf('   stored %.3f us', sv);
        if sv >= F.tcUpper
            fprintf(' -- already past the step; shortening saves %.0f ns of dead time per half-repeat.\n', ...
                    1e3*(sv - F.recommended));
        else
            fprintf(' -- SHORT of the step: the MW is still on when the window opens.\n');
        end
    end

    % ---- the background, said once ----
    if F.degree >= 1 && isfinite(F.slope)
        rel = F.slope / max(F.plateau, eps);
        fprintf('\n   Background slope in the window: %+.4f%% per us (%+.1f%% of the contrast per us).\n', ...
                100*F.slope, 100*rel);
        if opts.balancedSequence
            fprintf('   Small and not the MW clearing; a balanced scan should be close to flat.\n');
        else
            fprintf(['   UNBALANCED scan: this is the fluorescence-after-dark effect (see the header), not\n' ...
                     '   the MW. Use the default balancedSequence = 1.\n']);
        end
    end
    fprintf('   Nothing was written (level2.m cell 2.3 writes).\n\n');
end

% ======================================================================
%  PLOT
% ======================================================================

function plotFit(g, C, dC, F, sv, haveStored)
    COL.data = [0.15 0.42 0.69];  COL.fit = [0.10 0.10 0.10];
    COL.bg   = [0.55 0.55 0.55];  COL.rec = [0.00 0.45 0.35];  COL.band = [0.86 0.94 0.91];
    COL.old  = [0.85 0.33 0.10];

    figure('Color', 'w', 'Position', [100 100 1100 420]);
    if isfinite(F.tcUpper)
        zEnd = min(g(end), max(4*F.tcUpper, F.tcUpper + 0.15));
    else
        zEnd = g(end);
    end
    gMax = g(end); if isfield(F, 'gMax'); gMax = F.gMax; end
    gg1 = linspace(g(1), min(zEnd, gMax), 600)';
    gg2 = linspace(g(1), gMax, 800)';           % the fit is drawn only where it was fitted

    for panel = 1:2
        ax = subplot(1, 2, panel); hold(ax, 'on'); tidy(ax);
        if panel == 1; k = g <= zEnd; gg = gg1; else; k = true(size(g)); gg = gg2; end
        yv = 100*[C(k) - dC(k); C(k) + dC(k)];
        yl = [min(yv) max(yv)]; yl = yl + 0.08*max(diff(yl), eps)*[-1 1];
        hP = [];
        if isfinite(F.tcUpper) && isfinite(F.tcLower)
            hP = patch(ax, [F.tcLower F.tcUpper F.tcUpper F.tcLower], yl([1 1 2 2]), COL.band, ...
                       'EdgeColor', 'none');
        end
        hD = errorbar(ax, g(k), 100*C(k), 100*dC(k), '.');
        set(hD, 'Color', COL.data, 'MarkerSize', 12, 'LineWidth', 1);
        try; set(hD, 'CapSize', 0); catch; end
        hF = plot(ax, gg, 100*F.model(gg), '-', 'Color', COL.fit, 'LineWidth', 1.4);
        hB = plot(ax, gg, 100*(F.background(gg) + F.h), '--', 'Color', COL.bg, 'LineWidth', 1.1);
        hs = [hD hF hB]; ls = {'measured', 'fit', 'background + full step'};
        if ~isempty(hP); hs(end+1) = hP; ls{end+1} = '95% range of the 99% point'; end
        if isfinite(F.recommended)
            hR = plot(ax, F.recommended*[1 1], yl, '-', 'Color', COL.rec, 'LineWidth', 1.8);
            hs(end+1) = hR; ls{end+1} = sprintf('lastDelay %.3f us', F.recommended); %#ok<AGROW>
        end
        if haveStored && sv <= xr0(panel, g, zEnd)
            hO = plot(ax, sv*[1 1], yl, ':', 'Color', COL.old, 'LineWidth', 1.4);
            hs(end+1) = hO; ls{end+1} = sprintf('stored %.3f us', sv); %#ok<AGROW>
        end
        if panel == 1; xr = [g(1) zEnd]; else; xr = [g(1) g(end)]; end
        ylim(ax, yl); xlim(ax, xr);
        xlabel(ax, 'MW-to-detection gap (us)'); ylabel(ax, 'contrast (%)');
        if panel == 2
            legend(ax, hs, ls, 'Location', 'southeast', 'Box', 'off');
        end
        if panel == 1
            switch F.status
                case {'ok', 'atstart'}
                    title(ax, sprintf('MW gone by %.3f us (%.0f%% of step)  ->  lastDelay = %.3f us', ...
                          F.tc, 100*F.completion, F.recommended), 'FontWeight', 'bold');
                case 'nostep'
                    title(ax, sprintf('no step: no MW effect at any gap  ->  lastDelay = %.3f us', ...
                          F.recommended), 'FontWeight', 'bold');
                otherwise
                    title(ax, 'no step separable inside the scan -- extend the gaps', 'FontWeight', 'bold');
            end
        else
            title(ax, sprintf('whole scan  (fit only up to %.2f us)', gMax), 'FontWeight', 'bold');
        end
    end
end

function x = xr0(panel, g, zEnd)
    if panel == 1; x = zEnd; else; x = g(end); end
end

function tidy(ax)
    set(ax, 'Box', 'off', 'TickDir', 'out', 'Layer', 'top', ...
            'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.12, 'FontSize', 9);
end

% ======================================================================
%  HELPERS
% ======================================================================

function [good, bad, limit, limCh] = safeGaps(g)
% Gaps that survive fixSequenceDelays, and the floor that decides it: max over channels of
% (offDelay - onDelay). Fails OPEN -- if the delays cannot be read, nothing is dropped.
    good = g; bad = []; limit = 0; limCh = '';
    try
        j = JsonInfoReader.getJson();
    catch
        return
    end
    D = struct();
    try
        for L = asList(j.lasers)
            LL = L{1};
            if isstruct(LL) && isfield(LL, 'switch')
                D.(vname(LL.switch.switchChannelName)) = ...
                    [LL.switch.onDelay, LL.switch.offDelay];
            end
        end
    catch
    end
    try
        f = fieldnames(j.microwave);
        for i = 1:numel(f)
            m = j.microwave.(f{i});
            D.(vname(m.switchChannelName)) = [m.onDelay, m.offDelay];
        end
    catch
    end
    names = fieldnames(D);
    for i = 1:numel(names)
        d = D.(names{i})(2) - D.(names{i})(1);
        if d > limit; limit = d; limCh = names{i}; end
    end
    if limit <= 0; return; end
    keep = g > limit;
    good = g(keep); bad = g(~keep);
end

function n = vname(s)
    n = regexprep(char(s), '\W', '_');
    if isempty(regexp(n, '^[A-Za-z]', 'once')); n = ['x' n]; end
end

function c = asList(x)
% jsondecode gives a 1x1 STRUCT for a one-element array and a CELL for a heterogeneous one.
    if iscell(x); c = x; elseif isstruct(x); c = num2cell(x); else; c = {x}; end
end

function s = numList(v)
    p = cell(1, numel(v));
    for i = 1:numel(v); p{i} = sprintf('%.4f', v(i)); end
    s = strjoin(p, ', ');
end

function [v, ok] = storedVal(field, jsonLocation)
    v = NaN; ok = false;
    try
        if isempty(jsonLocation)
            j = JsonInfoReader.getJson();
        else
            j = jsondecode(fileread(jsonLocation));
        end
        if isfield(j, field) && isnumeric(j.(field)) && isscalar(j.(field))
            v = j.(field); ok = true;
        end
    catch
    end
end

function m = nanMean2(X)
% Row means ignoring NaN. Written out because Octave's std has no 'omitnan'.
    ok = isfinite(X); X(~ok) = 0;
    m = sum(X, 2) ./ max(sum(ok, 2), 1);
    m(sum(ok, 2) == 0) = NaN;
end

function s = nanStd2(X)
    ok = isfinite(X); n = sum(ok, 2);
    m = nanMean2(X); D = X - m; D(~ok) = 0;
    s = sqrt(sum(D.^2, 2) ./ max(n - 1, 1));
    s(n < 2) = NaN;
end

function s = rmfieldSafe(s, f)
    if isfield(s, f); s = rmfield(s, f); end
end

function o = withDefaults(o, d)
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(o, f{i}) || (isempty(o.(f{i})) && ~strcmp(f{i}, 'data'))
            o.(f{i}) = d.(f{i});
        end
    end
end
