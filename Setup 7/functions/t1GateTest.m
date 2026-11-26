function out = t1GateTest(opts)
%T1GATETEST Is the open microwave switch shortening T1? Measured three ways, no fit required.
%
%   t1GateTest                                  % the real thing, ~25 min at the defaults
%   t1GateTest(struct('selfTest', true))        % no hardware: grade the statistics
%   t1GateTest(struct('selfTest', true, 'leakRatePerUs', 2e-5))   % with a leak injected
%   t1GateTest(struct('repeats', 200, 'averages', 2))             % faster, less sensitive
%
% ===================================================================
% THE QUESTION
% ===================================================================
%
% AWGExpT1 held the microwave switch open for the WHOLE of tau. The AWG plays for 0.335 us
% after each trigger -- a 0.15 us zero pad, piTime = 0.035, and another 0.15 of pad -- and is
% idle after that, so at tau = 4900 us the switch was open for 4900 us where 0.3 is needed.
% echoGateTest bounded the gate as harmless out to 240 us of excess; this is twenty times that,
% and Ramsey's measured-harmless case was 4.05 us, which is twelve hundred times smaller.
%
% What makes it worth a dedicated test rather than an extrapolation is that on a T1 THE ARTEFACT
% HAS THE SHAPE OF THE SIGNAL. Anything leaking through an open switch removes contrast at some
% rate r per microsecond of open gate, so with the gate open for tau the measured contrast is
%
%       c(tau) = c0 * exp(-tau/T1) * exp(-r*tau)
%
% and all that does is add to the relaxation rate:  1/T1_measured = 1/T1_true + r.  There is no
% feature to notice, no kink, no bad residual -- just a T1 that is too short, with a small error
% bar. Echo's fixed gate at least put its worst excess at SHORT tau, against the direction of
% the decay, so it was separable in shape. Here it is not.
%
% And the mechanism is not exotic. An open switch passes broadband noise out of the amplifier
% chain, not only the carrier, and off-resonant noise driving spin flips is exactly the process
% T1 measures. The lab's own non-AWG ExpT1 closes the gate during the wait -- MW off, gate open
% only for piTime -- which is the strongest hint that leaving it open here was inherited from
% Ramsey (where the second pi/2 genuinely lands at the end of tau) rather than chosen.
%
% ===================================================================
% HOW YOU MEASURE A RATE THAT LOOKS LIKE THE SIGNAL
% ===================================================================
%
% You do not try to separate it inside one run. You vary the DOSE at fixed physics, and the
% patched AWGExpT1 has a knob for exactly that -- mwGateMode -- which changes when the switch is
% open without changing anything else:
%
%       arm   mwGateMode   gate-open time G        signature in log(c)
%       ----------------------------------------------------------------------------
%       A     'tau'        G = tau                 slope -r, the historical behaviour
%       B     'pulse'      G = hold (0.435 us)     the baseline, essentially G = 0
%       C     'all'        G = lastDelay + max(tau), CONSTANT -- offset, NO slope
%       D     'tau'        a repeat of A, for drift and the error-bar factor
%
% Because log c = log c0 - tau/T1 - r*G, taking RATIOS between arms cancels c0 and T1 exactly.
% Nothing about the relaxation has to be modelled, fitted, or assumed:
%
%       log(cA/cB) = -r*(tau - hold)       -> weighted straight line, SLOPE = -r
%       log(cC/cB) = -r*(Gc  - hold)       -> a CONSTANT; r = -intercept/Gc, and the slope of
%                                             this one must come back at ZERO
%
% Two independent handles on the same r, with different systematics: one lives in the slope
% against tau, the other in an offset with no tau dependence at all. A drift, a tracking step or
% a laser change cannot produce both with the same r -- that is the whole reason for a third arm
% rather than the obvious two.
%
% THE ARMS ARE TIMING-IDENTICAL, which is what makes this clean. With constantTime = true (this
% test forces it) every arm has the same total sequence duration at every tau: 'pulse' splits the
% wait into hold + (tau - hold) and 'all' only changes whether lastDelay carries MW. Same duty
% cycle, same laser load, same repetition rate, same everything -- the switch's open window is
% the single difference. That is rarer than it sounds and it is why this test can quote a rate
% rather than a bound on a bound.
%
% ===================================================================
% WHAT IT REPORTS
% ===================================================================
%
%   T1 per arm, model-free: the slope of log(contrast) against tau is -1/T1. Single-exponential
%   by assumption -- t1FitRobust is the place for a shape -- but the same assumption in all four
%   arms, so the COMPARISON survives it even where the absolute number does not.
%
%   r from the slope (A vs B) and r from the offset (C vs B), with errors, and whether they
%   agree. Plus the same r as a T1 shortening in per cent, which is the number you actually
%   want: "the T1 you have been quoting was X% short."
%
%   The drift control (A vs D), including the error-bar factor -- two identical runs should
%   differ by 0.674 sigma per point, the median of a folded standard normal. ESR came back at
%   12.3 on this setup and Echo at 2.63, so this is measured, not assumed.
%
%   A SENSITIVITY BOUND when nothing is found: the smallest r that would have shown at 3 sigma,
%   converted into the T1 shortening it corresponds to. "No leak" without that number is not a
%   result.
%
% ===================================================================
% BEFORE YOU RUN IT
% ===================================================================
%
%   >> testT1Gate      % ~2 s, no hardware. Grades the sequence arithmetic and these statistics.
%
% It needs the patched AWGExpT1 -- the one with mwGateMode. Without it this file stops with a
% message saying so rather than silently running four identical arms and reporting no leak,
% which is the one failure mode that would look like good news.
%
% Time: four arms, each a full T1 sweep with constantTime on. At the defaults (14 taus to
% 4900 us, 500 repeats, 3 averages) that is about 6 minutes per arm and 25 in total. The banner
% prints the estimate before anything starts, so you can stop and lower repeats.

    DEFAULTS = struct( ...
        'tauUs', [], 'tauMaxUs', 4900, 'nTau', 14, 'tauMinUs', 20, ...
        'repeats', 500, 'averages', 3, 'nChannels', 2, 'amplitude', -3, ...
        'jsonLocation', '', 'initUs', [], 'lastDelayUs', 1, 'mwHoldUs', [], ...
        'contrastSigma', 3, 'plotIt', true, 'verbose', true, ...
        'selfTest', false, 'leakRatePerUs', 0, 'synthT1Us', 2500, ...
        'synthContrast', 0.03, 'synthNoise', 4e-4, 'synthSeed', 1, ...
        'synthArmDriftFrac', 0, 'synthArmT1DriftFrac', 0);

    if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
    opts = withDefaults(opts, DEFAULTS, 't1GateTest');

    % ---- the tau list ----------------------------------------------------
    % Log-spaced, and starting well above the gate hold. Log because the decay knee is what
    % carries T1, and because the log-RATIO statistics below want tau spread over a wide range:
    % a slope of -r against tau is measured by the lever arm, and a linear grid spends most of
    % its points at the long end where the contrast has already gone.
    tau = opts.tauUs;
    if isempty(tau)
        tau = unique(roundTo(logspace(log10(opts.tauMinUs), log10(opts.tauMaxUs), opts.nTau), 3));
    end
    tau = sort(tau(:)');
    out = struct('ok', false, 'why', '', 'opts', opts, 'tau', tau);

    if ~opts.selfTest && ~hasGateKnob()
        out.why = ['this needs the patched AWGExpT1, the one with an mwGateMode property. ' ...
                   'Without it all four arms would be identical and the test would report ' ...
                   '"no leak" for the wrong reason.'];
        fprintf('\nt1GateTest STOPPED: %s\n\n', out.why);
        return
    end

    banner(tau, opts);

    % ================================================================
    % THE FOUR ARMS
    % ================================================================
    % D last, so that A and D bracket the whole test in time. A drift control that sits next to
    % the arm it controls measures nothing -- the point is to span every minute that B and C
    % were taken in.
    A = runOne(tau, 'tau',   opts, 'A  gate = tau     (as shipped)');
    B = runOne(tau, 'pulse', opts, 'B  gate = pulse   (closed during the wait)');
    C = runOne(tau, 'all',   opts, 'C  gate = all     (open through tau AND lastDelay)');
    D = runOne(tau, 'tau',   opts, 'D  gate = tau     (repeat of A, drift control)');

    for g = {A, B, C, D}
        if isempty(g{1}.contrast)
            out.why = sprintf('arm "%s" produced no usable contrast -- see the run output above.', g{1}.tag);
            fprintf('\nt1GateTest STOPPED: %s\n\n', out.why); return
        end
    end
    out.A = A; out.B = B; out.C = C; out.D = D;

    % ================================================================
    % DRIFT AND THE ERROR-BAR FACTOR, FIRST
    % ================================================================
    % Before any comparison between arms, because the answer decides how to read every number
    % below it. If A and D disagree, nothing that separates B from C means anything.
    out.drift = driftBetween(A, D);
    if opts.verbose; reportDrift(out.drift); end

    % ================================================================
    % T1 PER ARM, MODEL-FREE
    % ================================================================
    out.T1 = struct();
    for g = {A, B, C, D}
        s = logSlope(g{1}.tau, g{1}.contrast, g{1}.cerr, opts.contrastSigma);
        s.T1 = -1/s.slope; s.T1Err = abs(s.slopeErr / s.slope^2);
        out.T1.(g{1}.key) = s;
    end
    if opts.verbose; reportT1(out.T1, out); end

    % ================================================================
    % THE LEAK RATE, TWO WAYS
    % ================================================================
    % SLOPE: log(cA/cB) against tau. The gate-open time differs by (tau - hold) between these
    % two arms, so the slope IS -r and the intercept should be small and positive.
    out.slopeAB = logRatio(A, B, opts.contrastSigma);
    out.rFromSlope    = -out.slopeAB.slope;
    out.rFromSlopeErr =  out.slopeAB.slopeErr;

    % OFFSET: log(cC/cB) against tau. C's dose is CONSTANT, so this ratio must be flat, and its
    % height divided by that dose is r. The slope being zero is not a bonus -- it is the check
    % that the offset is a dose effect rather than a drift.
    out.slopeCB = logRatio(C, B, opts.contrastSigma);
    Gc = C.gateOpenUs - B.gateOpenUs;
    out.doseC = Gc;
    out.rFromOffset    = -out.slopeCB.intercept / max(Gc, eps);
    out.rFromOffsetErr =  out.slopeCB.interceptErr / max(Gc, eps);

    % And the controlled pair, for the record: A vs C differ only in whether lastDelay carries
    % MW, at identical total duration. Its slope must be +r.
    out.slopeAC = logRatio(A, C, opts.contrastSigma);

    % The agreement between the two estimates is computed HERE rather than in verdictGate,
    % because reportLeak prints it and used to run first -- which raised "structure has no
    % member 'agree'" on the very first self test. verdictGate reads this field rather than
    % recomputing it, so the number printed and the number judged are the same one.
    out.agree = NaN;
    if isfinite(out.rFromSlope) && isfinite(out.rFromOffset)
        out.agree = (out.rFromSlope - out.rFromOffset) / ...
                    max(sqrt(out.rFromSlopeErr^2 + out.rFromOffsetErr^2), eps);
    end

    if opts.verbose; reportLeak(out); end

    % ================================================================
    % VERDICT
    % ================================================================
    out = verdictGate(out, opts);
    if opts.verbose; reportVerdict(out); end
    if opts.plotIt;  drawGate(out); end
    out.ok = true;
end

% ======================================================================
% RUNNING ONE ARM
% ======================================================================

function g = runOne(tau, mode, opts, tag)
    if opts.verbose; fprintf('\n--- %s ---\n', tag); end
    key = tag(1);

    if opts.selfTest
        g = synthRun(tau, mode, opts, tag);
        g.key = key;
        return
    end

    t = AWGExpT1;
    if ~isempty(opts.jsonLocation)
        t.loadJsonParams('all', opts.jsonLocation);
    end
    t.nChannels = opts.nChannels;                       % BEFORE amplitude, always:
    t.amplitude = opts.amplitude + zeros(1, opts.nChannels);   % set.amplitude adds amplitudeDiff
    t.timeDelay = 0;                                    % channelPhase is the phase knob
    if ~isempty(opts.initUs); t.laserInitializationDuration = opts.initUs; end
    t.repeats           = opts.repeats;
    t.averages          = opts.averages;
    t.lastDelay         = opts.lastDelayUs;
    t.doubleMeasurement = true;                         % both projections: the |0> arm is the control
    t.isTracking        = 0;                            % a tracking step mid-arm is a systematic
    t.shouldAutosave    = false;
    t.constantTime      = true;                         % REQUIRED: see the header on timing-identical arms
    t.tau               = tau;
    t.mwGateMode        = mode;
    if ~isempty(opts.mwHoldUs); t.mwHoldUs = opts.mwHoldUs; end

    styleWas = '';
    try; styleWas = t.AWG.progressStyle; t.AWG.progressStyle = 'milestone'; catch; end
    tId = tic;                      % ID FORM: a bare tic is reset by any tic inside t.run
    try
        t.run;
    catch err
        fprintf('  *** the run threw: %s\n', err.message);
    end
    secs = toc(tId);
    try; t.AWG.progressStyle = styleWas; catch; end

    g = expGrab(t, tag);
    g.key = key; g.mode = mode; g.seconds = secs;
    g.gateOpenUs = NaN;
    if isfield(g.settings, 'mwGateOpenUs'); g.gateOpenUs = g.settings.mwGateOpenUs; end
    g.holdUs = NaN;
    if isfield(g.settings, 'mwHoldUsed'); g.holdUs = g.settings.mwHoldUsed; end
    g = addContrast(g, tau);
    if opts.verbose
        fprintf('  %.0f s, gate open %.3f us per repeat-half at the longest tau, %d/%d points with contrast\n', ...
                g.seconds, g.gateOpenUs, sum(isfinite(g.contrast)), numel(tau));
    end
end

function g = addContrast(g, tau)
% The |0> minus |1> contrast, from t1Contrast -- one implementation of the orientation rule,
% shared with t1FitRobust and testable on its own. See that file for why the orientation is
% measured rather than read off the field names.
    c = t1Contrast(struct('grab', g, 'tau', tau));
    g.tau = c.tau; g.contrast = c.contrast; g.cerr = c.cerr;
    g.flipped = c.flipped; g.confident = c.confident;
    if ~isempty(c.why)
        fprintf('  *** contrast: %s\n', c.why);
    end
end

% ======================================================================
% STATISTICS
% ======================================================================

function s = logSlope(tau, c, e, kSig)
%LOGSLOPE Weighted straight line through log(contrast) against tau.
%
% The slope is -1/T1 when the contrast is a single exponential. Points whose contrast is not
% significant are DROPPED rather than down-weighted: log of a quantity consistent with zero is
% not a number with an error bar, it is unbounded below, and one such point drags a weighted fit
% wherever its enormous nominal weight points. kSig = 3 by default and the count is reported, so
% a sweep that ran past the decay is visible as a shrinking n rather than as a wild slope.
    s = struct('slope', NaN, 'slopeErr', NaN, 'intercept', NaN, 'interceptErr', NaN, ...
               'n', 0, 'chi2', NaN, 'used', []);
    keep = isfinite(c) & isfinite(e) & (c > kSig*max(e, eps)) & isfinite(tau);
    s.used = keep; s.n = sum(keep);
    if s.n < 4; return; end
    x = tau(keep)'; L = log(c(keep))'; sL = (e(keep)./c(keep))';
    s = fillFit(s, wlsLine(x, L, sL));
end

function s = logRatio(g1, g2, kSig)
%LOGRATIO Weighted straight line through log(c1/c2) against tau.
%
% c0 and T1 cancel in the ratio, so what is left is entirely the dose difference between the two
% arms. That is the point: no model of the relaxation enters, and the two arms need not even have
% the same T1 for the SLOPE to be meaningful -- though if they do not, the drift control should
% already have said so.
    s = struct('slope', NaN, 'slopeErr', NaN, 'intercept', NaN, 'interceptErr', NaN, ...
               'n', 0, 'chi2', NaN, 'used', []);
    n = min(numel(g1.contrast), numel(g2.contrast));
    if n < 4; return; end
    c1 = g1.contrast(1:n); e1 = g1.cerr(1:n);
    c2 = g2.contrast(1:n); e2 = g2.cerr(1:n);
    tau = g1.tau(1:n);
    keep = isfinite(c1) & isfinite(c2) & (c1 > kSig*max(e1, eps)) & (c2 > kSig*max(e2, eps));
    s.used = keep; s.n = sum(keep);
    if s.n < 4; return; end
    x = tau(keep)';
    L = log(c1(keep)./c2(keep))';
    sL = sqrt((e1(keep)./c1(keep)).^2 + (e2(keep)./c2(keep)).^2)';
    s = fillFit(s, wlsLine(x, L, sL));
end

function s = fillFit(s, f)
    s.slope = f.b; s.slopeErr = f.bErr;
    s.intercept = f.a; s.interceptErr = f.aErr;
    s.chi2 = f.chi2;
end

function f = wlsLine(x, y, sy)
%WLSLINE Weighted least squares straight line, with Birge-scaled errors.
%
% The Birge scaling is not decoration. Every error bar on this setup is understated -- the drift
% control measures the factor and it has come back at 12.3 for ESR and 2.63 for Echo -- so an
% error propagated from sterr alone would make a 1-sigma disagreement look like 4. Scaling by
% sqrt(reduced chi2) where that exceeds 1 lets the SCATTER OF THE POINTS set the error bar
% instead, which is the only honest option when the quoted sigma is known to be wrong.
    f = struct('a', NaN, 'aErr', NaN, 'b', NaN, 'bErr', NaN, 'chi2', NaN);
    ok = isfinite(x) & isfinite(y) & isfinite(sy) & sy > 0;
    x = x(ok); y = y(ok); sy = sy(ok);
    n = numel(x);
    if n < 3; return; end
    w = 1./sy.^2;
    X = [ones(n,1), x(:)];
    M = X' * (X .* w(:));
    if rcond(M) < 1e-14; return; end
    th = M \ (X' * (y(:) .* w(:)));
    r = y(:) - X*th;
    chi2 = sum(w(:).*r.^2) / max(n - 2, 1);
    C = inv(M) * max(chi2, 1);
    f.a = th(1); f.b = th(2);
    f.aErr = sqrt(abs(C(1,1))); f.bErr = sqrt(abs(C(2,2)));
    f.chi2 = chi2;
end

function d = driftBetween(a, b)
%DRIFTBETWEEN Two runs of the SAME arm. Everything they disagree about is drift or understated
% error, and the two are told apart by which statistic moves.
    d = struct('sigmaPerPoint', NaN, 'factor', NaN, 'scale', NaN, 'scaleErr', NaN, 'n', 0);
    n = min(numel(a.contrast), numel(b.contrast));
    if n < 4; return; end
    ca = a.contrast(1:n); ea = a.cerr(1:n);
    cb = b.contrast(1:n); eb = b.cerr(1:n);
    s = sqrt(ea.^2 + eb.^2);
    ok = isfinite(ca) & isfinite(cb) & s > 0;
    if sum(ok) < 4; return; end
    d.n = sum(ok);
    d.sigmaPerPoint = median(abs(ca(ok) - cb(ok)) ./ s(ok));
    % 0.674 is the median of a folded standard normal: two runs that agree, with honest error
    % bars, differ by that much per point. Anything above it is either drift or an understated
    % sigma, and the factor is how much.
    d.factor = d.sigmaPerPoint / 0.674;
    % And the scale, which separates the two: a pure error-bar problem leaves the SHAPE alone, so
    % the regression of one on the other still comes back at 1.
    f = wlsLine(cb(ok)', ca(ok)', s(ok)');
    d.scale = f.b; d.scaleErr = f.bErr;
end

% ======================================================================
% VERDICT
% ======================================================================

function out = verdictGate(out, opts)
    r1 = out.rFromSlope;    s1 = out.rFromSlopeErr;
    r2 = out.rFromOffset;   s2 = out.rFromOffsetErr;
    out.leakDetected = false;
    if ~isfield(out, 'agree'); out.agree = NaN; end

    sig1 = r1 / max(s1, eps);
    sig2 = r2 / max(s2, eps);
    out.sigmaSlope = sig1; out.sigmaOffset = sig2;
    % BOTH have to fire, and they have to agree. One arm crossing 3 sigma on its own is what a
    % drift looks like; the slope and the offset carrying the same rate is what a dose effect
    % looks like, and nothing else does.
    out.leakDetected = isfinite(sig1) && isfinite(sig2) && sig1 > 3 && sig2 > 3 && ...
                       abs(out.agree) < 3;

    % The rate as a T1 shortening. 1/T1_measured = 1/T1_true + r, so with the gate open for the
    % whole of tau the fractional shortening at the T1 arm B measures is r*T1/(1 + r*T1).
    T1clean = NaN;
    if isfield(out.T1, 'B'); T1clean = out.T1.B.T1; end
    out.T1clean = T1clean;
    out.shortenFrac    = NaN; out.shortenFracErr = NaN;
    % TWO SENSITIVITY BOUNDS, because the two estimates have genuinely different systematics and
    % quoting only the tighter one would be overselling the test.
    %
    % The OFFSET arm is far more precise -- measured in the self test, 7 to 10 times tighter --
    % because a constant uses every point at full weight while a slope is limited by the lever
    % arm in tau, and a leak kills the long-tau contrast first, shortening exactly that lever.
    % But an offset is also what ANY constant multiplicative difference between two arms looks
    % like, from any cause, and the "C/B slope is zero" check does not exclude that.
    %
    % The SLOPE arm is immune to a constant offset by construction, which makes it the
    % conservative bound and the one the verdict quotes. Both are reported.
    out.sensitivityR       = 3*max(s1, eps);
    out.sensitivityROffset = 3*max(s2, eps);
    out.sensitivityFrac = NaN; out.sensitivityFracOffset = NaN;
    if isfinite(T1clean) && T1clean > 0
        out.shortenFrac     = r1*T1clean / (1 + r1*T1clean);
        out.shortenFracErr  = s1*T1clean / (1 + r1*T1clean)^2;
        rS = out.sensitivityR;
        out.sensitivityFrac = rS*T1clean / (1 + rS*T1clean);
        rO = out.sensitivityROffset;
        out.sensitivityFracOffset = rO*T1clean / (1 + rO*T1clean);
    end
    out.doseAtLongestTau = out.A.gateOpenUs;
end

% ======================================================================
% REPORTING
% ======================================================================

function banner(tau, opts)
    if ~opts.verbose; return; end
    fprintf('\n=== t1GateTest ===\n');
    if opts.selfTest
        fprintf('SELF TEST: no hardware. Injected leak rate r = %.3g per us.\n', opts.leakRatePerUs);
    end
    fprintf('tau: %d points, %.0f to %.0f us. constantTime is forced ON, so all four arms have\n', ...
            numel(tau), min(tau), max(tau));
    fprintf('the same total duration at every tau and differ ONLY in when the switch is open.\n');
    if ~opts.selfTest
        % Per repeat-half: tau (or the constant max) plus the fixed overhead, times two halves.
        overheadUs = 0.1 + opts.lastDelayUs + 0.5 + 9.5 + 0.5 + 1;
        perArm = numel(tau) * 2 * (max(tau) + overheadUs) * opts.repeats * opts.averages * 1e-6;
        fprintf('estimated %.0f s per arm, %.0f min for all four. Lower repeats if that is too long.\n', ...
                perArm, 4*perArm/60);
    end
end

function reportDrift(d)
    fprintf('\nDRIFT CONTROL (A vs D, the two ends of the test)\n');
    if ~isfinite(d.sigmaPerPoint)
        fprintf('  not computable -- one of the arms has too few usable points.\n'); return
    end
    fprintf('  median %.2f sigma per point (two honest runs give 0.674), scale %+.3f +/- %.3f\n', ...
            d.sigmaPerPoint, d.scale, d.scaleErr);
    fprintf('  ERROR BARS: sterr here is UNDERSTATED BY %.1fx.\n', d.factor);
    if d.factor > 3
        fprintf(['  Every error propagated from sterr below is scaled by the scatter of the points\n' ...
                 '  (Birge), which is the only honest option when the quoted sigma is known wrong.\n']);
    end
    if isfinite(d.scale) && abs(d.scale - 1) > 3*d.scaleErr
        fprintf(['  *** and the SHAPE moved too, not just the scatter: the scale is %.1f sigma from 1.\n' ...
                 '      Something drifted across the test. Read the arm comparisons as upper bounds\n' ...
                 '      on how well they can be separated, not as measurements.\n'], ...
                abs(d.scale - 1)/max(d.scaleErr, eps));
    end
end

function reportT1(T, out)
    fprintf('\nT1 PER ARM, from the slope of log(contrast) -- single exponential, no fit\n');
    keys = {'A', 'B', 'C', 'D'};
    labs = {'A gate = tau  ', 'B gate = pulse', 'C gate = all  ', 'D gate = tau  '};
    for i = 1:4
        if ~isfield(T, keys{i}); continue; end
        s = T.(keys{i});
        fprintf('  %s  T1 = %8.1f +/- %6.1f us   (%d points, chi2 %.2f)\n', ...
                labs{i}, s.T1, s.T1Err, s.n, s.chi2);
    end
    if isfield(T, 'A') && isfield(T, 'B') && isfinite(T.A.T1) && isfinite(T.B.T1)
        dd = T.B.T1 - T.A.T1;
        sd = sqrt(T.A.T1Err^2 + T.B.T1Err^2);
        fprintf(['  B minus A = %+.1f +/- %.1f us (%.1f sigma). This is the headline in the units\n' ...
                 '  you care about: how much T1 the open switch was costing.\n'], dd, sd, abs(dd)/max(sd,eps));
    end
    if isfield(out, 'A') && isfinite(out.A.gateOpenUs)
        fprintf('  doses: A %.1f us, B %.3f us, C %.1f us of open gate per repeat-half\n', ...
                out.A.gateOpenUs, out.B.gateOpenUs, out.C.gateOpenUs);
    end
end

function reportLeak(out)
    fprintf('\nTHE LEAK RATE, TWO INDEPENDENT WAYS\n');
    fprintf('  from the SLOPE  (A vs B): r = %+.4g +/- %.3g /us   (%.1f sigma, %d points, chi2 %.2f)\n', ...
            out.rFromSlope, out.rFromSlopeErr, out.rFromSlope/max(out.rFromSlopeErr,eps), ...
            out.slopeAB.n, out.slopeAB.chi2);
    fprintf('  from the OFFSET (C vs B): r = %+.4g +/- %.3g /us   (%.1f sigma, dose %.0f us)\n', ...
            out.rFromOffset, out.rFromOffsetErr, out.rFromOffset/max(out.rFromOffsetErr,eps), out.doseC);
    fprintf('    C vs B slope %+.3g +/- %.3g -- MUST be zero: C''s dose does not depend on tau,\n', ...
            out.slopeCB.slope, out.slopeCB.slopeErr);
    if isfinite(out.slopeCB.slope) && abs(out.slopeCB.slope) > 3*out.slopeCB.slopeErr
        fprintf(['    *** it is %.1f sigma from zero. C''s offset is then NOT purely a dose effect,\n' ...
                 '        so the second estimate of r is contaminated and the two-way agreement\n' ...
                 '        below is not the check it is meant to be.\n'], ...
                abs(out.slopeCB.slope)/max(out.slopeCB.slopeErr,eps));
    else
        fprintf('        it is, within errors. Good: the offset is a dose effect.\n');
    end
    % A and C have IDENTICAL total timing and differ only in whether lastDelay carries MW, so
    % this is the most controlled pair in the test. Its slope carries the same -r as A vs B
    % (C's dose does not depend on tau, so it drops out of the slope), which makes it a check on
    % the first estimate rather than a third one.
    fprintf('  the controlled pair (A vs C, identical timing): slope %+.3g +/- %.3g against the\n', ...
            out.slopeAC.slope, out.slopeAC.slopeErr);
    fprintf('    A vs B slope of %+.3g +/- %.3g -- the same quantity, measured against a different\n', ...
            out.slopeAB.slope, out.slopeAB.slopeErr);
    dAC = (out.slopeAC.slope - out.slopeAB.slope) / ...
          max(sqrt(out.slopeAC.slopeErr^2 + out.slopeAB.slopeErr^2), eps);
    fprintf('    reference. They differ by %.1f sigma.\n', abs(dAC));
    if isfinite(out.agree)
        fprintf('  the two estimates of r differ by %.1f sigma\n', abs(out.agree));
    end
end

function reportVerdict(out)
    fprintf('\nVERDICT\n');
    if out.leakDetected
        fprintf('  THE GATE LEAKS. r = %.4g +/- %.3g per us of open switch, confirmed by both the\n', ...
                out.rFromSlope, out.rFromSlopeErr);
        fprintf('  slope (%.1f sigma) and the offset (%.1f sigma), agreeing to %.1f sigma.\n', ...
                out.sigmaSlope, out.sigmaOffset, abs(out.agree));
        if isfinite(out.shortenFrac)
            fprintf(['  With the gate open for the whole of tau, that makes the measured T1 SHORT by\n' ...
                     '  %.1f%% +/- %.1f%%. Every T1 taken on this setup before now is affected.\n'], ...
                    100*out.shortenFrac, 100*out.shortenFracErr);
        end
        fprintf('  WHAT TO DO: set mwGateMode = ''pulse'' in the run script and leave it there.\n');
    else
        fprintf('  NO LEAK FOUND.\n');
        if isfinite(out.sensitivityFrac)
            fprintf(['  SENSITIVITY, conservatively (from the SLOPE, which no constant offset can\n' ...
                     '  fake): a rate above %.3g /us would have shown at 3 sigma, which at the %.0f us\n' ...
                     '  T1 measured with the gate closed is a %.1f%% shortening.\n'], ...
                    out.sensitivityR, out.T1clean, 100*out.sensitivityFrac);
            fprintf(['  and from the OFFSET arm, %.1fx tighter but assuming no constant arm-to-arm\n' ...
                     '  systematic: above %.3g /us, i.e. a %.1f%% shortening.\n'], ...
                    out.sensitivityR/max(out.sensitivityROffset, eps), ...
                    out.sensitivityROffset, 100*out.sensitivityFracOffset);
            fprintf(['  Smaller than that is NOT excluded -- "no leak" means "none above this", and\n' ...
                     '  the number is the result. More repeats and a longer tau reach both help.\n']);
        end
        fprintf(['  mwGateMode = ''pulse'' is still the better default: it costs nothing, it is what\n' ...
                 '  the non-AWG ExpT1 does, and an open switch passing amplifier noise for %.0f us is\n' ...
                 '  a risk with no upside even when it is currently below this bound.\n'], ...
                out.doseAtLongestTau);
    end
    if isfinite(out.drift.factor) && out.drift.factor > 3
        fprintf(['  Read all of the above with the drift control in mind: sterr on this setup is\n' ...
                 '  understated %.1fx, and only the Birge-scaled errors quoted here account for it.\n'], ...
                out.drift.factor);
    end
    fprintf('\n');
end

% ======================================================================
% SELF TEST
% ======================================================================

function g = synthRun(tau, mode, opts, tag)
%SYNTHRUN A T1 arm built from the physics, with a KNOWN leak rate.
%
% The dose per mode is computed here the same way AWGExpT1 computes it, and deliberately not by
% calling the class: a self test that gets its dose from the code under test cannot detect a dose
% that is wrong.
    holdUs = 2*0.15 + 0.035 + 0.1;                  % 0.435 us, as prepare derives it
    lastDelay = opts.lastDelayUs;
    maxLastDelay = lastDelay + max(tau);
    switch mode
        case 'pulse'; G = 0.1 + holdUs + zeros(size(tau));
        case 'all';   G = 0.1 + max(tau) + lastDelay + zeros(size(tau));   % constant, by design
        otherwise;    G = 0.1 + tau;
    end
    % ---- the two confounds this test has to survive ------------------------------
    % The arms run in sequence -- A, B, C, D -- so anything that drifts over the half hour the
    % test takes drifts BETWEEN them, and the two ways it can do so are exactly the two ways a
    % leak can look real. Both are injectable so that testT1Gate can check the design claim
    % rather than the header merely asserting it.
    %
    %   synthArmDriftFrac    the contrast AMPLITUDE drifts by this fraction per arm. A constant
    %                        scale factor is a pure OFFSET in the log-ratio and no slope, so it
    %                        fakes the C-vs-B estimate and cannot touch the A-vs-B one.
    %   synthArmT1DriftFrac  T1 itself drifts by this fraction per arm. That is a pure SLOPE and
    %                        no offset, so it fakes A-vs-B and cannot touch C-vs-B.
    %
    % Requiring the two estimates to AGREE is what rejects both. Neither confound can produce
    % the same r in a slope and in an offset, because each produces only one of them.
    armIdx = find(strcmp(tag(1), {'A', 'B', 'C', 'D'})) - 1;
    if isempty(armIdx); armIdx = 0; end
    amp = opts.synthContrast * (1 + opts.synthArmDriftFrac   * armIdx);
    T1a = opts.synthT1Us     * (1 + opts.synthArmT1DriftFrac * armIdx);
    c = amp * exp(-tau/T1a) .* exp(-opts.leakRatePerUs*G);
    e = opts.synthNoise * ones(size(tau));
    seed = 100*opts.synthSeed + double(tag(1));
    c = c + e .* synthNoise(numel(tau), seed);

    g = struct('tag', tag, 'mode', mode, 'seconds', 0, 'flipped', false, ...
               'tau', tau, 'contrast', c, 'cerr', e, ...
               'gateOpenUs', G(end), 'holdUs', holdUs, ...
               'value', [], 'sterr2', [], 'value2', [], 'settings', struct());
    if strcmp(mode, 'all'); g.gateOpenUs = G(1); end
    if opts.verbose
        fprintf('  synthetic: dose %.3f us at the longest tau, %d points\n', g.gateOpenUs, numel(tau));
    end
end

function r = synthNoise(n, seed)
% Written-out LCG plus Box-Muller: rng() does not exist in Octave, where this file is
% parse-checked, so a failure is reproducible on whichever machine saw it.
    m = 2^31 - 1;
    x = mod(1103515245*(seed + 12345) + 12345, m);
    u = zeros(1, 2*n);
    for i = 1:2*n
        x = mod(1103515245*x + 12345, m);
        u(i) = (x + 0.5)/m;
    end
    r = sqrt(-2*log(u(1:n))) .* cos(2*pi*u(n+1:2*n));
end

% ======================================================================
% PLOT
% ======================================================================

function drawGate(out)
% Two panels, and the second is the measurement. The top one shows the four contrast decays,
% which is what anyone looks at first and which -- if there is a leak -- shows almost nothing,
% because a leak is a change of slope on a log scale hiding inside a decay. The bottom shows the
% LOG RATIOS, where a leak is a straight line with a slope and no leak is a flat line at zero.
% That is the whole argument of this test in one axes.
    S = fitStyle();
    fig = S.newFigure('T1 gate test', true);

    axT = S.axes(subplot(2, 1, 1));
    H = []; L = {};
    keys = {'B', 'A', 'C'};
    labs = {'B  gate = pulse', 'A  gate = tau', 'C  gate = all'};
    for i = 1:3
        g = out.(keys{i});
        H(end+1) = S.data(axT, g.tau, g.contrast, g.cerr, S.slots{i}, '', true, []);  %#ok<AGROW>
        L{end+1} = labs{i};                                                            %#ok<AGROW>
    end
    set(axT, 'YScale', 'log');      % a leak is a slope change; on a linear axis it is invisible
    xlabel(axT, '\tau  (\mus)', 'Color', S.col.ink2);
    ylabel(axT, 'contrast   |0\rangle - |1\rangle', 'Color', S.col.ink2);
    S.tight(axT, out.A.tau);
    S.legend(axT, H, L);
    S.title(axT, gateTitle(out));

    axR = S.axes(subplot(2, 1, 2));
    S.zeroLine(axR, [min(out.A.tau) max(out.A.tau)]);
    H = []; L = {};
    pairs = {{out.A, out.B, 'log(A/B): slope = -r', 1}, ...
             {out.C, out.B, 'log(C/B): offset only, slope must be 0', 2}};
    for i = 1:numel(pairs)
        p = pairs{i};
        n = min(numel(p{1}.contrast), numel(p{2}.contrast));
        c1 = p{1}.contrast(1:n); e1 = p{1}.cerr(1:n);
        c2 = p{2}.contrast(1:n); e2 = p{2}.cerr(1:n);
        good = c1 > 0 & c2 > 0;
        Lr = nan(1, n); Er = nan(1, n);
        Lr(good) = log(c1(good)./c2(good));
        Er(good) = sqrt((e1(good)./c1(good)).^2 + (e2(good)./c2(good)).^2);
        H(end+1) = S.data(axR, p{1}.tau(1:n), Lr, Er, S.slots{p{4}}, '', true, []);   %#ok<AGROW>
        L{end+1} = p{3};                                                             %#ok<AGROW>
    end
    td = linspace(min(out.A.tau), max(out.A.tau), 200);
    if isfinite(out.slopeAB.slope)
        h = S.model(axR, td, out.slopeAB.intercept + out.slopeAB.slope*td);
        H(end+1) = h; L{end+1} = sprintf('fit, r = %.3g \\pm %.3g /\\mus', ...
                                         out.rFromSlope, out.rFromSlopeErr);
    end
    if isfinite(out.slopeCB.intercept)
        S.envelope(axR, td, out.slopeCB.intercept + zeros(size(td)));
    end
    xlabel(axR, '\tau  (\mus)', 'Color', S.col.ink2);
    ylabel(axR, 'log contrast ratio', 'Color', S.col.ink2);
    S.tight(axR, out.A.tau);
    S.legend(axR, H, L);
    S.title(axR, sprintf(['A/B slope %+.3g \\pm %.3g     |     C/B offset %+.3g \\pm %.3g over a ' ...
                          '%.0f \\mus dose     |     the two agree to %.1f\\sigma'], ...
                         out.slopeAB.slope, out.slopeAB.slopeErr, ...
                         out.slopeCB.intercept, out.slopeCB.interceptErr, out.doseC, ...
                         abs(out.agree)));
end

function T = gateTitle(out)
    a = 'T1 gate test';
    if out.leakDetected
        a = sprintf('THE GATE LEAKS: r = %.3g \\pm %.3g /\\mus, T_1 short by %.1f%%', ...
                    out.rFromSlope, out.rFromSlopeErr, 100*out.shortenFrac);
    elseif isfinite(out.sensitivityFrac)
        a = sprintf('no leak above %.3g /\\mus, i.e. above a %.1f%% T_1 shortening (3\\sigma)', ...
                    out.sensitivityR, 100*out.sensitivityFrac);
    end
    bits = {};
    for k = {'A','B','C','D'}
        if isfield(out.T1, k{1}) && isfinite(out.T1.(k{1}).T1)
            bits{end+1} = sprintf('T_1(%s) = %.0f \\mus', k{1}, out.T1.(k{1}).T1);   %#ok<AGROW>
        end
    end
    b = strjoin(bits, '     |     ');
    c = sprintf('drift %.2f\\sigma per point, sterr understated %.1f\\times', ...
                out.drift.sigmaPerPoint, out.drift.factor);
    T = {a, b, c};
end

% ======================================================================
% HELPERS
% ======================================================================

function tf = hasGateKnob()
% Does the installed AWGExpT1 have mwGateMode? Asked of the CLASS rather than of an instance, so
% nothing is constructed -- constructing an AWGExpT1 deletes any live one, and a capability check
% has no business destroying the object the user was working with.
    tf = false;
    try
        m = meta.class.fromName('AWGExpT1');
        tf = any(strcmp({m.PropertyList.Name}, 'mwGateMode'));
    catch
        try
            tf = ismember('mwGateMode', properties('AWGExpT1'));
        catch
        end
    end
end

function y = roundTo(x, n)
% round(x, n) is a MATLAB-only two-argument form; this file is parse-checked in Octave.
    y = round(x .* 10.^n) ./ 10.^n;
end

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d); given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('%s: unknown option "%s". Valid: %s.', who, bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}); o.(known{i}) = d.(known{i}); end
    end
end
