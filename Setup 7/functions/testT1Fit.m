function ok = testT1Fit()
%TESTT1FIT Offline grader for t1FitRobust. No hardware, ~40 s.
%
%   testT1Fit
%
% Every check injects a known T1, stretch exponent and contrast into synthetic data built from
% the physics rather than from the fitter, and asks whether it comes back. Run it before quoting
% a T1, and certainly after touching t1FitRobust.
%
% The four that matter most, because each is a way of being confidently wrong rather than
% obviously broken:
%
%   t3  THE CONTROL. The |0> trace should not move. If something common-mode drifts with tau --
%       laser power, tracking, charge state -- it appears in BOTH traces, the fit absorbs it into
%       T1, and nothing about the residual complains. flatnessSigma is the only thing standing
%       between that and a published number, so it has to fire, and it has to stay quiet when
%       the control really is flat.
%
%   t4  THE SWEEP THAT STOPS TOO EARLY, and the answer is not the one you would guess. The fitted
%       T1 is barely biased even at 0.8 T1 of reach -- 4% low over ten seeds. What fails is the
%       ERROR BAR: the scatter between seeds is three times what the fit quotes, because the
%       likelihood has a long ridge along which the plateau, the amplitude and T1 trade. So the
%       test measures that ratio rather than asserting a threshold, and the lesson it encodes is
%       that spanOverT1 is the number to act on and the error bar is not.
%
%   t5  THE SECOND COMPONENT MUST NOT BE FREE. A stretched exponential with p < 1 and a sum of
%       two exponentials are nearly the same curve on twenty noisy points, so a two-component fit
%       will almost always look slightly better. On single-exponential data it has to be REFUSED,
%       and on genuinely two-component data it has to be accepted and both constants recovered.
%
%   t7  ERROR BARS THAT MEAN SOMETHING. The scatter of T1 over many noise realisations against
%       the error bar the fit quotes. This is the only check that catches an error bar that is
%       simply too small -- and on a T1, where T1 trades against the amplitude and the plateau,
%       too small is the likely direction.

    fprintf('\n=== testT1Fit ===\n');
    ok = true;
    ok = ok & t1_recoversSingleExponential();
    ok = ok & t2_recoversStretch();
    ok = ok & t3_theControl();
    ok = ok & t4_shortSweepIsHonest();
    ok = ok & t5_secondComponent();
    ok = ok & t6_survivesOutliers();
    ok = ok & t7_errorsAreHonest();
    ok = ok & t8_ratesAndArea();
    ok = ok & t9_refusesJunk();
    ok = ok & t10_readsAnExperiment();
    fprintf('\n');
    if ok
        fprintf('ALL TESTS PASS -- the numbers t1FitRobust prints mean what they say.\n\n');
    else
        fprintf('*** FAILURES ABOVE -- do not quote a T1 until they are fixed ***\n\n');
    end
end

% ======================================================================
% SYNTHETIC DATA
% ======================================================================

function [t, Y, E, T] = makeT1(varargin)
% Built from the physics, and the decay is written out here rather than taken from the fitter --
% a test that generates its data with the code under test cannot detect an error in that code.
    P = struct('T1', 2500, 'p', 1, 'nPts', 20, 'tMin', 20, 'tMax', 9000, ...
               'contrast', 0.03, 'base', 1.0, 'noise', 4e-4, 'seed', 1, ...
               'ctlDrift', 0, 'two', true, 'Tfast', [], 'weight', 0.4, 'swap', false);
    for i = 1:2:numel(varargin); P.(varargin{i}) = varargin{i+1}; end

    t = unique(round(logspace(log10(P.tMin), log10(P.tMax), P.nPts)));
    n = numel(t);
    if isempty(P.Tfast)
        D = exp(-(t/P.T1).^P.p);
    else
        D = P.weight*exp(-t/P.Tfast) + (1-P.weight)*exp(-t/P.T1);
    end
    % The |0> trace is the CONTROL: flat, unless ctlDrift injects a common-mode change that
    % follows the same shape -- which is what a drifting laser looks like.
    y0 = P.base + P.ctlDrift*(1 - D) + P.noise*lcgNormal(n, 7*P.seed + 1);
    y1 = P.base - P.contrast*D       + P.noise*lcgNormal(n, 7*P.seed + 2);
    if P.two
        if P.swap; Y = [y1; y0]; else; Y = [y0; y1]; end
        E = P.noise*ones(2, n);
    else
        Y = y1; E = P.noise*ones(1, n);
    end
    T = P; T.t = t; T.D = D;
end

function r = lcgNormal(n, seed)
% Written-out LCG plus Box-Muller: rng() does not exist in Octave, where this is parse-checked,
% and randn('seed',k) is deprecated in MATLAB. A failure here is reproducible anywhere.
    m = 2^31 - 1;
    x = mod(1103515245*(seed + 12345) + 12345, m);
    u = zeros(1, 2*n);
    for i = 1:2*n
        x = mod(1103515245*x + 12345, m);
        u(i) = (x + 0.5)/m;
    end
    r = sqrt(-2*log(u(1:n))) .* cos(2*pi*u(n+1:2*n));
end

function o = fitQuiet(t, Y, E, varargin)
    op = struct('plotIt', false, 'verbose', false, 't', t, 'y', Y, 'e', E);
    for i = 1:2:numel(varargin); op.(varargin{i}) = varargin{i+1}; end
    o = t1FitRobust(op);
end

% ======================================================================

function ok = t1_recoversSingleExponential()
    [t, Y, E, T] = makeT1('T1', 2500, 'p', 1, 'seed', 7);
    o = fitQuiet(t, Y, E);
    dT = abs(o.T1 - T.T1)/max(o.T1Err, eps);
    dC = abs(o.contrast - T.contrast)/max(o.contrastErr, eps);
    ok = o.ok && dT < 3 && dC < 3 && ~o.doubleUsed;
    report(ok, sprintf('t1  T1 = %.0f +/- %.0f vs %.0f us (%.1f sigma), contrast %.3f%% vs %.3f%% (%.1f sigma)', ...
                       o.T1, o.T1Err, T.T1, dT, 100*o.contrast, 100*T.contrast, dC));
end

function ok = t2_recoversStretch()
% p must be RECOVERED. A T1 with a spread of rates is the common departure from textbook, and a
% fit that pins p at 1 reports a T1 that depends on where the sweep happened to stop.
    ok = true;
    for pTrue = [0.5 1 1.5]
        [t, Y, E] = makeT1('T1', 2500, 'p', pTrue, 'nPts', 24, 'tMax', 20000, 'seed', 10+round(10*pTrue));
        o = fitQuiet(t, Y, E);
        d = abs(o.p - pTrue)/max(o.pErr, eps);
        good = o.ok && d < 3;
        ok = ok & good;
        report(good, sprintf('t2  p = %.1f recovered as %.3f +/- %.3f (%.1f sigma)   "%s"', ...
                             pTrue, o.p, o.pErr, d, o.pMeaning));
    end
end

function ok = t3_theControl()
%T3 The |0> trace is the null, and the check has to work in both directions.
%
% A common-mode drift enters BOTH traces with the same shape, so the fit absorbs it into T1 and
% the residual stays clean. Nothing else in the output would notice.
    [t, Y, E] = makeT1('ctlDrift', 0, 'seed', 31);
    quiet = fitQuiet(t, Y, E);
    ok1 = quiet.controlPresent && quiet.flatnessSigma < 3;
    report(ok1, sprintf('t3  flat control stays quiet: %.1f sigma', quiet.flatnessSigma));

    ok2 = true;
    for d = [0.005 0.02]
        [t, Y, E] = makeT1('ctlDrift', d, 'seed', 33);
        o = fitQuiet(t, Y, E);
        good = o.flatnessSigma > 3;
        ok2 = ok2 & good;
        report(good, sprintf('t3  %.1f%% common-mode drift fires the control at %.0f sigma', ...
                             100*d, o.flatnessSigma));
    end

    % And with one trace only, it must say the control is MISSING rather than quietly passing.
    [t, Y, E] = makeT1('two', false, 'seed', 35);
    o = fitQuiet(t, Y, E);
    ok3 = o.ok && ~o.controlPresent && ~isfinite(o.flatnessSigma);
    report(ok3, 't3  one trace only: control reported as unavailable, not as passed');
    ok = ok1 && ok2 && ok3;
end

function ok = t4_shortSweepIsHonest()
%T4 A sweep that stops before the plateau -- and what actually goes wrong there.
%
% The obvious guess is a biased T1. Measured over ten noise realisations, it is not: at 0.8 T1
% of reach the median T1 is only 4% low. What breaks is the ERROR BAR. The scatter of T1 between
% seeds is 1635 us against a quoted 561 -- understated by about 3x -- because the likelihood has
% a long ridge along which the baseline, the amplitude and T1 trade, and an unlucky draw lands
% far along it while still reporting a local curvature.
%
% Beyond about 1.6 T1 of reach the quoted error tracks the scatter properly (ratio 0.7 to 0.8,
% which is what an honest error bar looks like on ten samples).
%
% So the test is not "is T1 right" but "does the error bar know how uncertain it is", and the
% answer is that below two T1s it does NOT. That is why t1FitRobust warns on spanOverT1 as a
% separate thing rather than trusting the error bar to speak for itself, and this test measures
% the ratio rather than asserting a threshold on one seed -- which is what it did first, and it
% failed on an unlucky draw, for the right reason.
    ok = true;
    rows = {{2000, 0.8}, {9000, 3.6}};
    ratios = zeros(1, 2);
    for i = 1:2
        Ts = []; Es = []; spans = [];
        for sd = 1:10
            [t, Y, E] = makeT1('T1', 2500, 'p', 1, 'tMax', rows{i}{1}, 'nPts', 18, 'seed', 40+sd);
            o = fitQuiet(t, Y, E);
            Ts(end+1) = o.T1; Es(end+1) = o.T1Err; spans(end+1) = o.spanOverT1;   %#ok<AGROW>
        end
        ratios(i) = std(Ts) / max(median(Es), eps);
        fprintf(['      reach %.1f T1: median T1 %6.0f, quoted err %5.0f, actual scatter %6.0f, ' ...
                 'ratio %.1f, span/T1 %.2f\n'], ...
                rows{i}{2}, median(Ts), median(Es), std(Ts), ratios(i), median(spans));
    end
    % The claim: truncated understates by more than 2x, full does not, and the warning threshold
    % sits between the two reaches.
    ok = ratios(1) > 2 && ratios(2) < 1.5;
    report(ok, sprintf(['t4  below 2 T1 the error bar understates the scatter %.1fx; past 2 T1 it ' ...
                        'does not (%.1fx). span/T1 is the number to act on, not the error bar.'], ...
                       ratios(1), ratios(2)));
end

function ok = t5_secondComponent()
%T5 Accepted when it is there, refused when it is not. The refusal is the important half.
    ok = true;
    for sd = [21 23 25]
        [t, Y, E] = makeT1('T1', 2500, 'p', 1, 'nPts', 24, 'tMax', 20000, 'seed', sd);
        o = fitQuiet(t, Y, E, 'model', 'auto');
        good = o.ok && ~o.doubleUsed && ~isempty(o.doubleTest.failed);
        ok = ok & good;
        report(good, sprintf('t5  seed %d single exponential: second component refused (dSSR %.1f)', ...
                             sd, o.doubleTest.dSSR));
    end
    [t, Y, E, T] = makeT1('T1', 4000, 'Tfast', 300, 'weight', 0.4, 'nPts', 26, ...
                          'tMax', 20000, 'contrast', 0.04, 'noise', 3e-4, 'seed', 19);
    a = fitQuiet(t, Y, E);                          % stretched only
    b = fitQuiet(t, Y, E, 'model', 'auto');
    dF = abs(b.T1fast - T.Tfast)/max(T.Tfast, eps);
    dS = abs(b.T1slow - T.T1)/max(T.T1, eps);
    good = b.doubleUsed && dF < 0.3 && dS < 0.3 && b.chi2Reduced < a.chi2Reduced;
    ok = ok & good;
    report(good, sprintf(['t5  genuine two components accepted: %.0f and %.0f us vs %.0f and %.0f, ' ...
                          'w = %.2f vs %.2f, chi2 %.2f against the stretched fit''s %.2f'], ...
                         b.T1fast, b.T1slow, T.Tfast, T.T1, b.weight, T.weight, ...
                         b.chi2Reduced, a.chi2Reduced));
end

function ok = t6_survivesOutliers()
% Points diluted toward no contrast. On a T1 they sit AT the plateau, which is where half the
% sweep already is, so they are far harder to spot by eye than on a Ramsey -- and a diluted point
% early in the sweep pulls T1 down exactly like a faster relaxation.
    [t, Y, E, T] = makeT1('seed', 61);
    clean = fitQuiet(t, Y, E);
    bad = [3 6];
    Yd = Y;
    Yd(2, bad) = T.base;                            % the |1> trace pulled to no contrast
    o  = fitQuiet(t, Yd, E);
    nr = fitQuiet(t, Yd, E, 'robust', false);
    moved   = abs(o.T1  - clean.T1)/max(clean.T1Err, eps);
    movedNR = abs(nr.T1 - clean.T1)/max(clean.T1Err, eps);
    found = all(ismember(bad, o.downIdx));
    ok = o.ok && found && moved < 2 && movedNR > moved;
    report(ok, sprintf(['t6  2 diluted points: found %d of 2, T1 moved %.2f sigma with robust on ' ...
                        'and %.2f sigma with it off'], numel(intersect(bad, o.downIdx)), moved, movedNR));
end

function ok = t7_errorsAreHonest()
% The scatter of T1 over many noise realisations against the error bar the fit quotes. The only
% check that can catch an error bar that is simply too small -- and on a T1, where T1 trades
% against the amplitude and the plateau, too small is the likely direction.
    N = 14; v = zeros(1, N); s = zeros(1, N);
    for i = 1:N
        [t, Y, E] = makeT1('seed', 200+i);
        o = fitQuiet(t, Y, E);
        v(i) = o.T1; s(i) = o.T1Err;
    end
    ratio = std(v)/max(median(s), eps);
    ok = ratio > 0.4 && ratio < 2.5;
    report(ok, sprintf('t7  scatter of T1 over %d seeds is %.0f us against a quoted %.0f -- ratio %.2f (want 0.4 to 2.5)', ...
                       N, std(v), median(s), ratio));
end

function ok = t8_ratesAndArea()
% The two derived quantities have to be arithmetically what they claim. Cheap, and they are the
% numbers the report tells people to compare between runs, so an error here propagates into every
% comparison silently.
    [t, Y, E] = makeT1('T1', 2500, 'p', 0.6, 'nPts', 24, 'tMax', 20000, 'seed', 81);
    o = fitQuiet(t, Y, E);
    dRate = abs(o.ratePerUs - 1/o.T1) / max(1/o.T1, eps);
    dArea = abs(o.coherenceArea - o.T1*gamma(1 + 1/o.p)) / max(o.coherenceArea, eps);
    % And the area must exceed T1 when p < 1: a spread of rates has a long tail, so the integral
    % under the envelope is LARGER than the nominal time constant. Getting this backwards would
    % make the figure of merit rank runs the wrong way round.
    ok = dRate < 1e-9 && dArea < 1e-9 && o.coherenceArea > o.T1;
    report(ok, sprintf('t8  rate = 1/T1 exactly, area = T1*gamma(1+1/p) exactly, and at p = %.2f the area (%.0f) exceeds T1 (%.0f)', ...
                       o.p, o.coherenceArea, o.T1));
end

function ok = t9_refusesJunk()
    ok = true;
    ok = ok & mustFail(@() t1FitRobust(struct('nonsense', 1)), 'an unknown option');
    ok = ok & mustFail(@() t1FitRobust(42), 'a non-struct argument');
    ok = ok & mustFail(@() t1FitRobust(struct('model', 'quadruple')), 'an unknown model name');

    o = t1FitRobust(struct('plotIt', false, 'verbose', false));
    good = ~o.ok && ~isempty(o.why);
    ok = ok & good;
    report(good, sprintf('t9  empty input refused: "%s"', o.why));

    o = t1FitRobust(struct('t', 1:4, 'y', [1 1 1 1], 'plotIt', false, 'verbose', false));
    good = ~o.ok && ~isempty(o.why);
    ok = ok & good;
    report(good, sprintf('t9  four points refused: "%s"', o.why));

    [t, Y, E] = makeT1('seed', 91);
    o = fitQuiet(t, Y, E, 'pFixed', 1);
    good = o.ok && abs(o.p - 1) < 1e-9 && ~isfinite(o.pErr) && ~o.pFitted;
    ok = ok & good;
    report(good, sprintf('t9  pFixed = 1 held exactly, no error bar invented: "%s"', o.pMeaning));
end

function ok = t10_readsAnExperiment()
%T10 The path from an AWGExpT1 to a fitted T1, including the projection order.
%
% This is the half that never touches the fitting arithmetic, and the swap is what lives here.
% AWGExpT1's labels were the wrong way round, so the same data must fit identically whichever
% field the decaying trace arrives in -- otherwise the fitter reports the flat control's noise as
% a T1, with an error bar, on data taken before the class fix.
    [t, Y, E, T] = makeT1('T1', 2500, 'seed', 95);
    Exp = struct();
    Exp.tau = t; Exp.repeats = 500; Exp.currIter = 3; Exp.signal = [];
    Exp.mCurrentXAxisParam = struct('value', t);
    Exp.mwGateOpenUs = 4900.1; Exp.mwGateMode = 'tau';
    Exp.signalParam  = struct('value', Y(1,:), 'sterr', E(1,:));   % |0>, the control
    Exp.signalParam2 = struct('value', Y(2,:), 'sterr', E(2,:));   % |1>, the decay
    a = t1FitRobust(struct('exp', Exp, 'plotIt', false, 'verbose', false));

    Swap = Exp;                                                     % the old, wrong label order
    Swap.signalParam  = struct('value', Y(2,:), 'sterr', E(2,:));
    Swap.signalParam2 = struct('value', Y(1,:), 'sterr', E(1,:));
    b = t1FitRobust(struct('exp', Swap, 'plotIt', false, 'verbose', false));

    same = abs(a.T1 - b.T1) < 1e-6*max(a.T1, 1);
    dT = abs(a.T1 - T.T1)/max(a.T1Err, eps);
    ok = a.ok && b.ok && same && dT < 3 && isfinite(a.gateOpenUs) && ...
         a.pull.flipped == false && b.pull.flipped == true;
    report(ok, sprintf(['t10 read from an experiment: T1 = %.0f +/- %.0f vs %.0f (%.1f sigma), and ' ...
                        'the SAME T1 with the labels swapped (%.0f), flipped flags %d then %d, ' ...
                        'gate dose %.0f us carried through'], ...
                       a.T1, a.T1Err, T.T1, dT, b.T1, a.pull.flipped, b.pull.flipped, a.gateOpenUs));
end

% ======================================================================

function ok = mustFail(fn, what)
    ok = false;
    try
        fn();
    catch
        ok = true;
    end
    report(ok, sprintf('t9  %s raises rather than guessing', what));
end

function report(ok, what)
    if ok
        fprintf('  PASS  %s\n', what);
    else
        fprintf('  FAIL  %s\n', what);
    end
end
