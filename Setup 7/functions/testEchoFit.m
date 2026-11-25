function ok = testEchoFit()
%TESTECHOFIT Offline grader for echoFitRobust. No hardware, ~20 s.
%
%   testEchoFit
%
% Every check injects a KNOWN T2, stretch exponent, contrast and -- where relevant -- a known
% 13-C revival period into synthetic data built from the physics rather than from the fitter,
% and asks whether it comes back. Run it before quoting a T2 from the setup, and certainly
% after touching echoFitRobust.
%
% The four that matter most, because each corresponds to a way of being CONFIDENTLY WRONG
% rather than obviously broken:
%
%   t3  REVIVAL DATA WITH THE COMB OFF. The echo's version of the single-cosine-on-a-triplet
%       mistake. A monotone stretched exponential fitted through a collapse-and-revival trace
%       returns a short T2 with a small error bar, because the first collapse looks exactly
%       like the end of the decay. The test asserts the size of the lie -- the monotone T2
%       must come back at less than a fifth of the truth -- and then that turning the comb on
%       recovers the real value. If the second half ever starts passing while the first half
%       stops, the model has quietly become a decay fitter again.
%
%   t4  THE COMB MUST NOT BE FREE. On monotone data the comb has to be REFUSED. A model that
%       is allowed to add two parameters whenever it likes will find revivals in noise, and a
%       fitted 13-C field is exactly the kind of number that gets believed because it has
%       units.
%
%   t6  OUTLIERS. The classic path dilutes isolated points toward no contrast, and on an echo
%       a point pulled toward the plateau looks like an honest early loss of coherence rather
%       than like a broken point. Two of them must not move T2.
%
%   t8  THE DRIFT DEGENERACY, measured rather than asserted. This is where the numbers quoted
%       in echoFitRobust's header come from, and the test prints them so a change in the
%       fitter shows up as a change in the table rather than as a stale comment.

    fprintf('\n=== testEchoFit ===\n');
    ok = true;
    ok = ok & t1_recoversMonotone();
    ok = ok & t2_recoversStretch();
    ok = ok & t3_revivalsOrLies();
    ok = ok & t4_combIsNotFree();
    ok = ok & t5_jointBeatsSingle();
    ok = ok & t6_survivesOutliers();
    ok = ok & t7_errorsAreHonest();
    ok = ok & t8_driftDegeneracy();
    ok = ok & t9_flatDataIsCalledFlat();
    ok = ok & t10_refusesJunk();
    ok = ok & t11_readsAnExperiment();
    ok = ok & t12_singleRevival();
    fprintf('\n');
    if ok
        fprintf('ALL TESTS PASS -- the numbers echoFitRobust prints mean what they say.\n\n');
    else
        fprintf('*** FAILURES ABOVE -- do not quote a T2 until they are fixed ***\n\n');
    end
end

% ======================================================================
% SYNTHETIC DATA, BUILT FROM THE PHYSICS RATHER THAN FROM THE FITTER
% ======================================================================

function [t, Y, E, truth] = makeEcho(varargin)
% The comb here is written out independently of echoFitRobust's revivalComb. That is
% deliberate: a test that generates data with the function under test cannot detect an error
% in that function, only in the code around it. The two agree because the formula is the same,
% not because they share a line.
    P = struct('T2', 40, 'p', 3, 'span', 120, 'nPts', 40, 'tMin', 2, ...
               'contrast', 0.0121, 'base', 0.99, 'drift', 0, ...
               'Trev', [], 'wFrac', 0.15, 'noise', 2.0e-4, 'seed', 1, 'two', true);
    for i = 1:2:numel(varargin); P.(varargin{i}) = varargin{i+1}; end

    t = linspace(P.tMin, P.span, P.nPts);
    n = numel(t);
    D = exp(-(t/P.T2).^P.p);
    if ~isempty(P.Trev)
        w = P.wFrac*P.Trev;
        K = floor(max(t)/P.Trev) + 2;
        R = zeros(1, n); R0 = 0;
        for k = 0:K
            R  = R  + exp(-((t - k*P.Trev)/w).^2);
            R0 = R0 + exp(-((0 - k*P.Trev)/w).^2);
        end
        D = D .* (R/R0);
    end
    a = -P.contrast * P.base;                  % |1> loses fluorescence as it decoheres

    y1 = P.base + P.drift*t + a*D;
    Y = y1; E = P.noise*ones(1, n);
    if P.two
        y0 = (P.base + 0.001) + P.drift*t - a*D;      % complementary: OPPOSITE sign
        Y = [y1; y0]; E = [E; P.noise*ones(1, n)];
    end
    for r = 1:size(Y,1)
        Y(r,:) = Y(r,:) + P.noise*lcgNormal(n, 7*P.seed + r);
    end
    truth = P; truth.amp = a; truth.t = t;
end

function r = lcgNormal(n, seed)
% Written-out LCG plus Box-Muller. rng() does not exist in Octave, where this is
% parse-checked, and randn('seed',k) is deprecated in MATLAB -- so a failure here is
% reproducible on whichever machine saw it.
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
% ONE struct in, one struct out -- the data rides in the same struct as the options.
    op = struct('plotIt', false, 'verbose', false, 't', t, 'y', Y, 'e', E);
    for i = 1:2:numel(varargin); op.(varargin{i}) = varargin{i+1}; end
    o = echoFitRobust(op);
end

% ======================================================================

function ok = t1_recoversMonotone()
    [t, Y, E, T] = makeEcho('T2', 40, 'p', 3, 'seed', 11);
    o = fitQuiet(t, Y, E);
    dT = abs(o.T2 - T.T2)/max(o.T2Err, eps);
    dC = abs(o.contrast - T.contrast)/max(o.contrastErr, eps);
    ok = o.ok && dT < 3 && dC < 3 && ~o.revivalsUsed;
    report(ok, sprintf('t1  monotone echo: T2 %.2f +/- %.2f vs %.1f (%.1f sigma), contrast %.3f%% vs %.3f%% (%.1f sigma)', ...
                       o.T2, o.T2Err, T.T2, dT, 100*o.contrast, 100*T.contrast, dC));
end

function ok = t2_recoversStretch()
% p must be RECOVERED, not assumed. Three baths, three exponents, one model.
    ok = true;
    for pTrue = [1 2 3]
        [t, Y, E] = makeEcho('T2', 40, 'p', pTrue, 'span', 140, 'nPts', 50, 'seed', 20+pTrue);
        o = fitQuiet(t, Y, E);
        d = abs(o.p - pTrue)/max(o.pErr, eps);
        good = o.ok && d < 3;
        ok = ok & good;
        report(good, sprintf('t2  p = %d recovered as %.3f +/- %.3f (%.1f sigma)   "%s"', ...
                             pTrue, o.p, o.pErr, d, o.pMeaning));
    end
end

function ok = t3_revivalsOrLies()
%T3 The check this whole file is built around.
%
% Part one: with the comb switched OFF, revival data must produce a BADLY WRONG T2 -- that is
% the trap being documented, and if it ever stops being wrong then the monotone model has
% quietly grown the ability to describe collapses and the comb is no longer needed. Asserting
% the size of the lie is what keeps the header honest.
%
% Part two: with the comb on, T2, the period and the implied field must all come back.
    [t, Y, E, T] = makeEcho('T2', 60, 'p', 2, 'span', 90, 'nPts', 120, 'tMin', 1, ...
                            'Trev', 18.7, 'wFrac', 0.15, 'contrast', 0.01414, ...
                            'noise', 1.5e-4, 'seed', 31);
    off = fitQuiet(t, Y, E, 'revivals', 'off');
    lied = off.T2 < 0.2*T.T2;
    report(lied, sprintf('t3a comb OFF on revival data: T2 = %.2f us against a true %.1f (%.0f%% of it) -- the trap is still there, chi2 %.1f, acf %+.2f', ...
                         off.T2, T.T2, 100*off.T2/T.T2, off.chi2Reduced, off.residAutocorr));

    on = fitQuiet(t, Y, E);
    dT  = abs(on.T2 - T.T2)/max(on.T2Err, eps);
    dTr = abs(on.revivalPeriodUs - T.Trev)/max(on.revivalPeriodErr, eps);
    got = on.ok && on.revivalsUsed && dT < 3 && dTr < 3;
    report(got, sprintf('t3b comb ON:  T2 = %.2f +/- %.2f (%.1f sigma), period %.3f +/- %.3f vs %.2f (%.1f sigma), B = %.1f G, dSSR = %.0f', ...
                        on.T2, on.T2Err, dT, on.revivalPeriodUs, on.revivalPeriodErr, T.Trev, dTr, ...
                        on.fieldGauss, on.revivalTest.dSSR));
    ok = lied && got;
end

function ok = t4_combIsNotFree()
% On monotone data the comb must be REFUSED, on all three grounds where they apply. A model
% allowed two free parameters whenever it likes will find revivals in noise, and a fitted
% 13-C field is exactly the sort of number that gets believed because it has units.
    ok = true;
    for sd = 1:4
        [t, Y, E] = makeEcho('T2', 40, 'p', 3, 'seed', 40+sd);
        o = fitQuiet(t, Y, E);
        good = o.ok && ~o.revivalsUsed && ~isempty(o.revivalTest.failed);
        ok = ok & good;
        report(good, sprintf('t4  seed %d: comb refused, dSSR %.1f against a threshold of %.2f', ...
                             sd, o.revivalTest.dSSR, o.revivalTest.threshold));
    end
    % ---- the two geometries that really do have to be refused ---------------------
    %
    % This replaces an earlier check that demanded 1.8 periods inside the sweep. That check
    % was wrong, and it was wrong in the expensive direction: revivalComb is anchored at t = 0,
    % so ONE revival away from the origin is a complete period, and refusing it threw away a
    % real 105 us revival on hardware and reported a T2 of 11.8 us instead of one over a
    % hundred. The test encoded a belief rather than a requirement, and then defended it.
    %
    % What actually has to be refused is a comb with NO measured revival in it. Two ways that
    % happens, and both are checked on the fitted values rather than guessed from a ratio.

    % (a) a period LONGER than the sweep: only the peak at the origin is in the data, which is
    %     a decay wearing a comb's clothes.
    [t, Y, E] = makeEcho('T2', 400, 'p', 2, 'span', 100, 'nPts', 60, 'tMin', 1, ...
                         'Trev', 150, 'wFrac', 0.15, 'noise', 1.5e-4, 'seed', 49);
    o = fitQuiet(t, Y, E);
    good = ~o.revivalsUsed || o.revivalTest.revivalsVisible >= 1;
    ok = ok & good;
    report(good, sprintf('t4  period longer than the sweep: comb %s -- %s', ...
                         eitherText(o.revivalsUsed, 'ACCEPTED', 'refused'), ...
                         shortReason(o)));

    % (b) a revival whose FALLING EDGE runs off the end. Its height is then whatever the last
    %     point happens to be, and the height is the only thing carrying the envelope.
    [t, Y, E] = makeEcho('T2', 300, 'p', 2, 'span', 100, 'nPts', 60, 'tMin', 1, ...
                         'Trev', 99, 'wFrac', 0.15, 'noise', 1.5e-4, 'seed', 50);
    o = fitQuiet(t, Y, E);
    edgeOk = ~o.revivalsUsed || ...
             (o.revivalsVisible*o.revivalPeriodUs + o.revivalWidthUs <= max(o.t) + 1e-9);
    ok = ok & edgeOk;
    report(edgeOk, sprintf('t4  revival cut off at the end of the sweep: comb %s -- %s', ...
                           eitherText(o.revivalsUsed, 'accepted, and its edge is inside', 'refused'), ...
                           shortReason(o)));
end

function ok = t5_jointBeatsSingle()
% Fitting both projections jointly must actually buy something. If it does not, the joint
% machinery is decoration and the code should be simpler.
    [t, Y, E] = makeEcho('T2', 40, 'p', 3, 'seed', 51);
    both = fitQuiet(t, Y, E);
    one  = fitQuiet(t, Y, E, 'useProjection', 'one');
    ok = both.ok && one.ok && both.T2Err < one.T2Err;
    report(ok, sprintf('t5  joint fit tightens T2: %.4f us against %.4f from one projection (x%.2f), |0>/|1> = %+.3f', ...
                       both.T2Err, one.T2Err, one.T2Err/max(both.T2Err,eps), both.ampRatio));
end

function ok = t6_survivesOutliers()
% Two points diluted toward the plateau -- the classic path's arm-dilution signature. On an
% echo they are more dangerous than on a Ramsey, because a point pulled toward the plateau sits
% INSIDE the data range and looks like coherence that was lost early rather than like a broken
% point.
    [t, Y, E, T] = makeEcho('T2', 40, 'p', 3, 'seed', 61);
    clean = fitQuiet(t, Y, E);
    bad = [7 15];
    Yd = Y;
    for r = 1:size(Yd,1)
        Yd(r, bad) = Yd(r, bad) + (T.base - Yd(r, bad));     % pulled all the way to no contrast
    end
    o  = fitQuiet(t, Yd, E);
    nr = fitQuiet(t, Yd, E, 'robust', false);
    moved   = abs(o.T2  - clean.T2)/max(clean.T2Err, eps);
    movedNR = abs(nr.T2 - clean.T2)/max(clean.T2Err, eps);
    found = all(ismember(bad, o.downIdx));
    ok = o.ok && found && moved < 2 && movedNR > moved;
    report(ok, sprintf(['t6  2 diluted points: found %d of 2, T2 moved %.2f sigma with robust on ' ...
                        'and %.2f sigma with it off'], numel(intersect(bad, o.downIdx)), moved, movedNR));
end

function ok = t7_errorsAreHonest()
% Repeat the same fit over many noise realisations and compare the SCATTER of T2 with the
% error bar the fit quotes. The ratio should be about 1. This is the only check that can catch
% an error bar that is simply too small, which is the failure mode that makes a fit look
% better than it is -- and on an echo it is the likely one, because T2 trades against the
% amplitude and the plateau and an error computed over T2 alone comes out far too tight.
    N = 14; v = zeros(1, N); s = zeros(1, N);
    for i = 1:N
        [t, Y, E] = makeEcho('T2', 40, 'p', 3, 'seed', 200+i);
        o = fitQuiet(t, Y, E);
        v(i) = o.T2; s(i) = o.T2Err;
    end
    ratio = std(v) / max(median(s), eps);
    ok = ratio > 0.4 && ratio < 2.5;
    report(ok, sprintf('t7  scatter of T2 over %d seeds is %.4f us against a quoted %.4f -- ratio %.2f (want 0.4 to 2.5)', ...
                       N, std(v), median(s), ratio));
end

function ok = t8_driftDegeneracy()
%T8 MEASURE the drift degeneracy rather than asserting it.
%
% This is where the table in echoFitRobust's header comes from. It is printed every run so
% that a change in the fitter shows up as a change in the numbers instead of as a comment that
% quietly stopped being true. The assertion is deliberately weak -- only that the ORDERING
% holds, a sharp decay being less affected than a shallow one -- because the exact percentages
% are a property of the data, not a specification.
    rows = {{40, 3, 120, 'p = 3, decay well inside the sweep'}, ...
            {40, 1, 120, 'p = 1, decay well inside the sweep'}, ...
            {60, 1,  80, 'p = 1, sweep only 1.3x T2'}};
    move = zeros(1, numel(rows));
    for ir = 1:numel(rows)
        r = rows{ir}; dT = []; dE = [];
        for sd = 1:5
            [t, Y, E] = makeEcho('T2', r{1}, 'p', r{2}, 'span', r{3}, 'seed', 300+sd);
            a = fitQuiet(t, Y, E);
            b = fitQuiet(t, Y, E, 'fitDrift', true);
            dT(end+1) = abs(b.T2 - a.T2)/max(a.T2, eps);                   %#ok<AGROW>
            dE(end+1) = b.T2Err/max(a.T2Err, eps);                         %#ok<AGROW>
        end
        move(ir) = median(dT);
        eTxt = sprintf('x%.2f', median(dE));
        if ~isfinite(median(dE)); eTxt = 'not computable'; end
        fprintf('      %-36s T2 moves %6.1f%%   error bar %s\n', r{4}, 100*move(ir), eTxt);
    end
    ok = move(1) < move(2) && move(2) < move(3);
    report(ok, sprintf('t8  drift degeneracy grows with a shallower decay (%.1f%% < %.1f%% < %.1f%%)', ...
                       100*move(1), 100*move(2), 100*move(3)));
end

function ok = t9_flatDataIsCalledFlat()
% A T2 fitted to a trace with no decay in it is a number with nothing under it, and it still
% comes back with an error bar. decayErrorBars is the statistic that says so, and it has to
% fire. Nothing here checks T2 itself -- there is no right answer for it, which is the point.
    [t, Y, E] = makeEcho('T2', 40, 'p', 3, 'contrast', 2e-5, 'noise', 2.0e-4, 'seed', 71);
    o = fitQuiet(t, Y, E);
    ok = o.decayErrorBars < 5 && (~isfinite(o.contrastErr) || o.contrastErr >= 0.5*o.contrast);
    report(ok, sprintf('t9  flat trace: decay spans %.1f error bars, contrast %.4f +/- %.4f %% -- not measured, and says so', ...
                       o.decayErrorBars, 100*o.contrast, 100*o.contrastErr));
end

function ok = t10_refusesJunk()
    ok = true;
    ok = ok & mustFail(@() echoFitRobust(struct('nonsense', 1)), 'an unknown option');
    ok = ok & mustFail(@() echoFitRobust(42), 'a non-struct argument');

    o = echoFitRobust(struct('plotIt', false, 'verbose', false));
    good = ~o.ok && ~isempty(o.why);
    ok = ok & good;
    report(good, sprintf('t10 empty input refused: "%s"', o.why));

    o = echoFitRobust(struct('t', 1:4, 'y', [1 1 1 1], 'plotIt', false, 'verbose', false));
    good = ~o.ok && ~isempty(o.why);
    ok = ok & good;
    report(good, sprintf('t10 four points refused: "%s"', o.why));

    % A fixed p must be reported as fixed rather than as a measurement with an error bar.
    [t, Y, E] = makeEcho('T2', 40, 'p', 3, 'seed', 81);
    o = fitQuiet(t, Y, E, 'pFixed', 3);
    good = o.ok && abs(o.p - 3) < 1e-9 && ~isfinite(o.pErr) && ~o.pFitted;
    ok = ok & good;
    report(good, sprintf('t10 pFixed = 3 held exactly, no error bar invented, T2 = %.3f +/- %.3f: "%s"', ...
                         o.T2, o.T2Err, o.pMeaning));
end

function ok = t11_readsAnExperiment()
%T11 The path from an AWGExpEcho to a fitted T2, on a stand-in for the object.
%
% This is the half of the file that never touches the fitting arithmetic, and it is where the
% quiet mistakes live -- the wrong x axis being the obvious one. AWGExpEcho publishes tau in
% topParam and 2*tau*cycles in mCurrentXAxisParam, and a T2 read off the wrong one is short by
% 2*cycles. That is a factor of 2 on a Hahn echo and 16 on an XY8, and it does not look like
% an error, it looks like a worse sample.
%
% So the mock below carries a KNOWN cycles = 4, and the test asserts that the default
% ('evolution') and xUnits = 'tau' differ by exactly 2*cycles = 8. A struct is enough: pullEcho
% reaches for named fields through tryGet, which accepts a struct or an object.
    T2 = 48; p = 2.4; cycles = 4;
    tau  = linspace(0.4, 18, 44);
    tEvo = 2*tau*cycles;
    D = exp(-(tEvo/T2).^p);
    base = 0.99; a = -0.012; noise = 2.0e-4;
    y1 = base + a*D + noise*lcgNormal(numel(tau), 91);
    y0 = base + 0.001 - a*D + noise*lcgNormal(numel(tau), 92);

    Exp = struct();
    Exp.cycles = cycles; Exp.halfPiTime = 0.018; Exp.piTime = 0.0367;
    Exp.tau = tau; Exp.doubleMeasurement = true; Exp.repeats = 2000; Exp.currIter = 12;
    Exp.mCurrentXAxisParam = struct('value', tEvo);
    Exp.topParam           = struct('value', tau);
    Exp.signalParam  = struct('value', y1, 'sterr', noise*ones(1, numel(tau)));
    Exp.signalParam2 = struct('value', y0, 'sterr', noise*ones(1, numel(tau)));
    Exp.signal = [];                              % no per-average cube: the sterr path

    oE = echoFitRobust(struct('exp', Exp, 'plotIt', false, 'verbose', false));
    oT = echoFitRobust(struct('exp', Exp, 'plotIt', false, 'verbose', false, 'xUnits', 'tau'));

    dT    = abs(oE.T2 - T2)/max(oE.T2Err, eps);
    ratio = oE.T2 / max(oT.T2, eps);
    ok = oE.ok && oT.ok && dT < 3 && abs(ratio - 2*cycles) < 0.05*2*cycles && ...
         oE.cycles == cycles;
    report(ok, sprintf(['t11 read from an experiment: T2 = %.2f +/- %.2f vs %.0f (%.1f sigma) on the ' ...
                        'evolution axis, %.3f on the tau axis -- ratio %.3f, expected 2*cycles = %d'], ...
                       oE.T2, oE.T2Err, T2, dT, oT.T2, ratio, 2*cycles));
end

function ok = t12_singleRevival()
%T12 ONE revival in the sweep -- the case a real run on this setup actually produced.
%
% A 60-point Hahn echo, tau linear 1 to 60 us at cycles = 1, so 2 to 120 us of evolution. It
% collapsed by 20 us, sat on the plateau for sixty microseconds, and then plainly revived
% around 105 us: both projections moving in OPPOSITE directions by three or four error bars
% each across fifteen points. The numbers below are matched to that trace, including its
% error bars -- it was taken at repeats = 200, so the whole decay spans only about seven of
% them.
%
% Two things have to be true and they pull in opposite directions.
%
%   THE COMB MUST BE ACCEPTED, and the period recovered. That is a real revival and refusing
%   it costs a factor of ten on T2.
%
%   T2 MUST NOT BE CLAIMED CONFIDENTLY with p free. With one revival the envelope rests on a
%   single number -- that revival's height relative to t = 0 -- and exp(-(T/T2)^p) has two
%   unknowns in it, with the comb WIDTH quietly competing for the same variance as a third.
%   So the fit has to come back either with p on a bound or with an error bar that is a large
%   fraction of T2. What it must not do is return a tight T2 and a plausible p, because that
%   would be a number with nothing behind it.
%
% And with p FIXED the period must still be recovered, which is what makes "set pFixed and
% quote T2 conditionally" a usable instruction rather than a shrug.
    T2 = 180; pTrue = 2; Trev = 105; wFrac = 0.11;
    [t, Y, E, T] = makeEcho('T2', T2, 'p', pTrue, 'span', 120, 'nPts', 60, 'tMin', 2, ...
                            'Trev', Trev, 'wFrac', wFrac, 'contrast', 0.00314, ...
                            'noise', 6.0e-4, 'seed', 61);
    o = fitQuiet(t, Y, E);
    seen = o.revivalsUsed && o.revivalsVisible == 1;
    dTr  = abs(o.revivalPeriodUs - Trev);
    gotP = seen && dTr < 0.1*Trev;
    report(gotP, sprintf('t12 single revival ACCEPTED: period %.1f vs %.0f us, %d peak in the data, dSSR %.0f, chi2 %.2f, acf %+.2f', ...
                         o.revivalPeriodUs, Trev, o.revivalsVisible, o.revivalTest.dSSR, ...
                         o.chi2Reduced, o.residAutocorr));

    onBound = o.p < 0.55 || o.p > 3.95;
    loose   = isfinite(o.T2Err) && o.T2Err > 0.25*o.T2;
    honest  = onBound || loose || ~isfinite(o.T2Err);
    report(honest, sprintf('t12 and T2 is NOT claimed confidently: T2 = %.0f +/- %.0f us (%.0f%%), p = %.2f%s', ...
                           o.T2, o.T2Err, 100*o.T2Err/max(o.T2,eps), o.p, ...
                           eitherText(onBound, ' -- on a bound', '')));

    % with p fixed, the period must still land, and T2 becomes conditional on the choice
    pf = [];
    for pv = [1 2 3]
        q = fitQuiet(t, Y, E, 'pFixed', pv);
        pf(end+1) = q.T2;                                                  %#ok<AGROW>
        good = q.revivalsUsed && abs(q.revivalPeriodUs - Trev) < 0.1*Trev;
        report(good, sprintf('t12   pFixed = %d: period %.1f us, T2 = %.0f +/- %.0f us, area = %.0f us, chi2 %.2f', ...
                             pv, q.revivalPeriodUs, q.T2, q.T2Err, q.coherenceArea, q.chi2Reduced));
        honest = honest && good;
    end
    % The point of printing all three: chi2 barely moves across them, so the data does not
    % choose. T2 sliding by this much between equally good fits IS the degeneracy, stated in
    % microseconds instead of in prose.
    fprintf('      T2 slides %.0f -> %.0f -> %.0f us across p = 1, 2, 3 on fits of equal quality (truth %g at p = %g)\n', ...
            pf(1), pf(2), pf(3), T2, pTrue);
    ok = gotP && honest;
end

% ======================================================================

function ok = mustFail(fn, what)
    ok = false;
    try
        fn();
    catch
        ok = true;
    end
    report(ok, sprintf('t10 %s raises rather than guessing', what));
end

function report(ok, what)
    if ok
        fprintf('  PASS  %s\n', what);
    else
        fprintf('  FAIL  %s\n', what);
    end
end

function s = eitherText(tf, a, b)
    if tf; s = a; else; s = b; end
end

function s = shortReason(o)
% The fitter's own words for why the comb was refused, trimmed to one clause. Better than the
% test re-deriving a reason and then agreeing with itself.
    s = o.revivalTest.failed;
    if isempty(s); s = 'accepted'; return; end
    k = strfind(s, ';');
    if ~isempty(k); s = s(1:k(1)-1); end
end
