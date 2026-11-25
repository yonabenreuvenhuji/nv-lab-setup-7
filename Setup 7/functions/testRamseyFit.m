function ok = testRamseyFit()
%TESTRAMSEYFIT Offline grader for ramseyFitRobust. No hardware, ~10 s.
%
%   testRamseyFit
%
% Every check injects a KNOWN T2*, detuning and contrast into synthetic data built from the
% same physics the fitter assumes, and asks whether it comes back. Run it before trusting a
% T2* from the setup, and certainly after touching ramseyFitRobust.
%
% The three that matter most, because each corresponds to a way of being confidently wrong:
%
%   t2  ZERO DETUNING. The case that motivated the whole file: on resonance the fringe is
%       pure hyperfine beat, and a single-cosine fit reports a 2.16 MHz detuning that is not
%       there. The fit must return delta consistent with zero AND the right T2*.
%
%   t6  OUTLIERS. Classic Ramsey runs on this setup show isolated points pulled toward
%       ratio = 1 by the per-point re-arm. Two of them must not move T2*.
%
%   t8  A SINGLE-COSINE MODEL ON TRIPLET DATA. The check that the hyperfine term is doing
%       something rather than being decorative -- nHyperfine = 1 on triplet data must fit
%       measurably WORSE than nHyperfine = 3, or the model has no content.

    fprintf('\n=== testRamseyFit ===\n');
    ok = true;
    ok = ok & t1_finiteDetuning();
    ok = ok & t2_zeroDetuning();
    ok = ok & t3_recoversT2();
    ok = ok & t4_recoversContrast();
    ok = ok & t5_jointBeatsSingle();
    ok = ok & t6_survivesOutliers();
    ok = ok & t7_errorsAreHonest();
    ok = ok & t8_hyperfineTermEarnsItsPlace();
    ok = ok & t9_refusesJunk();
    fprintf('\n');
    if ok
        fprintf('ALL TESTS PASS -- the numbers ramseyFitRobust prints mean what they say.\n\n');
    else
        fprintf('*** FAILURES ABOVE -- do not quote a T2* until they are fixed ***\n\n');
    end
end

% ======================================================================
% SYNTHETIC DATA, BUILT FROM THE PHYSICS RATHER THAN FROM THE FITTER
% ======================================================================

function [t, Y, E, truth] = makeRamsey(varargin)
    P = struct('T2', 1.4, 'delta', 0, 'A', 2.158, 'nH', 3, 'contrast', 0.012, ...
               'phi', 0, 'span', 2.5, 'perBeat', 15, 'noise', 2.2e-4, 'base', 0.995, ...
               'drift', 0, 'seed', 1, 'two', true, 'p', 2);
    for i = 1:2:numel(varargin); P.(varargin{i}) = varargin{i+1}; end

    dt = 1/(P.A * P.perBeat);
    t  = 0.005 : dt : P.span;
    n  = numel(t);

    switch P.nH
        case 3; m = [-1 0 1];
        case 2; m = [-0.5 0.5];
        otherwise; m = 0;
    end
    G = zeros(1, n);
    for k = 1:numel(m)
        G = G + cos(2*pi*(P.delta + m(k)*P.A)*t + P.phi);
    end
    G = G / numel(m);
    D = exp(-(t/P.T2).^P.p);
    amp = P.contrast * P.base / 2;                 % contrast = 2*amp/base

    y1 = P.base + P.drift*t + amp*D.*G;
    Y = y1; E = P.noise*ones(1, n);
    if P.two
        y0 = (P.base + 0.001) + P.drift*t - amp*D.*G;   % complementary: OPPOSITE sign
        Y = [y1; y0]; E = [E; P.noise*ones(1, n)];
    end
    for r = 1:size(Y,1)
        Y(r,:) = Y(r,:) + P.noise*lcgNormal(n, 7*P.seed + r);
    end
    truth = P; truth.amp = amp; truth.t = t;
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
    o = ramseyFitRobust(op);
end

% ======================================================================

function ok = t1_finiteDetuning()
    [t, Y, E, T] = makeRamsey('delta', 3.0, 'T2', 1.2, 'seed', 11);
    o = fitQuiet(t, Y, E);
    ok = o.ok && abs(o.delta - T.delta) < max(4*o.deltaErr, 0.15) && ...
                 abs(o.T2star - T.T2)  < max(4*o.T2starErr, 0.15);
    report(ok, sprintf('detuned 3.00 MHz: delta %.3f +/- %.3f, T2* %.3f +/- %.3f (want 3.000, 1.200)', ...
                       o.delta, o.deltaErr, o.T2star, o.T2starErr));
end

function ok = t2_zeroDetuning()
% THE ONE THIS FILE EXISTS FOR. On resonance the fringe is entirely the 2.158 MHz hyperfine
% beat. A single-cosine fitter reports 2.16 MHz of detuning that is not there; this must not.
    [t, Y, E, T] = makeRamsey('delta', 0, 'T2', 1.4, 'seed', 3);
    o = fitQuiet(t, Y, E);
    zeroish = o.delta < max(3*o.deltaErr, 0.25);
    t2ok = abs(o.T2star - T.T2) < max(4*o.T2starErr, 0.2);
    ok = o.ok && zeroish && t2ok;
    report(ok, sprintf(['ON RESONANCE: delta %.3f +/- %.3f (must be ~0, NOT 2.16), ' ...
                        'T2* %.3f (want 1.400)'], o.delta, o.deltaErr, o.T2star));
    % and the oscillation really is there, so this is not passing by fitting a flat line
    ok2 = o.contrast > 0.5*T.contrast;
    report(ok2, sprintf('...and it still finds the beat: contrast %.3f%% (injected %.3f%%)', ...
                        100*o.contrast, 100*T.contrast));
    ok = ok & ok2;
end

function ok = t3_recoversT2()
    ok = true;
    for T2 = [0.4 0.8 1.6]
        [t, Y, E, T] = makeRamsey('T2', T2, 'delta', 1.0, 'span', 3.0, 'seed', round(100*T2));
        o = fitQuiet(t, Y, E);
        good = o.ok && abs(o.T2star - T.T2) < max(4*o.T2starErr, 0.12*T.T2);
        report(good, sprintf('T2* = %.2f us recovered as %.3f +/- %.3f', ...
                             T.T2, o.T2star, o.T2starErr));
        ok = ok & good;
    end
end

function ok = t4_recoversContrast()
    ok = true;
    for c = [0.004 0.012 0.030]
        [t, Y, E, T] = makeRamsey('contrast', c, 'delta', 0.8, 'seed', round(1e4*c));
        o = fitQuiet(t, Y, E);
        good = o.ok && abs(o.contrast - T.contrast) < max(4*o.contrastErr, 0.15*T.contrast);
        report(good, sprintf('contrast %.2f%% recovered as %.3f +/- %.3f %%', ...
                             100*c, 100*o.contrast, 100*o.contrastErr));
        ok = ok & good;
    end
end

function ok = t5_jointBeatsSingle()
% Fitting both projections together should tighten T2* and delta, because it is twice the
% data on the same parameters. If it does not, the joint fit is not actually sharing them.
    [t, Y, E] = makeRamsey('delta', 1.5, 'T2', 1.1, 'seed', 21);
    oBoth = fitQuiet(t, Y, E, 'useProjection', 'both');
    oOne  = fitQuiet(t, Y, E, 'useProjection', 'one');
    ok = oBoth.ok && oOne.ok && oBoth.traces == 2 && oOne.traces == 1 && ...
         oBoth.T2starErr < oOne.T2starErr;
    report(ok, sprintf('joint fit tightens T2*: +/- %.4f with both vs +/- %.4f with one', ...
                       oBoth.T2starErr, oOne.T2starErr));
    % and the two projections must come back complementary
    ok2 = abs(oBoth.ampRatio + 1) < 0.15;
    report(ok2, sprintf('|0>/|1> amplitude ratio %+.3f (must be near -1)', oBoth.ampRatio));
    ok = ok & ok2;
end

function ok = t6_survivesOutliers()
% Two points dragged toward ratio = 1, which is what the classic path's per-point re-arm
% does on this setup. Least squares would chase them; the bisquare must not.
    [t, Y, E, T] = makeRamsey('delta', 0.5, 'T2', 1.3, 'seed', 31);
    clean = fitQuiet(t, Y, E);
    bad = [round(0.45*numel(t)), round(0.7*numel(t))];
    Yb = Y; Yb(:, bad) = 1.0;                       % pulled to "no contrast"
    dirty  = fitQuiet(t, Yb, E, 'robust', true);
    naive  = fitQuiet(t, Yb, E, 'robust', false);

    moved  = abs(dirty.T2star - clean.T2star) / max(clean.T2star, eps);
    movedN = abs(naive.T2star - clean.T2star) / max(clean.T2star, eps);
    found  = dirty.nDown >= 2;
    report(found, sprintf('the %d injected outliers are found (%d down-weighted)', ...
                          numel(bad), dirty.nDown));
    report(moved < 0.15, sprintf('robust T2* moves %.1f%% with them present (want < 15%%)', ...
                                 100*moved));
    % Reported, NOT asserted. The first version of this test demanded that the non-robust fit
    % move MORE, and it does not always: two outliers among eighty barely shift an envelope
    % that many points constrain. That was a made-up requirement, and a test that fails on
    % correct behaviour is worse than no test. What matters is that the robust fit finds them
    % and does not move -- which is what the two assertions above check.
    fprintf('        (for reference: non-robust moves %.1f%%, robust %.1f%%)\n', ...
            100*movedN, 100*moved);
    ok = dirty.ok && found && moved < 0.15;
end

function ok = t7_errorsAreHonest()
% The error bar has to match the actual spread of the parameter over noise realisations. This
% is the check that turns "T2* = 1.40 +/- 0.05" into a statement rather than a decoration.
    N = 15;
    T2s = zeros(1, N); errs = zeros(1, N);
    for k = 1:N
        [t, Y, E] = makeRamsey('T2', 1.2, 'delta', 1.0, 'seed', 500 + k);
        o = fitQuiet(t, Y, E);
        T2s(k) = o.T2star; errs(k) = o.T2starErr;
    end
    spread = std(T2s); rep = mean(errs);
    ratio = rep / max(spread, eps);
    ok = ratio > 0.5 && ratio < 2.5;
    report(ok, sprintf('reported T2* error %.4f vs actual spread %.4f (ratio %.2f, want 0.5-2.5)', ...
                       rep, spread, ratio));
    bias = mean(T2s) - 1.2;
    ok2 = abs(bias) < 0.1;
    report(ok2, sprintf('and T2* is unbiased: mean %.3f against a true 1.200 (bias %+.3f)', ...
                        mean(T2s), bias));
    ok = ok & ok2;
end

function ok = t8_hyperfineTermEarnsItsPlace()
% If nHyperfine = 1 fitted triplet data just as well, the hyperfine term would be decoration.
% It must fit measurably worse -- and, more importantly, it must report a WRONG detuning,
% which is the actual failure mode being guarded against.
    [t, Y, E] = makeRamsey('delta', 0, 'T2', 1.5, 'seed', 41, 'noise', 1.0e-4);
    o3 = fitQuiet(t, Y, E, 'nHyperfine', 3);
    o1 = fitQuiet(t, Y, E, 'nHyperfine', 1);
    better = o3.chi2Reduced < 0.8*o1.chi2Reduced;
    report(better, sprintf('triplet model fits better than a single cosine: chi2 %.2f vs %.2f', ...
                           o3.chi2Reduced, o1.chi2Reduced));
    % the single-cosine fit lands on the beat and calls it a detuning
    fooled = abs(o1.delta - 2.158) < 0.4;
    report(fooled, sprintf(['and the single cosine IS fooled: it reports %.3f MHz of detuning ' ...
                            'on resonance (the beat)'], o1.delta));
    ok = better && fooled;
end

function ok = t9_refusesJunk()
    ok = true;
    % A DATA problem returns ok = false with a reason; a CALLER problem raises. That split is
    % deliberate: five points is something the setup can hand you, and killing a script over
    % it is unhelpful, whereas a misspelt option is a bug in the line you just wrote.
    oFew = ramseyFitRobust(struct('t', 1:5, 'y', 1:5, 'e', 1:5, ...
                                  'plotIt', false, 'verbose', false));
    refused = ~oFew.ok && ~isempty(oFew.why);
    report(refused, 'too few points -> ok = false with a reason, not a crash');
    ok = ok & refused;
    ok = ok & mustFail(@() ramseyFitRobust(struct('t', 1:30, 'y', 1:30, 'e', 1:30, ...
                                                  'nonsense', 1)), ...
                       'unknown option refused rather than ignored');
    % A flat trace must not come back with a confident T2*.
    t = 0.005:0.03:2.5;
    y = 0.995 + 2e-4*lcgNormal(numel(t), 99);
    o = fitQuiet(t, y, 2e-4*ones(size(t)));
    weak = ~o.ok || o.contrast < 0.002 || o.contrastErr > 0.3*o.contrast;
    report(weak, sprintf(['a flat trace does not yield a confident fringe (contrast %.4f%% ' ...
                          '+/- %.4f%%)'], 100*o.contrast, 100*o.contrastErr));
    ok = ok & weak;
end

% ======================================================================

function ok = mustFail(fn, what)
    ok = false;
    try
        fn();
    catch
        ok = true;
    end
    report(ok, what);
end

function report(ok, what)
    if ok
        fprintf('  PASS  %s\n', what);
    else
        fprintf('  FAIL  %s\n', what);
    end
end
