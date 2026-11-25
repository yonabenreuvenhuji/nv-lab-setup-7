function ok = testWalkAudit()
%TESTWALKAUDIT Offline grader for walkAudit's difference decomposition. No hardware, ~1 s.
%
%   testWalkAudit
%
% walkAudit's whole conclusion rests on two numbers -- the contrast scale and the
% gate-tail tilt -- so those two numbers have to be shown to mean what the printout says
% they mean BEFORE the audit is run on the setup. Every check here injects a known effect
% into synthetic data and asks whether the decomposition recovers it.
%
% The dangerous failure is not a crash. It is a decomposition that returns a confident
% -0.5 contrast scale on data that has no contrast loss at all, because a tilt or a
% baseline shift leaked into that column. Two of these tests exist for exactly that.

    fprintf('\n=== testWalkAudit ===\n');
    ok = true;
    ok = ok & t1_recoversPureScale();
    ok = ok & t2_recoversPureTilt();
    ok = ok & t3_separatesBoth();
    ok = ok & t4_identicalIsZero();
    ok = ok & t4b_identicalWithRealNoise();
    ok = ok & t5_offsetDoesNotLeak();
    ok = ok & t6_noTiltColumnForEsr();
    ok = ok & t7_errorsAreHonest();
    fprintf('\n');
    if ok
        fprintf('ALL TESTS PASS -- the two reported numbers mean what walkAudit claims.\n\n');
    else
        fprintf('*** FAILURES ABOVE -- do not run walkAudit until they are fixed ***\n\n');
    end
end

% ======================================================================

function [t, c, cE] = classicRabi(n)
% A plausible classic Rabi: damped oscillation on a ratio near 1, with error bars the size
% the real thing has (median sterr ~1e-3 at repeats 500, averages 2).
    t = linspace(0.02, 1.0, n);
    c = 1 - 0.06 * exp(-t/1.2) .* (1 - cos(2*pi*5.5*t)) / 2;
    cE = 1.0e-3 * ones(1, n);
end

function ok = t1_recoversPureScale()
% Inject a known contrast scale and nothing else. This is the halved-contrast case.
    [t, c, cE] = classicRabi(51);
    kTrue = -0.45;
    w = mean(c) + (1 + kTrue) * (c - mean(c));
    g = (max(t) - t) / max(t);
    p = walkAuditDecompose(c, cE, w, cE, g);
    ok = p.ok && abs(p.scale - kTrue) < 0.01 && abs(p.tilt) < 3*p.tiltErr;
    report(ok, sprintf('pure scale %+.2f recovered as %+.3f, tilt stays flat (%+.4f +/- %.4f)', ...
                       kTrue, p.scale, p.tilt, p.tiltErr));
end

function ok = t2_recoversPureTilt()
% Inject a ramp along the gate tail and nothing else. The scale column must NOT absorb it
% -- if it does, a gate leak gets reported as a contrast loss and the wrong fix follows.
    [t, c, cE] = classicRabi(51);
    g = (max(t) - t) / max(t);
    mTrue = -0.004;
    w = c + mTrue * (g - mean(g));
    p = walkAuditDecompose(c, cE, w, cE, g);
    ok = p.ok && abs(p.tilt - mTrue) < 2e-4 && abs(p.scale) < 3*max(p.scaleErr, eps);
    report(ok, sprintf('pure tilt %+.4f recovered as %+.4f, scale stays zero (%+.4f +/- %.4f)', ...
                       mTrue, p.tilt, p.scale, p.scaleErr));
end

function ok = t3_separatesBoth()
% Both at once. The two columns are not orthogonal on real Rabi data -- a damped
% oscillation has a downward trend of its own -- so this is the check that the
% orthogonalisation actually works rather than just being claimed in a comment.
    [t, c, cE] = classicRabi(51);
    g = (max(t) - t) / max(t);
    kTrue = -0.30; mTrue = -0.003;
    w = mean(c) + (1 + kTrue)*(c - mean(c)) + mTrue*(g - mean(g));
    p = walkAuditDecompose(c, cE, w, cE, g);
    ok = p.ok && abs(p.scale - kTrue) < 0.02 && abs(p.tilt - mTrue) < 5e-4;
    report(ok, sprintf('scale AND tilt together: %+.3f (want %+.2f), %+.4f (want %+.4f)', ...
                       p.scale, kTrue, p.tilt, mTrue));
end

function ok = t4_identicalIsZero()
% Two identical traces, with the error declared as ZERO so there is nothing for the dilution
% correction to add back. This is the pure algebra check.
%
% The earlier version of this test passed identical NOISELESS traces while declaring
% cE = 1e-3, and after the dilution correction was added it started failing with scale
% +0.004. The test was wrong, not the code: a trace with no noise in it whose sterr claims
% 1e-3 is not a situation that occurs, and the correction correctly added back the dilution
% that sterr said was present. Real data has noise matching its error bars, which is what
% t4b checks.
    [t, c, ~] = classicRabi(41);
    g = (max(t) - t) / max(t);
    z = zeros(size(c));
    p = walkAuditDecompose(c, z, c, z, g);
    ok = p.ok && abs(p.scale) < 1e-9 && abs(p.tilt) < 1e-9;
    report(ok, sprintf('identical noiseless traces give scale %.2e and tilt %.2e (exactly zero)', ...
                       p.scale, p.tilt));
end

function ok = t4b_identicalWithRealNoise()
% Two runs of the SAME thing, each with real noise matching its declared error. The corrected
% scale must average to zero -- this is the case the dilution correction exists for, and the
% case that read -0.446 before it existed.
    [t, c, cE] = classicRabi(34);
    g = (max(t) - t) / max(t);
    N = 300;
    raw = zeros(1, N); cor = zeros(1, N);
    big = 4*cE;                      % noise heavy enough that dilution really bites
    for i = 1:N
        aa = c + big.*randnLike(numel(t), 5*i);
        bb = c + big.*randnLike(numel(t), 5*i + 1);
        p = walkAuditDecompose(aa, big, bb, big, g);
        raw(i) = p.scaleRaw; cor(i) = p.scale;
    end
    fprintf('     uncorrected mean %+.3f, corrected mean %+.3f\n', mean(raw), mean(cor));
    ok = mean(raw) < -0.05 && abs(mean(cor)) < 0.05;
    report(ok, sprintf(['dilution bias is real (%+.3f uncorrected) and the correction removes ' ...
                        'it (%+.3f)'], mean(raw), mean(cor)));
end

function ok = t5_offsetDoesNotLeak()
% A common-mode shift in the ratio -- laser power, reference level -- must land entirely
% in the offset. If it leaked into scale, every run with a slightly different laser would
% be reported as a contrast change and the audit would cry wolf on every pair.
    [t, c, cE] = classicRabi(51);
    g = (max(t) - t) / max(t);
    w = c + 0.01;                        % 1% baseline shift, ten times any real contrast
    z = zeros(size(c));                  % zero declared error: no dilution to correct
    p = walkAuditDecompose(c, z, w, z, g);
    ok = p.ok && abs(p.scale) < 1e-9 && abs(p.tilt) < 1e-9 && abs(p.offset - 0.01) < 1e-9;
    report(ok, sprintf('a 1%% baseline shift lands in the offset (%.4f), not in scale (%.2e)', ...
                       p.offset, p.scale));
end

function ok = t6_noTiltColumnForEsr()
% ESR passes g = [], and then tilt must come back NaN rather than 0. Zero would print as
% "gate-tail tilt +0.0000, CLEAN" and read as evidence about a gate ESR does not have.
    n = 41;
    f = linspace(2730, 2790, n);
    c = 1 - 0.03 ./ (1 + ((f - 2760)/5).^2);
    cE = zeros(1, n);                    % noiseless synthetic pair, so declare it
    w = mean(c) + 0.8*(c - mean(c));
    p = walkAuditDecompose(c, cE, w, cE, []);
    ok = p.ok && isnan(p.tilt) && abs(p.scale - (-0.2)) < 0.01;
    report(ok, sprintf('ESR (no tilt column): scale %+.3f, tilt is NaN not zero (%s)', ...
                       p.scale, mat2str(isnan(p.tilt))));
end

function ok = t7_errorsAreHonest()
% The error bar on the scale has to track the noise: a null case must sit inside 2 sigma
% about 95% of the time, and the reported error must match the ACTUAL spread of the
% coefficient over many noise realisations. That last check is the strongest one in this
% file and it caught a real mistake -- in the first version of this test the noise went
% into the walk trace only while both traces were declared to have error cE, so the
% combined sigma was overstated by sqrt(2), the reduced chi2 came back at 0.50 instead of
% 1, and the reported error was 1.38x the true spread. The decomposition was right and the
% test was wrong. Noise now goes into BOTH traces, which is also what happens on the
% setup: classic and walk each have their own photon statistics.
    [t, c, cE] = classicRabi(51);
    g = (max(t) - t) / max(t);
    nTrial = 400;
    s1 = zeros(1, nTrial); e1 = zeros(1, nTrial); ch = zeros(1, nTrial);
    for i = 1:nTrial
        cc = c + cE .* randnLike(numel(t), 2*i);
        ww = c + cE .* randnLike(numel(t), 2*i + 1);
        p = walkAuditDecompose(cc, cE, ww, cE, g);
        s1(i) = p.scale; e1(i) = p.scaleErr; ch(i) = p.chi2;
    end
    inside = mean(abs(s1) < 2*e1);
    spread = std(s1);

    ok0 = abs(median(ch) - 1) < 0.25;
    report(ok0, sprintf('reduced chi2 of a null decomposition is %.2f (want ~1)', median(ch)));

    ok1 = inside > 0.88 && inside < 1.0;
    report(ok1, sprintf('null case sits inside 2 sigma %.0f%% of the time (want 88-100%%)', ...
                        100*inside));

    ok2 = abs(mean(e1) - spread)/max(spread, eps) < 0.20;
    report(ok2, sprintf('mean reported error %.4f vs actual spread %.4f (%.0f%% apart)', ...
                        mean(e1), spread, 100*abs(mean(e1)-spread)/max(spread,eps)));

    % scaling: ten times the noise, ten times the error bar
    p10 = walkAuditDecompose(c + 10*cE.*randnLike(numel(t), 91), 10*cE, ...
                             c + 10*cE.*randnLike(numel(t), 92), 10*cE, g);
    p01 = walkAuditDecompose(c +    cE.*randnLike(numel(t), 91), cE, ...
                             c +    cE.*randnLike(numel(t), 92), cE, g);
    ratio = p10.scaleErr / max(p01.scaleErr, eps);
    ok3 = ratio > 7 && ratio < 14;
    report(ok3, sprintf('10x the noise gives %.1fx the error bar (want ~10)', ratio));

    ok = ok0 & ok1 & ok2 & ok3;
end

% ======================================================================

function r = randnLike(n, seed)
% Deterministic standard normals, so a failure here is reproducible on the machine that
% failed. rng() behaves differently across releases and does not exist in Octave, and
% randn('seed', k) is deprecated in one and absent in the other, so the generator is
% written out: a plain 32-bit LCG feeding Box-Muller.
%
% NOT mean-removed and NOT variance-normalised, deliberately. Normalising would remove
% the very fluctuation whose size test t7 is checking the reported error bar against, and
% the test would then pass on an error bar that is systematically too small.
    m = 2^31 - 1;
    x = mod(1103515245 * (seed + 12345) + 12345, m);
    u = zeros(1, 2*n);
    for i = 1:2*n
        x = mod(1103515245 * x + 12345, m);
        u(i) = (x + 0.5) / m;
    end
    r = sqrt(-2*log(u(1:n))) .* cos(2*pi*u(n+1:2*n));
end

function report(ok, what)
    if ok
        fprintf('  PASS  %s\n', what);
    else
        fprintf('  FAIL  %s\n', what);
    end
end
