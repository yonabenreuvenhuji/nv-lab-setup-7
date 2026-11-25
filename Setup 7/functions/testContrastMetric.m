function ok = testContrastMetric()
%TESTCONTRASTMETRIC Why ramseyWalkDiag's contrast number was not a measurement. ~2 s.
%
%   testContrastMetric
%
% ===================================================================
% WHAT THIS RETRACTS
% ===================================================================
%
% ramseyWalkDiag reported, on 34 taus spanning 2.5 us at 300 passes:
%
%       classic     0.001107
%       walk P=1    0.0004962      <-- "0.45 of classic"
%       walk P=50   0.0003453
%       walk P=228  0.0003453      <-- "0.31"
%
% and that 0.45 is what set off four nights of investigation into the Ramsey walk. It was
% not a measurement of contrast. It was the sampling scatter of a badly chosen statistic.
%
% THE STATISTIC. ramseyWalkDiag's contrastOf did
%
%       resid = v - movMedian(v, 5);
%       c     = 2 * 1.4826 * mad(resid);
%
% i.e. subtract a 5-point moving median and take the robust spread of what is left. The
% intent was "remove slow drift, measure the fringe". On that grid it does the opposite.
%
% THE ARITHMETIC NOBODY DID. 34 points over 2.5 us is a 0.075 us step. The hyperfine beat
% is 1/2.16 MHz = 0.463 us, so there are 6.17 points per beat, and a FIVE-point window spans
% 0.375 us -- 81% of one full oscillation. A moving median whose window is most of a period
% TRACKS that period. So the fringe is largely subtracted along with the drift, and what
% survives into the residual is a small remainder whose size depends on how the sampling
% grid happens to line up with the oscillation phase.
%
% Two consequences, both measured below rather than argued:
%
%   1. At FIXED true amplitude, with nothing changing but the noise draw, the statistic
%      scatters by about 28% and a high draw over a low draw is a factor of ~2.6. The
%      reported 0.001107 against 0.0004962 is a factor of 2.2 -- INSIDE the statistic's own
%      noise. There was never a signal to explain.
%
%   2. At fixed true amplitude and NO noise at all, moving the fringe frequency from 2.00 to
%      2.35 MHz changes the statistic by a factor of ~4.6, and at 2.16 MHz it returns only
%      about 31% of the true peak-to-peak. Its response is dominated by grid-versus-phase
%      alignment, not by contrast.
%
% So the "0.45, sliding to 0.31 with passesPerRun" was one statistic wandering inside its own
% 28% scatter across four runs, and the reason the error bars looked correct while the
% contrast looked halved is the plainest possible one: the error bars WERE real measurements
% and the contrast was not.
%
% WHAT REPLACED IT. walkAuditDecompose projects (walk - classic) onto (classic - mean
% classic) and returns a scale WITH an error bar, Birge-scaled. On the same question it gave
% Rabi 98% +/- 1.7%, ESR 99% +/- 1.9%, and Ramsey 102% +/- 4.1%. testWalkAudit grades it
% against injected effects. This file exists so that contrastOf is not quietly reused.
%
% KEEP THE LESSON, NOT JUST THE FIX: a summary statistic with no error bar is not a
% measurement, and the moment two of them get compared as a ratio the missing error bar
% becomes the whole result.

    fprintf('\n=== testContrastMetric ===\n');
    fprintf('Grading ramseyWalkDiag''s retired contrast statistic against KNOWN contrast.\n\n');
    ok = true;

    % ramseyWalkDiag's exact grid: step = min(1/(6*2.16), 0.075) = 0.075, tauMax 2.5
    t = 0.005 : 0.075 : 2.5;
    beat = 1/2.16;
    A = 5.5e-4;                 % amplitude -> 1.1e-3 peak to peak, matching the real fringe
    noise = 2.2e-4;             % matching the reported median sterr

    fprintf('ramseyWalkDiag grid: %d points, %.4f us step, %.2f points per beat.\n', ...
            numel(t), median(diff(t)), beat/median(diff(t)));
    fprintf('A 5-point median window spans %.4f us = %.0f%% of one beat.\n\n', ...
            5*0.075, 100*5*0.075/beat);

    ok = ok & t1_scatterAtFixedAmplitude(t, A, noise);
    ok = ok & t2_dependsOnFrequencyNotAmplitude(t, A);
    ok = ok & t3_thatGridHadNoPowerAtAll(t, A, noise);

    fprintf('\n');
    if ok
        fprintf(['CONCLUSION, in two parts:\n' ...
                 '  1. The retired statistic''s own scatter covers the effect it reported, and it\n' ...
                 '     responds to grid-versus-phase alignment more than to contrast.\n' ...
                 '  2. Even a correct metric has NO POWER on that grid and shot count. So the\n' ...
                 '     Ramsey walk deficit was never a measurement, and no reanalysis of that\n' ...
                 '     data can rescue one. The 30-point, 3000-repeat configuration in\n' ...
                 '     ramseyGateTest resolves 12%%, and it says walk is at 102%%.\n\n']);
    else
        fprintf('*** a check above did not behave as documented -- read it before trusting this\n\n');
    end
end

% ======================================================================

function ok = t1_scatterAtFixedAmplitude(t, A, noise)
% The decisive one. Fix the contrast, vary only the noise, and see how much the statistic
% moves. If its own scatter spans the reported effect, the effect was never measured.
    n = numel(t);
    N = 400;
    c = zeros(1, N);
    for k = 1:N
        v = 1 + A*exp(-t/1.8).*cos(2*pi*2.16*t) + noise*lcgNormal(n, k);
        c(k) = oldContrast(v);
    end
    rel = std(c)/mean(c);
    hiLo = pct(c, 95) / max(pct(c, 5), eps);
    fprintf('  1. SAME contrast, noise only:  mean %.6f, sd %.6f (%.0f%% relative)\n', ...
            mean(c), std(c), 100*rel);
    fprintf('     95th over 5th percentile = %.2f\n', hiLo);
    fprintf('     ramseyWalkDiag reported 0.001107 vs 0.0004962, a factor of 2.2.\n');
    ok = rel > 0.15 && hiLo > 1.8;
    report(ok, sprintf(['the statistic''s OWN scatter (%.0f%%, %.1fx high-to-low) covers the ' ...
                        '2.2x that was reported as a finding'], 100*rel, hiLo));
end

function ok = t2_dependsOnFrequencyNotAmplitude(t, A)
% No noise at all. Only the fringe FREQUENCY changes, by less than the uncertainty in the
% detuning. A contrast metric must not care.
    fs = linspace(2.00, 2.35, 40);
    c = zeros(size(fs));
    for k = 1:numel(fs)
        c(k) = oldContrast(1 + A*exp(-t/1.8).*cos(2*pi*fs(k)*t));
    end
    at216 = oldContrast(1 + A*exp(-t/1.8).*cos(2*pi*2.16*t));
    span = max(c)/max(min(c), eps);
    frac = at216 / (2*A);
    fprintf('\n  2. SAME contrast, NO noise, f from 2.00 to 2.35 MHz:\n');
    fprintf('     statistic ranges %.6f to %.6f -- a factor of %.2f\n', min(c), max(c), span);
    fprintf('     at 2.16 MHz it returns %.0f%% of the true peak-to-peak %.6f\n', ...
            100*frac, 2*A);
    ok = span > 2 && frac < 0.6;
    report(ok, sprintf(['the statistic tracks grid-vs-phase alignment (%.1fx over a small ' ...
                        'frequency change) and recovers only %.0f%% of the real contrast'], ...
                       span, 100*frac));
end

function ok = t3_thatGridHadNoPowerAtAll(t, A, noise)
% The conclusion this file exists for, and it is stronger than "the old metric was bad".
%
% Put the REPLACEMENT metric on ramseyWalkDiag's own grid and shot count and ask it the same
% question. It comes back unbiased -- 0 for an unchanged pair, -0.55 for a halved one -- and
% with an error bar so large that the two cases DO NOT SEPARATE. That is the honest verdict on
% the original run: 34 points, 300 passes and one average could not have distinguished a
% halved Ramsey contrast from an unchanged one by any means. There was no measurement to
% interpret, and four nights went into interpreting it.
%
% (The first version of this check asserted the replacement WOULD separate them here. It does
% not, and the assertion was wrong -- an estimator that reported a clean separation on this
% data would be the broken one.)
    n = numel(t);
    base = A*exp(-t/1.8).*cos(2*pi*2.16*t);
    e = noise*ones(1, n);
    N = 200;
    sSame = zeros(1, N); sHalf = zeros(1, N); eSame = zeros(1, N);
    for k = 1:N
        cRun  = 1 + base + noise*lcgNormal(n, 3*k);
        wSame = 1 + base + noise*lcgNormal(n, 3*k+1);
        wHalf = 1 + 0.45*base + noise*lcgNormal(n, 3*k+2);
        pS = walkAuditDecompose(cRun, e, wSame, e, [], noise^2);
        pH = walkAuditDecompose(cRun, e, wHalf, e, [], noise^2);
        sSame(k) = pS.scale; sHalf(k) = pH.scale; eSame(k) = pS.scaleErr;
    end
    sep = abs(mean(sSame) - mean(sHalf)) / sqrt(std(sSame)^2 + std(sHalf)^2);
    fprintf('\n  3. The REPLACEMENT metric, on ramseyWalkDiag''s grid and shot count:\n');
    fprintf('     unchanged pair : %+.3f, spread %.3f, reported error %.3f\n', ...
            mean(sSame), std(sSame), mean(eSame));
    fprintf('     55%% removed    : %+.3f, spread %.3f\n', mean(sHalf), std(sHalf));
    fprintf('     the two cases separate by only %.1f sigma\n', sep);
    unbiased = abs(mean(sSame)) < 0.25 && abs(mean(sHalf) + 0.55) < 0.25;
    powerless = sep < 3;
    report(unbiased, 'the replacement is UNBIASED even here (0 for unchanged, -0.55 for halved)');
    report(powerless, sprintf(['...but on THAT grid nothing separates them (%.1f sigma), so the ' ...
                               'original run had no power to answer its own question'], sep));
    ok = unbiased && powerless;
end

% ======================================================================

function c = oldContrast(v)
% ramseyWalkDiag's contrastOf, reproduced verbatim so this file grades the real thing.
    resid = v - movMedian(v, 5);
    c = 2 * 1.4826 * median(abs(resid - median(resid)));
end

function m = movMedian(v, w)
    n = numel(v); m = zeros(1, n); h = floor(w/2);
    for i = 1:n
        lo = max(1, i-h); hi = min(n, i+h);
        m(i) = median(v(lo:hi));
    end
end

function r = lcgNormal(n, seed)
% Written-out LCG plus Box-Muller: rng() is absent in Octave and randn('seed',k) is
% deprecated in MATLAB, and a test that cannot be reproduced is not much of a test.
    m = 2^31 - 1;
    x = mod(1103515245*(seed + 12345) + 12345, m);
    u = zeros(1, 2*n);
    for i = 1:2*n
        x = mod(1103515245*x + 12345, m);
        u(i) = (x + 0.5)/m;
    end
    r = sqrt(-2*log(u(1:n))) .* cos(2*pi*u(n+1:2*n));
end

function q = pct(v, p)
    s = sort(v(:));
    k = max(1, min(numel(s), round(p/100 * numel(s))));
    q = s(k);
end

function report(ok, what)
    if ok
        fprintf('     PASS  %s\n', what);
    else
        fprintf('     FAIL  %s\n', what);
    end
end
