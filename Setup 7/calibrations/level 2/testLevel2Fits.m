function ok = testLevel2Fits()
%TESTLEVEL2FITS Offline grader for the level-2 fits. No hardware, no json, ~1-2 minutes.
%
%   testLevel2Fits
%
% Synthetic data with KNOWN answers, run through the replay path (opts.data) of calibLastDelay
% and calibDetectionWindow. Expect "ALL CHECKS PASSED". Run it after any edit to either fit.
%
%   lastDelay     sharp step + creeping background (must ignore the creep)
%                 smooth exponential approach (no sharp edge; must not undershoot)
%                 no step at all (must say so and return minLastDelay, 50 ns)
%                 step past the end of the scan (must refuse to answer)
%                 a coarse 0-150 us scan with the slow hardware rise (must refuse, not invent)
%   dark time     calibDarkTime's verdicts on synthetic conditions (artefact / switch / fix)
%   window        3 exponentials + a 60 kHz artefact + reference contamination, with
%                 per-average traces and laser drift (must recover the TRUE amplitudes)
%                 1 exponential (must not invent components)

    fails = {};
    seed = 1;   % local generator below; the global RNG is not touched

    fprintf('\n================ testLevel2Fits ================\n');

    % ------------------------------------------------------------------
    fprintf('\n[1] lastDelay: step at 0.12 us, width 8 ns, plus a creep that turns over\n');
    g  = [0.005:0.01:0.3, 0.4:0.2:5]';
    g0 = 0.12; w = 0.008; tcTrue = g0 + 2.326*w;
    base = 0.0100 + 0.00035*sin(2*pi*(g + 1.5)/16.75) .* exp(-g/32);
    C  = base + 0.0003 * 0.5*(1 + erf((g - g0)/(sqrt(2)*w)));
    [Ck, seed] = perAvg(C, 20, 4e-5, 0.01, seed);
    r = calibLastDelay(struct('data', struct('gaps', g, 'contrast', mean(Ck, 2), ...
                       'contrastPerAvg', Ck), 'plot', false, 'jsonLocation', 'none'));
    fails = check(fails, strcmp(r.fit.status, 'ok'), 'step located');
    fails = check(fails, r.lastDelay >= tcTrue && r.lastDelay <= tcTrue + 0.05, ...
      sprintf('lastDelay %.3f in [%.3f, %.3f] (true 99%% point + 50 ns)', r.lastDelay, tcTrue, tcTrue + 0.05));
    fails = check(fails, abs(r.fit.g0 - g0) < 0.01, sprintf('centre %.4f vs %.4f', r.fit.g0, g0));

    fprintf('\n[2] lastDelay: smooth exponential approach, tau 50 ns, no sharp edge\n');
    g  = [0.005:0.01:0.3, 0.4:0.2:5]';
    C  = 0.0100 - 0.0003*exp(-g/0.05);
    [Ck, seed] = perAvg(C, 20, 3e-5, 0.0, seed);
    r = calibLastDelay(struct('data', struct('gaps', g, 'contrast', mean(Ck, 2), ...
                       'contrastPerAvg', Ck), 'plot', false, 'jsonLocation', 'none'));
    t99 = 0.05*log(100);
    fails = check(fails, r.lastDelay >= 0.6*t99 && r.lastDelay <= 1.8*t99, ...
      sprintf('lastDelay %.3f vs the 99%% point of the exponential %.3f (within -40%%/+80%%)', r.lastDelay, t99));

    fprintf('\n[3] lastDelay: no step\n');
    g  = (0.005:0.025:1.2)';
    C  = 0.0100 + 0*g;
    [Ck, seed] = perAvg(C, 20, 4e-5, 0.01, seed);
    r = calibLastDelay(struct('data', struct('gaps', g, 'contrast', mean(Ck, 2), ...
                       'contrastPerAvg', Ck), 'plot', false, 'jsonLocation', 'none'));
    fails = check(fails, strcmp(r.fit.status, 'nostep') && abs(r.lastDelay - 0.05) < 1e-9, ...
      sprintf('status %s, lastDelay %.3f (expect nostep, 0.05 = minLastDelay)', r.fit.status, r.lastDelay));

    fprintf('\n[4] lastDelay: step past the end of the scan\n');
    g  = (0.005:0.01:0.3)';
    C  = 0.0100 + 0.0003*0.5*(1 + erf((g - 0.33)/(sqrt(2)*0.02)));
    [Ck, seed] = perAvg(C, 20, 3e-5, 0.0, seed);
    r = calibLastDelay(struct('data', struct('gaps', g, 'contrast', mean(Ck, 2), ...
                       'contrastPerAvg', Ck), 'plot', false, 'jsonLocation', 'none'));
    fails = check(fails, isnan(r.lastDelay), sprintf('lastDelay %g (expect NaN: refuse)', r.lastDelay));

    fprintf('\n[4b] lastDelay: coarse 0-150 us scan, slow saturating rise, no resolved step\n');
    g  = (0.005:1.5:150)';
    C  = 0.0112 + 0.0065*(1 - exp(-g/25));
    [Ck, seed] = perAvg(C, 10, 5e-5, 0.0, seed);
    r = calibLastDelay(struct('data', struct('gaps', g, 'contrast', mean(Ck, 2), ...
                       'contrastPerAvg', Ck), 'plot', false, 'jsonLocation', 'none'));
    fails = check(fails, isnan(r.lastDelay) && strcmp(r.fit.status, 'coarse'), ...
      sprintf('status %s, lastDelay %g (expect coarse, NaN)', r.fit.status, r.lastDelay));

    fprintf('\n[4c] lastDelay: fine scan, step at 0.10 us on top of the slow rise\n');
    g  = (0.005:0.01:0.5)';
    C  = 0.0100 + 0.0065*(1 - 0.8*exp(-g/25)) + 0.0003*0.5*(1 + erf((g - 0.10)/(sqrt(2)*0.01)));
    [Ck, seed] = perAvg(C, 20, 3e-5, 0.0, seed);
    r = calibLastDelay(struct('data', struct('gaps', g, 'contrast', mean(Ck, 2), ...
                       'contrastPerAvg', Ck), 'plot', false, 'jsonLocation', 'none'));
    tcTrue = 0.10 + 2.326*0.01;
    fails = check(fails, r.lastDelay >= tcTrue && r.lastDelay <= tcTrue + 0.06, ...
      sprintf('lastDelay %.3f in [%.3f, %.3f]', r.lastDelay, tcTrue, tcTrue + 0.06));

    % ------------------------------------------------------------------
    fprintf('\n[5] window: 3 exponentials + 60 kHz artefact, contaminated reference, drift\n');
    aT = [0.022 0.009 0.006]; tT = [1.0 4.0 37.0]; P = 65;
    osc = struct('A', 0.0009, 'P', 16.75, 'td', 32, 'ph', 0.06);
    [d, seed] = synthWindow(aT, tT, P, osc, 30, seed);
    r = calibDetectionWindow(struct('data', d, 'plot', false, 'jsonLocation', 'none'));
    f = r.fit;
    fails = check(fails, f.n == 3 && ~isempty(f.osc), sprintf('model: %d exponentials, oscillation %d (expect 3, 1)', f.n, ~isempty(f.osc)));
    if f.n == 3
        fails = check(fails, all(abs(f.tau - tT) ./ tT < 0.08), sprintf('taus %s vs %s (8%%)', mat2str(f.tau, 3), mat2str(tT)));
        fails = check(fails, all(abs(f.aTrue - aT) ./ aT < 0.10), sprintf('true amplitudes %s vs %s (10%%)', ...
          mat2str(f.aTrue, 3), mat2str(aT)));
    end
    tpT = sum(aT .* tT) / sum(aT);
    fails = check(fails, abs(r.tauPol - tpT)/tpT < 0.05, sprintf('tauPol %.2f vs %.2f (5%%)', r.tauPol, tpT));
    if ~isempty(f.osc)
        fails = check(fails, abs(f.osc.P - osc.P)/osc.P < 0.05, sprintf('artefact period %.2f vs %.2f', f.osc.P, osc.P));
    end
    fails = check(fails, r.noiseExponent > 0.42 && r.noiseExponent < 0.62, sprintf('noise exponent %.2f (white: 0.5)', r.noiseExponent));
    % the choice, recomputed here from the TRUE parameters by brute force
    [dE, LE] = bruteChoice(aT, tT, 0.5, 2, 50, 0.05);
    fails = check(fails, abs(r.laserInitializationDuration - LE) <= 3 && abs(r.detectionDuration - dE) <= 0.5, ...
      sprintf('choice %.1f / %.2f us vs brute force on the true decay %.1f / %.2f us', ...
      r.laserInitializationDuration, r.detectionDuration, LE, dE));

    fprintf('\n[6] window: a single exponential, tau 1.5 us\n');
    [d, seed] = synthWindow(0.03, 1.5, P, [], 30, seed);
    r = calibDetectionWindow(struct('data', d, 'plot', false, 'jsonLocation', 'none'));
    fails = check(fails, r.fit.n == 1 && isempty(r.fit.osc), sprintf('model: %d exponentials, oscillation %d (expect 1, 0)', ...
      r.fit.n, ~isempty(r.fit.osc)));
    fails = check(fails, abs(r.tauPol - 1.5)/1.5 < 0.02, sprintf('tauPol %.3f vs 1.5 (2%%)', r.tauPol));

    % ------------------------------------------------------------------
    fprintf('\n[7] dark time: artefact off resonance, gone with the switch fixes, flat spin contrast\n');
    g = [0.2 0.5 1 2 3 5 7 10 15 20 30 40 60 80 100 130 160 200]';
    art = 0.0072*(1 - 0.81*exp(-g/30));
    [Z, seed] = randnLocal(numel(g), 5, seed); Z = 4e-5*Z;
    cs = {};
    cs{1} = mkCond('off', struct('frequency', 2600), g, art + Z(:,1));
    cs{2} = mkCond('off mwOff', struct('frequency', 2600, 'mwOff', true), g, art + Z(:,2));
    cs{3} = mkCond('off gap', struct('frequency', 2600, 'gapSwitch', false), g, 0.0003 + Z(:,3));
    cs{4} = mkCond('off ref', struct('frequency', 2600, 'refSwitch', true), g, Z(:,4));
    cs{5} = mkCond('on ref', struct('frequency', 2870, 'refSwitch', true), g, 0.0103 + Z(:,5));
    dd = struct('cond', [cs{:}], 'opts', struct('frequency', 2870));
    r = calibDarkTime(struct('data', dd, 'plot', false));
    fl = arrayfun(@(c) c.fit.flat, r.cond);
    fails = check(fails, isequal(fl, [false false true true true]), ...
      sprintf('flat flags %s (expect 0 0 1 1 1)', mat2str(fl)));
    fails = check(fails, abs(r.cond(1).fit.rise - 0.0072*0.81*(exp(-0.2/30) - exp(-200/30))) < 3e-4, ...
      sprintf('off-resonance rise %.4f%%', 100*r.cond(1).fit.rise));
    v = strjoin(r.verdict, ' ');
    fails = check(fails, ~isempty(strfind(v, 'UNBALANCED')) && ~isempty(strfind(v, 'refSwitch has no effect')), ...
      'unrecorded balance read as unbalanced (data before the fix), and refSwitch flagged as inert there');

    fprintf('\n[8] dark time: round 3 -- unbalanced rises, balanced flat on and off resonance\n');
    cs = {};
    cs{1} = mkCond('off UB', struct('frequency', 2600, 'balancedSequence', 0), g, art + Z(:,1));
    cs{2} = mkCond('off B',  struct('frequency', 2600, 'balancedSequence', 1), g, Z(:,2));
    cs{3} = mkCond('on B',   struct('frequency', 2870, 'balancedSequence', 1), g, 0.019 + Z(:,3));
    r = calibDarkTime(struct('data', struct('cond', [cs{:}], 'opts', struct('frequency', 2870, ...
                       'balancedSequence', 1)), 'plot', false));
    v = strjoin(r.verdict, ' ');
    fl = arrayfun(@(c) c.fit.flat, r.cond);
    fails = check(fails, isequal(fl, [false true true]), sprintf('flat flags %s (expect 0 1 1)', mat2str(fl)));
    fails = check(fails, ~isempty(strfind(v, 'was the unbalanced sequence')) && ...
      ~isempty(strfind(v, 'lastDelay is set by the MW step alone')), ...
      'verdicts: rise was the unbalanced sequence; spin contrast flat');

    fprintf('\n================================================\n');
    if isempty(fails)
        fprintf('ALL CHECKS PASSED\n\n');
    else
        fprintf('%d CHECK(S) FAILED:\n', numel(fails));
        for i = 1:numel(fails); fprintf('   %s\n', fails{i}); end
        fprintf('\n');
    end
    ok = isempty(fails);
end

% ======================================================================

function fails = check(fails, cond, msg)
    if cond
        fprintf('   pass  %s\n', msg);
    else
        fprintf('   FAIL  %s\n', msg); fails{end+1} = msg;
    end
end

function c = mkCond(label, ov, g, C)
    c = struct('label', label, 'overrides', ov, 'gaps', g, 'contrast', C, ...
               'contrastErr', 4e-5*ones(size(g)), 'contrastPerAvg', []);
end

function [Ck, seed] = perAvg(C, K, sigma, scaleJitter, seed)
% K per-average contrast curves: white noise per average (sigma is the error of the MEAN) and
% a per-average multiplicative jitter of the contrast (focus drift).
    [Z, seed] = randnLocal(numel(C), K, seed);
    [J, seed] = randnLocal(1, K, seed);
    Ck = C(:) .* (1 + scaleJitter*J) + sigma*sqrt(K)*Z;
end

function [d, seed] = synthWindow(aT, tT, P, osc, K, seed)
% Per-average signal and reference traces, 0..50 us at 10 ns, as the class records them.
    t  = (0:0.01:49.99)';
    c  = @(s) sum(aT(:)' .* exp(-s(:) ./ tT(:)'), 2);
    Cm = c(t) - c(t + P);                               % the reference carries c(t + P)
    if ~isempty(osc)
        Cm = Cm + osc.A * exp(-t/osc.td) .* sin(2*pi*t/osc.P + osc.ph);
    end
    V0 = 0.0225;
    rise = 1 - 0.08*exp(-t/0.012);                      % laser turn-on, mostly done at t = 0
    [Zs, seed] = randnLocal(numel(t), K, seed);
    [Zr, seed] = randnLocal(numel(t), K, seed);
    [Dk, seed] = randnLocal(1, K, seed);
    lvl = V0 * rise .* (1 + 0.01*Dk);                   % 1% laser drift between averages
    sig = 6e-5 * sqrt(K) / sqrt(2) * V0;                % -> 6e-5 contrast error on the mean
    S2 = lvl + sig*Zr;
    S1 = lvl .* (1 - Cm) + sig*Zs;
    d = struct('time', t, 'signal', mean(S1, 2), 'reference', mean(S2, 2), ...
               'signalPerAvg', S1, 'referencePerAvg', S2, 'nIter', K, ...
               'opts', struct('window', 50, 'extraInit', 15, 'mwBurst', 1, 'lastDelay', 0.15));
end

function [dBest, LBest] = bruteChoice(a, tau, beta, D, window, cap)
    I = @(s0, s1) sum(a .* tau .* (exp(-s0 ./ tau) - exp(-s1 ./ tau)));
    best = -Inf; dBest = NaN; LBest = NaN;
    for d = 0.05:0.05:min(window, 20)
        for L = 0.5:0.5:200
            if d > L; continue; end
            num = I(0, d) - I(L, L + d);
            if num <= 0 || I(2*L, 2*L + d) / num > cap; continue; end
            F = num / (d^beta * sqrt(D + L));
            if F > best; best = F; dBest = d; LBest = L; end
        end
    end
end

function [Z, seed] = randnLocal(m, n, seed)
% Normal deviates from a local LCG + Box-Muller, so the grader never touches the global RNG.
    N = m*n; U = zeros(2*ceil(N/2), 1);
    for i = 1:numel(U)
        seed = mod(1103515245*seed + 12345, 2^31);
        U(i) = (seed + 0.5) / 2^31;
    end
    U1 = U(1:2:end); U2 = U(2:2:end);
    R = sqrt(-2*log(U1));
    Z = [R .* cos(2*pi*U2); R .* sin(2*pi*U2)];
    Z = reshape(Z(1:N), m, n);
end
