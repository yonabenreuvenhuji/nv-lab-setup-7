function results = calibDetectionWindow(opts)
%CALIBDETECTIONWINDOW Polarisation decay -> tauPol, detectionDuration, laserInitializationDuration.
%
%   results = calibDetectionWindow()
%   results = calibDetectionWindow(opts)
%   results = calibDetectionWindow(struct('data', r))   % refit a saved result, no hardware
%
% LEVEL 2 CALIBRATION, run after calibLastDelay. Writes nothing: results.toWrite holds the four
% numbers, and level2.m cell 2.3 writes them with expPropsWrite.
%
% ===================================================================
% THE MEASUREMENT
% ===================================================================
%
% AWGExpDigitizerInitDuration: saturating MW burst, lastDelay, then ONE long laser window
% recorded sample by sample; then the same with no MW as the reference. The contrast
% C(t) = 1 - S1(t)/S2(t) is the spin polarisation deficit left after t of green light.
%
% A SATURATING burst rather than a pi pulse, so this needs no piTime (a level 4/5 number). It
% gives about half the contrast of a pi pulse but the same SHAPE, and every answer here comes
% from the shape. It also erases history: whatever deficit was left from the previous repeat,
% the burst equalises the driven pair again, so every signal window starts from the same state.
%
% ===================================================================
% THE DECAY MODEL, AND WHY
% ===================================================================
%
% The old analysis fitted a straight line to log C, i.e. one exponential. On 2026-11-03 that
% returned tau_pol = 11.3 us at chi2_nu = 1370 (against the true point noise) and the line
% missed the data everywhere: too slow for the first microseconds, too fast for the tail. The
% data are not one exponential, and fitting log C hid it because the 5000 tail samples outvote
% the few hundred that carry the fast decay.
%
% What the data do support, tested on that run:
%
%     model (from light-on, 0 to 50 us)       chi2_nu
%     1 exponential                            1370
%     stretched exponential                    119
%     2 exponentials                           56
%     3 exponentials                           24
%     3 exponentials + damped oscillation      2.9    (white residuals)
%
% So: a SUM OF EXPONENTIALS, as many as the data need (1 to 3, chosen by BIC), plus an
% optional damped sinusoid as a NUISANCE term. A sum of exponentials is the physical shape
% for an ensemble whose NVs see different laser intensity (centre vs wings of the spot, depth
% of focus): each sub-population re-polarises at its own rate. On that run the components were
% about 1.0 us (60%), 4 us (24%) and 37 us (17%).
%
% The oscillation is NOT polarisation: 16.8 us period (60 kHz), decaying in ~32 us, 0.09% of
% contrast, present in the SIGNAL window only -- so it is triggered by the MW burst (amplifier
% supply, pickup or thermo-mechanical ringing). It is fitted so it cannot bend the decay, and
% then left out of everything downstream. To confirm it is an artefact, re-run with
% opts.frequency far off resonance: the spin contrast goes, the oscillation should stay.
%
% Fitting: for fixed time constants the amplitudes are linear, so they are solved exactly
% (variable projection) and fminsearch only searches the 1-5 time constants, from several
% starts. No toolbox needed.
%
% THE REFERENCE IS CONTAMINATED, AND IT IS CORRECTED. The reference window opens
% laserInitializationDuration (window + extraInit) of light after the signal's MW, so it still
% carries that much-decayed deficit:  C_measured(t) = c(t) - c(t + P).  Because the burst erases
% history, that is exact to first order, and for an exponential component it only scales the
% amplitude by (1 - exp(-P/tau)). The time constants are unbiased; the amplitudes are divided by
% that factor. It matters only for the slow component (17% larger at P = 65 us, tau = 37 us).
%
% ===================================================================
% tauPol
% ===================================================================
%
% With several components there is no single time constant, so tauPol has to be DEFINED. It is
% the MEAN re-polarisation time, the area under c(t)/c(0):
%
%       tauPol = sum(a_i tau_i) / sum(a_i)
%
% That equals tau for a single exponential, so every rule written as "k x tauPol" keeps its
% meaning. It is also the quantity those rules need: the contrast a window loses to the start-
% of-window transient is the area of that transient, which is what esrWindowPick's
% 10 x tauPol floor is about. The 1/e time is printed too; it is much shorter (~2 us) because
% it ignores the tail.
%
% ===================================================================
% detectionDuration AND laserInitializationDuration ARE CHOSEN TOGETHER
% ===================================================================
%
% In every pulsed class (Rabi, Ramsey, Echo, T1, ...) laserInitializationDuration is the TOTAL
% light per half-repeat and the detection window sits inside it:
%       initDuration = laserInitializationDuration - detectionDuration
% So the old "SNR per unit time" panel, which charged 2 us of sequence per us of window, was
% optimising a cost that does not exist. The window's real cost is different: the reference
% window opens L = laserInitializationDuration of light after the signal's MW, so it carries
% the signal's residual deficit integrated over ITS OWN length -- and the slow component makes
% that contamination grow with the window. In a balanced half-repeat of dark D and light L:
%
%       signal      I(0, d) - I(L, L + d)        I(a,b) = integral of c from a to b
%       noise       d^beta                       (beta = 0.5 for white noise; measured if possible)
%       time        D + L                        per half-repeat
%       FOM         [I(0,d) - I(L,L+d)] / ( d^beta * sqrt(D + L) )
%
% and the constraint is CARRY-OVER: the next point's signal window opens 2L after this point's
% MW, so it inherits I(2L, 2L+d) of this point's deficit. In WALK mode consecutive repeats are
% different points, so that is cross-talk between neighbouring tau values. The default cap is
% 5% of the point's own contrast (opts.carryOverMax). In CLASSIC mode it only rescales the
% contrast, and the unconstrained optimum -- much shorter -- is printed for that case.
%
% The recommendation is the shortest L whose carry-over is under the cap, with the d that
% maximises the FOM there. The table it prints shows the price of a tighter or looser cap.
%
% ===================================================================
% ERROR BARS
% ===================================================================
%
% Contrast is formed PER AVERAGE and its error is the standard error across averages. Laser and
% focus drift are common to signal and reference and cancel in the ratio; the old propagation of
% separate S1 and S2 errors charged that drift as noise, ~4x too large. A replay of a result
% without per-average data falls back to the local point-to-point scatter.
%
% ===================================================================
% OPTIONS (struct fields, all optional)
%   data           previous result struct: refit, touch no hardware
%   jsonLocation   expProperties.json (default setup7Paths().json)
%   window         detection window in us (default 50)
%   extraInit      us of laser after the window (default 15)
%   mwBurst        saturating burst, us (default 1)
%   lastDelay      MW-to-window gap, us (default: the json value, else 0.2)
%   repeats, averages (100, 100), frequency (MHz, 2870), amplitude (dBm, -10), fixDelays (1)
%   carryOverMax   allowed carry-over into the next point (default 0.05)
%   darkPerHalf    us of dark per half-repeat in the target experiments: trigger + MW/free
%                  evolution + lastDelay (default 2). Set it larger for T1-like sequences.
%   maxExp         most exponentials tried (default 3)
%   fitOscillation try the damped-sinusoid nuisance term (default true)
%   plot           default true
%
% OUTPUT (main fields)
%   tauPol, tauPolErr                 mean re-polarisation time, us
%   detectionDuration                 recommended, us
%   laserInitializationDuration       recommended, us
%   toWrite                           struct of the four json fields
%   fit                               the decay model (components, oscillation, chi2, ...)
%   trade                             the carry-over table
%   time, signal, reference, contrast, contrastErr, contrastPerAvg (single), opts

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('data', [], 'jsonLocation', '', 'window', 50, ...
                        'extraInit', 15, 'mwBurst', 1, 'lastDelay', [], ...
                        'repeats', 100, 'averages', 100, 'frequency', 2870, ...
                        'amplitude', -10, 'fixDelays', 1, ...
                        'nChannels', 2, 'AWGorSRSswitch', 0, 'GIMeas', 0, ...
                        'timeDelay', 0, 'amplitudeDiff', [0 0], ...
                        'carryOverMax', 0.05, 'darkPerHalf', 2, 'maxExp', 3, ...
                        'fitOscillation', true, 'noiseExponent', [], 'plot', true));

    if isempty(opts.jsonLocation)
        try
            opts.jsonLocation = setup7Paths().json;
        catch
            opts.jsonLocation = '';
        end
    end

    S1 = []; S2 = []; Ck = [];
    if ~isempty(opts.data)
        % ================= REPLAY =========================================
        d = opts.data;
        if isfield(d, 'opts')
            keep = {'window', 'extraInit', 'mwBurst', 'lastDelay', 'repeats', 'averages', ...
                    'frequency', 'amplitude'};
            for i = 1:numel(keep)
                if isfield(d.opts, keep{i}) && ~isempty(d.opts.(keep{i}))
                    opts.(keep{i}) = d.opts.(keep{i});
                end
            end
        end
        t  = d.time(:);
        s1 = d.signal(:);  s2 = d.reference(:);
        if isfield(d, 'signalPerAvg') && ~isempty(d.signalPerAvg)
            S1 = double(d.signalPerAvg); S2 = double(d.referencePerAvg);
            Ck = 1 - S1 ./ S2;
        elseif isfield(d, 'contrastPerAvg') && ~isempty(d.contrastPerAvg)
            Ck = double(d.contrastPerAvg);
        end
        nIter = NaN; if isfield(d, 'nIter'); nIter = d.nIter; end
        fprintf('\n=== calibDetectionWindow (REPLAY of saved data, no hardware) ===\n');
    else
        % ================= MEASURE ========================================
        if isempty(opts.lastDelay)
            [opts.lastDelay, ok] = storedVal('lastDelay', opts.jsonLocation);
            if ~ok; opts.lastDelay = 0.2; end
        end
        fprintf('\n=== calibDetectionWindow ===\n');
        fprintf('window %g us, saturating MW burst %g us, gap %g us, %d x %d.  ZERO FIELD (%g MHz).\n', ...
                opts.window, opts.mwBurst, opts.lastDelay, opts.repeats, opts.averages, opts.frequency);

        c = AWGExpDigitizerInitDuration;
        if ~isempty(opts.jsonLocation)
            c.loadJsonParams('all', opts.jsonLocation);
        end
        % Setup 7 hardware, set explicitly: this class never overrides ExperimentAWG's
        % AWGorSRSswitch = 1 (a channel Setup 7 does not have), and nChannels defaults to 1.
        c.AWGorSRSswitch               = opts.AWGorSRSswitch;
        c.nChannels                    = opts.nChannels;
        c.GIMeas                       = opts.GIMeas;
        c.timeDelay                    = opts.timeDelay;
        c.amplitudeDiff                = opts.amplitudeDiff;
        c.repeats                      = opts.repeats;
        c.averages                     = opts.averages;
        c.detectionDuration            = opts.window;
        c.referenceDetectionDuration   = opts.window;
        c.laserInitializationDuration  = opts.window + opts.extraInit;
        c.frequency                    = opts.frequency;
        c.amplitude                    = opts.amplitude + zeros(1, opts.nChannels);
        c.piTime                       = opts.mwBurst;      % NOT a pi pulse -- see header
        c.lastDelay                    = opts.lastDelay;
        c.isTracking                   = 0;
        c.fixDelays                    = opts.fixDelays;
        % A full trace is thousands of samples, and the json's checkGlitchInData NaNs a whole
        % repeat when ANY sample is 4 sigma out -- on a 5000-sample trace that is most repeats.
        % Same fix as level 1.
        try; c.checkGlitchInData = false; catch; end

        c.run;

        if isempty(c.signal); error('The experiment produced no data.'); end
        nIter = c.currIter;
        S1 = reshape(c.signal(1, :, 1:nIter), [], nIter);
        S2 = reshape(c.signal(2, :, 1:nIter), [], nIter);
        t  = c.time(:);
        n  = min([numel(t), size(S1, 1)]);
        t  = t(1:n); S1 = S1(1:n, :); S2 = S2(1:n, :);
        bad = all(S1 == 0, 1) | all(~isfinite(S1), 1);
        S1(:, bad) = []; S2(:, bad) = [];
        s1 = nanMean2(S1); s2 = nanMean2(S2);
        Ck = 1 - S1 ./ S2;
    end

    dt = median(diff(t));
    Pl = opts.window + opts.extraInit;      % light between signal MW and reference window

    % ---- contrast and its errors -------------------------------------------
    if ~isempty(Ck) && size(Ck, 2) >= 3
        C  = nanMean2(Ck);
        dC = nanStd2(Ck) ./ sqrt(sum(isfinite(Ck), 2));
        noiseSrc = sprintf('standard error over %d averages of the per-average contrast', size(Ck, 2));
    else
        C  = 1 - s1 ./ s2;
        if ~isempty(opts.data) && isfield(opts.data, 'contrast'); C = opts.data.contrast(:); end
        dC = localScatter(C);
        noiseSrc = 'local point-to-point scatter (no per-average data)';
    end
    dC = max(dC, 0.25*median(dC));

    results = struct('time', t, 'signal', s1, 'reference', s2, 'contrast', C, ...
                     'contrastErr', dC, 'nIter', nIter, 'dt', dt, 'noiseSource', noiseSrc, ...
                     'opts', rmfieldSafe(opts, 'data'));
    if ~isempty(S1)
        % Kept so a replay can redo everything, including the noise exponent (~4 MB as single).
        results.signalPerAvg = single(S1); results.referencePerAvg = single(S2);
    elseif ~isempty(Ck)
        results.contrastPerAvg = single(Ck);
    end

    % ---- 1. where the light arrives ----------------------------------------
    tOn = lightOn(t, s2, dt);
    results.lightOn = tOn;

    % ---- 2. the decay fit ---------------------------------------------------
    use = t >= tOn & isfinite(C) & isfinite(dC);
    x = t(use) - tOn;  y = C(use);  sg = dC(use);
    fit = decayFit(x, y, sg, opts);
    fit.P = Pl;
    fit.aTrue = fit.a ./ (1 - exp(-Pl ./ fit.tau));         % undo the reference contamination
    fit.tauPol = sum(fit.aTrue .* fit.tau) / sum(fit.aTrue);
    fit.t1e = shapeCrossing(fit.aTrue, fit.tau, exp(-1));
    results.fit = fit;

    % ---- 3. noise growth with window length ----------------------------------
    beta = 0.5; betaSrc = 'white noise assumed';
    if ~isempty(opts.noiseExponent)
        beta = opts.noiseExponent; betaSrc = 'given';
    elseif ~isempty(S1) && size(S1, 2) >= 10
        [beta, betaSrc] = noiseExponent(t, tOn, S1, S2, dt, opts.window);
    end
    results.noiseExponent = beta;

    % ---- 4. window and init, together ----------------------------------------
    [stDet, haveDet]   = storedVal('detectionDuration', opts.jsonLocation);
    [stInit, haveInit] = storedVal('laserInitializationDuration', opts.jsonLocation);
    G = fomGrid(fit.aTrue, fit.tau, beta, opts.darkPerHalf, opts.window);
    pick = choose(G, opts.carryOverMax);
    trade = struct('cap', {}, 'det', {}, 'L', {}, 'fom', {}, 'carry', {});
    for cap = [0.10 0.05 0.02 0.01]
        pk = choose(G, cap);
        trade(end+1) = struct('cap', cap, 'det', pk.det, 'L', pk.L, 'fom', pk.fom, 'carry', pk.carry); %#ok<AGROW>
    end
    free = G.free;
    stored = struct('det', NaN, 'L', NaN, 'fom', NaN, 'carry', NaN);
    if haveDet && haveInit
        stored = evalPoint(fit.aTrue, fit.tau, beta, opts.darkPerHalf, stDet, stInit);
    end
    results.trade = trade;
    results.free  = free;
    results.stored = stored;

    % ---- 5. errors on the answers, by propagating the fit covariance -----------
    E = answerErrors(fit, Pl, beta, opts);
    fit.tauPolErr = E.tauPol;  fit.t1eErr = E.t1e;
    results.fit = fit;

    results.tauPol    = fit.tauPol;
    results.tauPolErr = E.tauPol;
    results.detectionDuration = pick.det;
    results.referenceDetectionDuration = pick.det;
    results.laserInitializationDuration = pick.L;
    results.carryOver = pick.carry;
    results.LErr = E.L;
    results.toWrite = struct();
    if isfinite(pick.det) && isfinite(pick.L)
        results.toWrite.detectionDuration          = pick.det;
        results.toWrite.referenceDetectionDuration = pick.det;
        results.toWrite.laserInitializationDuration = pick.L;
    end
    if isfinite(fit.tauPol)
        results.toWrite.tauPol = round(100*fit.tauPol)/100;
    end

    printReport(results, fit, opts, noiseSrc, beta, betaSrc, haveDet && haveInit);
    if opts.plot
        plotResult(t, C, tOn, fit, G, pick, stored, opts);
    end
end

% ======================================================================
%  LIGHT-ON
% ======================================================================

function tOn = lightOn(t, s2, dt)
% 90% point of the reference rise, the level-1 convention. 0 if no edge is visible.
    tOn = 0;
    base = s2(1);
    top  = median(s2(t > 0.5*t(end)));
    if abs(top - base) < 1e-9; return; end
    f = (s2 - base) / (top - base);
    k = find(f(1:end-1) < 0.9 & f(2:end) >= 0.9, 1);
    if isempty(k); return; end
    tOn = t(k) + (0.9 - f(k)) / (f(k+1) - f(k)) * (t(k+1) - t(k));
    tOn = max(tOn, 0);
end

% ======================================================================
%  THE DECAY FIT
% ======================================================================

function fit = decayFit(x, y, sg, opts)
    sw = 1 ./ sg(:);  x = x(:);  y = y(:);
    N  = numel(x);  xmax = x(end);  dtx = median(diff(x));
    lim = struct('tauMin', dtx/2, 'tauMax', 20*xmax, 'Pmin', max(1, 20*dtx), ...
                 'Pmax', 0.8*xmax, 'tdMin', 0.5, 'tdMax', 20*xmax);

    % rough 1/e time of the smoothed data, for starting values
    ys = movingMean(y, 11);
    k  = find(ys < ys(1)*exp(-1), 1);
    if isempty(k); t1 = xmax/3; else; t1 = max(x(k), 2*dtx); end

    ladder = struct('n', {}, 'osc', {}, 'chi2', {}, 'bic', {}, 'theta', {}, 'ok', {});
    bestByN = cell(1, opts.maxExp);
    for n = 1:opts.maxExp
        starts = expStarts(n, t1, bestByN, lim);
        [th, c2, ok] = multiStart(starts, x, y, sw, n, false, lim);
        bestByN{n} = th;
        ladder(end+1) = struct('n', n, 'osc', false, 'chi2', c2, ...
                               'bic', c2 + 2*n*log(N), 'theta', th, 'ok', ok); %#ok<AGROW>
        if opts.fitOscillation && N > 50
            % start the oscillation from the strongest periodic residual
            [~, cf] = vpChi(th, x, y, sw, n, false, lim);
            res = y - expBasis(x, exp(th(1:n))) * cf;
            P0  = periodogramStart(x, res, sw, lim);
            os  = {};
            for td0 = [0.5 2] * xmax
                for s = 1:numel(starts)
                    os{end+1} = [starts{s}(:)', log(P0), log(td0)]; %#ok<AGROW>
                end
                os{end+1} = [th(:)', log(P0), log(td0)]; %#ok<AGROW>
            end
            [th2, c22, ok2] = multiStart(os, x, y, sw, n, true, lim);
            ladder(end+1) = struct('n', n, 'osc', true, 'chi2', c22, ...
                                   'bic', c22 + (2*n + 4)*log(N), 'theta', th2, 'ok', ok2); %#ok<AGROW>
        end
    end
    okL = [ladder.ok] & isfinite([ladder.chi2]);
    if ~any(okL); error('calibDetectionWindow: no decay model could be fitted.'); end
    bics = [ladder.bic]; bics(~okL) = Inf;
    [~, ib] = min(bics);
    L = ladder(ib);

    [c2, cf] = vpChi(L.theta, x, y, sw, L.n, L.osc, lim);
    n   = L.n;
    tau = exp(L.theta(1:n)); tau = tau(:)';
    a   = cf(1:n)';
    [tau, o] = sort(tau); a = a(o);
    osc = [];
    if L.osc
        osc = struct('P', exp(L.theta(n+1)), 'td', exp(L.theta(n+2)), ...
                     'as', cf(n+1), 'ac', cf(n+2));
        osc.amp = hypot(osc.as, osc.ac);
    end
    k   = 2*n + 4*L.osc;
    dof = max(N - k, 1);
    chi2nu = c2 / dof;

    % covariance of [a, tau, (as, ac, P, td)], from the analytic Jacobian
    J = modelJacobian(x, a, tau, osc) .* sw;
    cov = pinv(J' * J) * max(chi2nu, 1);

    fit = struct('n', n, 'a', a, 'tau', tau, 'osc', osc, 'chi2nu', chi2nu, 'N', N, ...
                 'cov', cov, 'ladder', ladder, 'x', x);
    fit.model = @(xx) modelEval(xx, a, tau, osc);
    fit.decay = @(xx) modelEval(xx, a, tau, []);
    % chi2_nu of the single exponential, for the report
    oneExp = ladder([ladder.n] == 1 & ~[ladder.osc]);
    fit.chi2nuOneExp = oneExp.chi2 / max(N - 2, 1);
end

function S = expStarts(n, t1, bestByN, lim)
    S = {};
    switch n
        case 1
            S = {log(t1), log(3*t1), log(t1/3)};
        case 2
            S = {log([0.5 5]*t1), log([0.3 10]*t1), log([1 20]*t1)};
        otherwise
            S = {log([0.4 2 15]*t1), log([0.3 1.5 8]*t1), log([0.5 4 30]*t1)};
    end
    prev = bestByN{max(n-1, 1)};
    if n > 1 && ~isempty(prev)
        tp = exp(prev(:)');
        S{end+1} = log(sort([tp, 8*max(tp)]));
        S{end+1} = log(sort([tp, 0.25*min(tp)]));
    end
    for i = 1:numel(S)
        S{i} = min(max(S{i}, log(lim.tauMin*1.01)), log(lim.tauMax*0.99));
    end
end

function [thBest, cBest, ok] = multiStart(starts, x, y, sw, n, hasOsc, lim)
    thBest = starts{1}; cBest = Inf;
    o = optimset('Display', 'off', 'MaxFunEvals', 3000, 'MaxIter', 3000, ...
                 'TolX', 1e-7, 'TolFun', 1e-4);
    f = @(th) vpChi(th, x, y, sw, n, hasOsc, lim);
    for s = 1:numel(starts)
        th0 = starts{s};
        if ~isfinite(f(th0)); continue; end
        th = fminsearch(f, th0, o);
        th = fminsearch(f, th, o);              % one restart: the simplex often stalls once
        cv = f(th);
        if cv < cBest; cBest = cv; thBest = th; end
    end
    ok = isfinite(cBest);
end

function [c2, cf] = vpChi(th, x, y, sw, n, hasOsc, lim)
% Variable projection: amplitudes solved exactly for the given time constants.
    cf = [];
    tau = exp(th(1:n));
    if any(tau < lim.tauMin | tau > lim.tauMax); c2 = Inf; return; end
    B = expBasis(x, tau);
    if hasOsc
        P = exp(th(n+1)); td = exp(th(n+2));
        if P < lim.Pmin || P > lim.Pmax || td < max(lim.tdMin, 0.5*P) || td > lim.tdMax
            c2 = Inf; return
        end
        e = exp(-x/td);
        B = [B, e.*sin(2*pi*x/P), e.*cos(2*pi*x/P)];
    end
    Bw = B .* sw;  yw = y .* sw;
    cf = Bw \ yw;
    if any(~isfinite(cf)) || any(cf(1:n) <= 0); c2 = Inf; return; end
    r  = yw - Bw*cf;
    c2 = r' * r;
end

function B = expBasis(x, tau)
    B = exp(-x(:) ./ tau(:)');
end

function P0 = periodogramStart(x, res, sw, lim)
% The period whose damped sin/cos pair removes the most chi-square from the residuals.
    Ps = exp(linspace(log(lim.Pmin), log(lim.Pmax), 200));
    e  = exp(-x / (0.5*x(end)));
    rw = res .* sw;  base = rw' * rw;  gain = zeros(size(Ps));
    for i = 1:numel(Ps)
        B = [e.*sin(2*pi*x/Ps(i)), e.*cos(2*pi*x/Ps(i))] .* sw;
        cf = B \ rw;  r = rw - B*cf;
        gain(i) = base - r'*r;
    end
    [~, k] = max(gain);
    P0 = Ps(k);
end

function v = modelEval(xx, a, tau, osc)
    xx = xx(:);
    v = expBasis(xx, tau) * a(:);
    if ~isempty(osc)
        e = exp(-xx/osc.td);
        v = v + e .* (osc.as*sin(2*pi*xx/osc.P) + osc.ac*cos(2*pi*xx/osc.P));
    end
end

function J = modelJacobian(x, a, tau, osc)
    E  = expBasis(x, tau);
    J  = [E, E .* (x ./ tau.^2) .* a];
    if ~isempty(osc)
        e = exp(-x/osc.td); ph = 2*pi*x/osc.P;
        sn = sin(ph); cs = cos(ph);
        J = [J, e.*sn, e.*cs, ...
             e .* (osc.as*cs - osc.ac*sn) .* (-ph/osc.P), ...
             e .* (osc.as*sn + osc.ac*cs) .* (x/osc.td^2)];
    end
end

function tc = shapeCrossing(a, tau, level)
% Time at which c(t)/c(0) of the sum of exponentials falls to 'level' (monotone, so bisect).
    c = @(s) sum(a .* exp(-s ./ tau)) / sum(a);
    lo = 0; hi = max(tau) * 50;
    if c(hi) > level; tc = NaN; return; end
    for i = 1:80
        mid = 0.5*(lo + hi);
        if c(mid) > level; lo = mid; else; hi = mid; end
    end
    tc = 0.5*(lo + hi);
end

% ======================================================================
%  WINDOW AND INIT
% ======================================================================

function G = fomGrid(a, tau, beta, D, window)
% FOM and carry-over over detection window d and total light per half-repeat L (d <= L).
    G.d = (0.05 : 0.05 : min(window, 20))';
    G.L = 0.5 : 0.5 : min(max(200, 4*max(tau)), 2000);
    a = a(:)'; tau = tau(:)';
    open = 1 - exp(-G.d ./ tau);                         % nd x nc
    Isig = (open .* (a .* tau)) * ones(numel(tau), 1);   % nd x 1
    Iat  = @(L0) open * ((a .* tau)' .* exp(-L0(:)' ./ tau'));   % nd x nL
    Iref  = Iat(G.L);
    Inext = Iat(2*G.L);
    num   = Isig - Iref;
    F     = num ./ (G.d.^beta .* sqrt(D + G.L));
    X     = Inext ./ max(num, eps);
    valid = G.d <= G.L & num > 0;
    F(~valid) = NaN;  X(~valid) = NaN;
    G.F = F; G.X = X;
    Fz = F; Fz(isnan(Fz)) = -Inf;
    [fm, k] = max(Fz(:)); [i, j] = ind2sub(size(F), k);
    G.fomMax = fm;
    G.free = struct('det', G.d(i), 'L', G.L(j), 'fom', fm, 'carry', X(i, j));
end

function pk = choose(G, cap)
% Best FOM over (d, L) with carry-over under the cap. Joint, because at a given L a SHORTER
% window also carries less into the next point -- picking the best window per L first and the
% L afterwards lands on a worse point.
    Fz = G.F; Fz(~(G.X <= cap)) = -Inf;
    [fm, k] = max(Fz(:));
    if ~isfinite(fm)
        pk = struct('det', NaN, 'L', NaN, 'fom', NaN, 'carry', NaN); return
    end
    [i, j] = ind2sub(size(Fz), k);
    pk = struct('det', G.d(i), 'L', G.L(j), 'fom', fm, 'carry', G.X(i, j));
end

function p = evalPoint(a, tau, beta, D, d, L)
    I = @(s0, s1) sum(a .* tau .* (exp(-s0 ./ tau) - exp(-s1 ./ tau)));
    num = I(0, d) - I(L, L + d);
    p = struct('det', d, 'L', L, 'fom', num / (d^beta * sqrt(D + L)), ...
               'carry', I(2*L, 2*L + d) / max(num, eps));
end

function [beta, src] = noiseExponent(t, tOn, S1, S2, dt, window)
% How the noise of a window-integrated contrast grows with the window, from the scatter of
% the per-average values. White noise gives beta = 0.5; drift between averages pushes it up.
    i0 = find(t >= tOn, 1);
    ds = exp(linspace(log(max(1, 10*dt)), log(0.8*window), 20));
    sd = nan(size(ds));
    for i = 1:numel(ds)
        idx = i0 : min(i0 + round(ds(i)/dt) - 1, size(S1, 1));
        w   = 1 - sum(S1(idx, :), 1) ./ sum(S2(idx, :), 1);
        sd(i) = std(w(isfinite(w)));
    end
    p = polyfit(log(ds), log(sd), 1);
    beta = min(max(1 + p(1), 0.5), 1);
    src  = sprintf('measured from %d averages (window-contrast noise ~ d^%.2f)', size(S1, 2), p(1));
end

function E = answerErrors(fit, Pl, beta, opts)
% Propagate the fit covariance into tauPol, the 1/e time and the chosen L, by sampling.
% Quasi-random (Halton) draws, so the global random-number state is left alone -- the walk
% modes rely on rng('shuffle') for their visit order.
    E = struct('tauPol', NaN, 't1e', NaN, 'L', NaN);
    n = fit.n; k = size(fit.cov, 1);
    [V, Dg] = eig((fit.cov + fit.cov')/2);
    Dg = max(diag(Dg), 0);
    q0 = [fit.a(:); fit.tau(:)];
    if ~isempty(fit.osc); q0 = [q0; fit.osc.as; fit.osc.ac; fit.osc.P; fit.osc.td]; end
    if numel(q0) ~= k; return; end
    M = 200;
    Z = haltonNormal(M, k);
    tp = nan(M, 1); te = tp; Ls = tp;
    for m = 1:M
        q = q0 + V * (sqrt(Dg) .* Z(m, :)');
        a = q(1:n)'; tau = q(n+1:2*n)';
        if any(a <= 0) || any(tau <= 0); continue; end
        at = a ./ (1 - exp(-Pl ./ tau));
        tp(m) = sum(at .* tau) / sum(at);
        te(m) = shapeCrossing(at, tau, exp(-1));
        G = fomGrid(at, tau, beta, opts.darkPerHalf, opts.window);
        pk = choose(G, opts.carryOverMax);
        Ls(m) = pk.L;
    end
    E.tauPol = std(tp(isfinite(tp)));
    E.t1e    = std(te(isfinite(te)));
    E.L      = std(Ls(isfinite(Ls)));
end

function Z = haltonNormal(M, k)
    pr = [2 3 5 7 11 13 17 19 23 29 31 37 41 43 47];
    Z = zeros(M, k);
    for j = 1:k
        b = pr(j);
        for i = 1:M
            f = 1; r = 0; ii = i;
            while ii > 0
                f = f / b; r = r + f * mod(ii, b); ii = floor(ii / b);
            end
            Z(i, j) = -sqrt(2) * erfcinv(2*r);
        end
    end
end

% ======================================================================
%  REPORT
% ======================================================================

function printReport(R, fit, opts, noiseSrc, beta, betaSrc, haveStored)
    fprintf('%d samples at %g ns, %s.\n', numel(R.time), 1e3*R.dt, noiseSrc);
    if R.lightOn > 5*R.dt
        fprintf(['light on: 90%% at %.3f us, %.0f samples after the window opens; level 1 aims for 0.\n' ...
                 '   Those samples carry little light. Worth a look at the greenLaser delay.\n'], ...
                R.lightOn, R.lightOn/R.dt);
    else
        fprintf('light on: 90%% at %.3f us, consistent with level 1.\n', R.lightOn);
    end

    % ---- decay ----
    fprintf('\nPOLARISATION DECAY  (fit from light-on)\n');
    word = sprintf('%d exponential%s', fit.n, plural(fit.n));
    if ~isempty(fit.osc); word = [word ' + oscillation']; end
    fprintf('   model  %-30s chi2_nu %.2f   (one exponential: %.3g)\n', word, fit.chi2nu, fit.chi2nuOneExp);
    sh = fit.aTrue / sum(fit.aTrue);
    nq = numel(fit.a);
    tErr = sqrt(max(diag(fit.cov(nq+1:2*nq, nq+1:2*nq)), 0))';
    for i = 1:fit.n
        fprintf('   component %d   tau %8.3f +/- %-7.3f us   share %4.0f%%\n', i, fit.tau(i), tErr(i), 100*sh(i));
    end
    fprintf('   tauPol  %.2f +/- %.2f us (stat.)   mean re-polarisation time: area under c/c0, = tau for one exponential\n', ...
            fit.tauPol, fit.tauPolErr);
    fprintf('   1/e     %.2f us\n', fit.t1e);
    if fit.n > 1 && max(fit.tau) > 0.5*opts.window
        fprintf(['   NOTE the slowest component (%.0f us) is long against the %g us window, so it is\n' ...
                 '   only loosely pinned down. A longer window would tighten it.\n'], max(fit.tau), opts.window);
    end
    if ~isempty(fit.osc)
        fprintf(['   oscillation  %.1f us period (%.0f kHz), %.3f%% amplitude, decays in %.0f us.\n' ...
                 '                Excluded from everything below. Signal window only, so MW-triggered; to\n' ...
                 '                confirm, re-run with opts.frequency off resonance -- it should stay.\n'], ...
                fit.osc.P, 1e3/fit.osc.P, 100*fit.osc.amp, fit.osc.td);
    end

    % ---- window and init ----
    fprintf('\nWINDOW AND INIT  (best SNR per unit time with carry-over into the next point <= %.0f%%)\n', ...
            100*opts.carryOverMax);
    fprintf('   assumes %.1f us of dark per half-repeat (opts.darkPerHalf); noise exponent %.2f, %s\n', ...
            opts.darkPerHalf, beta, betaSrc);
    ref = R.free.fom; refWord = 'unconstrained optimum';
    if haveStored && isfinite(R.stored.fom); ref = R.stored.fom; refWord = 'stored values'; end
    fprintf('\n   %-10s %-20s %-30s %s\n', 'carry-over', 'detectionDuration', ...
            'laserInitializationDuration', sprintf('SNR/time vs %s', refWord));
    for i = 1:numel(R.trade)
        T = R.trade(i);
        mark = '  '; if abs(T.cap - opts.carryOverMax) < 1e-12; mark = '->'; end
        fprintf('%s  <= %-5s  %8.2f us %20.1f us %24s\n', mark, sprintf('%.0f%%', 100*T.cap), ...
                T.det, T.L, pct(T.fom, ref));
    end
    if haveStored && isfinite(R.stored.fom)
        fprintf('   stored     %8.2f us %20.1f us %24s   (carry-over %.1f%%)\n', ...
                R.stored.det, R.stored.L, pct(R.stored.fom, ref), 100*R.stored.carry);
    end
    fprintf(['   classic mode only: %.2f / %.1f us gives %s, but carries %.0f%% of each\n' ...
             '   point into the next -- fine when every repeat is the same point, not in walk mode.\n'], ...
            R.free.det, R.free.L, pct(R.free.fom, ref), 100*R.free.carry);
    if isfinite(R.LErr)
        fprintf('   (the recommended init is known to about +/- %.1f us from the fit)\n', R.LErr);
    end
    if 2*R.laserInitializationDuration > opts.window
        fprintf(['   carry-over at L = %.1f us is read at 2L = %.0f us, past the %g us window, so it\n' ...
                 '   leans on the fitted slow component.\n'], R.laserInitializationDuration, ...
                2*R.laserInitializationDuration, opts.window);
    end

    % ---- to write ----
    fprintf('\nTO WRITE  (level2.m cell 2.3 does this with expPropsWrite)\n');
    w = R.toWrite; f = fieldnames(w);
    for i = 1:numel(f)
        fprintf('   %-28s %g\n', f{i}, w.(f{i}));
    end
    fprintf('\n');
end

function s = pct(v, ref)
    if ~isfinite(v) || ~isfinite(ref) || ref == 0; s = '--'; return; end
    s = sprintf('%+.0f%%', 100*(v/ref - 1));
end

function s = plural(n)
    if n == 1; s = ''; else; s = 's'; end
end

% ======================================================================
%  PLOT
% ======================================================================

function plotResult(t, C, tOn, fit, G, pick, stored, opts)
% Three axes, each titled with the number it decides.
    COL.data = [0.62 0.72 0.84]; COL.fit = [0.10 0.10 0.10]; COL.dec = [0.85 0.33 0.10];
    COL.cmp  = [0.45 0.45 0.45]; COL.rec = [0.00 0.45 0.35]; COL.old = [0.85 0.33 0.10];
    COL.fom  = [0.49 0.18 0.56];

    figure('Color', 'w', 'Position', [80 80 1250 560]);

    % ---- decay ----
    ax1 = subplot(2, 2, [1 3]); hold(ax1, 'on'); tidy(ax1);
    k  = t >= tOn;  xx = t(k) - tOn;  Cx = C(k);
    pos = Cx > 0;
    hD = plot(ax1, xx(pos), 100*Cx(pos), '.', 'Color', COL.data, 'MarkerSize', 4);
    xf = linspace(0, xx(end), 1500)';
    hF = plot(ax1, xf, 100*fit.model(xf), '-', 'Color', COL.fit, 'LineWidth', 1.0);
    hE = plot(ax1, xf, 100*fit.decay(xf), '-', 'Color', COL.dec, 'LineWidth', 1.8);
    hc = [];
    for i = 1:fit.n
        hc = plot(ax1, xf, 100*fit.a(i)*exp(-xf/fit.tau(i)), '--', 'Color', COL.cmp, 'LineWidth', 0.9);
    end
    tk = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20];
    set(ax1, 'YScale', 'log', 'YTick', tk, 'YTickLabel', arrayfun(@(v) sprintf('%g', v), tk, ...
        'UniformOutput', false));
    yl = [max(min(100*C(k & C > 0)), 1e-3*100*max(C)) , 1.3*100*max(C)];
    ylim(ax1, yl); xlim(ax1, [0 xx(end)]);
    xlabel(ax1, 'time after light-on (us)'); ylabel(ax1, 'contrast (%)');
    lab = {'measured', 'fit', 'decay part (used below)'}; hh = [hD hF hE];
    if ~isempty(hc); hh(end+1) = hc; lab{end+1} = 'components'; end
    legend(ax1, hh, lab, 'Location', 'northeast', 'Box', 'off');
    ttl = sprintf('tauPol = %.2f us  (%d exp', fit.tauPol, fit.n);
    if ~isempty(fit.osc); ttl = [ttl sprintf(' + %.0f kHz artefact', 1e3/fit.osc.P)]; end
    title(ax1, [ttl sprintf(', \\chi^2_\\nu = %.1f)', fit.chi2nu)], 'FontWeight', 'bold');

    % ---- the map that decides both numbers ----
    ax2 = subplot(2, 2, [2 4]); hold(ax2, 'on'); tidy(ax2);
    Lmax = min(G.L(end), 2*max([pick.L, stored.L, 30]));
    dmax = min(G.d(end), 2.5*max([pick.det, stored.det, G.free.det, 2]));
    jL = G.L <= Lmax;  iD = G.d <= dmax;
    ref = pick.fom; if ~isfinite(ref); ref = G.fomMax; end
    Z = 100 * G.F(iD, jL) / ref;
    lev = [40 55 70 85 100 115 130 150];
    [~, hcf] = contourf(ax2, G.L(jL), G.d(iD), min(max(Z, 40), 160), lev);
    set(hcf, 'LineColor', 'none');
    colormap(ax2, purpleMap(64)); caxis(ax2, [40 160]);
    cb = colorbar(ax2); ylabel(cb, 'SNR per unit time, % of the chosen point');
    Xc = 100 * G.X(iD, jL);
    [c2, h2] = contour(ax2, G.L(jL), G.d(iD), Xc, [2 5 10 20], 'LineColor', [0 0 0], 'LineWidth', 1.1);
    try; clabel(c2, h2, 'FontSize', 8, 'Color', [0 0 0], 'LabelSpacing', 400); catch; end
    hX = plot(ax2, NaN, NaN, 'k-', 'LineWidth', 1.1);
    hs = hX; ls = {'carry-over into the next point, %'};
    if isfinite(pick.L)
        hR = plot(ax2, pick.L, pick.det, 'p', 'MarkerSize', 16, 'MarkerFaceColor', COL.rec, ...
                  'MarkerEdgeColor', 'w');
        hs(end+1) = hR; ls{end+1} = sprintf('USE  %.1f / %.2f us', pick.L, pick.det);
    end
    if isfinite(stored.L) && stored.L <= Lmax && stored.det <= dmax
        hO = plot(ax2, stored.L, stored.det, 'o', 'MarkerSize', 9, 'LineWidth', 1.8, ...
                  'MarkerEdgeColor', COL.old);
        hs(end+1) = hO; ls{end+1} = sprintf('stored  %.1f / %.2f us', stored.L, stored.det);
    end
    legend(ax2, hs, ls, 'Location', 'northwest', 'Box', 'on');
    xlim(ax2, [0 Lmax]); ylim(ax2, [0 dmax]);
    xlabel(ax2, 'laserInitializationDuration (us)');
    ylabel(ax2, 'detectionDuration (us)');
    title(ax2, sprintf('laserInitializationDuration %.1f us,  detectionDuration %.2f us  (carry-over <= %.0f%%)', ...
          pick.L, pick.det, 100*opts.carryOverMax), 'FontWeight', 'bold');
end

function m = purpleMap(n)
    a = [0.97 0.96 0.99]; b = [0.42 0.16 0.53];
    t = linspace(0, 1, n)';
    m = a + (b - a) .* t;
end

function tidy(ax)
    set(ax, 'Box', 'off', 'TickDir', 'out', 'Layer', 'top', ...
            'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.12, 'FontSize', 9);
end

% ======================================================================
%  HELPERS
% ======================================================================

function s = localScatter(y)
% Local point-to-point noise from second differences, smoothed over ~100 samples. Robust to
% the decay itself, which second differences remove.
    y = y(:);
    d = (y(3:end) - 2*y(2:end-1) + y(1:end-2)) / sqrt(6);
    d = [d(1); d; d(end)];
    s = sqrt(movingMean(d.^2, 101));
end

function m = movingMean(y, w)
    y = y(:); n = numel(y); h = floor(w/2);
    cs = [0; cumsum(y)];
    lo = max((1:n)' - h, 1); hi = min((1:n)' + h, n);
    m = (cs(hi+1) - cs(lo)) ./ (hi - lo + 1);
end

function m = nanMean2(X)
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
