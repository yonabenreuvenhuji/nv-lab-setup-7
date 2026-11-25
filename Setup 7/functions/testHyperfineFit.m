function allOk = testHyperfineFit(realFile)
%TESTHYPERFINEFIT Offline grader for esrHyperfineFit. No hardware.
%
%   testHyperfineFit                       % synthetic cases only
%   testHyperfineFit('G:\...\_AutoSave\AWG ESR Pulsed_20261110_172200.mat')   % + a real save
%
% Every synthetic spectrum is generated HERE with its own formulas (the pulse shape by
% brute-force numerical convolution, not the fitter's FFT table), so a shared bug in the
% model cannot grade itself. Prints one line per check and ALL CHECKS PASSED at the end.

    if ~exist('realFile', 'var'); realFile = ''; end
    rng(20261116);
    nFail = 0; nCheck = 0;
    q = struct('plotIt', false, 'verbose', false);
    f = (2804.84:0.5967:2828.71)';              % the 2026-11-10 grid
    sig = 2e-4;                                 % noise, = the error bars given

    % ---- 1. 14N triplet, Lorentzian, clean ---------------------------------------
    y = model(f, 2816.33, 2.16, [1 1 1]*1.8e-3, 1.0, 'lorentz', []) + sig*randn(size(f));
    o = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('14N: 3 dips', o.nDips == 3);
    check('14N: isotope 14N', strcmp(o.isotope, '14N'));
    check('14N: confident', any(strcmp(o.confidence, {'strong', 'decisive'})));
    check('14N: A within 3 sigma', abs(o.splitting - 2.16) < 3*o.splittingErr);
    check('14N: f0 within 3 sigma', abs(o.f0 - 2816.33) < 3*o.f0Err);
    check('14N: modelAt matches model', max(abs(o.modelAt(o.f) - o.model)) < 1e-12);

    % ---- 2. 15N doublet -------------------------------------------------------------
    y = model(f, 2815.9, 3.03, [1 1]*1.8e-3, 1.0, 'gauss', []) + sig*randn(size(f));
    o = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('15N: 2 dips', o.nDips == 2);
    check('15N: isotope 15N', strcmp(o.isotope, '15N'));
    check('15N: A within 3 sigma', abs(o.splitting - 3.03) < 3*o.splittingErr);

    % ---- 3. one broad line (hyperfine unresolved) -------------------------------------
    y = model(f, 2817.0, 0, 3e-3, 6.0, 'lorentz', []) + sig*randn(size(f));
    o = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('broad: 1 dip', o.nDips == 1);
    check('broad: unresolved', strcmp(o.isotope, 'unresolved'));
    check('broad: f0 within 3 sigma', abs(o.f0 - 2817.0) < 3*o.f0Err);

    % ---- 4. nothing there -----------------------------------------------------------
    y = 1 + 1e-6*(f - 2816) + sig*randn(size(f));
    o = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('flat: no line', o.nDips == 0 && strcmp(o.isotope, 'none'));

    % ---- 5. pulsed-ESR triplet like the real run (tau 0.83 us, inhom. 1.0 MHz) --------
    tau = 0.828;
    y = model(f, 2816.33, 2.16, [1 1 1]*2.4e-3, 1.0, 'pulse', tau) + sig*randn(size(f));
    o = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y, 'e', sig*ones(size(f)), 'tau', tau)));
    check('pulse: 3 dips 14N', o.nDips == 3 && strcmp(o.isotope, '14N'));
    check('pulse: A within 3 sigma', abs(o.splitting - 2.16) < 3*o.splittingErr);
    o2 = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y, 'e', sig*ones(size(f)), ...
                                             'tau', tau, 'shape', 'pulse', 'nLines', 3)));
    check('pulse: inhomogeneous width within 3 sigma', ...
          abs(o2.widthParam - 1.0) < 3*o2.widthParamErr);
    yNoiseless = model(f, 2816.33, 2.16, [1 1 1]*2.4e-3, 1.0, 'pulse', tau);
    o3 = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', yNoiseless, 'e', sig*ones(size(f)), ...
                                             'tau', tau, 'shape', 'pulse', 'nLines', 3)));
    check('pulse: noiseless shape reproduced to 1e-6', max(abs(o3.residual)) < 1e-6);

    % ---- 6. polarised nitrogen: unequal depths -----------------------------------------
    y = model(f, 2816.33, 2.16, [2.6 1.6 0.7]*1e-3, 0.9, 'lorentz', []) + 1e-4*randn(size(f));
    o = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y, 'e', 1e-4*ones(size(f)))));
    check('polarised: 3 dips', o.nDips == 3);
    check('polarised: depths freed', strcmp(o.depthTest.used, 'free'));
    d = [o.lines.depth];
    check('polarised: depths ordered as generated', d(1) > d(2) && d(2) > d(3));

    % ---- 7. two lines 4.5 MHz apart: not a nitrogen pattern ----------------------------
    y = model(f, 2816.0, 4.5, [1 1]*1.8e-3, 1.0, 'lorentz', []) + sig*randn(size(f));
    o = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('4.5 MHz pair: 2 dips, unknown', o.nDips == 2 && strcmp(o.isotope, 'unknown'));

    % ---- 8. no error bars, sloped baseline ------------------------------------------------
    y = model(f, 2816.33, 2.16, [1 1 1]*1.8e-3, 1.0, 'gauss', []) + 2e-5*(f - 2816) ...
        + sig*randn(size(f));
    o  = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y)));
    oE = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('no e: 3 dips 14N', o.nDips == 3 && strcmp(o.isotope, '14N'));
    check('no e: same A error as with the true error bars (within 10%)', ...
          abs(o.splittingErr/oE.splittingErr - 1) < 0.1);

    % ---- 9. forced count -------------------------------------------------------------------
    o = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y, 'nLines', 2)));
    check('forced: nDips = 2, confidence forced', o.nDips == 2 && strcmp(o.confidence, 'forced'));

    % ---- 10. low SNR: never confidently wrong --------------------------------------------
    nWrong = 0; nConf = 0;
    for k = 1:20
        y = model(f, 2816.33, 2.16, [1 1 1]*4e-4, 1.0, 'lorentz', []) + sig*randn(size(f));
        o = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y, 'e', sig*ones(size(f)), ...
                                               'shape', 'lorentz')));
        if any(strcmp(o.confidence, {'strong', 'decisive'}))
            nConf = nConf + 1;
            if o.nDips ~= 3; nWrong = nWrong + 1; end
        end
    end
    fprintf('       low SNR: %d/20 confident, %d of them wrong\n', nConf, nWrong);
    check('low SNR: at most 1 confident wrong answer in 20', nWrong <= 1);

    % ---- 11. the error bars are right: pulls of A and f0 -------------------------------------
    nR = 60; pa = nan(nR, 1); pf = nan(nR, 1); cls = 0;
    for k = 1:nR
        y = model(f, 2816.33, 2.16, [1 1 1]*1.5e-3, 1.0, 'lorentz', []) + sig*randn(size(f));
        o = esrHyperfineFit(mergeOpts(q, struct('f', f, 'y', y, 'e', sig*ones(size(f)), ...
                                               'shape', 'lorentz')));
        cls = cls + (o.nDips == 3);
        pa(k) = (o.splitting - 2.16) / o.splittingErr;
        pf(k) = (o.f0 - 2816.33) / o.f0Err;
    end
    fprintf('       pulls over %d: A mean %+.2f std %.2f;  f0 mean %+.2f std %.2f;  %d/%d called 3\n', ...
            nR, mean(pa), std(pa), mean(pf), std(pf), cls, nR);
    check('pulls: std of A pull in [0.75 1.3]', std(pa) > 0.75 && std(pa) < 1.3);
    check('pulls: std of f0 pull in [0.75 1.3]', std(pf) > 0.75 && std(pf) < 1.3);
    check('pulls: |mean| < 0.4', abs(mean(pa)) < 0.4 && abs(mean(pf)) < 0.4);
    check('pulls: all called 3 dips', cls == nR);

    % ---- 12. a real save, if given -------------------------------------------------------------
    if ~isempty(realFile)
        o = esrHyperfineFit(mergeOpts(q, struct('file', realFile)));
        fprintf('       real file: %s\n', o.verdict);
        check('real file: 3 dips 14N', o.nDips == 3 && strcmp(o.isotope, '14N'));
        check('real file: tau read from the save', ~isempty(o.tau) && o.tau > 0);
    end

    allOk = nFail == 0;
    if allOk
        fprintf('\nALL CHECKS PASSED (%d)\n', nCheck);
    else
        fprintf('\n%d of %d CHECKS FAILED\n', nFail, nCheck);
    end

    function check(name, ok)
        nCheck = nCheck + 1;
        if ok; tag = 'pass'; else; tag = 'FAIL'; nFail = nFail + 1; end
        fprintf('  %s  %s\n', tag, name);
    end
end

% ======================================================================

function y = model(f, f0, A, C, w, shp, tau)
% The generator, independent of the fitter: y = 1 - SUM C_k * S(f - f_k).
    n = numel(C);
    switch n
        case 1; m = 0;
        case 2; m = [-0.5 0.5];
        case 3; m = [-1 0 1];
    end
    y = ones(size(f));
    for k = 1:n
        u = f - (f0 + m(k)*A);
        switch shp
            case 'lorentz'; s = 1 ./ (1 + (2*u/w).^2);
            case 'gauss';   s = exp(-4*log(2)*(u/w).^2);
            case 'pulse';   s = pulseBrute(u, w, tau);
        end
        y = y - C(k)*s;
    end
end

function s = pulseBrute(u, w, tau)
% pi-pulse response convolved with a Gaussian of FWHM w, by direct trapezoid integration.
    Om = 1/(2*tau);
    sg = w/(2*sqrt(2*log(2)));
    v = linspace(-7*sg, 7*sg, 6001);
    g = exp(-0.5*(v/sg).^2); g = g / trapz(v, g);
    s = zeros(size(u));
    for i = 1:numel(u)
        d = u(i) - v;
        W2 = Om^2 + d.^2;
        s(i) = trapz(v, g .* (Om^2 ./ W2 .* sin(pi*tau*sqrt(W2)).^2));
    end
end

function o = mergeOpts(a, b)
    o = a; fn = fieldnames(b);
    for i = 1:numel(fn); o.(fn{i}) = b.(fn{i}); end
end
