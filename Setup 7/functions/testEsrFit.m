function allOk = testEsrFit(realFile)
%TESTESRFIT Offline grader for esrFitRobust: lines, resonance, field, hyperfine. No hardware.
%
%   testEsrFit                       % synthetic cases only
%   testEsrFit('G:\...\_AutoSave\AWG ESR Pulsed_20261110_172200.mat')   % + a real pulsed save
%
% Every synthetic spectrum is generated HERE with its own formulas (the pulse shape by
% brute-force numerical convolution, not the fitter's FFT table), so a shared bug in the
% model cannot grade itself. The field cases use esrSynthSpectrum (exact Hamiltonian). Prints
% one line per check and ALL CHECKS PASSED at the end. testFieldPipeline grades the
% calibFieldSurvey chain on top of this.
%
%   1  hyperfine       14N / 15N / unresolved / none, pulse shape, polarised depths, a
%                      non-nitrogen pair, no error bars, forced count, low SNR, pulls
%   2  hyperfine gate  'auto' tests any line on a fine grid -- a resolved pulsed triplet whose
%                      one-line FWHM is the envelope included -- and leaves a 1 MHz-step CW
%                      line alone; a broad line on a fine grid comes back not resolved
%   3  resonance       at zero field just outside and inside the sweep, nLines = 1 counting
%                      the NV line only, no false resonance on clean data, off, and a hint
%   4  field           45 G with and without a resonance; one pair; zero field; a window
%   5  interface       unknown options, the positional form, struct input, modelAt

    if ~exist('realFile', 'var'); realFile = ''; end
    rng(20261126, 'twister');     % named: a bare rng(seed) errors in a legacy-mode session
    nFail = 0; nCheck = 0;
    q = struct('plotIt', false, 'verbose', false);
    lor = @(f, c, w) 1 ./ (1 + (2*(f - c)/w).^2);

    % =====================================================================
    fprintf('\n1. HYPERFINE\n');
    f = (2804.84:0.5967:2828.71)';              % the 2026-11-10 pulsed grid
    sig = 2e-4;                                 % noise, = the error bars given
    hq = mergeOpts(q, struct('hyperfine', true));

    y = model(f, 2816.33, 2.16, [1 1 1]*1.8e-3, 1.0, 'lorentz', []) + sig*randn(size(f));
    h = hfOf(mergeOpts(hq, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('14N: 3 dips', h.nDips == 3);
    check('14N: isotope 14N, resolved', strcmp(h.isotope, '14N') && h.resolved);
    check('14N: confident', any(strcmp(h.confidence, {'strong', 'decisive'})));
    check('14N: A within 3 sigma', abs(h.splitting - 2.16) < 3*h.splittingErr);
    check('14N: f0 within 3 sigma', abs(h.f0 - 2816.33) < 3*h.f0Err);
    check('14N: modelAt matches model', max(abs(h.modelAt(h.f) - h.model)) < 1e-12);

    y = model(f, 2815.9, 3.03, [1 1]*1.8e-3, 1.0, 'gauss', []) + sig*randn(size(f));
    h = hfOf(mergeOpts(hq, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('15N: 2 dips, isotope 15N', h.nDips == 2 && strcmp(h.isotope, '15N'));
    check('15N: A within 3 sigma', abs(h.splitting - 3.03) < 3*h.splittingErr);

    y = model(f, 2817.0, 0, 3e-3, 6.0, 'lorentz', []) + sig*randn(size(f));
    h = hfOf(mergeOpts(hq, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('broad: 1 dip, unresolved', h.nDips == 1 && strcmp(h.isotope, 'unresolved') && ~h.resolved);
    check('broad: f0 within 3 sigma', abs(h.f0 - 2817.0) < 3*h.f0Err);

    y = 1 + 1e-6*(f - 2816) + sig*randn(size(f));
    o = esrFitRobust(mergeOpts(hq, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('flat: no line reported', ~o.ok || o.hyperfine.nDips == 0);

    tau = 0.828;
    y = model(f, 2816.33, 2.16, [1 1 1]*2.4e-3, 1.0, 'pulse', tau) + sig*randn(size(f));
    h = hfOf(mergeOpts(hq, struct('f', f, 'y', y, 'e', sig*ones(size(f)), 'tau', tau)));
    check('pulse: 3 dips 14N', h.nDips == 3 && strcmp(h.isotope, '14N'));
    check('pulse: A within 3 sigma', abs(h.splitting - 2.16) < 3*h.splittingErr);
    h2 = hfOf(mergeOpts(hq, struct('f', f, 'y', y, 'e', sig*ones(size(f)), 'tau', tau, ...
                                   'hfShape', 'pulse', 'hfCount', 3)));
    check('pulse: inhomogeneous width within 3 sigma', abs(h2.widthParam - 1.0) < 3*h2.widthParamErr);
    yNoiseless = model(f, 2816.33, 2.16, [1 1 1]*2.4e-3, 1.0, 'pulse', tau);
    h3 = hfOf(mergeOpts(hq, struct('f', f, 'y', yNoiseless, 'e', sig*ones(size(f)), 'tau', tau, ...
                                   'hfShape', 'pulse', 'hfCount', 3)));
    check('pulse: noiseless shape reproduced to 1e-6', max(abs(h3.residual)) < 1e-6);

    y = model(f, 2816.33, 2.16, [2.6 1.6 0.7]*1e-3, 0.9, 'lorentz', []) + 1e-4*randn(size(f));
    h = hfOf(mergeOpts(hq, struct('f', f, 'y', y, 'e', 1e-4*ones(size(f)))));
    check('polarised: 3 dips, depths freed', h.nDips == 3 && strcmp(h.depthTest.used, 'free'));
    d = [h.lines.depth];
    check('polarised: depths ordered as generated', numel(d) == 3 && d(1) > d(2) && d(2) > d(3));

    y = model(f, 2816.0, 4.5, [1 1]*1.8e-3, 1.0, 'lorentz', []) + sig*randn(size(f));
    h = hfOf(mergeOpts(hq, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('4.5 MHz pair: 2 dips, unknown, not resolved', h.nDips == 2 && ...
          strcmp(h.isotope, 'unknown') && ~h.resolved);

    y = model(f, 2816.33, 2.16, [1 1 1]*1.8e-3, 1.0, 'gauss', []) + 2e-5*(f - 2816) ...
        + sig*randn(size(f));
    h  = hfOf(mergeOpts(hq, struct('f', f, 'y', y)));
    hE = hfOf(mergeOpts(hq, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('no e: 3 dips 14N', h.nDips == 3 && strcmp(h.isotope, '14N'));
    check('no e: same A error as with the true error bars (within 10%)', ...
          abs(h.splittingErr/hE.splittingErr - 1) < 0.1);
    h = hfOf(mergeOpts(hq, struct('f', f, 'y', y, 'hfCount', 2)));
    check('forced: nDips = 2, confidence forced', h.nDips == 2 && strcmp(h.confidence, 'forced'));

    nWrong = 0; nConf = 0;
    for k = 1:20
        y = model(f, 2816.33, 2.16, [1 1 1]*4e-4, 1.0, 'lorentz', []) + sig*randn(size(f));
        o = esrFitRobust(mergeOpts(hq, struct('f', f, 'y', y, 'e', sig*ones(size(f)), ...
                                              'hfShape', 'lorentz')));
        if o.ok && o.hyperfine.tested && any(strcmp(o.hyperfine.confidence, {'strong', 'decisive'}))
            nConf = nConf + 1;
            if o.hyperfine.nDips ~= 3; nWrong = nWrong + 1; end
        end
    end
    fprintf('       low SNR: %d/20 confident, %d of them wrong\n', nConf, nWrong);
    check('low SNR: at most 1 confident wrong answer in 20', nWrong <= 1);

    nR = 60; pa = nan(nR, 1); pf = nan(nR, 1); cls = 0;
    for k = 1:nR
        y = model(f, 2816.33, 2.16, [1 1 1]*1.5e-3, 1.0, 'lorentz', []) + sig*randn(size(f));
        h = hfOf(mergeOpts(hq, struct('f', f, 'y', y, 'e', sig*ones(size(f)), 'hfShape', 'lorentz')));
        cls = cls + (h.nDips == 3);
        pa(k) = (h.splitting - 2.16) / h.splittingErr;
        pf(k) = (h.f0 - 2816.33) / h.f0Err;
    end
    fprintf('       pulls over %d: A mean %+.2f std %.2f;  f0 mean %+.2f std %.2f;  %d/%d called 3\n', ...
            nR, mean(pa), std(pa), mean(pf), std(pf), cls, nR);
    check('pulls: std of A pull in [0.75 1.3]', std(pa) > 0.75 && std(pa) < 1.3);
    check('pulls: std of f0 pull in [0.75 1.3]', std(pf) > 0.75 && std(pf) < 1.3);
    check('pulls: |mean| < 0.4', abs(mean(pa)) < 0.4 && abs(mean(pf)) < 0.4);
    check('pulls: all called 3 dips', cls == nR);

    % =====================================================================
    fprintf('\n2. HYPERFINE GATE (hyperfine = ''auto'')\n');
    y = model(f, 2816.33, 2.16, [1 1 1]*1.8e-3, 1.0, 'lorentz', []) + sig*randn(size(f));
    o = esrFitRobust(mergeOpts(q, struct('f', f, 'y', y, 'e', sig*ones(size(f)))));
    check('auto: a narrow line on a fine grid is tested, and 14N', o.hyperfine.tested && ...
          strcmp(o.hyperfine.isotope, '14N'));
    check('auto: out.lines stays one line per transition', o.nLines == 1);
    fw = (2840:1:2900)';
    y = 1 - 0.012*lor(fw, 2870.4, 10) + 5e-4*randn(size(fw));
    o = esrFitRobust(mergeOpts(q, struct('f', fw, 'y', y, 'e', 5e-4*ones(size(fw)))));
    check('auto: a 10 MHz CW line at 1 MHz steps is NOT tested', ~o.hyperfine.tested);
    % The 2026-11-26 pulsed sweep: three pi-pulse lines (tau 0.572 us) on a 0.2 MHz grid. Fitted
    % as ONE line its FWHM is the envelope, above the 3A = 6.5 MHz the gate once demanded.
    fp = (2752:0.2:2772)'; tp = 0.572; sp = 4e-5;
    y = 0.99975 + model(fp, 2762.85, 2.16, [1 1 1]*2.7e-3, 0.8, 'pulse', tp) - 1 + sp*randn(size(fp));
    o = esrFitRobust(mergeOpts(q, struct('f', fp, 'y', y, 'e', sp*ones(size(fp)), 'tau', tp)));
    check(sprintf('auto: a resolved pulsed triplet whose one-line FWHM is %.1f MHz is tested', ...
          o.lines(1).fwhm), o.lines(1).fwhm > 3*2.16 && o.hyperfine.tested);
    check('   ... and called 14N with A within 3 sigma', strcmp(o.hyperfine.isotope, '14N') && ...
          o.hyperfine.resolved && abs(o.hyperfine.splitting - 2.16) < 3*o.hyperfine.splittingErr);
    fw = (2840:0.5:2900)';
    y = 1 - 0.012*lor(fw, 2870.4, 10) + 5e-4*randn(size(fw));
    o = esrFitRobust(mergeOpts(q, struct('f', fw, 'y', y, 'e', 5e-4*ones(size(fw)))));
    check('auto: a 10 MHz line at 0.5 MHz steps is tested and NOT resolved', ...
          o.hyperfine.tested && o.hyperfine.nDips == 1 && ~o.hyperfine.resolved);

    % =====================================================================
    fprintf('\n3. RESONANCE\n');
    fz = (2770:1:2970)'; s = 1e-3;
    y = 1 + 2e-5*(fz - 2870) - 0.02*lor(fz, 2870.4, 12) - 0.012*lor(fz, 2763, 33) + s*randn(size(fz));
    o = esrFitRobust(mergeOpts(q, struct('f', fz, 'y', y, 'e', s*ones(size(fz)), 'nLines', 1)));
    check('zero field, resonance BELOW the sweep: the one NV line is at D', ...
          o.ok && o.nLines == 1 && abs(o.lines(1).centre - 2870.4) < 0.5);
    check('   ... and the resonance is found near the low edge', o.resonance.found && ...
          o.resonance.centre < 2780);
    check('   ... and the field reads zero', strcmp(o.field.kind, 'zero'));

    y = 1 - 0.02*lor(fz, 2870.4, 12) - 0.012*lor(fz, 2800, 33) + s*randn(size(fz));
    o = esrFitRobust(mergeOpts(q, struct('f', fz, 'y', y, 'e', s*ones(size(fz)))));
    check('zero field, resonance INSIDE (2800, FWHM 33): found', o.resonance.found && ...
          abs(o.resonance.centre - 2800) < 3 && abs(o.resonance.fwhm/33 - 1) < 0.3);
    check('   ... and only the NV line is in out.lines', o.nLines == 1 && ...
          abs(o.lines(1).centre - 2870.4) < 0.5);

    y = 1 - 0.02*lor(fz, 2870.4, 12) - 0.02*lor(fz, 2780, 30) + s*randn(size(fz));
    o = esrFitRobust(mergeOpts(q, struct('f', fz, 'y', y, 'e', s*ones(size(fz)), 'nLines', 1)));
    check('zero field, resonance AS DEEP as the line, nLines = 1: the slot goes to D', ...
          o.ok && abs(o.lines(1).centre - 2870.4) < 0.5 && o.resonance.found);

    nFalse = 0;
    for k = 1:10
        y = 1 + 2e-5*(fz - 2870) - 0.02*lor(fz, 2870.4, 12) + s*randn(size(fz));
        o = esrFitRobust(mergeOpts(q, struct('f', fz, 'y', y, 'e', s*ones(size(fz)))));
        nFalse = nFalse + o.resonance.found;
    end
    check(sprintf('no resonance in the data: none reported (%d/10 false)', nFalse), nFalse == 0);

    y = 1 - 0.02*lor(fz, 2870.4, 12) - 0.012*lor(fz, 2800, 33) + s*randn(size(fz));
    o = esrFitRobust(mergeOpts(q, struct('f', fz, 'y', y, 'e', s*ones(size(fz)), 'resonance', false)));
    check('resonance = false: none reported', ~o.resonance.found);
    o = esrFitRobust(mergeOpts(q, struct('f', fz, 'y', y, 'e', s*ones(size(fz)), 'resonance', 2805)));
    check('resonance = 2805 (a hint): found at 2800', o.resonance.found && ...
          abs(o.resonance.centre - 2800) < 3);

    % =====================================================================
    fprintf('\n4. FIELD\n');
    sy = esrSynthSpectrum(struct('B', 45, 'tiltDeg', 7, 'verbose', false));
    o = esrFitRobust(mergeOpts(q, struct('exp', sy)));
    check(sprintf('45 G, 7 deg: |B| %.3f, tilt %.2f, D %.2f', o.B, o.tiltDeg, o.D), ...
          strcmp(o.field.kind, 'vector') && abs(o.B - 45) < 0.2 && abs(o.tiltDeg - 7) < 0.5);
    check('   ... the two aligned-axis lines are flagged', nnz(o.field.aligned) == 2 && ...
          all(abs(sort([o.lines(o.field.aligned).centre]) - [2745.45 2995.94]) < 1));
    yr = sy.signal(:) - 0.008*lor(sy.freqs(:), 2763, 33);
    o = esrFitRobust(mergeOpts(q, struct('f', sy.freqs, 'y', yr, 'e', sy.sterr)));
    check(sprintf('45 G + a resonance at 2763: |B| %.3f, tilt %.2f', o.B, o.tiltDeg), ...
          abs(o.B - 45) < 0.2 && abs(o.tiltDeg - 7) < 0.5);
    check(sprintf('   ... resonance found at %.1f MHz', o.resonance.centre), ...
          o.resonance.found && abs(o.resonance.centre - 2763) < 5);
    check('   ... and no line in out.lines within 10 MHz of it', ...
          ~any(abs([o.lines.centre] - 2763) < 10));

    fp = (2800:0.5:2940)';
    y = 1 - 0.012*lor(fp, 2870.5 - 28.0, 9) - 0.012*lor(fp, 2870.5 + 28.0, 9) + 5e-4*randn(size(fp));
    o = esrFitRobust(mergeOpts(q, struct('f', fp, 'y', y, 'e', 5e-4*ones(size(fp)))));
    check(sprintf('one pair, 56 MHz apart: B_par %.3f G (truth %.3f)', o.field.Bpar, 28/2.8025), ...
          strcmp(o.field.kind, 'pair') && abs(o.field.Bpar - 28/2.8025) < 0.05);

    fn = (2771:0.5:2801)';
    y = 1 - 0.01*lor(fn, 2786, 10) + 1e-3*randn(size(fn));
    o = esrFitRobust(mergeOpts(q, struct('f', fn, 'y', y, 'e', 1e-3*ones(size(fn)), 'baselineOrder', 1)));
    check('a 30 MHz window on one line: one line, no resonance, no field', o.ok && ...
          o.nLines == 1 && ~o.resonance.found && strcmp(o.field.kind, 'none'));

    % =====================================================================
    fprintf('\n5. INTERFACE\n');
    threw = false;
    try; esrFitRobust(struct('f', fz, 'y', y, 'nosuch', 1)); catch; threw = true; end
    check('an unknown option is an error', threw);
    threw = false;
    try; esrFitRobust(fz, y); catch; threw = true; end                           %#ok<NOCOM>
    check('the positional form is an error', threw);
    o = esrFitRobust(mergeOpts(q, struct('exp', [fp, 1 - 0.01*lor(fp, 2850, 9) + 5e-4*randn(size(fp))])));
    check('an [f y] matrix as exp works', o.ok && abs(o.lines(1).centre - 2850) < 1);
    check('modelAt(f) equals model', max(abs(o.modelAt(o.f) - o.model)) < 1e-12);

    % =====================================================================
    if ~isempty(realFile)
        fprintf('\n6. A REAL SAVE\n');
        o = esrFitRobust(mergeOpts(hq, struct('file', realFile)));
        fprintf('       real file: %s\n', o.hyperfine.verdict);
        check('real file: 3 dips 14N', o.hyperfine.nDips == 3 && strcmp(o.hyperfine.isotope, '14N'));
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

function h = hfOf(args)
% The hyperfine part of a fit, with f/model/residual over its own window for the checks.
    o = esrFitRobust(args);
    h = o.hyperfine;
    if ~h.tested
        fprintf('       (hyperfine not tested: %s%s)\n', o.why, h.why);
        h.nDips = NaN; h.isotope = ''; h.confidence = ''; h.resolved = false;
    end
end

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
