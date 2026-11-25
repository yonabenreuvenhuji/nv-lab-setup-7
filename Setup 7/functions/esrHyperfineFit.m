function out = esrHyperfineFit(args)
%ESRHYPERFINEFIT One model for the whole ESR sweep, and the decision: 1, 2 or 3 dips.
%
%   out = esrHyperfineFit(struct('exp', p))        % a live AWGExpESR / AWGExpESRPulsed
%   out = esrHyperfineFit(struct('file', 'G:\...\_AutoSave\AWG ESR Pulsed_20261110_172200.mat'))
%   out = esrHyperfineFit(struct('f', f, 'y', y, 'e', e))
%   out = esrHyperfineFit(struct('exp', p, 'nLines', 3))     % skip the choice, force 3
%   yy  = out.modelAt(ff);                                   % the fitted curve, any ff
%
% ONE STRUCT IN, ONE STRUCT OUT, like esrFitRobust. out.ok says whether to believe it.
%
% esrMultiLineFit is NOT the tool for this: it keeps candidate lines >= minSepMHz (6 MHz)
% apart and gives every line its own free centre and width, so a 2.16 MHz triplet comes
% back as one wide line. This fitter instead builds the hyperfine structure INTO the model.
%
% ===================================================================
% THE MODEL (one function for the whole plot)
% ===================================================================
%
%   y(f) = P(f) - SUM_k C_k * S(f - f0 - m_k*A; w)
%
%   3 dips   14N, I = 1     m = [-1 0 +1]      A ~ 2.16 MHz,  f0 = the mI = 0 line
%   2 dips   15N, I = 1/2   m = [-1/2 +1/2]    A ~ 3.03 MHz,  f0 = their midpoint
%   1 dip    hyperfine not resolved            m = 0
%   0        no line: P only
%
% P is a polynomial baseline in (f - xref) of order baselineOrder (default 1). All lines
% share ONE width w and ONE shape S. The spacing is EXACTLY equal: the quadrupole term
% Q*mI^2 is the same in ms = 0 and ms = +-1 and cancels in the transition. Second-order
% terms (A_perp^2/D ~ 3 kHz) are far below anything this setup resolves.
%
% Depths: equal is the physics away from the ~500 G level anticrossings (unpolarised
% nitrogen). depths = 'auto' (default) fits every line count BOTH ways, equal and free,
% and BIC picks, charging the N-1 extra depths. That matters for the COUNT, not just the
% depths: a polarised triplet (e.g. 2.6 : 1.6 : 0.7) fitted with equal depths loses to a
% 2-dip model. 'equal' or 'free' forces one variant.
%
% A is FITTED inside splitRange (default 1-5 MHz), not fixed to 2.16 or 3.03. That keeps
% the 2-vs-3 decision honest: both models get the same freedom, and the fitted A is then
% compared with the nitrogen values as a separate check. Fix it with 'splitting'.
%
% ===================================================================
% LINE SHAPES
% ===================================================================
%
%   'lorentz'  1/(1 + (2u/w)^2)                    CW, power broadened
%   'gauss'    exp(-4 ln2 u^2/w^2)                 inhomogeneous
%   'voigt'    pseudo-Voigt eta*L + (1-eta)*G, same FWHM, eta fitted in [0 1]
%   'pulse'    square pulse, length tau, on-resonance rotation theta (default 180 deg):
%                 Pr(u) = W^2/(W^2+u^2) sin^2(pi tau sqrt(W^2+u^2)),   W = theta/(2 pi tau)
%              convolved with a Gaussian of FWHM w. Here w is the INHOMOGENEOUS width, C
%              is the contrast at full inversion, and the line's own FWHM (pulse and
%              inhomogeneous together) is what out.fwhm reports.
%   'auto'     (default) lorentz and gauss, plus pulse when tau is known (read from a
%              pulsed experiment or save, or opts.tau). Same parameter count, so the lowest
%              chi2 wins for each line count.
%
% ===================================================================
% HOW IT DECIDES 0 / 1 / 2 / 3
% ===================================================================
%
% Every candidate is fitted to the WHOLE sweep. Start: a grid over (f0, A, w) with the
% baseline and the depth solved linearly at every grid point, so the start is the global
% minimum of that grid and nothing depends on peak finding. Then a bounded Levenberg-
% Marquardt on the nonlinear parameters with the linear ones re-solved at every step
% (variable projection), depths kept >= 0.
%
% The candidates are compared by BIC = chi2/s^2 + k ln(n). s^2 is the reduced chi2 of the
% best-fitting candidate (two-sided Birge, as in esrMultiLineFit: the error bars give the
% shape, the residuals give the scale). With equal depths the 2- and 3-dip models have the
% SAME number of parameters, so between those two it is plain delta-chi2. Strength from
% delta-BIC to the runner-up (Kass & Raftery): < 2 ambiguous, 2-6 positive, 6-10 strong,
% > 10 decisive. out.prob holds the BIC weights.
%
% Label: 3 dips with A within max(3 sigma, 0.3 MHz) of 2.16 -> '14N'; 2 dips within that
% of 3.03 -> '15N'; 1 dip -> 'unresolved'; 0 -> 'none'; anything else 'unknown', with the
% fitted spacing printed (e.g. two orientations a few MHz apart).
%
% ===================================================================
% ERRORS
% ===================================================================
%
% Full covariance (J'WJ)^-1 of the final fit, scaled by its Birge factor
% max(sqrt(chi2_nu), 0.05) in both directions -- the esrMultiLineFit convention.
% Line centres carry the f0-A covariance: var(f0 + m A) = s0^2 + m^2 sA^2 + 2 m c0A.
%
% ===================================================================
% OPTIONS (struct fields, all optional)
% ===================================================================
%   exp            live experiment (read with esrPullData), or a loaded save struct
%   file           path to a SaveLoad .mat: uses its summaryData, reads tau from it
%   f, y, e        plain vectors (e optional; without it the residuals set the scale)
%   range          [fmin fmax] MHz: fit only this part of the sweep (default [] = all).
%                  The model is ONE multiplet, so crop a wide sweep to one line group.
%   nLines         force 0, 1, 2 or 3 and skip the choice (default [])
%   nLinesList     candidates (default [1 2 3]; "no line" is always added)
%   shape          'auto' | 'lorentz' | 'gauss' | 'voigt' | 'pulse'      (default 'auto')
%   depths         'auto' | 'equal' | 'free'                             (default 'auto')
%   splitting      MHz, fix A instead of fitting it (default [])
%   splitRange     [min max] MHz allowed for A (default [1 5])
%   baselineOrder  0, 1 or 2 (default 1)
%   maxFwhmMHz     widest line accepted (default 40)
%   tau            us, pulse length for 'pulse' (default: read from the experiment)
%   pulseRotationDeg   on-resonance rotation of that pulse (default 180)
%   plotIt (true), showResid (false), verbose (true), name, yLabel,
%   avgOutlierMad (4, passed to esrPullData), maxIter (200)
%
% OUTPUT (struct)
%   ok, why, nDips, isotope, confidence, verdict, warnings
%   f0, f0Err, splitting, splittingErr, fwhm, fwhmErr, depth, depthErr, shape, eta,
%   widthParam (w as fitted: the inhomogeneous FWHM for 'pulse')
%   lines        struct array: centre, centreErr, fwhm, fwhmErr, depth, depthErr
%                (the same fields as esrMultiLineFit / esrFitRobust)
%   centreCov    covariance of the line centres, MHz^2
%   candidates   every model tried: nDips, shape, chi2, chi2Reduced, k, bic, dBic, prob,
%                params, modelAt
%   prob         [P(0) P(1) P(2) P(3)], BIC weights (NaN where not tried)
%   depthTest    pEqual (chi-square p-value that the depths are equal, for the chosen
%                count) and used ('equal' | 'free')
%   model, residual, modelAt, formula, params, cov, chi2Reduced, birge, f, y, e, opts
%
% Base MATLAB only: esrPullData for live experiments, fitStyle for the figure when it is on
% the path (a plain figure otherwise). Runs in Octave. Offline grader: testHyperfineFit.

    A14N = 2.16;   % |A_par| of 14N, MHz
    A15N = 3.03;   % |A_par| of 15N, MHz

    DEFAULTS = struct('exp', [], 'file', '', 'f', [], 'y', [], 'e', [], ...
        'nLines', [], 'nLinesList', [1 2 3], 'shape', 'auto', 'depths', 'auto', ...
        'splitting', [], 'splitRange', [1 5], 'baselineOrder', 1, 'maxFwhmMHz', 40, ...
        'tau', [], 'pulseRotationDeg', 180, 'avgOutlierMad', 4, 'maxIter', 200, ...
        'range', [], 'plotIt', true, 'showResid', false, 'verbose', true, 'name', '', ...
        'yLabel', '');

    if ~exist('args', 'var') || isempty(args); args = struct(); end
    if ~isstruct(args)
        error(['esrHyperfineFit takes ONE struct: esrHyperfineFit(struct(''exp'', p)), ' ...
               'esrHyperfineFit(struct(''file'', path)) or ' ...
               'esrHyperfineFit(struct(''f'', f, ''y'', y, ''e'', e)).']);
    end
    opts = withDefaults(args, DEFAULTS, 'esrHyperfineFit');
    opts.shape = lower(opts.shape);
    opts.depths = lower(opts.depths);
    if ~any(strcmp(opts.depths, {'auto', 'equal', 'free'}))
        error('esrHyperfineFit: depths must be ''auto'', ''equal'' or ''free''.');
    end
    if ~isempty(opts.nLines) && ~any(opts.nLines == [0 1 2 3])
        error('esrHyperfineFit: nLines must be 0, 1, 2 or 3.');
    end
    if any(~ismember(opts.nLinesList, [1 2 3]))
        error('esrHyperfineFit: nLinesList may only contain 1, 2 and 3.');
    end

    % Rank-deficient solves are expected and handled (a line pushed between samples), so
    % their warnings are noise here. Restored on exit, whatever happens.
    ids = {'MATLAB:rankDeficientMatrix', 'MATLAB:singularMatrix', ...
           'MATLAB:nearlySingularMatrix', 'Octave:singular-matrix', ...
           'Octave:nearly-singular-matrix'};
    wst = warning;
    for i = 1:numel(ids); warning('off', ids{i}); end
    restoreWarn = onCleanup(@() warning(wst));                                   %#ok<NASGU>

    out = emptyOut(opts);

    % ---- the data -------------------------------------------------------
    try
        [x, y, e, tau, src] = getData(opts);
    catch err
        out.why = err.message; out = failed(out, opts); return
    end
    if ~isempty(opts.tau); tau = opts.tau; end
    x = double(x(:)); y = double(y(:)); e = double(e(:));
    if numel(e) ~= numel(x); e = nan(size(x)); end
    good = isfinite(x) & isfinite(y);
    if ~isempty(opts.range)
        good = good & x >= min(opts.range) & x <= max(opts.range);
    end
    x = x(good); y = y(good); e = e(good);
    [x, ix] = sort(x); y = y(ix); e = e(ix);
    e(~isfinite(e) | e <= 0) = NaN;
    n = numel(x);
    out.f = x; out.y = y; out.e = e; out.source = src;
    if n < opts.baselineOrder + 8
        out.why = sprintf('only %d usable points', n); out = failed(out, opts); return
    end
    ee = e;
    if any(isfinite(ee))
        ee(~isfinite(ee)) = median(ee(isfinite(ee)));
    else
        % No error bars: a robust point-to-point noise estimate, so chi2 stays in sigma
        % units (the Birge floor and the BIC scale assume that). Successive differences
        % cancel smooth structure; median|diff| = 0.6745*sqrt(2)*sigma for Gaussian noise.
        % Birge still corrects the scale from the residuals afterwards.
        sg = 1.4826 * median(abs(diff(y))) / sqrt(2);
        if ~isfinite(sg) || sg <= 0; sg = std(y); end
        if ~isfinite(sg) || sg <= 0; sg = 1; end
        ee = sg * ones(n, 1);
    end

    D = struct('x', x, 'y', y, 'sw', 1 ./ ee, 'xr', mean(x), ...
               'nB', opts.baselineOrder + 1, 'dx', median(diff(x)), ...
               'span', max(x) - min(x));

    shapes = shapeList(opts.shape, tau);
    P = pulseParams(tau, opts.pulseRotationDeg);
    out.tau = tau;

    % ---- fit every candidate ---------------------------------------------
    if isempty(opts.nLines)
        Ns = unique([0, opts.nLinesList(:)']);
    else
        Ns = unique([0, opts.nLines]);
    end
    freeFromStart = strcmp(opts.depths, 'free');
    tryFree = strcmp(opts.depths, 'auto');
    vars = [];                       % every fit: per line count, best shape, equal (+ free)
    for N = Ns
        best = [];
        shp = shapes;
        if N == 0; shp = {'none'}; end
        for s = 1:numel(shp)
            F = fitCandidate(D, N, shp{s}, freeFromStart, opts, P, []);
            if isempty(best) || F.chi2 < best.chi2; best = F; end
        end
        vars = [vars, best];                                                  %#ok<AGROW>
        if tryFree && N >= 2
            % Polarised nitrogen (or a drive that differs across the lines) makes the depths
            % unequal, and an equal-depth model can then lose to the WRONG line count. So
            % each count also competes with free depths; BIC charges the N-1 extra depths.
            vars = [vars, fitCandidate(D, N, best.shape, true, opts, P, best.theta)]; %#ok<AGROW>
        end
    end

    % ---- choose ------------------------------------------------------------
    k    = [vars.k];
    chi2 = [vars.chi2];
    nu   = chi2 ./ max(n - k, 1);
    isLine = [vars.nDips] > 0;
    s2 = min(nu(isLine));
    if isempty(s2); s2 = nu(1); end
    s2 = max(s2, 0.05^2);
    bicV = chi2 / s2 + k * log(n);
    for i = 1:numel(vars)
        vars(i).chi2Reduced = nu(i); vars(i).bic = bicV(i);
    end
    % one candidate per line count: its better depth variant
    cands = [];
    for N = Ns
        idx = find([vars.nDips] == N);
        [~, j] = min(bicV(idx));
        cands = [cands, vars(idx(j))];                                       %#ok<AGROW>
    end
    bic  = [cands.bic];
    dBic = bic - min(bic);
    pr   = exp(-dBic/2); pr = pr / sum(pr);
    for i = 1:numel(cands)
        cands(i).dBic = dBic(i); cands(i).prob = pr(i);
    end
    out.prob = nan(1, 4);
    out.prob([cands.nDips] + 1) = pr;

    if isempty(opts.nLines)
        [~, iw] = min(bic);
    else
        iw = find([cands.nDips] == opts.nLines, 1);
    end
    final = cands(iw);
    others = setdiff(1:numel(cands), iw);
    if isempty(others)
        gap = Inf; iRun = [];
    else
        [gap, j] = min(bic(others) - bic(iw));
        iRun = others(j);
    end
    out.candidates = cands;

    % ---- equal depths? (reported; the choice above already made it) -----------
    out.depthTest = struct('pEqual', NaN, 'used', ternary(final.freeDepths, 'free', 'equal'));
    if tryFree && final.nDips >= 2
        idx = find([vars.nDips] == final.nDips);
        eq = vars(idx(~[vars(idx).freeDepths]));
        fr = vars(idx([vars(idx).freeDepths]));
        out.depthTest.pEqual = chi2sf(max(eq.chi2 - fr.chi2, 0)/s2, fr.k - eq.k);
    end

    % ---- covariance and the numbers ----------------------------------------
    R = finalNumbers(D, final, opts, P);
    if ~R.ok
        out.why = R.why; out = failed(out, opts); return
    end
    out = copyFields(out, R);
    out.nDips = final.nDips;
    out.shape = final.shape;
    out.modelAt = final.modelAt;
    out.model = final.yhat;
    out.residual = y - final.yhat;

    % ---- verdict ------------------------------------------------------------
    nd = final.nDips;
    tolA = max(3*out.splittingErr, 0.3);     % how close A must be to the nitrogen value
    switch nd
        case 0
            out.isotope = 'none';
        case 1
            out.isotope = 'unresolved';
        case 2
            if abs(out.splitting - A15N) <= tolA; out.isotope = '15N';
            else; out.isotope = 'unknown'; end
        case 3
            if abs(out.splitting - A14N) <= tolA; out.isotope = '14N';
            else; out.isotope = 'unknown'; end
    end
    if ~isempty(opts.nLines)
        out.confidence = 'forced';
    elseif ~isfinite(gap) || gap > 10
        out.confidence = 'decisive';
    elseif gap > 6
        out.confidence = 'strong';
    elseif gap > 2
        out.confidence = 'positive';
    else
        out.confidence = 'ambiguous';
    end
    out.runnerUp = [];
    if ~isempty(iRun); out.runnerUp = cands(iRun); end
    out.verdict = verdictText(out, gap, A14N, A15N);

    out.warnings = collectWarnings(out, final, D, opts, P, A14N, A15N);
    out.ok = true;
    out.opts = opts;

    if opts.verbose; report(out, A14N, A15N); end
    if opts.plotIt;  drawIt(out); end
end

% ======================================================================
% FITTING
% ======================================================================

function F = fitCandidate(D, N, shp, freeDepths, opts, P, th0)
% One line count, one shape: grid start, bounded VarPro Levenberg-Marquardt.
    F = struct('nDips', N, 'shape', shp, 'freeDepths', freeDepths && N >= 2, ...
               'm', hfOffsets(N), 'names', {{}}, 'theta', [], 'lo', [], 'hi', [], ...
               'beta', [], 'chi2', Inf, 'yhat', [], 'k', D.nB, 'chi2Reduced', NaN, ...
               'bic', NaN, 'dBic', NaN, 'prob', NaN, 'params', struct(), ...
               'modelAt', []);
    [F.names, F.lo, F.hi, hs] = thetaSpec(D, N, shp, opts, P);
    if N > 0
        if isempty(th0)
            th0 = gridStart(D, F, opts, P);
        end
        th0 = min(max(th0(:), F.lo(:)), F.hi(:));
        rfun = @(th) vpResid(th, D, F, opts, P);
        F.theta = lmBounded(rfun, th0, F.lo(:), F.hi(:), hs(:), opts.maxIter);
    end
    [F.chi2, F.beta, F.yhat] = vpSolve(F.theta, D, F, opts, P);
    nDepth = 0;
    if N > 0; nDepth = ternary(F.freeDepths, N, 1); end
    F.k = numel(F.theta) + D.nB + nDepth;
    F.params = paramStruct(F, D, opts);
    Fc = F;      % the handle keeps its own copy of the parameters
    F.modelAt = @(ff) curveAt(ff, Fc, D, opts, P);
end

function [names, lo, hi, hs] = thetaSpec(D, N, shp, opts, P)
    names = {}; lo = []; hi = []; hs = [];
    if N == 0; return; end
    wLo = 0.5*D.dx;
    if strcmp(shp, 'pulse'); wLo = P.wMin; end
    wHi = max(min(D.span, opts.maxFwhmMHz), 2*wLo);
    hx  = 1e-4 * max(D.dx, 0.01);
    names = {'f0', 'w'}; lo = [min(D.x), wLo]; hi = [max(D.x), wHi]; hs = [hx, hx];
    if N >= 2 && isempty(opts.splitting)
        names{end+1} = 'A';
        lo(end+1) = opts.splitRange(1); hi(end+1) = opts.splitRange(2); hs(end+1) = hx;
    end
    if strcmp(shp, 'voigt')
        names{end+1} = 'eta'; lo(end+1) = 0; hi(end+1) = 1; hs(end+1) = 1e-4;
    end
end

function th0 = gridStart(D, F, opts, P)
% Global start on a grid. For every (w, eta, A) the baseline and the (common) depth are
% solved in closed form for ALL f0 at once, by projecting out the baseline.
    x = D.x; sw = D.sw; names = F.names;
    [Q, ~] = qr(baseCols(x, D.xr, D.nB) .* sw, 0);
    yw = D.y .* sw;
    yp = yw - Q*(Q'*yw);
    chi0 = yp'*yp;
    fcs = linspace(F.lo(1), F.hi(1), max(ceil(3*D.span/D.dx), 30) + 1);
    wTop = min(F.hi(2), max(D.span/2, 2*F.lo(2)));
    ws = exp(linspace(log(F.lo(2)), log(wTop), 16));
    iA = find(strcmp(names, 'A'));
    iE = find(strcmp(names, 'eta'));
    As = NaN;  if ~isempty(iA); As = F.lo(iA):0.05:F.hi(iA); end
    ets = NaN; if ~isempty(iE); ets = [0 0.5 1]; end
    best = Inf;
    th0 = (F.lo + F.hi)/2;
    for w = ws
        for et = ets
            for A = As
                th = th0; th(2) = w;
                if ~isempty(iA); th(iA) = A; end
                if ~isempty(iE); th(iE) = et; end
                [~, ~, Av, eta] = unpackTheta(th, names, opts);
                S = zeros(numel(x), numel(fcs));
                for k = 1:numel(F.m)
                    S = S + shapeAt(x - (fcs + F.m(k)*Av), w, eta, F.shape, P);
                end
                Sw = S .* sw;
                Sp = Sw - Q*(Q'*Sw);
                num = -(yp' * Sp);
                den = sum(Sp.^2, 1);
                ok  = den > 0 & num > 0;          % a dip, not a peak
                chi = chi0 * ones(size(num));
                chi(ok) = chi0 - num(ok).^2 ./ den(ok);
                [cm, i] = min(chi);
                if cm < best
                    best = cm; th0 = th; th0(1) = fcs(i);
                end
            end
        end
    end
end

function th = lmBounded(rfun, th, lo, hi, hs, maxIter)
% Levenberg-Marquardt with Marquardt column scaling; every trial step is clamped to the box.
    r = rfun(th); c = r'*r;
    lambda = 1e-2;
    np = numel(th);
    for it = 1:maxIter
        J = numJac(rfun, th, hs, lo, hi);
        H = J'*J; g = J'*r;
        sc = sqrt(max(diag(H), realmin));
        Hs = H ./ (sc*sc'); gs = g ./ sc;
        took = false; rel = 0;
        for tries = 1:25
            d  = -((Hs + lambda*eye(np)) \ gs) ./ sc;
            tn = min(max(th + d, lo), hi);
            if all(tn == th); lambda = lambda*5; continue; end
            rn = rfun(tn); cn = rn'*rn;
            if isfinite(cn) && cn < c
                rel = (c - cn) / max(c, realmin);
                th = tn; r = rn; c = cn;
                lambda = max(lambda/5, 1e-10);
                took = true;
                break
            end
            lambda = lambda*5;
        end
        if ~took || rel < 1e-11; break; end
    end
end

function J = numJac(fun, th, hs, lo, hi)
% Central differences, one-sided against a bound.
    f0 = fun(th);
    J = zeros(numel(f0), numel(th));
    for j = 1:numel(th)
        tp = th; tm = th;
        tp(j) = min(th(j) + hs(j), hi(j));
        tm(j) = max(th(j) - hs(j), lo(j));
        J(:, j) = (fun(tp) - fun(tm)) / (tp(j) - tm(j));
    end
end

function r = vpResid(th, D, F, opts, P)
    [~, ~, yhat] = vpSolve(th, D, F, opts, P);
    r = D.sw .* (D.y - yhat);
end

function [chi2, beta, yhat, X] = vpSolve(th, D, F, opts, P)
% The linear part for given nonlinear parameters: baseline free, depths >= 0.
    B = baseCols(D.x, D.xr, D.nB);
    Xd = depthCols(D.x, th, F, opts, P);
    [beta, chi2] = nnDepthLS(B, Xd, D.y, D.sw);
    X = [B, Xd];
    yhat = X * beta;
end

function Xd = depthCols(x, th, F, opts, P)
    if F.nDips == 0
        Xd = zeros(numel(x), 0); return
    end
    L = lineMatrix(x, th, F, opts, P);
    if F.freeDepths; Xd = -L; else; Xd = -sum(L, 2); end
end

function L = lineMatrix(x, th, F, opts, P)
    [f0, w, A, eta] = unpackTheta(th, F.names, opts);
    U = x(:) - (f0 + F.m(:)' * A);
    L = shapeAt(U, w, eta, F.shape, P);
end

function [beta, chi2] = nnDepthLS(B, Xd, y, sw)
% Least squares with the depth columns constrained >= 0. At most 3 depths, so the exact
% answer is the best feasible subset; the unconstrained solve is tried first and is the
% answer whenever it is already feasible.
    nB = size(B, 2); nd = size(Xd, 2);
    X = [B, Xd];
    b = (X .* sw) \ (y .* sw);
    if all(b(nB+1:end) >= 0)
        beta = b; r = sw .* (y - X*b); chi2 = r'*r; return
    end
    chi2 = Inf; beta = zeros(nB + nd, 1);
    for mask = 0:2^nd - 2
        on = logical(bitget(mask, 1:nd));
        Xs = [B, Xd(:, on)];
        bs = (Xs .* sw) \ (y .* sw);
        if any(bs(nB+1:end) < 0); continue; end
        r = sw .* (y - Xs*bs); ch = r'*r;
        if ch < chi2
            chi2 = ch;
            beta = zeros(nB + nd, 1); beta(1:nB) = bs(1:nB);
            beta(nB + find(on)) = bs(nB+1:end);
        end
    end
end

function R = finalNumbers(D, F, opts, P)
% Full covariance over [theta; baseline; depths], Birge-scaled, and the derived numbers.
    R = struct('ok', false, 'why', '');
    n  = numel(D.x);
    nB = D.nB;
    nth = numel(F.theta);
    p  = [F.theta(:); F.beta(:)];
    np = numel(p);
    dof = max(n - F.k, 1);
    birge = max(sqrt(F.chi2 / dof), 0.05);

    % Jacobian of the model: numeric in theta, exact in the linear parameters.
    B  = baseCols(D.x, D.xr, nB);
    Xd = depthCols(D.x, F.theta, F, opts, P);
    J  = [zeros(n, nth), B, Xd];
    if nth > 0
        [~, ~, ~, hs] = thetaSpec(D, F.nDips, F.shape, opts, P);
        mf = @(th) [B, depthCols(D.x, th, F, opts, P)] * F.beta;
        J(:, 1:nth) = numJac(mf, F.theta(:), hs(:), F.lo(:), F.hi(:));
    end
    % A depth clamped at 0 is fixed, not fitted: leave it out of the covariance.
    free = true(np, 1);
    dIdx = nth + nB + (1:size(Xd, 2));
    free(dIdx(F.beta(nB+1:end) <= 0)) = false;
    Jw = J(:, free) .* D.sw;
    H  = Jw' * Jw;
    sc = sqrt(max(diag(H), realmin));
    Hs = H ./ (sc*sc');
    if rcond(Hs) < 1e-13
        Cf = pinv(Hs) ./ (sc*sc');
        R.degenerate = true;
    else
        Cf = (Hs \ eye(size(Hs))) ./ (sc*sc');
        R.degenerate = false;
    end
    C = zeros(np); C(free, free) = Cf * birge^2;
    if ~all(isfinite(C(:)))
        R.why = 'the covariance came out non-finite'; return
    end
    se = sqrt(max(diag(C), 0));

    R.params = p; R.cov = C; R.birge = birge;
    R.chi2Reduced = F.chi2 / dof;
    R.baseline = F.beta(1:nB)'; R.xref = D.xr;
    R.atBound = false(1, nth);
    for j = 1:nth
        R.atBound(j) = abs(F.theta(j) - F.lo(j)) < 1e-9*max(1, abs(F.lo(j))) || ...
                       abs(F.theta(j) - F.hi(j)) < 1e-9*max(1, abs(F.hi(j)));
    end
    R.thetaNames = F.names;

    R.f0 = NaN; R.f0Err = NaN; R.splitting = NaN; R.splittingErr = NaN;
    R.fwhm = NaN; R.fwhmErr = NaN; R.widthParam = NaN; R.widthParamErr = NaN;
    R.eta = NaN; R.depth = NaN; R.depthErr = NaN; R.lines = emptyLines();
    R.centreCov = []; R.formula = baselineText(R.baseline, D.xr);
    if F.nDips == 0
        R.ok = true; return
    end

    iF = find(strcmp(F.names, 'f0')); iW = find(strcmp(F.names, 'w'));
    iA = find(strcmp(F.names, 'A'));  iE = find(strcmp(F.names, 'eta'));
    [f0, w, A, eta] = unpackTheta(F.theta, F.names, opts);
    R.f0 = f0; R.f0Err = se(iF);
    R.widthParam = w; R.widthParamErr = se(iW);
    if F.nDips >= 2
        R.splitting = A;
        if isempty(iA); R.splittingErr = 0; else; R.splittingErr = se(iA); end
    end
    if ~isempty(iE); R.eta = eta; end

    % The single-line FWHM and peak, numerically for 'pulse'. Its error comes through w.
    [fw, pk] = effWidth(w, eta, F.shape, P);
    dw = max(1e-4*w, 1e-6);
    fw2 = effWidth(w + dw, eta, F.shape, P);
    R.fwhm = fw; R.fwhmErr = abs((fw2 - fw)/dw) * se(iW);
    R.peak = pk;

    % Line centres with the f0-A covariance.
    m = F.m(:);
    G = zeros(numel(m), np);
    G(:, iF) = 1;
    if ~isempty(iA); G(:, iA) = m; end
    R.centreCov = G * C * G';
    cen = f0 + m*A;
    cErr = sqrt(max(diag(R.centreCov), 0));
    dep = F.beta(nB+1:end) * pk;
    dErr = se(dIdx) * pk;
    if ~F.freeDepths
        dep = repmat(dep, numel(m), 1); dErr = repmat(dErr, numel(m), 1);
    end
    R.depth = mean(dep);
    if F.freeDepths
        dv = dIdx;
        R.depthErr = sqrt(sum(sum(C(dv, dv)))) * pk / numel(dv);
    else
        R.depthErr = dErr(1);
    end
    L = emptyLines();
    for kk = 1:numel(m)
        L(kk) = struct('centre', cen(kk), 'centreErr', cErr(kk), 'fwhm', fw, ...
                       'fwhmErr', R.fwhmErr, 'depth', dep(kk), 'depthErr', dErr(kk));
    end
    R.lines = L;
    R.formula = modelText(R, F);
    R.ok = true;
end

% ======================================================================
% MODEL PIECES
% ======================================================================

function m = hfOffsets(N)
    switch N
        case 0; m = zeros(1, 0);
        case 1; m = 0;
        case 2; m = [-0.5 0.5];
        case 3; m = [-1 0 1];
    end
end

function [f0, w, A, eta] = unpackTheta(th, names, opts)
    f0 = NaN; w = NaN; eta = 0.5;
    A = opts.splitting; if isempty(A); A = 0; end
    for j = 1:numel(names)
        switch names{j}
            case 'f0';  f0 = th(j);
            case 'w';   w = th(j);
            case 'A';   A = th(j);
            case 'eta'; eta = th(j);
        end
    end
end

function B = baseCols(x, xr, nB)
    x = x(:) - xr;
    B = ones(numel(x), nB);
    for j = 2:nB; B(:, j) = x.^(j-1); end
end

function S = shapeAt(U, w, eta, shp, P)
% Peak-1 line shapes (for 'pulse': peak 1 at zero inhomogeneous width and a true pi pulse).
    switch shp
        case 'lorentz'
            S = 1 ./ (1 + (2*U/w).^2);
        case 'gauss'
            S = exp(-4*log(2)*(U/w).^2);
        case 'voigt'
            S = eta ./ (1 + (2*U/w).^2) + (1 - eta)*exp(-4*log(2)*(U/w).^2);
        case 'pulse'
            S = pulseShape(U, w, P);
        otherwise
            error('esrHyperfineFit: unknown shape %s', shp);
    end
end

function S = pulseShape(U, w, P)
% Square-pulse response convolved with a Gaussian of FWHM w, tabulated on a fine grid and
% splined. The table is cached: the grid search and every LM step at the same w reuse it.
    persistent key coefs u0 reach
    sig = w / (2*sqrt(2*log(2)));
    h = P.h;
    need = max(abs(U(:))) + 2*h;
    k = [w, P.tau, P.Om, h];
    if isempty(key) || ~isequal(key, k) || reach < need
        nU = ceil(max(1.3*need, 5)/h);
        kh = max(ceil(5*sig/h), 1);
        ue = (-(nU+kh):(nU+kh))' * h;
        Pe = rabiP(ue, P);
        g  = exp(-0.5*(((-kh:kh)')*h/sig).^2);
        g  = g / sum(g);
        nf = 2^nextpow2(numel(Pe) + numel(g) - 1);
        c  = real(ifft(fft(Pe, nf) .* fft(g, nf)));
        same = c(kh+1 : kh+numel(Pe));
        Sg = same(kh+1 : kh+2*nU+1);
        ug = (-nU:nU)' * h;
        pp = spline(ug, Sg);
        coefs = pp.coefs; u0 = ug(1);
        key = k; reach = nU*h;
    end
    % The knots are uniform, so the spline piece is found by arithmetic instead of
    % ppval's search -- this is most of the run time of a pulse-shape fit.
    t = (U(:) - u0) / h;
    i = min(max(floor(t), 0), size(coefs, 1) - 1);
    t = (t - i) * h;
    i = i + 1;
    S = reshape(((coefs(i,1).*t + coefs(i,2)).*t + coefs(i,3)).*t + coefs(i,4), size(U));
end

function v = rabiP(u, P)
    W2 = P.Om^2 + u.^2;
    v = P.Om^2 ./ W2 .* sin(pi*P.tau*sqrt(W2)).^2;
end

function P = pulseParams(tau, rotDeg)
    P = struct('tau', NaN, 'Om', NaN, 'h', NaN, 'wMin', NaN);
    if isempty(tau) || ~isfinite(tau) || tau <= 0; return; end
    P.tau = tau;
    P.Om = rotDeg/360/tau;                 % MHz: rotation = 2 pi Om tau
    P.h = min(0.01, 1/(40*tau));           % MHz, table step
    P.wMin = 3*P.h*2*sqrt(2*log(2));       % keeps the kernel >= 7 samples wide
end

function shapes = shapeList(shape, tau)
    switch shape
        case 'auto'
            shapes = {'lorentz', 'gauss'};
            if ~isempty(tau); shapes{end+1} = 'pulse'; end
        case {'lorentz', 'gauss', 'voigt'}
            shapes = {shape};
        case 'pulse'
            if isempty(tau)
                error(['esrHyperfineFit: shape ''pulse'' needs the pulse length. Pass tau ' ...
                       '(us), or a pulsed experiment or save it can be read from.']);
            end
            shapes = {'pulse'};
        otherwise
            error('esrHyperfineFit: shape must be auto, lorentz, gauss, voigt or pulse.');
    end
end

function [fw, pk] = effWidth(w, eta, shp, P)
    if ~strcmp(shp, 'pulse')
        fw = w; pk = 1; return
    end
    span = max(4*w, 4/P.tau) + 1;
    uu = linspace(0, span, 40001)';
    s = pulseShape(uu, w, P);
    pk = s(1);
    i = find(s < pk/2, 1);
    if isempty(i) || i < 2; fw = NaN; return; end
    fw = 2*interp1(s([i-1 i]), uu([i-1 i]), pk/2);
end

function v = curveAt(ff, F, D, opts, P)
    sz = size(ff);
    xx = ff(:);
    X = [baseCols(xx, D.xr, D.nB), depthCols(xx, F.theta, F, opts, P)];
    v = reshape(X * F.beta, sz);
end

function s = paramStruct(F, D, opts)
    s = struct('baseline', F.beta(1:D.nB)', 'depths', F.beta(D.nB+1:end)');
    [f0, w, A, eta] = unpackTheta(F.theta, F.names, opts);
    s.f0 = f0; s.w = w; s.A = A; s.eta = eta;
end

% ======================================================================
% DATA
% ======================================================================

function [x, y, e, tau, src] = getData(opts)
    tau = [];
    if ~isempty(opts.file)
        if exist(opts.file, 'file') ~= 2
            error('esrHyperfineFit: no such file: %s', opts.file);
        end
        S = load(opts.file);
        s = findSaved(S);
        if isempty(s)
            error('esrHyperfineFit: no experiment with summaryData in %s', opts.file);
        end
        [x, y, e] = fromSummary(s.summaryData);
        tau = pulseFrom(s); src = 'file';
        return
    end
    X = opts.exp;
    if ~isempty(X)
        if isstruct(X)
            s = findSaved(X);
            if ~isempty(s)
                [x, y, e] = fromSummary(s.summaryData);
                tau = pulseFrom(s); src = 'saved struct';
            else
                [x, y, e] = fromVectorStruct(X); src = 'struct';
            end
            return
        end
        if isnumeric(X)
            x = X(:, 1); y = X(:, 2);
            if size(X, 2) >= 3; e = X(:, 3); else; e = nan(size(x)); end
            src = 'matrix'; return
        end
        [x, y, e] = esrPullData(X, struct('avgOutlierMad', opts.avgOutlierMad, ...
                                          'verbose', opts.verbose, 'name', 'esrHyperfineFit'));
        tau = pulseFrom(X); src = 'experiment';
        return
    end
    if isempty(opts.f) || isempty(opts.y)
        error(['nothing to fit: give exp, file, or f and y. ' ...
               'esrHyperfineFit(struct(''exp'', p)).']);
    end
    x = opts.f(:); y = opts.y(:); e = opts.e(:);
    if isempty(e); e = nan(size(x)); end
    if numel(y) ~= numel(x) || numel(e) ~= numel(x)
        error('esrHyperfineFit: f, y and e have different lengths (%d, %d, %d).', ...
              numel(x), numel(y), numel(e));
    end
    src = 'vectors';
end

function s = findSaved(S)
% The experiment inside a SaveLoad file: load(file) gives myStruct.<ExpName>.summaryData.
    s = [];
    if isfield(S, 'summaryData'); s = S; return; end
    if isfield(S, 'myStruct') && isstruct(S.myStruct)
        S = S.myStruct;
        if isfield(S, 'summaryData'); s = S; return; end
    end
    fn = fieldnames(S);
    for i = 1:numel(fn)
        v = S.(fn{i});
        if isstruct(v) && isscalar(v) && isfield(v, 'summaryData')
            s = v; return
        end
    end
end

function [x, y, e] = fromSummary(sd)
    x = sd.x(:); y = sd.y(:);
    if isfield(sd, 'err') && numel(sd.err) == numel(x); e = sd.err(:);
    else; e = nan(size(x)); end
end

function [x, y, e] = fromVectorStruct(s)
    x = pick(s, {'f', 'freqs', 'freq', 'frequency', 'x'});
    y = pick(s, {'y', 'signal', 'value', 'values'});
    e = pick(s, {'e', 'err', 'sterr', 'error', 'errors'});
    if isempty(x) || isempty(y)
        error('esrHyperfineFit: this struct has no recognisable sweep (f/x and y/signal).');
    end
    if isempty(e); e = nan(size(x)); end
end

function v = pick(s, names)
    v = [];
    for k = 1:numel(names)
        if isfield(s, names{k}) && isnumeric(s.(names{k})) && ~isempty(s.(names{k}))
            v = s.(names{k})(:); return
        end
    end
end

function tau = pulseFrom(s)
% Pulse length of a pulsed ESR (AWGExpESRPulsed). Deliberately NOT piTime: CW classes can
% carry a piTime too, and a CW sweep must not be offered the pulse shape.
    tau = [];
    nm = {'pulseDurationActual', 'mwPulseDuration'};
    for k = 1:numel(nm)
        try; v = s.(nm{k}); catch; v = []; end
        if isnumeric(v) && isscalar(v) && isfinite(v) && v > 0
            tau = double(v); return
        end
    end
end

% ======================================================================
% REPORTING
% ======================================================================

function out = emptyOut(opts)
    out = struct('ok', false, 'why', '', 'nDips', NaN, 'isotope', '', ...
        'confidence', '', 'verdict', '', 'warnings', {{}}, ...
        'f0', NaN, 'f0Err', NaN, 'splitting', NaN, 'splittingErr', NaN, ...
        'fwhm', NaN, 'fwhmErr', NaN, 'widthParam', NaN, 'widthParamErr', NaN, ...
        'depth', NaN, 'depthErr', NaN, 'shape', '', 'eta', NaN, 'tau', [], ...
        'lines', emptyLines(), 'centreCov', [], 'candidates', [], 'runnerUp', [], ...
        'prob', nan(1, 4), 'depthTest', struct('pEqual', NaN, 'used', ''), ...
        'model', [], 'residual', [], 'modelAt', [], 'formula', '', 'params', [], ...
        'cov', [], 'chi2Reduced', NaN, 'birge', NaN, 'baseline', [], 'xref', NaN, ...
        'f', [], 'y', [], 'e', [], 'source', '', 'opts', opts);
end

function L = emptyLines()
    L = struct('centre', {}, 'centreErr', {}, 'fwhm', {}, 'fwhmErr', {}, ...
               'depth', {}, 'depthErr', {});
end

function out = copyFields(out, R)
    fn = fieldnames(R);
    for i = 1:numel(fn)
        if ~any(strcmp(fn{i}, {'ok', 'why'})); out.(fn{i}) = R.(fn{i}); end
    end
end

function out = failed(out, opts)
    if opts.verbose; fprintf('\nesrHyperfineFit FAILED: %s\n', out.why); end
end

function t = verdictText(o, gap, A14N, A15N)
    switch o.nDips
        case 0
            t = 'no line in this sweep';
        case 1
            t = sprintf('1 dip (hyperfine not resolved), FWHM %.2f MHz', o.fwhm);
        otherwise
            nom = ternary(o.nDips == 3, A14N, A15N);
            lab = ternary(o.nDips == 3, '14N triplet', '15N doublet');
            if strcmp(o.isotope, 'unknown')
                lab = sprintf('NOT a %s (expected %.2f MHz)', lab, nom);
            end
            t = sprintf('%d dips: %s, A = %.3f +- %.3f MHz', o.nDips, lab, ...
                        o.splitting, o.splittingErr);
    end
    if isfinite(gap) && ~strcmp(o.confidence, 'forced')
        t = sprintf('%s  [%s: dBIC %.1f over %s]', t, o.confidence, gap, ...
                    dipName(o.runnerUp.nDips));
    elseif strcmp(o.confidence, 'forced')
        t = sprintf('%s  [forced]', t);
    end
end

function s = dipName(N)
    switch N
        case 0; s = 'no line';
        case 1; s = '1 dip';
        otherwise; s = sprintf('%d dips', N);
    end
end

function W = collectWarnings(o, F, D, opts, P, A14N, A15N)
    W = {};
    if F.nDips == 0; return; end
    nm = F.names;
    for j = 1:numel(nm)
        if o.atBound(j)
            switch nm{j}
                case 'A'
                    W{end+1} = sprintf(['the splitting sits on its bound (%.2f MHz): widen ' ...
                        'splitRange, or the pattern is not hyperfine'], F.theta(j));     %#ok<AGROW>
                case 'w'
                    if F.theta(j) <= F.lo(j)*1.0001
                        if strcmp(F.shape, 'pulse')
                            W{end+1} = ['the inhomogeneous width went to ~0: the pulse alone ' ...
                                'sets the line width'];                                   %#ok<AGROW>
                        else
                            W{end+1} = sprintf(['the width hit its floor (%.2f MHz, half the ' ...
                                'frequency step): the lines are not resolved by this grid'], ...
                                F.theta(j));                                              %#ok<AGROW>
                        end
                    else
                        W{end+1} = 'the width hit its ceiling: no resolved line';         %#ok<AGROW>
                    end
                case 'f0'
                    W{end+1} = 'the pattern centre sits on the edge of the sweep';        %#ok<AGROW>
            end
        end
    end
    c = [o.lines.centre];
    if any(c < min(D.x)) || any(c > max(D.x))
        W{end+1} = 'a line centre lies outside the sweep: extend it for a reliable count'; %#ok<AGROW>
    end
    if o.depth < 3*o.depthErr
        W{end+1} = sprintf('the dips are only %.1f sigma deep', o.depth/max(o.depthErr, eps)); %#ok<AGROW>
    end
    if strcmp(o.depthTest.used, 'free')
        W{end+1} = sprintf(['the depths are NOT equal (p = %.2g): nitrogen polarisation, ' ...
            'or the drive differs across the lines'], o.depthTest.pEqual);             %#ok<AGROW>
        if any([o.lines.depth] <= 0)
            W{end+1} = ['a line came out with zero depth: the pattern is not a full ' ...
                'multiplet, check the count'];                                         %#ok<AGROW>
        end
    end
    if F.nDips >= 2 && strcmp(o.isotope, 'unknown')
        W{end+1} = sprintf(['spacing %.3f MHz matches neither 14N (%.2f) nor 15N (%.2f); ' ...
            'two orientations or sites a few MHz apart look like this'], ...
            o.splitting, A14N, A15N);                                                  %#ok<AGROW>
    end
    if ~strcmp(F.shape, 'pulse') && isfinite(P.tau) && F.nDips > 0
        % nothing: the pulse shape was offered and lost, which is a legitimate outcome
    end
    if 2*D.dx > o.fwhm && F.nDips > 0
        W{end+1} = sprintf(['the frequency step (%.2f MHz) is more than half the line FWHM ' ...
            '(%.2f): the width is barely sampled, use a finer step for line shapes'], ...
            D.dx, o.fwhm);                                                             %#ok<AGROW>
    end
    if strcmp(o.confidence, 'ambiguous')
        W{end+1} = 'the line count is ambiguous: more averages, or a finer/longer sweep'; %#ok<AGROW>
    end
end

function report(o, A14N, A15N)
    fprintf('\n--- esrHyperfineFit: %s ---\n', o.verdict);
    fprintf('  %-10s %-8s %9s %9s %7s\n', 'candidate', 'shape', 'chi2_nu', 'dBIC', 'P');
    c = o.candidates;
    for i = 1:numel(c)
        mark = ''; if c(i).nDips == o.nDips; mark = '  <--'; end
        fprintf('  %-10s %-8s %9.2f %9.1f %7.3f%s\n', dipName(c(i).nDips), ...
                ternary(c(i).nDips == 0, '-', c(i).shape), c(i).chi2Reduced, ...
                c(i).dBic, c(i).prob, mark);
    end
    if o.nDips == 0
        fprintf('  reduced chi2 %.2f\n', o.chi2Reduced); return
    end
    lab = ternary(o.nDips == 3, '(the mI = 0 line)', ternary(o.nDips == 2, '(midpoint)', ''));
    fprintf('  f0          %10.3f +- %.3f MHz  %s\n', o.f0, o.f0Err, lab);
    if o.nDips >= 2
        nom = ternary(o.nDips == 3, A14N, A15N);
        fprintf('  splitting   %10.3f +- %.3f MHz  (%s %.2f: %+.1f sigma)\n', ...
                o.splitting, o.splittingErr, ternary(o.nDips == 3, '14N', '15N'), nom, ...
                (o.splitting - nom)/max(o.splittingErr, eps));
        fprintf('  lines       %s MHz\n', numList([o.lines.centre], '%.3f', '  '));
    end
    if strcmp(o.shape, 'pulse')
        fprintf('  FWHM        %10.3f +- %.3f MHz  (pulse %.3f us, inhomogeneous FWHM %.3f +- %.3f)\n', ...
                o.fwhm, o.fwhmErr, o.tau, o.widthParam, o.widthParamErr);
    else
        fprintf('  FWHM        %10.3f +- %.3f MHz  (%s)\n', o.fwhm, o.fwhmErr, o.shape);
    end
    if strcmp(o.depthTest.used, 'free')
        fprintf('  depths      %s %%\n', numList(100*[o.lines.depth], '%.3f', '  '));
    else
        fprintf('  depth       %10.3f +- %.3f %% per line\n', 100*o.depth, 100*o.depthErr);
    end
    if isfinite(o.depthTest.pEqual)
        fprintf('  equal depths: p = %.2g -> %s\n', o.depthTest.pEqual, o.depthTest.used);
    end
    fprintf('  reduced chi2 %.2f (errors scaled %.2fx)\n', o.chi2Reduced, o.birge);
    fprintf('  %s\n', o.formula);
    for i = 1:numel(o.warnings)
        fprintf('  NOTE: %s\n', o.warnings{i});
    end
end

function s = baselineText(b, xr)
    s = sprintf('%.6f', b(1));
    for j = 2:numel(b)
        s = sprintf('%s %+.3g*(f-%.3f)%s', s, b(j), xr, ternary(j > 2, sprintf('^%d', j-1), ''));
    end
end

function s = modelText(R, F)
    switch F.shape
        case 'lorentz'; sh = 'L';
        case 'gauss';   sh = 'G';
        case 'voigt';   sh = sprintf('V[eta %.2f]', R.eta);
        case 'pulse';   sh = 'Pulse*G';
    end
    mStr = numList(F.m, '%g', ',');
    C = R.params(numel(F.theta) + numel(R.baseline) + 1:end);
    if numel(C) == 1; cs = sprintf('%.4g', C); else; cs = ['[' numList(C, '%.4g', ' ') ']_m']; end
    s = sprintf('y(f) = %s - %s * SUM_{m=%s} %s(f - %.3f - m*%.3f; w = %.3f)', ...
                baselineText(R.baseline, R.xref), cs, mStr, sh, R.f0, ...
                ternary(F.nDips >= 2, R.splitting, 0), R.widthParam);
end

function drawIt(o)
    haveStyle = exist('fitStyle', 'file') == 2;
    ff = linspace(min(o.f), max(o.f), 3000)';
    ttl = titleLines(o);
    if haveStyle
        S = fitStyle();
        fig = S.newFigure('ESR hyperfine', o.opts.showResid);
        if o.opts.showResid; ax = S.axes(subplot(4, 1, 1:3)); else; ax = S.axes(fig); end
        H = []; L = {};
        H(end+1) = S.data(ax, o.f, o.y, fillErr(o.e), S.col.s1); L{end+1} = 'measured';
        if ~isempty(o.runnerUp) && o.runnerUp.nDips > 0
            H(end+1) = plot(ax, ff, o.runnerUp.modelAt(ff), '--', 'LineWidth', 1.3, ...
                            'Color', S.col.s3);
            L{end+1} = sprintf('%s (dBIC %.1f)', dipName(o.runnerUp.nDips), o.runnerUp.dBic);
        end
        H(end+1) = S.model(ax, ff, o.modelAt(ff)); L{end+1} = sprintf('fit: %s', dipName(o.nDips));
        S.tight(ax, o.f); S.headroom(ax, 0.30);
        yl = ylim(ax); top = yl(1) + 0.05*(yl(2) - yl(1));
        for i = 1:numel(o.lines)
            plot(ax, [1 1]*o.lines(i).centre, [yl(1) top], '-', 'LineWidth', 1.6, ...
                 'Color', S.col.flag);
        end
        ylim(ax, yl);
        xlabel(ax, 'frequency  [MHz]', 'Color', S.col.ink2);
        ylab = S.yLabel; if ~isempty(o.opts.yLabel); ylab = o.opts.yLabel; end
        ylabel(ax, ylab, 'Color', S.col.ink2);
        S.legend(ax, H, L);
        S.title(ax, ttl);
        if o.opts.showResid
            axR = S.axes(subplot(4, 1, 4));
            S.zeroLine(axR, [min(o.f) max(o.f)]);
            S.resid(axR, o.f, o.residual ./ max(fillErr(o.e)*o.birge, eps), S.col.s1);
            xlabel(axR, 'frequency  [MHz]', 'Color', S.col.ink2);
            ylabel(axR, 'residual / \sigma', 'Color', S.col.ink2);
            S.tight(axR, o.f);
        end
    else
        figure('Name', 'ESR hyperfine'); hold on; box on;
        errorbar(o.f, o.y, fillErr(o.e), 'o-');
        if ~isempty(o.runnerUp) && o.runnerUp.nDips > 0
            plot(ff, o.runnerUp.modelAt(ff), '--');
        end
        plot(ff, o.modelAt(ff), 'k-', 'LineWidth', 2);
        xlabel('frequency [MHz]'); ylabel('normalised fluorescence');
        title(ttl);
    end
end

function e = fillErr(e)
    if isempty(e) || all(~isfinite(e)); e = zeros(size(e)); return; end
    e(~isfinite(e)) = median(e(isfinite(e)));
end

function T = titleLines(o)
% Two lines, like esrFitRobust: a third is clipped above a residual panel.
    if o.nDips == 0
        T = {'no line found', sprintf('\\chi^2_\\nu = %.2f', o.chi2Reduced)}; return
    end
    if o.nDips >= 2
        a = sprintf(['%d dips (%s):   A = %.3f \\pm %.3f MHz     f_0 = %.3f \\pm %.3f MHz' ...
                     '     FWHM = %.2f \\pm %.2f MHz (%s)'], o.nDips, o.isotope, ...
                    o.splitting, o.splittingErr, o.f0, o.f0Err, o.fwhm, o.fwhmErr, o.shape);
    else
        a = sprintf('1 dip (unresolved):   f_0 = %.3f \\pm %.3f MHz     FWHM = %.2f \\pm %.2f MHz (%s)', ...
                    o.f0, o.f0Err, o.fwhm, o.fwhmErr, o.shape);
    end
    b = sprintf(['depth %.3f \\pm %.3f %% per line     |     P(%d dips) = %.3f, %s     |     ' ...
                 '\\chi^2_\\nu = %.2f     |     %d points'], 100*o.depth, 100*o.depthErr, ...
                o.nDips, o.prob(o.nDips + 1), o.confidence, o.chi2Reduced, numel(o.f));
    T = {a, b};
end

% ======================================================================
% SMALL HELPERS
% ======================================================================

function p = chi2sf(x, k)
% Upper tail of chi-square with k dof, base MATLAB (no Statistics Toolbox).
    if k <= 0; p = NaN; return; end
    p = gammainc(x/2, k/2, 'upper');
end

function s = numList(v, fmt, sep)
    parts = cell(1, numel(v));
    for i = 1:numel(v); parts{i} = sprintf(fmt, v(i)); end
    s = strjoin(parts, sep);
end

function v = ternary(c, a, b)
    if c; v = a; else; v = b; end
end

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
    known = fieldnames(d); given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('%s: unknown option "%s". Valid: %s.', who, bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || isempty(o.(known{i}))
            o.(known{i}) = d.(known{i});
        end
    end
end
