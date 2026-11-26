function out = coilMatrixFit(data, opts)
%COILMATRIXFIT The 3x3 coil matrix, from ESR lines measured at known coil currents.
%
%   out = coilMatrixFit(data)
%   out = coilMatrixFit(data, struct('gaugeRef', h.coilMatrix, 'nominalM', Mnom))
%
% THE MODEL
%
%       B = M * I + b0                 B in Gauss, I in Ampere, M in G/A, b0 in Gauss
%
% with B in the NV CRYSTAL FRAME: x, y, z along the diamond's cubic <100> axes, so the four
% NV axes are (+-1, +-1, +-1)/sqrt(3) and an NV-aligned field has theta = acosd(1/sqrt(3))
% = 54.7356 deg and phi = 45 deg (+ k*90). M is a full 3x3: column k is the field coil k
% makes per Ampere, so the off-diagonal terms are exactly the cross-talk between coil axes
% and crystal axes -- a turned diamond, tilted coils, coils not quite perpendicular, all in
% one matrix. b0 is the field at zero current (the earth, the lab, anything magnetised
% nearby). Nothing about where the diamond sits is entered by hand: the NVs define the frame.
%
% Given (M, b0, D) every spectrum's eight lines follow from the exact spin-1 Hamiltonian per
% orientation, so the fit is 13 numbers against every line of every spectrum at once.
%
% ===================================================================
% HOW, AND WHY IN THIS ORDER
% ===================================================================
%
% The final answer is a direct fit to the line centres. That cost has many local minima --
% the right field with the crystal turned 90 degrees fits three of four orientations -- so
% it needs a start that is already right. The start is built from steps that cannot get
% stuck:
%
%   1. EVERY SPECTRUM'S OWN FIELD, BY EXHAUSTIVE SEARCH. The lines only see |B . n_k|, and
%      every signed permutation of (x,y,z) maps the NV axes onto themselves, so B only needs
%      searching over the 1/48 of directions with x >= y >= z >= 0. That wedge is gridded
%      at ~1.5 deg and |B| at ~1%, the eight lines of every grid point come from the
%      closed-form cubic (no eig in the loop), and the best distinct candidates are
%      polished. esrFieldSolve is NOT used here: it seeds within 25 deg of one NV axis,
%      which is right for an alignment and wrong for a calibration, whose fields point
%      anywhere -- tested, it converged to a wrong field on 3 of 16 generic spectra even
%      with exact line positions. A spectrum whose second-best DISTINCT field fits almost
%      as well is marked ambiguous and kept out of the start (it still counts in step 4).
%
%   2. |B| IS LINEAR IN THE UNKNOWNS. |B|^2 = I'*(M'*M)*I + 2*(M'*b0)'*I + |b0|^2, so the
%      ten numbers G = M'*M, h = M'*b0, c = |b0|^2 come from an ordinary (robustly
%      reweighted) least-squares fit, with no start point. G fixes every coil's
%      sensitivity and the angles between the coil fields; h fixes b0 against the coils.
%      Only the crystal's ORIENTATION is left: M = Q*chol(G).
%
%   3. THE ORIENTATION IS THREE NUMBERS, SO IT CAN BE SEARCHED. Thousands of random
%      rotations are scored against every spectrum's sorted |projections| (sorted, so the
%      score does not care which line is which), and the best few are polished.
%
%   4. THE REAL FIT, from each of the best few starts: Levenberg-Marquardt on all 13
%      parameters against the line centres, each measured line assigned to its nearest
%      predicted transition, assignment redone until it stops changing. A measured line
%      further than huberMHz from every prediction is an outlier (the resonator dip, a bad
%      refit), reported and not fitted. BLENDS: a predicted line that no measured line
%      claims, sitting within ~blendMHz of one that is claimed, is averaged into that
%      line's prediction (weight 1/(1+(dist/blendMHz)^4)), because a line fit to two
%      overlapping transitions returns their centroid, not either one.
%
% THE 48 EQUIVALENT ANSWERS. (g*M, g*b0) fits exactly as well as (M, b0) for every signed
% permutation g. That is the freedom to NAME the crystal axes, and it is fixed by choosing
% the g that puts the NV frame closest to gaugeRef: the previous calibration when there is
% one (so "the [1 1 1] NV" keeps meaning the same physical NV across recalibrations -- B1
% depends on which one it is), otherwise the nominal coil matrix (so NV x ~ coil x, etc.).
%
% Base MATLAB only: no Optimization Toolbox on this setup, so fminsearch for the small
% searches and a hand-written Levenberg-Marquardt for the final fit.
%
% ===================================================================
% INPUT
% ===================================================================
%   data       struct array, one element per spectrum:
%                I          1x3 currents, A (as SET -- the supplies have no readback)
%                centres    fitted NV line centres, MHz (esrFitRobust out.lines.centre)
%                centreErr  1-sigma, MHz (optional; floored at centreFloorMHz)
%                sweep      [fmin fmax] MHz (optional; a predicted line outside it is not
%                           expected to be seen)
%
% OPTIONS
%   gaugeRef        3x3: the answer's frame is put closest to this (default nominalM)
%   nominalM        3x3 G/A, the model in use before calibrating: gauge fallback, fallback
%                   start, and the comparison printout
%   nominalB0       1x3 G
%   D0, dTol        MHz (2870, 15)
%   perSpectrumD    fit D separately per spectrum (default false). Only if a long run drifts.
%   huberMHz        outlier distance (default 5)
%   blendMHz        blend scale (default 4, about half a narrow linewidth)
%   centreFloorMHz  added in quadrature to every line's error (default 0.1). esrFitRobust
%                   errors can be tiny on deep lines; one 0.01 MHz line must not run the fit.
%   nRotations      random rotations in step 3 (default 6000)
%   nStarts         distinct starts taken into step 4 (default 3)
%   seed            for steps 3 and 3b (default 7; private RandStreams, the global RNG is
%                   untouched, so the same data always gives the same answer)
%   ransacSubsets   four-spectrum draws in step 3b (default 30)
%   ransacTolG      a spectrum agrees with a candidate if its field is this close, Gauss
%                   (default 0.5)
%   verbose         print the report (default true)
%
% ===================================================================
% OUTPUT
% ===================================================================
%   ok, why
%   M, b0, D             the calibration (D is a vector if perSpectrumD)
%   Merr, b0err, Derr    1-sigma, Birge-scaled (two-sided, as esrFitRobust)
%   cov                  of [M(:); b0(:); D(:)]
%   rmsMHz, chi2Red, nLines, nOutliers, nSpectra, nUsed
%   spectra(j)           I, B (model), Bmag, BmagSolve (its own solve), ambiguous, centres,
%                        predicted, resid, outlier, rms, allLines
%   coilSensitivity      |M(:,k)|, G/A -- compare magnetic_field_ratio in setupInfo
%   coilDirections       M(:,k)/|M(:,k)|, in the NV frame
%   coilThetaPhi         3x2, each coil axis as (theta, phi) in the NV frame, degrees
%   coilAnglesDeg        3x3, angles between the coil fields (90 = orthogonal)
%   gauge                the signed permutation applied
%   init                 how the start was made, and how well it fitted
%   calibration          the struct Helmholtz.setCoilCalibration takes and the json keeps
%   fieldAt              @(I) B, for 1x3 or Nx3 currents
%   linesAt              @(I) the 8 transitions at those currents (mean D)

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('gaugeRef', [], 'nominalM', [], 'nominalB0', [], ...
        'D0', 2870, 'dTol', 15, 'perSpectrumD', false, 'huberMHz', 5, 'blendMHz', 4, ...
        'centreFloorMHz', 0.1, 'nRotations', 6000, 'nStarts', 3, 'seed', 7, ...
        'ransacSubsets', 30, 'ransacTolG', 0.5, ...
        'verbose', true));

    out = struct('ok', false, 'why', '');
    AX = nvAxes();

    % ---- tidy the input -------------------------------------------------------
    data = data(:)';
    N = numel(data);
    I = zeros(N, 3);
    hasSweep = isfield(data, 'sweep');
    for j = 1:N
        I(j, :) = reshape(data(j).I, 1, 3);
        c = data(j).centres(:)';
        if isfield(data, 'centreErr') && numel(data(j).centreErr) == numel(c)
            ce = data(j).centreErr(:)';
        else
            ce = 0.5*ones(size(c));
        end
        ce(~isfinite(ce) | ce <= 0) = 0.5;
        keep = isfinite(c);
        data(j).centres = c(keep);
        data(j).centreErr = sqrt(ce(keep).^2 + opts.centreFloorMHz^2);
        if ~hasSweep || numel(data(j).sweep) ~= 2
            data(j).sweep = [-inf inf];
        end
    end
    out.nSpectra = N;
    if N < 5
        out.why = sprintf(['%d spectra. The model has 13 parameters and a spectrum pins about ' ...
            '4, so at least 5 are needed -- 12 or more to be any good.'], N);
        if opts.verbose; fprintf('\ncoilMatrixFit: %s\n', out.why); end
        return
    end

    % ---- 1. every spectrum's own field, by exhaustive search -------------------
    dirs = wedgeDirections(1.5);
    Bmag = NaN(N, 1); P = NaN(4, N); Dj = NaN(N, 1); Bhat = NaN(N, 3);
    solved = false(N, 1); ambiguous = false(N, 1);
    for j = 1:N
        s = solveSpectrum(data(j), dirs, opts);
        if s.ok
            Bmag(j) = norm(s.Bvec);
            Bhat(j, :) = canon(s.Bvec);
            P(:, j) = sort(abs(AX * s.Bvec(:)), 'descend');
            Dj(j) = s.D;
            solved(j) = true;
            ambiguous(j) = s.ambiguous;
        end
    end
    usable = solved & ~ambiguous;
    out.nUsed = nnz(usable);

    % ---- 2. M'*M, M'*b0, |b0|^2 from |B|^2: robust linear least squares ---------
    init = struct('source', '', 'G', [], 'h', [], 'c', NaN, 'cCheck', NaN, ...
                  'rotRmsG', NaN, 'nGood', 0, 'downweighted', [], 'startUsed', '', ...
                  'nStarts', 0);
    metrics = struct('L', {}, 'h', {}, 'w', {}, 'label', {});
    wBase = double(solved(:)'); wBase(ambiguous(:)') = 0.25;
    if nnz(usable) >= 11
        uIdx = find(usable);
        [G, h, cc, inl] = robustMetricFit(I(usable, :), Bmag(usable));
        [L, notPD] = chol(G);
        if ~notPD && rcond(G) > 1e-8
            init.source = 'linear |B|^2 fit (trimmed)';
            init.G = G; init.h = h; init.c = cc;
            init.cCheck = h' * (G \ h);
            init.downweighted = uIdx(~inl)';
            w = wBase; w(uIdx(~inl)) = 0.1;
            metrics(end+1) = struct('L', L, 'h', h, 'w', w, 'label', 'data |B|^2 metric');
        end
    end
    % The nominal coil model is ALWAYS a second source of starts. It is what saves the day
    % when too few spectra are unambiguous for the |B|^2 fit, and costs a couple of seconds
    % when it is not needed -- step 4 picks whichever start fits the lines best.
    Mn = opts.nominalM;
    if isempty(Mn); Mn = eye(3) * nanMedian(Bmag ./ max(sqrt(sum(I.^2, 2)), 0.1)); end
    b0n = opts.nominalB0; if isempty(b0n); b0n = [0 0 0]; end
    if all(isfinite(Mn(:))) && rcond(Mn' * Mn) > 1e-10
        metrics(end+1) = struct('L', chol(Mn' * Mn), 'h', Mn' * b0n(:), 'w', wBase, ...
                                'label', 'nominal metric');
    end
    if isempty(init.source)
        init.source = 'nominal coil model only (too few unambiguous spectra for the |B|^2 fit)';
    end

    % ---- 3. the crystal orientation: rotation search, per metric --------------
    % B_j = Q*u_j, u_j = L*I_j + L'\h. Proper rotations suffice: -eye(3) is a symmetry.
    starts = struct('R', {}, 'cost', {}, 'p0', {}, 'label', {}, 'nInliers', {});
    for k = 1:numel(metrics)
        starts = [starts, orientationStarts(metrics(k).L, metrics(k).h, I, P, ...
                  metrics(k).w, AX, opts.nRotations, opts.seed + k, opts.nStarts, ...
                  metrics(k).label)];                                         %#ok<AGROW>
    end
    % ---- 3b. RANSAC on the per-spectrum fields: needs no metric at all ------------
    % Steps 2-3 lean on the |B|^2 regression, and ten parameters from ~15 spectra cannot
    % tell one wrong spectrum from a good one (tested: with 3 to spare, ANY 11 of them fit
    % almost exactly, so trimming picked the wrong ones). This start uses the vectors
    % instead: four spectra fix M and b0 exactly once each one's symmetry label is chosen,
    % all 48^3 labellings are tried at once, and each is scored by how many OTHER spectra
    % it predicts -- compared through canon(), which no labelling changes.
    rsR = RandStream('mt19937ar', 'Seed', opts.seed + 100);
    rStarts = ransacStarts(I, Bhat, find(usable), find(solved), opts.ransacSubsets, ...
                           opts.ransacTolG, rsR, opts.nStarts);
    starts = [starts, rStarts];
    init.ransacInliers = [];
    if ~isempty(rStarts); init.ransacInliers = rStarts(1).nInliers; end
    if isempty(starts)
        out.why = 'fewer than 3 spectra have a field solution; the orientation cannot be found';
        if opts.verbose; fprintf('\ncoilMatrixFit: %s\n', out.why); end
        return
    end
    init.nGood = nnz(usable);
    D00 = nanMedian(Dj(usable)); if ~isfinite(D00); D00 = nanMedian(Dj); end
    if ~isfinite(D00); D00 = opts.D0; end

    % ---- 4. the real fit, to every line centre, from each start ------------------
    nD = 1; if opts.perSpectrumD; nD = N; end
    best = []; bestScore = inf;
    for s = 1:numel(starts)
        p0 = [starts(s).p0; D00 * ones(nD, 1)];
        [p, asg] = fitLines(p0, data, I, nD, opts);
        sc = robustScore(p, data, I, nD, opts);
        if sc < bestScore
            bestScore = sc;
            best = struct('p', p, 'asg', {asg}, 'start', s);
        end
    end
    p = best.p; asg = best.asg;
    init.rotRmsG = starts(best.start).cost;
    init.startUsed = sprintf('%s, start %d of %d', starts(best.start).label, best.start, ...
                             numel(starts));
    init.nStarts = numel(starts);
    out.init = init;

    r = residuals(p, data, I, nD, asg, opts);
    J = numJac(@(q) residuals(q, data, I, nD, asg, opts), p);
    nRes = numel(r); nPar = numel(p);
    if nRes <= nPar
        out.why = sprintf(['only %d lines were matched to a prediction, for %d parameters. ' ...
            'The spectra do not determine the coils.'], nRes, nPar);
        if opts.verbose; fprintf('\ncoilMatrixFit: %s\n', out.why); end
        return
    end
    chi2Red = sum(r.^2) / (nRes - nPar);
    JtJ = J' * J;
    if rcond(JtJ) < 1e-14; cov = pinv(JtJ); else; cov = JtJ \ eye(nPar); end
    cov = cov * max(chi2Red, 0.05^2);                    % two-sided Birge

    % ---- 5. name the crystal axes ----------------------------------------------
    M = reshape(p(1:9), 3, 3);
    ref = opts.gaugeRef;
    if isempty(ref); ref = opts.nominalM; end
    if isempty(ref); ref = eye(3); end
    g = pickGauge(M, ref);
    T = blkdiag(kron(eye(3), g), g, eye(nD));
    p = T * p;
    cov = T * cov * T';
    M = reshape(p(1:9), 3, 3);
    b0 = p(10:12);
    D = p(13:end);
    sd = sqrt(max(diag(cov), 0));

    out.M = M; out.b0 = b0(:)'; out.D = D(:)';
    out.Merr = reshape(sd(1:9), 3, 3); out.b0err = sd(10:12)'; out.Derr = sd(13:end)';
    out.cov = cov; out.gauge = g; out.chi2Red = chi2Red;

    % ---- per-spectrum table ------------------------------------------------------
    Sx = struct('I', {}, 'B', {}, 'Bmag', {}, 'BmagSolve', {}, 'ambiguous', {}, ...
                'centres', {}, 'predicted', {}, 'resid', {}, 'outlier', {}, 'rms', {}, ...
                'allLines', {});
    allRes = []; nOut = 0;
    asgFinal = assignLines(p, data, I, nD, opts.huberMHz);
    for j = 1:N
        Bj = M * I(j, :)' + b0;
        pred = nvLinesVec(Bj', D(min(j, nD)));
        [pr, ~] = blendPrediction(pred, asgFinal{j}, opts.blendMHz);
        c = data(j).centres;
        isOut = asgFinal{j} == 0;
        res = c - pr;
        Sx(j) = struct('I', I(j, :), 'B', Bj', 'Bmag', norm(Bj), 'BmagSolve', Bmag(j), ...
                       'ambiguous', ambiguous(j) || ~solved(j), 'centres', c, ...
                       'predicted', pr, 'resid', res, 'outlier', isOut, ...
                       'rms', sqrt(mean(res(~isOut).^2)), 'allLines', sort(pred));
        allRes = [allRes, res(~isOut)];                                       %#ok<AGROW>
        nOut = nOut + nnz(isOut);
    end
    out.spectra = Sx;
    out.rmsMHz = sqrt(mean(allRes.^2));
    out.nLines = numel(allRes);
    out.nOutliers = nOut;

    % ---- what the matrix says about the coils -------------------------------------
    sens = sqrt(sum(M.^2, 1));
    dirsC = M ./ sens;
    out.coilSensitivity = sens;
    out.coilDirections = dirsC;
    out.coilThetaPhi = [acosd(max(min(dirsC(3, :)', 1), -1)), atan2d(dirsC(2, :)', dirsC(1, :)')];
    out.coilAnglesDeg = acosd(max(min(dirsC' * dirsC, 1), -1));
    out.condM = cond(M);

    out.fieldAt = @(Iq) (M * reshape(Iq, [], 3)' + b0)';
    Dm = mean(D);
    out.linesAt = @(Iq) nvLinesVec((M * Iq(:) + b0)', Dm);

    out.calibration = struct( ...
        'date',       datestr(now, 'yyyy-mm-dd HH:MM'), ...                  %#ok<TNOW1,DATST>
        'frame',      'NV', ...
        'matrix',     M, ...
        'offset',     out.b0, ...
        'matrixErr',  out.Merr, ...
        'offsetErr',  out.b0err, ...
        'D',          Dm, ...
        'rmsMHz',     out.rmsMHz, ...
        'nSpectra',   N, ...
        'currentMin', min(I, [], 1), ...
        'currentMax', max(I, [], 1));

    out.ok = true;
    fracOut = out.nOutliers / max(out.nLines + out.nOutliers, 1);
    if fracOut > 0.25
        % Not noise: a quarter of the lines disagreeing with every prediction means the
        % model put the crystal the wrong way round. Do not hand this to the coils.
        out.ok = false;
        out.why = sprintf(['%d of %d lines are outliers -- more than a quarter. The fit did ' ...
            'not find a consistent crystal orientation. Look at the per-spectrum table: ' ...
            'too few lines per spectrum (lower the MW amplitude to narrow them), or the ' ...
            'sweep missing the outer lines, are the usual causes.'], out.nOutliers, ...
            out.nLines + out.nOutliers);
    elseif out.rmsMHz > 1.0
        out.why = sprintf(['line residual rms %.2f MHz is large for fitted centres; look ' ...
            'for one bad spectrum in the table before using this'], out.rmsMHz);
    end
    if opts.verbose; report(out, opts); end
end

% ======================================================================
%   the physics
% ======================================================================

function AX = nvAxes()
% The four NV axes in the crystal frame, in esrFieldSolve's order.
    AX = [ 1  1  1; 1 -1 -1; -1  1 -1; -1 -1  1] / sqrt(3);
end

function f = nvLinesVec(B, D)
% All eight transitions, MHz, for K fields at once (B is K x 3). Lines 2k-1, 2k belong to
% orientation k. Exact: the eigenvalues of H = D*Sz^2 + gamma*(Bpar*Sz + Bperp*Sx) from
% the closed-form cubic,
%     lambda^3 - 2D lambda^2 + (D^2 - z^2 - 2p^2) lambda + 2D p^2 = 0,
% z = gamma*Bpar, p = gamma*Bperp/sqrt(2). Checked against eig to 3e-9 MHz, including
% zero field and exact alignment (the degenerate cases).
    GAMMA = 2.8025;
    AX = nvAxes();
    proj = B * AX';                                   % K x 4
    B2 = sum(B.^2, 2);
    z = GAMMA * proj;
    p2 = GAMMA^2 * max(B2 - proj.^2, 0) / 2;
    A = -2*D;
    Bc = D.^2 - z.^2 - 2*p2;
    Cc = 2*D .* p2;
    Pq = Bc - A.^2/3;
    Qq = 2*A.^3/27 - A.*Bc/3 + Cc;
    m = 2*sqrt(max(-Pq/3, 0));
    arg = (3*Qq ./ (2*Pq)) .* sqrt(max(-3 ./ Pq, 0));
    arg(~isfinite(arg)) = 1;
    arg = max(min(arg, 1), -1);
    th = acos(arg) / 3;
    l1 = m .* cos(th) - A/3;
    l2 = m .* cos(th - 2*pi/3) - A/3;
    l3 = m .* cos(th - 4*pi/3) - A/3;
    lo = min(min(l1, l2), l3);
    hi = max(max(l1, l2), l3);
    mid = l1 + l2 + l3 - lo - hi;
    K = size(B, 1);
    f = zeros(K, 8);
    f(:, 1:2:end) = mid - lo;
    f(:, 2:2:end) = hi - lo;
end

function dirs = wedgeDirections(stepDeg)
% Unit vectors with x >= y >= z >= 0 -- one copy of every field direction up to the 48
% crystal symmetries -- at about stepDeg spacing (a Fibonacci sphere, cut to the wedge).
    n = ceil(4*pi / (stepDeg*pi/180)^2);
    k = (0:n-1)' + 0.5;
    zc = 1 - 2*k/n;
    r = sqrt(1 - zc.^2);
    ph = k * pi * (3 - sqrt(5));
    v = [r.*cos(ph), r.*sin(ph), zc];
    keep = v(:,1) >= v(:,2) - 1e-12 & v(:,2) >= v(:,3) - 1e-12 & v(:,3) >= -1e-12;
    dirs = [v(keep, :); 1 0 0; [1 1 0]/sqrt(2); [1 1 1]/sqrt(3)];
end

function s = solveSpectrum(d, dirs, opts)
% One spectrum's field, up to the crystal symmetries: grid, then polish distinct winners.
    GAMMA = 2.8025;
    s = struct('ok', false, 'Bvec', [NaN NaN NaN], 'D', NaN, 'cost', NaN, ...
               'altCost', NaN, 'ambiguous', true);
    c = d.centres; ce = max(d.centreErr, 0.2);
    if numel(c) < 3; return; end
    span = max(c) - min(c);
    Bspan = span / (2*GAMMA);
    if Bspan < 1; return; end
    % The widest pair belongs to the largest projection, and max|B.n| >= |B|/sqrt(3), so
    % |B| is within [Bspan, sqrt(3)*Bspan] when the outer pair is in the sweep. Wider when
    % it may not be.
    mags = Bspan * exp(linspace(log(0.85), log(2.2), 90));
    nd = size(dirs, 1);
    Bc = kron(mags(:), dirs);                         % (90*nd) x 3, magnitude-major
    cost = specCost(Bc, opts.D0, c, ce, d.sweep, opts.huberMHz);
    [~, ord] = sort(cost);
    fo = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-8, ...
                  'MaxFunEvals', 2000, 'MaxIter', 2000);
    sols = zeros(0, 4); scs = zeros(0, 1);
    for k = ord(1:min(200, numel(ord)))'
        cand = Bc(k, :);
        if ~isempty(sols) && any(sqrt(sum((canon(sols(:, 1:3)) - canon(cand)).^2, 2)) < 0.03*norm(cand))
            continue
        end
        f = @(q) specCost(q(1:3), q(4), c, ce, d.sweep, opts.huberMHz) + ...
                 1e3*max(abs(q(4) - opts.D0) - opts.dTol, 0)^2;
        [q, cv] = fminsearch(f, [cand, opts.D0], fo);
        sols(end+1, :) = q; scs(end+1) = cv;                                 %#ok<AGROW>
        if size(sols, 1) >= 6; break; end
    end
    [scs, o2] = sort(scs); sols = sols(o2, :);
    % the best and the best DISTINCT alternative
    alt = inf;
    for k = 2:numel(scs)
        if norm(canon(sols(k, 1:3)) - canon(sols(1, 1:3))) > max(0.5, 0.02*norm(sols(1, 1:3)))
            alt = scs(k); break
        end
    end
    s.ok = true;
    s.Bvec = sols(1, 1:3);
    s.D = sols(1, 4);
    s.cost = scs(1);
    s.altCost = alt;
    % ambiguous when the alternative is within ~3 sigma of one line in total cost
    s.ambiguous = (alt - scs(1)) * numel(c) < 9;
end

function J = specCost(B, D, c, ce, sweep, huber)
% esrFieldSolve's cost, vectorized over K candidate fields: each measured line to its
% nearest prediction (capped at huber), plus a charge for a predicted line inside the sweep
% with nothing measured near it -- without that, half the true |B| fits a clipped sweep.
    pred = nvLinesVec(B, D);                          % K x 8
    K = size(pred, 1);
    J = zeros(K, 1);
    cap = (huber ./ ce).^2;
    for l = 1:numel(c)
        dl = min(abs(pred - c(l)), [], 2);
        J = J + min((dl / ce(l)).^2, cap(l));
    end
    inSweep = pred >= sweep(1) + 2 & pred <= sweep(2) - 2;
    dmin = inf(K, 8);
    for l = 1:numel(c)
        dmin = min(dmin, abs(pred - c(l)));
    end
    miss = inSweep & dmin > huber;
    J = (J + 0.5 * median(cap) * sum(miss, 2)) / numel(c);
end

function v = canon(B)
% The representative of B in the x >= y >= z >= 0 wedge.
    v = sort(abs(B), 2, 'descend');
end

% ======================================================================
%   the global fit
% ======================================================================

function [M, b0, D] = unpack(p, nD)
    M = reshape(p(1:9), 3, 3);
    b0 = p(10:12);
    D = p(13:12+nD);
end

function asg = assignLines(p, data, I, nD, huber)
% For every measured line, the nearest predicted transition (1..8), or 0 if none is within
% huber MHz.
    [M, b0, D] = unpack(p, nD);
    asg = cell(1, numel(data));
    for j = 1:numel(data)
        pred = nvLinesVec((M * I(j, :)' + b0)', D(min(j, nD)));
        c = data(j).centres;
        a = zeros(size(c));
        for l = 1:numel(c)
            [dd, k] = min(abs(pred - c(l)));
            if dd <= huber; a(l) = k; end
        end
        asg{j} = a;
    end
end

function [pr, extraVar] = blendPrediction(pred, a, blendMHz)
% Each measured line's prediction: its assigned transition, averaged with any UNCLAIMED
% transition close to it (a single line fitted to two overlapping transitions lands near
% their centroid). extraVar inflates that line's variance by how unsure the blend is.
    pr = NaN(size(a)); extraVar = zeros(size(a));
    claimed = false(1, numel(pred));
    claimed(a(a > 0)) = true;
    for l = 1:numel(a)
        k = a(l);
        if k == 0; continue; end
        dd = pred - pred(k);
        w = 1 ./ (1 + (dd / blendMHz).^4);
        w(claimed) = 0; w(k) = 1;
        pr(l) = sum(w .* pred) / sum(w);
        extraVar(l) = sum(w .* (1 - w) .* dd.^2) / 4;
    end
end

function r = residuals(p, data, I, nD, asg, opts)
    [M, b0, D] = unpack(p, nD);
    r = [];
    for j = 1:numel(data)
        a = asg{j};
        if ~any(a); continue; end
        pred = nvLinesVec((M * I(j, :)' + b0)', D(min(j, nD)));
        [pr, ev] = blendPrediction(pred, a, opts.blendMHz);
        m = a > 0;
        sig = sqrt(data(j).centreErr(m).^2 + ev(m));
        r = [r, (data(j).centres(m) - pr(m)) ./ sig];                        %#ok<AGROW>
    end
    r = r(:);
end

function sc = robustScore(p, data, I, nD, opts)
% For choosing between starts: every measured line counts, an outlier at the cap.
    asg = assignLines(p, data, I, nD, inf);
    [M, b0, D] = unpack(p, nD);
    sc = 0;
    for j = 1:numel(data)
        pred = nvLinesVec((M * I(j, :)' + b0)', D(min(j, nD)));
        [pr, ~] = blendPrediction(pred, asg{j}, opts.blendMHz);
        dd = data(j).centres - pr;
        sc = sc + sum(min((dd ./ data(j).centreErr).^2, (opts.huberMHz ./ data(j).centreErr).^2));
    end
end

function [p, asg] = fitLines(p, data, I, nD, opts)
% Assign, fit, re-assign until stable. The outlier window starts wide and tightens, so a
% start that is a few MHz off can still pull its lines in.
    sched = opts.huberMHz * [3 2 1.5 1 1 1 1 1 1 1 1 1];
    asg = assignLines(p, data, I, nD, sched(1));
    for it = 1:numel(sched)
        p = levmar(@(q) residuals(q, data, I, nD, asg, opts), p);
        nxt = sched(min(it + 1, numel(sched)));
        asg2 = assignLines(p, data, I, nD, nxt);
        if it >= 4 && isequal(asg2, asg); break; end
        asg = asg2;
    end
end

function J = numJac(f, p)
    n = numel(p);
    r0 = f(p);
    J = zeros(numel(r0), n);
    for i = 1:n
        hstep = 1e-6 * max(1, abs(p(i)));
        e = zeros(n, 1); e(i) = hstep;
        J(:, i) = (f(p + e) - f(p - e)) / (2*hstep);
    end
end

function p = levmar(f, p)
% Levenberg-Marquardt with Marquardt scaling; 13 parameters and a few hundred residuals,
% so a dense solve per step is fine.
    r = f(p); cost = r' * r;
    lambda = 1e-3;
    for it = 1:200
        J = numJac(f, p);
        A = J' * J; gvec = J' * r;
        dA = diag(A); dA(dA <= 0) = 1;
        improved = false; rel = 0; step = zeros(size(p));
        for tries = 1:20
            step = -(A + lambda * diag(dA)) \ gvec;
            pn = p + step;
            rn = f(pn); cn = rn' * rn;
            if isfinite(cn) && cn < cost
                rel = (cost - cn) / max(cost, eps);
                p = pn; r = rn; cost = cn;
                lambda = max(lambda / 5, 1e-12);
                improved = true;
                break
            end
            lambda = lambda * 10;
        end
        if ~improved || rel < 1e-12 || norm(step) < 1e-12 * (norm(p) + 1e-12)
            break
        end
    end
end

function [G, h, c, inl, resG] = robustMetricFit(I, Bm)
% |B|^2 = I'*G*I + 2*h'*I + c, linear in (G, h, c) -- by EXHAUSTIVE least trimmed squares.
%
% Ten parameters from ~15-20 spectra leaves few spare equations, and one wrong per-spectrum
% field (a spurious line that happens to fit a different |B|) is several Gauss off. Tested:
% a single such spectrum among 15 dragged an ordinary fit far enough that bisquare
% reweighting, which starts FROM that fit, converged to the wrong G and threw out four good
% spectra instead. So every way of leaving out up to 3 spectra is tried (a few hundred
% 10x10 solves, milliseconds), each scored by the sum of its N-3 smallest squared
% residuals over ALL spectra, and the winner's inliers are refitted.
    N = numel(Bm);
    A = [I(:,1).^2, I(:,2).^2, I(:,3).^2, 2*I(:,1).*I(:,2), 2*I(:,1).*I(:,3), ...
         2*I(:,2).*I(:,3), 2*I, ones(N, 1)];
    y = Bm(:).^2;
    base = 1 ./ max(Bm(:), 1);                         % residuals in ~Gauss of |B|
    kMax = max(0, min(3, N - 11));
    hKeep = N - kMax;
    bestScore = inf; xBest = [];
    for k = 0:kMax
        if k == 0; combos = zeros(1, 0); else; combos = nchoosek(1:N, k); end
        for ci = 1:size(combos, 1)
            keep = true(N, 1); keep(combos(ci, :)) = false;
            x = (A(keep, :) .* base(keep)) \ (y(keep) .* base(keep));
            res = (A * x - y) .* base / 2;
            r2 = sort(res.^2);
            sc = sum(r2(1:hKeep));
            if sc < bestScore; bestScore = sc; xBest = x; end
        end
    end
    res = (A * xBest - y) .* base / 2;
    s = 1.4826 * sqrt(median(res.^2));
    inl = abs(res) <= max(4*s, 0.3);
    x = (A(inl, :) .* base(inl)) \ (y(inl) .* base(inl));
    resG = (A * x - y) .* base / 2;
    G = [x(1) x(4) x(5); x(4) x(2) x(6); x(5) x(6) x(3)];
    h = x(7:9);
    c = x(10);
end

function starts = orientationStarts(L, h, I, P, w, AX, nRot, seed, nKeep, label)
% Step 3 for one metric (L, h): score random rotations against the measured sorted
% |projections|, polish the best, keep the best nKeep DISTINCT orientations.
    N = size(I, 1);
    U = L * I' + repmat(L' \ h, 1, N);
    sel = find(all(isfinite(P), 1) & w(:)' > 0);
    starts = struct('R', {}, 'cost', {}, 'p0', {}, 'label', {}, 'nInliers', {});
    if numel(sel) < 3; return; end
    Ps = P(:, sel); Us = U(:, sel); ws = w(sel); ws = ws(:)';
    capG2 = 4;                    % a spectrum off by more than 2 G costs no more than 2 G
    rotCost = @(R) sum(ws .* min(sum((sort(abs(AX * R * Us), 1, 'descend') - Ps).^2, 1), capG2));
    rs = RandStream('mt19937ar', 'Seed', seed);
    qs = randn(rs, 4, nRot);
    qs = qs ./ sqrt(sum(qs.^2, 1));
    costs = zeros(1, nRot);
    for k = 1:nRot
        costs(k) = rotCost(quatToRot(qs(:, k)));
    end
    [~, ord] = sort(costs);
    fo = optimset('Display', 'off', 'TolX', 1e-9, 'TolFun', 1e-10, ...
                  'MaxFunEvals', 3000, 'MaxIter', 3000);
    cands = struct('R', {}, 'cost', {});
    for k = ord(1:min(16, numel(ord)))
        R0 = quatToRot(qs(:, k));
        [v, cv] = fminsearch(@(v) rotCost(rotvecToRot(v) * R0), [0 0 0], fo);
        cands(end+1) = struct('R', rotvecToRot(v) * R0, 'cost', cv);       %#ok<AGROW>
    end
    [~, o2] = sort([cands.cost]);
    cands = cands(o2);
    kept = cands(1);
    for k = 2:numel(cands)
        if numel(kept) >= nKeep; break; end
        isNew = true;
        for m = 1:numel(kept)
            if sameModuloCube(cands(k).R, kept(m).R); isNew = false; break; end
        end
        if isNew; kept(end+1) = cands(k); end                                  %#ok<AGROW>
    end
    for k = 1:numel(kept)
        R = kept(k).R;
        starts(end+1) = struct('R', R, 'cost', sqrt(kept(k).cost / max(4*sum(ws), 1)), ...
                               'p0', [reshape(R * L, [], 1); R * (L' \ h)], ...
                               'label', label, 'nInliers', NaN);            %#ok<AGROW>
    end
end

function starts = ransacStarts(I, Bhat, goodIdx, solvedIdx, nSubsets, tau, rs, nKeep)
% M and b0 from four spectra at a time, every symmetry labelling at once.
%
% With spectrum r's label fixed (that fixes the gauge) and labels g_a, g_b, g_c for three
% more, B_k - B_r = M*(I_k - I_r) gives M = [g_a*Bh_a - Bh_r, ...] / [dI_a dI_b dI_c] exactly,
% and b0 = Bh_r - M*I_r. Every other spectrum k is then predicted, and agrees if
% canon(predicted) is within tau Gauss of canon(Bh_k) -- canon is the same for all 48
% labels of a vector, so the others need no label. 48^3 = 110592 combinations per four,
% vectorized. The four are drawn from the unambiguous spectra, with well-spread currents.
    starts = struct('R', {}, 'cost', {}, 'p0', {}, 'label', {}, 'nInliers', {});
    if numel(goodIdx) < 4; return; end
    G48 = symmetryGroup();
    N = size(I, 1);
    best = struct('M', {}, 'b0', {}, 'n', {}, 'res', {});
    for s = 1:nSubsets
        pick = [];
        for tries = 1:200
            cand = goodIdx(randperm(rs, numel(goodIdx), 4));
            dI = (I(cand(2:4), :) - I(cand(1), :))';
            if rcond(dI) > 0.02; pick = cand; break; end
        end
        if isempty(pick); continue; end
        r = pick(1); abc = pick(2:4);
        Br = Bhat(r, :)';
        dI = (I(abc, :) - I(r, :))';
        dIinv = dI \ eye(3);
        % all 48 images of each of a, b, c, minus Br: 3 x 48
        Dm = cell(1, 3);
        for q = 1:3
            Dm{q} = reshape(pagemtimes(G48, Bhat(abc(q), :)'), 3, 48) - Br;
        end
        % prediction for spectrum k at every (ia, ib, ic):
        %   Br + w1*Da(:,ia) + w2*Db(:,ib) + w3*Dc(:,ic),  w = dIinv*(I_k - I_r)
        nIn = zeros(48, 48, 48); res = zeros(48, 48, 48);
        for k = solvedIdx(:)'
            w = dIinv * (I(k, :) - I(r, :))';
            X = Br(1) + w(1)*reshape(Dm{1}(1,:),48,1,1) + w(2)*reshape(Dm{2}(1,:),1,48,1) + w(3)*reshape(Dm{3}(1,:),1,1,48);
            Y = Br(2) + w(1)*reshape(Dm{1}(2,:),48,1,1) + w(2)*reshape(Dm{2}(2,:),1,48,1) + w(3)*reshape(Dm{3}(2,:),1,1,48);
            Z = Br(3) + w(1)*reshape(Dm{1}(3,:),48,1,1) + w(2)*reshape(Dm{2}(3,:),1,48,1) + w(3)*reshape(Dm{3}(3,:),1,1,48);
            X = abs(X); Y = abs(Y); Z = abs(Z);
            hi = max(max(X, Y), Z); lo = min(min(X, Y), Z); mid = X + Y + Z - hi - lo;
            d2 = (hi - Bhat(k, 1)).^2 + (mid - Bhat(k, 2)).^2 + (lo - Bhat(k, 3)).^2;
            nIn = nIn + (d2 < tau^2);
            res = res + min(d2, tau^2);
        end
        sc = nIn * 1e3 - res / tau^2;           % most inliers, then the tightest
        [~, ord] = sort(sc(:), 'descend');
        for t = ord(1:3)'
            [ia, ib, ic] = ind2sub([48 48 48], t);
            M = [Dm{1}(:, ia), Dm{2}(:, ib), Dm{3}(:, ic)] * dIinv;
            b0 = Br - M * I(r, :)';
            best(end+1) = struct('M', M, 'b0', b0, 'n', nIn(t), 'res', res(t));  %#ok<AGROW>
        end
    end
    if isempty(best); return; end
    % refine each candidate by linear least squares on its inliers (each with its own
    % best label), twice, then keep the best few DISTINCT ones
    for c = 1:numel(best)
        M = best(c).M; b0 = best(c).b0;
        for it = 1:3
            [inl, Bl] = labelAndInliers(M, b0, I, Bhat, solvedIdx, tau, G48);
            if nnz(inl) < 5; break; end
            X = [I(inl, :), ones(nnz(inl), 1)];
            sol = X \ Bl(inl, :);
            M = sol(1:3, :)'; b0 = sol(4, :)';
        end
        [inl, ~] = labelAndInliers(M, b0, I, Bhat, solvedIdx, tau, G48);
        best(c).M = M; best(c).b0 = b0; best(c).n = nnz(inl);
    end
    [~, ord] = sort([best.n], 'descend');
    best = best(ord);
    kept = best(1);
    for c = 2:numel(best)
        if numel(kept) >= nKeep; break; end
        isNew = true;
        for m = 1:numel(kept)
            g = pickGauge(best(c).M, kept(m).M);
            if norm(g * best(c).M - kept(m).M) < 0.02 * norm(kept(m).M); isNew = false; break; end
        end
        if isNew; kept(end+1) = best(c); end                                   %#ok<AGROW>
    end
    for c = 1:numel(kept)
        starts(end+1) = struct('R', eye(3), 'cost', NaN, ...
            'p0', [kept(c).M(:); kept(c).b0(:)], 'label', 'RANSAC on per-spectrum fields', ...
            'nInliers', kept(c).n);                                             %#ok<AGROW>
    end
end

function [inl, Bl] = labelAndInliers(M, b0, I, Bhat, idx, tau, G48)
% For each spectrum, the label that brings its own field closest to the model's.
    N = size(I, 1);
    inl = false(N, 1); Bl = NaN(N, 3);
    for k = idx(:)'
        pr = M * I(k, :)' + b0;
        imgs = reshape(pagemtimes(G48, Bhat(k, :)'), 3, 48);
        [d2, j] = min(sum((imgs - pr).^2, 1));
        Bl(k, :) = imgs(:, j)';
        inl(k) = d2 < tau^2;
    end
end

function G = symmetryGroup()
% The 48 signed permutation matrices, 3 x 3 x 48.
    P = perms(1:3);
    G = zeros(3, 3, 48); n = 0;
    for i = 1:6
        Pm = eye(3); Pm = Pm(P(i, :), :);
        for s = 0:7
            n = n + 1;
            G(:, :, n) = diag(1 - 2*[bitget(s, 1), bitget(s, 2), bitget(s, 3)]) * Pm;
        end
    end
end

% ======================================================================
%   rotations and the gauge
% ======================================================================

function R = quatToRot(q)
    q = q / norm(q);
    w = q(1); x = q(2); y = q(3); z = q(4);
    R = [1-2*(y^2+z^2), 2*(x*y-z*w),   2*(x*z+y*w);
         2*(x*y+z*w),   1-2*(x^2+z^2), 2*(y*z-x*w);
         2*(x*z-y*w),   2*(y*z+x*w),   1-2*(x^2+y^2)];
end

function R = rotvecToRot(v)
    th = norm(v);
    if th < 1e-12; R = eye(3); return; end
    k = v(:) / th;
    K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    R = eye(3) + sin(th) * K + (1 - cos(th)) * (K * K);
end

function same = sameModuloCube(R1, R2)
% True when R1 = g*R2 for one of the 24 proper cube rotations: then the two starts are
% the same orientation with the axes named differently.
    X = R1 * R2';
    same = all(abs(abs(X(:)) - round(abs(X(:)))) < 0.05) && ...
           all(sum(abs(round(X)), 1) == 1);
end

function g = pickGauge(M, ref)
% The signed permutation g (one of 48) that puts g*M closest to ref, by column DIRECTION,
% so a nominal matrix with the wrong scale still picks the right frame.
    Mn = M ./ sqrt(sum(M.^2, 1));
    Rn = ref ./ max(sqrt(sum(ref.^2, 1)), eps);
    P = perms(1:3);
    best = -inf; g = eye(3);
    for i = 1:size(P, 1)
        Pm = eye(3); Pm = Pm(P(i, :), :);
        for s = 0:7
            gi = diag(1 - 2*[bitget(s, 1), bitget(s, 2), bitget(s, 3)]) * Pm;
            sc = sum(sum((gi * Mn) .* Rn));
            if sc > best + 1e-12; best = sc; g = gi; end
        end
    end
end

% ======================================================================

function report(out, opts)
    fprintf('\n============================================================\n');
    fprintf('coilMatrixFit: B = M*I + b0, in the NV crystal frame\n');
    fprintf('============================================================\n');
    fprintf('  %d spectra (%d unambiguous on their own), %d lines matched, %d outliers\n', ...
            out.nSpectra, out.nUsed, out.nLines, out.nOutliers);
    fprintf('  metric: %s\n  winning start: %s (orientation search rms %.3f G)\n', ...
            out.init.source, out.init.startUsed, out.init.rotRmsG);
    if ~isempty(out.init.downweighted)
        fprintf('  spectra the |B|^2 fit down-weighted as inconsistent: %s\n', ...
                mat2str(out.init.downweighted));
    end
    if isfinite(out.init.c) && startsWith(out.init.source, 'linear')
        fprintf(['  |b0|^2: free fit %.3f G^2, implied by the rest %.3f G^2 (agree within ' ...
                 'noise = the linear model holds)\n'], out.init.c, out.init.cCheck);
    end
    fprintf('\n  line residual rms %.3f MHz, reduced chi2 %.2f\n', out.rmsMHz, out.chi2Red);

    ax = 'xyz';
    fprintf('\n  M (G/A)          coil x               coil y               coil z\n');
    for i = 1:3
        fprintf('    B%s (NV)', ax(i));
        for k = 1:3
            fprintf('  %8.4f +/- %6.4f', out.M(i, k), out.Merr(i, k));
        end
        fprintf('\n');
    end
    fprintf('\n  b0, the field at zero current = [%s] +/- [%s] G\n', ...
            sprintf(' %+.3f', out.b0), sprintf(' %.3f', out.b0err));
    if numel(out.D) == 1
        fprintf('  D = %.3f +/- %.3f MHz\n', out.D, out.Derr);
    else
        fprintf('  D per spectrum: %.3f to %.3f MHz\n', min(out.D), max(out.D));
    end

    nom = [];
    if ~isempty(opts.nominalM); nom = sqrt(sum(opts.nominalM.^2, 1)); end
    fprintf('\n  coil   G/A (nominal)     axis in the NV frame: theta (deg)  phi (deg)\n');
    for k = 1:3
        if isempty(nom); ns = '      '; else; ns = sprintf('%6.3f', nom(k)); end
        fprintf('    %s   %7.3f (%s)                       %7.2f    %7.2f\n', ax(k), ...
                out.coilSensitivity(k), ns, out.coilThetaPhi(k, 1), out.coilThetaPhi(k, 2));
    end
    A = out.coilAnglesDeg;
    fprintf('  angles between the coil fields: xy %.2f, xz %.2f, yz %.2f deg (90 = orthogonal)\n', ...
            A(1,2), A(1,3), A(2,3));

    fprintf('\n  %3s  %21s  %8s  %8s  %6s  %s\n', '#', 'I (A)', '|B| fit', '|B| own', ...
            'rms', 'outliers (MHz)');
    for j = 1:numel(out.spectra)
        s = out.spectra(j);
        oStr = '';
        if any(s.outlier); oStr = sprintf('%.1f ', s.centres(s.outlier)); end
        own = sprintf('%8.3f', s.BmagSolve);
        if s.ambiguous; own = [own(1:end-1) '?']; end
        fprintf('  %3d  [%5.3f %5.3f %5.3f]  %8.3f  %s  %6.3f  %s\n', j, s.I, s.Bmag, ...
                own, s.rms, oStr);
    end
    fprintf('  (? = that spectrum alone allows more than one field; the global fit decides)\n');
    if ~isempty(out.why); fprintf('\n  NOTE %s\n', out.why); end
    fprintf('============================================================\n\n');
end

function m = nanMedian(x)
    x = x(isfinite(x));
    if isempty(x); m = NaN; else; m = median(x); end
end

function o = withDefaults(o, d)
    if isempty(o); o = struct(); end
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(o, fn{i}) || (isempty(o.(fn{i})) && ~isempty(d.(fn{i})))
            o.(fn{i}) = d.(fn{i});
        end
    end
end
