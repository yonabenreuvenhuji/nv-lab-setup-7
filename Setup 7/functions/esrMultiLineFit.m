function out = esrMultiLineFit(freqs, signal, sterr, opts)
%ESRMULTILINEFIT Fit every line in an ESR sweep at once, with a shared baseline.
%
%   out = esrMultiLineFit(f, y, e)
%   out = esrMultiLineFit(f, y, e, struct('nLines', 4, 'baselineOrder', 2))
%
% Replaces the hand-transcribed myESRFit buried in yyMagneticFieldAlignment. A field
% splits the four NV orientations into up to eight lines, and the quantities that
% matter -- splittings, their ratios, their means -- are DIFFERENCES AND SUMS of line
% centres. So all lines have to be fitted together, and the covariance between their
% centres has to come out with them.
%
% ===================================================================
% WHY NOT JUST CALL esrLineFit REPEATEDLY
% ===================================================================
%
% Three reasons, in increasing order of importance.
%
% The baseline is shared. Fitting each line in its own window gives each one its own
% independent baseline, and the resonator envelope is then absorbed differently by each
% -- which puts a frequency-dependent offset into exactly the centres being compared.
%
% Lines overlap. At small field, or when misalignment splits the three-fold line, the
% components sit within a linewidth of each other. A single-line fit to overlapping
% lines reports a centre somewhere between them, biased and with a plausible error bar.
%
% AND THE ERRORS ON THE THINGS YOU ACTUALLY WANT NEED THE COVARIANCE. The splitting is
% c2 - c1 and its variance is C11 + C22 - 2*C12. The cross term is not optional: the
% centres are correlated through the shared baseline, and ignoring C12 can be wrong in
% either direction. confint, which esrLineFit uses, does not expose covariance at all.
% That is why this fits by hand rather than through the Curve Fitting Toolbox: an
% explicit (J'WJ)^-1 gives the full matrix, and out.centreCov carries it.
%
% ===================================================================
% THE MODEL
% ===================================================================
%
%   y(x) = P(x) - SUM_k  A_k * g_k^2 / ((x - c_k)^2 + g_k^2)
%
% P is a polynomial baseline of order opts.baselineOrder in (x - xref).
%
% ORDER 1 IS THE DEFAULT AND ORDER 2 EXISTS FOR A SPECIFIC REASON. The plan's section
% 3.4 identified a real systematic: `amplitude` is a scalar for a whole sweep, so B1
% varies across it by whatever the resonator does, and an asymmetric envelope under a
% symmetric line pulls the fitted centre -- a direct error in the field alignment. Its
% proposed fix was per-frequency amplitude flattening, and that fix is NOT AVAILABLE:
% calibResonator established there is no reflection standard on this setup, so the
% flattening table produced |Gamma| > 1 and +156 dB entries and was removed.
%
% So the envelope has to be absorbed by the baseline instead. Over a narrow window it
% is essentially linear and order 1 is enough. Over the 200+ MHz a field survey needs,
% it curves, and order 1 would leave that curvature to be soaked up by the line
% parameters. Use order 2 for wide sweeps, order 1 for narrow refits.
%
% ===================================================================
% ERRORS
% ===================================================================
%
% Two-sided Birge scaling, the same convention and for the same measured reason as
% esrLineFit: AWGExpESR's sterr is dominated by drift BETWEEN averages, common to every
% frequency in the sweep and measured at 42x the real point-to-point scatter. So the
% input errors carry the right SHAPE and the wrong SCALE, and parameter errors are
% scaled by sqrt(reduced chi2) in BOTH directions, floored at 0.05. Refusing to scale
% down is not conservative, it is useless -- it once rejected a whole channel of
% calibMwDrive on dips that were plainly there.
%
% ===================================================================
% OPTIONS (struct fields, all optional)
%   nLines         how many lines to fit ([] = find them, default [])
%   seedCentres    MHz, force the fit to start from these centres (default [] = find
%                  them). This is how a line the finder MERGED with its neighbour gets
%                  recovered: the caller knows the spectrum must be symmetric about D,
%                  so it can predict a missing partner at 2D - f and hand it in here.
%                  Physics supplying a start point, not a free parameter.
%   maxLines       cap when finding automatically (default 8 -- four orientations,
%                  two transitions each)
%   baselineOrder  1 or 2 (default 1; use 2 for sweeps wider than ~100 MHz)
%   snr            a candidate must be this many noise units deep (default 4)
%   minSepMHz      minimum separation between candidate centres (default 6)
%   minDepthSigma  drop a fitted line shallower than this (default 3)
%   maxIter        Gauss-Newton iterations (default 200)
%   plot           draw data, fit, and residuals (default false)
%   name           label for the plot
%
% OUTPUT (struct)
%   lines          struct array, sorted by centre: centre, centreErr, fwhm, fwhmErr,
%                  depth, depthErr
%   centreCov      nLines x nLines covariance of the centres, MHz^2. The splitting
%                  and pairMean tables below are derived from it.
%   baseline       polynomial coefficients, and xref they are about
%   chi2Reduced, birge, ok, why, residual
%
%   splitting, splittingErr    n x n, |c_j - c_i| and its error WITH the covariance
%                              term: var = Cii + Cjj - 2*Cij
%   pairMean, pairMeanErr      n x n, (c_i + c_j)/2 and its error, var uses +2*Cij.
%                              For a pair from ONE orientation this must equal D,
%                              so it is a free check on both the field and the fit.

    if ~exist('sterr', 'var'); sterr = []; end
    if ~exist('opts', 'var'); opts = struct(); end
    % maxFwhmMHz: the widest FWHM that can still be an NV line rather than baseline. Lines
    % on this setup run 7 MHz unsaturated to about 25 MHz hard-driven, so 40 leaves headroom
    % and still rejects a Lorentzian standing in for the baseline. See the tooWide test.
    opts = withDefaults(opts, struct('nLines', [], 'maxLines', 8, ...
        'baselineOrder', 1, 'snr', 4, 'minSepMHz', 6, 'minDepthSigma', 3, ...
        'maxFwhmMHz', 40, ...
        'maxIter', 200, 'plot', false, 'name', '', 'seedCentres', []));

    x = freqs(:); y = signal(:);
    if isempty(sterr); e = nan(size(x)); else; e = sterr(:); end
    e(~isfinite(e) | e <= 0) = NaN;
    good = isfinite(x) & isfinite(y);
    x = x(good); y = y(good); e = e(good);
    [x, ix] = sort(x); y = y(ix); e = e(ix);

    out = struct('ok', false, 'why', '', 'lines', [], 'centreCov', [], ...
                 'splitting', [], 'splittingErr', [], 'pairMean', [], 'pairMeanErr', [], ...
                 'baseline', [], 'xref', NaN, 'chi2Reduced', NaN, 'birge', NaN, ...
                 'residual', [], 'model', [], 'modelAt', [], 'nLines', 0);

    nB = opts.baselineOrder + 1;
    if numel(x) < nB + 6
        out.why = sprintf('only %d usable points', numel(x));
        return
    end
    out.xref = mean(x);

    % ---- find the lines ----------------------------------------------------
    [c0, a0, g0, noise, env] = findLines(x, y, opts);
    if ~isempty(opts.seedCentres)
        % Caller-supplied centres replace the finder's. Amplitude and width start points
        % still come from the data AT those positions, so a seeded line that is not
        % really there simply fails the depth test later and gets dropped.
        %
        % THE DEPTH START POINT COMES FROM THE GLOBAL ENVELOPE, NOT A LOCAL WINDOW.
        %
        % It used to be interp1(x, movmaxLocal(y, 7), cs) - y(cs): the running maximum
        % over 7 points, which on a 1.5 MHz grid spans 10.5 MHz. Every line here is 10 to
        % 44 MHz wide, so that window sits ENTIRELY INSIDE the line and the "local
        % maximum" is itself a point on the dip. base - yv then came out near zero and
        % every seeded line started at the 4*noise floor -- a fraction of its true depth.
        %
        % That is why seeded refits could not beat the unseeded fit. The first fit gets
        % its amplitudes from the finder, which measures the deficit under a robust
        % GLOBAL envelope and is therefore right; every refit after it started from
        % near-zero depths, and Levenberg-Marquardt with a shrinking trust region does
        % not always climb back out. The shallow component of a blend is exactly the one
        % that then fails its 2.5-sigma depth test and gets dropped, which returned the
        % refit to the same line count and made the caller reject it. Same envelope as
        % the finder, so a refit now starts where the first fit started.
        cs = sort(opts.seedCentres(:));
        cs = cs(cs >= min(x) & cs <= max(x));
        if ~isempty(cs)
            dEnv = env - y;                     % deficit: positive inside a dip
            c0 = cs;
            a0 = max(interp1(x, dEnv, cs, 'linear', 'extrap'), 4*noise);
            % Width start, capped by the distance to the nearest neighbouring seed. Two
            % seeds placed inside one 44 MHz blend must not both start as wide as the
            % blend, or they sit on top of each other and one is redundant from the first
            % iteration -- which is the other way a seeded component gets dropped.
            step = median(diff(x));
            gap  = inf(size(cs));
            if numel(cs) > 1
                d2 = diff(cs);
                gap(1:end-1) = min(gap(1:end-1), d2);
                gap(2:end)   = min(gap(2:end),   d2);
            end
            g0 = min(max(3*step, 2), max(gap/2, 2*step));
        end
    end
    if ~isempty(opts.nLines)
        % Asked for a specific count: keep the deepest that many, or pad by splitting
        % the sweep evenly if the finder came up short.
        if numel(c0) > opts.nLines
            [~, ord] = sort(a0, 'descend');
            keep = sort(ord(1:opts.nLines));
            c0 = c0(keep); a0 = a0(keep); g0 = g0(keep);
        elseif numel(c0) < opts.nLines
            need = opts.nLines - numel(c0);
            extra = linspace(min(x), max(x), need + 2);
            c0 = [c0(:); extra(2:end-1)'];
            a0 = [a0(:); repmat(max([a0(:); 4*noise]), need, 1)];
            g0 = [g0(:); repmat(max([median(g0), 5]), need, 1)];
        end
    end
    if isempty(c0)
        out.why = sprintf(['no line found deeper than %g noise units (noise %.2g). ' ...
                           'Either there is no resonance here or the drive is off.'], ...
                          opts.snr, noise);
        return
    end

    % ---- fit, dropping lines that do not survive ---------------------------
    % A line the finder proposed may fit to nothing, and one bad line distorts every
    % other parameter through the shared baseline. So refit without it rather than
    % reporting a fit that contains it.
    for attempt = 1:max(numel(c0), 1)
        [p, C, chi2, birge, resid, yhat, fitOk, why] = ...
            gnFit(x, y, e, c0, a0, g0, nB, out.xref, opts);
        if ~fitOk
            out.why = why; return
        end
        [A, c, g] = unpack(p, nB);
        eA = sqrt(max(diag(C(nB+1:3:end, nB+1:3:end)), 0));
        % Reject: too shallow to believe, centre pushed outside the sweep, or a width
        % narrower than the sampling can resolve.
        % Reject too-WIDE as well as too-narrow. A component far broader than a real line is
        % not a line, it is the baseline being fitted twice: the polynomial and a 100 MHz
        % Lorentzian are nearly degenerate, so the fit is free to put a huge shallow "line"
        % anywhere it likes and chi2 does not object. On a synthetic 45 G spectrum at 2
        % degrees of tilt this produced a 112 MHz component that was accepted, and it
        % distorts every real width in the fit -- which matters because the blend width IS
        % the alignment metric near alignment.
        %
        % THE THRESHOLD IS NOT A FRACTION OF THE SWEEP, AND THAT WAS A REAL BUG. It used to
        % be 0.25*span: fine on a 130 MHz survey, where it allows 32 MHz, and catastrophic
        % on a narrow window, where a 30 MHz sweep makes the limit 7.5 MHz and this setup's
        % lines are 10 to 13 MHz FWHM. A visibly clean single dip -- centre 2826.9 MHz,
        % depth 0.8 %, FWHM 10.2 -- was thrown out with "the only line left does not survive
        % its checks", and the caller reported "ESR fitted no lines". The width of a real
        % line has nothing to do with how wide a window you chose to look through, so the
        % threshold must not depend on it either.
        %
        % Two independent statements, and a component is rejected if either holds:
        %
        %   FWHM > the sweep span    you cannot measure a width you have not spanned
        %   FWHM > maxFwhmMHz        no NV line here is that broad, saturated or not
        %
        % which together is FWHM > min(span, maxFwhmMHz). The 112 MHz synthetic component is
        % still caught, by the second test -- its sweep was 350 MHz wide, so the first would
        % not have caught it -- and the 10.2 MHz real line now passes on a 30 MHz window.
        tooWide = 2*g > min(max(x) - min(x), opts.maxFwhmMHz);
        bad = (A < opts.minDepthSigma*eA) | (c < min(x)) | (c > max(x)) | ...
              (2*g < 1.5*median(diff(x))) | tooWide | ~isfinite(A) | ~isfinite(c);
        if ~any(bad)
            break
        end
        if numel(c0) == 1
            % The last remaining line failed its own test. Reporting it anyway is how a
            % diverged fit reached the caller and got treated as a measurement.
            out.why = sprintf(['the only line left does not survive its checks ' ...
                '(centre %.1f MHz, depth %.3g, FWHM %.2f). The fit did not converge on ' ...
                'anything usable.'], c(1), A(1), 2*g(1));
            return
        end
        drop = find(bad, 1);
        c0(drop) = []; a0(drop) = []; g0(drop) = [];
        if isempty(c0)
            out.why = 'every candidate line failed its depth or width test';
            return
        end
    end

    % ---- assemble ----------------------------------------------------------
    [A, c, g] = unpack(p, nB);
    n = numel(c);
    idxA = nB + (1:3:3*n);
    idxC = nB + (2:3:3*n);
    idxG = nB + (3:3:3*n);
    se = sqrt(max(diag(C), 0));

    [c, ord] = sort(c);
    A = A(ord); g = g(ord);
    idxA = idxA(ord); idxC = idxC(ord); idxG = idxG(ord);

    L = struct('centre', {}, 'centreErr', {}, 'fwhm', {}, 'fwhmErr', {}, ...
               'depth', {}, 'depthErr', {});
    for k = 1:n
        L(k) = struct('centre', c(k), 'centreErr', se(idxC(k)), ...
                      'fwhm', 2*g(k), 'fwhmErr', 2*se(idxG(k)), ...
                      'depth', A(k), 'depthErr', se(idxA(k)));
    end
    out.lines = L;
    out.nLines = n;
    out.centreCov = C(idxC, idxC);       % the whole point: sums and differences
    [out.splitting, out.splittingErr, out.pairMean, out.pairMeanErr] = ...
        pairTables(c, out.centreCov);
    out.baseline = p(1:nB)';
    out.chi2Reduced = chi2; out.birge = birge;
    out.residual = resid; out.model = yhat;

    % out.model is the fitted vector AT THE DATA POINTS. out.modelAt is the same model as a
    % function handle, so a caller can draw a smooth curve or evaluate between points
    % without re-deriving the formula. Added rather than replacing out.model, because that
    % field is what the residual plots already use.
    %
    % The two conventions worth stating, since getting either wrong gives a plausible curve
    % in the wrong place: the baseline coefficients are ASCENDING powers of (x - xref), not
    % MATLAB's polyval order, and the lines are SUBTRACTED because they are dips.
    out.modelAt = @(xx) esrModelAt(xx, out.baseline, out.xref, out.lines);

    % Final sanity, cheap and absolute. A depth above 1 in normalised fluorescence, or a
    % centre outside the swept range, means the fit diverged whatever chi2 says.
    if any(A > 1) || any(c < min(x) - 1e-6) || any(c > max(x) + 1e-6) || ...
       any(~isfinite([A(:); c(:); g(:)]))
        out.ok = false;
        out.why = sprintf(['the fit diverged: depths up to %.3g, centres %.1f to %.1f MHz ' ...
            'against a sweep of %.1f to %.1f.'], max(A), min(c), max(c), min(x), max(x));
        return
    end
    out.ok = true;

    if opts.plot; drawFit(x, y, e, yhat, resid, out, opts); end
end

% ======================================================================

function v = esrModelAt(xx, baseline, xref, L)
% y(x) = P(x) - SUM_k A_k * g_k^2 / ((x - c_k)^2 + g_k^2),  P ascending in (x - xref).
    xx = xx(:);
    v = zeros(size(xx));
    b = baseline(:);
    for j = 1:numel(b)
        v = v + b(j) * (xx - xref).^(j-1);
    end
    for k = 1:numel(L)
        g = L(k).fwhm/2;
        v = v - L(k).depth * g^2 ./ ((xx - L(k).centre).^2 + g^2);
    end
end

function [c0, a0, g0, noise, env] = findLines(x, y, opts)
% Candidate centres from the DEFICIT below a robust upper envelope.
%
% The envelope is fitted iteratively to the points that are NOT in a dip: start from
% all points, then keep only those at or above the current line, and repeat. Using all
% points once would drag the baseline down into the dips and the deficit would
% understate every line -- which is how a shallow line gets missed entirely.
    x = x(:); y = y(:);
    keep = true(size(x));
    b = [0; median(y)];
    for it = 1:6
        X = [x(keep) - mean(x), ones(sum(keep), 1)];
        b = X \ y(keep);
        base = b(1)*(x - mean(x)) + b(2);
        r = y - base;
        s = 1.4826*median(abs(r - median(r)));
        if ~isfinite(s) || s == 0; break; end
        newKeep = r > -0.5*s;
        if isequal(newKeep, keep) || sum(newKeep) < 4; break; end
        keep = newKeep;
    end
    base = b(1)*(x - mean(x)) + b(2);
    d = base - y;                       % deficit: positive inside a dip
    env = base;                         % returned so SEEDED refits use the same envelope

    % NOISE FROM POINT-TO-POINT SCATTER, not from the spread of the deficit.
    %
    % The first version took the MAD of the deficit itself after trimming the top. That
    % fails exactly when it matters: with several lines present the deficit's spread is
    % dominated by the LINES, so the "noise" comes out at line depth, the threshold ends
    % up above the shallow lines, and they are silently dropped. On a synthetic <111>
    % pattern it found the two deep three-fold lines and missed both aligned lines --
    % which are the ones the field magnitude is measured from.
    %
    % Successive differences cancel any smooth structure, so median|diff(y)| sees only
    % noise wherever the spectrum is locally flat, and the median is robust to the
    % minority of points sitting on a line flank. For Gaussian noise
    % median|diff| = 0.6745 * sigma * sqrt(2), hence the 1.4826/sqrt(2).
    noise = 1.4826 * median(abs(diff(y))) / sqrt(2);
    if ~isfinite(noise) || noise <= 0
        noise = 1.4826*median(abs(d - median(d)));
    end
    if ~isfinite(noise) || noise <= 0; noise = eps; end

    ds = med3local(d);                   % kill single-point spikes before peak finding
    thresh = opts.snr * noise;

    % Local maxima above threshold, then greedy selection by depth with a minimum
    % separation, so one broad line does not register as three.
    isPeak = [false; ds(2:end-1) >= ds(1:end-2) & ds(2:end-1) >= ds(3:end); false];
    cand = find(isPeak & ds > thresh);
    if isempty(cand)
        c0 = []; a0 = []; g0 = []; return
    end
    [~, ord] = sort(ds(cand), 'descend');
    cand = cand(ord);
    chosen = [];
    for k = cand(:)'
        if all(abs(x(k) - x(chosen)) >= opts.minSepMHz)
            chosen(end+1) = k;                                   %#ok<AGROW>
        end
        if numel(chosen) >= opts.maxLines; break; end
    end
    chosen = sort(chosen);

    c0 = x(chosen);
    a0 = ds(chosen);
    % Width start point from where the deficit falls to half, measured outward from
    % each peak; falls back to a few frequency steps if the flank runs into a neighbour.
    g0 = zeros(size(c0));
    for k = 1:numel(chosen)
        i0 = chosen(k);
        half = ds(i0)/2;
        iR = i0; while iR < numel(ds) && ds(iR) > half; iR = iR + 1; end
        iL = i0; while iL > 1        && ds(iL) > half; iL = iL - 1; end
        w = x(iR) - x(iL);
        if ~isfinite(w) || w <= 0; w = 4*median(diff(x)); end
        g0(k) = max(w/2, 1.5*median(diff(x)));
    end
end

function [p, C, chi2, birge, resid, yhat, ok, why] = gnFit(x, y, e, c0, a0, g0, nB, xr, opts)
% Weighted Levenberg-Marquardt on the explicit model, with an analytic Jacobian.
%
% Hand-rolled rather than through fit(): the covariance between line centres is the
% output this function exists for, and confint does not provide it. Analytic
% derivatives because a Lorentzian's are trivial and finite differences on a
% (2 + 3N)-parameter fit near a narrow line are not reliable.
    ok = false; why = ''; C = []; chi2 = NaN; birge = 1; resid = []; yhat = [];
    n = numel(c0);
    p = zeros(nB + 3*n, 1);
    % Baseline start: robust straight line, higher orders start at zero.
    X0 = [ones(numel(x),1), (x - xr)];
    b0 = X0 \ y;
    p(1) = b0(1);
    if nB >= 2; p(2) = b0(2); end
    for k = 1:n
        p(nB + 3*k - 2) = max(a0(k), eps);
        p(nB + 3*k - 1) = c0(k);
        p(nB + 3*k    ) = max(g0(k), eps);
    end

    w = 1 ./ e.^2;
    if all(~isfinite(w)); w = ones(size(x)); end
    w(~isfinite(w)) = median(w(isfinite(w)));
    w = w(:);

    % Bounds for the clamp above, all taken from the data.
    xLo = min(x); xHi = max(x);
    dxm = median(diff(x));
    gMin = 0.5*dxm;                 % a HWHM under half a sample is not a resolved line
    gMax = (xHi - xLo)/2;           % nor is one as wide as the whole sweep

    lambda = 1e-3;
    converged = false;
    cost = @(pp) sum(w .* (y - modelOf(pp, x, nB, xr)).^2);
    fPrev = cost(p);
    for it = 1:opts.maxIter
        [f, J] = modelAndJac(p, x, nB, xr);
        r = y - f;
        JW = J .* w;
        H = J' * JW;
        gr = JW' * r;
        % COLUMN SCALING, Marquardt's own. Parameters differ by ten orders of magnitude
        % here -- centres are ~2870 MHz, depths ~0.01 -- so the normal matrix is wildly
        % badly scaled and the raw solve produced rcond down to 1e-37 and a stream of
        % "matrix singular to machine precision" warnings. Scaling each column of J to
        % unit norm makes the solve well conditioned without changing the step it
        % represents: solve in scaled variables, then unscale.
        sc = sqrt(max(diag(H), eps));
        Ds = diag(1 ./ sc);
        Hs = Ds * H * Ds;
        gs = Ds * gr;
        stepTaken = false;
        for tries = 1:12
            try
                dps = (Hs + lambda*eye(size(Hs))) \ gs;
                dp  = Ds * dps;
            catch
                dp = Ds * (pinv(Hs + lambda*eye(size(Hs))) * gs);
            end
            pTry = p + dp;
            % BOUND EVERY LINE PARAMETER ON EVERY STEP.
            %
            % Amplitude and width were already forced positive. The CENTRE was not, and
            % that was a real divergence: a fit came back with a centre at -22,034,007 MHz
            % and a depth of 2492%. The mechanism is that once a centre is far enough
            % outside the data its Lorentzian is zero everywhere, so the model collapses to
            % baseline-only, the cost DECREASES, and the line-search happily accepts the
            % step. The line has escaped, and with a quadratic baseline able to absorb a
            % good deal of the spectrum on its own, escaping can genuinely look cheaper
            % than sitting in a slightly wrong place.
            %
            % Clamping is the fix rather than a post-hoc check: a diverged parameter
            % corrupts the shared baseline and every other line before anyone inspects it.
            pTry(nB+1:3:end) = abs(pTry(nB+1:3:end));                    % depth >= 0
            gg = abs(pTry(nB+3:3:end));
            gg = min(max(gg, gMin), gMax);
            pTry(nB+3:3:end) = gg;                                       % width in range
            cc = pTry(nB+2:3:end);
            pTry(nB+2:3:end) = min(max(cc, xLo), xHi);                   % centre in range
            fTry = cost(pTry);
            if isfinite(fTry) && fTry <= fPrev
                rel = abs(fPrev - fTry) / max(fPrev, eps);
                p = pTry; fPrev = fTry; lambda = max(lambda/3, 1e-9);
                stepTaken = true;
                converged = rel < 1e-12;
                break
            end
            lambda = lambda * 4;
        end
        if ~stepTaken || converged; break; end
    end

    [yhat, J] = modelAndJac(p, x, nB, xr);
    resid = y - yhat;
    dof = max(numel(x) - numel(p), 1);
    ee = e; ee(~isfinite(ee)) = median(e(isfinite(e)));
    if all(~isfinite(ee)); ee = ones(size(x))*std(resid); end
    chi2 = sum((resid./ee).^2) / dof;

    % TWO-SIDED Birge, floored -- see the header. The sterr scale is not trusted; the
    % residuals are.
    if isfinite(chi2) && chi2 > 0
        birge = max(sqrt(chi2), 0.05);
    end

    JW = J .* w;
    H = J' * JW;
    % Same scaling for the covariance, so a badly scaled problem does not produce a
    % meaningless matrix (or a warning) here either.
    sc = sqrt(max(diag(H), eps));
    Ds = diag(1 ./ sc);
    Hs = Ds * H * Ds;
    if rcond(Hs) < 1e-12
        why = ['the fit is degenerate (singular normal matrix even after scaling) -- ' ...
               'too many lines, or a line with no data under it'];
        C = Ds * pinv(Hs) * Ds * birge^2;
    else
        C = Ds * inv(Hs) * Ds * birge^2;
    end
    % inv(J'WJ) is the covariance only when the weights are 1/sigma^2. They are, up to
    % the scale factor Birge corrects for, which is why it multiplies here.
    if ~all(isfinite(C(:)))
        why = 'covariance came out non-finite'; return
    end
    ok = true;
end

function f = modelOf(p, x, nB, xr)
    f = polyBase(p, x, nB, xr);
    n = (numel(p) - nB)/3;
    for k = 1:n
        A = p(nB+3*k-2); c = p(nB+3*k-1); g = p(nB+3*k);
        f = f - A * g^2 ./ ((x - c).^2 + g^2);
    end
end

function [f, J] = modelAndJac(p, x, nB, xr)
    n = (numel(p) - nB)/3;
    f = polyBase(p, x, nB, xr);
    J = zeros(numel(x), numel(p));
    for j = 1:nB
        J(:, j) = (x - xr).^(j-1);
    end
    for k = 1:n
        A = p(nB+3*k-2); c = p(nB+3*k-1); g = p(nB+3*k);
        u  = x - c;
        dn = u.^2 + g^2;
        L  = g^2 ./ dn;
        f  = f - A * L;
        J(:, nB+3*k-2) = -L;                          % d/dA
        J(:, nB+3*k-1) = -A * 2*g^2 .* u ./ dn.^2;    % d/dc
        J(:, nB+3*k  ) = -A * 2*g   .* u.^2 ./ dn.^2; % d/dg
    end
end

function b = polyBase(p, x, nB, xr)
    b = zeros(size(x));
    for j = 1:nB
        b = b + p(j) * (x - xr).^(j-1);
    end
end

function [A, c, g] = unpack(p, nB)
    n = (numel(p) - nB)/3;
    A = p(nB + (1:3:3*n)); c = p(nB + (2:3:3*n)); g = p(nB + (3:3:3*n));
end

function [S, SE, M, ME] = pairTables(c, C)
% Every pairwise splitting and pair-mean, WITH the covariance term. These are the
% quantities the field work actually consumes, and getting their errors right is the
% reason this file fits all lines together instead of one at a time:
%
%   splitting  c_j - c_i        var = Cii + Cjj - 2*Cij
%   pair mean  (c_i + c_j)/2    var = (Cii + Cjj + 2*Cij)/4
%
% The cross term is not a refinement. The centres are correlated through the shared
% baseline, so dropping Cij can be wrong in either direction, and for a well-determined
% pair on a common baseline it is usually a large negative correlation -- which makes
% the splitting MORE precise than the individual centres, not less.
    n = numel(c);
    S = nan(n); SE = nan(n); M = nan(n); ME = nan(n);
    for i = 1:n
        for j = 1:n
            S(i,j)  = c(j) - c(i);
            SE(i,j) = sqrt(max(C(i,i) + C(j,j) - 2*C(i,j), 0));
            M(i,j)  = (c(i) + c(j))/2;
            ME(i,j) = sqrt(max(C(i,i) + C(j,j) + 2*C(i,j), 0)) / 2;
        end
    end
end

function y = movmaxLocal(x, k)
% Local maximum over a window, used only as a baseline estimate at a seeded centre.
    n = numel(x); y = x; h = floor(k/2);
    for i = 1:n
        y(i) = max(x(max(1,i-h):min(n,i+h)));
    end
end

function y = med3local(x)
    y = x;
    for i = 2:numel(x)-1
        y(i) = median(x(i-1:i+1));
    end
end

function drawFit(x, y, e, yhat, resid, out, opts)
    figure('Name', ['esrMultiLineFit ' opts.name], 'Color', 'w');
    subplot(3,1,[1 2]); hold on; box on;
    if any(isfinite(e)); errorbar(x, y, e, '.', 'Color', [.4 .6 .9]);
    else; plot(x, y, '.', 'Color', [.4 .6 .9]); end
    plot(x, yhat, 'r-', 'LineWidth', 1.2);
    for k = 1:out.nLines
        xline(out.lines(k).centre, 'k:');
        text(out.lines(k).centre, min(y), sprintf(' %.2f', out.lines(k).centre), ...
             'Rotation', 90, 'FontSize', 8, 'VerticalAlignment', 'bottom');
    end
    ylabel('normalized FL');
    title(sprintf('%s%d lines, reduced chi2 %.2f (errors scaled %.2fx)', ...
                  ternary(isempty(opts.name), '', [opts.name '  ']), ...
                  out.nLines, out.chi2Reduced, out.birge));
    subplot(3,1,3); hold on; box on;
    plot(x, resid, '.-'); yline(0, 'k:');
    xlabel('MHz'); ylabel('residual');
end

function v = ternary(c, a, b)
    if c; v = a; else; v = b; end
end

function o = withDefaults(o, d)
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(o, fn{i}) || isempty(o.(fn{i})); o.(fn{i}) = d.(fn{i}); end
    end
end
