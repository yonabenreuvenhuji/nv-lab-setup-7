function out = resonatorLorentzianFit(f, w, wErr, opts)
%RESONATORLORENTZIANFIT Fit the antenna response from Rabi frequencies measured across it.
%
%   out = resonatorLorentzianFit(fLine, fRabi, fRabiErr)
%   out = resonatorLorentzianFit(fLine, fRabi, fRabiErr, struct('baseline', false))
%
% fLine   the NV line frequency at each field setting, MHz
% fRabi   the Rabi frequency measured at that line, MHz
% fRabiErr its error, MHz. Pass [] and every point is weighted equally.
%
% ===================================================================
% WHY THE FIT IS IN OMEGA SQUARED, AND WHY AGAINST FREQUENCY
% ===================================================================
%
% The Rabi frequency is proportional to the driving field:
%
%       Omega = gamma * B1 * <matrix element>
%
% and the matrix element is constant across a field scan of a few gauss to about one
% part in 5000 -- the state mixing that would change it goes as (gamma*Bperp/D)^2, which
% at 20 G is (56/2870)^2 = 4e-4 even for a fully off-axis orientation. So Omega IS B1,
% up to one fixed scale factor, and nothing else in it moves.
%
% A resonator delivers POWER as a Lorentzian in frequency, and B1 goes as the square root
% of the power. So B1 itself is a square root of a Lorentzian, which is not a Lorentzian
% and must not be fitted as one. Omega SQUARED is the Lorentzian:
%
%       Omega^2(f) = A / (1 + 4 (f - f0)^2 / W^2)  +  c
%
% with f0 the centre, W the FWHM in power, A the on-resonance height and c a floor from
% whatever couples to the NVs without going through the resonator. Q = f0/W.
%
% Fitting against FREQUENCY rather than against the field setting matters just as much.
% The frequency is what the resonator responds to; the field is only the knob that moves
% the NV line through it. Fitting a parabola in Brho would use three points near the top,
% discard the wings, return no physical parameters, and hide any nonlinearity in the coil
% calibration inside the answer. Fitting in frequency uses every point, and the wings are
% where W and therefore Q come from.
%
% ===================================================================
% THERE IS NO START POINT TO GET WRONG
% ===================================================================
%
% At fixed (f0, W) the model is LINEAR in A and c, so those two are solved exactly by
% weighted least squares rather than searched. That leaves a two-dimensional surface,
% which is scanned exhaustively on a grid before any iteration begins. A Lorentzian fit
% started from a guess routinely converges onto a wing or onto the baseline, and returns
% a plausible centre that is simply wrong; scanning cannot, because it evaluates the
% whole surface. This is the same variable-projection structure as rabiFitRobust and it
% is used here for the same reason.
%
% SELF-CONTAINED: base MATLAB only. No Curve Fitting Toolbox, no lab code, no findpeaks.
% fminsearch is a base function and is used only to polish a grid answer, never to find
% one.
%
% ===================================================================
% OUTPUT
% ===================================================================
%
%   ok, why            did it fit, and if not, why not
%   f0, f0Err          resonator centre, MHz
%   fwhm, fwhmErr      power FWHM, MHz
%   Q, QErr            f0/fwhm
%   A, AErr, c, cErr   Lorentzian height and floor, in (MHz)^2 of Omega^2
%   wMax, wMaxErr      the Rabi frequency at the centre, MHz = sqrt(A + c)
%   piMinNs, piMinErrNs the shortest pi pulse that implies, ns = 500/wMax
%   contrastOfFit      (A)/(A+c): how much of the drive comes through the resonator
%   chi2Reduced, birge, r2
%   centreInRange      false when f0 fell outside the measured frequencies -- the single
%                      most important flag on this page, see below
%   resolved           false when the fitted FWHM is wider than the span measured, in
%                      which case f0 is an extrapolation and Q is meaningless
%   model              function handle f -> Omega^2, for plotting
%
% A fitted centre outside the measured range is NOT a result. It is the fit telling you
% the scan did not bracket the peak, and the honest response is to extend the field scan,
% not to quote the number. The same is true of a FWHM wider than the span: a Lorentzian
% shoulder is very nearly a parabola, so any short arc admits a family of (f0, W) that
% fit equally well and trade off against each other. Both cases are flagged rather than
% left for you to notice.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('baseline', true, 'nGridF', 601, 'nGridW', 81, ...
                                     'polish', true, 'verbose', true, 'retryNoBaseline', true), ...
                        'resonatorLorentzianFit');

    f = f(:); w = w(:);
    if isempty(wErr); wErr = ones(size(w)); else; wErr = wErr(:); end

    % Every field a caller might read is initialised HERE, so an early return still returns
    % a struct with the same shape as a successful one. Otherwise every call site needs an
    % isfield guard around every field, and the one that gets forgotten fails at the end of
    % a twenty-minute scan.
    out = struct('ok', false, 'usable', false, 'why', '', 'opts', opts, ...
                 'f0', NaN, 'f0Err', NaN, 'fwhm', NaN, 'fwhmErr', NaN, ...
                 'Q', NaN, 'QErr', NaN, 'wMax', NaN, 'wMaxErr', NaN, ...
                 'piMinNs', NaN, 'piMinErrNs', NaN, 'baselineDropped', false, ...
                 'centreInRange', false, 'resolved', false, 'pedestal', false, ...
                 'crossCheckSigma', NaN, ...
                 'parabola', struct('ok', false, 'inRange', false, 'f0', NaN, ...
                                    'f0Err', NaN, 'why', 'the fit did not get that far', ...
                                    'curvature', NaN));

    keep = isfinite(f) & isfinite(w) & isfinite(wErr) & w > 0 & wErr > 0;
    f = f(keep); w = w(keep); wErr = wErr(keep);
    nPar = 3 + double(opts.baseline);
    if numel(f) < nPar + 1
        out.why = sprintf('%d usable points, need at least %d for a %d-parameter fit', ...
                          numel(f), nPar + 1, nPar);
        return
    end

    % Omega^2 and its error. The square is exact; the error is first order, which is
    % right to better than a percent whenever the Rabi frequency is known to better than
    % 10% -- and a Rabi frequency known worse than that should have been gated out
    % upstream rather than propagated carefully here.
    y  = w.^2;
    sy = 2*w.*wErr;
    out.f = f; out.y = y; out.sy = sy; out.w = w; out.wErr = wErr;

    span = max(f) - min(f);
    if span <= 0
        out.why = 'all the line frequencies are the same -- the field scan did not move the line';
        return
    end

    % ---- exhaustive grid over the two nonlinear parameters ----------------
    % The centre is allowed well outside the measured range on purpose. Clamping it to
    % the data would guarantee a centre inside the scan whether or not one exists, which
    % is exactly the failure that has to be detectable.
    f0g = linspace(min(f) - 0.6*span, max(f) + 0.6*span, opts.nGridF);
    Wmin = 0.15*span; Wmax = 12*span;
    Wg  = exp(linspace(log(Wmin), log(Wmax), opts.nGridW));

    best = gridScan(f0g, Wg, f, y, sy, opts.baseline);
    if ~isfinite(best.cost)
        out.why = 'no positive-amplitude Lorentzian fits these points at all';
        return
    end

    % ---- polish. Search log(W) so the width cannot walk negative --------
    % AND KEEP IT INSIDE THE GRID'S BOUNDS. Left unbounded, fminsearch will happily run
    % the width to zero on a data set that only samples one wing: as W -> 0 with A -> inf
    % at fixed A*W^2 the model becomes a smooth power law through the points and the cost
    % keeps falling. One test scan over a pure wing came back with FWHM = 0.00 MHz, A =
    % 1.5e15 and a Rabi frequency of 38 million MHz, which is a perfectly good fit to
    % those six points and complete nonsense. The bounds turn that into a detectable
    % pinned parameter instead of a silent absurdity.
    f0 = best.f0; W = best.W;
    if opts.polish
        obj = @(p) costBounded(p, f, y, sy, opts.baseline, Wmin, Wmax);
        try
            o1 = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-9, 'MaxFunEvals', 4000);
            p = fminsearch(obj, [f0, log(W)], o1);
            if isfinite(obj(p)) && obj(p) <= best.cost
                f0 = p(1); W = min(max(exp(p(2)), Wmin), Wmax);
            end
        catch
            % Grid answer stands. A failed polish is not a failed fit -- the grid is the
            % measurement and the polish is a refinement of the last digit.
        end
    end
    [cost, lin] = costAt(f0, W, f, y, sy, opts.baseline);
    A = lin(1);
    if opts.baseline; c = lin(2); else; c = 0; end

    % ---- covariance over the full parameter vector -----------------------
    % Numerically, and over ALL FOUR parameters including the two that were solved
    % linearly. Reporting an error on f0 from the 2-D search alone would ignore that A
    % and c are correlated with it -- a taller Lorentzian with a floor can mimic a shifted
    % narrower one, and that correlation is most of the uncertainty on the centre.
    if opts.baseline; p = [f0, W, A, c]; else; p = [f0, W, A]; end
    dof = max(numel(f) - numel(p), 1);
    chi2r = cost/dof;
    birge = max(sqrt(chi2r), 0.05);     % floored: refusing to scale errors down by more
                                        % than 20x guards against a fluke-good chi2
    J = zeros(numel(f), numel(p));
    for k = 1:numel(p)
        step = 1e-6*max(abs(p(k)), 1);
        pp = p; pp(k) = pp(k) + step;
        pm = p; pm(k) = pm(k) - step;
        J(:,k) = (evalModel(f, pp, opts.baseline) - evalModel(f, pm, opts.baseline)) ./ (2*step*sy);
    end
    C = nan(numel(p));
    out.covOk = false;
    try
        sc = sqrt(sum(J.^2, 1));                % column scaling: f0 in MHz and A in MHz^2
        sc(sc == 0 | ~isfinite(sc)) = 1;        % differ by orders of magnitude and an
        Js = J ./ sc;                           % unscaled J'J inverts badly
        JtJ = Js'*Js;
        % CHECK THE CONDITIONING RATHER THAN INVERTING AND HOPING. A degenerate fit makes
        % J'J singular, MATLAB returns Inf with a warning, and the warning scrolls past
        % while the Inf becomes an error bar of 0.00 after max(diag,0). An unusable
        % covariance is reported as NaN.
        if rcond(JtJ) > 1e-12
            C = (inv(JtJ) ./ (sc'*sc)) * birge^2;
            out.covOk = all(isfinite(C(:)));
        end
    catch
        % Leave C as NaN. An error bar that could not be computed is reported as NaN, not
        % as zero and not as the unscaled version.
    end
    if ~out.covOk; C = nan(numel(p)); end
    err = sqrt(max(diag(C), 0))';

    out.ok = true; out.baselineDropped = false;
    out.f0 = f0;   out.f0Err = err(1);
    out.fwhm = W;  out.fwhmErr = err(2);
    out.A = A;     out.AErr = err(3);
    out.c = c;     out.cErr = 0; if opts.baseline; out.cErr = err(4); end
    out.cov = C;
    out.chi2Reduced = chi2r; out.birge = birge; out.dof = dof;

    out.Q = f0/W;
    % Q = f0/W with both correlated, so the ratio needs the cross term, not just the two
    % variances added in quadrature.
    if all(isfinite([C(1,1) C(2,2) C(1,2)]))
        relVar = C(1,1)/f0^2 + C(2,2)/W^2 - 2*C(1,2)/(f0*W);
        out.QErr = out.Q * sqrt(max(relVar, 0));
    else
        out.QErr = NaN;
    end

    peak = A + c;
    out.wMax = sqrt(max(peak, 0));
    if opts.baseline && all(isfinite([C(3,3) C(4,4) C(3,4)]))
        vPeak = C(3,3) + C(4,4) + 2*C(3,4);
    elseif isfinite(C(3,3))
        vPeak = C(3,3);
    else
        vPeak = NaN;
    end
    out.wMaxErr = sqrt(max(vPeak, 0)) / max(2*out.wMax, eps);
    out.piMinNs = 500/max(out.wMax, eps);                       % 1/(2f) in us -> ns
    out.piMinErrNs = out.piMinNs * out.wMaxErr/max(out.wMax, eps);
    out.contrastOfFit = A/max(peak, eps);

    resid = (y - evalModel(f, p, opts.baseline))./sy;
    out.resid = resid;
    ssRes = sum((y - evalModel(f, p, opts.baseline)).^2);
    ssTot = sum((y - mean(y)).^2);
    out.r2 = 1 - ssRes/max(ssTot, eps);
    out.model = @(ff) evalModel(ff(:), p, opts.baseline);

    out.centreInRange = (f0 >= min(f)) && (f0 <= max(f));
    out.resolved      = W <= 1.5*span;
    out.pointsBelow   = sum(f < f0);
    out.pointsAbove   = sum(f > f0);
    out.span          = span;

    % ---- is this fit worth quoting? --------------------------------------
    % THREE SEPARATE QUESTIONS, kept separate because they have different answers:
    %
    %   ok        the arithmetic worked and no parameter is pinned or absurd
    %   usable    ok, AND the data actually constrain the peak
    %   pedestal  the fit put most of the drive in the constant and made the "resonance"
    %             a small ripple on top of it
    %
    % A short scan across the top of a broad response is the case that needs pedestal:
    % over a 9 MHz span a real 30 MHz-wide Lorentzian is nearly a parabola, and a NARROW
    % Lorentzian sitting on a big floor fits that parabola just as well. Tested on
    % synthetic data with the truth known, that alternative won: it returned FWHM 5.5 MHz
    % against a true 30, with a fictitious 18.7 MHz of "direct coupling" holding the
    % pedestal up. The centre survived, so the peak POSITION is still recoverable from
    % such a scan -- but Q and the on-resonance Rabi frequency are not, and nothing in
    % chi2 or R^2 says so.
    out.pinnedWidth = (W <= 1.02*Wmin) || (W >= 0.98*Wmax);
    out.absurd      = (A > 50*max(y)) || ~isfinite(A) || ~isfinite(W);
    out.pedestal    = opts.baseline && (c > A);
    out.negFloor    = opts.baseline && (c < -0.25*max(y));

    out.ok = ~(out.pinnedWidth || out.absurd);
    if ~out.ok
        if out.pinnedWidth
            out.why = sprintf(['the width ran to the edge of its allowed range (%.2f MHz, ' ...
                'bounds %.2f to %.0f). These points do not contain a peak -- they are an ' ...
                'arc, and any width fits it.'], W, Wmin, Wmax);
        else
            out.why = sprintf('the fitted height (%.3g) is not physical against data peaking at %.3g', A, max(y));
        end
    end
    out.usable = out.ok && out.centreInRange && out.resolved && ~out.pedestal && ...
                 ~out.negFloor && out.pointsBelow >= 2 && out.pointsAbove >= 2;

    % ---- the independent cross-check: a parabola in the same variables ---
    % A Lorentzian near its top IS a parabola -- expand it and the quadratic term is
    % -4A/W^2 -- so a weighted quadratic in Omega^2 against frequency finds the same
    % centre with two fewer parameters and no width to be degenerate with. It cannot give
    % Q or the on-resonance amplitude, and it is biased once the wings are included, so it
    % is not the answer. It is the check on the answer.
    %
    % It earns its place on the short scans. On synthetic data with a true centre of
    % 2823.0 the Lorentzian was REJECTED at a 12 MHz span (its width pinned at the bound)
    % while its centre was still 2823.16 -- the peak position survives a scan that cannot
    % determine anything else. The parabola recovers that case without needing the
    % rejected fit to be trusted, which is why calibB1Rabi can still hand you a field
    % setting from a scan too short for Q.
    out.parabola = parabolaVertex(f, y, sy);
    % Only compare when the vertex is INSIDE the scan. A parabola through an arc puts its
    % vertex hundreds of MHz away with an error bar to match, and "0.2 sigma agreement"
    % between two numbers 350 MHz apart is arithmetically true and completely useless.
    if out.parabola.ok && out.parabola.inRange && isfinite(out.f0Err)
        sep = abs(out.f0 - out.parabola.f0);
        comb = sqrt(out.f0Err^2 + out.parabola.f0Err^2);
        out.crossCheckSigma = sep/max(comb, eps);
    else
        out.crossCheckSigma = NaN;
    end

    % ---- automatic retry without the floor ------------------------------
    % A pedestal fit has an obvious alternative -- the same model with c forced to zero --
    % and trying it costs one more grid scan. If that one is clean, it is the answer, and
    % saying so here is better than printing a warning and leaving you to re-run by hand.
    if out.pedestal && opts.retryNoBaseline
        o2 = opts; o2.baseline = false; o2.retryNoBaseline = false; o2.verbose = false;
        alt = resonatorLorentzianFit(f, w, wErr, o2);
        if isstruct(alt) && isfield(alt, 'usable') && alt.usable
            alt.baselineDropped = true;
            alt.pedestalFit = rmfield(out, 'opts');
            if opts.verbose
                fprintf(['\n  The fit with a free floor came out as a narrow peak on a large\n' ...
                         '  pedestal (%.0f %% of the drive in the floor), which is the known\n' ...
                         '  degeneracy of a short scan. Refitted with the floor forced to zero;\n' ...
                         '  that fit is clean and is the one reported. out.pedestalFit holds the\n' ...
                         '  other one if you want to compare.\n'], 100*(1 - out.contrastOfFit));
                report(alt);
            end
            out = alt; return
        end
    end

    if opts.verbose; report(out); end
end

% ======================================================================

function p = parabolaVertex(f, y, sy)
% Weighted quadratic in Omega^2 vs frequency, vertex reported with its error.
%
% Centred on the mean frequency before fitting. Uncentred, the design matrix has columns
% of order 1, 2.8e3 and 8e6, and the normal equations lose most of their significant
% figures -- the vertex of a parabola through three points at 2820 MHz is a small
% difference of large numbers. Centring makes it conditioned.
    p = struct('ok', false, 'f0', NaN, 'f0Err', NaN, 'why', '', 'curvature', NaN);
    if numel(f) < 4; p.why = 'fewer than 4 points'; return; end
    fbar = mean(f);
    x = f - fbar;
    X = [x.^2, x, ones(size(x))];
    Xw = X ./ sy; yw = y ./ sy;
    b = Xw \ yw;
    if any(~isfinite(b)); p.why = 'the weighted solve failed'; return; end
    p.curvature = b(1);
    if b(1) >= 0
        p.why = 'the parabola opens UPWARD -- these points have a minimum, not a maximum';
        return
    end
    dx = -b(2)/(2*b(1));
    p.f0 = fbar + dx;

    dof = max(numel(f) - 3, 1);
    r = Xw*b - yw;
    p.chi2Reduced = sum(r.^2)/dof;
    birge = max(sqrt(p.chi2Reduced), 0.05);
    try
        Cb = inv(Xw'*Xw) * birge^2;
        % dx = -b2/(2 b1), so the vertex error needs the b1-b2 covariance: a steeper
        % parabola and a shifted one are strongly correlated, and adding the two variances
        % in quadrature would overstate the error by a factor of a few.
        g = [b(2)/(2*b(1)^2), -1/(2*b(1)), 0];
        p.f0Err = sqrt(max(g*Cb*g', 0));
    catch
        p.f0Err = NaN;
    end
    p.ok = isfinite(p.f0);
    p.inRange = p.ok && p.f0 >= min(f) && p.f0 <= max(f);
    if p.ok && ~p.inRange
        p.why = 'the vertex is outside the measured range';
    end
end

function best = gridScan(f0g, Wg, f, y, sy, useBaseline)
% The exhaustive scan, one width at a time and every centre at once.
%
% VECTORISED, and not as a flourish. The obvious nested loop calls a 2-column backslash
% 601*81 = 48681 times per fit, which took long enough that the 120-realisation error-bar
% test could not finish. The linear step is a 2x2 weighted least squares, so it has a
% closed form in five sums, and those sums are column reductions over an
% nPoints-by-nCentres matrix. Same arithmetic, same answer, about three orders of magnitude
% faster -- which is what makes it affordable to VERIFY the error bars by Monte Carlo
% instead of asserting them.
%
% For columns X = [g 1] with weights 1/sy, writing u = g/sy, v = 1/sy, z = y/sy:
%
%       A = (Suz*Svv - Svz*Suv)/det,   c = (Svz*Suu - Suz*Suv)/det,   det = Suu*Svv - Suv^2
%       chi2 = Szz - A*Suz - c*Svz
%
% the last line being the standard result that at the least-squares solution the residual
% sum is z'z - beta'X'z.
    best = struct('cost', inf, 'f0', NaN, 'W', NaN, 'lin', []);
    f = f(:); y = y(:); sy = sy(:);
    v   = 1 ./ sy;
    z   = y ./ sy;
    Svv = sum(v.^2);
    Szz = sum(z.^2);
    Svz = sum(v.*z);
    f0r = f0g(:)';                      % row, so f - f0r broadcasts to nPoints x nCentres

    for iW = 1:numel(Wg)
        W = Wg(iW);
        G = 1 ./ (1 + 4*(f - f0r).^2 / W^2);
        U = G .* v;                                     % nPoints x nCentres
        Suu = sum(U.^2, 1);
        Suz = sum(U .* z, 1);
        if useBaseline
            Suv = sum(U .* v, 1);
            det = Suu*Svv - Suv.^2;
            bad = ~isfinite(det) | abs(det) < eps*max(abs(Suu*Svv), eps);
            A = (Suz*Svv - Svz*Suv) ./ det;
            c = (Svz*Suu - Suz.*Suv) ./ det;
            cost = Szz - A.*Suz - c.*Svz;
            cost(bad | A <= 0) = inf;                   % A <= 0 is a dip, not a resonance
        else
            A = Suz ./ Suu;
            c = zeros(size(A));
            cost = Szz - A.*Suz;
            cost(~isfinite(Suu) | Suu <= 0 | A <= 0) = inf;
        end
        [cMin, iF] = min(cost);
        if isfinite(cMin) && cMin < best.cost
            best.cost = cMin; best.f0 = f0r(iF); best.W = W;
            if useBaseline; best.lin = [A(iF); c(iF)]; else; best.lin = A(iF); end
        end
    end
end

function cost = costBounded(p, f, y, sy, useBaseline, Wmin, Wmax)
% costAt with the width held inside the grid's range. Outside it the cost is Inf rather
% than a penalty that grows: fminsearch handles a hard wall by reflecting off it, and a
% soft penalty would let the width sit just outside while the reported value is clamped.
    W = exp(p(2));
    if ~isfinite(W) || W < Wmin || W > Wmax; cost = inf; return; end
    cost = costAt(p(1), W, f, y, sy, useBaseline);
end

function [cost, lin] = costAt(f0, W, f, y, sy, useBaseline)
% Weighted chi2 at fixed (f0, W), with the linear parameters solved exactly.
    lin = [];
    if ~isfinite(f0) || ~isfinite(W) || W <= 0; cost = inf; return; end
    g = 1 ./ (1 + 4*(f - f0).^2 / W^2);
    if useBaseline; X = [g, ones(size(g))]; else; X = g; end
    Xw = X ./ sy; yw = y ./ sy;
    lastwarn('');
    lin = Xw \ yw;
    if any(~isfinite(lin)); cost = inf; lin = []; return; end
    % A <= 0 is not a Lorentzian resonance, it is a dip. Rejecting it here rather than
    % letting the search wander into it is what keeps the grid answer physical.
    if lin(1) <= 0; cost = inf; return; end
    r = Xw*lin - yw;
    cost = sum(r.^2);
end

function m = evalModel(f, p, useBaseline)
    f0 = p(1); W = p(2); A = p(3);
    if useBaseline; c = p(4); else; c = 0; end
    m = A ./ (1 + 4*(f - f0).^2 / max(W, eps)^2) + c;
end

function report(o)
    fprintf('\n--- resonator response, from the Rabi frequencies ---\n');
    if ~o.ok
        fprintf(['  THE FIT DID NOT PRODUCE A USABLE RESONANCE.\n  %s\n' ...
                 '  Nothing below is worth quoting; it is printed only so you can see what\n' ...
                 '  the search settled on. The field scan needs to be extended.\n'], o.why);
    end
    fprintf('  centre f0        %.2f +/- %.2f MHz\n', o.f0, o.f0Err);
    fprintf('  power FWHM       %.2f +/- %.2f MHz     -> Q = %.0f +/- %.0f\n', ...
            o.fwhm, o.fwhmErr, o.Q, o.QErr);
    fprintf('  Rabi at f0       %.3f +/- %.3f MHz     -> pi pulse %.1f +/- %.1f ns\n', ...
            o.wMax, o.wMaxErr, o.piMinNs, o.piMinErrNs);
    if o.opts.baseline
        fprintf('  through-resonator share %.1f %%  (floor c = %.3f MHz^2, i.e. %.2f MHz of Rabi\n', ...
                100*o.contrastOfFit, o.c, sqrt(max(o.c, 0)));
        fprintf('                           off resonance from direct coupling)\n');
        if o.pedestal
            fprintf(['    *** MOST OF THE DRIVE (%.0f %%) WENT INTO THE FLOOR, not the resonance.\n' ...
                     '    Over a short span a narrow peak on a big pedestal imitates a broad peak\n' ...
                     '    exactly, and the two cannot be told apart by chi2. The centre is\n' ...
                     '    probably still right; the width, Q and the on-resonance Rabi frequency\n' ...
                     '    are not. Compare against baseline = false.\n'], ...
                    100*(1 - o.contrastOfFit));
        end
        if o.c < 0 && abs(o.c) > 2*o.cErr
            fprintf(['    *** the floor came out NEGATIVE by %.1f sigma. A Lorentzian plus a\n' ...
                     '    negative constant is not a physical response; it usually means the\n' ...
                     '    shape is not Lorentzian over this span (two coupled modes, or the\n' ...
                     '    amplifier response tilting across the scan). Try baseline = false\n' ...
                     '    and compare f0.\n'], abs(o.c)/max(o.cErr, eps));
        end
    end
    fprintf('  R^2 %.4f, reduced chi2 %.2f on %d dof (errors scaled %.2fx)\n', ...
            o.r2, o.chi2Reduced, o.dof, o.birge);
    fprintf('  %d points below f0, %d above, over a %.1f MHz span\n', ...
            o.pointsBelow, o.pointsAbove, o.span);

    if isfield(o, 'parabola') && o.parabola.ok && o.parabola.inRange
        fprintf('  cross-check: parabola vertex at %.2f +/- %.2f MHz', ...
                o.parabola.f0, o.parabola.f0Err);
        if isfinite(o.crossCheckSigma)
            fprintf('   (%.1f sigma from the Lorentzian)\n', o.crossCheckSigma);
        else
            fprintf('\n');
        end
    elseif isfield(o, 'parabola') && o.parabola.ok
        fprintf(['  cross-check: the parabola puts its vertex at %.0f MHz, outside the scan --\n' ...
                 '  so it has no opinion here and is not compared. Both fits are being asked\n' ...
                 '  to extrapolate past the end of the data.\n'], o.parabola.f0);
        if isfinite(o.crossCheckSigma) && o.crossCheckSigma > 3
            fprintf(['    *** the two centres disagree by more than 3 sigma. The parabola is\n' ...
                     '    biased when the wings are included, so a wide clean scan can disagree\n' ...
                     '    honestly -- but on a scan that is not wide they should agree, and a\n' ...
                     '    disagreement then means the response is not a single Lorentzian.\n']);
        end
    elseif isfield(o, 'parabola')
        fprintf('  cross-check: no parabola -- %s\n', o.parabola.why);
    end

    if ~o.centreInRange
        fprintf(['\n  *** THE CENTRE IS OUTSIDE THE MEASURED RANGE (%.1f MHz, scan covered\n' ...
                 '  %.1f to %.1f). This is not a measurement of f0, it is an extrapolation\n' ...
                 '  off the end of an arc, and the quoted error bar does not cover being on\n' ...
                 '  the wrong side of the peak. EXTEND THE FIELD SCAN in that direction and\n' ...
                 '  re-run before using any of these numbers.\n'], ...
                o.f0, min(o.f), max(o.f));
    elseif o.pointsBelow < 2 || o.pointsAbove < 2
        fprintf(['\n  *** ONLY %d point(s) below the centre and %d above. The centre is inside\n' ...
                 '  the scan but barely bracketed, so f0 rests on one or two points and their\n' ...
                 '  noise. Add points on the thin side.\n'], o.pointsBelow, o.pointsAbove);
    end
    if ~o.resolved
        fprintf(['\n  *** THE FITTED FWHM (%.0f MHz) IS WIDER THAN THE SPAN MEASURED (%.0f MHz).\n' ...
                 '  Over a short arc a Lorentzian shoulder is indistinguishable from a\n' ...
                 '  parabola, so (f0, W, A) trade off against each other and a whole family\n' ...
                 '  of widths fits equally well. Q is meaningless here and f0 is much less\n' ...
                 '  certain than its error bar suggests. Widen the field scan, or accept the\n' ...
                 '  peak position and ignore Q.\n'], o.fwhm, o.span);
    end
    if o.chi2Reduced > 4
        fprintf(['\n  Reduced chi2 %.1f: the points scatter well beyond their error bars. Either\n' ...
                 '  the response is not a single Lorentzian, or the Rabi errors are\n' ...
                 '  underestimated, or the drive level moved during the scan. Check that the\n' ...
                 '  amplitude was identical at every field.\n'], o.chi2Reduced);
    end
end

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d);
    given = fieldnames(o);
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
