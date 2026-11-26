function out = esrFitRobust(args)
%ESRFITROBUST The ESR fit: lines, field, resonance, hyperfine. ONE STRUCT IN, ONE STRUCT OUT.
%
%   out = esrFitRobust(struct('exp', a))                   % a live AWGExpESR / AWGExpESRPulsed
%   out = esrFitRobust(struct('f', f, 'y', y, 'e', e))     % plain vectors (e optional)
%   out = esrFitRobust(struct('file', 'G:\...\_AutoSave\AWG ESR Pulsed_20261110_172200.mat'))
%   out = esrFitRobust(struct('exp', a, 'nLines', 1))      % force the number of NV lines
%   yy  = out.modelAt(ff);                                 % the fitted curve, any ff
%
% THE ONE ESR FITTER ON THIS SETUP. It absorbed esrMultiLineFit (the N-line engine),
% esrHyperfineFit (the nitrogen multiplet), esrLineFit, and the line-finding half of
% calibFieldSurvey (its mirror, split and physics-seeded refits). The run script, every
% calibration and the diagnostics all come here, so there is one definition of "the line is
% at f" and one place to fix it.
%
% One call gives, as far as the sweep supports it:
%
%   out.lines       every NV line, fitted TOGETHER on one shared baseline, with the full
%                   centre covariance (a splitting's variance is C11 + C22 - 2*C12)
%   out.field       |B|, the tilt of B from the ALIGNED <111> axis, D, and which line
%                   belongs to which orientation -- from the exact spin Hamiltonian
%   out.resonance   the microwave resonator's own dip, when the sweep shows one, fitted as
%                   its OWN component and kept out of the NV lines and the field
%   out.hyperfine   the 14N triplet / 15N doublet, when the data can resolve it
%
% and one figure with all of it on it. out.ok says whether to believe any of it. The old
% positional outputs are still named fields -- out.f0, out.width, out.contrast, out.SNR,
% out.coeffError -- and out.f is the FREQUENCY AXIS, not a line centre.
%
% ===================================================================
% THE MODEL
% ===================================================================
%
%   y(f) = P(f) - SUM_k A_k L(f; c_k, w_k)  [ - A_r L(f; c_r, w_r) ]
%
% P is a polynomial baseline in (f - xref), ascending powers, order baselineOrder. By default
% it follows the span: 2 for a sweep of 60 MHz or more (over 200 MHz the resonator envelope
% curves, and order 1 leaves that curvature to be soaked up by the line centres), 1 below
% that (in a narrow window a quadratic trades off against the line itself: on a resolved
% 14N triplet in 24 MHz the line's depth error grew until the line was dropped). L is a
% unit-depth Lorentzian. The sum is the NV lines; the bracket is the resonance, when there
% is one.
%
% Weighted Levenberg-Marquardt with an analytic Jacobian and Marquardt column scaling, every
% parameter CLAMPED to its range on every step -- a centre once escaped to -22,034,007 MHz
% because a line far outside the data costs nothing. Start points come from the deficit
% under a robust upper envelope, so a shallow line is not missed because the deep ones
% dragged the baseline down. A component that fits nothing -- shallower than
% minDepthSigma, narrower than the sampling, or wider than any NV line here (FWHM above
% min(span, maxFwhmMHz)) -- is dropped and the rest refitted without it.
%
% ===================================================================
% THE RESONANCE
% ===================================================================
%
% The microwave resonator puts a dip of its own into the ratio: broad (~30 MHz FWHM seen),
% not an NV transition, and NOT AT A FIXED FREQUENCY -- it sat near 2820 MHz and has been
% measured near 2763 since, so nothing here assumes where it is. At zero field it is the
% second dip that used to throw the fit: taken for an NV line it bent the baseline, stole
% the slot when nLines was forced, and once paired with the real line at D to report
% |B| = 8.26 G with the field off.
%
% It is recognised by what it is NOT:
%
%   not an NV line   with >= 3 lines, the field solution has no transition there; with
%                    fewer, it neither sits at D nor pairs with another line about D
%   not NV-narrow    its FWHM is at least resonanceWidthRatio (1.3) x the NV lines'
%
% and is then refitted as its own component, with a wider width range and a centre allowed
% up to 15% of the span OUTSIDE the sweep -- a resonance just below the first point still
% bends the edge. If no fitted line qualifies, the residual is smoothed and its deepest
% broad deficit is tried as the resonance, kept only if BIC improves by more than 10 AND the
% result passes the same two tests -- so a weak NV line the finder missed is not mistaken
% for it. Only on sweeps wider than resonanceMinSpanMHz (80): in a narrow window a
% resonance cannot be told from the baseline. resonance = false switches it off; a number
% (MHz) seeds it there when you already know.
%
% ===================================================================
% THE FIELD
% ===================================================================
%
%   >= 3 NV lines   esrFieldSolve fits Bx, By, Bz and D to every line through the exact
%                   Hamiltonian -- second-order transverse shift included, no pairing
%                   assumed. field.tiltDeg is the angle between B and the nearest <111>
%                   axis: the ALIGNED orientation, the one Rabi and Ramsey use, and
%                   field.lineAxis says which orientation each line belongs to. Below ~2.5
%                   degrees the off-axis lines merge and their positions stop carrying the
%                   tilt; calibFieldAlign's span metric takes over there.
%   2 lines about D one pair: B_par = splitting / (2 gamma). |B| >= B_par, direction open.
%   1 line at D     zero field (any splitting is under the linewidth).
%
% A wide sweep across D is REFINED first, as calibFieldSurvey used to: mirror-seed every
% line about D, split a line much wider than the narrowest, and seed the solved field's
% predicted positions -- each refit kept only if it gains a line, narrows the widest blend,
% or drops chi2 materially. That is what resolves the 44 MHz off-axis blob at 45 G.
%
% ===================================================================
% THE HYPERFINE
% ===================================================================
%
% Tried automatically when the data can show it: ONE NV line and a step <= 0.6 MHz -- a
% pulsed sweep (CW sweeps here step 1 to 2 MHz and never qualify). The one-line FWHM is NOT
% a gate: fitted to a resolved 14N triplet it measures the envelope, about 2A plus the
% pulse's sidelobes -- 7.1 MHz on three ~1.4 MHz lines at tau = 0.572 us, which a 3A =
% 6.5 MHz ceiling once skipped. A broad line on a fine grid comes back "1 dip, not
% resolved" and the figure is unchanged; on a wide fine grid the test costs ~25 s (301
% points), against 1 to 4 s for a 100-point pulsed sweep. The model is one multiplet, shared
% width and shape, EXACTLY equal spacing (the quadrupole term cancels in the transition):
%
%   3 dips  14N, m = -1 0 +1    A ~ 2.16 MHz       2 dips  15N, m = -1/2 +1/2   A ~ 3.03 MHz
%   1 dip   not resolved                           0       no line
%
% Every count is fitted with lorentz and gauss shapes -- plus, when the pulse length is
% known, the square-pulse response convolved with a Gaussian -- with equal and with free
% depths, and BIC picks. A is FITTED inside hfSplitRange and only then compared with the
% nitrogen values, so the 2-versus-3 decision is honest. out.lines stays one line per NV
% transition either way; the multiplet is out.hyperfine, and the figure draws it.
% hyperfine = true forces the test (still one line: crop with range), false skips it.
%
% ===================================================================
% ERRORS
% ===================================================================
%
% Two-sided Birge: parameter errors scaled by max(sqrt(reduced chi2), 0.05) in BOTH
% directions. AWGExpESR's sterr is drift between averages, common to every point, so it has
% the right SHAPE and the wrong SCALE (measured at 42x the point-to-point scatter). The
% error bars themselves come from esrPullData -- the larger of the propagated error and the
% scatter across averages, nothing discarded -- and NOT from extract1Ddata, whose 5, 3 and
% ONE sigma cuts throw away a third of good data and shrink every error bar downstream.
%
% ===================================================================
% OPTIONS (struct fields, all optional)
% ===================================================================
%   exp, file, f/y/e   the data: a live experiment (esrPullData), a saved struct or
%                  [f y e] matrix, a SaveLoad .mat, or plain vectors
%   range          [fmin fmax] MHz: fit only this part of the sweep
%   nLines         force the number of NV LINES (the resonance is not counted) ([] = find)
%   maxLines       cap when finding (default 10)
%   seedCentres    MHz, start from these centres instead of searching
%   baselineOrder  ([] = 2 for spans of 60 MHz or more, 1 below)
%   snr, minSepMHz candidate depth in noise units (4) and spacing (6 MHz)
%   minDepthSigma  drop a component shallower than this (default 3)
%   maxFwhmMHz     widest NV line (default 40)
%   field          'auto' | false                                        (default 'auto')
%   D0, dTol       D and how far it may move, MHz (2870, 15)
%   huberMHz       a line this far from every predicted transition is not NV (default 6)
%   refine         'auto' | true | false: the mirror/split/physics refits (default 'auto':
%                  sweeps of 60+ MHz across D0, nLines not forced)
%   pairTol, refineDepthSigma, splitRatio, splitGain, widthGain, physicsSeed,
%   verboseRefine  the refinement's knobs (3, 2.5, 1.8, 0.9, 0.75, true, false)
%   resonance      'auto' | true | false | MHz                           (default 'auto')
%   resonanceMinSpanMHz (80), resonanceMaxFwhmMHz ([] = the span),
%   resonanceWidthRatio (1.3), resonanceSnr (5)
%   hyperfine      'auto' | true | false                                 (default 'auto')
%   hfCount        force 0-3 dips ([] = choose)
%   hfShape        'auto' | 'lorentz' | 'gauss' | 'voigt' | 'pulse'      (default 'auto')
%   hfDepths       'auto' | 'equal' | 'free'                             (default 'auto')
%   hfSplitting    MHz, fix A instead of fitting it; hfSplitRange [1 5] MHz otherwise
%   tau            us, pulse length for 'pulse' (default: read from the experiment)
%   pulseRotationDeg   on-resonance rotation of that pulse (180)
%   avgOutlierMad  passed to esrPullData (4)
%   plotIt (true), verbose (true), showResid (false), yLabel, name, maxIter (200)
%
% ===================================================================
% OUTPUT
% ===================================================================
%   ok, why
%   lines          NV lines, by frequency: centre, centreErr, fwhm, fwhmErr, depth, depthErr
%   nLines, centreCov, splitting(+Err), pairMean(+Err)   all from the joint covariance
%   f0, width, contrast, coeffError, SNR, r2, r2Ceiling  the old positional outputs
%   resonance      found, centre, centreErr, fwhm, fwhmErr, depth, depthErr, outside, why
%   field          kind ('vector' | 'pair' | 'zero' | 'none' | 'off'), ok, B, tiltDeg, D,
%                  Bpar, projections, alignedAxisIdx, lineAxis, aligned, note, solve
%   B, tiltDeg, D, solve   the same, at the top level
%   hyperfine      tested, why, resolved, nDips, isotope, confidence, verdict, f0, splitting,
%                  fwhm, depth (+Err), shape, lines, prob, candidates, modelAt, warnings
%   model, modelAt, residual, baseline, xref, chi2Reduced, birge, refineLog
%   f, y, e, repeats, averages, dropped, source, tau, opts
%
% Base MATLAB: esrPullData for experiments, esrFieldSolve for the field, fitStyle for the
% figure. No Curve Fitting, Optimization or Statistics toolbox. Offline grader: testEsrFit.

    K = struct('gamma', 2.8025, ...     % MHz/G, NV electron (g = 2.0028)
               'A14N', 2.16, ...        % |A_par| of 14N, MHz
               'A15N', 3.03);           % |A_par| of 15N, MHz

    DEFAULTS = struct( ...
        'exp', [], 'file', '', 'f', [], 'y', [], 'e', [], 'range', [], ...
        'nLines', [], 'maxLines', 10, 'seedCentres', [], 'baselineOrder', NaN, ...
        'snr', 4, 'minSepMHz', 6, 'minDepthSigma', 3, 'maxFwhmMHz', 40, 'maxIter', 200, ...
        'field', 'auto', 'D0', 2870, 'dTol', 15, 'huberMHz', 6, ...
        'refine', 'auto', 'pairTol', 3, 'refineDepthSigma', 2.5, 'splitRatio', 1.8, ...
        'splitGain', 0.9, 'widthGain', 0.75, 'physicsSeed', true, 'verboseRefine', false, ...
        'resonance', 'auto', 'resonanceMinSpanMHz', 80, 'resonanceMaxFwhmMHz', [], ...
        'resonanceWidthRatio', 1.3, 'resonanceSnr', 5, ...
        'hyperfine', 'auto', 'hfCount', [], 'hfShape', 'auto', 'hfDepths', 'auto', ...
        'hfSplitting', [], 'hfSplitRange', [1 5], 'tau', [], 'pulseRotationDeg', 180, ...
        'avgOutlierMad', 4, 'plotIt', true, 'verbose', true, 'showResid', false, ...
        'yLabel', '', 'name', '');

    if ~exist('args', 'var') || isempty(args); args = struct(); end
    if ~isstruct(args)
        error(['esrFitRobust takes ONE struct: esrFitRobust(struct(''exp'', a)) or ' ...
               'esrFitRobust(struct(''f'', f, ''y'', y, ''e'', e)).']);
    end
    opts = withDefaults(args, DEFAULTS, 'esrFitRobust');
    checkOptions(opts);

    % Rank-deficient solves are expected and handled (a line pushed between samples, a
    % degenerate multiplet candidate), so their warnings are noise. Restored on exit. Each is
    % switched off on its own: newer MATLAB rejects an identifier it does not know, and the
    % Octave ones killed the old hyperfine fitter outright on R2025b.
    wst = warning;
    ids = {'MATLAB:rankDeficientMatrix', 'MATLAB:singularMatrix', ...
           'MATLAB:nearlySingularMatrix', 'Octave:singular-matrix', ...
           'Octave:nearly-singular-matrix'};
    for i = 1:numel(ids)
        try; warning('off', ids{i}); catch; end
    end
    restoreWarn = onCleanup(@() warning(wst));                                   %#ok<NASGU>

    out = emptyOut(opts);

    % ---- the data ----------------------------------------------------------
    try
        [x, y, e, info] = getData(opts);
    catch err
        out.why = err.message; out = failed(out); return
    end
    tau = info.tau;
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
    haveErr = any(isfinite(e));

    out.f = x; out.y = y; out.e = e;
    if ~haveErr; out.e = []; end
    out.repeats = info.reps; out.averages = info.nAvg; out.dropped = info.dropped;
    out.source = info.source; out.tau = tau;
    if ~isfinite(opts.baselineOrder)
        opts.baselineOrder = 1 + (numel(x) >= 2 && (max(x) - min(x)) >= 60);  % see THE MODEL
        out.opts = opts;
    end
    if numel(x) < opts.baselineOrder + 8
        out.why = sprintf('only %d usable points', numel(x)); out = failed(out); return
    end

    % Weights. Without error bars, a robust point-to-point noise stands in, so chi2 stays in
    % sigma units -- the resonance's BIC test and the hyperfine choice both assume that.
    % Birge rescales from the residuals either way.
    noise = pointNoise(y);
    ee = e;
    if haveErr
        ee(~isfinite(ee)) = median(ee(isfinite(ee)));
    else
        ee = noise * ones(size(x));
    end

    % ---- lines, and the resonance ----------------------------------------------
    [F, msgs] = fitAll(x, y, ee, noise, opts);
    out.refineLog = msgs;
    if opts.verbose && opts.verboseRefine && ~isempty(msgs)
        fprintf('  %s\n', msgs{:});
    end
    hf = [];
    if ~F.ok
        % A RESOLVED multiplet is exactly where the one-Lorentzian line search is weakest:
        % two 1 MHz lines sampled every 0.6 MHz are one or two points each, and the finder's
        % spike filter takes them for noise. The multiplet model does not need that search
        % (its start is a grid, and "no line" is one of its candidates), so where the
        % hyperfine is asked for -- or the step could resolve it -- it gets the sweep itself.
        [F2, hf2] = multipletFallback(x, y, ee, tau, opts, K);
        if F2.ok
            msgs{end+1} = sprintf(['the one-Lorentzian line search failed (%s); the multiplet ' ...
                                   'model found %s'], F.why, hf2.verdict);
            out.refineLog = msgs;
            F = F2; hf = hf2;
        end
    end
    if ~F.ok
        out.why = F.why; out = failed(out); return
    end

    % ---- the field, and the hyperfine ------------------------------------------
    fld = fieldFrom(F, x, opts, K);
    if isempty(hf)
        hf = hyperfineStage(x, y, ee, F, tau, opts, K);
    end

    % ---- assemble ----------------------------------------------------------------
    L = F.lines;
    out.lines = L;            out.nLines = F.nLines;
    out.centreCov = F.centreCov;
    out.splitting = F.splitting;  out.splittingErr = F.splittingErr;
    out.pairMean  = F.pairMean;   out.pairMeanErr  = F.pairMeanErr;
    out.baseline = F.baseline; out.xref = F.xref;
    out.chi2Reduced = F.chi2Reduced; out.birge = F.birge;
    out.residual = F.residual; out.model = F.model; out.modelAt = F.modelAt;
    out.resonance = F.resonance;
    out.field = fld; out.solve = fld.solve;
    out.B = fld.B; out.tiltDeg = fld.tiltDeg; out.D = fld.D;
    if hf.tested && hf.ok && ~isempty(hf.model)
        % The multiplet's own R^2 and ceiling over its window, for the title when it is drawn.
        ssH = max(sum((hf.y - mean(hf.y)).^2), eps);
        hf.r2 = 1 - sum((hf.y - hf.model).^2) / ssH;
        if haveErr
            inW = x >= hf.window(1) & x <= hf.window(2);
            hf.r2Ceiling = 1 - sum(ee(inW).^2) / ssH;
        end
    end
    out.hyperfine = hf;

    out.f0 = [L.centre]; out.width = [L.fwhm]; out.contrast = [L.depth];
    out.coeffError = [[L.depthErr]; [L.centreErr]; [L.fwhmErr]]';
    % Same definition as ESRfit's printed SNR, so the two stay comparable: depth over the
    % mean error scaled by the total shots. A figure of merit, not a fitted quantity.
    if haveErr
        out.SNR = out.contrast ./ (mean(ee) * sqrt(max(info.nAvg, 1) * max(info.reps, 1)));
    else
        out.SNR = nan(size(out.contrast));
    end
    ss = max(sum((y - mean(y)).^2), eps);
    out.r2 = 1 - sum((y - F.model).^2) / ss;
    % The ceiling the noise imposes: the best R^2 any correct model could reach on this
    % data. 0.82 against a ceiling of 0.84 is as good as the photons allow; the same 0.82
    % against 0.99 is a real failure, and one number cannot tell those apart.
    if haveErr
        out.r2Ceiling = 1 - sum(ee.^2) / ss;
    else
        out.r2Ceiling = NaN;
    end
    out.ok = true;

    if opts.verbose; report(out, K); end
    if opts.plotIt;  drawIt(out); end
end

% ======================================================================
% THE LINE FIT
% ======================================================================

function [F, msgs] = fitAll(x, y, ee, noise, opts)
% Lines first, then the resonance, then the forced count, then the refinement.
    msgs = {};
    span = max(x) - min(x);
    [resMode, resHint] = resonanceMode(opts.resonance);
    resOn = ~strcmp(resMode, 'off') && (strcmp(resMode, 'on') || span >= opts.resonanceMinSpanMHz);

    [c0, a0, g0, env] = findLines(x, y, opts, noise);
    if ~isempty(opts.seedCentres)
        cs = sort(opts.seedCentres(:));
        cs = cs(cs >= min(x) & cs <= max(x));
        if ~isempty(cs); [c0, a0, g0] = seedStarts(x, y, cs, env, noise); end
    end
    res0 = [];
    if resOn && ~isempty(resHint)
        res0 = resonanceSeedAt(x, y, env, noise, resHint, 15);
        near = abs(c0 - resHint) < 10;           % that dip IS the resonance, not a line
        c0(near) = []; a0(near) = []; g0(near) = [];
    end
    if ~resOn && ~isempty(opts.nLines)
        % Nothing to separate out, so the count can be applied to the candidates directly.
        [c0, a0, g0] = forceStarts(c0, a0, g0, opts.nLines, x, noise);
    end
    if isempty(c0)
        F = blankFit();
        F.why = sprintf(['no line found deeper than %g noise units (noise %.2g). Either ' ...
                         'there is no ESR line in this sweep or the drive is off.'], ...
                        opts.snr, noise);
        return
    end

    F = coreFit(x, y, ee, starts(c0, a0, g0), res0, opts, opts.minDepthSigma);
    if ~F.ok; return; end

    if resOn
        if ~F.resonance.found
            [F, msgs] = classifyResonance(F, x, y, ee, opts, msgs);
        end
        if ~F.resonance.found
            [F, msgs] = residualResonance(F, x, y, ee, noise, opts, msgs, resHint);
        end
        if ~isempty(opts.nLines) && F.nLines ~= opts.nLines
            % Counted AFTER the resonance is out, so nLines = 1 at zero field means the
            % line at D -- not whichever of it and the resonance happened to be deeper.
            F = forceCountFit(F, x, y, ee, noise, opts);
            if ~F.ok; return; end
        end
    end

    if refineApplies(F, x, opts)
        [F, msgs] = refineLines(F, x, y, ee, noise, env, opts, msgs);
        if resOn && ~F.resonance.found
            [F, msgs] = classifyResonance(F, x, y, ee, opts, msgs);
        end
    end
end

function F = coreFit(x, y, ee, nv, res, opts, minDepthSigma, testRes)
% One joint fit: NV lines (nv.c/a/g starts) plus, optionally, the resonance (res.c/a/g).
% Components that fail their checks are dropped and the rest refitted -- one bad line
% distorts every other parameter through the shared baseline.
%
% testRes (default true) applies the in-sweep depth test to the resonance. The callers that
% decide whether there IS a resonance pass false and compare the fit with and without it
% instead (dBIC > 10): a resonance centred outside the sweep is a degenerate set of
% parameters, and a significance test on them can fail while the dip in the data is plain.
% Refits that carry an already-accepted resonance pass false too, so it is not lost to that
% same degeneracy halfway through the refinement.
    if ~exist('testRes', 'var'); testRes = true; end
    F = blankFit();
    nB = opts.baselineOrder + 1;
    xLo = min(x); xHi = max(x); span = xHi - xLo; dx = median(diff(x));
    F.xref = mean(x);
    if numel(x) < nB + 6
        F.why = sprintf('only %d usable points', numel(x)); return
    end
    % The resonance's own ranges: centre up to 15% of the span outside the sweep, width up
    % to the span (or resonanceMaxFwhmMHz). Its depth test is the guard against it simply
    % duplicating the baseline: a component the polynomial can stand in for comes out with
    % a huge depth error and fails it.
    resMax = opts.resonanceMaxFwhmMHz;
    if isempty(resMax); resMax = span; end
    rcLo = xLo - 0.15*span; rcHi = xHi + 0.15*span;
    rgLo = max(dx, 1);      rgHi = max(resMax/2, 2*rgLo);

    c0 = nv.c(:); a0 = nv.a(:); g0 = nv.g(:);
    haveRes = ~isempty(res);
    resWhy = '';
    settled = false;
    for attempt = 1:(numel(c0) + 3)                                         %#ok<FXUP>
        nNV = numel(c0);
        if nNV == 0
            F.why = 'every candidate line failed its depth or width test'; return
        end
        cs = c0; as = a0; gs = g0;
        bd = struct('cLo', repmat(xLo, nNV, 1), 'cHi', repmat(xHi, nNV, 1), ...
                    'gLo', repmat(0.5*dx, nNV, 1), 'gHi', repmat(span/2, nNV, 1));
        if haveRes
            cs(end+1, 1) = min(max(res.c, rcLo), rcHi);                     %#ok<AGROW>
            as(end+1, 1) = max(res.a, eps);                                 %#ok<AGROW>
            gs(end+1, 1) = min(max(res.g, rgLo), rgHi);                     %#ok<AGROW>
            bd.cLo(end+1, 1) = rcLo; bd.cHi(end+1, 1) = rcHi;
            bd.gLo(end+1, 1) = rgLo; bd.gHi(end+1, 1) = rgHi;
        end
        [p, C, chi2r, chi2, birge, resid, yhat, ok, why] = ...
            gnFit(x, y, ee, cs, as, gs, nB, F.xref, bd, opts.maxIter);
        if ~ok; F.why = why; return; end
        [A, c, g] = unpack(p, nB);
        eA = sqrt(max(diag(C(nB+1:3:end, nB+1:3:end)), 0));
        if haveRes
            % Judged on the dip it puts INTO THE SWEEP, not on its fitted peak. A resonance
            % centred just past the first point shows only its flank, so its peak depth is an
            % extrapolation strongly correlated with its centre and width -- poorly determined
            % even when the dip in the data is plain (it failed a peak-depth test on half of
            % the noise draws of a 2763 MHz resonance on a 2770 MHz sweep). The depth at the
            % nearest swept point, with its error propagated through the covariance, is what
            % the data actually says. Centre inside the sweep: identical to the peak test.
            k = nNV + 1; j = nB + 3*nNV;
            xv = min(max(c(k), xLo), xHi);
            u = xv - c(k); dn = u^2 + g(k)^2;
            v = A(k) * g(k)^2 / dn;
            G = zeros(1, numel(p));
            G(j+1) = g(k)^2 / dn;                                  % d/dA
            G(j+2) = A(k) * 2*g(k)^2 * u / dn^2;                   % d/dc
            G(j+3) = A(k) * 2*g(k) * u^2 / dn^2;                   % d/dg
            ev = sqrt(max(G * C * G', 0));
            if ~(isfinite(v) && isfinite(c(k))) || (testRes && v < minDepthSigma*ev)
                resWhy = sprintf(['the resonance candidate (%.1f MHz, FWHM %.1f) is only %.1f ' ...
                    'sigma deep in the sweep'], c(k), 2*g(k), v/max(ev, eps));
                haveRes = false;
                continue
            end
        end
        % Reject: too shallow to believe, centre pushed onto the sweep edge, narrower than the
        % sampling resolves, or WIDER than an NV line can be. Too wide has two independent
        % limits, and a component is rejected if either holds: wider than the span (you
        % cannot measure a width you have not spanned) or wider than maxFwhmMHz (no NV line
        % here is that broad, saturated or not). It must NOT be a fraction of the span: at
        % 0.25*span a clean 10.2 MHz line in a 30 MHz window was thrown out as baseline.
        i = 1:nNV;
        tooWide = 2*g(i) > min(span, opts.maxFwhmMHz);
        bad = (A(i) < minDepthSigma*eA(i)) | (c(i) < xLo) | (c(i) > xHi) | ...
              (2*g(i) < 1.5*dx) | tooWide | ~isfinite(A(i)) | ~isfinite(c(i));
        if ~any(bad); settled = true; break; end
        if nNV == 1
            F.why = sprintf(['the only line left does not survive its checks (centre %.1f ' ...
                'MHz, depth %.3g, FWHM %.2f). The fit did not converge on anything usable.'], ...
                c(1), A(1), 2*g(1));
            return
        end
        drop = find(bad, 1);
        c0(drop) = []; a0(drop) = []; g0(drop) = [];
    end
    if ~settled
        F.why = 'the line set did not settle'; return
    end

    % ---- assemble ----------------------------------------------------------------
    se = sqrt(max(diag(C), 0));
    n = nNV;
    idxA = nB + (1:3:3*n); idxC = nB + (2:3:3*n); idxG = nB + (3:3:3*n);
    An = A(1:n); cn = c(1:n); gn = g(1:n);
    [cn, ord] = sort(cn);
    An = An(ord); gn = gn(ord);
    idxA = idxA(ord); idxC = idxC(ord); idxG = idxG(ord);
    L = emptyLines();
    for k = 1:n
        L(k) = struct('centre', cn(k), 'centreErr', se(idxC(k)), 'fwhm', 2*gn(k), ...
                      'fwhmErr', 2*se(idxG(k)), 'depth', An(k), 'depthErr', se(idxA(k)));
    end
    R = noResonance(resWhy);
    if haveRes
        k = n + 1; j = nB + 3*n;
        R = struct('found', true, 'centre', c(k), 'centreErr', se(j+2), ...
                   'fwhm', 2*g(k), 'fwhmErr', 2*se(j+3), 'depth', A(k), 'depthErr', se(j+1), ...
                   'outside', c(k) < xLo || c(k) > xHi, 'why', '');
    end
    F.lines = L; F.nLines = n;
    F.centreCov = C(idxC, idxC);        % the whole point: sums and differences of centres
    [F.splitting, F.splittingErr, F.pairMean, F.pairMeanErr] = pairTables(cn, F.centreCov);
    F.baseline = p(1:nB)';
    F.chi2Reduced = chi2r; F.chi2 = chi2; F.birge = birge;
    F.residual = resid; F.model = yhat; F.nPar = numel(p);
    F.resonance = R;
    % The model as a function, so a caller can draw it smooth or evaluate between points.
    % Baseline coefficients are ASCENDING powers of (f - xref), not polyval order, and every
    % component is SUBTRACTED, because they are dips.
    F.modelAt = @(xx) esrModelAt(xx, F.baseline, F.xref, L, R);

    % Final sanity, cheap and absolute: a depth above 1 in normalised fluorescence, or an NV
    % centre outside the sweep, means the fit diverged whatever chi2 says.
    if any(A > 1) || any(cn < xLo - 1e-6) || any(cn > xHi + 1e-6) || ...
       any(~isfinite([A(:); c(:); g(:)]))
        F.why = sprintf(['the fit diverged: depths up to %.3g, centres %.1f to %.1f MHz ' ...
            'against a sweep of %.1f to %.1f.'], max(A), min(cn), max(cn), xLo, xHi);
        return
    end
    F.ok = true;
end

function [p, C, chi2r, chi2, birge, resid, yhat, ok, why] = ...
            gnFit(x, y, ee, c0, a0, g0, nB, xr, bd, maxIter)
% Weighted Levenberg-Marquardt on the explicit model, analytic Jacobian. Hand-rolled rather
% than through fit(): the covariance between line centres is an output this function
% exists for, and confint does not provide it.
    ok = false; why = ''; C = []; chi2r = NaN; chi2 = NaN; birge = 1; resid = []; yhat = [];
    n = numel(c0);
    p = zeros(nB + 3*n, 1);
    b0 = [ones(numel(x), 1), (x - xr)] \ y;         % baseline start: a straight line
    p(1) = b0(1);
    if nB >= 2; p(2) = b0(2); end
    for k = 1:n
        p(nB + 3*k - 2) = max(a0(k), eps);
        p(nB + 3*k - 1) = c0(k);
        p(nB + 3*k    ) = max(g0(k), eps);
    end
    p = clampP(p, nB, bd);

    w = 1 ./ ee(:).^2;
    cost = @(pp) sum(w .* (y - modelOf(pp, x, nB, xr)).^2);
    fPrev = cost(p);
    lambda = 1e-3;
    for it = 1:maxIter                                                      %#ok<FXUP>
        [f, J] = modelAndJac(p, x, nB, xr);
        r = y - f;
        JW = J .* w;
        H = J' * JW;
        gr = JW' * r;
        % Marquardt's column scaling: centres are ~2870 MHz and depths ~0.01, so the raw
        % normal matrix ran to rcond 1e-37. Solve in scaled variables, then unscale.
        sc = sqrt(max(diag(H), eps));
        Ds = diag(1 ./ sc);
        Hs = Ds * H * Ds;
        gs = Ds * gr;
        stepTaken = false; converged = false;
        for tries = 1:12                                                    %#ok<FXUP>
            try
                dp = Ds * ((Hs + lambda*eye(size(Hs))) \ gs);
            catch
                dp = Ds * (pinv(Hs + lambda*eye(size(Hs))) * gs);
            end
            % EVERY parameter bounded on EVERY step. A centre left free escaped to
            % -22,034,007 MHz: outside the data its Lorentzian is zero, the cost drops, and
            % the step is accepted. Clamping beats a post-hoc check, because a diverged
            % parameter corrupts the shared baseline before anyone inspects it.
            pTry = clampP(p + dp, nB, bd);
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
    chi2 = sum((resid ./ ee(:)).^2);
    chi2r = chi2 / dof;
    % TWO-SIDED Birge, floored: the sterr gives the shape, the residuals give the scale.
    if isfinite(chi2r) && chi2r > 0
        birge = max(sqrt(chi2r), 0.05);
    end
    JW = J .* w;
    H = J' * JW;
    sc = sqrt(max(diag(H), eps));
    Ds = diag(1 ./ sc);
    Hs = Ds * H * Ds;
    if rcond(Hs) < 1e-12
        C = Ds * pinv(Hs) * Ds * birge^2;
    else
        C = Ds * (Hs \ eye(size(Hs))) * Ds * birge^2;
    end
    if ~all(isfinite(C(:)))
        why = 'covariance came out non-finite'; return
    end
    ok = true;
end

function p = clampP(p, nB, bd)
    iA = nB+1:3:numel(p); iC = nB+2:3:numel(p); iG = nB+3:3:numel(p);
    p(iA) = abs(p(iA));                                            % depth >= 0
    p(iG) = min(max(abs(p(iG)), bd.gLo), bd.gHi);                  % width in range
    p(iC) = min(max(p(iC), bd.cLo), bd.cHi);                       % centre in range
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

function v = esrModelAt(xx, baseline, xref, L, R)
% y(f) = P(f) - SUM_k A_k g_k^2/((f - c_k)^2 + g_k^2), P ascending in (f - xref), the
% resonance (when found) one more term of the same form.
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
    if R.found
        g = R.fwhm/2;
        v = v - R.depth * g^2 ./ ((xx - R.centre).^2 + g^2);
    end
end

function [S, SE, M, ME] = pairTables(c, C)
% Every pairwise splitting and pair mean, WITH the covariance term:
%   splitting  c_j - c_i        var = Cii + Cjj - 2*Cij
%   pair mean  (c_i + c_j)/2    var = (Cii + Cjj + 2*Cij)/4
% The centres are correlated through the shared baseline, usually strongly negatively for
% a well-determined pair, which makes the splitting MORE precise than either centre.
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

function [c0, a0, g0, env] = findLines(x, y, opts, noise)
% Candidate centres from the DEFICIT below a robust upper envelope. The envelope is fitted
% to the points NOT in a dip (start from all, keep those at or above the line, repeat): one
% pass over all points drags it into the dips and a shallow line then gets missed.
%
% The noise is point-to-point scatter, NOT the spread of the deficit: with several lines
% the deficit's spread IS the lines, the threshold lands above the shallow ones, and on a
% synthetic <111> pattern exactly the aligned lines -- the ones |B| comes from -- vanished.
    x = x(:); y = y(:);
    keep = true(size(x));
    b = [0; median(y)];
    for it = 1:6                                                            %#ok<FXUP>
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
    env = base;                         % SEEDED refits start from the same envelope

    ds = med3local(d);                  % kill single-point spikes before peak finding
    thresh = opts.snr * noise;
    isPeak = [false; ds(2:end-1) >= ds(1:end-2) & ds(2:end-1) >= ds(3:end); false];
    cand = find(isPeak & ds > thresh);
    c0 = []; a0 = []; g0 = [];
    if isempty(cand); return; end
    % Greedy by depth with a minimum separation, so one broad line does not register as three.
    [~, ord] = sort(ds(cand), 'descend');
    cand = cand(ord);
    chosen = [];
    for k = cand(:)'
        if all(abs(x(k) - x(chosen)) >= opts.minSepMHz)
            chosen(end+1) = k;                                               %#ok<AGROW>
        end
        if numel(chosen) >= opts.maxLines; break; end
    end
    chosen = sort(chosen);
    c0 = x(chosen);
    a0 = ds(chosen);
    % Width start from where the deficit falls to half, outward from each peak.
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

function [c0, a0, g0] = seedStarts(x, y, cs, env, noise)
% Start points AT caller-supplied centres. The depth comes from the GLOBAL envelope, not a
% local window: a 7-point running maximum spans 10 MHz on a 1.5 MHz grid, sits entirely
% inside a 10-44 MHz line, and started every seeded line at the noise floor -- which is why
% seeded refits once could not beat the unseeded fit. The width start is capped by the
% distance to the neighbouring seed, so two seeds inside one blend do not start on top of
% each other.
    cs = sort(cs(:));
    cs = cs(cs >= min(x) & cs <= max(x));
    c0 = cs;
    a0 = max(interp1(x, env - y, cs, 'linear', 'extrap'), 4*noise);
    step = median(diff(x));
    gap = inf(size(cs));
    if numel(cs) > 1
        d2 = diff(cs);
        gap(1:end-1) = min(gap(1:end-1), d2);
        gap(2:end)   = min(gap(2:end),   d2);
    end
    g0 = min(max(3*step, 2), max(gap/2, 2*step));
end

function [c0, a0, g0] = forceStarts(c0, a0, g0, n, x, noise)
% A specific count: keep the deepest n, or pad by spreading extra starts across the sweep.
    c0 = c0(:); a0 = a0(:); g0 = g0(:);
    if numel(c0) > n
        [~, ord] = sort(a0, 'descend');
        keep = sort(ord(1:n));
        c0 = c0(keep); a0 = a0(keep); g0 = g0(keep);
    elseif numel(c0) < n
        need = n - numel(c0);
        extra = linspace(min(x), max(x), need + 2);
        c0 = [c0; extra(2:end-1)'];
        a0 = [a0; repmat(max([a0; 4*noise]), need, 1)];
        g0 = [g0; repmat(max([median(g0), 5]), need, 1)];
    end
end

function F = forceCountFit(F, x, y, ee, noise, opts)
% nLines applied to the NV lines once the resonance has been separated. Lines that belong
% to the NV pattern are kept before lines that do not, then the deepest.
    L = F.lines;
    c = [L.centre]; a = [L.depth]; g = [L.fwhm]/2;
    if numel(L) > opts.nLines
        isNV = nvByD(c, opts);
        [~, ord] = sortrows([-double(isNV(:)), -a(:)]);
        keep = sort(ord(1:opts.nLines));
        c = c(keep); a = a(keep); g = g(keep);
    else
        [c, a, g] = forceStarts(c, a, g, opts.nLines, x, noise);
    end
    F = coreFit(x, y, ee, starts(c, a, g), resonanceAsSeed(F.resonance), opts, ...
                opts.minDepthSigma, false);
end

function s = starts(c, a, g)
    s = struct('c', c(:), 'a', a(:), 'g', g(:));
end

function y = med3local(x)
    y = x;
    for i = 2:numel(x)-1
        y(i) = median(x(i-1:i+1));
    end
end

function s = pointNoise(y)
% Successive differences cancel smooth structure, so median|diff| sees only noise wherever
% the spectrum is locally flat; for Gaussian noise median|diff| = 0.6745*sqrt(2)*sigma.
    s = 1.4826 * median(abs(diff(y))) / sqrt(2);
    if ~isfinite(s) || s <= 0; s = std(y); end
    if ~isfinite(s) || s <= 0; s = eps; end
end

% ======================================================================
% THE RESONANCE
% ======================================================================

function [F, msgs] = classifyResonance(F, x, y, ee, opts, msgs)
% Is one of the fitted lines the resonance? Convert it and refit if so.
    [k, how] = resonanceCandidate(F, x, opts);
    if k == 0; return; end
    L = F.lines;
    others = setdiff(1:numel(L), k);
    nv = starts([L(others).centre], [L(others).depth], [L(others).fwhm]/2);
    res = struct('c', L(k).centre, 'a', L(k).depth, 'g', L(k).fwhm/2);
    F2 = coreFit(x, y, ee, nv, res, opts, opts.minDepthSigma, false);
    % Is the dip really there? The same fit WITHOUT it, and BIC charging its 3 parameters.
    dBic = NaN;
    if F2.ok && F2.resonance.found
        F0 = coreFit(x, y, ee, starts([F2.lines.centre], [F2.lines.depth], ...
                     [F2.lines.fwhm]/2), [], opts, opts.minDepthSigma);
        dBic = Inf;
        if F0.ok
            n = numel(x);
            s2 = max(min(F0.chi2Reduced, F2.chi2Reduced), 0.05^2);
            dBic = (F0.chi2/s2 + F0.nPar*log(n)) - (F2.chi2/s2 + F2.nPar*log(n));
        end
    end
    if F2.ok && F2.resonance.found && F2.chi2Reduced <= 1.05*F.chi2Reduced && dBic > 10
        msgs{end+1} = sprintf(['resonance: the %.1f MHz line (FWHM %.1f) is not an NV line ' ...
                               '-- %s -- refitted as the resonance, now %.1f MHz (dBIC %.0f)'], ...
                              L(k).centre, L(k).fwhm, how, F2.resonance.centre, dBic);
        F = F2;
    else
        if ~F2.ok
            why = F2.why;
        elseif ~F2.resonance.found
            why = F2.resonance.why;
        elseif ~(dBic > 10)
            why = sprintf('dBIC only %.1f against no resonance at all', dBic);
        else
            why = sprintf('chi2 %.3g -> %.3g', F.chi2Reduced, F2.chi2Reduced);
        end
        msgs{end+1} = sprintf(['resonance: the %.1f MHz line looked like it (%s) but did ' ...
                               'not refit as one (%s) -- kept as a line'], L(k).centre, how, why);
    end
end

function [k, how] = resonanceCandidate(F, x, opts)
% The widest fitted line that belongs to no NV pattern and is broader than the NV lines.
% Returns 0 when there is no NV structure to compare against: a lone line, or lines none
% of which are NV-consistent, cannot be told apart from the resonance by position.
    k = 0; how = '';
    L = F.lines;
    if numel(L) < 2; return; end
    c = [L.centre]; ce = [L.centreErr]; w = [L.fwhm];
    [isNV, how0] = nvMembership(c, ce, x, opts);
    cand = find(~isNV);
    if isempty(cand) || ~any(isNV); return; end
    [~, j] = max(w(cand));
    kk = cand(j);
    wNV = median(w(isNV));
    if w(kk) < opts.resonanceWidthRatio * wNV; return; end
    k = kk;
    how = sprintf('%s, and %.1f MHz wide against %.1f for the NV lines', how0, w(kk), wNV);
end

function [isNV, how] = nvMembership(c, ce, x, opts)
% Which lines does NV physics account for? With 3+ lines the field solution decides; with
% fewer (or when it fails), a line is NV if it sits at D or pairs with another about D.
    isNV = true(1, numel(c));
    if numel(c) >= 3
        sol = solveField(c, ce, x, opts);
        if sol.ok
            u = sol.unmatched(:)';
            for i = 1:numel(c)
                isNV(i) = ~any(abs(c(i) - u) < 1e-9);
            end
            how = sprintf('the %.1f G field solution has no transition there', sol.B);
            return
        end
    end
    [isNV, how] = nvByD(c, opts);
end

function [isNV, how] = nvByD(c, opts)
    n = numel(c);
    isNV = false(1, n);
    for i = 1:n
        if abs(c(i) - opts.D0) <= opts.dTol; isNV(i) = true; continue; end
        others = c([1:i-1, i+1:n]);
        isNV(i) = any(abs((c(i) + others)/2 - opts.D0) <= opts.dTol);
    end
    how = sprintf('neither at D nor paired with another line about D = %g MHz', opts.D0);
end

function [ok, how] = resonancePosition(cr, F, x, opts)
% Could a dip at cr be the resonance, as far as position goes? Not if an NV line belongs
% there (a predicted transition within huberMHz, or at D / mirroring a line about D when
% there are too few lines for a field), and not if it sits on a fitted line (that is a
% line-shape residual, not a separate dip).
    ok = false;
    L = F.lines;
    c = [L.centre]; ce = [L.centreErr]; w = [L.fwhm];
    if isempty(c); how = 'no NV line to compare it with'; return; end
    if any(abs(c - cr) <= w/2)
        how = 'it sits on a fitted line'; return
    end
    if numel(c) >= 3
        sol = solveField(c, ce, x, opts);
        if sol.ok
            d = min(abs(sol.predicted - cr));
            if d <= opts.huberMHz
                how = sprintf('the %.1f G field puts a transition %.1f MHz away', sol.B, d);
                return
            end
            ok = true; how = sprintf('%.1f MHz from the nearest predicted transition', d);
            return
        end
    end
    if abs(cr - opts.D0) <= opts.dTol
        how = 'it sits at D'; return
    end
    if any(abs((cr + c)/2 - opts.D0) <= opts.dTol)
        how = 'it pairs with a line about D'; return
    end
    ok = true; how = 'not at D and pairs with nothing about D';
end

function [F, msgs] = residualResonance(F, x, y, ee, noise, opts, msgs, hint)
% A resonance too broad, or too far off the edge, for the line search to have found: the
% deepest BROAD deficit left in the residual, tried as one, kept only on a BIC gain of 10
% AND if it passes the same position and width tests as a converted line.
    r = y - F.model;
    dx = median(diff(x)); span = max(x) - min(x);
    W = max(5, 2*floor(max(10, 0.04*span)/dx/2) + 1);    % ~10 MHz, odd
    W = min(W, 2*floor(numel(x)/4) + 1);
    rs = boxSmooth(r, W);
    sigS = noise / sqrt(W);
    if isempty(hint)
        [mn, i] = min(rs);
        if -mn < opts.resonanceSnr * sigS; return; end
        cSeed = x(i);
    else
        [~, i] = min(abs(x - hint));
        cSeed = hint;
    end
    d0 = max(-rs(i), 3*sigS);
    iL = i; while iL > 1 && -rs(iL) > d0/2; iL = iL - 1; end
    iR = i; while iR < numel(x) && -rs(iR) > d0/2; iR = iR + 1; end
    wSeed = max(x(iR) - x(iL), max(10, 4*dx));
    if iL == 1 || iR == numel(x); wSeed = 2*wSeed; end     % runs off the edge: half is outside

    [posOk, how] = resonancePosition(cSeed, F, x, opts);
    if ~posOk && isempty(hint)
        msgs{end+1} = sprintf(['resonance: the residual dip at %.1f MHz (%.1f sigma) is not ' ...
                               'it -- %s'], cSeed, -rs(i)/sigS, how);
        return
    end
    L = F.lines;
    nv = starts([L.centre], [L.depth], [L.fwhm]/2);
    F2 = coreFit(x, y, ee, nv, struct('c', cSeed, 'a', d0, 'g', wSeed/2), opts, ...
                 opts.minDepthSigma, false);           % the dBIC below decides, not its depth
    if ~F2.ok || ~F2.resonance.found
        why = F2.why;
        if F2.ok; why = F2.resonance.why; end
        msgs{end+1} = sprintf('resonance: the residual dip at %.1f MHz did not fit as one (%s)', ...
                              cSeed, why);
        return
    end
    n = numel(x);
    s2 = max(min(F.chi2Reduced, F2.chi2Reduced), 0.05^2);
    dBic = (F.chi2/s2 + F.nPar*log(n)) - (F2.chi2/s2 + F2.nPar*log(n));
    R = F2.resonance;
    [posOk2, how2] = resonancePosition(R.centre, F2, x, opts);
    wNV = median([F2.lines.fwhm]);
    broad = R.fwhm >= opts.resonanceWidthRatio * wNV;
    if dBic > 10 && ((posOk2 && broad) || ~isempty(hint))
        msgs{end+1} = sprintf(['resonance: found in the residual at %.1f MHz, FWHM %.1f ' ...
                               '(dBIC %.0f; %s)'], R.centre, R.fwhm, dBic, how2);
        F = F2;
    else
        if ~posOk2; why = how2;
        elseif ~broad; why = sprintf('%.1f MHz wide, no broader than the NV lines (%.1f)', R.fwhm, wNV);
        else; why = sprintf('dBIC only %.1f', dBic);
        end
        msgs{end+1} = sprintf('resonance: tried %.1f MHz from the residual, rejected (%s)', ...
                              R.centre, why);
    end
end

function rs = boxSmooth(r, W)
    k = ones(W, 1);
    rs = conv(r(:), k, 'same') ./ conv(ones(numel(r), 1), k, 'same');
end

function res = resonanceSeedAt(x, y, env, noise, c, g)
    cc = min(max(c, min(x)), max(x));
    a = max(interp1(x, env - y, cc, 'linear', 'extrap'), 3*noise);
    res = struct('c', c, 'a', a, 'g', g);
end

function res = resonanceAsSeed(R)
    res = [];
    if R.found
        res = struct('c', R.centre, 'a', R.depth, 'g', R.fwhm/2);
    end
end

function [m, hint] = resonanceMode(v)
    hint = [];
    if isnumeric(v) && isscalar(v) && v > 100
        m = 'on'; hint = double(v); return
    end
    m = triMode(v, 'resonance');
end

% ======================================================================
% REFINEMENT (mirror, split, physics seeding) -- formerly inside calibFieldSurvey
% ======================================================================

function tf = refineApplies(F, x, opts)
    tf = false;
    m = triMode(opts.refine, 'refine');
    if strcmp(m, 'off') || ~isempty(opts.nLines) || F.nLines < 2; return; end
    if strcmp(m, 'on'); tf = true; return; end
    tf = min(x) < opts.D0 && max(x) > opts.D0 && (max(x) - min(x)) >= 60;
end

function [fit, msgs] = refineLines(fit, x, y, ee, noise, env, opts, msgs)
% Two failures, two repairs, interleaved:
%   ASYMMETRIC OMISSION -- a line with no partner. Every orientation gives a pair
%       straddling D, so its mirror 2D - f is a start point physics supplies.
%   SYMMETRIC BLEND -- several lines merged on BOTH sides, so the spectrum still looks
%       symmetric. Only a width test finds it: split the widest line into 3 and 2.
% Splitting one side breaks the symmetry, which the next mirror pass repairs. Then the
% solved FIELD seeds every predicted transition -- the only source of seeds inside a blob
% that neither side has resolved. Every refit is judged by acceptRefit; a seeded line that
% is not there fails its depth test and costs a refit, nothing else.
    for pass = 1:8
        improved = false;
        c = [fit.lines.centre];
        % D from the OUTERMOST pair (largest splitting, least blended), not assumed: the
        % nominal 2870 put every predicted partner 1.9 MHz out on a 2870.94 run.
        Dest = (min(c) + max(c))/2;
        if abs(Dest - opts.D0) > opts.dTol; Dest = opts.D0; end
        % Mirror EVERY line, not only unpaired ones: at 45 G a 43 MHz blob sat within
        % tolerance of one high-side line's mirror, counted as paired, and the other
        % high-side line's mirror was never proposed.
        allMir = 2*Dest - c;
        allMir = allMir(allMir >= min(x) & allMir <= max(x));
        cand = sort([c(:); allMir(:)]);
        cand = cand([true; diff(cand) > opts.pairTol]);
        if numel(cand) > numel(c)
            f2 = refitSeeds(fit, cand, x, y, ee, noise, env, opts);
            [take, why] = acceptRefit(fit, f2, opts);
            msgs{end+1} = sprintf('pass %d mirror about D = %.2f, %d seeds (%s): %s', pass, ...
                                  Dest, numel(cand), strtrim(sprintf('%.1f ', cand)), ...
                                  describeRefit(f2, take, why));             %#ok<AGROW>
            if take; fit = f2; improved = true; end
        end

        c = [fit.lines.centre]; w = [fit.lines.fwhm];
        [wMax, iw] = max(w);
        if numel(w) >= 2 && wMax > opts.splitRatio * min(w)
            % Two and three components: at 45 G the three off-axis orientations span ~22 MHz
            % a side, and cutting that in two leaves a residual a 2-way split cannot remove.
            % But the 3-way is taken over an accepted 2-way only if the extra component
            % EARNS it -- a chi2 drop to splitGain, or a narrower widest line. On a tie the
            % spare component otherwise wanders off to fit baseline curvature at the sweep
            % edge (a 3011.8 MHz "line" on a 45 G synthetic).
            best = []; bestN = NaN;
            for nSub = [2 3]
                off = linspace(-wMax/2.5, wMax/2.5, nSub);
                seeds = c; seeds(iw) = [];
                seeds = sort([seeds(:); c(iw) + off(:)]);
                fc = refitSeeds(fit, seeds, x, y, ee, noise, env, opts);
                [tk, wy] = acceptRefit(fit, fc, opts);
                msgs{end+1} = sprintf('pass %d split %.1f MHz line at %.1f into %d: %s', ...
                                      pass, wMax, c(iw), nSub, describeRefit(fc, tk, wy)); %#ok<AGROW>
                if ~tk; continue; end
                if isempty(best) || fc.chi2Reduced < opts.splitGain*best.chi2Reduced || ...
                        max([fc.lines.fwhm]) < opts.widthGain*max([best.lines.fwhm])
                    best = fc; bestN = nSub;
                end
            end
            if ~isempty(best)
                msgs{end+1} = sprintf(['pass %d: accepted the %d-way split of the %.1f MHz ' ...
                    'line at %.1f -> %d lines, chi2 %.2f -> %.2f'], pass, bestN, wMax, c(iw), ...
                    best.nLines, fit.chi2Reduced, best.chi2Reduced);          %#ok<AGROW>
                fit = best; improved = true;
            end
        end
        if ~improved; break; end
    end

    if opts.physicsSeed
        for pass = 1:3
            c = [fit.lines.centre];
            if numel(c) < 3; break; end
            sol = solveField(c, [fit.lines.centreErr], x, opts);
            if ~sol.ok
                msgs{end+1} = sprintf('pass %d physics seed: field solve failed (%s)', ...
                                      pass, sol.why);                        %#ok<AGROW>
                break
            end
            seeds = sol.predicted(sol.predicted >= min(x) & sol.predicted <= max(x));
            seeds = sort([seeds(:); sol.unmatched(:)]);
            seeds = seeds([true; diff(seeds) > opts.pairTol]);
            f2 = refitSeeds(fit, seeds, x, y, ee, noise, env, opts);
            [take, why] = acceptRefit(fit, f2, opts);
            msgs{end+1} = sprintf(['pass %d physics seed: |B| = %.2f G, tilt %.2f deg -> ' ...
                '%d seeds: %s'], pass, sol.B, sol.tiltDeg, numel(seeds), ...
                describeRefit(f2, take, why));                               %#ok<AGROW>
            if ~take; break; end
            fit = f2;
        end
    end
end

function f2 = refitSeeds(fit, seeds, x, y, ee, noise, env, opts)
    [c0, a0, g0] = seedStarts(x, y, seeds, env, noise);
    f2 = coreFit(x, y, ee, starts(c0, a0, g0), resonanceAsSeed(fit.resonance), opts, ...
                 opts.refineDepthSigma, false);
end

function [take, why] = acceptRefit(fit, f2, opts)
% Does a seeded refit replace the current fit? Three ways in, and it SAYS WHY either way:
%   gained a line, chi2 not materially worse
%   the WIDEST line got narrower (widthGain), chi2 not worse -- the one that matters for a
%       blend: splitting a 44 MHz blob into 12 MHz lines barely moves chi2, because the blob
%       already absorbed most of the deficit, so a chi2-only test refuses an obvious split
%   same lines, chi2 down to splitGain of its value
    take = false;
    if isempty(f2) || ~isstruct(f2) || ~f2.ok
        why = 'fit failed';
        if isstruct(f2) && ~isempty(f2.why); why = sprintf('fit failed (%s)', f2.why); end
        return
    end
    chi0 = fit.chi2Reduced; chi1 = f2.chi2Reduced;
    w0 = max([fit.lines.fwhm]); w1 = max([f2.lines.fwhm]);
    if f2.nLines > fit.nLines && chi1 < 1.05*chi0
        take = true; why = sprintf('gained %d line(s)', f2.nLines - fit.nLines);
    elseif w1 < opts.widthGain*w0 && chi1 < 1.05*chi0
        take = true; why = sprintf('widest line %.1f -> %.1f MHz', w0, w1);
    elseif f2.nLines == fit.nLines && chi1 < opts.splitGain*chi0
        take = true; why = sprintf('chi2 %.2f -> %.2f', chi0, chi1);
    elseif f2.nLines > fit.nLines
        why = sprintf('gained %d line(s) but chi2 %.2f -> %.2f is worse', ...
                      f2.nLines - fit.nLines, chi0, chi1);
    elseif f2.nLines < fit.nLines
        why = sprintf(['LOST %d line(s) -- the seeded component did not survive its ' ...
                       '%.1f-sigma depth test'], fit.nLines - f2.nLines, opts.refineDepthSigma);
    else
        why = sprintf('no gain: %d lines, chi2 %.2f -> %.2f, widest %.1f -> %.1f MHz', ...
                      f2.nLines, chi0, chi1, w0, w1);
    end
end

function s = describeRefit(f2, take, why)
    if isstruct(f2) && f2.ok
        if take; v = 'ACCEPTED'; else; v = 'rejected'; end
        s = sprintf('%d lines, chi2 %.2f, widest %.1f MHz -- %s (%s)', ...
                    f2.nLines, f2.chi2Reduced, max([f2.lines.fwhm]), v, why);
    else
        s = sprintf('REJECTED (%s)', why);
    end
end

% ======================================================================
% THE FIELD
% ======================================================================

function sol = solveField(c, ce, x, opts)
    sol = esrFieldSolve(c, ce, struct('D0', opts.D0, 'dTol', opts.dTol, ...
                        'huberMHz', opts.huberMHz, 'sweepRange', [min(x) max(x)], ...
                        'verbose', false));
end

function fld = blankField(n)
    fld = struct('kind', 'none', 'ok', false, 'B', NaN, 'tiltDeg', NaN, 'D', NaN, ...
                 'Bpar', NaN, 'projections', [], 'alignedAxisIdx', [], ...
                 'lineAxis', nan(1, n), 'aligned', false(1, n), 'note', '', ...
                 'solve', struct('ok', false, 'why', 'not attempted', 'B', NaN, ...
                                 'tiltDeg', NaN, 'D', NaN));
end

function fld = fieldFrom(F, x, opts, K)
    L = F.lines; n = numel(L);
    fld = blankField(n);
    if strcmp(triMode(opts.field, 'field'), 'off')
        fld.kind = 'off'; fld.note = 'field analysis switched off'; return
    end
    if n == 0; return; end
    c = [L.centre]; ce = [L.centreErr]; w = [L.fwhm];
    if n >= 3
        sol = solveField(c, ce, x, opts);
        fld.solve = sol;
        if sol.ok
            fld.kind = 'vector'; fld.ok = true;
            fld.B = sol.B; fld.tiltDeg = sol.tiltDeg; fld.D = sol.D;
            fld.projections = sol.projections; fld.alignedAxisIdx = sol.alignedAxisIdx;
            if isfield(sol, 'matched') && ~isempty(sol.matched)
                mc = [sol.matched.centre]; ma = [sol.matched.axis];
                for k = 1:n
                    [dd, j] = min(abs(mc - c(k)));
                    if dd < 1e-6; fld.lineAxis(k) = ma(j); end
                end
                fld.aligned = fld.lineAxis == sol.alignedAxisIdx;
            end
            notes = {};
            if isfield(sol, 'outerVisible') && ~sol.outerVisible
                notes{end+1} = 'the aligned pair is outside the sweep, so |B| is weakly constrained';
            end
            if sol.tiltDeg < 2.5
                notes{end+1} = 'below ~2.5 deg the line positions stop resolving the tilt';
            end
            if ~isempty(sol.unmatched)
                notes{end+1} = sprintf('%d line(s) match no NV transition (%s MHz)', ...
                                       numel(sol.unmatched), numList(sol.unmatched, '%.1f', ', '));
            end
            fld.note = strjoin(notes, '; ');
            return
        end
        fld.note = ['field solve failed: ' sol.why];
    end
    if n == 2 && abs(mean(c) - opts.D0) <= opts.dTol
        fld.kind = 'pair'; fld.ok = true;
        fld.Bpar = abs(c(2) - c(1)) / (2*K.gamma);
        fld.D = mean(c);
        fld.note = 'one pair: B_par along that orientation; |B| >= B_par, direction open';
        if abs(c(2) - c(1)) < max(w)
            fld.note = [fld.note '. Split by less than a linewidth: strain or a small field'];
        end
    elseif n == 1 && abs(c(1) - opts.D0) <= opts.dTol
        fld.kind = 'zero'; fld.ok = true; fld.D = c(1); fld.Bpar = 0;
        fld.note = sprintf(['one line at D: zero field (any splitting is under the %.1f MHz ' ...
                            'linewidth, |B_par| < ~%.1f G)'], w(1), w(1)/(2*K.gamma));
    elseif isempty(fld.note)
        fld.note = 'the lines do not determine a field';
    end
end

% ======================================================================
% THE HYPERFINE
% ======================================================================

function hf = hyperfineStage(x, y, ee, F, tau, opts, K)
    hf = hfBlank();
    m = triMode(opts.hyperfine, 'hyperfine');
    if strcmp(m, 'off'); hf.why = 'switched off'; return; end
    if F.nLines ~= 1
        hf.why = sprintf(['%d NV lines in the sweep -- the multiplet test needs exactly one ' ...
                          '(crop with range)'], F.nLines);
        return
    end
    L = F.lines(1);
    dx = median(diff(x));
    % Only the STEP gates it, not L.fwhm: on a resolved triplet that is the envelope (see THE
    % HYPERFINE in the header).
    if strcmp(m, 'auto') && dx > 0.6
        hf.why = sprintf('not resolvable here: step %.2f MHz (needs <= 0.6)', dx);
        return
    end
    % A window around the line: the multiplet model has one baseline and no resonance.
    half = max([5*L.fwhm, 4*K.A15N + 2*L.fwhm, 12]);
    in = x >= L.centre - half & x <= L.centre + half;
    hf = hfRun(x(in), y(in), ee(in), tau, opts, K);
end

function hf = hfRun(xs, ys, es, tau, opts, K)
    hf = hfBlank();
    if numel(xs) < 12
        hf.why = sprintf('only %d points around the line', numel(xs)); return
    end
    bo = opts.baselineOrder;
    if max(xs) - min(xs) < 60; bo = min(bo, 1); end     % order 2 eats a line in a narrow window
    ho = struct('nLines', opts.hfCount, 'nLinesList', [1 2 3], 'shape', lower(opts.hfShape), ...
                'depths', lower(opts.hfDepths), 'splitting', opts.hfSplitting, ...
                'splitRange', opts.hfSplitRange, 'baselineOrder', bo, ...
                'maxFwhmMHz', opts.maxFwhmMHz, 'pulseRotationDeg', opts.pulseRotationDeg, ...
                'maxIter', opts.maxIter);
    R = hfSelect(xs, ys, es, tau, ho, K);
    hf = copyFields(hf, R);
    hf.tested = true; hf.ok = R.ok; hf.why = R.why;
    hf.window = [min(xs) max(xs)];
    hf.f = xs; hf.y = ys;
    hf.resolved = R.ok && R.nDips >= 2 && ~strcmp(R.isotope, 'unknown') && ...
                  ~strcmp(R.confidence, 'ambiguous');
end

function [F, hf] = multipletFallback(x, y, ee, tau, opts, K)
% When the line search failed: the multiplet model over the whole sweep, if the hyperfine
% is asked for, or on 'auto' if the step could resolve it. A multiplet it finds becomes ONE
% transition line -- out.lines keeps its meaning -- with the centre and depth of the
% multiplet and the FWHM of its envelope, and the model, residual and chi2 are the
% multiplet's.
    F = blankFit(); hf = hfBlank();
    m = triMode(opts.hyperfine, 'hyperfine');
    if strcmp(m, 'off') || (strcmp(m, 'auto') && median(diff(x)) > 0.6); return; end
    if max(x) - min(x) > 80; return; end      % a wide sweep has other lines: not one multiplet
    hf = hfRun(x, y, ee, tau, opts, K);
    if ~hf.ok || ~(hf.nDips >= 1); return; end

    ff = linspace(min(x), max(x), 8000)';
    base = zeros(size(ff));
    for j = 1:numel(hf.baseline)
        base = base + hf.baseline(j) * (ff - hf.xref).^(j-1);
    end
    dip = base - hf.modelAt(ff);
    [depth, ip] = max(dip);
    in = find(dip >= depth/2);
    fwhm = ff(in(end)) - ff(in(1));
    L = struct('centre', hf.f0, 'centreErr', hf.f0Err, 'fwhm', fwhm, 'fwhmErr', NaN, ...
               'depth', depth, 'depthErr', hf.depthErr);
    if ~isfinite(L.centre); L.centre = ff(ip); end
    F.lines = L; F.nLines = 1;
    F.centreCov = hf.f0Err^2;
    [F.splitting, F.splittingErr, F.pairMean, F.pairMeanErr] = pairTables(L.centre, F.centreCov);
    F.baseline = hf.baseline; F.xref = hf.xref;
    F.model = hf.model; F.residual = hf.residual; F.modelAt = hf.modelAt;
    F.chi2Reduced = hf.chi2Reduced; F.birge = hf.birge;
    F.ok = true;
end

function R = hfSelect(x, y, ee, tau, ho, K)
% Fit 0, 1, 2 and 3 dips (each with every shape, equal and free depths) and let BIC choose.
%
% Start: a grid over (f0, A, w) with the baseline and depth solved linearly at every grid
% point, so the start is that grid's global minimum and nothing depends on peak finding.
% Then bounded Levenberg-Marquardt on the nonlinear parameters with the linear ones
% re-solved at every step (variable projection), depths >= 0. BIC = chi2/s^2 + k ln(n), s^2
% the reduced chi2 of the best candidate (Birge: the error bars give the shape, the
% residuals the scale). Strength from dBIC to the runner-up: < 2 ambiguous, 2-6 positive,
% 6-10 strong, > 10 decisive. A 14N label needs 3 dips with A within max(3 sigma, 0.3 MHz)
% of 2.16; 15N needs 2 dips within that of 3.03; anything else with 2+ dips is 'unknown'.
    R = hfBlank();
    n = numel(x);
    D = struct('x', x, 'y', y, 'sw', 1 ./ ee, 'xr', mean(x), 'nB', ho.baselineOrder + 1, ...
               'dx', median(diff(x)), 'span', max(x) - min(x));
    try
        shapes = shapeList(ho.shape, tau);
    catch err
        R.why = err.message; return
    end
    P = pulseParams(tau, ho.pulseRotationDeg);
    R.tau = tau;

    if isempty(ho.nLines)
        Ns = unique([0, ho.nLinesList(:)']);
    else
        Ns = unique([0, ho.nLines]);
    end
    freeFromStart = strcmp(ho.depths, 'free');
    tryFree = strcmp(ho.depths, 'auto');
    vars = [];
    for N = Ns
        best = [];
        shp = shapes;
        if N == 0; shp = {'none'}; end
        for s = 1:numel(shp)
            Fc = fitCandidate(D, N, shp{s}, freeFromStart, ho, P, []);
            if isempty(best) || Fc.chi2 < best.chi2; best = Fc; end
        end
        vars = [vars, best];                                                  %#ok<AGROW>
        if tryFree && N >= 2
            % Polarised nitrogen makes the depths unequal, and an equal-depth model can then
            % lose to the WRONG count. So each count also competes with free depths.
            vars = [vars, fitCandidate(D, N, best.shape, true, ho, P, best.theta)]; %#ok<AGROW>
        end
    end

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
    R.prob = nan(1, 4);
    R.prob([cands.nDips] + 1) = pr;

    if isempty(ho.nLines)
        [~, iw] = min(bic);
    else
        iw = find([cands.nDips] == ho.nLines, 1);
    end
    final = cands(iw);
    others = setdiff(1:numel(cands), iw);
    if isempty(others)
        gap = Inf; iRun = [];
    else
        [gap, j] = min(bic(others) - bic(iw));
        iRun = others(j);
    end
    R.candidates = cands;

    R.depthTest = struct('pEqual', NaN, 'used', ternary(final.freeDepths, 'free', 'equal'));
    if tryFree && final.nDips >= 2
        idx = find([vars.nDips] == final.nDips);
        eq = vars(idx(~[vars(idx).freeDepths]));
        fr = vars(idx([vars(idx).freeDepths]));
        R.depthTest.pEqual = chi2sf(max(eq.chi2 - fr.chi2, 0)/s2, fr.k - eq.k);
    end

    Rn = finalNumbers(D, final, ho, P);
    if ~Rn.ok
        R.why = Rn.why; return
    end
    R = copyFields(R, Rn);
    R.nDips = final.nDips;
    R.shape = final.shape;
    R.modelAt = final.modelAt;
    R.model = final.yhat;
    R.residual = y - final.yhat;

    nd = final.nDips;
    tolA = max(3*R.splittingErr, 0.3);
    switch nd
        case 0; R.isotope = 'none';
        case 1; R.isotope = 'unresolved';
        case 2; R.isotope = ternary(abs(R.splitting - K.A15N) <= tolA, '15N', 'unknown');
        case 3; R.isotope = ternary(abs(R.splitting - K.A14N) <= tolA, '14N', 'unknown');
    end
    if ~isempty(ho.nLines)
        R.confidence = 'forced';
    elseif ~isfinite(gap) || gap > 10
        R.confidence = 'decisive';
    elseif gap > 6
        R.confidence = 'strong';
    elseif gap > 2
        R.confidence = 'positive';
    else
        R.confidence = 'ambiguous';
    end
    R.runnerUp = [];
    if ~isempty(iRun); R.runnerUp = cands(iRun); end
    R.verdict = verdictText(R, gap, K.A14N, K.A15N);
    R.warnings = collectWarnings(R, final, D, K.A14N, K.A15N);
    R.ok = true;
end

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
% For every (w, eta, A) the baseline and the common depth are solved in closed form for ALL
% f0 at once, by projecting out the baseline.
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
% Levenberg-Marquardt with Marquardt column scaling; every trial step clamped to the box.
    r = rfun(th); c = r'*r;
    lambda = 1e-2;
    np = numel(th);
    for it = 1:maxIter                                                      %#ok<FXUP>
        J = numJac(rfun, th, hs, lo, hi);
        H = J'*J; g = J'*r;
        sc = sqrt(max(diag(H), realmin));
        Hs = H ./ (sc*sc'); gs = g ./ sc;
        took = false; rel = 0;
        for tries = 1:25                                                    %#ok<FXUP>
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
% Least squares with the depth columns >= 0. At most 3 depths, so the exact answer is the
% best feasible subset; the unconstrained solve is the answer whenever it is feasible.
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
% Line centres carry the f0-A covariance: var(f0 + m A) = s0^2 + m^2 sA^2 + 2 m c0A.
    R = struct('ok', false, 'why', '');
    n  = numel(D.x);
    nB = D.nB;
    nth = numel(F.theta);
    p  = [F.theta(:); F.beta(:)];
    np = numel(p);
    dof = max(n - F.k, 1);
    birge = max(sqrt(F.chi2 / dof), 0.05);

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
            error('esrFitRobust: unknown hyperfine shape %s', shp);
    end
end

function S = pulseShape(U, w, P)
% Square-pulse response convolved with a Gaussian of FWHM w:
%   Pr(u) = W^2/(W^2+u^2) sin^2(pi tau sqrt(W^2+u^2)),   W = theta/(2 pi tau)
% tabulated on a fine grid and splined. Cached: the grid search and every LM step at the
% same w reuse the table.
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
    % Uniform knots, so the spline piece is found by arithmetic instead of ppval's search.
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
                error(['esrFitRobust: hfShape ''pulse'' needs the pulse length. Pass tau ' ...
                       '(us), or a pulsed experiment or save it can be read from.']);
            end
            shapes = {'pulse'};
        otherwise
            error('esrFitRobust: hfShape must be auto, lorentz, gauss, voigt or pulse.');
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

function t = verdictText(o, gap, A14N, A15N)
    switch o.nDips
        case 0
            t = 'no line in this window';
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

function W = collectWarnings(o, F, D, A14N, A15N)
    W = {};
    if F.nDips == 0; return; end
    nm = F.names;
    for j = 1:numel(nm)
        if o.atBound(j)
            switch nm{j}
                case 'A'
                    W{end+1} = sprintf(['the splitting sits on its bound (%.2f MHz): widen ' ...
                        'hfSplitRange, or the pattern is not hyperfine'], F.theta(j));   %#ok<AGROW>
                case 'w'
                    if F.theta(j) <= F.lo(j)*1.0001
                        if strcmp(F.shape, 'pulse')
                            W{end+1} = ['the inhomogeneous width went to ~0: the pulse alone ' ...
                                'sets the line width'];                                 %#ok<AGROW>
                        else
                            W{end+1} = sprintf(['the width hit its floor (%.2f MHz, half the ' ...
                                'frequency step): the lines are not resolved by this grid'], ...
                                F.theta(j));                                            %#ok<AGROW>
                        end
                    else
                        W{end+1} = 'the width hit its ceiling: no resolved line';       %#ok<AGROW>
                    end
                case 'f0'
                    W{end+1} = 'the pattern centre sits on the edge of the window';     %#ok<AGROW>
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
            'or the drive differs across the lines'], o.depthTest.pEqual);           %#ok<AGROW>
        if any([o.lines.depth] <= 0)
            W{end+1} = ['a line came out with zero depth: the pattern is not a full ' ...
                'multiplet, check the count'];                                       %#ok<AGROW>
        end
    end
    if F.nDips >= 2 && strcmp(o.isotope, 'unknown')
        W{end+1} = sprintf(['spacing %.3f MHz matches neither 14N (%.2f) nor 15N (%.2f); ' ...
            'two orientations or sites a few MHz apart look like this'], ...
            o.splitting, A14N, A15N);                                                %#ok<AGROW>
    end
    if 2*D.dx > o.fwhm
        W{end+1} = sprintf(['the frequency step (%.2f MHz) is more than half the line FWHM ' ...
            '(%.2f): the width is barely sampled, use a finer step for line shapes'], ...
            D.dx, o.fwhm);                                                           %#ok<AGROW>
    end
    if strcmp(o.confidence, 'ambiguous')
        W{end+1} = 'the line count is ambiguous: more averages, or a finer/longer sweep'; %#ok<AGROW>
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

function p = chi2sf(x, k)
% Upper tail of chi-square with k dof, base MATLAB (no Statistics Toolbox).
    if k <= 0; p = NaN; return; end
    p = gammainc(x/2, k/2, 'upper');
end

function hf = hfBlank()
    hf = struct('tested', false, 'why', '', 'resolved', false, 'ok', false, ...
        'nDips', NaN, 'isotope', '', 'confidence', '', 'verdict', '', 'warnings', {{}}, ...
        'f0', NaN, 'f0Err', NaN, 'splitting', NaN, 'splittingErr', NaN, ...
        'fwhm', NaN, 'fwhmErr', NaN, 'widthParam', NaN, 'widthParamErr', NaN, ...
        'depth', NaN, 'depthErr', NaN, 'shape', '', 'eta', NaN, 'tau', [], ...
        'lines', emptyLines(), 'centreCov', [], 'candidates', [], 'runnerUp', [], ...
        'prob', nan(1, 4), 'depthTest', struct('pEqual', NaN, 'used', ''), ...
        'model', [], 'residual', [], 'modelAt', [], 'formula', '', 'params', [], ...
        'cov', [], 'chi2Reduced', NaN, 'birge', NaN, 'baseline', [], 'xref', NaN, ...
        'window', [NaN NaN], 'f', [], 'y', [], 'r2', NaN, 'r2Ceiling', NaN);
end

function s = hfLabel(hf)
    switch hf.isotope
        case '14N'; s = '14N triplet';
        case '15N'; s = '15N doublet';
        otherwise;  s = sprintf('%d dips', hf.nDips);
    end
end

% ======================================================================
% DATA
% ======================================================================

function [x, y, e, info] = getData(opts)
    info = struct('reps', 1, 'nAvg', 1, 'dropped', [], 'tau', [], 'source', '');
    if ~isempty(opts.file)
        if exist(opts.file, 'file') ~= 2
            error('esrFitRobust: no such file: %s', opts.file);
        end
        S = load(opts.file);
        s = findSaved(S);
        if isempty(s)
            error('esrFitRobust: no experiment with summaryData in %s', opts.file);
        end
        [x, y, e] = fromSummary(s.summaryData);
        info.tau = pulseFrom(s); info.source = 'file';
        return
    end
    X = opts.exp;
    if ~isempty(X)
        if isstruct(X)
            s = findSaved(X);
            if ~isempty(s)
                [x, y, e] = fromSummary(s.summaryData);
                info.tau = pulseFrom(s); info.source = 'saved struct';
            else
                [x, y, e] = fromVectorStruct(X); info.source = 'struct';
            end
            return
        end
        if isnumeric(X)
            x = X(:, 1); y = X(:, 2);
            if size(X, 2) >= 3; e = X(:, 3); else; e = nan(size(x)); end
            info.source = 'matrix'; return
        end
        [x, y, e, reps, nAvg, pinfo] = esrPullData(X, struct('avgOutlierMad', ...
            opts.avgOutlierMad, 'verbose', opts.verbose, 'name', 'esrFitRobust'));
        info.reps = reps; info.nAvg = nAvg; info.dropped = pinfo.dropped;
        info.tau = pulseFrom(X); info.source = 'experiment';
        return
    end
    if isempty(opts.f) || isempty(opts.y)
        error(['nothing to fit: give exp, file, or f and y. ' ...
               'esrFitRobust(struct(''exp'', a)).']);
    end
    x = opts.f(:); y = opts.y(:); e = opts.e(:);
    if isempty(e); e = nan(size(x)); end
    if numel(y) ~= numel(x) || numel(e) ~= numel(x)
        error('esrFitRobust: f, y and e have different lengths (%d, %d, %d).', ...
              numel(x), numel(y), numel(e));
    end
    info.source = 'vectors';
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
% Any struct carrying a sweep: plain vectors, a previous esrFitRobust or calibFieldSurvey
% output, an esrSynthSpectrum.
    x = pick(s, {'f', 'freqs', 'freq', 'frequency', 'frequencies', 'x'});
    y = pick(s, {'y', 'signal', 'value', 'values'});
    e = pick(s, {'e', 'sterr', 'err', 'error', 'errors'});
    if isempty(x) || isempty(y)
        error('esrFitRobust: this struct has no recognisable sweep (f/freqs and y/signal).');
    end
    if numel(e) ~= numel(x); e = nan(size(x)); end
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
    out = struct('ok', false, 'why', '', 'lines', emptyLines(), 'nLines', 0, ...
        'f0', NaN, 'width', NaN, 'contrast', NaN, 'SNR', NaN, 'coeffError', nan(1, 3), ...
        'r2', NaN, 'r2Ceiling', NaN, 'resonance', noResonance(''), ...
        'field', blankField(0), 'solve', struct('ok', false, 'why', 'not attempted'), ...
        'B', NaN, 'tiltDeg', NaN, 'D', NaN, 'hyperfine', hfBlank(), 'refineLog', {{}}, ...
        'centreCov', [], 'splitting', [], 'splittingErr', [], 'pairMean', [], ...
        'pairMeanErr', [], 'baseline', [], 'xref', NaN, 'chi2Reduced', NaN, 'birge', NaN, ...
        'residual', [], 'model', [], 'modelAt', [], 'f', [], 'y', [], 'e', [], ...
        'repeats', 1, 'averages', 1, 'dropped', [], 'source', '', 'tau', [], 'opts', opts);
end

function F = blankFit()
    F = struct('ok', false, 'why', '', 'lines', emptyLines(), 'nLines', 0, ...
        'centreCov', [], 'splitting', [], 'splittingErr', [], 'pairMean', [], ...
        'pairMeanErr', [], 'baseline', [], 'xref', NaN, 'chi2Reduced', NaN, 'chi2', NaN, ...
        'birge', NaN, 'residual', [], 'model', [], 'modelAt', [], 'nPar', 0, ...
        'resonance', noResonance(''));
end

function R = noResonance(why)
    R = struct('found', false, 'centre', NaN, 'centreErr', NaN, 'fwhm', NaN, ...
               'fwhmErr', NaN, 'depth', NaN, 'depthErr', NaN, 'outside', false, 'why', why);
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

function out = failed(out)
    if out.opts.verbose; fprintf('\nesrFitRobust FAILED: %s\n', out.why); end
end

function report(o, K)
    L = o.lines; fld = o.field; R = o.resonance; hf = o.hyperfine;
    hdr = sprintf('%d line(s)', numel(L));
    if R.found; hdr = [hdr ' + the resonance']; end
    fprintf('\n--- esrFitRobust: %s ---\n', hdr);
    fprintf('%6s %17s %14s %15s %8s\n', 'line', 'centre (MHz)', 'FWHM (MHz)', 'depth (%)', 'SNR');
    for i = 1:numel(L)
        sn = NaN; if numel(o.SNR) >= i; sn = o.SNR(i); end
        fprintf('%6d %9.3f+-%-7.3f %6.2f+-%-5.2f %7.3f+-%-5.3f %8.1f  %s\n', i, ...
                L(i).centre, L(i).centreErr, L(i).fwhm, L(i).fwhmErr, ...
                100*L(i).depth, 100*L(i).depthErr, sn, lineTag(fld, i));
    end
    if R.found
        where = '';
        if R.outside; where = ', centre OUTSIDE the sweep (extrapolated from the edge)'; end
        fprintf('%6s %9.3f+-%-7.3f %6.2f+-%-5.2f %7.3f+-%-5.3f %8s  the resonance, not NV%s\n', ...
                'res', R.centre, R.centreErr, R.fwhm, R.fwhmErr, 100*R.depth, 100*R.depthErr, ...
                '', where);
    end
    fprintf('R^2 %.5f (noise ceiling %.5f), reduced chi2 %.2f\n', ...
            o.r2, o.r2Ceiling, o.chi2Reduced);
    if isfinite(o.r2Ceiling) && o.r2 < o.r2Ceiling - 0.05
        fprintf(['  R^2 sits %.3f below what the noise allows, so this is a MODEL problem and ' ...
                 'not a\n  quiet spectrum. Try nLines, or a higher baselineOrder if the ' ...
                 'background is curved.\n'], o.r2Ceiling - o.r2);
    end
    if numel(L) > 1
        fprintf('splittings: %s MHz\n', numList(diff(sort([L.centre])), '%.3f', ', '));
    end

    switch fld.kind
        case 'vector'
            fprintf(['field: |B| = %.3f G, tilt %.2f deg from the aligned <111> axis, ' ...
                     'D = %.3f MHz\n'], fld.B, fld.tiltDeg, fld.D);
            fprintf('       projections %s G, line residual rms %.2f MHz\n', ...
                    numList(fld.projections, '%+.2f', ' '), fld.solve.rms);
        case 'pair'
            fprintf('field: one pair, B_par = %.3f G, D = %.3f MHz\n', fld.Bpar, fld.D);
        case 'zero'
            fprintf('field: zero -- one line, at D = %.3f MHz\n', fld.D);
    end
    if ~isempty(fld.note) && ~strcmp(fld.kind, 'off')
        fprintf('       %s\n', fld.note);
    end

    if hf.tested
        if ~hf.ok
            fprintf('hyperfine: test failed (%s)\n', hf.why);
        else
            fprintf('hyperfine: %s\n', hf.verdict);
            if hf.nDips >= 1
                fprintf('       f0 %.3f +- %.3f MHz, FWHM %.3f +- %.3f MHz (%s)', ...
                        hf.f0, hf.f0Err, hf.fwhm, hf.fwhmErr, hf.shape);
                if strcmp(hf.shape, 'pulse')
                    fprintf(', pulse %.3f us, inhomogeneous %.3f MHz', hf.tau, hf.widthParam);
                end
                fprintf('\n');
            end
            for i = 1:numel(hf.warnings)
                fprintf('       NOTE: %s\n', hf.warnings{i});
            end
        end
    elseif strcmp(triMode(o.opts.hyperfine, 'hyperfine'), 'on')
        fprintf('hyperfine: not tested -- %s\n', hf.why);
    end
end

function t = lineTag(fld, i)
    t = '';
    if ~strcmp(fld.kind, 'vector') || numel(fld.lineAxis) < i; return; end
    if ~isfinite(fld.lineAxis(i))
        t = 'matches no NV transition';
    elseif fld.aligned(i)
        t = 'ALIGNED axis';
    else
        t = sprintf('orientation %d', fld.lineAxis(i));
    end
end

function drawIt(o)
% ONE axes, in fitStyle's look (shared with the Rabi, Ramsey, echo and T1 fitters). What
% it adds over a bare fit: the resonance drawn as the background the NV lines hang from,
% the aligned-axis lines picked out, and the field in the title.
    S = fitStyle();
    COL_ALIGNED = [0.00 0.45 0.35];            % the axis every later calibration uses
    nm = 'ESR'; if ~isempty(o.opts.name); nm = ['ESR ' o.opts.name]; end
    fig = S.newFigure(nm, o.opts.showResid);
    if o.opts.showResid
        ax = S.axes(subplot(4, 1, 1:3));
    else
        ax = S.axes(fig);
    end

    H = []; Lg = {};
    H(end+1) = S.data(ax, o.f, o.y, o.e, S.col.s1);  Lg{end+1} = 'measured';
    ff = linspace(min(o.f), max(o.f), 4000)';
    yy = o.modelAt(ff);
    fitLabel = 'fit';
    hf = o.hyperfine;
    if hf.resolved
        inW = ff >= hf.window(1) & ff <= hf.window(2);
        yy(inW) = hf.modelAt(ff(inW));
        yy(find(diff(inW)) + 1) = NaN;            % no vertical connector at the window edge
        fitLabel = sprintf('fit (%s)', hfLabel(hf));
    end
    R = o.resonance;
    if R.found
        bg = esrModelAt(ff, o.baseline, o.xref, emptyLines(), R);
        H(end+1) = plot(ax, ff, bg, '--', 'LineWidth', 1.4, 'Color', S.col.s3);
        Lg{end+1} = sprintf('resonance %.1f MHz, FWHM %.0f MHz', R.centre, R.fwhm);
    end
    H(end+1) = S.model(ax, ff, yy);   Lg{end+1} = fitLabel;

    % A tick UNDER each line, not a line through the data: the positions are the result,
    % and drawing them over the points they came from makes both harder to read.
    S.tight(ax, o.f); S.headroom(ax, 0.26);
    yl = ylim(ax);
    tickTop = yl(1) + 0.05*(yl(2) - yl(1));
    hAl = [];
    if hf.resolved
        pos = [hf.lines.centre]; al = false(size(pos));
    else
        pos = [o.lines.centre];  al = o.field.aligned;
        if numel(al) ~= numel(pos); al = false(size(pos)); end
    end
    for i = 1:numel(pos)
        if al(i)
            h = plot(ax, [pos(i) pos(i)], [yl(1) tickTop], '-', 'LineWidth', 2.6, ...
                     'Color', COL_ALIGNED);
            if isempty(hAl); hAl = h; end
        else
            plot(ax, [pos(i) pos(i)], [yl(1) tickTop], '-', 'LineWidth', 1.6, ...
                 'Color', S.col.flag);
        end
    end
    if R.found && R.centre >= min(o.f) && R.centre <= max(o.f)
        plot(ax, [R.centre R.centre], [yl(1) tickTop], '-', 'LineWidth', 1.6, 'Color', S.col.s3);
    end
    if ~isempty(hAl)
        H(end+1) = hAl; Lg{end+1} = 'aligned-axis lines';
    end
    ylim(ax, yl);

    xlabel(ax, 'frequency  (MHz)', 'Color', S.col.ink2);
    ylab = S.yLabel;
    if ~isempty(o.opts.yLabel); ylab = o.opts.yLabel; end
    ylabel(ax, ylab, 'Color', S.col.ink2);
    S.legend(ax, H, Lg);
    S.title(ax, esrTitle(o));

    if o.opts.showResid
        axR = S.axes(subplot(4, 1, 4));
        S.zeroLine(axR, [min(o.f) max(o.f)]);
        sg = o.e; if isempty(sg); sg = ones(size(o.residual)); end
        S.resid(axR, o.f, o.residual ./ max(sg, eps), S.col.s1);
        xlabel(axR, 'frequency  (MHz)', 'Color', S.col.ink2);
        ylabel(axR, 'residual / \sigma', 'Color', S.col.ink2);
        S.tight(axR, o.f);
        S.title(axR, sprintf('\\chi^2_\\nu = %.2f', o.chi2Reduced));
    end
end

function T = esrTitle(o)
% Two lines: a third is clipped above a residual panel. The first carries the RESULT --
% the multiplet, the field, or the line -- the second how well it fits, and the resonance.
    fld = o.field; hf = o.hyperfine; L = o.lines;
    if hf.resolved
        a = sprintf(['%s:   A = %.3f \\pm %.3f MHz     f_0 = %.3f \\pm %.3f MHz     ' ...
                     'FWHM = %.2f \\pm %.2f MHz'], hfLabel(hf), hf.splitting, ...
                    hf.splittingErr, hf.f0, hf.f0Err, hf.fwhm, hf.fwhmErr);
    elseif strcmp(fld.kind, 'vector')
        a = sprintf(['|B| = %.2f G     tilt from the aligned axis = %.2f^\\circ     ' ...
                     'D = %.2f MHz     %d lines'], fld.B, fld.tiltDeg, fld.D, numel(L));
    elseif numel(L) == 1
        a = sprintf(['f_0 = %.3f \\pm %.3f MHz     FWHM = %.2f \\pm %.2f MHz     ' ...
                     'depth = %.3f \\pm %.3f %%'], L.centre, L.centreErr, L.fwhm, ...
                    L.fwhmErr, 100*L.depth, 100*L.depthErr);
        if strcmp(fld.kind, 'zero'); a = ['zero field:   ' a]; end
    elseif strcmp(fld.kind, 'pair')
        a = sprintf('one pair:   B_{||} = %.2f G     D = %.2f MHz     lines at %s MHz', ...
                    fld.Bpar, fld.D, numList([L.centre], '%.2f', ',  '));
    else
        a = sprintf('%d lines at %s MHz', numel(L), numList([L.centre], '%.2f', ',  '));
    end
    if hf.resolved          % the stats of the curve that is drawn
        b = sprintf(['R^2 = %.4f of %.4f     |     \\chi^2_\\nu = %.2f     |     %d points     ' ...
                     '|     P(%d dips) = %.3f, %s'], hf.r2, hf.r2Ceiling, hf.chi2Reduced, ...
                    numel(hf.f), hf.nDips, hf.prob(hf.nDips + 1), hf.confidence);
    else
        b = sprintf('R^2 = %.4f of %.4f     |     \\chi^2_\\nu = %.2f     |     %d points', ...
                    o.r2, o.r2Ceiling, o.chi2Reduced, numel(o.f));
    end
    if o.resonance.found
        b = sprintf('%s     |     resonance %.1f MHz', b, o.resonance.centre);
    end
    T = {a, b};
end

% ======================================================================
% SMALL HELPERS
% ======================================================================

function checkOptions(opts)
    triMode(opts.field, 'field'); triMode(opts.refine, 'refine');
    triMode(opts.hyperfine, 'hyperfine'); resonanceMode(opts.resonance);
    if ~any(strcmpi(opts.hfDepths, {'auto', 'equal', 'free'}))
        error('esrFitRobust: hfDepths must be ''auto'', ''equal'' or ''free''.');
    end
    if ~isempty(opts.hfCount) && ~any(opts.hfCount == [0 1 2 3])
        error('esrFitRobust: hfCount must be 0, 1, 2 or 3.');
    end
    if ~isempty(opts.nLines) && (opts.nLines < 1 || opts.nLines ~= round(opts.nLines))
        error('esrFitRobust: nLines must be a positive integer.');
    end
end

function m = triMode(v, name)
    if ischar(v) || (isstring(v) && isscalar(v))
        v = lower(char(v));
        if strcmp(v, 'auto')
            m = 'auto';
        elseif any(strcmp(v, {'on', 'true', 'yes'}))
            m = 'on';
        elseif any(strcmp(v, {'off', 'false', 'no', 'none'}))
            m = 'off';
        else
            error('esrFitRobust: %s must be ''auto'', true or false (got ''%s'').', name, v);
        end
    elseif (islogical(v) || isnumeric(v)) && isscalar(v)
        if v; m = 'on'; else; m = 'off'; end
    else
        error('esrFitRobust: %s must be ''auto'', true or false.', name);
    end
end

function s = numList(v, fmt, sep)
% strjoin(compose(...)) needs R2016b and compose does not exist in Octave.
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
