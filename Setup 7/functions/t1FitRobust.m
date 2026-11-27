function out = t1FitRobust(args)
%T1FITROBUST Robust NV T1 fit: stretched exponential with the exponent fitted, and a second
%component only if the data pays for it.
%
%   out = t1FitRobust(struct('exp', t))
%   out = t1FitRobust(struct('exp', t, 'model', 'double'))
%   out = t1FitRobust(struct('t', tau, 'y', y, 'e', e))        % plain vectors
%   out = t1FitRobust(struct('t', tau, 'y', Y, 'e', E))        % Y is 2xN: |0> and |1>
%
% ONE STRUCT IN, ONE STRUCT OUT -- the fifth fitter on the convention shared with
% ramseyFitRobust, rabiFitRobust, esrFitRobust and echoFitRobust. The data rides in the same
% struct as the options and everything comes back in `out`, with out.ok saying whether to
% believe it. Self-contained: no Curve Fitting Toolbox, no lab fitting helpers. The figure comes
% from fitStyle so all five look the same.
%
% ===================================================================
% WHAT A T1 IS, AND WHY IT NEEDS ITS OWN FITTER
% ===================================================================
%
% AWGExpT1 publishes TWO traces and they play different roles, which is not true of any of the
% other four experiments. The |1> trace relaxes; the |0> trace is a CONTROL that should not move
% at all. So the fit is not two views of one signal -- it is a signal and its null.
%
%       y_0(t) = c_0 + a_0 * D(t)        a_0 should come back at ZERO
%       y_1(t) = c_1 + a_1 * D(t)        a_1 is the contrast
%
% with D(t) shared. That structure turns the control into a number: out.flatnessSigma is a_0
% divided by its own error, and it should be under about 3. A |0> trace that decays as well means
% something common-mode is moving -- laser power, tracking, charge state -- and the T1 from the
% same fit is then measuring that too. No amount of care in the fitting fixes it, which is why
% the fit reports it rather than absorbing it.
%
% (If only one trace is given, that control is simply unavailable. The fit still works; the
% report says the control is missing, and it is worth knowing that it is.)
%
% ===================================================================
% THE MODEL, AND WHY p IS FITTED
% ===================================================================
%
%       D(t) = exp(-(t/T1)^p)                                   'stretched'  (default)
%       D(t) = w*exp(-t/T1fast) + (1-w)*exp(-t/T1slow)          'double'
%
% p is FITTED, bounds 0.3 to 3, not pinned at 1. A single exponential is what T1 relaxation is
% for one two-level system with one bath, and p is how the data says otherwise:
%
%   p = 1     one rate. The textbook case, and what a clean single NV should give.
%   p < 1     a DISTRIBUTION of rates -- several relaxation channels, or one that fluctuates
%             over the run. p = 0.5 is the classic signature of a spin bath with a spread of
%             couplings. This is the common departure and it is worth knowing about, because a
%             single exponential forced onto it returns a T1 that depends on where the sweep
%             stopped.
%   p > 1     relaxation that is SLOWER than exponential at short times, which for a T1 is
%             hard to arrange physically. Suspect a common-mode drift, or a plateau the sweep
%             never reached, before believing it.
%
% RATES, NOT TIMES, ARE WHAT ADD. Relaxation channels combine as 1/T1 = sum of rates, so when
% comparing runs -- or asking what a gate leak cost -- the rate is the quantity to work in. The
% report gives both, and out.ratePerUs is directly comparable with the leak rate t1GateTest
% measures: that whole test exists because an open microwave switch adds to exactly this number.
%
% AND T1 ALONE IS NOT COMPARABLE ACROSS RUNS WHEN p HAS MOVED, for the same reason as in
% echoFitRobust. out.coherenceArea = T1*gamma(1+1/p) is the shape-independent figure of merit;
% compare that when p changes.
%
% ===================================================================
% THE SECOND COMPONENT, AND WHY IT IS GUARDED
% ===================================================================
%
% A stretched exponential with p < 1 and a sum of two exponentials describe the same physics --
% more than one rate -- and on a sweep of twenty noisy points they are very nearly the same
% curve. Which means a two-component fit will almost always improve the residual slightly, and
% will almost always return two time constants with confident-looking error bars.
%
% So with model = 'auto' the stretched fit is the null and the double has to beat it on all
% three of:
%
%   1. the weighted SSR drops by more than 9 -- 3 sigma for the one extra parameter (double has
%      three nonlinear parameters against stretched's two);
%   2. the two time constants differ by more than a factor of 2. Closer than that and they are
%      not two components, they are one component and a fitting artefact;
%   3. the fast component's weight w is more than 3 sigma from zero AND more than 3 sigma from
%      one. Either bound means the second component is not there.
%
% All three are reported pass or fail, so a refused second component is visible. The default is
% 'stretched' rather than 'auto', because a two-component T1 is a claim worth making
% deliberately.
%
% ===================================================================
% THE LONG END IS WHERE A T1 FIT GOES WRONG
% ===================================================================
%
% Not the short end. A T1 sweep that stops before the plateau leaves the baseline c_i and the
% amplitude a_i free to trade against T1 along a ridge, and the fit will pick a point on that
% ridge and quote an error bar for it. The symptoms are a T1 comparable to or larger than the
% sweep and an error bar that looks fine.
%
% Two numbers are reported for it. out.spanOverT1 is how many T1s the sweep covers -- under
% about 2 and the plateau is an extrapolation. out.plateauSigma is how far the LAST quarter of the
% POINTS sits from the fitted plateau, in sigma of that mean.
%
% AND THE ERROR BAR DOES NOT KNOW. This is the part worth stating precisely, because the obvious
% assumption is wrong in a way that makes a truncated T1 look better than it is. Measured in
% testT1Fit's t4 over ten noise realisations of a 2500 us T1:
%
%       reach       median T1    quoted error    actual scatter    ratio
%       0.8 T1        2399           561              1635          2.9
%       1.6 T1        2473           378               294          0.8
%       3.6 T1        2465           107                75          0.7
%       8.0 T1        2517            71                56          0.8
%
% The fitted VALUE is barely biased even at 0.8 T1 -- 4% low. What fails is the uncertainty: the
% quoted error understates the real scatter by about three times, because the likelihood has a
% long ridge along which the plateau, the amplitude and T1 trade, and an unlucky draw lands far
% along that ridge while still reporting the local curvature it happens to sit in. Past about
% 1.6 T1 of reach the quoted error tracks the scatter properly.
%
% So spanOverT1 is the number to act on, and it is warned about separately for exactly this
% reason. A T1 with a 20% error bar on a sweep reaching 0.8 T1 is not a 20% measurement.
%
% ===================================================================
% ROBUSTNESS
% ===================================================================
%
% Tukey bisquare re-weighting (IRLS, 4.685 sigma), same as the other four. On a T1 the points
% that need it are the ones the classic path dilutes toward no contrast, and on this experiment
% they are especially hard to see by eye: a diluted point sits at the PLATEAU, which is where
% half the sweep already is, so it hides in the flat part instead of standing out. The report
% names how many were down-weighted and at which tau.
%
% ===================================================================
% OPTIONS  (defaults in brackets)
% ===================================================================
%
%   exp []                  an AWGExpT1, or
%   t, y, e []              vectors; y may be 1xN or 2xN
%   model ['stretched']     'stretched' | 'double' | 'auto'
%   pFixed []               fix the stretch exponent (1 for a plain exponential)
%   useProjection ['both']  'both' | 'one' | 'zero'
%   robust [true]           Tukey bisquare re-weighting
%   tMin [], tMax []        window in tau
%   nGridT1 [30]            T1 grid points
%   avgOutlierMad [4]       MAD cut on whole averages, applied ONCE
%   plotIt [true], verbose [true], showResid [false], logX [true], yLabel ['']
%
% ===================================================================
% OUT
% ===================================================================
%
%   ok, why                 did it work, and if not why
%   T1, T1Err               microseconds
%   ratePerUs, rateErr      1/T1. THE additive quantity -- comparable with t1GateTest's r
%   p, pErr, pMeaning       stretch exponent and what it implies
%   coherenceArea           T1*gamma(1+1/p), the shape-independent figure of merit
%   contrast, contrastErr   |a| over the plateau, for the |1> trace
%   flatnessSigma           a(|0>) over its error. THE CONTROL. Should be under 3.
%   doubleUsed, T1fast, T1slow, weight, doubleTest
%   spanOverT1, plateauSigma, decayErrorBars
%   residAutocorr, chi2Reduced, r2, r2Ceiling, r2Kept, r2CeilingKept
%   plateau, amp                per projection
%   t, y, e, model, resid, weights, nDown, downIdx, downT, traces, labels
%   gateOpenUs              the gate dose the run was taken at, if the object carried it
%   coeffError              [T1Err, contrastErr, pErr]

    DEFAULTS = struct('exp', [], 't', [], 'y', [], 'e', [], ...
        'model', 'stretched', 'pFixed', [], 'useProjection', 'both', 'robust', true, ...
        'tMin', [], 'tMax', [], 'nGridT1', 30, 'avgOutlierMad', 4, ...
        'plotIt', true, 'verbose', true, 'showResid', false, 'logX', true, 'yLabel', '');

    if ~exist('args', 'var') || isempty(args); args = struct(); end
    if ~isstruct(args)
        error(['t1FitRobust takes ONE struct: t1FitRobust(struct(''exp'', t)) or ' ...
               't1FitRobust(struct(''t'', tau, ''y'', y, ''e'', e)).']);
    end
    opts = withDefaults(args, DEFAULTS, 't1FitRobust');
    if ~ismember(lower(opts.model), {'stretched', 'double', 'auto'})
        error('t1FitRobust: model must be ''stretched'', ''double'' or ''auto'', not ''%s''.', opts.model);
    end
    opts.model = lower(opts.model);

    out = struct('ok', false, 'why', '', 'opts', opts);
    Exp = opts.exp; fromVectors = isempty(Exp);
    gateOpenUs = NaN; pullInfo = [];

    if fromVectors
        if isempty(opts.t) || isempty(opts.y)
            out.why = 'nothing to fit: give exp, or t and y.';
            fprintf('\nt1FitRobust FAILED: %s\n', out.why); return
        end
        t = opts.t; Y = opts.y; E = opts.e;
    else
        try
            [t, Y, E, pullInfo] = pullT1(Exp, opts);
            gateOpenUs = pullInfo.gateOpenUs;
        catch err
            describeExp(Exp);
            rethrow(err);
        end
    end

    t = t(:)';
    if isvector(Y); Y = Y(:)'; end
    if isvector(E); E = E(:)'; end
    if size(Y, 2) ~= numel(t)
        if size(Y, 1) == numel(t); Y = Y.'; E = E.'; else
            out.why = sprintf('%d taus against a %dx%d signal', numel(t), size(Y,1), size(Y,2));
            fprintf('\nt1FitRobust FAILED: %s\n', out.why); return
        end
    end
    if isempty(E) || all(~isfinite(E(:))); E = ones(size(Y)); end
    if ~isequal(size(E), size(Y)); E = repmat(median(E(isfinite(E) & E>0)), size(Y)); end

    [Y, E, labels, ctlRow] = pickProjections(Y, E, opts.useProjection);

    tMin = opts.tMin; if isempty(tMin); tMin = 0; end
    tMax = opts.tMax; if isempty(tMax); tMax = Inf; end
    keep = isfinite(t) & t >= tMin & t <= tMax & all(isfinite(Y), 1);
    nDrop = sum(~keep);
    t = t(keep); Y = Y(:, keep); E = E(:, keep);
    n = numel(t); nT = size(Y, 1);
    if n < 6
        out.why = sprintf(['only %d usable points. A T1 needs enough to separate a decay from a ' ...
                           'plateau -- 12 or more, reaching past the decay.'], n);
        fprintf('\nt1FitRobust FAILED: %s\n', out.why); return
    end
    E(~isfinite(E) | E <= 0) = median(E(isfinite(E) & E > 0));
    span = max(t) - min(t);

    % ================================================================
    % THE FIT
    % ================================================================
    fitP = isempty(opts.pFixed); pFix = opts.pFixed;
    w0 = 1 ./ E.^2;

    S = fitStretched(t, Y, w0, span, opts, fitP, pFix);
    P = S.P; kind = 'stretched';
    test = struct('tried', false, 'dSSR', NaN, 'threshold', 9, 'ratio', NaN, ...
                  'wSigma', NaN, 'failed', '', 'accepted', false);

    % ---- the second component, on PLAIN weights ---------------------------
    % Before any robust re-weighting, for the same reason echoFitRobust compares its comb on
    % plain weights: Tukey weights derived from one of two competing models delete the points
    % that distinguish them, and the null hypothesis should not get to choose the data.
    if ~strcmp(opts.model, 'stretched')
        test.tried = true;
        [~, ~, ssrS] = chiOf(P, t, Y, w0, 'stretched');
        Dbl = fitDouble(t, Y, w0, span, opts, P);
        if isempty(Dbl)
            test.failed = 'no two-component fit converged on this span';
        else
            [~, linD, ssrD] = chiOf(Dbl.P, t, Y, w0, 'double');
            test.dSSR = ssrS - ssrD;
            test.ratio = max(Dbl.P(1), Dbl.P(2)) / max(min(Dbl.P(1), Dbl.P(2)), eps);
            [errD, ~] = paramErrors(Dbl.P, linD, t, Y, w0, 'double', true, 1);
            test.wSigma = min(abs(Dbl.P(3)), abs(1 - Dbl.P(3))) / max(errD.w, eps);
            fails = {};
            if ~(test.dSSR > test.threshold)
                if test.dSSR <= 0
                    fails{end+1} = sprintf('the two-component fit is WORSE (SSR up by %.1f)', -test.dSSR);
                else
                    fails{end+1} = sprintf('SSR dropped by only %.1f, under the %.0f needed for one extra parameter', ...
                                           test.dSSR, test.threshold);
                end
            end
            if ~(test.ratio > 2)
                fails{end+1} = sprintf(['the two time constants differ by only %.2fx (needs 2) -- ' ...
                                        'that is one component and a fitting artefact'], test.ratio);
            end
            if ~(test.wSigma > 3)
                fails{end+1} = sprintf(['the weight is %.1f sigma from a bound (needs 3) -- one of ' ...
                                        'the two components is not there'], test.wSigma);
            end
            if isempty(fails) || strcmp(opts.model, 'double')
                P = Dbl.P; kind = 'double'; test.accepted = true;
                if ~isempty(fails)
                    test.failed = ['FORCED by model = ''double'' despite: ' strjoin(fails, '; ')];
                end
            else
                test.failed = strjoin(fails, '; ');
            end
        end
    end

    % ---- robust re-weighting, ITERATED, on the model that was chosen -------
    % Three IRLS passes rather than one, and the reason is worth writing down because one pass
    % looked fine until it was inspected.
    %
    % The weights come from the residuals of the CURRENT fit, and with a single pass that is the
    % unweighted fit -- which is the one the outlier has already distorted. On a 24-point log
    % sweep with one point diluted to no contrast, that first fit is pulled hard enough that the
    % whole short-tau end sits more than 3.3 sigma out, so ten points get flagged instead of one.
    % The refit was fine (chi2 0.90, R^2 0.999 on the kept points) but the REPORT said ten points
    % were down-weighted, which reads as a broken run.
    %
    % Iterating fixes it: after the first refit the model is right, the residuals are honest, and
    % the second pass flags the one point that is actually bad. Three passes is the standard
    % number for a bisquare and it converges well before that here.
    w = w0; nDown = 0; downIdx = [];
    if opts.robust
        for pass = 1:3
            [~, lin] = chiOf(P, t, Y, w, kind);
            R = (Y - modelOf(P, t, lin, kind)) ./ E;
            s = 1.4826 * median(abs(R(:) - median(R(:))));
            if ~(s > 0); break; end
            u = R / (4.685 * s);
            bw = (abs(u) < 1) .* (1 - u.^2).^2;
            % bw < 0.25 is |residual| > 3.31 sigma. The WEIGHTING tapers smoothly from zero
            % either way; this number only decides what gets NAMED an outlier and ringed.
            newIdx = find(any(bw < 0.25, 1));
            w = w0 .* max(bw, 1e-3);
            if strcmp(kind, 'stretched')
                S = fitStretched(t, Y, w, span, opts, fitP, pFix);
                P = S.P;
            else
                P = nelderMead(@(Q) chiOf(Q, t, Y, w, 'double'), P, boundsFor('double', span, fitP, pFix));
            end
            converged = isequal(newIdx, downIdx);
            downIdx = newIdx; nDown = numel(downIdx);
            if converged; break; end
        end
    end

    [chi2, lin] = chiOf(P, t, Y, w, kind);
    M = modelOf(P, t, lin, kind);
    resid = Y - M;

    % ================================================================
    % PARAMETERS
    % ================================================================
    o = struct();
    o.kind = kind; o.doubleUsed = strcmp(kind, 'double'); o.doubleTest = test;
    o.T1fast = NaN; o.T1slow = NaN; o.weight = NaN;
    if o.doubleUsed
        Tf = min(P(1), P(2)); Ts = max(P(1), P(2));
        wf = P(3); if P(1) > P(2); wf = 1 - P(3); end
        o.T1fast = Tf; o.T1slow = Ts; o.weight = wf;
        % The reported T1 for a two-component decay is the AMPLITUDE-WEIGHTED mean lifetime,
        % which is the integral of D and therefore the only single number that reproduces the
        % total area. Quoting the slow constant alone is the usual shortcut and it overstates
        % the coherence by however much weight sits in the fast one.
        o.T1 = wf*Tf + (1-wf)*Ts;
        o.p = 1;
    else
        o.T1 = P(1); o.p = P(2);
    end
    o.plateau = zeros(1, nT); o.amp = zeros(1, nT);
    for i = 1:nT
        o.plateau(i) = lin(i, 1);
        o.amp(i)     = lin(i, 2);
    end
    % The |1> trace is the one that moves. With both projections present that is the second row
    % once pickProjections has ordered them; with one, it is the only row.
    sigRow = nT; if nT == 1; sigRow = 1; end
    o.contrast = abs(o.amp(sigRow)) / max(abs(o.plateau(sigRow)), eps);

    [err, cov] = paramErrors(P, lin, t, Y, w, kind, fitP, chi2);
    B = boundsFor(kind, span, fitP, pFix);
    o.T1Err = pickErr(err.T1, profileError(1, P, t, Y, w, kind, B), 20*span);
    o.pErr  = NaN;
    if ~o.doubleUsed && fitP
        o.pErr = pickErr(err.p, profileError(2, P, t, Y, w, kind, B), 3);
    end
    o.contrastErr = err.contrast;
    o.cov = cov;

    % RATES. 1/T1 is what adds when relaxation channels combine, so it is the quantity to
    % compare between runs and the one t1GateTest's leak rate is measured in.
    o.ratePerUs = 1/max(o.T1, eps);
    o.rateErr   = o.T1Err / max(o.T1^2, eps);
    o.coherenceArea = o.T1 * gamma(1 + 1/max(o.p, 0.05));
    o.pMeaning = pMeaningOf(o.p, o.pErr, fitP, o.doubleUsed);

    % ---- THE CONTROL: is the |0> trace flat -------------------------------
    o.flatnessSigma = NaN; o.controlPresent = false; o.controlTestable = false;
    if nT == 2 && ~isempty(ctlRow)
        o.controlPresent = true;
        ampErr = err.ampPerTrace(ctlRow);
        % GUARD THE DIVISION. max(NaN, eps) returns eps -- max ignores NaN -- so a
        % non-computable amplitude error turned this into abs(amp)/2.2e-16 and printed a
        % flatness of 1.07e17 sigma, with "*** THE CONTROL IS MOVING ***" under it. On a run
        % whose fit was degenerate in every direction. That is the worst kind of wrong: a
        % fabricated alarm on data that could not support any statement at all, and it looked
        % like the most confident line in the report.
        %
        % An amplitude error that is not computable is itself the finding -- it means the fit
        % has no curvature in that direction -- so the control is reported as NOT TESTABLE
        % rather than as passing or failing.
        if isfinite(ampErr) && ampErr > 0
            o.flatnessSigma = abs(o.amp(ctlRow)) / ampErr;
            o.controlTestable = true;
        end
    end

    % ---- did the sweep reach the plateau ---------------------------------
    o.spanOverT1 = span / max(o.T1, eps);
    % The last quarter of the POINTS, not of the range. A T1 sweep is log-spaced over two or
    % three decades, so the last quarter of the range holds one or two points and the check
    % silently never ran -- which is worse than a wrong answer, because nothing was printed.
    nLate = max(3, ceil(n/4));
    late = false(1, n); late(max(1, n-nLate+1):n) = true;
    o.plateauSigma = NaN;
    if sum(late) >= 3
        r1 = resid(sigRow, late); e1 = E(sigRow, late);
        wv = 1./e1.^2;
        mu = sum(wv.*r1)/sum(wv);
        o.plateauSigma = mu / max(sqrt(1/sum(wv)), eps);
    end
    o.decayErrorBars = (max(M(sigRow,:)) - min(M(sigRow,:))) / max(median(E(sigRow,:)), eps);

    o.chi2Reduced = chi2;
    yAll = Y(:); eAll = E(:); rAll = resid(:);
    o.r2 = 1 - sum(rAll.^2)/max(sum((yAll - mean(yAll)).^2), eps);
    o.r2Ceiling = 1 - sum(eAll.^2)/max(sum((yAll - mean(yAll)).^2), eps);
    o.r2Kept = o.r2; o.r2CeilingKept = o.r2Ceiling;
    if nDown > 0
        m = true(1, n); m(downIdx) = false;
        rk = resid(:, m); yk = Y(:, m); ek = E(:, m);
        o.r2Kept = 1 - sum(rk(:).^2)/max(sum((yk(:) - mean(yk(:))).^2), eps);
        o.r2CeilingKept = 1 - sum(ek(:).^2)/max(sum((yk(:) - mean(yk(:))).^2), eps);
    end
    o.residAutocorr = lag1(resid ./ max(E, eps));

    % IS THERE A MEASUREMENT HERE AT ALL. Three independent ways for the answer to be no, and
    % any one of them is enough: the fitted decay does not clear the noise, the contrast is not
    % separable from zero, or T1 has no computable error bar. When this is false the report and
    % the figure LEAD with it, because the first version led with "T1 = 9.798e+04 +/- NaN us"
    % on a run that measured nothing -- putting the least trustworthy number on the page in the
    % largest type.
    o.measured = o.decayErrorBars >= 5 && isfinite(o.T1Err) && ...
                 isfinite(o.contrastErr) && o.contrastErr < 0.5*o.contrast;
    % And the single most informative statement about a run like that, which a negative R^2 says
    % only obliquely: sum(sterr^2) exceeding the total variance of the data means the quoted
    % error bars are BIGGER than everything that happens in the trace.
    o.errorBarsExceedSpread = o.r2Ceiling < 0;

    o.t = t; o.y = Y; o.e = E; o.model = M; o.resid = resid;
    o.weights = w; o.nDown = nDown; o.downIdx = downIdx; o.downT = t(downIdx);
    o.labels = {labels{:}};                                              %#ok<CCAT1>
    o.traces = nT; o.spanUs = span; o.tMin = tMin; o.tMax = tMax; o.nDropped = nDrop;
    o.ctlRow = ctlRow; o.sigRow = sigRow; o.gateOpenUs = gateOpenUs;
    o.pFitted = fitP && ~o.doubleUsed;
    o.ok = true; o.why = '';
    o.coeffError = [o.T1Err, o.contrastErr, o.pErr];

    out = o; out.opts = opts;
    if ~isempty(pullInfo); out.pull = pullInfo; end

    % For the experiment window (ExpFit): the curves the figure draws (one row per
    % projection) and the figure's title, so the window shows the same fit as this figure.
    oc = out; oc.opts.exp = [];
    out.curveAt = @(tt) t1CurveAt(oc, tt);
    try; out.titleLines = titleLines(out); catch; out.titleLines = {}; end

    if opts.verbose; report(out); end
    if opts.plotIt;  drawIt(out); end
end

function M = t1CurveAt(o, tt)
% The fitted model at any tau, one row per projection -- what drawIt draws.
    M = modelOf(paramVec(o), reshape(tt, 1, []), [o.plateau(:), o.amp(:)], o.kind);
end

% ======================================================================
% MODEL
% ======================================================================

function D = decayOf(P, t, kind)
    tt = max(t, 0);
    if strcmp(kind, 'double')
        Ta = max(P(1), eps); Tb = max(P(2), eps); wf = min(max(P(3), 0), 1);
        D = wf*exp(-tt/Ta) + (1-wf)*exp(-tt/Tb);
    else
        D = exp(-(tt/max(P(1), eps)).^max(P(2), 0.05));
    end
end

function X = designAt(P, t, kind)
% [1, D]. The plateau and the amplitude are linear in the model, so they are solved exactly at
% every step of the search rather than being given a starting guess -- which is what stops the
% amplitude and the baseline settling into one of the many pairs that describe a truncated decay
% equally well.
    D = decayOf(P, t, kind);
    X = [ones(numel(t),1), D(:)];
end

function [chi2, lin, ssr] = chiOf(P, t, Y, w, kind)
    X = designAt(P, t, kind);
    nT = size(Y, 1);
    lin = zeros(nT, 2); s = 0; ndof = 0;
    for i = 1:nT
        wi = w(i,:)';
        [th, nUsed] = wls(X, Y(i,:)', wi);
        lin(i,:) = th';
        r = Y(i,:)' - X*th;
        s = s + sum(wi .* r.^2);
        ndof = ndof + numel(t) - nUsed;
    end
    chi2 = s / max(ndof, 1);
    ssr = s;
end

function M = modelOf(P, t, lin, kind)
    X = designAt(P, t, kind);
    M = zeros(size(lin,1), numel(t));
    for i = 1:size(lin,1)
        M(i,:) = (X * lin(i,:)')';
    end
end

function r = lag1(R)
    r = 0; num = 0; den = 0;
    for i = 1:size(R,1)
        v = R(i,:) - mean(R(i,:));
        if numel(v) < 3; continue; end
        num = num + sum(v(1:end-1).*v(2:end));
        den = den + sum(v.^2);
    end
    if den > 0; r = num/den; end
end

function s = pMeaningOf(p, pErr, fitP, isDouble)
    if isDouble
        s = 'not applicable: a two-component fit was used instead'; return
    end
    if ~fitP
        s = 'held fixed -- an assumption, not a measurement'; return
    end
    if p < 0.35;  s = 'ON ITS LOWER BOUND -- a very wide spread of rates, or the model is wrong'; return; end
    if p > 2.95;  s = 'ON ITS UPPER BOUND -- for a T1 that is hard to arrange physically; suspect a common-mode drift'; return; end
    if ~isfinite(pErr) || pErr <= 0
        s = 'no usable error bar, so it is not distinguishable from anything'; return
    end
    if abs(p - 1) < 2*pErr
        s = 'consistent with 1: a single relaxation rate, the clean case';
    elseif p < 1
        s = sprintf('%.1f sigma BELOW 1: a distribution of rates, not one rate', (1-p)/pErr);
    else
        s = sprintf('%.1f sigma ABOVE 1: slower than exponential early, which for a T1 wants explaining', (p-1)/pErr);
    end
end

% ======================================================================
% SEARCH
% ======================================================================

function S = fitStretched(t, Y, w, span, opts, fitP, pFix)
% Grid then polish over (T1, p). The T1 grid is LOG spaced across two and a half decades around
% the span, because T1 is a scale and a linear grid spends nearly all its points at the long end
% where the model barely changes.
    T1Grid = logspace(log10(0.03*span), log10(10*span), max(opts.nGridT1, 10));
    if fitP; pGrid = [0.4 0.6 0.8 1 1.3 1.7 2.2]; else; pGrid = pFix; end
    best = struct('chi2', Inf, 'P', [T1Grid(1), pGrid(1)]);
    for ip = 1:numel(pGrid)
        for it = 1:numel(T1Grid)
            c = chiOf([T1Grid(it), pGrid(ip)], t, Y, w, 'stretched');
            if c < best.chi2; best.chi2 = c; best.P = [T1Grid(it), pGrid(ip)]; end
        end
    end
    B = boundsFor('stretched', span, fitP, pFix);
    S = struct('P', nelderMead(@(Q) chiOf(Q, t, Y, w, 'stretched'), best.P, B), 'grid', best);
end

function D = fitDouble(t, Y, w, span, opts, P0)
% Two time constants and a weight. Seeded from a COARSE two-dimensional grid over the pair,
% restricted to Ta < Tb, plus three weights -- rather than from the stretched fit, because the
% stretched T1 sits between the two constants it is standing in for and starting both there is
% the one place Nelder-Mead reliably collapses them onto each other.
    D = [];
    g = logspace(log10(0.03*span), log10(10*span), 12);
    wGrid = [0.25 0.5 0.75];
    best = struct('chi2', Inf, 'P', []);
    for ia = 1:numel(g)
        for ib = ia+1:numel(g)
            for iw = 1:numel(wGrid)
                Q = [g(ia), g(ib), wGrid(iw)];
                c = chiOf(Q, t, Y, w, 'double');
                if c < best.chi2; best.chi2 = c; best.P = Q; end
            end
        end
    end
    if isempty(best.P); return; end
    B = boundsFor('double', span, true, []);
    D = struct('P', nelderMead(@(Q) chiOf(Q, t, Y, w, 'double'), best.P, B), 'grid', best);
end

function B = boundsFor(kind, span, fitP, pFix)
    if strcmp(kind, 'double')
        B = [0.01*span, 20*span; 0.01*span, 20*span; 0.01, 0.99];
    else
        B = [0.01*span, 20*span; 0.3, 3];
        if ~fitP
            if ~exist('pFix', 'var') || isempty(pFix); pFix = 1; end
            B(2,:) = [pFix pFix];
        end
    end
end

function [th, nUsed] = wls(X, y, wi)
%WLS Weighted least squares that survives a degenerate column.
%
% The degeneracy here is specific and common: when T1 wanders far past the sweep, exp(-(t/T1)^p)
% is 1 to working precision at every point, so the decay column becomes a copy of the constant
% column. The normal equations are then singular in that direction and the solve returns an
% enormous amplitude with an equally enormous opposite plateau -- a fit with several hundred per
% cent contrast that reproduces the data perfectly and means nothing. On a T1 that is not an edge
% case, it is what happens whenever the sweep stops early.
    sc = sqrt(sum((X.^2) .* wi, 1));
    keep = sc > 1e-9 * max(sc);
    th = zeros(size(X,2), 1);
    nUsed = sum(keep);
    if nUsed == 0; return; end
    Xk = X(:, keep);
    M = Xk' * (Xk .* wi);
    b = Xk' * (y .* wi);
    if rcond(M) < 1e-14; return; end
    th(keep) = M \ b;
end

function P = nelderMead(fun, P0, B)
% Written out rather than calling fminsearch, so this file behaves identically in Octave, where
% it is parse-checked, and so there is no question about which toolbox anything came from.
    n = numel(P0);
    P0 = clip(P0, B);
    S = repmat(P0, n+1, 1);
    for i = 1:n
        step = 0.08*max(abs(P0(i)), 0.05*(B(i,2) - B(i,1)));
        S(i+1,i) = clipOne(P0(i) + step, B(i,:));
    end
    F = zeros(n+1,1);
    for i = 1:n+1; F(i) = fun(S(i,:)); end
    for iter = 1:500
        [F, ord] = sort(F); S = S(ord,:);
        if abs(F(end) - F(1)) <= 1e-10*(abs(F(1)) + 1e-12); break; end
        m = mean(S(1:n,:), 1);
        xr = clip(m + (m - S(end,:)), B); fr = fun(xr);
        if fr < F(1)
            xe = clip(m + 2*(m - S(end,:)), B); fe = fun(xe);
            if fe < fr; S(end,:) = xe; F(end) = fe; else; S(end,:) = xr; F(end) = fr; end
        elseif fr < F(n)
            S(end,:) = xr; F(end) = fr;
        else
            xc = clip(m + 0.5*(S(end,:) - m), B); fc = fun(xc);
            if fc < F(end)
                S(end,:) = xc; F(end) = fc;
            else
                for i = 2:n+1
                    S(i,:) = clip(S(1,:) + 0.5*(S(i,:) - S(1,:)), B);
                    F(i) = fun(S(i,:));
                end
            end
        end
    end
    [~, k] = min(F);
    P = S(k,:);
end

function x = clip(x, B)
    for i = 1:numel(x); x(i) = clipOne(x(i), B(i,:)); end
end

function v = clipOne(v, b)
    v = min(max(v, b(1)), b(2));
end

% ======================================================================
% ERRORS
% ======================================================================

function [err, C] = paramErrors(P, lin, t, Y, w, kind, fitP, chi2)
% Numerical Jacobian over the nonlinear parameters AND every linear one, then C = inv(J'WJ),
% Birge scaled where the reduced chi2 exceeds 1.
%
% The linear block has to be in there. On a T1 the correlation between T1, the amplitude and the
% plateau is not a nuisance, it is the dominant feature of the likelihood -- a smaller amplitude
% with a longer T1 describes a truncated decay almost as well -- so an error bar computed over
% the nonlinear parameters alone comes out far too small. That is the shape of mistake that makes
% a fit look better than it is.
    err = struct('T1', NaN, 'p', NaN, 'w', NaN, 'contrast', NaN, 'ampPerTrace', []);
    C = [];
    nT = size(Y,1); nP = numel(P); nLin = 2;
    theta = [P(:); reshape(lin', [], 1)];
    nTot = numel(theta);
    f = @(th) flatModel(th, nP, nT, nLin, t, kind);
    y0 = f(theta);
    J = zeros(numel(y0), nTot);
    for k = 1:nTot
        h = 1e-6*max(abs(theta(k)), 1e-6);
        tp = theta; tp(k) = tp(k) + h;
        J(:,k) = (f(tp) - y0)/h;
    end
    wv = reshape(w', [], 1);
    sc = sqrt(sum((J.^2) .* wv, 1));
    good = sc > 1e-9*max(sc);
    % A FIXED p is not a free direction. Leaving it in would hand its variance to T1, which is
    % the wrong sign of error: fixing a parameter must make the others' bars smaller.
    if strcmp(kind, 'stretched') && ~fitP && nP >= 2; good(2) = false; end
    s = nan(nTot, 1);
    if any(good)
        Jk = J(:, good);
        Mk = Jk' * (Jk .* wv);
        if rcond(Mk) > 1e-14
            Ck = inv(Mk) * max(chi2, 1);
            C = nan(nTot); C(good, good) = Ck;
            s(good) = sqrt(abs(diag(Ck)));
        end
    end
    if strcmp(kind, 'double')
        % The reported T1 is the amplitude-weighted mean w*Ta + (1-w)*Tb, so its error has to be
        % propagated through that rather than taken from either constant.
        if ~isempty(C)
            g = zeros(1, nTot);
            g(1) = P(3); g(2) = 1 - P(3); g(3) = P(1) - P(2);
            Cg = C; Cg(~isfinite(Cg)) = 0;
            err.T1 = sqrt(abs(g * Cg * g'));
        end
        err.w = s(3);
    else
        err.T1 = s(1);
        if fitP; err.p = s(2); end
    end
    err.ampPerTrace = nan(1, nT);
    for i = 1:nT
        err.ampPerTrace(i) = s(nP + nLin*(i-1) + 2);
    end
    % contrast = |a|/|c0| for the signal trace, which is the LAST one.
    i0 = nP + nLin*(nT-1);
    c0 = lin(nT,1); a = lin(nT,2);
    if abs(c0) > 0 && ~isempty(C)
        g = zeros(1, nTot);
        g(i0+1) = -abs(a)*sign(c0)/c0^2;
        g(i0+2) =  sign(a)/abs(c0);
        Cg = C; Cg(~isfinite(Cg)) = 0;
        err.contrast = sqrt(abs(g * Cg * g'));
    end
end

function v = flatModel(theta, nP, nT, nLin, t, kind)
    P = theta(1:nP)';
    lin = reshape(theta(nP+1:end), nLin, nT)';
    M = modelOf(P, t, lin, kind);
    v = reshape(M', [], 1);
end

function e = pickErr(jac, prof, range)
    cand = [];
    if isfinite(jac)  && jac  > 0 && jac  < 0.5*range; cand(end+1) = jac;  end
    if isfinite(prof) && prof > 0 && prof < 0.5*range; cand(end+1) = prof; end
    if isempty(cand); e = NaN; else; e = max(cand); end
end

function e = profileError(k, P, t, Y, w, kind, B)
%PROFILEERROR 1-sigma from scanning one parameter: where the weighted SSR rises by 1.
%
% EXPAND TO BRACKET, THEN BISECT. The obvious implementation -- step out on a geometric ladder
% and return the first rung past the crossing -- is wrong in one direction only, which is the
% dangerous kind: the answer can only be too big. echoFitRobust shipped with that version and it
% returned the same 0.9431 us for every seed and for one projection and two, because the first
% rung was already past the crossing and the "error bar" was simply 2% of the box width. The true
% scatter was 0.2559. A constant error bar does not shrink with more averages and does not shrink
% when a second trace is added, which destroys the two things it is used for.
%
% ramseyFitRobust, rabiFitRobust and esrFitRobust still carry the ladder version.
    e = NaN;
    if k > numel(P); return; end
    if B(k,1) == B(k,2); return; end
    [~, ~, ssr0] = chiOf(P, t, Y, w, kind);
    if ~isfinite(ssr0); return; end
    lo = B(k,1); hi = B(k,2); p0 = P(k);
    for side = [1 -1]
        if side > 0; reach = hi - p0; else; reach = p0 - lo; end
        if ~(reach > 0); continue; end
        step = max(1e-4*max(abs(p0), eps), 1e-7*(hi - lo));
        below = 0; above = NaN;
        for it = 1:90
            d = min(step, reach);
            Q = P; Q(k) = p0 + side*d;
            [~, ~, ssr] = chiOf(Q, t, Y, w, kind);
            if ssr - ssr0 >= 1; above = d; break; end
            below = d;
            if d >= reach; break; end
            step = step*1.6;
        end
        if ~isfinite(above); continue; end
        a = below; b = above;
        for it = 1:50
            mid = 0.5*(a+b);
            Q = P; Q(k) = p0 + side*mid;
            [~, ~, ssr] = chiOf(Q, t, Y, w, kind);
            if ssr - ssr0 >= 1; b = mid; else; a = mid; end
            if (b - a) <= 1e-4*max(b, eps); break; end
        end
        v = 0.5*(a+b);
        if v > 0
            if isnan(e); e = v; else; e = max(e, v); end
        end
    end
end

% ======================================================================
% DATA
% ======================================================================

function [t, Y, E, info] = pullT1(Exp, opts)
%PULLT1 Get tau and both traces out of an AWGExpT1, in a fixed order: |0> first, |1> second.
%
% THE ORDER IS DECIDED BY t1Contrast, not by the field names. AWGExpT1's projection labels were
% the wrong way round until they were fixed, and this fitter has to give the same answer on data
% taken before the fix as after it. Reading the roles off signalParam and signalParam2 would put
% the flat control where the decay belongs, and the fit would then report the control's noise as
% a T1 with an error bar.
%
% Deliberately NOT extract1Ddata: that applies three successive cuts across averages at 5, then
% 3, then ONE standard deviation, and a one-sigma cut discards about a third of perfectly good
% Gaussian data after two harsher passes have already shrunk the spread it measures against. All
% four of the other robust fits refuse it and so does this.
    info = struct('nAvg', NaN, 'repeats', NaN, 'usedScatter', false, 'flipped', false, ...
                  'gateOpenUs', NaN, 'mwGateMode', '', 'confident', false, 'why', '');

    t = tryGet(Exp, 'tau', []);
    ax = tryGet(Exp, 'mCurrentXAxisParam', []);
    if ~isempty(ax)
        v = tryGet(ax, 'value', []);
        if ~isempty(v); t = v; end
    end
    if isempty(t); error('t1FitRobust: the experiment has no tau axis.'); end
    t = t(:)';

    info.repeats = tryGet(Exp, 'repeats', NaN);
    info.nAvg    = tryGet(Exp, 'currIter', tryGet(Exp, 'averages', NaN));
    info.gateOpenUs = tryGet(Exp, 'mwGateOpenUs', NaN);
    info.mwGateMode = tryGet(Exp, 'mwGateMode', '');

    [y1, e1] = paramData(Exp, 'signalParam');
    [y2, e2] = paramData(Exp, 'signalParam2');
    if isempty(y1); error('t1FitRobust: signalParam is empty -- did the run finish?'); end

    if isempty(y2) || numel(y2) ~= numel(y1)
        % One trace only. That happens with doubleMeasurement off, where AWGExpT1 runs the pi
        % half alone -- so what is here IS the |1> relaxation, and there is no control.
        Y = y1; E = e1;
    else
        c = t1Contrast(struct('value', y1, 'sterr', e1, 'value2', y2, 'sterr2', e2, 'tau', t));
        info.flipped = c.flipped; info.confident = c.confident; info.why = c.why;
        if c.flipped
            Y = [y2; y1]; E = [e2; e1];       % |0> first, |1> second
        else
            Y = [y1; y2]; E = [e1; e2];
        end
    end

    % Widen the error bars using the scatter across averages where it is bigger. Two different
    % questions -- photon statistics within an average, and everything that drifts between them
    % -- and on this setup the second is routinely the bigger. t1GateTest measures the factor.
    S = tryGet(Exp, 'signal', []);
    nIter = tryGet(Exp, 'currIter', NaN);
    if ~isempty(S) && ndims(S) == 3 && isfinite(nIter) && nIter >= 3
        try
            for row = 1:size(Y,1)
                src = row; if info.flipped; src = 3 - row; end
                dS = 2*src - 1; dR = 2*src;
                if size(S,1) >= dR
                    r = squeeze(S(dS,:,1:nIter)) ./ max(squeeze(S(dR,:,1:nIter)), eps);
                    sc = std(r, 0, 2)' / sqrt(nIter);
                    bigger = sc > E(row,:);
                    if any(bigger)
                        E(row, bigger) = sc(bigger);
                        info.usedScatter = true;
                    end
                end
            end
        catch
        end
    end
end

function [v, s] = paramData(e, which)
    v = []; s = [];
    p = tryGet(e, which, []);
    if isempty(p); return; end
    v = tryGet(p, 'value', []); s = tryGet(p, 'sterr', []);
    v = v(:)'; s = s(:)';
end

function v = tryGet(o, name, dflt)
    v = dflt;
    try
        if isobject(o) || isstruct(o); v = o.(name); end
    catch
    end
    if isempty(v) && ~isempty(dflt); v = dflt; end
end

function [Y, E, labels, ctlRow] = pickProjections(Y, E, which)
% Rows arrive as |0> then |1>. The CONTROL is the |0> row, and knowing which row that is is what
% makes the flatness check possible.
    labels = {'|1\rangle'}; ctlRow = [];
    if size(Y,1) >= 2
        labels = {'|0\rangle  (control)', '|1\rangle'};
        ctlRow = 1;
    end
    switch lower(which)
        case 'one'
            if size(Y,1) >= 2; Y = Y(2,:); E = E(2,:); end
            labels = {'|1\rangle'}; ctlRow = [];
        case 'zero'
            if size(Y,1) >= 2; Y = Y(1,:); E = E(1,:); end
            labels = {'|0\rangle  (control)'}; ctlRow = [];
        otherwise
    end
end

% ======================================================================
% REPORT
% ======================================================================

function report(o)
    fprintf('\n=== t1FitRobust ===\n');
    if o.doubleUsed
        fprintf('model: TWO exponentials, %d trace(s)\n', o.traces);
    else
        fprintf('model: stretched exponential, %d trace(s)\n', o.traces);
    end

    if ~o.measured
        fprintf('\n*** NO MEASURABLE DECAY IN THIS RUN. ***\n');
        fprintf(['    Read the reasons below and nothing else. T1, p, the rate and the area are\n' ...
                 '    printed because the fit produced them, not because they mean anything: they\n' ...
                 '    describe wherever the optimiser stopped on a flat likelihood.\n']);
        why = {};
        if o.decayErrorBars < 5
            why{end+1} = sprintf('the fitted decay spans %.1f error bars (needs 5)', o.decayErrorBars);
        end
        if ~isfinite(o.T1Err)
            why{end+1} = 'T1 has no computable error bar';
        end
        if ~isfinite(o.contrastErr) || o.contrastErr >= 0.5*o.contrast
            why{end+1} = 'the contrast is not separable from zero';
        end
        for i = 1:numel(why); fprintf('      - %s\n', why{i}); end
        if o.errorBarsExceedSpread
            fprintf(['    AND THE ERROR BARS ARE BIGGER THAN THE WHOLE TRACE. The noise ceiling on\n' ...
                     '    R^2 came out NEGATIVE (%.3f), which means sum(sterr^2) exceeds the total\n' ...
                     '    variance of the data -- so the quoted uncertainty on a single point is\n' ...
                     '    larger than everything that happens across the sweep. Nothing can be\n' ...
                     '    fitted through that. Raise repeats or averages, or find the contrast\n' ...
                     '    first with a Rabi.\n'], o.r2Ceiling);
        end
    end

    fprintf('\n  T1          %10.1f +/- %.1f us\n', o.T1, o.T1Err);
    if ~isfinite(o.T1Err)
        fprintf(['              *** NO ERROR BAR. Both the Jacobian and the profile scan came back\n' ...
                 '              wider than half the search range, and an "error" that big is a\n' ...
                 '              degenerate direction rather than an uncertainty. Almost always this\n' ...
                 '              means the sweep stopped before the plateau -- see spanOverT1 below.\n']);
    end
    fprintf('  rate        %10.4g +/- %.3g per us    RATES ADD: this is the number to compare,\n', ...
            o.ratePerUs, o.rateErr);
    fprintf('              and the one t1GateTest measures a gate leak in.\n');
    if o.doubleUsed
        fprintf('  components  %.1f us (weight %.2f) and %.1f us -- T1 above is the weighted mean\n', ...
                o.T1fast, o.weight, o.T1slow);
    else
        fprintf('  stretch p   %10.3f +/- %.3f   %s\n', o.p, o.pErr, o.pMeaning);
        fprintf('  area        %10.1f us        T1*gamma(1+1/p): compare THIS when p moves\n', o.coherenceArea);
    end
    fprintf('  contrast    %10.3f +/- %.3f %%\n', 100*o.contrast, 100*o.contrastErr);
    if isfinite(o.contrastErr) && o.contrastErr >= o.contrast
        fprintf(['              *** the error is as big as the value: the decay is NOT MEASURED and\n' ...
                 '              everything above describes noise.\n']);
    end

    % ---- THE CONTROL ------------------------------------------------------
    fprintf('\n  CONTROL\n');
    if o.controlPresent && ~o.controlTestable
        fprintf(['  NOT TESTABLE. The |0> amplitude has no computable error, so there is no sigma\n' ...
                 '              to divide by -- which is itself the finding: the fit has no\n' ...
                 '              curvature in that direction. The control is neither passing nor\n' ...
                 '              failing here, and any flatness number would be fabricated.\n']);
    elseif o.controlPresent
        fprintf('  |0> flatness %9.1f sigma    the |0> trace''s amplitude over its own error\n', o.flatnessSigma);
        if o.flatnessSigma > 3
            fprintf(['              *** THE CONTROL IS MOVING. The |0> trace should not decay -- it was\n' ...
                     '              never excited. That it does means something common-mode changes\n' ...
                     '              with tau: laser power, tracking, charge state, or the readout\n' ...
                     '              windows themselves. The T1 above is measuring that as well, and\n' ...
                     '              no choice of model separates them. Fix this before quoting a T1.\n']);
        else
            fprintf('              flat, as it should be. The T1 above is a relaxation and not a drift.\n');
        end
    else
        fprintf(['              NOT AVAILABLE -- only one trace was fitted. The |0> trace is the only\n' ...
                 '              check that the decay is relaxation rather than something common-mode,\n' ...
                 '              so a T1 without it rests on trust. Run with doubleMeasurement on.\n']);
    end

    % ---- did the sweep reach the plateau ---------------------------------
    fprintf('\n  REACH\n');
    fprintf('  span/T1     %10.2f            the sweep covers this many T1s\n', o.spanOverT1);
    if o.spanOverT1 < 2
        fprintf(['              *** under 2. The plateau is an EXTRAPOLATION, so the baseline and the\n' ...
                 '              amplitude are trading against T1 along a ridge and the fit has picked\n' ...
                 '              a point on it. Extend tau to at least %.0f us.\n'], 3*o.T1);
    end
    if isfinite(o.plateauSigma)
        fprintf('  plateau     %+10.2f sigma      how far the last quarter sits from the fitted plateau\n', ...
                o.plateauSigma);
        if abs(o.plateauSigma) > 3
            fprintf(['              *** the end of the sweep does not sit on the fitted plateau. Either\n' ...
                     '              the decay has not finished or the shape is wrong -- try\n' ...
                     '              model = ''auto'' before extending the sweep.\n']);
        end
    end
    fprintf('  decay span  %10.1f error bars\n', o.decayErrorBars);
    if o.decayErrorBars < 5
        fprintf('              *** barely above the noise. A T1 from this is a number with nothing under it.\n');
    end

    % ---- the second component --------------------------------------------
    if o.doubleTest.tried
        fprintf('\n  SECOND COMPONENT\n');
        if o.doubleUsed
            fprintf('  accepted: dSSR %.1f (needs %.0f), constants differ %.2fx, weight %.1f sigma from a bound\n', ...
                    o.doubleTest.dSSR, o.doubleTest.threshold, o.doubleTest.ratio, o.doubleTest.wSigma);
            if ~isempty(o.doubleTest.failed); fprintf('  NOTE: %s\n', o.doubleTest.failed); end
        else
            fprintf('  refused: %s\n', o.doubleTest.failed);
            fprintf(['              A stretched exponential with p < 1 and a sum of two exponentials\n' ...
                     '              describe the same physics, so on a short noisy sweep they are nearly\n' ...
                     '              the same curve. Being refused here is the normal outcome and not a\n' ...
                     '              statement that there is only one rate -- read p for that.\n']);
        end
    end

    % ---- the gate dose ---------------------------------------------------
    if isfinite(o.gateOpenUs)
        fprintf('\n  GATE\n');
        fprintf('  the run was taken with the microwave switch open for %.1f us per repeat-half', o.gateOpenUs);
        if ~isempty(o.opts) && isfield(o, 'pull') && ~isempty(o.pull.mwGateMode)
            fprintf(' (mwGateMode = ''%s'')', o.pull.mwGateMode);
        end
        fprintf('.\n');
        if o.gateOpenUs > 100
            fprintf(['              Any leak through it ADDS to the rate above, and cannot be separated\n' ...
                     '              from it inside one run. If t1GateTest has not been run at this dose,\n' ...
                     '              the T1 quoted here is a lower bound on the true one.\n']);
        end
    end

    % ---- fit quality ------------------------------------------------------
    fprintf('\n  QUALITY\n');
    fprintf('  chi2 (red)  %10.2f\n', o.chi2Reduced);
    fprintf('  R^2         %10.4f   (noise ceiling %.4f)\n', o.r2, o.r2Ceiling);
    if o.nDown > 0
        fprintf('  R^2 kept    %10.4f   (ceiling %.4f)\n', o.r2Kept, o.r2CeilingKept);
    end
    fprintf('  resid acf   %+10.2f\n', o.residAutocorr);
    if abs(o.residAutocorr) > 0.4 && o.chi2Reduced > 0.05
        fprintf(['              *** STRUCTURED. The residual walks rather than scatters. On a T1 the\n' ...
                 '              usual causes in order: more than one rate (try model = ''auto''), a\n' ...
                 '              sweep that stopped before the plateau, and a control that is moving.\n']);
    end
    if o.nDown > 0
        fprintf('\n  robust      %d point(s) down-weighted, at tau = %s us\n', o.nDown, numList(o.downT, '%.0f'));
        fprintf(['              on a T1 these are easy to miss by eye: a point diluted toward no\n' ...
                 '              contrast sits AT the plateau, which is where half the sweep already is.\n']);
    end
    fprintf('\n');
end

function s = numList(v, fmt)
    if isempty(v); s = '(none)'; return; end
    parts = cell(1, numel(v));
    for i = 1:numel(v); parts{i} = sprintf(fmt, v(i)); end
    s = strjoin(parts, ', ');
end

function s = eitherText(tf, a, b)
    if tf; s = a; else; s = b; end
end

% ======================================================================
% PLOT
% ======================================================================

function drawIt(o)
    S = fitStyle();
    nExtra = double(o.opts.showResid);
    fig = S.newFigure('T1', nExtra > 0);
    if nExtra == 0
        axM = S.axes(fig);
    else
        axM = S.axes(subplot(3,1,1:2));
    end

    % LOG x by default. A T1 sweep spans two or three decades in tau, and on a linear axis the
    % whole decay is squeezed into the first tenth of the width while the plateau gets the rest.
    % The decay is the measurement; it should get the space.
    tdLo = min(o.t); if o.opts.logX; tdLo = max(min(o.t), 1e-3); end
    if o.opts.logX
        td = logspace(log10(tdLo), log10(max(o.t)), 1500);
    else
        td = linspace(min(o.t), max(o.t), 1500);
    end
    lin = [o.plateau(:), o.amp(:)];
    Md = t1CurveAt(o, td);

    % COLOUR BY ROLE, NOT BY ROW. The |1> trace is the measurement and the |0> trace is a
    % control, and the control is exactly the thing that should recede -- so the signal takes the
    % primary orange and the control the desaturated slate, whichever row each happens to be in.
    % Taking the slots in row order gave the control the primary colour and the measurement the
    % secondary one, which inverts the figure's emphasis on the one experiment where the two
    % traces are not equals.
    %
    % Orange against slate measures dE 19.7 normal and 14.7 worst-case CVD on this palette, so
    % the pair is still separable; what changes is which one the eye goes to first.
    cols = cell(1, o.traces);
    for i = 1:o.traces
        if o.traces == 2 && i == o.ctlRow
            cols{i} = S.col.s3;      % the control: desaturated, recessive
        else
            cols{i} = S.col.s1;      % the measurement
        end
    end

    H = []; L = {};
    for i = 1:o.traces
        H(end+1) = S.data(axM, o.t, o.y(i,:), o.e(i,:), cols{i}, '', true, o.downIdx); %#ok<AGROW>
        L{end+1} = o.labels{i};                                                        %#ok<AGROW>
    end
    for i = 1:o.traces
        h = S.model(axM, td, Md(i,:));
        if i == 1; H(end+1) = h; L{end+1} = 'fit'; end                                  %#ok<AGROW>
    end
    % The plateau of the SIGNAL trace, as a reading aid: the distance from the data to this line
    % is the contrast, and where the data meets it is where the spin has forgotten.
    if isfinite(o.contrastErr) && o.contrastErr < 0.5*o.contrast
        H(end+1) = S.envelope(axM, td, o.plateau(o.sigRow) + zeros(size(td)));
        L{end+1} = sprintf('plateau = %.4f', o.plateau(o.sigRow));
    end
    % T1 itself, marked. On a log axis this is a readable position rather than a squeezed one.
    if isfinite(o.T1) && o.T1 > min(o.t) && o.T1 < max(o.t)
        yl = ylim(axM);
        hk = plot(axM, [o.T1 o.T1], yl, '--', 'Color', S.col.flag, 'LineWidth', 1.0);
        H(end+1) = hk; L{end+1} = sprintf('T_1 = %.3g \\mus', o.T1);
    end
    if o.nDown > 0
        H(end+1) = S.flag(axM, o.downT, o.y(o.sigRow, o.downIdx));
        L{end+1} = sprintf('%d down-weighted', o.nDown);
    end

    if o.opts.logX; set(axM, 'XScale', 'log'); end
    xlabel(axM, '\tau  (\mus)', 'Color', S.col.ink2);
    ylabel(axM, yLabelOf(o, S), 'Color', S.col.ink2);
    keepMask = true(1, numel(o.t)); keepMask(o.downIdx) = false;
    S.focus(axM, [Md(:); reshape(o.y(:, keepMask), [], 1)]);
    if ~o.opts.logX; S.tight(axM, o.t); end
    S.headroom(axM, 0.30);
    S.legend(axM, H, L);
    S.title(axM, titleLines(o));

    if o.opts.showResid
        axR = S.axes(subplot(3,1,3));
        S.zeroLine(axR, [min(o.t) max(o.t)]);
        Rn = o.resid ./ max(o.e, eps);
        for i = 1:o.traces
            S.resid(axR, o.t, Rn(i,:), cols{i});
        end
        S.focus(axR, reshape(Rn(:, keepMask), [], 1), 0.15);
        if o.opts.logX; set(axR, 'XScale', 'log'); end
        xlabel(axR, '\tau  (\mus)', 'Color', S.col.ink2);
        ylabel(axR, 'residual / \sigma', 'Color', S.col.ink2);
        S.title(axR, sprintf('\\chi^2_\\nu = %.2f     lag-1 autocorr %+.2f%s', ...
                             o.chi2Reduced, o.residAutocorr, ...
                             eitherText(abs(o.residAutocorr) > 0.4, ' (structured)', '')));
    end
end

function P = paramVec(o)
    if o.doubleUsed
        P = [o.T1fast, o.T1slow, o.weight];
    else
        P = [o.T1, o.p];
    end
end

function s = yLabelOf(o, S)
    s = S.yLabel;
    if isfield(o.opts, 'yLabel') && ~isempty(o.opts.yLabel); s = o.opts.yLabel; end
end

function T = titleLines(o)
    if ~o.measured
        % Lead with the verdict, not with T1. The numbers still appear -- they are what the
        % optimiser returned and hiding them would be its own kind of dishonesty -- but on the
        % second line and labelled as unquotable.
        bits = {};
        if o.decayErrorBars < 5
            bits{end+1} = sprintf('decay spans %.1f error bars', o.decayErrorBars);
        end
        if o.errorBarsExceedSpread
            bits{end+1} = 'error bars exceed the whole trace';
        end
        if ~isfinite(o.T1Err); bits{end+1} = 'no error bar on T_1'; end
        L1 = ['NO MEASURABLE DECAY  --  ' strjoin(bits, ', ')];
        L2 = sprintf(['unquotable: T_1 = %.3g \\mus, p = %.2f, contrast = %.2f %% -- where the ' ...
                      'optimiser stopped, not a measurement'], o.T1, o.p, 100*o.contrast);
        b3 = {sprintf('\\chi^2_\\nu = %.2f', o.chi2Reduced), ...
              sprintf('R^2 = %.3f of %.3f', rKept(o), rCeil(o))};
        if o.controlPresent && ~o.controlTestable
            b3{end+1} = '|0\rangle control not testable';
        elseif o.controlPresent
            b3{end+1} = sprintf('|0\\rangle control %.1f\\sigma', o.flatnessSigma);
        end
        b3{end+1} = sprintf('sweep = %.2f T_1', o.spanOverT1);
        if isfinite(o.gateOpenUs); b3{end+1} = sprintf('gate %.0f \\mus', o.gateOpenUs); end
        if o.nDown > 0; b3{end+1} = sprintf('%d down-weighted', o.nDown); end
        T = {L1, L2, strjoin(b3, '     |     ')};
        return
    end
    if o.doubleUsed
        L1 = sprintf('T_1 = %.4g \\pm %.3g \\mus (weighted mean of %.3g and %.3g)     contrast = %.3f \\pm %.3f %%', ...
                     o.T1, o.T1Err, o.T1fast, o.T1slow, 100*o.contrast, 100*o.contrastErr);
    else
        L1 = sprintf('T_1 = %.4g \\pm %.3g \\mus     p = %.2f \\pm %.2f     contrast = %.3f \\pm %.3f %%', ...
                     o.T1, o.T1Err, o.p, o.pErr, 100*o.contrast, 100*o.contrastErr);
    end
    b2 = {sprintf('rate = %.3g /\\mus', o.ratePerUs)};
    if ~o.doubleUsed; b2{end+1} = sprintf('area = %.3g \\mus', o.coherenceArea); end
    b2{end+1} = sprintf('sweep = %.1f T_1', o.spanOverT1);
    b2{end+1} = sprintf('decay spans %.0f error bars', o.decayErrorBars);
    if isfinite(o.gateOpenUs); b2{end+1} = sprintf('gate %.0f \\mus', o.gateOpenUs); end

    b3 = {};
    if o.controlPresent && o.controlTestable
        b3{end+1} = sprintf('|0\\rangle control %.1f\\sigma%s', o.flatnessSigma, ...
                            eitherText(o.flatnessSigma > 3, ' -- MOVING', ' (flat)'));
    elseif o.controlPresent
        b3{end+1} = '|0\rangle control not testable';
    else
        b3{end+1} = 'no control trace';
    end
    b3{end+1} = sprintf('\\chi^2_\\nu = %.2f', o.chi2Reduced);
    b3{end+1} = sprintf('R^2 = %.3f of %.3f', rKept(o), rCeil(o));
    b3{end+1} = sprintf('resid acf %+.2f', o.residAutocorr);
    if isfinite(o.plateauSigma); b3{end+1} = sprintf('plateau %+.1f\\sigma', o.plateauSigma); end
    if o.nDown > 0; b3{end+1} = sprintf('%d down-weighted', o.nDown); end

    T = {L1, strjoin(b2, '     |     '), strjoin(b3, '     |     ')};
end

function r = rKept(o)
    r = o.r2; if o.nDown > 0; r = o.r2Kept; end
end

function r = rCeil(o)
    r = o.r2Ceiling; if o.nDown > 0; r = o.r2CeilingKept; end
end

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d); given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('%s: unknown option "%s". Valid: %s.', who, bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}); o.(known{i}) = d.(known{i}); end
    end
end
