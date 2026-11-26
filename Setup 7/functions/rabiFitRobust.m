function out = rabiFitRobust(args)
%RABIFITROBUST Robust Rabi fit. ONE STRUCT IN, ONE STRUCT OUT.
%
%   out = rabiFitRobust(struct('exp', b))
%   out = rabiFitRobust(struct('exp', b, 'nFreq', 2, 'showResid', true))
%   out = rabiFitRobust(struct('t', t, 'y', y, 'e', e))       % plain vectors
%
% The data goes in the same struct as the options -- `exp` for an experiment object, or
% `t`/`y`/`e` for vectors -- and everything comes back in `out`, with out.ok saying whether
% to believe it. The old six positional outputs are gone: the struct was the SIXTH, and
% asking for one output bound `pulses`, a 1x3 double, so `o.ok` raised "Dot indexing is not
% supported for variables of this type". That error was blamed on the data extraction for
% three rounds before anyone read the signature.
%
% The old positional values are all still there as NAMED fields, so nothing that read them
% has to change: out.pulses, out.f, out.T, out.coeff, out.coeffError.
%
%RABIFITROBUST Rabi fit with no start point to get wrong, and honest error bars.
%
%   (see the usage block at the top -- one struct in, one struct out)
%
% THE DIAGNOSTICS STRUCT IS THE SIXTH OUTPUT. `out = rabiFitRobust(b)` binds `pulses`,
% a 1x3 double, and any out.field then raises "Dot indexing is not supported for
% variables of this type". The order is fixed at rabiFit's for drop-in compatibility,
% so the struct has to come last.
%
% ONE OPTIONS STRUCT, not positional flags. Everything optional lives in opts and gets
% filled in by withDefaults, the same convention as calibFieldSurvey, calibFieldAlign
% and esrMultiLineFit. That is deliberate and it is not only taste: rabiFit's signature
% is (Exp, toSave, path, toPlot, expDecayOrder, varargin), so the common call reads
%
%       rabiFit(b, 0, 0, 1, '', 'f', 9)
%
% in which the two zeros mean different things, the empty string is a decay order, and
% the last pair is the start point that pinned the frequency at its bound. Named fields
% cannot be silently transposed, and a typo in a field name is caught by withDefaults
% rather than becoming a value for whatever argument happens to sit in that slot.
%
% NOTE this therefore does NOT accept rabiFit's positional arguments. Calls have to be
% updated -- see the mapping at the bottom of this help.
%
% SELF-CONTAINED. It calls nothing from the lab codebase -- not extract1Ddata, not
% fitStd, not findpeaks, not the Curve Fitting Toolbox. Base MATLAB only. The data
% extraction, the averaging, the error propagation, the search and the covariance are
% all in this file, so there is one place to read and one place to fix.
%
% Same first five outputs as rabiFit, same meanings and units, so it is a drop-in.
% The sixth is a diagnostics struct rather than a cfit object.
%
% ===================================================================
% WHY A NEW ONE. THE OLD FIT DID NOT FAIL BY BAD LUCK.
% ===================================================================
%
% On the 13.5 MHz / T = 0.19 us result: the data oscillates at full contrast all the
% way to 0.5 us, with peaks near 0.07, 0.14, 0.20, 0.27, 0.33, 0.40 and 0.46 us. That
% is a period of 0.0655 us -- 15.3 MHz -- and no measurable decay. Four structural
% causes, not bad luck:
%
% 1. THE FREQUENCY WAS PINNED AT A BOUND. rabiFit sets lower and upper to 0.5x and
%    1.5x the start point. The call passed 'f', 9, so the fit could not leave
%    [4.5, 13.5] MHz, and it returned 13.50000 -- exactly 1.5x9. The true value was
%    outside the box. Any fitter that reports a parameter resting on its own bound as
%    a measurement will do this every time.
%
% 2. THE AUTOMATIC GUESS IS BROKEN ANYWAY. It calls findpeaks with
%    'MinPeakDistance', 3 -- three SAMPLES, not microseconds -- and then does
%       pks = [signal(1); pks];   locs = [1; locs];
%    injecting the first data point as a peak whether or not it is one. A first point
%    near a minimum inflates the mean peak spacing and guesses the frequency low.
%
% 3. THE AMPLITUDE START IS 2x TOO BIG AND BOUNDED BELOW BY THE TRUTH. contrast is
%    max - min, peak to peak, but the model's 'a' is the half-amplitude -- hence the
%    2*a in the printed contrast. So the start point is 2x the true value and
%    lower = 0.5*contrast is exactly the true value: 'a' is pinned at its lower bound
%    even when everything else goes right.
%
% 4. THE DECAY GUESS DIVIDES BY THAT FAKE PEAK. find(pks/pks(1) < exp(-1), 1) is
%    meaningless when pks(1) - centreLine is negative, which it is whenever the first
%    point sits below the centre line.
%
% With the frequency wrong, the decay is the parameter that absorbs the phase error.
% That is the flat line through the second half of the data.
%
% ===================================================================
% THE METHOD: VARIABLE PROJECTION, SO THERE IS NO FREQUENCY GUESS
% ===================================================================
%
%   y(t) = c + m*(t - tbar) + E(t) * SUM_k [ a_k cos(2 pi f_k t) + b_k sin(2 pi f_k t) ]
%   E(t) = exp(-(t/T)^q)
%
% The model is LINEAR in c, m and every a_k, b_k once (f_k, T, q) are fixed. So those
% are never searched: at each trial (f, T, q) they are solved in closed form by
% weighted least squares, leaving a two or three parameter problem. The frequency is
% therefore found by SCANNING every frequency the sampling can resolve -- there is no
% start point to be wrong and no bound to be pinned against. In order:
%
%   (a) a weighted periodogram over f with the linear block solved exactly at each
%       point. A Lomb-Scargle in all but name, and it copes with the uneven tau grid
%       that AWGExpRabi's LoadAWG leaves behind when it snaps taus onto existing
%       residue waveforms.
%   (b) the best few minima become candidates; each is gridded over decay rate 1/T,
%       INCLUDING zero, and over q.
%   (c) Levenberg-Marquardt polish with a numerical Jacobian, weights 1/sterr^2, and
%       the covariance matrix retained.
%
% Step (a) on its own gets the 15.3 MHz dataset right.
%
% ===================================================================
% WHAT IT REPORTS THAT THE OLD ONE COULD NOT
% ===================================================================
%
% A DECAY LOWER BOUND WHEN THE DECAY IS NOT RESOLVED. If the fitted rate is not at
% least two sigma above zero, a number for T is invented information. The honest
% output is "T > x", with decayResolved = false.
%
% ERRORS ON THE PULSE TIMES. piTime = 0.5/f + x0 depends on the frequency and the
% phase, which are correlated, so its variance needs the cross term. Errors are
% Birge-scaled in both directions and floored at 0.05, the same convention as
% esrMultiLineFit and for the same measured reason: sterr on this setup is
% drift-dominated and carries the right shape with the wrong scale.
%
% A RESIDUAL AUTOCORRELATION. The number that says "one frequency is the wrong
% model". A single frequency fitted to a beating pair leaves strongly correlated
% residuals while chi2 can still look survivable. 14N hyperfine sits 2.16 MHz away,
% and a field off the <111> axis drives more than one orientation, so this is not a
% theoretical concern. opts.nFreq = 2 then fits both.
%
% ===================================================================
% OPTIONS (opts struct, all optional)
%   nFreq        1 or 2 (default 1)
%   q            decay exponent: [] searches {1, 2}, or fix it (default [])
%   linearDrift  allow a sloped baseline (default true). One parameter, and it stops
%                slow laser drift being absorbed by the decay -- over half a
%                microsecond the two look alike.
%   fMax, fMin   MHz, ends of the frequency scan (defaults: just under Nyquist, and
%                half an oscillation across the sweep)
%   nGrid        periodogram resolution (default 4000)
%   robust       one 4-sigma outlier pass AFTER a fit that is already right, then
%                refit (default true). rabiFit's version runs against a model that
%                may be wrong, which removes data rather than glitches.
%   avgOutlierMad  reject a single average's point if it sits this many MAD from the
%                median across averages, and only when there are at least 5 averages
%                (default 4). See pullData for why the existing 1-sigma cut is worse.
%   plotIt       draw the data and the fit, one axes, nothing else (default true)
%   showResid    add a residuals panel below it (default false)
%   showScan     add the frequency-scan figure (default false). Worth turning on when a
%                fit looks wrong: a single sharp minimum means the frequency is not in
%                doubt, several comparable ones mean the sampling cannot separate the
%                aliases and the tau grid needs changing.
%
% OUTPUT (out) -- f, fErr, T, TErr, Tlower, decayResolved, q, contrast, contrastErr,
%   x0, x0Err, pulses, pulsesErr, centre, slope, chi2Reduced, birge, residAutocorr,
%   pointsPerPeriod, periodsSpanned, fNyquist, periodogram, SNR, B1equiv, t, y, e,
%   model, resid, ok, why
%
% MIGRATING A rabiFit CALL
%   rabiFit(b, 0, 0, 1, '', 'f', 9)      ->  rabiFitRobust(struct('exp', b))
%   rabiFit(b, 1, path, 1)               ->  rabiFitRobust(struct('exp', b, 'save', true, 'path', path))
%   rabiFit(b, 0, 0, 1, 2)               ->  rabiFitRobust(struct('exp', b, 'q', 2))
%   rabiFit(..., 'f', 9)                 ->  no equivalent, and none needed: the
%                                            frequency is scanned, not started from

    % ---- one options struct, filled in by withDefaults --------------------
    % Two accepted call shapes and nothing else: (Exp, opts) or (t, y, e, opts). An
    % unrecognised field raises rather than being ignored, because a silently dropped
    % 'nfreq' would produce a perfectly plausible single-frequency fit.
    DEFAULTS = struct('nFreq', 1, 'q', [], 'linearDrift', true, ...
        'fMax', [], 'fMin', [], 'nGrid', 4000, 'robust', true, ...
        'avgOutlierMad', 4, 'plotIt', true, 'save', false, 'path', '', 'verbose', true, ...
        'showResid', false, 'showScan', false);

    if ~exist('args', 'var') || isempty(args); args = struct(); end
    if ~isstruct(args)
        error(['rabiFitRobust takes ONE struct: rabiFitRobust(struct(''exp'', b)) or ' ...
               'rabiFitRobust(struct(''t'', t, ''y'', y, ''e'', e)). The positional forms ' ...
               '(Exp, opts) and (t, y, e, opts) are gone, and so are the six outputs.']);
    end
    DEFAULTS.exp = []; DEFAULTS.t = []; DEFAULTS.y = []; DEFAULTS.e = [];
    DEFAULTS.yLabel = '';
    opts = args;
    reps = 1; nAvg = 1;
    opts = withDefaults(opts, DEFAULTS, 'rabiFitRobust');
    if isempty(opts.path); opts.path = pwd; end
    Exp = opts.exp;
    fromVectors = isempty(Exp);
    if fromVectors
        if isempty(opts.t) || isempty(opts.y)
            out = struct('ok', false, 'why', ['nothing to fit: give exp, or t and y. ' ...
                         'rabiFitRobust(struct(''exp'', b)).'], 'opts', opts);
            fprintf('\nrabiFitRobust FAILED: %s\n', out.why); return
        end
        t = opts.t; y = opts.y; e = opts.e;
    end
    if ~fromVectors
        try
            [t, y, e, reps, nAvg] = pullData(Exp, opts.avgOutlierMad, opts.verbose);
        catch err
            % SAY WHICH ACCESS FAILED. "Dot indexing is not supported for variables of
            % this type" names neither the variable nor the field, so on an unfamiliar
            % experiment class it is nearly useless -- it cost two runs on the setup
            % before this diagnostic existed. Walk the properties the fitter needs and
            % report the class of each, so one run identifies the problem.
            describeExp(Exp);
            rethrow(err);
        end
    end

    % ANY failure past this point reports WHERE. Three rounds of guessing at a
    % "Dot indexing is not supported for variables of this type" have cost more of your
    % time than the fit is worth: MATLAB names neither the variable nor the field, and
    % without the stack there is nothing to go on. err.stack has the file, the function
    % and the line, so one run pins it.
    try
        out = fitCore(t, y, e, opts);
    catch err
        fprintf('\n--- rabiFitRobust failed INSIDE the fit ---\n');
        fprintf('  %s\n  %s\n', err.identifier, err.message);
        for k = 1:numel(err.stack)
            fprintf('    at %s > %s, line %d\n', err.stack(k).file, ...
                    err.stack(k).name, err.stack(k).line);
        end
        fprintf('  sizes handed in: t %s, y %s, e %s\n', ...
                mat2str(size(t)), mat2str(size(y)), mat2str(size(e)));
        fprintf('  classes: t %s, y %s, e %s\n', class(t), class(y), class(e));
        rethrow(err);
    end
    out.repeats = reps; out.averages = nAvg;
    if ~out.ok
        fprintf('\nrabiFitRobust FAILED: %s\n', out.why);
        out.pulses = nan(1,3); out.f = NaN; out.T = NaN;
        out.coeff = nan(1,5); out.coeffError = nan(1,5);
        out.opts = opts; return
    end

    % Same definition as the old function's printed SNR, so the two are comparable.
    out.SNR = 2*out.amplitude / (mean(out.e)*sqrt(max(nAvg,1)*max(reps,1)));

    % The five-element coeff vectors the old positional outputs carried, kept as named
    % fields so nothing that read them has to change.
    out.coeff      = [out.amplitude,    out.x0,    out.T,    out.f,    out.centre];
    out.coeffError = [out.amplitudeErr, out.x0Err, out.TErr, out.fErr, out.centreErr];
    out.opts = opts;

    if opts.verbose; showReport(out); end
    if opts.plotIt; showPlots(out); end
    if opts.save
        stamp = datestr(now, 'yyyymmdd_HHMMSS');                          %#ok<TNOW1,DATST>
        fileName = fullfile(opts.path, ['RabiFitRobust_' stamp]);
        if opts.plotIt
            try; savefig([fileName '.fig']); saveas(gcf, [fileName '.png']); catch; end
        end
        save([fileName '.mat'], 'out');
    end
end

% ======================================================================
%  DATA
% ======================================================================

function [t, y, e, reps, nAvg] = pullData(Exp, madCut, verbose)
% Signal, error and tau off the experiment object, in TWO layers.
%
% PRIMARY: Exp.signalParam.value and .sterr, which perform() has already filled in, and
% Exp.mCurrentXAxisParam.value for the tau axis. These always exist after a run and
% cannot be got wrong.
%
% ENHANCEMENT, attempted second: rebuild the point errors from the individual averages,
% because the reduced sterr understates them. Per-average, the error is taken as the
% LARGER of two estimates -- propagated from each average's own sterr, and the observed
% scatter of the averages about their mean. Those answer different questions: photon
% statistics within an average, versus everything that drifts between them. On this
% setup the second is routinely an order of magnitude larger, which is the same drift
% that made AWGExpESR's sterr read 42x the point-to-point scatter. Taking the larger is
% the only choice that does not understate the uncertainty.
%
% If the enhancement fails for any reason the primary numbers are used and a line says
% so. An earlier version had only the enhanced path and died with "Dot indexing is not
% supported for variables of this type" -- a fitter that cannot fit because of how it
% reached for the data is worse than one with slightly optimistic error bars.
%
% NOT extract1Ddata, and not only for independence: that function applies three
% successive cuts to the ratio across averages, at 5, then 3, then ONE standard
% deviation --
%       s(abs(s - mean(s,2)) > 1*mean(std(s,1,2))) = NaN;
% A one-sigma cut discards about a third of perfectly good Gaussian data, after two
% harsher passes have already shrunk the spread it measures against, and it biases the
% survivors toward the mean so the scatter looks smaller than it is. Here nothing is cut
% unless there are at least 5 averages to judge with, and then on a median-absolute-
% deviation rule.
    reps = 1; nAvg = 1;

    % ---- tau axis --------------------------------------------------------
    t = [];
    try
        ax = Exp.mCurrentXAxisParam;
        if ~isempty(ax); t = ax.value(:); end
    catch
    end
    if isempty(t)
        try; t = Exp.tau(:); catch; end
    end

    % ---- primary: the reduced result -------------------------------------
    y = []; e = [];
    try
        y = Exp.signalParam.value(:);
        e = Exp.signalParam.sterr(:);
    catch err
        error(['rabiFitRobust: could not read Exp.signalParam (%s). Has the experiment ' ...
               'actually run?'], err.message);
    end
    if isempty(y)
        error('rabiFitRobust: Exp.signalParam.value is empty -- the experiment has not run.');
    end
    try; reps = Exp.repeats; catch; reps = 1; end
    usedEnhanced = false;

    % ---- enhancement: rebuild the errors from the averages ---------------
    try
        sig = Exp.signal; ste = Exp.sterr;
        if ~isempty(sig) && size(sig, 1) >= 2 && numel(size(sig)) == 3
            filled = squeeze(all(all(sig ~= 0, 1), 2));
            nA = find(~filled(:)', 1) - 1;
            if isempty(nA); nA = size(sig, 3); end
            nAvg = max(nA, 1);

            s1 = reshape(sig(1, :, 1:nAvg), [], nAvg);
            s2 = reshape(sig(2, :, 1:nAvg), [], nAvg);
            d1 = reshape(ste(1, :, 1:nAvg), [], nAvg);
            d2 = reshape(ste(2, :, 1:nAvg), [], nAvg);

            r  = s1 ./ s2;
            dr = abs(r) .* sqrt((d1./max(abs(s1), eps)).^2 + (d2./max(abs(s2), eps)).^2);

            bad = ~isfinite(r);
            if nAvg >= 5
                med = median(r, 2);
                mad = 1.4826 * median(abs(r - med), 2);
                mad(mad <= 0) = inf;
                bad = bad | (abs(r - med) > madCut * mad);
            end
            r(bad) = NaN; dr(bad) = NaN;

            nPt = size(r, 1);
            yE = nan(nPt, 1); eP = nan(nPt, 1); eS = nan(nPt, 1);
            for k = 1:nPt
                v = r(k, ~isnan(r(k,:)));
                d = dr(k, ~isnan(dr(k,:)));
                if isempty(v); continue; end
                yE(k) = mean(v);
                eP(k) = sqrt(sum(d.^2))/numel(v);
                if numel(v) >= 2; eS(k) = std(v)/sqrt(numel(v)); end
            end
            eBoth = max(eP, eS);
            eBoth(~isfinite(eBoth)) = eP(~isfinite(eBoth));
            if numel(yE) == numel(y) && any(isfinite(eBoth))
                cmp = eS > eP;
                scatterWon = sum(cmp(isfinite(eS) & isfinite(eP)));
                y = yE; e = eBoth; usedEnhanced = true;
                if verbose
                    fprintf(['rabiFitRobust: %d taus, %d averages, %d repeats. Error bars ' ...
                             'from the scatter\n               between averages on %d of %d ' ...
                             'points, from photon statistics on the rest.\n'], ...
                            numel(y), nAvg, reps, scatterWon, numel(y));
                end
            end
        end
    catch err
        if verbose
            fprintf(['rabiFitRobust: using the reduced signalParam errors -- could not ' ...
                     'rebuild them\n               from the individual averages (%s). The fit ' ...
                     'is unaffected; the\n               error bars may be optimistic.\n'], ...
                    err.message);
        end
    end
    if ~usedEnhanced && verbose
        fprintf('rabiFitRobust: %d taus, %d repeats, errors from signalParam.sterr.\n', ...
                numel(y), reps);
    end

    % ---- match up --------------------------------------------------------
    if isempty(t) || numel(t) ~= numel(y)
        error(['rabiFitRobust: a tau axis of %d points does not match %d data points. ' ...
               'AWGExpRabi rewrites tau inside LoadAWG, so read it from ' ...
               'mCurrentXAxisParam rather than the vector you assigned.'], numel(t), numel(y));
    end
    keep = isfinite(t) & isfinite(y);
    t = t(keep); y = y(keep); e = e(keep);
end

% ======================================================================
%  THE FIT
% ======================================================================

function o = fitCore(t, y, e, opts)
    [t, y, e] = tidy(t, y, e);
    o = struct('ok', false, 'why', '', 't', t, 'y', y, 'e', e, 'opts', opts);
    n = numel(t);
    if n < 12
        o.why = sprintf('only %d usable points; a frequency and a decay need more', n);
        return
    end
    w = 1./e.^2;
    dt = median(diff(t));
    fNyq = 1/(2*dt);
    fMax = opts.fMax; if isempty(fMax); fMax = 0.9*fNyq; end
    fMin = opts.fMin; if isempty(fMin); fMin = 0.5/(max(t) - min(t)); end

    % ---- (a) frequency scan ----------------------------------------------
    fGrid = linspace(fMin, fMax, opts.nGrid);
    sse = inf(size(fGrid));
    for k = 1:numel(fGrid)
        sse(k) = sseOf(t, y, w, fGrid(k), 0, 1, opts);
    end
    o.periodogram = struct('f', fGrid, 'sse', sse);

    isMin = [false, sse(2:end-1) < sse(1:end-2) & sse(2:end-1) < sse(3:end), false];
    cand = find(isMin);
    if isempty(cand); [~, cand] = min(sse); end
    [~, ord] = sort(sse(cand));
    fCand = fGrid(cand(ord(1:min(5, numel(cand)))));

    % ---- (b) decay rate and exponent -------------------------------------
    if isempty(opts.q); qList = [1 2]; else; qList = opts.q; end
    lamMax = 6/(max(t) - min(t));
    lamGrid = [0, logspace(log10(0.05*lamMax), log10(lamMax), 40)];
    best = struct('sse', inf);
    for f0 = fCand
        for q0 = qList
            for lam = lamGrid
                s = sseOf(t, y, w, f0, lam, q0, opts);
                if s < best.sse
                    best = struct('sse', s, 'f', f0, 'lam', lam, 'q', q0);
                end
            end
        end
    end
    if ~isfinite(best.sse)
        o.why = 'no trial model produced a finite residual'; return
    end

    % second frequency, found the same way the first was: from the residual scan
    fs = best.f;
    if opts.nFreq == 2
        r = y - modelAt(t, y, w, best.f, best.lam, best.q, opts);
        sse2 = inf(size(fGrid));
        for k = 1:numel(fGrid)
            sse2(k) = sseOf(t, r, w, fGrid(k), 0, 1, opts);
        end
        sse2(abs(fGrid - best.f) <= 0.15*best.f) = inf;
        [~, i2] = min(sse2);
        fs = [best.f, fGrid(i2)];
    end

    % ---- (c) polish -------------------------------------------------------
    % WHEN THERE IS NO DECAY, q DOES NOT EXIST. exp(-(0*t)^q) = 1 for every q, so the
    % derivative with respect to q is identically zero; and with q = 2 the derivative
    % with respect to lam also vanishes at lam = 0, since it goes as 2*lam*t^2. Both
    % parameters disappear together and the covariance is singular. So pin q = 1,
    % where d/dlam = -t is well behaved, and report q as undetermined.
    qFixed = false;
    if best.lam <= 1e-6*lamMax
        best.q = 1; qFixed = true;
    end
    nF = numel(fs);
    p0 = [fs, best.lam, best.q];
    lb = [repmat(fMin, 1, nF), 0,        0.3];
    ub = [repmat(fMax, 1, nF), 4*lamMax, 3];
    if ~isempty(opts.q); lb(end) = opts.q; ub(end) = opts.q; end
    if qFixed;            lb(end) = 1;     ub(end) = 1;     end
    opts.nFreq = nF;

    res = @(p) sqrt(w) .* (y - modelFromP(t, y, w, p, opts));
    [p, C, chi2r, birge, ok, why] = lmFit(res, p0, lb, ub, n);
    if ~ok; o.why = why; return; end

    % ---- one outlier pass, against a model that is already right ----------
    if opts.robust
        rr = (y - modelFromP(t, y, w, p, opts)) ./ e;
        cut = 4 * 1.4826 * median(abs(rr - median(rr)));
        keep = abs(rr - median(rr)) < max(cut, 3);
        if any(~keep) && sum(keep) > numel(p) + 8
            o.dropped = t(~keep);
            t = t(keep); y = y(keep); e = e(keep); w = w(keep);
            res = @(p) sqrt(w) .* (y - modelFromP(t, y, w, p, opts));
            [p2, C2, chi2r2, birge2, ok2] = lmFit(res, p, lb, ub, numel(y));
            if ok2; p = p2; C = C2; chi2r = chi2r2; birge = birge2; end
        end
    end

    % ---- assemble ---------------------------------------------------------
    [yhat, lin, X] = modelFromP(t, y, w, p, opts);
    fFit = p(1:nF); lam = p(nF+1); q = p(nF+2);

    Clin = pinv(X' * (w .* X)) * birge^2;
    if opts.linearDrift; iC = 1; iM = 2; iA = 3; else; iC = 1; iM = 0; iA = 2; end
    a1 = lin(iA); b1 = lin(iA+1);
    A = hypot(a1, b1);
    phi = mod(atan2(b1, a1) + pi, 2*pi) - pi;      % y ~ A cos(2 pi f t - phi)

    Cab = Clin(iA:iA+1, iA:iA+1);
    dA   = [a1, b1] / max(A, eps);
    dphi = [-b1, a1] / max(A^2, eps);
    AErr   = sqrt(max(dA*Cab*dA', 0));
    phiErr = sqrt(max(dphi*Cab*dphi', 0));

    f1 = fFit(1);
    fErr = sqrt(max(C(1,1), 0));
    x0 = phi/(2*pi*f1);
    x0Err = sqrt((phiErr/(2*pi*f1))^2 + (phi*fErr/(2*pi*f1^2))^2);

    kFrac = [0.25 0.5 0.75];
    pulses = kFrac/f1 + x0;
    pulsesErr = zeros(1,3);
    for k = 1:3
        dfd = -kFrac(k)/f1^2 - phi/(2*pi*f1^2);
        dpd = 1/(2*pi*f1);
        pulsesErr(k) = sqrt((dfd*fErr)^2 + (dpd*phiErr)^2);
    end

    lamErr = sqrt(max(C(nF+1, nF+1), 0));
    resolved = isfinite(lamErr) && lam > 2*lamErr && lam > 0;
    o.Tlower = 1/max(lam + 2*max(lamErr, 0), eps);

    o.ok = true;
    o.f = f1; o.fErr = fErr;
    o.decayResolved = resolved;
    if resolved
        o.T = 1/lam; o.TErr = lamErr/lam^2;
    else
        o.T = o.Tlower; o.TErr = Inf;
    end
    o.lambda = lam; o.lambdaErr = lamErr;
    o.q = q; o.qFixed = qFixed;
    o.amplitude = A;   o.amplitudeErr = AErr;
    o.contrast = 2*A;  o.contrastErr = 2*AErr;
    o.x0 = x0;         o.x0Err = x0Err;
    o.centre = lin(iC); o.centreErr = sqrt(max(Clin(iC,iC), 0));
    o.slope = 0; if iM; o.slope = lin(iM); end
    o.pulses = pulses; o.pulsesErr = pulsesErr;
    o.chi2Reduced = chi2r; o.birge = birge;
    o.linCoeffs = lin; o.fAll = fFit;
    o.model = yhat; o.resid = y - yhat;
    o.t = t; o.y = y; o.e = e;
    if nF == 2
        o.f2 = fFit(2); o.f2Err = sqrt(max(C(2,2), 0));
        o.amplitude2 = hypot(lin(iA+2), lin(iA+3));
    end

    o.pointsPerPeriod = (1/f1)/median(diff(t));
    o.periodsSpanned  = (max(t) - min(t))*f1;
    o.fNyquist = fNyq;
    rw = o.resid ./ e;
    o.residAutocorr = sum(rw(1:end-1).*rw(2:end)) / max(sum(rw.^2), eps);

    % R^2, the plain coefficient of determination, because it is the number people read
    % at a glance. It is reported alongside reduced chi2 rather than instead of it: R^2
    % says how much of the variation the model captures, chi2 says whether the residuals
    % are consistent with the error bars. A Rabi with good contrast lands at R^2 > 0.99
    % while a model that is the wrong SHAPE can still score well there -- so the title
    % carries both, and residAutocorr is the one that catches a wrong shape.
    ssRes = sum(o.resid.^2);
    ssTot = sum((y - mean(y)).^2);
    o.r2 = 1 - ssRes/max(ssTot, eps);
    wR = 1./e.^2;
    o.r2Weighted = 1 - sum(wR.*o.resid.^2) / max(sum(wR.*(y - sum(wR.*y)/sum(wR)).^2), eps);
    o.B1equiv = f1/2.8025;
end

function [X, E] = designMatrix(t, f, lam, q, opts)
    t = t(:);
    E = exp(-(max(lam, 0)*max(t, 0)).^q);
    X = ones(numel(t), 1);
    if opts.linearDrift; X = [X, t - mean(t)]; end
    for k = 1:numel(f)
        X = [X, E.*cos(2*pi*f(k)*t), E.*sin(2*pi*f(k)*t)];                  %#ok<AGROW>
    end
end

function s = sseOf(t, y, w, f, lam, q, opts)
% Weighted residual with the linear block PROFILED OUT -- solved, not searched. This
% is what makes an exhaustive frequency scan affordable.
    X = designMatrix(t, f, lam, q, opts);
    sw = sqrt(w(:));
    A = X .* sw; b = y(:) .* sw;
    M = A'*A;
    if ~all(isfinite(M(:))) || rcond(M) < 1e-14; s = inf; return; end
    r = b - A*(M \ (A'*b));
    s = sum(r.^2);
end

function [yhat, lin, X] = modelAt(t, y, w, f, lam, q, opts)
    X = designMatrix(t, f, lam, q, opts);
    sw = sqrt(w(:));
    lin = (X .* sw) \ (y(:) .* sw);
    yhat = X*lin;
end

function [yhat, lin, X] = modelFromP(t, y, w, p, opts)
    nF = opts.nFreq;
    [yhat, lin, X] = modelAt(t, y, w, p(1:nF), p(nF+1), p(nF+2), opts);
end

function [p, C, chi2r, birge, ok, why] = lmFit(resFun, p0, lb, ub, nData)
% Levenberg-Marquardt, numerical Jacobian, column-scaled, with the covariance kept.
% Hand-rolled for the same reason as esrMultiLineFit: the covariance is the output
% that matters and neither confint nor the Optimization Toolbox provides it here.
    ok = false; why = ''; chi2r = NaN; birge = 1;
    p = p0(:)';  nP = numel(p);
    C = inf(nP);
    r = resFun(p);
    if any(~isfinite(r)); why = 'the initial model is not finite'; return; end
    cost = sum(r.^2);
    lambda = 1e-3; step = zeros(nP, 1); J = zeros(numel(r), nP);
    for iter = 1:200
        for j = 1:nP
            h = max(1e-6*abs(p(j)), 1e-9);
            pp = p; pp(j) = pp(j) + h;
            rj = resFun(pp);
            if any(~isfinite(rj)); rj = r; end
            J(:, j) = (rj(:) - r(:))/h;
        end
        H = J'*J; g = J'*r(:);
        % A parameter the model does not depend on gives a ZERO Jacobian column, and
        % then the step solve fails. That is not a numerical accident, it is the data
        % saying the parameter is not measurable -- so step only in the subspace that
        % IS constrained and let the rest keep their values and infinite errors.
        dh = sqrt(diag(H))';
        act = dh > 1e-10*max(dh);
        if ~any(act); break; end
        sc = sqrt(max(diag(H), eps));
        Ds = diag(1./sc);
        improved = false;
        for tryLambda = 1:12
            Hs = Ds(act,act)*H(act,act)*Ds(act,act) + lambda*eye(sum(act));
            stepA = -Ds(act,act)*(Hs \ (Ds(act,act)*g(act)));
            step = zeros(nP, 1); step(act) = stepA;
            pTry = min(max(p + step(:)', lb), ub);
            rTry = resFun(pTry);
            if all(isfinite(rTry)) && sum(rTry.^2) < cost
                p = pTry; r = rTry; cost = sum(r.^2);
                lambda = max(lambda/3, 1e-9); improved = true; break
            end
            lambda = lambda*4;
        end
        if ~improved || norm(step) < 1e-12; break; end
    end
    dof = max(nData - nP, 1);
    chi2r = cost/dof;
    birge = max(sqrt(chi2r), 0.05);         % two-sided, floored: refusing to scale
                                            % DOWN is not conservative, it is useless
    H = J'*J;
    dh = sqrt(diag(H))';
    act = dh > 1e-10*max(dh);
    if any(act)
        Ha = H(act, act);
        if rcond(Ha + 1e-14*eye(sum(act))) < 1e-14
            why = 'the covariance is singular on the constrained subspace -- two parameters are doing one job';
            return
        end
        C(act, act) = inv(Ha) * birge^2;
    end
    ok = true;
end

function [t, y, e] = tidy(t, y, e)
    t = t(:); y = y(:);
    if isempty(e); e = ones(size(y)); else; e = e(:); end
    if numel(e) ~= numel(y); e = ones(size(y)); end
    good = isfinite(t) & isfinite(y);
    t = t(good); y = y(good); e = e(good);
    ok = isfinite(e) & e > 0;
    if any(ok); e(~ok) = median(e(ok)); else; e = ones(size(y)); end
    [t, i] = sort(t); y = y(i); e = e(i);
end

% ======================================================================
%  REPORTING
% ======================================================================

function showReport(o)
    fprintf('\n=== rabiFitRobust ===\n');
    fprintf('  Rabi frequency   %.4f +/- %.4f MHz\n', o.f, o.fErr);
    fprintf('  contrast         %.3f +/- %.3f %%\n', 100*o.contrast, 100*o.contrastErr);
    fprintf('  x0 (RabiX0)      %+.2f +/- %.2f ns\n', 1e3*o.x0, 1e3*o.x0Err);
    fprintf('  pi/2, pi, 3pi/2  [%.1f %.1f %.1f] +/- [%.1f %.1f %.1f] ns\n', ...
            1e3*o.pulses, 1e3*o.pulsesErr);
    if o.decayResolved
        fprintf('  decay T          %.3f +/- %.3f us,  exponent q = %.2f\n', o.T, o.TErr, o.q);
    else
        fprintf('  decay T          NOT RESOLVED. T > %.2f us at 2 sigma.\n', o.Tlower);
        fprintf('                   The oscillation does not measurably decay across this\n');
        fprintf('                   sweep, so any single value of T would be invented. The\n');
        fprintf('                   envelope needs a longer tau to be seen at all.\n');
        if o.qFixed
            fprintf('                   q was pinned to 1: with no decay the envelope shape is\n');
            fprintf('                   unmeasurable, so q is not a parameter here.\n');
        end
    end
    fprintf('  R^2              %.5f   (weighted %.5f)\n', o.r2, o.r2Weighted);
    fprintf('  reduced chi2     %.2f  (parameter errors scaled %.2fx)\n', o.chi2Reduced, o.birge);
    fprintf('  sampling         %.1f points/period, %.1f periods spanned, Nyquist %.1f MHz\n', ...
            o.pointsPerPeriod, o.periodsSpanned, o.fNyquist);
    if o.pointsPerPeriod < 5
        fprintf('  *** fewer than 5 points per period -- the frequency may be an alias.\n');
    end
    if o.periodsSpanned < 2
        fprintf('  *** fewer than 2 periods spanned -- a frequency is not really measured here.\n');
    end
    if o.f > 0.7*o.fNyquist
        fprintf('  *** within 30%% of Nyquist. Halve the tau step before believing this.\n');
    end
    fprintf('  resid autocorr   %+.2f  ', o.residAutocorr);
    if abs(o.residAutocorr) > 0.4
        fprintf('<- STRUCTURED, one frequency is the wrong model.\n');
        fprintf('                   Try opts.nFreq = 2. 14N hyperfine sits 2.16 MHz away, and\n');
        fprintf('                   a field off the <111> axis drives more than one orientation.\n');
    else
        fprintf('(structureless, so one frequency describes the data)\n');
    end
    if isfield(o, 'f2')
        fprintf('  second frequency %.4f +/- %.4f MHz, amplitude ratio %.2f\n', ...
                o.f2, o.f2Err, o.amplitude2/o.amplitude);
    end
    if isfield(o, 'dropped') && ~isempty(o.dropped)
        fprintf('  dropped %d point(s) at 4 sigma from the fitted curve: %s us\n', ...
                numel(o.dropped), num2str(o.dropped(:)', '%.3f '));
    end
    if isfield(o, 'SNR')
        fprintf('  SNR/sqrt(n)      %.3f\n', o.SNR);
    end
    fprintf('  B1 equivalent    %.3f G  (f/gamma; a convention factor of order 1 applies\n', o.B1equiv);
    fprintf('                   and cancels when B1 is compared across frequencies)\n');
end

function showPlots(o)
% ONE axes by default. showResid and showScan add panels, both off unless asked. The look
% comes from fitStyle, shared with esrFitRobust and ramseyFitRobust, so the three figures
% cannot drift apart -- which they had: this file drew its fit in ORANGE and Ramsey drew its
% in black, and orange is needed for a second data series the moment there is one.
    S = fitStyle();
    tf = linspace(min(o.t), max(o.t), 2000)';
    op = o.opts; op.nFreq = numel(o.fAll);
    Xf = designMatrix(tf, o.fAll, o.lambda, o.q, op);
    yf = Xf*o.linCoeffs;

    fig = S.newFigure('Rabi', o.opts.showResid);
    if o.opts.showResid
        axM = S.axes(subplot(4, 1, 1:3));
    else
        axM = S.axes(fig);
    end

    H = []; L = {};
    bad = [];
    if isfield(o, 'dropped') && ~isempty(o.dropped)
        [~, bad] = min(abs(o.t(:) - o.dropped(:)'), [], 1);   % nearest index per dropped tau
    end
    H(end+1) = S.data(axM, o.t, o.y, o.e, S.col.s1, '', true, bad);
    L{end+1} = 'measured';
    H(end+1) = S.model(axM, tf, yf);                   L{end+1} = 'fit';
    if isfield(o, 'dropped') && ~isempty(o.dropped)
        yd = interp1(o.t, o.y, o.dropped, 'nearest', 'extrap');
        H(end+1) = S.flag(axM, o.dropped, yd);
        L{end+1} = sprintf('%d dropped', numel(o.dropped));
    end

    ylabel(axM, yLabelOf(o, S), 'Color', S.col.ink2);
    if ~o.opts.showResid
        xlabel(axM, 'pulse duration \tau  (\mus)', 'Color', S.col.ink2);
    end
    keepMask = true(1, numel(o.t)); keepMask(bad) = false;
    S.focus(axM, [yf(:); reshape(o.y(keepMask), [], 1)]);
    S.tight(axM, o.t); S.headroom(axM, 0.28);
    S.legend(axM, H, L);
    S.title(axM, rabiTitle(o));

    if o.opts.showResid
        axR = S.axes(subplot(4, 1, 4));
        S.zeroLine(axR, [min(o.t) max(o.t)]);
        S.resid(axR, o.t, o.resid ./ max(o.e, eps), S.col.s1);
        xlabel(axR, 'pulse duration \tau  (\mus)', 'Color', S.col.ink2);
        ylabel(axR, 'residual / \sigma', 'Color', S.col.ink2);
        S.tight(axR, o.t);
        S.title(axR, sprintf('\\chi^2_\\nu = %.2f, lag-1 autocorr %+.2f', ...
                             o.chi2Reduced, o.residAutocorr));
    end

    if o.opts.showScan
        f2 = S.newFigure('Rabi: frequency scan');
        ax2 = S.axes(f2);
        % o.periodogram is a STRUCT with .f and .sse, not a two-column matrix. Indexing it
        % as o.periodogram(:,1) raised "out of bound 1 (dimensions are 1x1)" -- caught by the
        % stubbed-graphics smoke test rather than on the setup, which is the point of having
        % one: showScan is off by default, so this would have sat here until the first time
        % somebody turned it on to debug something else.
        plot(ax2, o.periodogram.f, o.periodogram.sse, '-', 'LineWidth', 1.4, ...
             'Color', S.col.s1);
        for k = 1:numel(o.fAll)
            plot(ax2, [o.fAll(k) o.fAll(k)], ylim(ax2), '--', 'Color', S.col.flag, ...
                 'LineWidth', 1.2);
        end
        xlabel(ax2, 'trial frequency  (MHz)', 'Color', S.col.ink2);
        ylabel(ax2, 'weighted residual', 'Color', S.col.ink2);
        S.title(ax2, sprintf('dashed = fitted %s MHz', ...
                             strjoin(cellstr(num2str(o.fAll(:), '%.3f'))', ', ')));
    end
end

function s = yLabelOf(o, S)
    s = S.yLabel;
    if isfield(o.opts, 'yLabel') && ~isempty(o.opts.yLabel); s = o.opts.yLabel; end
end

function T = rabiTitle(o)
    a = sprintf(['\\pi = %.1f \\pm %.1f ns     f = %.3f \\pm %.3f MHz     ' ...
                 'contrast = %.2f \\pm %.2f %%'], ...
                1e3*o.pulses(2), 1e3*o.pulsesErr(2), o.f, o.fErr, ...
                100*o.contrast, 100*o.contrastErr);
    ceil0 = 1 - sum(o.e.^2)/max(sum((o.y - mean(o.y)).^2), eps);
    bits = {sprintf('T = %.2f \\mus', o.T), ...
            sprintf('\\chi^2_\\nu = %.2f', o.chi2Reduced), ...
            sprintf('R^2 = %.3f of %.3f', o.r2, ceil0), ...
            sprintf('resid acf %+.2f', o.residAutocorr)};
    T = {a, strjoin(bits, '     |     ')};
end

function o = withDefaults(o, d, who)
% Fill in the defaults, and REFUSE an unknown field name.
%
% The usual version silently ignores anything it does not recognise, which turns
% struct('nfreq', 2) into a single-frequency fit that looks entirely reasonable and is
% answering a different question. A misspelled option should be an error, not a
% surprise, so this one lists the valid names when it finds a stray.
%
% A field present but EMPTY takes the default, so struct('q', []) means "search q"
% rather than "q is nothing".
    if ~exist('who', 'var'); who = 'options'; end
    if isempty(o); o = struct(); end
    known = fieldnames(d);
    given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('%s: unknown option "%s". Valid options are: %s.', ...
              who, bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || isempty(o.(known{i}))
            o.(known{i}) = d.(known{i});
        end
    end
end
