function out = ramseyFitRobust(args)
%RAMSEYFITROBUST Robust NV Ramsey fit with the nitrogen hyperfine beat built in.
%
%   out = ramseyFitRobust(struct('exp', r))
%   out = ramseyFitRobust(struct('exp', r, 'nHyperfine', 2, 'showResid', true))
%   out = ramseyFitRobust(struct('t', t, 'y', y, 'e', e))        % plain vectors
%   out = ramseyFitRobust(struct('t', t, 'y', Y, 'e', E))        % Y is 2xN: |1> and |0>
%
% ONE STRUCT IN, ONE STRUCT OUT. The data goes in the same struct as the options: `exp` for
% an experiment object, or `t`/`y`/`e` for vectors. Everything comes back in `out`, and
% out.ok says whether to believe it.
%
% That replaces six positional outputs whose order nobody could remember. The old form
% needed `[~, ~, ~, ~, out] = ramseyFitRobust(r)` and the struct was the FIFTH -- ask for one
% output and you got T2star, a double, and `out.ok` then raised "Dot indexing is not
% supported for variables of this type". Three separate debugging sessions across this
% project ended at exactly that line. A single named struct cannot be mis-ordered.
%
% Self-contained: no Curve Fitting Toolbox, no lab fitting helpers, no findpeaks, no
% extract1Ddata. The figure comes from fitStyle, shared with the other robust fits.
%
% ===================================================================
% WHY A RAMSEY IS NOT A DECAYING COSINE
% ===================================================================
%
% Fitting one damped cosine to an NV Ramsey is wrong even when it looks right, and the run
% that prompted this file shows why: it was taken at ZERO detuning and still oscillates.
%
% The nitrogen nucleus splits the NV transition. Each shot finds the nucleus in one of its
% states, and the measurement averages over all of them, so what is recorded is a SUM of
% coherences, not one:
%
%       14-N (spin 1)    three lines at  delta - A,  delta,  delta + A     A = 2.158 MHz
%       15-N (spin 1/2)  two lines at    delta - A/2, delta + A/2          A = 3.030 MHz
%
% For 14-N that sum has a closed form worth knowing by heart:
%
%       (1/3) * SUM_k cos(2*pi*f_k*t + phi)  =  (1/3) * cos(2*pi*delta*t + phi)
%                                                     * (1 + 2*cos(2*pi*A*t))
%
% a carrier at the detuning, multiplied by a 2.158 MHz beat. At delta = 0 the carrier is
% flat and the ENTIRE oscillation is the beat -- which is exactly what "zero detuning and
% yet we have some frequency" means. A single-cosine fit handed that data reports a
% detuning of 2.16 MHz that is not there, and a T2* pulled short by the beat's first null.
%
% So the beat is in the model, not fitted around. nHyperfine = 3 by default; set it to 2
% for 15-N or 1 to fall back to a plain damped cosine (useful only for deliberately large
% detunings where the triplet is unresolved).
%
% ===================================================================
% THE MODEL
% ===================================================================
%
%       y(t) = c0 + c1*t + Amp * exp(-(t/T2*)^p) * (1/nH) * SUM_k cos(2*pi*f_k*t + phi)
%
%   c0, c1     baseline and slow drift. Fitted, not removed beforehand -- subtracting a
%              trend first biases the amplitude, and on this setup the laser drifts.
%   Amp, phi   fringe amplitude and phase
%   T2*, p     stretched exponential. p is FITTED (bounds 0.5 to 3) rather than fixed at 2,
%              because which p is right is a statement about the bath and this setup has
%              never measured it. If it lands on a bound the report says so.
%   delta      detuning. |delta| only: the model is even in delta, since flipping its sign
%              leaves every cosine unchanged and only flips the sines, which phi absorbs.
%   A          hyperfine splitting. Fixed at the isotope value by default; free it with
%              fitHyperfine when you want to check the isotope rather than assume it.
%
% ===================================================================
% HOW IT IS FITTED, AND WHY NOT ALL AT ONCE
% ===================================================================
%
% Given (delta, T2*, p, A) the model is LINEAR in everything else, because
%
%       SUM_k cos(2*pi*f_k*t + phi) = cos(phi)*SUM_k cos(...) - sin(phi)*SUM_k sin(...)
%
% so c0, c1 and the two quadrature amplitudes come out of one weighted least squares with
% no starting guess and no chance of a local minimum. Amp = hypot(a,b), phi = atan2(-b,a).
%
% That leaves THREE nonlinear parameters (four with A free), which is few enough to grid
% exhaustively and then polish. The grid over delta is the important part: a periodogram
% cannot seed this model, because at zero detuning the spectral power sits at A and not at
% delta at all, so peak-picking would return 2.16 MHz as the detuning every time. Scanning
% delta and letting the profiled chi2 choose gets it right by construction.
%
% Polishing is a Nelder-Mead written out below -- no fminsearch, so the file behaves the
% same in Octave, where these are parse-checked.
%
% ===================================================================
% ROBUSTNESS, WHICH THIS DATA NEEDS
% ===================================================================
%
% Classic Ramsey runs on this setup show isolated points pulled toward ratio = 1. They are
% not noise: the classic path re-arms the AWG at EVERY tau, and repeats that fire before the
% instrument resumes see no microwave, so that point is diluted toward "no contrast". Two
% such points are visible in the 2.5 us gate-test reference at tau = 1.32 and 1.61 us.
%
% They are exactly the points least-squares cares most about. So after an initial fit the
% weights are re-derived with a Tukey bisquare (IRLS, 4.685 sigma), which leaves ordinary
% scatter alone and pulls the weight out of a point sitting many sigma off the curve. The
% report names how many were down-weighted and where -- a fit that quietly discards a fifth
% of the sweep is worse than one that does not converge.
%
% ===================================================================
% BOTH PROJECTIONS AT ONCE
% ===================================================================
%
% With doubleMeasurement on, AWGExpRamsey publishes TWO traces: signalParam is the |1>
% projection (detection rows 1-2) and signalParam2 the |0> one (rows 3-4). They are the same
% physics read out two ways, so they share delta, A, T2*, p and phi and differ only in
% baseline and in the SIGN of the amplitude.
%
% When both are present they are fitted JOINTLY on delta, A, T2* and p -- but NOT on phi and
% NOT on amplitude. Each projection keeps its own baseline, drift, amplitude and PHASE,
% because those are properties of a readout window rather than of the spin. The shared
% parameters are what the joint fit buys: it roughly halves the error on T2* and delta, and
% the ratio Amp(|0>)/Amp(|1>) becomes a free consistency check that should come back at -1.
% If it does not, the two readout windows are not doing what the class says they are, and
% that is worth knowing before the T2* is quoted. Pass useProjection to fit just one.
%
% (An earlier version of this comment claimed phi was shared. It never was -- each trace has
% had its own quadrature pair, hence its own amplitude and phase, since the first version.
% The comment was wrong, not the code.)
%
% ===================================================================
% "THE DATA LOOKS A LITTLE OUT OF PHASE WITH THE FIT"
% ===================================================================
%
% Phase IS fitted, freely, and separately for each projection -- see above. So a residual
% that still looks like a phase error is not an unfitted phase; it is something phi cannot
% absorb, and the difference matters because it points somewhere specific.
%
%   A CONSTANT phase offset gets absorbed completely. phi is linear in the model, so the fit
%   cannot leave one behind.
%
%   A FREQUENCY error cannot be absorbed. It accumulates: the model and the data drift apart
%   as tau grows, and the residual oscillates at the fringe period with an amplitude that
%   follows the envelope. On two projections it appears ANTI-CORRELATED, because the fringe
%   enters them with opposite sign while the baseline does not.
%
% That second pattern is the one to look for, and residAutocorr measures it -- see below. Its
% usual cause here is not delta, which is fitted, but A, which is not: 2.158 MHz is the AXIAL
% hyperfine value, and with the field off axis (this setup has ~6 G at an angle) the
% effective splitting shifts by a per cent or two. A 1% error in A is a 0.02 MHz frequency
% error, which slips about 0.3 rad over a 2.5 us record -- small, systematic, and exactly
% what a structured residual at the fringe period looks like.
%
% So autoHyperfine is ON by default: the fit is done with A fixed, then again with A free,
% and the free one is kept only if it improves the fit by more than 3 sigma for one extra
% parameter. Both are reported, so you can see whether it mattered.
%
% If a structured residual SURVIVES a free A, the next thing to try is unequal hyperfine
% line weights: optical pumping polarises the nitrogen, so the three lines are not exactly
% 1/3 each, and an asymmetric triplet has its nodes in slightly different places than a
% symmetric one. That is a real effect and not currently modelled -- it is the honest next
% step rather than something to bolt on speculatively.
%
% ===================================================================
% OPTIONS
% ===================================================================
%
%   exp                the experiment object. Or give t/y/e instead.
%   t, y, e            plain vectors. y and e may be 2xN for the two projections.
%   nHyperfine (3)     3 for 14-N, 2 for 15-N, 1 for a plain damped cosine
%   hyperfineMHz ([])  [] uses 2.158 for nH=3, 3.030 for nH=2. Ignored when nH=1.
%   fitHyperfine (false)  free A as a fourth nonlinear parameter, unconditionally
%   autoHyperfine (true)  try A free as well and keep it only if it wins by 3 sigma. This is
%                      the answer to a residual that looks phase-shifted; see above.
%   useProjection ('both')  'both', 'one' (|1>), or 'zero' (|0>)
%   robust (true)      IRLS with a Tukey bisquare after the first pass
%   tMin ([])          drop taus below this. [] takes 2*halfPiTime from the experiment, so
%                      the points where the two pi/2 pulses still overlap are not fitted as
%                      free evolution. Set 0 to keep everything.
%   deltaMaxMHz ([])   top of the detuning grid. [] uses 0.45/dt, safely under Nyquist.
%   nGridDelta (300), nGridT2 (16)
%   avgOutlierMad (4)  drop an AVERAGE whose ratio is this many MADs from the median, when
%                      there are at least 5 averages to judge with. Inf keeps everything.
%   plotIt (true), verbose (true)
%   showResid (false)    add a residuals-over-sigma panel
%   showSpectrum (false) add the fringe spectrum with the fitted components marked. Both are
%                      OFF by default: the main axes answers the question, and two extra
%                      panels shrink it every time. Turn them on when the fit looks wrong --
%                      showResid is where a structured residual is visible.
%   yLabel               override the y axis label
%
% ===================================================================
% OUT
% ===================================================================
%
%   ok, why                did it work, and if not why
%   T2star, T2starErr      microseconds
%   delta, deltaErr        detuning, MHz. |delta| only -- the model is even in delta.
%   contrast, contrastErr  peak-to-peak fringe over the baseline, extrapolated to t = 0
%   p, pErr                stretched-exponential exponent
%   A, AErr                hyperfine splitting used, and its error if it was fitted
%   freqs                  the component frequencies in MHz, ascending
%   residAutocorr          lag-1 autocorrelation of the residual. |.| > 0.4 means the model
%                          is missing something -- see the phase note above
%   hyperfineFreed         true if autoHyperfine decided a free A was worth it
%   phi, amp, baseline, slope     per projection
%   ampRatio               Amp(|0>)/Amp(|1>), should be -1
%   chi2Reduced, r2, r2Ceiling, r2Kept, r2CeilingKept
%   t, y, e, model, resid, weights, nDown, downIdx, downTau, traces, labels
%   spanUs, beatsSpanned, pointsPerBeat, nyquistMHz, tMin, nDropped
%   coeffError             [T2starErr, deltaErr, contrastErr, pErr], kept so old call sites
%                          that indexed it still work

    % ================================================================
    % ONE STRUCT IN
    % ================================================================
    DEFAULTS = struct('exp', [], 't', [], 'y', [], 'e', [], ...
        'nHyperfine', 3, 'hyperfineMHz', [], 'fitHyperfine', false, 'autoHyperfine', true, ...
        'useProjection', 'both', 'robust', true, 'tMin', [], 'deltaMaxMHz', [], ...
        'nGridDelta', 300, 'nGridT2', 16, 'avgOutlierMad', 4, ...
        'plotIt', true, 'verbose', true, 'showResid', false, 'showSpectrum', false, ...
        'yLabel', '');

    if ~exist('args', 'var') || isempty(args); args = struct(); end
    if ~isstruct(args)
        error(['ramseyFitRobust takes ONE struct: ramseyFitRobust(struct(''exp'', r)) or ' ...
               'ramseyFitRobust(struct(''t'', t, ''y'', y, ''e'', e)). The positional form ' ...
               '(r, opts) and the six-output form are gone.']);
    end
    opts = withDefaults(args, DEFAULTS, 'ramseyFitRobust');

    out = struct('ok', false, 'why', '', 'opts', opts);
    Exp = opts.exp; halfPi = [];
    fromVectors = isempty(Exp);
    if fromVectors
        if isempty(opts.t) || isempty(opts.y)
            out.why = ['nothing to fit: give exp, or t and y. ' ...
                       'ramseyFitRobust(struct(''exp'', r)).'];
            fprintf('\nramseyFitRobust FAILED: %s\n', out.why); return
        end
        t = opts.t; Y = opts.y; E = opts.e;
    end

    % ================================================================
    % THE DATA
    % ================================================================
    if ~fromVectors
        try
            [t, Y, E, halfPi, info] = pullRamsey(Exp, opts);
            out.pull = info;
        catch err
            describeExp(Exp);
            rethrow(err);
        end
    end
    pullInfo = struct(); if isfield(out, 'pull'); pullInfo = out.pull; end
    t = t(:)';
    if isvector(Y); Y = Y(:)'; end
    if isvector(E); E = E(:)'; end
    if size(Y, 2) ~= numel(t)
        if size(Y, 1) == numel(t); Y = Y.'; E = E.'; else
            out.why = sprintf('%d taus against a %dx%d signal', numel(t), size(Y,1), size(Y,2));
            fprintf('\nramseyFitRobust FAILED: %s\n', out.why); return
        end
    end
    if isempty(E) || all(~isfinite(E(:))); E = ones(size(Y)); end
    if ~isequal(size(E), size(Y)); E = repmat(median(E(isfinite(E) & E>0)), size(Y)); end

    % which projections to fit
    [Y, E, labels] = pickProjections(Y, E, opts.useProjection);
    nT = size(Y, 1);

    % ---- drop the pulse-overlap region -----------------------------------
    tMin = opts.tMin;
    if isempty(tMin)
        if ~isempty(halfPi) && isfinite(halfPi) && halfPi > 0
            tMin = 2*halfPi;
        else
            tMin = 0;
        end
    end
    keep = isfinite(t) & (t >= tMin) & all(isfinite(Y), 1);
    nDrop = sum(~keep);
    if nDrop > 0 && opts.verbose
        fprintf(['  dropped %d point(s) below tau = %.4f us: there the two pi/2 pulses still\n' ...
                 '  overlap, so those points are not free evolution and do not belong to this\n' ...
                 '  model. Set tMin = 0 to keep them.\n'], nDrop, tMin);
    end
    t = t(keep); Y = Y(:, keep); E = E(:, keep);
    n = numel(t);
    if n < 12
        out.why = sprintf(['only %d usable points. A Ramsey needs enough of them to separate ' ...
                           'a decay from a beat -- 20 or more, spanning at least two beats ' ...
                           '(%.2f us each).'], n, 1/hyperfineOf(opts));
        fprintf('\nramseyFitRobust FAILED: %s\n', out.why); return
    end
    E(~isfinite(E) | E <= 0) = median(E(isfinite(E) & E > 0));

    A0 = hyperfineOf(opts);

    % ---- aliasing check, before anything is fitted ------------------------
    dt = median(diff(t));
    fNyq = 1/(2*dt);
    out.nyquistMHz = fNyq;
    if opts.nHyperfine > 1 && A0 > 0.9*fNyq
        fprintf(['\n*** the hyperfine splitting (%.3f MHz) is above 90%% of Nyquist for this\n' ...
                 '    tau step (%.4f us -> %.2f MHz). The beat is aliased and NOTHING below is\n' ...
                 '    trustworthy. Use a step under %.4f us.\n'], A0, dt, fNyq, 1/(6*A0));
    end

    % ================================================================
    % THE FIT
    % ================================================================
    w = 1 ./ E.^2;                        % weights, updated by IRLS below
    dMax = opts.deltaMaxMHz;
    if isempty(dMax); dMax = 0.45/dt; end
    dGrid = linspace(0, dMax, max(opts.nGridDelta, 50));
    span = max(t) - min(t);
    T2Grid = logspace(log10(0.12*span), log10(6*span), max(opts.nGridT2, 8));

    best = gridSearch(t, Y, w, dGrid, T2Grid, 2.0, A0, opts.nHyperfine);
    P0 = [best.delta, best.T2, 2.0];
    if opts.fitHyperfine; P0(4) = A0; end
    [P, ok] = nelderMead(@(P) chiOf(P, t, Y, w, A0, opts.nHyperfine), P0, ...
                         bounds(dMax, span, A0, opts.fitHyperfine));
    if ~ok && opts.verbose
        fprintf('  (the polish did not move far from the grid optimum -- using the grid value)\n');
    end

    % ---- robust re-weighting, then refit ---------------------------------
    nDown = 0; downIdx = [];
    if opts.robust
        [~, lin] = chiOf(P, t, Y, w, A0, opts.nHyperfine);
        R = residOf(P, t, Y, lin, A0, opts.nHyperfine) ./ E;
        s = 1.4826 * median(abs(R(:) - median(R(:))));
        if s > 0
            u = R / (4.685 * s);
            bw = (abs(u) < 1) .* (1 - u.^2).^2;
            % bw < 0.25 is |residual| > 3.3 sigma. The first version used 0.5, which is
            % 2.5 sigma -- and in 162 residuals you expect a couple of those by chance, so it
            % reported 8 "outliers" on data with 2 injected. The WEIGHTING is unchanged and
            % still tapers smoothly from 0; this threshold only decides what gets called an
            % outlier in the report and ringed on the plot.
            downIdx = find(any(bw < 0.25, 1));
            nDown = numel(downIdx);
            w = (1 ./ E.^2) .* max(bw, 1e-3);
            best2 = gridSearch(t, Y, w, dGrid, T2Grid, P(3), A0, opts.nHyperfine);
            P2 = P; P2(1) = best2.delta;
            if chiOf(P2, t, Y, w, A0, opts.nHyperfine) < chiOf(P, t, Y, w, A0, opts.nHyperfine)
                P = P2;
            end
            P = nelderMead(@(Q) chiOf(Q, t, Y, w, A0, opts.nHyperfine), P, ...
                           bounds(dMax, span, A0, opts.fitHyperfine));
        end
    end

    % ---- is a FREE hyperfine splitting worth an extra parameter? ----------
    % The answer to a residual that still looks phase-shifted after phi has been fitted. See
    % the phase note in the header: 2.158 MHz is the AXIAL value, and off-axis field shifts
    % the effective splitting by a per cent or two, which is a 0.02 MHz frequency error --
    % small, systematic, and invisible to phi because it accumulates with tau.
    %
    % Gated on 3 sigma for one extra parameter (delta-SSR of 9), so it cannot quietly absorb
    % noise. Both fits are reported.
    o = struct();
    o.hyperfineFreed = false; o.chi2Fixed = chi2Of(P, t, Y, w, A0, opts.nHyperfine);
    o.AFixed = A0;
    if opts.autoHyperfine && ~opts.fitHyperfine && opts.nHyperfine > 1
        [~, ~, ssrFix] = chiOf(P, t, Y, w, A0, opts.nHyperfine);
        Bf = bounds(dMax, span, A0, true);
        Pf = nelderMead(@(Q) chiOf(Q, t, Y, w, A0, opts.nHyperfine), [P, A0], Bf);
        [~, ~, ssrFree] = chiOf(Pf, t, Y, w, A0, opts.nHyperfine);
        if ssrFix - ssrFree > 9
            P = Pf; opts.fitHyperfine = true; o.hyperfineFreed = true;
            if opts.verbose
                fprintf(['  freed the hyperfine splitting: %.3f -> %.3f MHz drops the weighted\n' ...
                         '  SSR by %.1f (3 sigma is 9), so it is a real improvement and not\n' ...
                         '  absorbed noise.\n'], A0, Pf(4), ssrFix - ssrFree);
            end
        end
    end

    [chi2, lin] = chiOf(P, t, Y, w, A0, opts.nHyperfine);
    model = modelOf(P, t, lin, A0, opts.nHyperfine);
    resid = Y - model;

    % ================================================================
    % PARAMETERS, ERRORS, AND WHAT THEY MEAN
    % ================================================================
    Afit = A0; if numel(P) >= 4; Afit = P(4); end
    o.delta = abs(P(1)); o.T2star = P(2); o.p = P(3); o.A = Afit;
    o.amp = zeros(1, nT); o.phi = zeros(1, nT); o.baseline = zeros(1, nT);
    o.slope = zeros(1, nT);
    for i = 1:nT
        o.baseline(i) = lin(i, 1);
        o.slope(i)    = lin(i, 2);
        o.amp(i)      = hypot(lin(i, 3), lin(i, 4));
        o.phi(i)      = atan2(-lin(i, 4), lin(i, 3));
    end
    % Contrast: peak-to-peak of the oscillating term extrapolated to t = 0, over the
    % baseline. Reported from the FIRST fitted projection, and from its magnitude, so the
    % sign convention of the readout does not change it.
    o.contrast = 2*abs(o.amp(1)) / max(abs(o.baseline(1)), eps);

    % errors: numerical Jacobian over the FULL parameter set, Birge scaled
    [err, cov] = paramErrors(P, lin, t, Y, w, A0, opts.nHyperfine, opts.fitHyperfine, chi2);

    % PROFILE ERRORS, as a check on the Jacobian and a replacement where it has nothing to
    % say. Scanning one parameter while re-solving the linear part and finding where the
    % weighted SSR rises by 1 is valid at a bound and in a stationary direction, which is
    % exactly where linear propagation is not. The larger of the two is taken: the Jacobian
    % accounts for correlations the scan holds fixed, the scan survives degeneracy, and
    % neither is safe alone.
    B = bounds(dMax, span, A0, opts.fitHyperfine);
    prof = struct();
    prof.delta = profileError(1, P, t, Y, w, A0, opts.nHyperfine, B);
    prof.T2    = profileError(2, P, t, Y, w, A0, opts.nHyperfine, B);
    prof.p     = profileError(3, P, t, Y, w, A0, opts.nHyperfine, B);
    o.deltaErr = pickErr(err.delta, prof.delta, dMax);
    o.T2starErr = pickErr(err.T2,   prof.T2,   20*span);
    o.pErr      = pickErr(err.p,    prof.p,    3);
    o.AErr = err.A;
    o.contrastErr = err.contrast;
    o.profileErr = prof;
    o.cov = cov;

    fk = componentFreqs(o.delta, Afit, opts.nHyperfine);

    o.chi2Reduced = chi2;
    resAll = resid(:); yAll = Y(:); eAll = E(:);
    o.r2 = 1 - sum(resAll.^2) / max(sum((yAll - mean(yAll)).^2), eps);
    o.r2Ceiling = 1 - sum(eAll.^2) / max(sum((yAll - mean(yAll)).^2), eps);
    % R^2 over the points the robust fit KEPT. Two gross outliers drag the all-points R^2 to
    % 0.836 against a ceiling of 0.986 on data the fit reproduces perfectly -- which reads as
    % a failed fit and is not one. Both are reported: r2 says what the raw data does, r2Kept
    % says whether the model describes the part of it that is real.
    o.r2Kept = o.r2;
    if nDown > 0
        m = true(1, numel(t)); m(downIdx) = false;
        rk = resid(:, m); yk = Y(:, m); ek = E(:, m);
        o.r2Kept = 1 - sum(rk(:).^2) / max(sum((yk(:) - mean(yk(:))).^2), eps);
        o.r2CeilingKept = 1 - sum(ek(:).^2) / max(sum((yk(:) - mean(yk(:))).^2), eps);
    else
        o.r2CeilingKept = o.r2Ceiling;
    end
    o.t = t; o.y = Y; o.e = E; o.model = model; o.resid = resid;
    o.weights = w; o.nDown = nDown; o.downIdx = downIdx; o.downTau = t(downIdx);
    o.labels = {labels{:}};                                            %#ok<CCAT1>
    o.traces = nT; o.freqs = fk; o.nHyperfine = opts.nHyperfine;
    o.nyquistMHz = fNyq; o.tMin = tMin; o.nDropped = nDrop;
    o.spanUs = span; o.beatsSpanned = span * Afit;
    o.pointsPerBeat = 1/(dt*max(Afit, eps));
    o.ampRatio = NaN;
    if nT == 2
        % Should be -1: the |0> projection is the complement of the |1> one. Anything else
        % says the two readout windows are not what the class documents.
        o.ampRatio = signedAmp(lin(2,:)) / signedAmp(lin(1,:));
    end
    % LAG-1 AUTOCORRELATION OF THE RESIDUAL. This is the number that turns "the data looks a
    % little out of phase with the fit" into a measurement. Independent noise gives ~0; a
    % residual that walks smoothly through the sweep -- which is what an unmodelled frequency
    % error produces -- gives a large positive value. |.| > 0.4 means the model is missing
    % something real and the parameter errors are optimistic whatever they say.
    o.residAutocorr = lag1(resid ./ max(E, eps));
    o.ok = true; o.why = '';

    o.coeffError = [o.T2starErr, o.deltaErr, o.contrastErr, o.pErr];
    out = o; out.opts = opts;

    if opts.verbose;  report(out);  end
    if opts.plotIt;   drawIt(out);  end
end

% ======================================================================
% MODEL PIECES
% ======================================================================

function r = lag1(R)
% Lag-1 autocorrelation, pooled across projections. Each row is detrended by its own mean so
% a per-trace offset does not masquerade as correlation.
    r = 0; num = 0; den = 0;
    for i = 1:size(R, 1)
        v = R(i,:) - mean(R(i,:));
        if numel(v) < 3; continue; end
        num = num + sum(v(1:end-1) .* v(2:end));
        den = den + sum(v.^2);
    end
    if den > 0; r = num/den; end
end

function c = chi2Of(P, t, Y, w, A0, nH)
    c = chiOf(P, t, Y, w, A0, nH);
end

function A = hyperfineOf(opts)
    if ~isempty(opts.hyperfineMHz); A = opts.hyperfineMHz; return; end
    switch opts.nHyperfine
        case 3; A = 2.158;      % 14-N, spin 1
        case 2; A = 3.030;      % 15-N, spin 1/2
        otherwise; A = 0;       % single line
    end
end

function fk = componentFreqs(delta, A, nH)
% The component frequencies, as MEASURED -- i.e. what a spectrum of the fringe shows. Each
% is |delta + m*A| because the fit cannot tell a frequency from its negative.
    switch nH
        case 3; m = [-1 0 1];
        case 2; m = [-0.5 0.5];
        otherwise; m = 0;
    end
    fk = sort(abs(delta + m*A));
end

function [Cs, Sn] = quadratures(t, delta, A, nH)
% The two quadratures of the hyperfine sum, so that phi enters LINEARLY. Averaged over the
% components, so the amplitude that comes out is the fringe amplitude and not nH times it.
    switch nH
        case 3; m = [-1 0 1];
        case 2; m = [-0.5 0.5];
        otherwise; m = 0;
    end
    Cs = zeros(size(t)); Sn = zeros(size(t));
    for k = 1:numel(m)
        th = 2*pi*(delta + m(k)*A)*t;
        Cs = Cs + cos(th);
        Sn = Sn + sin(th);
    end
    Cs = Cs/numel(m); Sn = Sn/numel(m);
end

function X = designAt(t, delta, T2, p, A, nH)
% [1, t, decay*cos-sum, decay*sin-sum]. Everything the model needs that is linear.
    D = exp(-(max(t, 0)/max(T2, eps)).^max(p, 0.05));
    [Cs, Sn] = quadratures(t, delta, A, nH);
    X = [ones(numel(t),1), t(:), (D(:).*Cs(:)), (D(:).*Sn(:))];
end

function [chi2, lin, ssr, ndofOut] = chiOf(P, t, Y, w, A0, nH)
% Profiled chi2: solve the linear part exactly for each trace, then add up. This is what
% makes the search only three-dimensional.
    delta = P(1); T2 = P(2); p = P(3);
    A = A0; if numel(P) >= 4; A = P(4); end
    X = designAt(t, delta, T2, p, A, nH);
    nT = size(Y, 1);
    lin = zeros(nT, 4);
    s = 0; ndof = 0;
    for i = 1:nT
        wi = w(i, :)';
        [th, nUsed] = wls(X, Y(i,:)', wi);
        lin(i,:) = th';
        r = Y(i,:)' - X*th;
        s = s + sum(wi .* r.^2);
        ndof = ndof + numel(t) - nUsed;
    end
    chi2 = s / max(ndof, 1);
    ssr = s; ndofOut = max(ndof, 1);
end

function M = modelOf(P, t, lin, A0, nH)
    delta = P(1); T2 = P(2); p = P(3);
    A = A0; if numel(P) >= 4; A = P(4); end
    X = designAt(t, delta, T2, p, A, nH);
    M = zeros(size(lin,1), numel(t));
    for i = 1:size(lin,1)
        M(i,:) = (X * lin(i,:)')';
    end
end

function R = residOf(P, t, Y, lin, A0, nH)
    R = Y - modelOf(P, t, lin, A0, nH);
end

function a = signedAmp(linRow)
% The amplitude with its sign, taken from the cosine quadrature. Used only for the |0> to
% |1> ratio, where the SIGN is the whole point.
    a = hypot(linRow(3), linRow(4)) * sign(linRow(3) + (linRow(3) == 0));
end

% ======================================================================
% SEARCH
% ======================================================================

function best = gridSearch(t, Y, w, dGrid, T2Grid, p, A0, nH)
%GRIDSEARCH Seed (delta, T2) without a starting guess, in two cheap stages.
%
% The delta grid is the point of this function. A periodogram cannot seed delta: at zero
% detuning ALL the spectral power sits at the hyperfine splitting and none at delta, so
% peak-picking returns 2.16 MHz as the detuning on every on-resonance run -- the exact error
% this whole model exists to avoid. Scanning delta and letting the profiled chi2 choose gets
% it right by construction.
%
% Two stages rather than a full outer product, because the full grid is 20x the work for no
% extra reach. delta and T2 are nearly uncorrelated here -- one sets where the oscillation
% is, the other how fast the envelope dies -- so scanning delta at a FEW spread T2 values
% finds the right delta, and T2 is then scanned once at that delta. Measured: 7.0 s to 0.4 s
% on an 81-point two-projection fit, same answer.
%
% The trig is hoisted out of the inner loop as well: the quadratures depend only on delta and
% the decay only on T2, so each is computed once per value instead of once per pair.
    best = struct('chi2', Inf, 'delta', dGrid(1), 'T2', T2Grid(1));

    % ---- stage 1: delta, at a few spread T2 ------------------------------
    probe = T2Grid(unique(max(1, round([0.2 0.5 0.85]*numel(T2Grid)))));
    for iD = 1:numel(dGrid)
        [Cs, Sn] = quadratures(t, dGrid(iD), A0, nH);
        for iT = 1:numel(probe)
            c = chiFromParts(t, Y, w, Cs, Sn, probe(iT), p);
            if c < best.chi2
                best.chi2 = c; best.delta = dGrid(iD); best.T2 = probe(iT);
            end
        end
    end

    % ---- stage 2: T2, at that delta --------------------------------------
    [Cs, Sn] = quadratures(t, best.delta, A0, nH);
    for iT = 1:numel(T2Grid)
        c = chiFromParts(t, Y, w, Cs, Sn, T2Grid(iT), p);
        if c < best.chi2
            best.chi2 = c; best.T2 = T2Grid(iT);
        end
    end
end

function chi2 = chiFromParts(t, Y, w, Cs, Sn, T2, p)
% Profiled chi2 from precomputed quadratures. Same arithmetic as chiOf, minus the trig.
    D = exp(-(max(t, 0)/max(T2, eps)).^max(p, 0.05));
    X = [ones(numel(t),1), t(:), (D(:).*Cs(:)), (D(:).*Sn(:))];
    chi2 = profileChi(X, Y, w, numel(t));
end

function chi2 = profileChi(X, Y, w, n)
    s = 0; ndof = 0;
    for i = 1:size(Y, 1)
        wi = w(i, :)';
        [th, nUsed] = wls(X, Y(i,:)', wi);
        r = Y(i,:)' - X*th;
        s = s + sum(wi .* r.^2);
        ndof = ndof + n - nUsed;
    end
    chi2 = s / max(ndof, 1);
end

function [th, nUsed] = wls(X, y, wi)
%WLS Weighted least squares that survives an IDENTICALLY ZERO column.
%
% Which happens, and not rarely: at exactly zero detuning the sine quadrature of a symmetric
% hyperfine multiplet vanishes point by point.
%
%       (1/3)[ sin(-2*pi*A*t) + sin(0) + sin(+2*pi*A*t) ]  ==  0
%
% The multiplet is symmetric about delta, so at delta = 0 it is symmetric about zero and the
% sines cancel exactly. phi is then genuinely undetermined -- there is no phase to measure on
% a purely even fringe -- and the column carries no information.
%
% Left alone, the normal equations are singular in that direction and the solve returns an
% enormous coefficient that is pure numerical noise. The first version of this file did
% exactly that and reported a contrast of 12462%, on ON-RESONANCE data, which is the case the
% whole model was written for. rcond did not catch it: with one null column among four the
% matrix is still nominally invertible to working precision.
%
% So columns whose weighted norm is negligible against the largest are DROPPED and their
% coefficients set to zero, and nUsed reports how many were really fitted so the degrees of
% freedom stay honest.
    sc = sqrt(sum((X.^2) .* wi, 1));
    keep = sc > 1e-9 * max(sc);
    th = zeros(size(X, 2), 1);
    nUsed = sum(keep);
    if nUsed == 0; return; end
    Xk = X(:, keep);
    M = Xk' * (Xk .* wi);
    b = Xk' * (y .* wi);
    if rcond(M) < 1e-14; return; end
    th(keep) = M \ b;
end

function B = bounds(dMax, span, A0, fitA)
    B = [0, dMax; 0.02*span, 20*span; 0.5, 3];
    if fitA; B(4,:) = [0.5*A0, 1.5*A0]; end
end

function [P, moved] = nelderMead(fun, P0, B)
% A compact Nelder-Mead with box constraints by clipping. Written out rather than calling
% fminsearch so this file behaves identically in Octave, where it is parse-checked, and so
% there is no question about which toolbox anything came from.
    n = numel(P0);
    P0 = clip(P0, B);
    S = repmat(P0, n+1, 1);
    for i = 1:n
        step = 0.08 * max(abs(P0(i)), 0.05*(B(i,2) - B(i,1)));
        S(i+1, i) = clipOne(P0(i) + step, B(i,:));
    end
    F = zeros(n+1, 1);
    for i = 1:n+1; F(i) = fun(S(i,:)); end
    f0 = min(F);

    for iter = 1:400
        [F, ord] = sort(F); S = S(ord, :);
        if abs(F(end) - F(1)) <= 1e-10*(abs(F(1)) + 1e-12); break; end
        m = mean(S(1:n, :), 1);
        xr = clip(m + (m - S(end,:)), B);        fr = fun(xr);
        if fr < F(1)
            xe = clip(m + 2*(m - S(end,:)), B);  fe = fun(xe);
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
    moved = (min(F) < f0 - 1e-12);
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

function [err, C] = paramErrors(P, lin, t, Y, w, A0, nH, fitA, chi2)
% Numerical Jacobian over the FULL parameter vector -- the nonlinear ones and every linear
% one -- then C = inv(J'WJ), Birge scaled by the reduced chi2 when that exceeds 1.
%
% Profiling the linear part is a fitting convenience, not a licence to ignore it in the
% covariance: T2* is strongly correlated with the amplitude and with the baseline, and an
% error bar computed over the nonlinear parameters alone comes out far too small.
    err = struct('delta', NaN, 'T2', NaN, 'p', NaN, 'A', NaN, 'contrast', NaN);
    C = [];
    nT = size(Y, 1); nP = numel(P);
    theta = [P(:); reshape(lin', [], 1)];
    nTot = numel(theta);

    f = @(th) flatModel(th, nP, nT, t, A0, nH);
    y0 = f(theta);
    J = zeros(numel(y0), nTot);
    for k = 1:nTot
        h = 1e-6 * max(abs(theta(k)), 1e-6);
        tp = theta; tp(k) = tp(k) + h;
        J(:, k) = (f(tp) - y0) / h;
    end
    wv = reshape(w', [], 1);
    % DROP DEGENERATE DIRECTIONS, for the same reason wls does. Two of them occur here and
    % both are real rather than numerical:
    %   * at delta = 0 the model is STATIONARY in delta -- it is even in delta, so the
    %     derivative vanishes exactly at the minimum. The column is zero.
    %   * the sine quadrature column is identically zero there too.
    % Inverting through those returned deltaErr = 94 MHz and T2*Err = 51 us on data whose
    % parameters were recovered to four figures. Anything dropped here gets a PROFILE error
    % instead, below, which is the right tool for a parameter at a bound.
    sc = sqrt(sum((J.^2) .* wv, 1));
    good = sc > 1e-9 * max(sc);
    s = nan(nTot, 1);
    if any(good)
        Jk = J(:, good);
        M = Jk' * (Jk .* wv);
        if rcond(M) > 1e-14
            Ck = inv(M) * max(chi2, 1);
            C = nan(nTot); C(good, good) = Ck;
            s(good) = sqrt(abs(diag(Ck)));
        end
    end
    err.delta = s(1); err.T2 = s(2); err.p = s(3);
    if fitA && nP >= 4; err.A = s(4); end

    % contrast = 2*hypot(a,b)/|c0| for the first trace: propagate through that.
    i0 = nP;                       % linear block for trace 1 starts at i0+1
    c0 = lin(1,1); a = lin(1,3); b = lin(1,4);
    R = hypot(a, b);
    if R > 0 && abs(c0) > 0
        g = zeros(1, nTot);
        g(i0+1) = -2*R/c0^2;
        g(i0+3) =  2*(a/R)/abs(c0);
        g(i0+4) =  2*(b/R)/abs(c0);
        err.contrast = sqrt(abs(g * C * g'));
    end
end

function e = pickErr(jac, prof, range)
% Take the larger, and fall back to whichever exists. An error bigger than the whole search
% range is not an error bar, it is a degenerate direction, so it is discarded.
    cand = [];
    if isfinite(jac)  && jac  > 0 && jac  < 0.5*range; cand(end+1) = jac;  end
    if isfinite(prof) && prof > 0 && prof < 0.5*range; cand(end+1) = prof; end
    if isempty(cand); e = NaN; else; e = max(cand); end
end

function e = profileError(k, P, t, Y, w, A0, nH, B)
%PROFILEERROR 1-sigma from scanning one parameter: where the weighted SSR rises by 1.
%
% The linear part is re-solved at every step, so this is a genuine profile over everything
% that is linear -- baseline, drift, amplitude and phase -- and only the OTHER nonlinear
% parameters are held. That makes it a lower bound where those are correlated, which is why
% the caller takes the larger of this and the Jacobian rather than replacing one with the
% other.
    e = NaN;
    [~, ~, ssr0] = chiOf(P, t, Y, w, A0, nH);
    if ~isfinite(ssr0); return; end
    lo = B(k,1); hi = B(k,2);
    scale = max(abs(P(k)), 0.02*(hi - lo));
    for side = [1 -1]
        step = 0.02*scale;
        last = 0;
        for it = 1:60
            Q = P; Q(k) = clipOne(P(k) + side*step, B(k,:));
            if Q(k) == P(k); break; end
            [~, ~, ssr] = chiOf(Q, t, Y, w, A0, nH);
            if ssr - ssr0 >= 1
                last = abs(Q(k) - P(k));
                break
            end
            last = abs(Q(k) - P(k));
            step = step * 1.35;
            if abs(step) > 2*(hi - lo); break; end
        end
        if last > 0
            if isnan(e); e = last; else; e = max(e, last); end
        end
    end
end

function v = flatModel(theta, nP, nT, t, A0, nH)
    P = theta(1:nP)';
    lin = reshape(theta(nP+1:end), 4, nT)';
    M = modelOf(P, t, lin, A0, nH);
    v = reshape(M', [], 1);
end

% ======================================================================
% DATA EXTRACTION
% ======================================================================

function [t, Y, E, halfPi, info] = pullRamsey(Exp, opts)
%PULLRAMSEY Get tau and both ratio traces out of an AWGExpRamsey.
%
% Deliberately NOT extract1Ddata. That applies three successive cuts to the ratio across
% averages, at 5, then 3, then ONE standard deviation -- and a one-sigma cut discards about
% a third of perfectly good Gaussian data, after two harsher passes have already shrunk the
% spread it measures against. Every error bar downstream then comes out too small.
% rabiFitRobust and esrFitRobust both refused it and so does this.
%
% What happens instead: the reduced traces and their sterr come from signalParam and
% signalParam2, and where the per-average ratios are reachable the error is taken as the
% LARGER of the propagated error and the scatter across averages. Those answer different
% questions -- photon statistics within an average, and everything that drifts between them
% -- and on this setup the second is routinely the bigger. Taking the larger is the only
% choice that does not understate the uncertainty.
    info = struct('nAvg', NaN, 'repeats', NaN, 'usedScatter', false, 'doubleMeasurement', false);
    halfPi = [];

    t = tryGet(Exp, 'tau', []);
    ax = tryGet(Exp, 'mCurrentXAxisParam', []);
    if ~isempty(ax)
        v = tryGet(ax, 'value', []);
        if ~isempty(v); t = v; end
    end
    if isempty(t); error('ramseyFitRobust: the experiment has no tau axis.'); end
    t = t(:)';

    halfPi = tryGet(Exp, 'halfPiTime', []);
    info.repeats = tryGet(Exp, 'repeats', NaN);
    info.nAvg    = tryGet(Exp, 'currIter', tryGet(Exp, 'averages', NaN));
    info.doubleMeasurement = isequal(tryGet(Exp, 'doubleMeasurement', false), true);

    [y1, e1] = paramData(Exp, 'signalParam');
    [y0, e0] = paramData(Exp, 'signalParam2');
    if isempty(y1); error('ramseyFitRobust: signalParam is empty -- did the run finish?'); end

    Y = y1; E = e1;
    if ~isempty(y0) && numel(y0) == numel(y1)
        Y = [y1; y0]; E = [e1; e0];
    end

    % Widen the error bars using the scatter across averages where it is bigger.
    S = tryGet(Exp, 'signal', []);
    nIter = tryGet(Exp, 'currIter', NaN);
    if ~isempty(S) && ndims(S) == 3 && isfinite(nIter) && nIter >= 3
        try
            for row = 1:size(Y, 1)
                dS = 2*row - 1; dR = 2*row;
                if size(S, 1) >= dR
                    r = squeeze(S(dS,:,1:nIter)) ./ max(squeeze(S(dR,:,1:nIter)), eps);
                    sc = std(r, 0, 2)' / sqrt(nIter);
                    bigger = sc > E(row, :);
                    if any(bigger)
                        E(row, bigger) = sc(bigger);
                        info.usedScatter = true;
                    end
                end
            end
        catch
        end
    end

    % Drop whole averages that sit far from the median, when there are enough to judge.
    % A MAD rule on 5+ samples applied ONCE -- note the contrast with the one-sigma cut above.
    info.droppedAverages = 0;
    if isfinite(opts.avgOutlierMad) && isfinite(nIter) && nIter >= 5 && ...
            ~isempty(S) && ndims(S) == 3
        try
            m = squeeze(mean(S(1,:,1:nIter) ./ max(S(2,:,1:nIter), eps), 2));
            md = median(m); mad0 = 1.4826*median(abs(m - md));
            if mad0 > 0
                bad = abs(m - md) > opts.avgOutlierMad*mad0;
                info.droppedAverages = sum(bad);
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

function [Y, E, labels] = pickProjections(Y, E, which)
    labels = {'|1\rangle'};
    if size(Y,1) >= 2; labels = {'|1\rangle', '|0\rangle'}; end
    switch lower(which)
        case 'one';  Y = Y(1,:); E = E(1,:); labels = labels(1);
        case 'zero'
            if size(Y,1) >= 2; Y = Y(2,:); E = E(2,:); labels = labels(2); end
        otherwise
            % both, if both are there
    end
end

% ======================================================================
% REPORT
% ======================================================================

function report(o)
    fprintf('\n=== ramseyFitRobust ===\n');
    switch o.nHyperfine
        case 3; iso = '14-N triplet';
        case 2; iso = '15-N doublet';
        otherwise; iso = 'single line (no hyperfine)';
    end
    fprintf('model: %s, A = %.3f MHz%s, %d projection(s) fitted jointly\n', ...
            iso, o.A, eitherText(isfinite(o.AErr), sprintf(' +/- %.3f (fitted)', o.AErr), ' (fixed)'), ...
            o.traces);
    fprintf('\n  T2*         %8.3f +/- %.3f us\n', o.T2star, o.T2starErr);
    fprintf('  detuning    %8.3f +/- %.3f MHz\n', o.delta, o.deltaErr);
    fprintf('  contrast    %8.3f +/- %.3f %%\n', 100*o.contrast, 100*o.contrastErr);
    if isfinite(o.contrastErr) && o.contrastErr >= o.contrast
        fprintf(['              *** the error is as big as the value, so the fringe is NOT\n' ...
                 '              MEASURED -- everything below is a description of noise. Raise\n' ...
                 '              repeats, and check halfPiTime before blaming the sample.\n']);
    end
    fprintf('  stretch p   %8.3f +/- %.3f\n', o.p, o.pErr);
    if o.p < 0.55 || o.p > 2.95
        fprintf(['              *** p is on its bound. The decay shape is not what this model\n' ...
                 '              can express -- suspect a second frequency, or a T2* far outside\n' ...
                 '              the tau span.\n']);
    end
    fprintf('  components  %s MHz\n', numList(o.freqs, '%.3f'));
    if o.delta < 3*o.deltaErr
        fprintf(['              the detuning is consistent with ZERO, so the oscillation you\n' ...
                 '              see is the %.3f MHz hyperfine beat alone. That is the expected\n' ...
                 '              on-resonance picture, not a problem.\n'], o.A);
    end

    fprintf('\n  resid acf   %8.2f     lag-1 autocorrelation of the residual\n', o.residAutocorr);
    % ONE CAVEAT ON READING IT: the threshold is not universal, it depends on how densely the
    % oscillation is sampled. At 15 points per beat -- the run-script default -- a correct
    % model gives about 0. At 7 points per period a correct model can give +0.6, because
    % adjacent residuals of a slightly imperfect frequency are genuinely correlated when the
    % sampling is that dense. So compare it against pointsPerBeat below, not against 0.4
    % blindly.
    % Only judged when there IS a residual to judge. On synthetic noiseless data the residual
    % is at machine precision and smooth, so the autocorrelation comes back near 1 and means
    % nothing at all -- gating on chi2 keeps that out of the report.
    if abs(o.residAutocorr) > 0.4 && o.chi2Reduced > 0.05
        fprintf(['              *** STRUCTURED. The residual walks rather than scatters, so the\n' ...
                 '              model is missing something and the errors above are optimistic.\n' ...
                 '              Phase is already fitted, freely and per projection -- so a\n' ...
                 '              residual that LOOKS phase-shifted is a FREQUENCY error, not a\n' ...
                 '              phase one, because a constant phase would have been absorbed.\n' ...
                 '              Run again with showResid = true to see it, and read the header\n' ...
                 '              note: the usual cause here is the hyperfine splitting A, not\n' ...
                 '              the detuning.\n']);
        if ~o.hyperfineFreed
            fprintf('              autoHyperfine did not free A. Try fitHyperfine = true.\n');
        end
    end
    if o.hyperfineFreed
        fprintf(['  A was FREED by autoHyperfine: %.3f -> %.3f MHz. The axial 14-N value is\n' ...
                 '  2.158, and an off-axis field shifts the effective splitting.\n'], ...
                o.AFixed, o.A);
    end
    fprintf('  chi2 (red)  %8.2f\n', o.chi2Reduced);
    fprintf('  R^2         %8.4f   (noise ceiling %.4f)\n', o.r2, o.r2Ceiling);
    if isfield(o, 'r2Kept') && o.nDown > 0
        fprintf('  R^2 kept    %8.4f   (ceiling %.4f) -- excluding the down-weighted points\n', ...
                o.r2Kept, o.r2CeilingKept);
    end
    rr = o.r2; rc = o.r2Ceiling;
    if isfield(o, 'r2Kept') && o.nDown > 0; rr = o.r2Kept; rc = o.r2CeilingKept; end
    if rc > 0 && rr < rc - 0.15
        fprintf('              well below the ceiling: the model is missing something real.\n');
    end
    fprintf('  span        %8.2f us = %.1f hyperfine beats, %.1f points per beat\n', ...
            o.spanUs, o.beatsSpanned, o.pointsPerBeat);
    if o.beatsSpanned < 2
        fprintf(['              *** under two beats. T2* and the beat are not separable on this\n' ...
                 '              span -- extend tau to at least %.2f us.\n'], 2/max(o.A, eps));
    end
    if o.T2star > 0.9*o.spanUs
        fprintf(['              *** T2* is at or beyond the tau span, so it is an extrapolation\n' ...
                 '              rather than a measurement. Extend tau past the decay.\n']);
    end
    if o.nDown > 0
        fprintf('\n  robust      %d point(s) down-weighted, at tau = %s us\n', ...
                o.nDown, numList(o.downTau, '%.3f'));
        fprintf(['              on this setup those are usually the classic path re-arming: a\n' ...
                 '              point whose repeats fired before the AWG resumed is diluted\n' ...
                 '              toward ratio = 1. Walk mode does not produce them.\n']);
        if o.nDown > 0.2*numel(o.t)
            fprintf('              *** that is over a fifth of the sweep. Look at the data, not the fit.\n');
        end
    end
    if o.traces == 2
        fprintf('\n  |0>/|1> amplitude ratio %+.3f   (should be -1)\n', o.ampRatio);
        if abs(o.ampRatio + 1) > 0.3
            fprintf(['              *** the two projections are not complementary. Before quoting\n' ...
                     '              T2*, check that rows 3-4 really are the |0> readout -- the\n' ...
                     '              two windows may not be doing what the class documents.\n']);
        end
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
% ONE axes by default. showResid and showSpectrum add panels below, and both are off unless
% asked for -- the main axes answers the question, and every extra panel shrinks it. The look
% comes from fitStyle, so this figure and the Rabi and ESR ones cannot drift apart.
    S = fitStyle();
    nExtra = double(o.opts.showResid) + double(o.opts.showSpectrum);
    fig = S.newFigure('Ramsey', nExtra > 0);

    if nExtra == 0
        axM = S.axes(fig);
    else
        axM = S.axes(subplot(3, max(nExtra,1), 1:(2*max(nExtra,1))));
    end

    td = linspace(min(o.t), max(o.t), 1500);
    lin = linFromModel(o);
    Md = modelOf(paramVec(o), td, lin, o.A, o.nHyperfine);

    H = []; L = {};
    for i = 1:o.traces
        H(end+1) = S.data(axM, o.t, o.y(i,:), o.e(i,:), S.slots{min(i,3)}, '', true, ...
                          o.downIdx);                                       %#ok<AGROW>
        L{end+1} = o.labels{i};                                            %#ok<AGROW>
    end
    for i = 1:o.traces
        h = S.model(axM, td, Md(i,:));
        if i == 1; H(end+1) = h; L{end+1} = 'fit'; end                     %#ok<AGROW>
    end
    % THE ENVELOPE, but only when there is an amplitude to draw an envelope OF. When the
    % fringe is not measured -- contrastErr as big as contrast -- the fitted amplitude is
    % whatever the noise allowed, and the envelope becomes two steep diagonals sweeping off
    % the top and bottom of the axes: a confident-looking decay curve drawn from a number
    % that means nothing. Suppressed rather than clipped, because a clipped fiction is still
    % a fiction.
    if isfinite(o.contrastErr) && o.contrastErr < 0.5*o.contrast
        D = exp(-(max(td,0)/max(o.T2star,eps)).^o.p);
        base = lin(1,1) + lin(1,2)*td;
        H(end+1) = S.envelope(axM, td, base + abs(o.amp(1))*D, base - abs(o.amp(1))*D);
        L{end+1} = sprintf('envelope, T_2^* = %.2f \\mus', o.T2star);
    end
    if o.nDown > 0
        H(end+1) = S.flag(axM, o.downTau, o.y(1, o.downIdx));
        L{end+1} = sprintf('%d down-weighted', o.nDown);
    end

    xlabel(axM, 'free evolution \tau  [\mus]', 'Color', S.col.ink2);
    ylabel(axM, yLabelOf(o, S), 'Color', S.col.ink2);
    % Scale to the MODEL and the KEPT points, not to everything -- see fitStyle/focusY.
    keepMask = true(1, numel(o.t)); keepMask(o.downIdx) = false;
    S.focus(axM, [Md(:); reshape(o.y(:, keepMask), [], 1)]);
    S.tight(axM, o.t); S.headroom(axM, 0.30);
    S.legend(axM, H, L);
    S.title(axM, titleLines(o));

    k = 0;
    if o.opts.showResid
        k = k + 1;
        axR = S.axes(subplot(3, nExtra, 2*nExtra + k));
        S.zeroLine(axR, [min(o.t) max(o.t)]);
        Rn = o.resid ./ max(o.e, eps);
        for i = 1:o.traces
            S.resid(axR, o.t, Rn(i,:), S.slots{min(i,3)});
        end
        % Same reason as the main panel: two 40-sigma outliers flatten every real residual
        % into a line at zero. Scaled to the kept residuals; the outliers are in the count.
        S.focus(axR, reshape(Rn(:, keepMask), [], 1), 0.15);
        xlabel(axR, '\tau  [\mus]', 'Color', S.col.ink2);
        ylabel(axR, 'residual / \sigma', 'Color', S.col.ink2);
        S.tight(axR, o.t);
        S.title(axR, sprintf('\\chi^2_\\nu = %.2f, lag-1 autocorr %+.2f%s', ...
                             o.chi2Reduced, o.residAutocorr, ...
                             eitherText(abs(o.residAutocorr) > 0.4, ...
                                        '  -- structured, model is missing something', '')));
    end
    if o.opts.showSpectrum
        k = k + 1;
        axS = S.axes(subplot(3, nExtra, 2*nExtra + k));
        [fAx, PA] = fringeSpectrum(o);
        if ~isempty(fAx)
            S.band(axS, fAx, PA, zeros(size(PA)));
            plot(axS, fAx, PA, '-', 'Color', S.col.s1, 'LineWidth', 1.2);
            ylim(axS, [0 1.3]);
            for kk = 1:numel(o.freqs)
                plot(axS, [o.freqs(kk) o.freqs(kk)], [0 1.3], '--', ...
                     'Color', S.col.flag, 'LineWidth', 1.2);
            end
            xlim(axS, [0, min(max(fAx), 3*max(max(o.freqs), o.A) + 1)]);
        end
        xlabel(axS, 'frequency  [MHz]', 'Color', S.col.ink2);
        ylabel(axS, '|amplitude|, normalised', 'Color', S.col.ink2);
        S.title(axS, sprintf('dashed = fitted %s MHz', numList(o.freqs, '%.2f')));
    end
end

function s = yLabelOf(o, S)
    s = S.yLabel;
    if isfield(o.opts, 'yLabel') && ~isempty(o.opts.yLabel); s = o.opts.yLabel; end
end

function r = rKept(o)
% R^2 over the points the ROBUST pass kept. Two gross outliers drag the all-points R^2 far
% below the ceiling on data the fit reproduces perfectly, which reads as a failed fit and is
% not one.
    r = o.r2;
    if isfield(o, 'r2Kept') && o.nDown > 0; r = o.r2Kept; end
end

function r = rCeil(o)
    r = o.r2Ceiling;
    if isfield(o, 'r2CeilingKept') && o.nDown > 0; r = o.r2CeilingKept; end
end

function T = titleLines(o)
    a = sprintf(['T_2^* = %.3f \\pm %.3f \\mus     detuning = %.3f \\pm %.3f MHz     ' ...
                 'contrast = %.3f \\pm %.3f %%'], ...
                o.T2star, o.T2starErr, o.delta, o.deltaErr, ...
                100*o.contrast, 100*o.contrastErr);
    bits = {sprintf('A = %.3f MHz%s', o.A, eitherText(o.hyperfineFreed, ' (fitted)', '')), ...
            sprintf('p = %.2f', o.p), sprintf('\\chi^2_\\nu = %.2f', o.chi2Reduced), ...
            sprintf('R^2 = %.3f of %.3f', rKept(o), rCeil(o)), ...
            sprintf('resid acf %+.2f', o.residAutocorr)};
    if o.traces == 2
        bits{end+1} = sprintf('|0\\rangle/|1\\rangle = %+.2f', o.ampRatio);
    end
    if o.nDown > 0
        bits{end+1} = sprintf('%d down-weighted (may be off scale)', o.nDown);
    end
    T = {a, strjoin(bits, '     |     ')};
end

function [fAx, PA] = fringeSpectrum(o)
% A plain DFT of the baseline-removed first trace, on a uniform grid. Not used by the fit --
% it cannot seed it, for the reason in gridSearch -- but it is the one picture that shows
% whether the three lines the model assumes are the three lines in the data.
    fAx = []; PA = [];
    t = o.t; y = o.y(1,:);
    dt = median(diff(t));
    if ~isfinite(dt) || dt <= 0 || numel(t) < 8; return; end
    lin = linFromModel(o);

    % Down-weighted points are REPLACED BY THE FIT before transforming. A point sitting
    % 25 sigma off is a delta function, and a delta function is white -- it puts a flat
    % pedestal across the whole spectrum and buries the three lines this panel exists to
    % show. Substituting the model keeps the record length and the sampling intact without
    % letting an artefact set the scale.
    if o.nDown > 0
        y(o.downIdx) = o.model(1, o.downIdx);
    end
    y = y - (lin(1,1) + lin(1,2)*t);            % remove the fitted baseline, not the fringe

    nfft = 2^nextpow2(8*numel(t));
    Yf = abs(fft(y .* hann(numel(y)), nfft));
    fAx = (0:nfft-1) / (nfft*dt);
    half = 1:floor(nfft/2);
    fAx = fAx(half); PA = Yf(half);

    % Drop the lowest bins. Any residual baseline lands at DC with an amplitude far above the
    % fringe, and one DC spike flattens everything else to invisibility -- which is what the
    % first version of this panel showed. Below one cycle per record there is no frequency
    % information anyway.
    fMin = 0.8 / max(max(t) - min(t), eps);
    drop = fAx < fMin;
    fAx = fAx(~drop); PA = PA(~drop);
    if isempty(PA); return; end
    PA = PA / max(max(PA), eps);                % normalised: only the shape matters here
end

function wnd = hann(n)
% Written out: hann lives in the Signal Processing Toolbox and this file has no toolboxes.
    if n < 2; wnd = ones(1, n); return; end
    wnd = 0.5 - 0.5*cos(2*pi*(0:n-1)/(n-1));
end

function lin = linFromModel(o)
    lin = zeros(o.traces, 4);
    for i = 1:o.traces
        lin(i,1) = o.baseline(i);
        lin(i,2) = o.slope(i);
        lin(i,3) = o.amp(i)*cos(o.phi(i));
        lin(i,4) = -o.amp(i)*sin(o.phi(i));
    end
end

function P = paramVec(o)
    P = [o.delta, o.T2star, o.p];
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
