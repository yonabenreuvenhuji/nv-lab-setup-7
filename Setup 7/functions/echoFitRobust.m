function out = echoFitRobust(args)
%ECHOFITROBUST Robust NV Hahn-echo / XY-N fit: stretched exponential, T2, and the 13-C revivals.
%
%   out = echoFitRobust(struct('exp', e))
%   out = echoFitRobust(struct('exp', e, 'revivals', 'on'))
%   out = echoFitRobust(struct('t', t, 'y', y, 'e', e))        % plain vectors
%   out = echoFitRobust(struct('t', t, 'y', Y, 'e', E))        % Y is 2xN: |1> and |0>
%
% ONE STRUCT IN, ONE STRUCT OUT, the same convention as ramseyFitRobust, rabiFitRobust and
% esrFitRobust: the data rides in the same struct as the options -- `exp` for an experiment
% object, or `t`/`y`/`e` for vectors -- and everything comes back in `out`, with out.ok saying
% whether to believe it.
%
% Self-contained. No Curve Fitting Toolbox, no Signal Processing Toolbox, no lab fitting
% helpers, no EchoFit, no extract1Ddata, no findpeaks. The figure comes from fitStyle, shared
% with the other three robust fits, so they cannot drift apart.
%
% ===================================================================
% WHICH TIME AXIS, AND WHY IT MATTERS
% ===================================================================
%
% AWGExpEcho publishes TWO x axes and they differ by a factor that is easy to lose:
%
%       obj.topParam.value           = tau                      the per-arm free evolution
%       obj.mCurrentXAxisParam.value = 2 * tau * cycles         the TOTAL evolution time
%
% The decay is a function of the TOTAL time under the bath, so that is what is fitted, and
% out.t is total evolution time in microseconds. A T2 quoted against tau instead is short by
% 2*cycles -- a factor of 2 for a plain Hahn echo and 16 for an XY8. Set xUnits = 'tau' if you
% deliberately want tau on the axis; out.T2 then means a tau, and the report says so rather
% than letting the label carry the whole difference.
%
% ===================================================================
% THE MODEL
% ===================================================================
%
%       y(t) = c0 + [c1*t] + a * exp(-(t/T2)^p) * R(t)
%
%   c0      the PLATEAU: where the echo ends up once the coherence is gone. Fitted, not
%           assumed to be 1. On this setup it is not 1.
%   a       the coherence amplitude at t = 0, WITH ITS SIGN. |1> and |0> carry opposite signs.
%   T2, p   stretched exponential, p FITTED (bounds 0.5 to 4) rather than pinned. See below --
%           on an echo p is not a nuisance parameter, it is a measurement.
%   R(t)    the 13-C revival comb, off unless the data asks for it. See below.
%   c1      a linear drift. OFF BY DEFAULT, and that is not laziness -- see the degeneracy
%           note below, which is the one thing about this fit that is unlike the Ramsey.
%
% ===================================================================
% p IS THE POINT, NOT A NUISANCE
% ===================================================================
%
% Everyone fits an echo with exp(-(t/T2)^p) and then fixes p at 1, or at 2, or at 3, according
% to what they were taught. Fixing it throws away the one number in the fit that says what the
% coherence is limited BY:
%
%   p ~ 1     the bath fluctuates fast compared with T2 -- motional narrowing, or a
%             relaxation-limited echo. If T1 is known, check T2 against 2*T1 (pass T1us and
%             the report does it): an echo at p ~ 1 sitting near the 2*T1 ceiling is not
%             measuring dephasing at all, it is measuring T1.
%   p ~ 2     slow Gaussian noise -- the Ornstein-Uhlenbeck bath in its slow limit, and the
%             usual shape when a technical noise source dominates.
%   p ~ 3     a dipolar-coupled nuclear spin bath. This is the textbook Hahn-echo exponent
%             for an NV in natural-abundance 13-C, and it is what a clean single-NV echo on
%             this setup should return.
%
% So p is reported with an error bar and an interpretation, and if it lands on a bound the
% report says the decay shape is outside what the model can express rather than quoting the
% bound as a result. A p that comes back at 1.2 when you expected 3 is a real finding about
% the sample or the noise floor; a p that was fixed at 3 hides it.
%
% Because p is free, T2 ALONE IS NOT COMPARABLE BETWEEN RUNS. T2 = 25 us at p = 3 and
% T2 = 25 us at p = 1 are different amounts of coherence. The shape-independent figure of
% merit is the area under the envelope,
%
%       coherenceArea = integral_0^inf exp(-(t/T2)^p) dt = T2 * gamma(1 + 1/p)
%
% which is reported alongside T2 and is the number to compare when p has moved.
%
% ===================================================================
% THE 13-C REVIVALS, WHICH ARE THIS FIT'S HYPERFINE BEAT
% ===================================================================
%
% ramseyFitRobust exists because a single damped cosine fitted to a Ramsey reports a detuning
% that is not there. The echo has the exact same trap in a different costume.
%
% In a field, the 13-C nuclei precess at their Larmor frequency, and the echo periodically
% loses and recovers coherence in step with them: COLLAPSES between, REVIVALS at multiples of
% the Larmor period. A monotone stretched exponential fitted through a collapse-and-revival
% trace does not fail loudly. It returns a short T2 with a small error bar and a plausible p,
% because the first collapse looks exactly like the end of a decay. That is the same failure
% mode as the 2.16 MHz detuning on an on-resonance Ramsey: confidently wrong, not obviously
% broken.
%
% So the modulation is IN the model rather than fitted around:
%
%       R(t) = SUM_k exp(-((t - k*Trev)/w)^2)   normalised so R(0) = 1
%
% a comb of Gaussian revivals of period Trev and width w. Not the smooth 1 - m*sin^2(pi*t/Trev)
% form, which is easier to write and wrong in shape: a real 13-C collapse is BROAD and near
% total while the revival is NARROW, and a sinusoid gives them equal width by construction.
% The width is parametrised as a fraction of the period (wFrac = w/Trev) so that the two
% parameters have independent box bounds and the search stays rectangular.
%
% REVIVALS ARE OFF UNLESS THE DATA PAYS FOR THEM. With revivals = 'auto' (the default) the
% monotone model is fitted first, then the comb, and the comb is accepted only if all three
% hold:
%
%   1. the weighted SSR drops by more than 11.83 -- the 3-sigma joint threshold for the TWO
%      extra parameters, so noise cannot buy it (one parameter would be 9; this is chi2inv
%      at 0.9973 with 2 degrees of freedom);
%   2. wFrac < 0.4, i.e. the revivals are actually narrow. At wFrac -> 0.5 the Gaussians
%      overlap into a constant and "revivals" would be a two-parameter way of writing 1;
%   3. THE REVIVAL IS IN THE DATA AND SAMPLED: the last revival peak lands inside the sweep
%      together with its falling edge (k*Trev + w <= max(t)), and at least four points sit
%      within one width of it.
%
% All three are reported, pass or fail, so a refused comb is visible rather than silent.
%
% Test 3 used to read "at least 1.8 periods fit inside the sweep", and that was a mistake worth
% recording. The comb has a peak PINNED AT t = 0, so a single revival away from the origin is a
% COMPLETE period, anchored at both ends -- not half of one. Demanding two threw away a real
% measurement: a 60-point echo out to 120 us of evolution collapsed by 20 us and plainly
% revived at 105 us, both projections moving in opposite directions across fifteen points, and
% the fit reported T2 = 11.8 us -- the width of the first collapse -- because 105 was outside
% both the grid and the box. The guard meant to stop the model inventing revivals stopped it
% seeing one.
%
% ONE REVIVAL IS ENOUGH TO MEASURE THE PERIOD, AND NOT ENOUGH TO MEASURE T2 AND p BOTH. With a
% single revival the envelope rests on one number -- that revival's height relative to t = 0 --
% and exp(-(T/T2)^p) has two unknowns in it, with the comb width competing for the same
% variance as a third. T2 and p then slide along a curve of equally good fits. The report says
% so explicitly and names the two ways out (extend the sweep to catch a second revival, or fix
% p and quote T2 conditionally); testEchoFit's t12 measures the slide, and on that trace T2
% runs 371 -> 202 -> 163 us across p = 1, 2, 3 at a reduced chi2 that never moves past 0.72.
%
% WHAT THE PERIOD MEANS. Trev comes out in whatever units out.t is in -- total evolution time
% by default. The standard two-pulse ESEEM result is periodic in the PER-ARM evolution with
% period 1/f_Larmor, so on a total-time axis
%
%       f_Larmor = 2*cycles / Trev        and     B = f_Larmor / 1.0705e-3   gauss
%
% using gamma(13-C)/2pi = 10.705 MHz/T. Both are reported, and the honest use of them is as a
% CHECK: compare the implied field against the coil setting. If it comes back at half or twice
% the field you set, the factor above is the wrong one for this pulse sequence, and Trev
% itself -- which is what was measured -- is still right. Nothing downstream of T2 depends on
% the conversion.
%
% ===================================================================
% WHY THE DRIFT TERM IS OFF BY DEFAULT
% ===================================================================
%
% This is the one place where copying the Ramsey fit would be a mistake.
%
% In a Ramsey the signal is an oscillation, so a linear drift c1*t is very nearly orthogonal
% to it: the fringe and the trend live in different places and freeing the drift costs almost
% nothing. In an echo the signal IS a monotone decay. c1*t and a*exp(-(t/T2)^p) are then
% competing to explain the same downward slope, and the fit has no way to tell a slow decay
% from a laser that dimmed during the sweep.
%
% How bad it is depends entirely on the SHAPE, which is why a blanket "the drift term is
% harmless" is wrong. Measured in testEchoFit, on synthetic data with NO drift injected,
% comparing fitDrift false against true (median over 5 noise seeds):
%
%       decay shape                          T2 moves    T2 error bar
%       p = 3, decay well inside the sweep       0.3%          x1.3
%       p = 1, decay well inside the sweep       2.5%          x3.9
%       p = 1, sweep only 1.3x T2              145.2%      not computable
%
% Read the second COLUMN as much as the first. A steep p = 3 fall really is distinguishable
% from a straight line, so on a clean 13-C-limited echo the drift term costs almost nothing.
% A p = 1 decay is very nearly a straight line over the first part of the sweep, and there the
% fitted T2 barely moves -- 2.5% -- while its error bar nearly QUADRUPLES: the fit has not
% become wrong, it has become uninformative, and only the error bar says so. That is the case
% to watch for, because a 2.5% shift looks like nothing if you only compare the values.
%
% The third row is an echo whose sweep barely reaches past T2, which is the normal state of a
% first measurement. There the fit moves T2 by more than T2 and the error bar goes degenerate:
% pickErr discards it rather than print a reassuring number, so the report says NO ERROR BAR,
% which is the correct answer.
%
% The numbers above are printed by testEchoFit's t8 on every run, so they cannot go stale
% while this comment claims otherwise.
%
% So fitDrift defaults to FALSE, the plateau is a constant, and turning it on is a deliberate
% act for when there is independent reason to think the laser moved -- in which case read the
% widened error bar as the honest price rather than a regression.
%
% ===================================================================
% ROBUSTNESS, WHICH THIS DATA NEEDS
% ===================================================================
%
% The classic path re-arms the AWG at every point, and repeats that fire before the instrument
% resumes see no microwave, so that point is diluted toward "no contrast". They show up as
% isolated points pulled toward the plateau -- and on an echo, a point pulled toward the
% plateau is indistinguishable by eye from an honest early loss of coherence, which makes them
% more dangerous here than in a Ramsey where they break an obvious fringe. The last point of
% the classic reference run in the echoGateTest output sat about 11e-4 off, which is roughly
% 9 error bars.
%
% So after an initial fit the weights are re-derived with a Tukey bisquare (IRLS, 4.685 sigma):
% ordinary scatter is left alone, a point many sigma off the curve loses its weight. The report
% names how many were down-weighted and at which t, and the figure rings them -- a fit that
% quietly discards a fifth of the sweep is worse than one that does not converge. Walk mode
% does not produce them at all, which is the whole reason walk mode was built.
%
% ===================================================================
% BOTH PROJECTIONS AT ONCE
% ===================================================================
%
% With doubleMeasurement on, AWGExpEcho.publishEcho puts the |1> projection (detection rows
% 1-2) in signalParam and the |0> one (rows 3-4) in signalParam2. Same physics, two readout
% windows.
%
% They are fitted JOINTLY on T2, p, Trev and wFrac -- and NOT on the plateau or the amplitude,
% because those are properties of a readout window rather than of the spin. That roughly halves
% the error on T2, and it turns the amplitude ratio a(|0>)/a(|1>) into a free consistency check
% which should come back at -1. If it does not, the two windows are not doing what the class
% documents, and that is worth knowing before a T2 is quoted. Pass useProjection to fit one.
%
% ===================================================================
% THE SHORT-TIME DEFICIT, AND WHY IT IS IN HERE
% ===================================================================
%
% Everything the gate tests were about lands on this fit. A microwave gate that leaks through
% an open switch removes contrast at SHORT evolution times, and lost contrast at short t is
% arithmetically indistinguishable from a shorter T2 -- the failure mode is a wrong number,
% not a broken run.
%
% So out.shortDeficitSigma reports the weighted mean residual over the first quarter of the
% sweep, in sigma of that mean. Negative means the data sits BELOW the model early, which is
% the leak signature. echoGateTest has since measured the gate as harmless out to 240 us of
% excess at the 17% level, so a deficit that shows up here is more likely a pulse-calibration
% error, an unmodelled fast decay component, or a p the single stretched exponential cannot
% reach -- all three are named in the report rather than left to guesswork.
%
% ===================================================================
% OPTIONS  (defaults in brackets)
% ===================================================================
%
%   exp []                  an AWGExpEcho, or
%   t, y, e []              vectors. y may be 1xN or 2xN.
%   xUnits ['evolution']    'evolution' = 2*tau*cycles, or 'tau'. Only used when reading an
%                           experiment object; with vectors, t is whatever you passed.
%   cycles []               overrides the object's cycles, used for the Larmor conversion only
%   useProjection ['both']  'both' | 'one' (|1>) | 'zero' (|0>)
%   revivals ['auto']       'auto' | 'on' | 'off'
%   revivalPeriodUs []      fix Trev instead of fitting it
%   revivalWidthFrac []     fix w/Trev instead of fitting it
%   pFixed []               fix the stretch exponent, e.g. 3. [] fits it.
%   fitDrift [false]        free the linear drift. Read the degeneracy note first.
%   robust [true]           Tukey bisquare re-weighting
%   tMin []                 drop points below this t. [] derives it from the pulse times.
%   tMax []                 drop points above this t
%   T1us []                 known T1, for the T2 <= 2*T1 check
%   nGridT2 [28]            T2 grid points
%   nGridRev [240]          revival-period grid points
%   avgOutlierMad [4]       MAD cut on whole averages, applied ONCE
%   plotIt [true], verbose [true], showResid [false], yLabel ['']
%
% ===================================================================
% OUT
% ===================================================================
%
%   ok, why                 did it work, and if not why
%   T2, T2Err               microseconds, on the axis named by tUnits
%   p, pErr                 stretch exponent, and what it implies -- see pMeaning
%   pMeaning                one line of text
%   coherenceArea           T2*gamma(1+1/p) us, the shape-independent figure of merit
%   contrast, contrastErr   |a|/|plateau| at t = 0, from the first projection. NOTE: this is
%                           peak-to-PLATEAU, where ramseyFitRobust's contrast is
%                           peak-to-PEAK. The two are not the same quantity; an echo has no
%                           fringe to take a peak-to-peak of. A Ramsey contrast of 1.2% and an
%                           echo contrast of 1.2% mean different amounts of signal.
%   revivalsUsed            did the comb earn its place
%   revivalPeriodUs, revivalPeriodErr, revivalWidthUs, revivalWidthFrac
%   larmorMHz, fieldGauss   implied 13-C Larmor and field. A CHECK, not a measurement.
%   revivalTest             struct: dSSR, threshold, wFrac, periodsInSpan, and which failed
%   plateau, amp, drift     per projection
%   ampRatio                a(|0>)/a(|1>), should be -1
%   decayErrorBars          how many error bars the fitted decay spans -- is there a decay
%   shortDeficit, shortDeficitSigma
%   T2overTwoT1             T2/(2*T1) if T1us was given
%   residAutocorr           lag-1 autocorrelation of the residual
%   chi2Reduced, r2, r2Ceiling, r2Kept, r2CeilingKept
%   t, y, e, model, envelope, resid, weights, nDown, downIdx, downT, traces, labels
%   spanUs, tUnits, tMin, tMax, nDropped, nyquistMHz, pointsPerRevival
%   coeffError              [T2Err, contrastErr, pErr], for anything that indexes it

    % ================================================================
    % ONE STRUCT IN
    % ================================================================
    DEFAULTS = struct('exp', [], 't', [], 'y', [], 'e', [], ...
        'xUnits', 'evolution', 'cycles', [], 'useProjection', 'both', ...
        'revivals', 'auto', 'revivalPeriodUs', [], 'revivalWidthFrac', [], ...
        'pFixed', [], 'fitDrift', false, 'robust', true, ...
        'tMin', [], 'tMax', [], 'T1us', [], ...
        'nGridT2', 28, 'nGridRev', 240, 'avgOutlierMad', 4, ...
        'plotIt', true, 'verbose', true, 'showResid', false, 'yLabel', '');

    if ~exist('args', 'var') || isempty(args); args = struct(); end
    if ~isstruct(args)
        error(['echoFitRobust takes ONE struct: echoFitRobust(struct(''exp'', e)) or ' ...
               'echoFitRobust(struct(''t'', t, ''y'', y, ''e'', e)). There is no positional ' ...
               'form -- see ramseyFitRobust for why.']);
    end
    opts = withDefaults(args, DEFAULTS, 'echoFitRobust');

    out = struct('ok', false, 'why', '', 'opts', opts);
    Exp = opts.exp;
    fromVectors = isempty(Exp);
    tUnits = 'total evolution';
    cyc = opts.cycles;
    pulseUs = [];

    pullInfo = [];
    if fromVectors
        if isempty(opts.t) || isempty(opts.y)
            out.why = ['nothing to fit: give exp, or t and y. ' ...
                       'echoFitRobust(struct(''exp'', e)).'];
            fprintf('\nechoFitRobust FAILED: %s\n', out.why); return
        end
        t = opts.t; Y = opts.y; E = opts.e;
    else
        try
            [t, Y, E, info] = pullEcho(Exp, opts);
            pullInfo = info;
            if isempty(cyc); cyc = info.cycles; end
            pulseUs = info.pulseUs;
            tUnits = info.tUnits;
        catch err
            describeExp(Exp);
            rethrow(err);
        end
    end
    if isempty(cyc) || ~isfinite(cyc); cyc = 1; end

    t = t(:)';
    if isvector(Y); Y = Y(:)'; end
    if isvector(E); E = E(:)'; end
    if size(Y, 2) ~= numel(t)
        if size(Y, 1) == numel(t); Y = Y.'; E = E.'; else
            out.why = sprintf('%d times against a %dx%d signal', numel(t), size(Y,1), size(Y,2));
            fprintf('\nechoFitRobust FAILED: %s\n', out.why); return
        end
    end
    if isempty(E) || all(~isfinite(E(:))); E = ones(size(Y)); end
    if ~isequal(size(E), size(Y)); E = repmat(median(E(isfinite(E) & E>0)), size(Y)); end

    [Y, E, labels] = pickProjections(Y, E, opts.useProjection);

    % ---- window ----------------------------------------------------------
    tMin = opts.tMin;
    if isempty(tMin)
        if ~isempty(pulseUs) && isfinite(pulseUs) && pulseUs > 0
            tMin = pulseUs;
        else
            tMin = 0;
        end
    end
    tMax = opts.tMax; if isempty(tMax); tMax = Inf; end
    keep = isfinite(t) & (t >= tMin) & (t <= tMax) & all(isfinite(Y), 1);
    nDrop = sum(~keep);
    if nDrop > 0 && opts.verbose
        fprintf(['  dropped %d point(s) outside t = [%.4f, %.4g] us: below tMin the pi/2 and pi\n' ...
                 '  pulses still overlap, so those points are not free evolution and do not\n' ...
                 '  belong to this model. Set tMin = 0 to keep them.\n'], nDrop, tMin, tMax);
    end
    t = t(keep); Y = Y(:, keep); E = E(:, keep);
    n = numel(t); nT = size(Y, 1);
    if n < 8
        out.why = sprintf(['only %d usable points. A T2 needs enough of them to separate a ' ...
                           'decay from a plateau -- 15 or more, reaching past the decay.'], n);
        fprintf('\nechoFitRobust FAILED: %s\n', out.why); return
    end
    E(~isfinite(E) | E <= 0) = median(E(isfinite(E) & E > 0));

    dt = median(diff(sort(t)));
    span = max(t) - min(t);
    out.nyquistMHz = 1/(2*max(dt, eps));

    % ================================================================
    % THE FIT -- MONOTONE FIRST, ALWAYS
    % ================================================================
    % The monotone fit is not a starting guess for the revival fit. It is the NULL HYPOTHESIS
    % the revival fit has to beat, and both are kept and reported.
    %
    % ORDER MATTERS HERE. The two models are compared under the PLAIN inverse-variance
    % weights, BEFORE any robust re-weighting, because weights derived from one of two
    % competing models rig the comparison in its favour. On collapse-and-revival data the
    % monotone fit is badly wrong, so its Tukey weights land on the points around the first
    % revival -- exactly the points that distinguish the two models -- and the comb is then
    % asked to justify itself on data the null hypothesis has already been allowed to delete.
    %
    % Measured on the synthetic 18.7 us revival trace in testEchoFit: the monotone fit
    % down-weights 10 of the 120 points, all of them on the revivals, and the comb's SSR
    % improvement falls from 90093 to 19904 -- a factor of 4.5 of the evidence thrown away by
    % nothing but the order of two steps. It did NOT flip the decision on that trace, and this
    % comment does not claim it would; 19904 is still far past the 11.83 threshold. It would
    % take a shallower modulation to flip, and the point of doing it in the right order is
    % that nobody has to work out where that boundary is.
    %
    % So: fit both on plain weights, choose, and only THEN re-weight robustly on the model
    % that was chosen.
    nLin = 2 + double(opts.fitDrift);          % [plateau, (drift), envelope]
    w0 = 1 ./ E.^2;
    fitP = isempty(opts.pFixed);
    pFix = opts.pFixed;

    mono = fitMonotone(t, Y, w0, span, opts, nLin, fitP, pFix);
    P = mono.P; useRev = false; rev = [];

    % ================================================================
    % DO THE REVIVALS EARN THEIR PLACE?
    % ================================================================
    wantRev = ~strcmpi(opts.revivals, 'off');
    forceRev = strcmpi(opts.revivals, 'on') || strcmpi(opts.revivals, 'true');
    test = struct('tried', false, 'dSSR', NaN, 'threshold', 11.83, 'wFrac', NaN, ...
                  'periodsInSpan', NaN, 'revivalsVisible', NaN, 'pointsOnPeak', NaN, ...
                  'failed', '', 'accepted', false);
    if wantRev
        test.tried = true;
        [~, ~, ssrMono] = chiOf(P, t, Y, w0, nLin, false);
        rev = fitRevival(t, Y, w0, span, dt, opts, nLin, P, fitP, pFix);
        if isempty(rev)
            test.failed = 'no revival period could be gridded on this span';
        else
            [~, ~, ssrRev] = chiOf(rev.P, t, Y, w0, nLin, true);
            test.dSSR = ssrMono - ssrRev;
            test.wFrac = rev.P(4);
            test.periodsInSpan = span / max(rev.P(3), eps);
            % HOW MANY REVIVALS ARE ACTUALLY IN THE DATA, and is the last one sampled. These
            % replace the old "1.8 periods in the span" count, which was the wrong quantity:
            % the comb is anchored at t = 0, so what matters is not how many periods fit but
            % whether a revival PEAK lands inside the sweep with enough points on it to be
            % more than a rumour.
            Trv = rev.P(3); wUs = rev.P(4)*Trv;
            test.revivalsVisible = floor(max(t)/max(Trv, eps));
            if test.revivalsVisible >= 1
                tPk = test.revivalsVisible * Trv;
                test.pointsOnPeak = sum(abs(t - tPk) <= wUs);
            else
                tPk = NaN; test.pointsOnPeak = 0;
            end
            fails = {};
            if test.dSSR <= 0
                % Not a near miss, and worth wording carefully. The comb is NOT a strict
                % generalisation of the monotone model -- a normalised sum of Gaussians can
                % approach a constant but never equal one -- so on an unmodulated trace the
                % best comb comes out measurably WORSE rather than equal. Printing "SSR
                % dropped by -136" invites reading a negative number as a near miss. It is
                % not; there is no modulation in the trace at all.
                fails{end+1} = sprintf(['the comb made the fit WORSE (weighted SSR up by %.1f) -- ' ...
                                        'there is no modulation in this trace'], -test.dSSR);
            elseif ~(test.dSSR > test.threshold)
                fails{end+1} = sprintf('SSR dropped by only %.1f, under the %.2f needed for two extra parameters', ...
                                       test.dSSR, test.threshold);
            end
            if ~(test.wFrac < 0.4)
                fails{end+1} = sprintf('the revivals came out wide (w/Trev = %.2f, needs < 0.40) -- that is a constant, not a comb', ...
                                       test.wFrac);
            end
            if test.revivalsVisible < 1
                fails{end+1} = sprintf(['the fitted period (%.3g us) is longer than the sweep, so no ' ...
                                        'revival peak is in the data at all -- that is a shoulder, ' ...
                                        'not a comb'], Trv);
            elseif ~(tPk + wUs <= max(t))
                fails{end+1} = sprintf(['the last revival (%.3g us) plus its width (%.3g us) runs off ' ...
                                        'the end of the sweep, so its falling edge was never measured'], ...
                                       tPk, wUs);
            elseif test.pointsOnPeak < 4
                fails{end+1} = sprintf(['only %d point(s) sit within one width of the revival at %.3g us ' ...
                                        '(needs 4) -- it is not sampled'], test.pointsOnPeak, tPk);
            end
            if isempty(fails) || forceRev
                useRev = true; P = rev.P;
                test.accepted = true;
                if ~isempty(fails)
                    test.failed = ['FORCED ON by revivals = ''on'' despite: ' strjoin(fails, '; ')];
                end
            else
                test.failed = strjoin(fails, '; ');
            end
        end
    end

    % ================================================================
    % ROBUST RE-WEIGHTING, ON THE MODEL THAT WAS CHOSEN
    % ================================================================
    w = w0; nDown = 0; downIdx = [];
    if opts.robust
        [~, lin] = chiOf(P, t, Y, w0, nLin, useRev);
        R = (Y - modelOf(P, t, lin, nLin, useRev)) ./ E;
        s = 1.4826 * median(abs(R(:) - median(R(:))));
        if s > 0
            u = R / (4.685 * s);
            bw = (abs(u) < 1) .* (1 - u.^2).^2;
            % bw < 0.25 is |residual| > 3.3 sigma. Same threshold as ramseyFitRobust, and for
            % the same reason: 2.5 sigma calls a couple of honest points outliers in any sweep
            % of this length, so the first version reported 8 "outliers" on data with 2. The
            % WEIGHTING tapers smoothly from zero either way; this number only decides what
            % gets NAMED an outlier and ringed on the figure.
            downIdx = find(any(bw < 0.25, 1));
            nDown = numel(downIdx);
            w = w0 .* max(bw, 1e-3);
            if useRev
                P = nelderMead(@(Q) chiOf(Q, t, Y, w, nLin, true), P, ...
                               revBounds(span, dt, fitP, pFix, opts, max(t)));
            else
                mono = fitMonotone(t, Y, w, span, opts, nLin, fitP, pFix);
                P = mono.P;
            end
        end
    end

    % WHAT THE MONOTONE FIT WOULD HAVE SAID, under the FINAL weights. Refitted rather than
    % carried over from the plain-weight pass so that the T2 printed as "without the comb" is
    % comparable point for point with the T2 printed above it. On revival data the gap between
    % them is the whole argument for the comb being in the model.
    monoFinal = mono;
    if useRev
        monoFinal = fitMonotone(t, Y, w, span, opts, nLin, fitP, pFix);
    end

    [chi2, lin] = chiOf(P, t, Y, w, nLin, useRev);
    model = modelOf(P, t, lin, nLin, useRev);
    resid = Y - model;

    % ================================================================
    % PARAMETERS AND WHAT THEY MEAN
    % ================================================================
    o = struct();
    o.T2 = P(1); o.p = P(2);
    o.tUnits = tUnits; o.cycles = cyc;
    o.revivalsUsed = useRev;
    o.revivalTest = test;
    o.revivalPeriodUs = NaN; o.revivalWidthFrac = NaN; o.revivalWidthUs = NaN;
    o.larmorMHz = NaN; o.fieldGauss = NaN; o.pointsPerRevival = NaN;
    o.revivalsVisible = test.revivalsVisible;
    if useRev
        o.revivalPeriodUs = P(3); o.revivalWidthFrac = P(4);
        o.revivalWidthUs  = P(3)*P(4);
        % gamma(13-C)/2pi = 10.705 MHz/T = 1.0705e-3 MHz/G. The 2*cycles is the total-time to
        % per-arm conversion; see the header on why this is a CHECK and not a measurement.
        o.larmorMHz  = (2*cyc) / max(P(3), eps);
        if strcmpi(tUnits, 'tau'); o.larmorMHz = 1 / max(P(3), eps); end
        o.fieldGauss = o.larmorMHz / 1.0705e-3;
        o.pointsPerRevival = P(3) / max(dt, eps);
    end

    o.plateau = zeros(1, nT); o.amp = zeros(1, nT); o.drift = zeros(1, nT);
    for i = 1:nT
        o.plateau(i) = lin(i, 1);
        if opts.fitDrift
            o.drift(i) = lin(i, 2); o.amp(i) = lin(i, 3);
        else
            o.drift(i) = 0;         o.amp(i) = lin(i, 2);
        end
    end
    % Contrast: coherence amplitude at t = 0 over the plateau. Peak-to-PLATEAU, not
    % peak-to-peak -- see the note in OUT above. Magnitude, so the readout's sign convention
    % does not change it.
    o.contrast = abs(o.amp(1)) / max(abs(o.plateau(1)), eps);

    % ---- errors: Jacobian over the FULL parameter set, plus profiles ------
    [err, cov] = paramErrors(P, lin, t, Y, w, nLin, useRev, fitP, chi2);
    B = bounds(span, dt, fitP, useRev, opts.pFixed, max(t));
    prof = struct();
    prof.T2 = profileError(1, P, t, Y, w, nLin, useRev, B);
    prof.p  = NaN; if fitP; prof.p = profileError(2, P, t, Y, w, nLin, useRev, B); end
    prof.Trev = NaN;
    if useRev; prof.Trev = profileError(3, P, t, Y, w, nLin, useRev, B); end
    o.T2Err = pickErr(err.T2, prof.T2, 20*span);
    o.pErr  = pickErr(err.p,  prof.p,  4);
    o.revivalPeriodErr = NaN;
    if useRev; o.revivalPeriodErr = pickErr(err.Trev, prof.Trev, span); end
    o.contrastErr = err.contrast;
    o.profileErr = prof; o.cov = cov;

    % ---- derived, and the things that decide whether to believe it --------
    o.coherenceArea = o.T2 * gamma(1 + 1/max(o.p, 0.05));
    o.pMeaning = pMeaningOf(o.p, o.pErr, fitP);
    o.pFitted = fitP;

    o.chi2Reduced = chi2;
    resAll = resid(:); yAll = Y(:); eAll = E(:);
    o.r2 = 1 - sum(resAll.^2) / max(sum((yAll - mean(yAll)).^2), eps);
    o.r2Ceiling = 1 - sum(eAll.^2) / max(sum((yAll - mean(yAll)).^2), eps);
    o.r2Kept = o.r2; o.r2CeilingKept = o.r2Ceiling;
    if nDown > 0
        m = true(1, n); m(downIdx) = false;
        rk = resid(:, m); yk = Y(:, m); ek = E(:, m);
        o.r2Kept = 1 - sum(rk(:).^2) / max(sum((yk(:) - mean(yk(:))).^2), eps);
        o.r2CeilingKept = 1 - sum(ek(:).^2) / max(sum((yk(:) - mean(yk(:))).^2), eps);
    end

    % IS THERE A DECAY AT ALL. echoGateTest reports the same statistic for the same reason: a
    % T2 fitted to a trace that does not measurably fall is a number with no measurement under
    % it, and it will still come back with an error bar.
    o.decayErrorBars = (max(model(1,:)) - min(model(1,:))) / max(median(E(1,:)), eps);

    % THE SHORT-TIME DEFICIT. Weighted mean residual over the first quarter of the sweep, in
    % sigma of that mean. Negative = data below model early = the gate-leak signature.
    early = t <= (min(t) + 0.25*span);
    if sum(early) >= 3
        r1 = resid(1, early); e1 = E(1, early);
        wv = 1./e1.^2;
        mu = sum(wv.*r1)/sum(wv);
        se = sqrt(1/sum(wv));
        o.shortDeficit = mu;
        o.shortDeficitSigma = mu/max(se, eps);
    else
        o.shortDeficit = NaN; o.shortDeficitSigma = NaN;
    end

    o.T2overTwoT1 = NaN;
    if ~isempty(opts.T1us) && isfinite(opts.T1us) && opts.T1us > 0
        o.T2overTwoT1 = o.T2 / (2*opts.T1us);
    end

    o.ampRatio = NaN;
    if nT == 2
        o.ampRatio = o.amp(2) / (o.amp(1) + (o.amp(1) == 0));
    end
    o.residAutocorr = lag1(resid ./ max(E, eps));

    o.t = t; o.y = Y; o.e = E; o.model = model; o.resid = resid;
    o.envelope = envelopeOf(P, t, lin, nLin);
    o.weights = w; o.nDown = nDown; o.downIdx = downIdx; o.downT = t(downIdx);
    o.labels = {labels{:}};                                              %#ok<CCAT1>
    o.traces = nT; o.spanUs = span; o.tMin = tMin; o.tMax = tMax;
    o.nDropped = nDrop; o.nyquistMHz = out.nyquistMHz;
    o.monoT2 = monoFinal.P(1); o.monoP = monoFinal.P(2);
    o.ok = true; o.why = '';
    o.coeffError = [o.T2Err, o.contrastErr, o.pErr];

    out = o; out.opts = opts;
    if ~isempty(pullInfo); out.pull = pullInfo; end

    if opts.verbose;  report(out);  end
    if opts.plotIt;   drawIt(out);  end
end

% ======================================================================
% MODEL PIECES
% ======================================================================

function R = revivalComb(t, Trev, wFrac)
%REVIVALCOMB A comb of Gaussian revivals, normalised so R(0) = 1.
%
% Normalised at zero rather than at its own maximum so that the fitted amplitude `a` means the
% coherence at t = 0 and nothing else. Without that, `a` absorbs whatever the comb happens to
% be worth at the origin and the contrast becomes uninterpretable.
%
% K runs two periods past the data so the sum does not change discontinuously as Trev moves
% during the search -- a truncation that jumps by a whole Gaussian puts a step in the chi2
% surface, and Nelder-Mead walks straight into it.
    if ~isfinite(Trev) || Trev <= 0 || ~isfinite(wFrac) || wFrac <= 0
        R = ones(size(t)); return
    end
    w = max(wFrac*Trev, eps);
    K = floor(max(t)/Trev) + 2;
    R  = zeros(size(t));
    R0 = 0;
    for k = 0:K
        R  = R  + exp(-((t - k*Trev)/w).^2);
        R0 = R0 + exp(-((0 - k*Trev)/w).^2);
    end
    R = R / max(R0, eps);
end

function X = designAt(t, P, nLin, useRev)
% [1, (t), decay*comb]. Everything the model needs that is linear in it.
    D = exp(-(max(t, 0)/max(P(1), eps)).^max(P(2), 0.05));
    if useRev
        D = D .* revivalComb(t, P(3), P(4));
    end
    if nLin == 3
        X = [ones(numel(t),1), t(:), D(:)];
    else
        X = [ones(numel(t),1), D(:)];
    end
end

function [chi2, lin, ssr, ndofOut] = chiOf(P, t, Y, w, nLin, useRev)
% Profiled chi2: the linear part is solved EXACTLY for each trace, so the search is only over
% T2 and p (and Trev, wFrac when the comb is in). Same trick as ramseyFitRobust, same reason:
% no starting guess for the plateau or the amplitude, and no local minimum in them.
    X = designAt(t, P, nLin, useRev);
    nT = size(Y, 1);
    lin = zeros(nT, nLin);
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

function M = modelOf(P, t, lin, nLin, useRev)
    X = designAt(t, P, nLin, useRev);
    M = zeros(size(lin,1), numel(t));
    for i = 1:size(lin,1)
        M(i,:) = (X * lin(i,:)')';
    end
end

function Env = envelopeOf(P, t, lin, nLin)
% The MONOTONE envelope -- plateau plus amplitude times the stretched exponential, with the
% comb taken out. This is the curve T2 actually describes, and on a trace with revivals it is
% the one worth looking at: the comb is 13-C bookkeeping, the envelope is the coherence.
    D = exp(-(max(t, 0)/max(P(1), eps)).^max(P(2), 0.05));
    Env = zeros(size(lin,1), numel(t));
    for i = 1:size(lin,1)
        base = lin(i,1);
        if nLin == 3; base = base + lin(i,2)*t; end
        Env(i,:) = base + lin(i,end)*D;
    end
end

function r = lag1(R)
% Lag-1 autocorrelation, pooled across projections, each row detrended by its own mean.
    r = 0; num = 0; den = 0;
    for i = 1:size(R, 1)
        v = R(i,:) - mean(R(i,:));
        if numel(v) < 3; continue; end
        num = num + sum(v(1:end-1) .* v(2:end));
        den = den + sum(v.^2);
    end
    if den > 0; r = num/den; end
end

function s = pMeaningOf(p, pErr, fitP)
    if ~fitP
        s = 'held fixed -- this is an assumption, not a measurement'; return
    end
    if p < 0.55
        s = 'ON ITS LOWER BOUND -- the decay is faster than a single exponential can be'; return
    end
    if p > 3.95
        s = 'ON ITS UPPER BOUND -- suspect an unmodelled collapse rather than a real p'; return
    end
    if ~isfinite(pErr) || pErr <= 0
        s = 'no usable error bar, so it is not distinguishable from anything'; return
    end
    d1 = abs(p - 1)/pErr; d2 = abs(p - 2)/pErr; d3 = abs(p - 3)/pErr;
    if d1 < 2
        s = 'consistent with 1: a fast bath, or a relaxation-limited echo (check T2 against 2*T1)';
    elseif d2 < 2
        s = 'consistent with 2: slow Gaussian noise, the usual shape when something technical dominates';
    elseif d3 < 2
        s = 'consistent with 3: a dipolar nuclear spin bath, the textbook 13-C-limited echo';
    else
        s = sprintf('%.1f/%.1f/%.1f sigma from 1/2/3 -- none of the standard baths', d1, d2, d3);
    end
end

% ======================================================================
% SEARCH
% ======================================================================

function mono = fitMonotone(t, Y, w, span, opts, nLin, fitP, pFix)
%FITMONOTONE Grid then polish, over T2 and p only.
%
% The T2 grid is LOG spaced over two and a half decades around the span, because T2 is a scale
% and a linear grid wastes almost all its points at the long end where the model barely moves.
% p is gridded over the four values anyone would fix it at, which is enough of a seed for the
% polish to find the rest.
    T2Grid = logspace(log10(0.04*span), log10(8*span), max(opts.nGridT2, 8));
    if fitP; pGrid = [0.8 1 1.5 2 2.5 3 3.5]; else; pGrid = pFix; end
    best = struct('chi2', Inf, 'T2', T2Grid(1), 'p', pGrid(1));
    for iP = 1:numel(pGrid)
        for iT = 1:numel(T2Grid)
            c = chiOf([T2Grid(iT), pGrid(iP)], t, Y, w, nLin, false);
            if c < best.chi2
                best.chi2 = c; best.T2 = T2Grid(iT); best.p = pGrid(iP);
            end
        end
    end
    P0 = [best.T2, best.p];
    B = bounds(span, NaN, fitP, false, pFix);
    P = nelderMead(@(Q) chiOf(Q, t, Y, w, nLin, false), P0, B);
    mono = struct('P', P, 'grid', best);
end

function rev = fitRevival(t, Y, w, span, dt, opts, nLin, P0mono, fitP, pFix)
%FITREVIVAL Grid the revival PERIOD jointly with T2, then polish everything.
%
% The period has to be gridded rather than seeded from a transform, for the same reason
% ramseyFitRobust grids the detuning: a periodogram of an echo is dominated by the decay
% itself, so the largest spectral feature sits at DC and the comb is a ripple on the side of
% it. Scanning Trev and letting the profiled chi2 choose gets it right by construction and
% costs a fraction of a second.
%
% IT IS GRIDDED WITH T2, NOT AT THE MONOTONE T2. That was the first version and it does not
% work, for a reason worth writing down. On collapse-and-revival data the monotone fit reads
% the FIRST COLLAPSE as the end of the decay, so it comes back with a T2 of a few
% microseconds -- on the synthetic 60 us trace in testEchoFit it returned 2.6 us. At that T2
% the envelope is dead before the second point and NO comb improves anything, so every period
% in the grid scores the same and the search picks noise. Probing three spread T2 values
% alongside the period fixes it: T2 and Trev are nearly uncorrelated -- one sets how fast the
% envelope dies, the other where the modulation sits -- so a few T2 values are enough to find
% the right period, and T2 is then refined once at that period.
%
% The trig is hoisted the same way ramseyFitRobust hoists its quadratures: the comb depends
% only on (Trev, wFrac) and the envelope only on (T2, p), so the comb is built once per
% period-width pair instead of once per combination.
%
% The grid runs to the END of the sweep. Whether a period that long is admissible is decided by
% the acceptance test in the caller, on the fitted values, rather than approximated here by a
% ratio -- see the note at the top of the file for what that ratio cost.
    rev = [];
    loT = max(3*dt, 0.02*span);
    % UP TO THE END OF THE SWEEP, not half of it. The first version stopped at 0.55*span on the
    % reasoning that a shorter grid could not hold two periods -- which double-counted the
    % constraint, because revivalComb has a peak PINNED AT t = 0. A single revival away from
    % the origin is therefore a COMPLETE period, anchored at both ends, not half of one.
    %
    % That mistake cost a real measurement. A 60-point Hahn echo out to 120 us of evolution
    % came back fully collapsed by 20 us and then plainly revived around 105 us, both
    % projections moving in opposite directions by 3-4 error bars each across fifteen points.
    % The period it needed was 105 us and the grid stopped at 65, so no period in the grid
    % could describe it, and the fit reported T2 = 11.8 us -- the width of the first collapse
    % -- against a real T2 well over a hundred. Precisely the failure this comb was written to
    % prevent, produced by the guard meant to keep it honest.
    %
    % What actually has to hold is that the revival be INSIDE the data and SAMPLED, and both
    % are now checked on the fitted values by the acceptance test in the caller rather than
    % approximated by a period count here.
    hiT = 0.98*max(t);
    if ~isfinite(loT) || loT >= hiT; return; end

    if ~isempty(opts.revivalPeriodUs)
        Tgrid = opts.revivalPeriodUs;
    else
        Tgrid = linspace(loT, hiT, max(opts.nGridRev, 40));
    end
    if ~isempty(opts.revivalWidthFrac)
        wGrid = opts.revivalWidthFrac;
    else
        wGrid = [0.06 0.12 0.20 0.30];
    end
    T2Grid = logspace(log10(0.05*span), log10(10*span), max(opts.nGridT2, 8));
    probe  = T2Grid(unique(max(1, round([0.25 0.55 0.85]*numel(T2Grid)))));
    pUse   = P0mono(2); if ~fitP; pUse = pFix; end

    % ---- stage 1: period and width, at a few spread T2 -------------------
    best = struct('chi2', Inf, 'P', []);
    for iT = 1:numel(Tgrid)
        for iW = 1:numel(wGrid)
            R = revivalComb(t, Tgrid(iT), wGrid(iW));
            for iP = 1:numel(probe)
                c = chiFromComb(t, Y, w, R, probe(iP), pUse, nLin);
                if c < best.chi2
                    best.chi2 = c;
                    best.P = [probe(iP), pUse, Tgrid(iT), wGrid(iW)];
                end
            end
        end
    end
    if isempty(best.P); return; end

    % ---- stage 2: T2, at that period -------------------------------------
    R = revivalComb(t, best.P(3), best.P(4));
    for iT = 1:numel(T2Grid)
        c = chiFromComb(t, Y, w, R, T2Grid(iT), best.P(2), nLin);
        if c < best.chi2; best.chi2 = c; best.P(1) = T2Grid(iT); end
    end

    % ---- stage 3: polish all four ----------------------------------------
    B = revBounds(span, dt, fitP, pFix, opts, max(t));
    P = nelderMead(@(Q) chiOf(Q, t, Y, w, nLin, true), best.P, B);
    rev = struct('P', P, 'grid', best);
end

function chi2 = chiFromComb(t, Y, w, R, T2, p, nLin)
% Profiled chi2 from a precomputed comb. Same arithmetic as chiOf, minus rebuilding the comb.
    D = exp(-(max(t, 0)/max(T2, eps)).^max(p, 0.05)) .* R;
    if nLin == 3
        X = [ones(numel(t),1), t(:), D(:)];
    else
        X = [ones(numel(t),1), D(:)];
    end
    s = 0; ndof = 0;
    for i = 1:size(Y, 1)
        wi = w(i, :)';
        [th, nUsed] = wls(X, Y(i,:)', wi);
        r = Y(i,:)' - X*th;
        s = s + sum(wi .* r.^2);
        ndof = ndof + numel(t) - nUsed;
    end
    chi2 = s / max(ndof, 1);
end

function B = revBounds(span, dt, fitP, pFix, opts, tMax)
%REVBOUNDS bounds() with the caller's FIXED period and width pinned into the box.
%
% It exists because the overrides were being applied in fitRevival and then thrown away by the
% robust refit, which built its own plain bounds() and let Nelder-Mead walk the period off
% again. Handing revivalPeriodUs = 105 produced a fitted period of 70.8 -- the option looked
% like it did nothing. One function, called from both places, cannot drift like that.
    B = bounds(span, dt, fitP, true, pFix, tMax);
    if ~isempty(opts.revivalPeriodUs);  B(3,:) = opts.revivalPeriodUs  * [1 1]; end
    if ~isempty(opts.revivalWidthFrac); B(4,:) = opts.revivalWidthFrac * [1 1]; end
end

function B = bounds(span, dt, fitP, useRev, pFix, tMax)
% Box bounds for the nonlinear block, in the order [T2, p, (Trev, wFrac)].
%
% A FIXED p is expressed as a degenerate box [pFix pFix] rather than by shortening the
% parameter vector. That keeps one parameter layout everywhere -- the same P indexes the same
% quantity whether or not p is free -- and profileError already returns NaN for a parameter
% whose box has zero width, so "fixed" propagates into the error reporting for free.
    B = [0.02*span, 20*span; 0.5, 4];
    if ~fitP
        if ~exist('pFix', 'var') || isempty(pFix); pFix = 1; end
        B(2,:) = [pFix pFix];
    end
    if useRev
        lo = 0.02*span; if isfinite(dt); lo = max(3*dt, lo); end
        % The ceiling is the END OF THE SWEEP, not a fraction of it, and it is expressed against
        % max(t) rather than the span because the comb is anchored at t = 0: a period of 105 us
        % is perfectly measurable on data running from 2 to 120 us even though nothing repeats
        % within it. This box read 0.60*span in the first version and quietly capped a real
        % 105 us revival at 70.8 -- the fitted period sat exactly on the wall and the comb was
        % then refused for "no modulation", which is what a model pinned away from the answer
        % always reports.
        if ~exist('tMax', 'var') || isempty(tMax) || ~isfinite(tMax); tMax = 1.02*span; end
        B(3,:) = [lo, max(0.98*tMax, 1.5*lo)];
        % wFrac up to 0.5, where the Gaussians merge into a constant. The acceptance test
        % refuses anything past 0.4, but the SEARCH is allowed to walk into the merged region
        % -- a box that stops at 0.4 would pile probability on the wall and make a wide comb
        % look like a narrow one that happened to hit the edge.
        B(4,:) = [0.02, 0.5];
    end
end

function [th, nUsed] = wls(X, y, wi)
%WLS Weighted least squares that survives a degenerate column.
%
% Which happens here in a way it does not in the Ramsey. If T2 wanders far past the span
% during the search, exp(-(t/T2)^p) is 1 to working precision across every point, so the
% envelope column becomes a copy of the constant column and the two are collinear. The normal
% equations are then singular in that direction and the solve returns a huge amplitude and an
% equally huge opposite plateau -- a "fit" with a contrast of several hundred per cent that
% reproduces the data perfectly and means nothing.
%
% Columns whose weighted norm is negligible against the largest are DROPPED and their
% coefficients set to zero; nUsed reports how many were really fitted so the degrees of
% freedom stay honest. rcond alone does not catch it: with one near-duplicate column among
% three the matrix is still nominally invertible.
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

    for iter = 1:500
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

function [err, C] = paramErrors(P, lin, t, Y, w, nLin, useRev, fitP, chi2)
% Numerical Jacobian over the FULL parameter vector -- the nonlinear ones AND every linear one
% -- then C = inv(J'WJ), Birge scaled by the reduced chi2 where that exceeds 1.
%
% The linear part has to be in there. Profiling it out is a fitting convenience, not a licence
% to ignore it in the covariance, and on an echo the correlation is severe: T2 trades directly
% against the amplitude and against the plateau, because a smaller amplitude with a longer T2
% describes almost the same curve. An error bar computed over T2 and p alone comes out far too
% small -- which is exactly the shape of mistake that makes a fit look better than it is.
    err = struct('T2', NaN, 'p', NaN, 'Trev', NaN, 'wFrac', NaN, 'contrast', NaN);
    C = [];
    nT = size(Y, 1); nP = numel(P);
    theta = [P(:); reshape(lin', [], 1)];
    nTot = numel(theta);

    f = @(th) flatModel(th, nP, nT, nLin, t, useRev);
    y0 = f(theta);
    J = zeros(numel(y0), nTot);
    for k = 1:nTot
        h = 1e-6 * max(abs(theta(k)), 1e-6);
        tp = theta; tp(k) = tp(k) + h;
        J(:, k) = (f(tp) - y0) / h;
    end
    wv = reshape(w', [], 1);
    % DROP DEGENERATE DIRECTIONS, for the same reason wls does. Two occur here and both are
    % real rather than numerical: p is stationary when the decay is barely resolved, and with
    % pFixed the p column is identically zero by construction. Anything dropped gets a PROFILE
    % error instead, below, which is the right tool at a bound.
    sc = sqrt(sum((J.^2) .* wv, 1));
    good = sc > 1e-9 * max(sc);
    % A FIXED p is not a free direction, and leaving it in the inversion would hand its
    % variance to T2 -- which is the wrong sign of error: fixing a parameter should make the
    % others' error bars SMALLER, not larger. The model still depends on p, so its Jacobian
    % column is perfectly well conditioned and nothing else would have dropped it.
    if ~fitP && nP >= 2; good(2) = false; end
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
    err.T2 = s(1);
    if fitP; err.p = s(2); end
    if useRev && nP >= 4; err.Trev = s(3); err.wFrac = s(4); end

    % contrast = |a|/|plateau| for trace 1: propagate through that.
    i0 = nP;
    c0 = lin(1,1); a = lin(1,end);
    if abs(c0) > 0 && ~isempty(C)
        g = zeros(1, nTot);
        g(i0+1)    = -abs(a)/c0^2 * sign(c0);
        g(i0+nLin) =  sign(a)/abs(c0);
        if all(isfinite(C(:))) || any(isfinite(C(:)))
            Cg = C; Cg(~isfinite(Cg)) = 0;
            err.contrast = sqrt(abs(g * Cg * g'));
        end
    end
end

function e = pickErr(jac, prof, range)
% Take the larger, and fall back to whichever exists. An "error" bigger than half the search
% range is not an error bar, it is a degenerate direction, so it is discarded rather than
% printed as a reassuringly large number.
    cand = [];
    if isfinite(jac)  && jac  > 0 && jac  < 0.5*range; cand(end+1) = jac;  end
    if isfinite(prof) && prof > 0 && prof < 0.5*range; cand(end+1) = prof; end
    if isempty(cand); e = NaN; else; e = max(cand); end
end

function e = profileError(k, P, t, Y, w, nLin, useRev, B)
%PROFILEERROR 1-sigma from scanning one parameter: where the weighted SSR rises by 1.
%
% The linear part is re-solved at every step, so this is a genuine profile over the plateau,
% the drift and the amplitude, and only the OTHER nonlinear parameters are held. That makes it
% a lower bound where those are correlated, which is why the caller takes the larger of this
% and the Jacobian rather than replacing one with the other.
%
% ===================================================================
% EXPAND, THEN BISECT -- AND WHY THAT IS NOT A DETAIL
% ===================================================================
%
% The obvious way to write this is a geometric ladder: step out by 2% of a scale, multiply the
% step by 1.35, stop when the SSR has risen by 1, return that distance. It is what the first
% version of this function did, and it is WRONG IN ONE DIRECTION ONLY, which is the dangerous
% kind of wrong: the returned value is the first ladder rung PAST the crossing, so it can only
% ever be too big, never too small.
%
% Caught by testEchoFit's t7, which compares the quoted error bar against the actual
% scatter of T2 over 14 noise realisations. The ladder returned 0.9431 us -- for one
% projection and for two, and for every seed, because the first rung was already past the
% crossing and the "error bar" was simply 2% of the box width. The true scatter was 0.2559 us.
% A factor of 3.7 too LARGE, and constant, so nothing about it looked like noise.
%
% Too large is the safe direction and this was never going to make a fit look better than it
% is. But a constant error bar is not an error bar: it does not shrink when you take more
% averages, it does not shrink when you fit both projections instead of one, and it therefore
% destroys the two things the number is actually used for.
%
% So: expand geometrically only to BRACKET the crossing, then bisect the bracket to 0.01%.
% If the scan reaches its bound without the SSR ever rising by 1, that side is unbounded and
% contributes nothing rather than contributing the distance to the bound.
%
% NOTE FOR ramseyFitRobust AND rabiFitRobust: they carry the ladder version of this function.
% Same mechanism, same one-sided error, and their quoted parameter errors will be too large by
% however much their first rung overshoots. Worth porting this across.
    e = NaN;
    if k > numel(P); return; end
    if B(k,1) == B(k,2); return; end          % held fixed: no error to report
    [~, ~, ssr0] = chiOf(P, t, Y, w, nLin, useRev);
    if ~isfinite(ssr0); return; end
    lo = B(k,1); hi = B(k,2); p0 = P(k);

    for side = [1 -1]
        if side > 0; reach = hi - p0; else; reach = p0 - lo; end
        if ~(reach > 0); continue; end

        % ---- expand until the SSR has risen by 1, or the bound is reached -----
        step = max(1e-4*max(abs(p0), eps), 1e-7*(hi - lo));
        below = 0; above = NaN;
        for it = 1:90
            d = min(step, reach);
            Q = P; Q(k) = p0 + side*d;
            [~, ~, ssr] = chiOf(Q, t, Y, w, nLin, useRev);
            if ssr - ssr0 >= 1
                above = d; break
            end
            below = d;
            if d >= reach; break; end          % ran into the bound without crossing
            step = step * 1.6;
        end
        if ~isfinite(above); continue; end      % unbounded on this side

        % ---- bisect the bracket ----------------------------------------------
        a = below; b = above;
        for it = 1:50
            mid = 0.5*(a + b);
            Q = P; Q(k) = p0 + side*mid;
            [~, ~, ssr] = chiOf(Q, t, Y, w, nLin, useRev);
            if ssr - ssr0 >= 1; b = mid; else; a = mid; end
            if (b - a) <= 1e-4*max(b, eps); break; end
        end
        v = 0.5*(a + b);
        if v > 0
            if isnan(e); e = v; else; e = max(e, v); end
        end
    end
end

function v = flatModel(theta, nP, nT, nLin, t, useRev)
    P = theta(1:nP)';
    lin = reshape(theta(nP+1:end), nLin, nT)';
    M = modelOf(P, t, lin, nLin, useRev);
    v = reshape(M', [], 1);
end

% ======================================================================
% DATA EXTRACTION
% ======================================================================

function [t, Y, E, info] = pullEcho(Exp, opts)
%PULLECHO Get the time axis and both ratio traces out of an AWGExpEcho.
%
% Deliberately NOT extract1Ddata. That applies three successive cuts to the ratio across
% averages, at 5, then 3, then ONE standard deviation -- and a one-sigma cut discards about a
% third of perfectly good Gaussian data, after two harsher passes have already shrunk the
% spread it measures against. Every error bar downstream then comes out too small.
% ramseyFitRobust, rabiFitRobust and esrFitRobust all refused it and so does this.
%
% What happens instead: the reduced traces and their sterr come from signalParam and
% signalParam2, and where the per-average ratios are reachable the error is taken as the
% LARGER of the propagated error and the scatter across averages. Those answer different
% questions -- photon statistics within an average, and everything that drifts between them --
% and on this setup the second is routinely bigger. echoGateTest measured the sterr from a
% repeated Echo run as UNDERSTATED BY 4x, so this is not a theoretical concern here; taking
% the larger is the only choice that does not quietly understate the uncertainty.
    info = struct('nAvg', NaN, 'repeats', NaN, 'usedScatter', false, ...
                  'doubleMeasurement', false, 'cycles', NaN, 'pulseUs', [], ...
                  'tUnits', 'total evolution', 'droppedAverages', 0);

    cycles = tryGet(Exp, 'cycles', NaN);
    info.cycles = cycles;
    halfPi = tryGet(Exp, 'halfPiTime', []);
    piT    = tryGet(Exp, 'piTime', []);
    if ~isempty(halfPi) && ~isempty(piT) && isfinite(cycles)
        % The pulses inside one sequence: two pi/2 and `cycles` pi. Below this much total time
        % there is no free evolution left to decay, so tMin defaults here.
        info.pulseUs = 2*halfPi + cycles*piT;
    end

    tau = tryGet(Exp, 'tau', []);
    tEvo = [];
    ax = tryGet(Exp, 'mCurrentXAxisParam', []);
    if ~isempty(ax); tEvo = tryGet(ax, 'value', []); end
    if isempty(tEvo) && ~isempty(tau) && isfinite(cycles)
        tEvo = 2*tau*cycles;
    end

    switch lower(opts.xUnits)
        case 'tau'
            t = tau;
            if isempty(t); t = tEvo / max(2*cycles, 1); end
            info.tUnits = 'tau';
        otherwise
            t = tEvo;
            if isempty(t); t = tau; info.tUnits = 'tau (no evolution axis found)'; end
    end
    if isempty(t); error('echoFitRobust: the experiment has no tau or evolution axis.'); end
    t = t(:)';
    if strcmpi(info.tUnits, 'tau') && ~isempty(info.pulseUs)
        info.pulseUs = info.pulseUs / max(2*cycles, 1);
    end

    info.repeats = tryGet(Exp, 'repeats', NaN);
    info.nAvg    = tryGet(Exp, 'currIter', tryGet(Exp, 'averages', NaN));
    info.doubleMeasurement = isequal(tryGet(Exp, 'doubleMeasurement', false), true);

    [y1, e1] = paramData(Exp, 'signalParam');
    [y0, e0] = paramData(Exp, 'signalParam2');
    if isempty(y1); error('echoFitRobust: signalParam is empty -- did the run finish?'); end

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

    % Count whole averages that sit far from the median. A MAD rule on 5+ samples applied
    % ONCE -- note the contrast with the one-sigma cut above.
    if isfinite(opts.avgOutlierMad) && isfinite(nIter) && nIter >= 5 && ...
            ~isempty(S) && ndims(S) == 3
        try
            m = squeeze(mean(S(1,:,1:nIter) ./ max(S(2,:,1:nIter), eps), 2));
            md = median(m); mad0 = 1.4826*median(abs(m - md));
            if mad0 > 0
                info.droppedAverages = sum(abs(m - md) > opts.avgOutlierMad*mad0);
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
    fprintf('\n=== echoFitRobust ===\n');
    fprintf('model: stretched exponential%s, %d projection(s)%s, axis = %s time\n', ...
            eitherText(o.revivalsUsed, ' x 13-C revival comb', ''), o.traces, ...
            eitherText(o.traces > 1, ' fitted jointly', ''), o.tUnits);

    fprintf('\n  T2          %8.3f +/- %.3f us\n', o.T2, o.T2Err);
    if ~isfinite(o.T2Err)
        fprintf(['              *** NO ERROR BAR. Both the Jacobian and the profile scan came\n' ...
                 '              back either non-finite or wider than half the search range, and\n' ...
                 '              an "error" that big is a degenerate direction rather than an\n' ...
                 '              uncertainty, so it is discarded instead of printed. T2 above is\n' ...
                 '              a number the fit landed on, not a measurement. The usual causes,\n' ...
                 '              in order: the sweep does not reach past the decay (compare T2\n' ...
                 '              with the span below), fitDrift is on and competing with the\n' ...
                 '              decay, or there is no decay to measure.\n']);
    end
    fprintf('  stretch p   %8.3f +/- %.3f   %s\n', o.p, o.pErr, o.pMeaning);
    fprintf('  coh. area   %8.3f us        T2*gamma(1+1/p): compare THIS across runs when p moves\n', ...
            o.coherenceArea);
    fprintf('  contrast    %8.3f +/- %.3f %%   (peak-to-PLATEAU, not peak-to-peak)\n', ...
            100*o.contrast, 100*o.contrastErr);
    if isfinite(o.contrastErr) && o.contrastErr >= o.contrast
        fprintf(['              *** the error is as big as the value, so the coherence is NOT\n' ...
                 '              MEASURED and everything above describes noise. Raise repeats,\n' ...
                 '              and check the pulse calibration before blaming the sample.\n']);
    end

    % ---- is there a decay at all -----------------------------------------
    fprintf('\n  decay span  %8.1f error bars\n', o.decayErrorBars);
    if o.decayErrorBars < 5
        fprintf(['              *** the fitted decay barely rises above the noise. A T2 from a\n' ...
                 '              trace that does not measurably fall is a number with nothing\n' ...
                 '              under it -- and it still comes with an error bar. Extend tau\n' ...
                 '              past the decay, or raise repeats.\n']);
    end
    if o.T2 > 0.9*o.spanUs
        fprintf(['              *** T2 = %.1f us against a %.1f us sweep, so it is an\n' ...
                 '              EXTRAPOLATION and not a measurement. The plateau and the\n' ...
                 '              amplitude are trading freely. Extend the sweep past the decay.\n'], ...
                o.T2, o.spanUs);
    end

    % ---- revivals ---------------------------------------------------------
    if o.revivalsUsed
        fprintf('\n  revivals    period %.3f +/- %.3f us, width %.3f us (w/T = %.2f)\n', ...
                o.revivalPeriodUs, o.revivalPeriodErr, o.revivalWidthUs, o.revivalWidthFrac);
        fprintf(['              implies 13-C Larmor %.4f MHz, i.e. B = %.1f G. That is a CHECK,\n' ...
                 '              not a measurement: compare it against the coil setting. Half or\n' ...
                 '              twice the field you set means the total-time-to-per-arm factor\n' ...
                 '              (2*cycles = %d) is wrong for this sequence -- the PERIOD above is\n' ...
                 '              still what was measured, and T2 does not depend on the factor.\n'], ...
                o.larmorMHz, o.fieldGauss, 2*o.cycles);
        fprintf('              %d revival(s) in the sweep, %.1f points per period, %d on the last peak\n', ...
                o.revivalsVisible, o.pointsPerRevival, o.revivalTest.pointsOnPeak);
        if o.pointsPerRevival < 6
            fprintf(['              *** under 6 points per revival. The comb is barely sampled and\n' ...
                     '              the period may be aliased. Nyquist for this step is %.4f MHz.\n'], ...
                    o.nyquistMHz);
        end
        if o.revivalsVisible < 2
            % The single-revival case is worth its own paragraph, because the fit will hand back
            % a T2 with an error bar and the error bar is the only thing saying how thin the
            % measurement is. With one revival there is exactly ONE number carrying the
            % envelope -- the ratio of that revival's height to the coherence at t = 0 -- and
            % exp(-(Trev/T2)^p) has two unknowns in it. T2 and p are then degenerate along a
            % curve, and which point on that curve you land on is decided by noise.
            fprintf(['              *** ONE revival only. The envelope rests on a single number: the\n' ...
                     '              height of that revival relative to t = 0. exp(-(T/T2)^p) has TWO\n' ...
                     '              unknowns in it, so T2 and p are degenerate here and the split\n' ...
                     '              between them is set by noise, not by data. Two ways out, and the\n' ...
                     '              first is better:\n' ...
                     '                1. extend tau so a SECOND revival lands in range (double the\n' ...
                     '                   sweep: tau to %.0f us at cycles = %d). Two revivals measure\n' ...
                     '                   T2 and p separately.\n' ...
                     '                2. or set pFixed (3 for a 13-C spin bath, 1 for a fast one) and\n' ...
                     '                   quote T2 as conditional on that assumption. The report will\n' ...
                     '                   then say p was held rather than measured.\n'], ...
                    2.05*o.revivalPeriodUs/max(2*o.cycles,1), o.cycles);
        end
        fprintf(['              WITHOUT the comb this fit returned T2 = %.3f us. The monotone fit\n' ...
                 '              has to explain the first collapse as the end of the decay, which\n' ...
                 '              is why it comes out short -- that gap is the whole reason the\n' ...
                 '              comb is in the model.\n'], o.monoT2);
        if ~isempty(o.revivalTest.failed)
            fprintf('              NOTE: %s\n', o.revivalTest.failed);
        end
    elseif o.revivalTest.tried
        fprintf('\n  revivals    NOT used.\n');
        if isempty(o.revivalTest.failed)
            fprintf('              (no reason recorded -- that is a bug, please report it)\n');
        else
            fprintf('              %s\n', o.revivalTest.failed);
        end
        fprintf(['              A monotone echo is the normal result at low field or in isotopically\n' ...
                 '              pure diamond. If you EXPECT collapses, the sweep has to reach past\n' ...
                 '              ONE full Larmor period -- the comb is anchored at t = 0, so a single\n' ...
                 '              revival is enough to fit a period. And look at the two diagnostics\n' ...
                 '              above before believing this line: a residual autocorrelation past\n' ...
                 '              0.4 or an R^2 well under its ceiling means something IS there.\n']);
    else
        fprintf('\n  revivals    switched off by the caller (revivals = ''off'').\n');
    end

    % ---- the short-time deficit -------------------------------------------
    if isfinite(o.shortDeficitSigma)
        fprintf('\n  short-t     %+.2f sigma over the first quarter of the sweep\n', o.shortDeficitSigma);
        if o.shortDeficitSigma < -3
            fprintf(['              *** the data sits BELOW the model early. That is the signature\n' ...
                     '              of contrast lost at short evolution times, which is\n' ...
                     '              arithmetically indistinguishable from a shorter T2 -- so this\n' ...
                     '              number, not the fit, is what tells you about it.\n' ...
                     '              echoGateTest has measured the fixed microwave gate as harmless\n' ...
                     '              out to 240 us of excess at the 17%% level, so the gate is the\n' ...
                     '              LEAST likely cause. Check, in this order: the pi and pi/2\n' ...
                     '              calibration, a fast second decay component the single\n' ...
                     '              stretched exponential cannot reach, and tMin.\n']);
        elseif o.shortDeficitSigma > 3
            fprintf(['              the data sits ABOVE the model early, i.e. the model decays too\n' ...
                     '              fast at the start. Usually a p the fit could not reach, or an\n' ...
                     '              unmodelled collapse -- try revivals = ''on''.\n']);
        end
    end

    % ---- T1 ceiling -------------------------------------------------------
    if isfinite(o.T2overTwoT1)
        fprintf('\n  T2/(2*T1)   %8.3f\n', o.T2overTwoT1);
        if o.T2overTwoT1 > 0.8
            fprintf(['              *** at or above the 2*T1 ceiling. This echo is not measuring\n' ...
                     '              dephasing, it is measuring relaxation -- a longer T2 is not\n' ...
                     '              available without a better T1.\n']);
        end
    end

    % ---- fit quality ------------------------------------------------------
    fprintf('\n  chi2 (red)  %8.2f\n', o.chi2Reduced);
    fprintf('  R^2         %8.4f   (noise ceiling %.4f)\n', o.r2, o.r2Ceiling);
    if o.nDown > 0
        fprintf('  R^2 kept    %8.4f   (ceiling %.4f) -- excluding the down-weighted points\n', ...
                o.r2Kept, o.r2CeilingKept);
    end
    rr = o.r2; rc = o.r2Ceiling;
    if o.nDown > 0; rr = o.r2Kept; rc = o.r2CeilingKept; end
    if rc > 0 && rr < rc - 0.15
        fprintf('              well below the ceiling: the model is missing something real.\n');
    end
    fprintf('  resid acf   %+8.2f     lag-1 autocorrelation of the residual\n', o.residAutocorr);
    if abs(o.residAutocorr) > 0.4 && o.chi2Reduced > 0.05
        fprintf(['              *** STRUCTURED. The residual walks rather than scatters, so the\n' ...
                 '              model is missing something and the errors above are optimistic.\n' ...
                 '              On an echo the usual causes, in order: an unmodelled collapse\n' ...
                 '              (try revivals = ''on''), a decay that is not a single stretched\n' ...
                 '              exponential, and a drifting laser (fitDrift = true -- but read\n' ...
                 '              the degeneracy note in the header first).\n']);
    end
    if o.nDown > 0
        fprintf('\n  robust      %d point(s) down-weighted, at t = %s us\n', ...
                o.nDown, numList(o.downT, '%.3f'));
        fprintf(['              on this setup those are usually the classic path re-arming: a\n' ...
                 '              point whose repeats fired before the AWG resumed is diluted\n' ...
                 '              toward no contrast, which on an echo looks exactly like an\n' ...
                 '              honest early loss of coherence. Walk mode does not produce them.\n']);
        if o.nDown > 0.2*numel(o.t)
            fprintf('              *** that is over a fifth of the sweep. Look at the data, not the fit.\n');
        end
    end
    if o.traces == 2
        fprintf('\n  |0>/|1> amplitude ratio %+.3f   (should be -1)\n', o.ampRatio);
        if abs(o.ampRatio + 1) > 0.3
            fprintf(['              *** the two projections are not complementary. Before quoting\n' ...
                     '              T2, check that rows 3-4 really are the |0> readout -- the two\n' ...
                     '              windows may not be doing what the class documents.\n']);
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
% ONE axes by default. showResid adds a panel below, and it is off unless asked for -- the main
% axes answers the question and every extra panel shrinks it. The look comes from fitStyle, so
% this figure and the Ramsey, Rabi and ESR ones cannot drift apart.
    S = fitStyle();
    nExtra = double(o.opts.showResid);
    fig = S.newFigure('Echo', nExtra > 0);
    if nExtra == 0
        axM = S.axes(fig);
    else
        axM = S.axes(subplot(3, 1, 1:2));
    end

    td = linspace(min(o.t), max(o.t), 2000);
    lin = linFromModel(o);
    nLin = 2 + double(o.opts.fitDrift);
    P = paramVec(o);
    Md = modelOf(P, td, lin, nLin, o.revivalsUsed);
    Ed = envelopeOf(P, td, lin, nLin);

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

    % THE PLATEAU, always -- it is where the echo ends up once the coherence is gone, it is
    % FITTED rather than assumed to be 1, and the distance from it to the data is the contrast.
    % A horizontal dashed line is the cheapest way to make that distance readable.
    %
    % THE ENVELOPE, only when the comb is in the model. Without a comb the envelope IS the fit
    % and drawing it puts a dashed line exactly underneath the ink one, which is clutter with a
    % misleading legend entry. With a comb the ink curve dives into every collapse and the
    % envelope becomes the most informative line on the figure: it is the curve T2 actually
    % describes, and it passes through the revival extrema.
    %
    % Both suppressed when the coherence is not measured (contrastErr as big as contrast),
    % because then the amplitude is whatever the noise allowed and an "envelope" is a
    % confident-looking decay drawn from a number that means nothing. A clipped fiction is
    % still a fiction.
    if isfinite(o.contrastErr) && o.contrastErr < 0.5*o.contrast
        plat = lin(1,1) + eitherNum(o.opts.fitDrift, lin(1,2), 0)*td;
        % Name WHOSE plateau when there are two of them. Each projection has its own -- they
        % are properties of a readout window, not of the spin -- and an unlabelled dashed line
        % sitting between two traces invites the reading that it belongs to both.
        who = '';
        if o.traces > 1; who = sprintf(' (%s)', o.labels{1}); end
        if o.revivalsUsed
            H(end+1) = S.envelope(axM, td, Ed(1,:), plat);
            L{end+1} = sprintf('envelope and plateau%s, T_2 = %.3g \\mus', who, o.T2);
        else
            H(end+1) = S.envelope(axM, td, plat);
            L{end+1} = sprintf('plateau%s = %.4f', who, lin(1,1));
        end
    end
    % Revival positions, so the period can be checked by eye against the data rather than
    % taken on trust from the number in the title.
    if o.revivalsUsed && isfinite(o.revivalPeriodUs) && o.revivalPeriodUs > 0
        yl = ylim(axM);
        kMax = floor(max(o.t)/o.revivalPeriodUs);
        for k = 1:kMax
            xk = k*o.revivalPeriodUs;
            hk = plot(axM, [xk xk], yl, '--', 'Color', S.col.flag, 'LineWidth', 1.0);
            if k == 1
                H(end+1) = hk;                                             %#ok<AGROW>
                L{end+1} = sprintf('revivals, %.3g \\mus apart', o.revivalPeriodUs); %#ok<AGROW>
            end
        end
    end
    if o.nDown > 0
        H(end+1) = S.flag(axM, o.downT, o.y(1, o.downIdx));
        L{end+1} = sprintf('%d down-weighted', o.nDown);
    end

    xlabel(axM, xLabelOf(o), 'Color', S.col.ink2);
    ylabel(axM, yLabelOf(o, S), 'Color', S.col.ink2);
    % Scale to the MODEL and the KEPT points, not to everything -- see fitStyle/focusY. On an
    % echo this matters more than on a Ramsey: a point diluted to the plateau sits inside the
    % data range rather than outside it, so it does not stretch the axis, but the two or three
    % that do sit outside would flatten a 1% decay into a band.
    keepMask = true(1, numel(o.t)); keepMask(o.downIdx) = false;
    S.focus(axM, [Md(:); Ed(1,:)'; reshape(o.y(:, keepMask), [], 1)]);
    S.tight(axM, o.t); S.headroom(axM, 0.30);
    S.legend(axM, H, L);
    S.title(axM, titleLines(o));

    if o.opts.showResid
        axR = S.axes(subplot(3, 1, 3));
        S.zeroLine(axR, [min(o.t) max(o.t)]);
        Rn = o.resid ./ max(o.e, eps);
        for i = 1:o.traces
            S.resid(axR, o.t, Rn(i,:), S.slots{min(i,3)});
        end
        S.focus(axR, reshape(Rn(:, keepMask), [], 1), 0.15);
        xlabel(axR, xLabelOf(o), 'Color', S.col.ink2);
        ylabel(axR, 'residual / \sigma', 'Color', S.col.ink2);
        S.tight(axR, o.t);
        S.title(axR, sprintf(['\\chi^2_\\nu = %.2f     lag-1 autocorr %+.2f%s     ' ...
                              'short-t %+.2f\\sigma'], ...
                             o.chi2Reduced, o.residAutocorr, ...
                             eitherText(abs(o.residAutocorr) > 0.4, ' (structured)', ''), ...
                             o.shortDeficitSigma));
    end
end

function v = eitherNum(tf, a, b)
    if tf; v = a; else; v = b; end
end

function s = xLabelOf(o)
    if strncmpi(o.tUnits, 'tau', 3)
        s = 'free evolution \tau  [\mus]';
    else
        s = 'total evolution time  2\tau\timescycles  [\mus]';
    end
end

function s = yLabelOf(o, S)
    s = S.yLabel;
    if isfield(o.opts, 'yLabel') && ~isempty(o.opts.yLabel); s = o.opts.yLabel; end
end

function T = titleLines(o)
%TITLELINES Three lines, grouped by the job each number does.
%
% One long line was the first version and it ran off the right-hand edge: nine items with the
% revival period and the down-weight count in it is about 180 characters, and MATLAB does not
% wrap a title, it clips it. The count -- the item most likely to change how you read the
% figure -- was the one that fell off.
%
% So they are grouped rather than shortened, because none of them is worth dropping:
%
%   line 1  THE ANSWER            T2, p, contrast
%   line 2  WHAT KIND OF DECAY    revivals and the field they imply, coherence area,
%                                 how many error bars the decay spans, T2 against 2*T1
%   line 3  WHETHER TO BELIEVE IT chi2, R^2 against its ceiling, residual autocorrelation,
%                                 the short-time deficit, the projection ratio, down-weights
%
% Reading down the three lines is then reading the fit in the order you should: what it says,
% what shape it is, and whether any of it holds up.
    L1 = sprintf('T_2 = %.4g \\pm %.3g \\mus     p = %.2f \\pm %.2f     contrast = %.3f \\pm %.3f %%', ...
                 o.T2, o.T2Err, o.p, o.pErr, 100*o.contrast, 100*o.contrastErr);

    b2 = {};
    if o.revivalsUsed
        b2{end+1} = sprintf('revivals %.4g \\mus \\rightarrow %.1f G', ...
                            o.revivalPeriodUs, o.fieldGauss);
    elseif o.revivalTest.tried
        b2{end+1} = 'no revivals';
    end
    b2{end+1} = sprintf('coherence area %.3g \\mus', o.coherenceArea);
    b2{end+1} = sprintf('decay spans %.0f error bars', o.decayErrorBars);
    if isfinite(o.T2overTwoT1)
        b2{end+1} = sprintf('T_2/2T_1 = %.2f', o.T2overTwoT1);
    end

    b3 = {sprintf('\\chi^2_\\nu = %.2f', o.chi2Reduced), ...
          sprintf('R^2 = %.3f of %.3f', rKept(o), rCeil(o)), ...
          sprintf('resid acf %+.2f', o.residAutocorr)};
    if isfinite(o.shortDeficitSigma)
        b3{end+1} = sprintf('short-t %+.1f\\sigma', o.shortDeficitSigma);
    end
    if o.traces == 2
        b3{end+1} = sprintf('|0\\rangle/|1\\rangle = %+.2f', o.ampRatio);
    end
    if o.nDown > 0
        b3{end+1} = sprintf('%d down-weighted (may be off scale)', o.nDown);
    end

    T = {L1, strjoin(b2, '     |     '), strjoin(b3, '     |     ')};
end

function r = rKept(o)
    r = o.r2; if o.nDown > 0; r = o.r2Kept; end
end

function r = rCeil(o)
    r = o.r2Ceiling; if o.nDown > 0; r = o.r2CeilingKept; end
end

function lin = linFromModel(o)
    nLin = 2 + double(o.opts.fitDrift);
    lin = zeros(o.traces, nLin);
    for i = 1:o.traces
        lin(i,1) = o.plateau(i);
        if nLin == 3; lin(i,2) = o.drift(i); end
        lin(i,end) = o.amp(i);
    end
end

function P = paramVec(o)
    P = [o.T2, o.p];
    if o.revivalsUsed; P = [P, o.revivalPeriodUs, o.revivalWidthFrac]; end
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
