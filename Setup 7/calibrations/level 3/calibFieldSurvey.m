function out = calibFieldSurvey(opts)
%CALIBFIELDSURVEY One wide ESR at the current field, and everything it can tell you.
%
%   out = calibFieldSurvey
%   out = calibFieldSurvey(struct('freqs', 2780:2:2960, 'targetMHz', 2820))
%
% LEVEL 3, first step. Measures the field you already have rather than changing it, and
% answers in one sweep: how well aligned you are, what |B| is in Gauss, whether the fit
% is being biased, and what coil setting would put a chosen line at a chosen frequency.
%
% Changes NO field setting. It reads the field that is currently applied and reports
% what to do. That is deliberate -- setting the field is a feedback loop and belongs in
% its own script, and this one has to be safe to run first.
%
% ===================================================================
% WHAT A <111> FIELD LOOKS LIKE, AND WHY THAT IS THE ONE TO WANT
% ===================================================================
%
% The four NV axes are the <111> directions, 109.47 degrees apart, cos = -1/3 exactly.
% Put B along one of them and:
%
%   the ALIGNED orientation      splits by 2*gamma*B,   depth 1 unit
%   the other THREE, degenerate  split by 2*gamma*B/3,  depth 3 units
%
% So four lines, and the ratio of the splittings is exactly 3.000. That is the target,
% and it is the right target for THIS setup because Rabi and Ramsey need a single clean
% orientation. A <100> field would make all four degenerate -- two lines at 4x the
% signal, tempting for ESR -- but then all four orientations are driven at four
% different Rabi rates and the Rabi beats.
%
% ===================================================================
% ALIGN ON THE SHARP FEATURE, NOT THE FLAT ONE
% ===================================================================
%
% The aligned splitting goes as cos(tilt), which is SECOND order: at 5 degrees off it
% has changed by 0.4%, about 0.4 MHz on a 100 MHz splitting, comparable to fit
% precision. Maximising it is a flat optimum.
%
% The three-fold line is degenerate ONLY at perfect alignment, and its components spread
% FIRST order in the tilt. Geometry: the three off-axis <111> directions have
% perpendicular components of magnitude sin(109.47) = 0.943 at 120 degrees of azimuth to
% each other, so a tilt theta spreads their splittings by 1.5*0.943*2*gamma*B*sin(theta)
% = 1.414 * s_aligned * sin(theta). At s_aligned = 100 MHz a 2 degree tilt spreads them
% by 3.5 MHz against a ~9 MHz linewidth -- visible. So the alignment metric is the
% EXCESS WIDTH of the three-fold line, and it gets more sensitive as B grows: align at a
% large field, work at a small one.
%
% This is the same lesson as the microwave phase calibration, where the antiphase null
% gave +/-0.4 degrees and the in-phase peak gave +/-8.5.
%
% ===================================================================
% THE PERPENDICULAR FIELD COMPONENT: WHERE IT MATTERS AND WHERE IT CANCELS
% ===================================================================
%
% For an orientation with B_par along its axis and B_perp across it, second-order
% perturbation on H = D*Sz^2 + gamma*(B_par*Sz + B_perp*Sx) gives, exactly to that order:
%
%   f+ - f-     = 2*g*B_par * [ 1 - (g*B_perp)^2 / (2*(D^2 - (g*B_par)^2)) ]
%   (f+ + f-)/2 = D + 1.5*(g*B_perp)^2 * D / (D^2 - (g*B_par)^2)
%
% Read those two lines together, because they say opposite things.
%
% THE MEAN IS SHIFTED, A LOT. Both transitions move UP by the same amount, so the pair
% mean rises. At 45 G the three off-axis orientations have B_perp = 41.2 G and their means
% sit 6.95 MHz above D, while the aligned orientation has B_perp = 0 and its mean IS D.
% That is why the pair means are not expected to agree, why the physical D is taken from
% the LEAST-shifted pair, and why the second-order check below predicts each mean rather
% than demanding they match.
%
% THE SPLITTING ALMOST CANCELS IT. The shift is common to both transitions, so it drops
% out of the difference and only the bracketed correction survives. At 45 G that is
% -0.081% on the off-axis splittings: 0.066 MHz on 81.5 MHz. Propagated through
% SUM p^2 = (4/3) B^2 it moves |B| by 0.009 G and the tilt by under 0.001 deg. Measured
% splitting errors are 0.3-0.8 MHz, so the systematic is an order of magnitude below the
% noise and is deliberately NOT corrected.
%
% This is the reason the metric is built from splittings rather than from line positions:
% splittings are immune to the largest systematic in the problem.
%
% ===================================================================
% THE BIAS THE PLAN WARNED ABOUT, AND WHAT REPLACES ITS FIX
% ===================================================================
%
% Section 3.4: `amplitude` is a scalar for a whole sweep, so B1 varies across it by
% whatever the resonator does, and an asymmetric envelope under a symmetric line pulls
% the fitted centre -- straight into the alignment numbers. The plan's fix was
% per-frequency amplitude flattening, and calibResonator showed that is NOT AVAILABLE:
% no reflection standard, so the table came out with |Gamma| > 1 and +156 dB entries.
%
% What replaces it: esrMultiLineFit with baselineOrder 2, so the envelope's curvature
% across a 200 MHz sweep is absorbed by the baseline rather than by the line centres,
% plus the pair-mean check above as an independent readout of whatever bias is left.
%
% ===================================================================
% OPTIONS (struct fields, all optional)
%   freqs        MHz (default 2780:2:2960 -- wide enough for the four-line <111>
%                pattern at 15-20 G, stepped at 2 MHz for >4 points per linewidth)
%   amplitude    dBm per channel (default [] = leave the json's value). -15 is what
%                worked at 45 G: it deepens the weak outer lines, and the inner ones are
%                already saturated so they barely change.
%   repeats      (default 50)
%   averages     (default 3 -- NOT 1. With a single average sterr is not a standard
%                error at all, which is how a 483-sigma dip got reported)
%   detectionMHz CW detection and reference window, us. Default [] = resolve it through
%                esrWindowPick, which is what you want. It was 10, justified in this
%                comment as "already many polarisation times" -- and that was simply
%                wrong: tau_pol is 13.6 us on this setup, so a 10 us window is less than
%                ONE polarisation time. Laser and microwave both start at the top of every
%                CW window, so its first ~tau_pol carries less than steady-state contrast;
%                at 10 us the dip came back about a third of its available depth. The
%                floor is ~10*tau_pol, and there is a hard AWG-trigger limit underneath it
%                as well -- see esrWindowPick.
%   jsonLocation expProperties.json, loaded with loadJsonParams('all', ...) before
%                anything else is set
%   isTracking   (default 0)
%   freshExperiment  construct a NEW AWGExpESR for every run (default true). See makeExp
%                -- this is what stops the "parameters have changed" dialog appearing.
%   newOrderEveryAverage  walk mode: draw a new random frequency order before every
%                average (true/false). Default [] = leave AWGExpESR's own setting (false
%                on a fresh object unless the json sets it). Each average then pays one
%                sequence upload; see AWGExpESR.newOrderEveryAverage.
%   targetMHz    the frequency you want a line moved to (default 2820, the
%                reflection minimum -- treat it as a hypothesis, see below)
%   whichLine    'aligned' or 'triple': which pair to move onto targetMHz
%                (default 'aligned')
%   baselineOrder passed to esrMultiLineFit (default 2)
%   plotIt       (default true)
%   D0, dTol     the zero-field splitting is known: D = 2870 MHz to within a few MHz
%                of temperature drift. Any pairing whose mean lands further than dTol
%                (default 15 MHz) from D0 is WRONG, and this refuses it rather than
%                reporting a field. The first run of this script did not have that
%                guard: with the field OFF it paired a real 2871 line against a broad
%                artefact at 2825 and reported |B| = 8.26 G and D = 2848.
%   pairTol      MHz. A line is considered to already have a partner if one sits within
%                this of the mirror position 2*D - f (default 3, and D is estimated
%                from the outermost pair rather than assumed). Otherwise the mirror
%                position is seeded into the fit as a missing line.
%   minDepthSigma passed to esrMultiLineFit for the seeded refits (default 2.5)
%   splitRatio   try splitting a line wider than this times the narrowest (default 1.8).
%                Two lines of one orientation share an intrinsic width, so a line much
%                wider than the narrowest is either a blend or envelope distortion.
%   splitGain    accept a split only if reduced chi2 falls to this fraction of its
%                previous value (default 0.9 -- was 0.75, which silently rejected the
%                genuine three-way split of a 43 MHz blend at 45 G). Adding parameters always improves a fit
%                slightly; requiring a real drop stops it inventing lines.
%   maxChi2      above this reduced chi2 the model is wrong, not merely noisy, and
%                centres are not trustworthy (default 10)
%   seedCentres  MHz, start the fit from these positions instead of finding them
%                (default [] = find them). Used by calibFieldAlign to carry each scan
%                point's lines into the next, where they have moved only a few MHz.
%   widthGain    accept a refit whose WIDEST line falls to this fraction of its previous
%                FWHM, provided chi2 does not degrade (default 0.75). This is the
%                acceptance rule that resolves a blend: splitting a 44 MHz blob into two
%                12 MHz lines barely moves reduced chi2, because the blob already absorbs
%                most of the deficit -- so a chi2-only test rejects a split that the
%                widths make obvious. No NV line here is wider than about 15 MHz.
%   verboseRefine print every mirror and split attempt and why it was taken or refused
%                (default true). Silence is how a rejected refit looked like an
%                unattempted one.
%   replay       a previous result struct (or one with freqs/signal/sterr, or an N-by-3
%                matrix): re-run the ENTIRE analysis on a sweep already taken, with no
%                experiment constructed and no hardware touched.
%
%                    sv = calibFieldSurvey(opts);            % once, on the setup
%                    save('field45.mat', 'sv')
%                    calibFieldSurvey(struct('replay', sv)); % as often as you like
%
%                Every fitter change so far has cost a live run on the magnet to test.
%                It should not.
%
% Writes nothing, changes no field setting.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('freqs', 2780:2:2960, 'amplitude', [], ...
        'repeats', 50, 'averages', 3, 'detectionMHz', [], 'targetMHz', 2820, ...
        'whichLine', 'aligned', 'baselineOrder', 2, 'plotIt', true, ...
        'D0', 2870, 'dTol', 15, 'maxChi2', 10, 'pairTol', 3, 'minDepthSigma', 2.5, 'splitRatio', 1.8, 'splitGain', 0.9, 'seedCentres', [], ...
        'jsonLocation', '', ...
        'isTracking', 0, 'freshExperiment', true, 'newOrderEveryAverage', [], ...
        'replay', [], 'verboseRefine', true, 'widthGain', 0.75, ...
        'physicsSeed', true, 'huberMHz', 6));

    % ---- where the calibration lives ------------------------------------
    % Resolved by setup7Paths rather than written down here -- see setup7Paths for why a
    % path repeated in twelve files is worse than N times the work to change.
    %
    % TOLERANT ON PURPOSE. This function can do useful work with no json at all: a replay
    % or dry run never reads it, and the offline graders run with no Setup 7 tree on the
    % path. So a failure to resolve leaves jsonLocation EMPTY, which every use site below
    % already guards with ~isempty, and the calibration falls back to class defaults. An
    % eager error here broke testFieldPipeline, which fits saved spectra and touches no
    % hardware and no json.
    if isempty(opts.jsonLocation)
        try
            opts.jsonLocation = setup7Paths().json;
        catch
            opts.jsonLocation = '';
        end
    end

    GAMMA = 2.8025;         % MHz/G, electron gyromagnetic ratio for NV (g = 2.0028).
                            % A physical constant, not a fitted quantity -- which is what
                            % makes the splitting an ABSOLUTE field measurement and lets
                            % it calibrate the coil setting rather than the reverse.

    fprintf('\n=== calibFieldSurvey ===\n');
    fprintf('%d points, %.0f to %.0f MHz, step %.1f MHz\n', numel(opts.freqs), ...
            min(opts.freqs), max(opts.freqs), median(diff(opts.freqs)));

    out = struct('opts', opts, 'gamma', GAMMA);
    a = [];
    saved = struct('frequency', [], 'repeats', [], 'averages', [], 'amplitude', []);

    try
        if isempty(opts.replay)
            a = makeExp(opts);
            if ~opts.freshExperiment
                saved = struct('frequency', a.frequency, 'repeats', a.repeats, ...
                               'averages', a.averages, 'amplitude', a.amplitude);
            end

            fprintf('\nrunning...\n');
            a.run;
            v = a.signalParam.value(:);
            e = a.signalParam.sterr(:);
            f = a.frequency(:);
        else
            % REPLAY: re-run the whole analysis on a sweep already taken, with no
            % hardware touched at all. Every iteration of the fitter to date has cost a
            % live run on the setup, which is both slow and a waste of the magnet's time
            % when the question being asked -- "why did the blend not split?" -- is
            % answered entirely by the numbers already on disk.
            %
            %   sv = calibFieldSurvey(opts);          % once, on the hardware
            %   save('field45.mat', 'sv')
            %   sv2 = calibFieldSurvey(struct('replay', sv, ...));   % as often as you like
            %
            % Accepts the output struct of a previous call, or any struct with freqs /
            % signal / sterr, or an N-by-2 or N-by-3 matrix.
            [f, v, e] = unpackReplay(opts.replay);
            fprintf('\nREPLAY: %d points from a previous sweep, %.0f to %.0f MHz.\n', ...
                    numel(f), min(f), max(f));
            fprintf('        No experiment was constructed and no hardware was touched.\n');
        end
        out.freqs = f; out.signal = v; out.sterr = e;

        % A sanity check on the data before anything is inferred from it. Point-to-point
        % scatter against the claimed error bars: if they disagree by a lot, every
        % uncertainty below inherits it. esrMultiLineFit's Birge scaling corrects the
        % SCALE, but a factor this large is worth seeing.
        scatter = 1.4826*median(abs(diff(v)))/sqrt(2);
        % median(e(isfinite(e))) rather than median(e,'omitnan') so the same file runs
%         under Octave, which is where the fitter is tested offline against synthetic
%         spectra before it is ever pointed at the magnet.
        eFin = e(isfinite(e));
        if isempty(eFin); eFin = NaN; end
        out.scatterRatio = scatter / max(median(eFin), eps);
        fprintf('  contrast %.3f%%, point-to-point scatter %.2fx the error bars\n', ...
                100*(1 - min(v)), out.scatterRatio);
        if out.scatterRatio > 3 || out.scatterRatio < 0.33
            fprintf(['  NOTE the error bars and the real scatter disagree by %.1fx. The fit\n' ...
                     '       rescales parameter errors by the residuals, so the numbers below\n' ...
                     '       are still usable, but averages = 1 is the usual cause.\n'], ...
                    max(out.scatterRatio, 1/out.scatterRatio));
        end

        % ---- fit every line ------------------------------------------------
        % seedCentres, when given, starts the fit from known line positions instead of
        % searching for them. During an angle scan the previous point's lines are only a
        % few MHz away, so every fit becomes a small refinement rather than a fresh search
        % from a peak-finder that can and does fail on a marginal line.
        fit = esrMultiLineFit(f, v, e, struct('baselineOrder', opts.baselineOrder, ...
                    'plot', opts.plotIt, 'name', 'field survey', ...
                    'seedCentres', opts.seedCentres, 'maxLines', 10));
        out.fit = fit;
        if ~fit.ok
            fprintf('\nFIT FAILED: %s\n', fit.why);
            restore(a, saved); return
        end
        c  = [fit.lines.centre];  ce = [fit.lines.centreErr];
        w  = [fit.lines.fwhm];    we = [fit.lines.fwhmErr];
        dp = [fit.lines.depth];
        n  = fit.nLines;

        % ---- is the MODEL right, before trusting any number from it ----------
        % A reduced chi2 of 20 does not mean noisy data, it means the model is the wrong
        % shape -- and Birge scaling then inflates the error bars to cover the model
        % error, which makes them look honest while the CENTRES stay biased.
        if isfinite(fit.chi2Reduced) && fit.chi2Reduced > opts.maxChi2
            fprintf(['\n  reduced chi2 = %.1f, so the model is imperfect (errors scaled %.2fx\n' ...
                     '  to cover it). Usually the resonator envelope, or a line the finder\n' ...
                     '  merged with its neighbour. The mirror-seeding below is the fix for the\n' ...
                     '  second; narrow-window refits are the fix for the first.\n'], ...
                    fit.chi2Reduced, fit.birge);
        end

        % ---- MIRROR-SEEDING: use the symmetry to find lines the fit missed ----
        % Every orientation gives a PAIR straddling D, so a line at f REQUIRES a partner
        % at 2D - f. D is known, so this adds no free parameter -- it is physics supplying
        % a start point. It exists because the first run of this script found one line
        % near 2894 where the data plainly showed two, and separately threw away a real
        % line at 2803 whose partner at 2938 sat right where the symmetry says it should.
        % Two distinct failures need two distinct repairs, and they interleave:
        %
        %   ASYMMETRIC OMISSION -- a line found with no partner. Mirror-seeding fixes it,
        %       because D is known and the partner position is therefore known.
        %   SYMMETRIC BLEND -- several lines merged into one, on BOTH sides, so the
        %       spectrum still looks symmetric and mirror-seeding sees nothing wrong.
        %       Only a width test finds it: a blend is much wider than a single line.
        %       Split the widest line in two, refit, and keep it only if chi2 improves
        %       materially. Splitting one side then breaks the symmetry, which the next
        %       mirror pass repairs on the other side -- hence the interleaving.
        for pass = 1:8
            improved = false;
            c = [fit.lines.centre]; w = [fit.lines.fwhm];

            % --- mirror pass. D is ESTIMATED from the outermost pair, not assumed: that
            % pair has the largest splitting so its centres are best determined and least
            % blended. On the run this was built for the true D is 2870.94, and using the
            % nominal 2870 puts every predicted partner 1.9 MHz out -- enough to make a
            % genuinely missing line look like it is already there.
            Dest = (min(c) + max(c))/2;
            if abs(Dest - opts.D0) > opts.dTol; Dest = opts.D0; end
            % SEED THE MIRROR OF EVERY LINE, not only of the unpaired ones.
            %
            % The narrower version only added a mirror when a line had NO partner, which
            % fails on an asymmetric blend -- exactly what 45 G produced. The high side
            % resolved into two lines at 2910.7 and 2924.8 while the low side merged into
            % one 43 MHz blob at 2832.8. Because 2832.8 sat within tolerance of 2910.7's
            % mirror it counted as "paired", so no seed was proposed for 2924.8's mirror at
            % 2817.1 and the blob was never given the components it needed.
            %
            % Seeding every mirror lets the better-resolved side teach the worse one. A
            % seeded line that is not really there fails its depth test and is dropped, so
            % over-seeding costs a refit and nothing else.
            missing = [];
            for k = 1:numel(c)
                mir = 2*Dest - c(k);
                if mir >= min(f) && mir <= max(f) && min(abs(c - mir)) > opts.pairTol
                    missing(end+1) = mir;                                   %#ok<AGROW>
                end
            end
            allMir = 2*Dest - c;
            allMir = allMir(allMir >= min(f) & allMir <= max(f));
            cand = sort([c(:); allMir(:)]);
            keepC = [true; diff(cand) > opts.pairTol];      % dedupe
            cand = cand(keepC);
            if numel(cand) > numel(c)
                f2 = refit(f, v, e, cand, opts);
                [take, why] = acceptRefit(fit, f2, opts);
                if opts.verboseRefine
                    fprintf('  pass %d mirror about D = %.2f, %d seeds (%s): %s\n', ...
                            pass, Dest, numel(cand), sprintf('%.1f ', cand), ...
                            describeRefit(fit, f2, take, why));
                end
                if take
                    fit = f2; improved = true;
                end
            end

            % --- split pass, on the widest line
            c = [fit.lines.centre]; w = [fit.lines.fwhm];
            [wMax, iw] = max(w);
            if wMax > opts.splitRatio * min(w)
                % Try TWO and THREE components. At 45 G the three off-axis orientations
                % span ~22 MHz on each side, so a blend of all three is common and cutting
                % it in two leaves a residual that a two-way split cannot remove -- the
                % chi2 gain test then rejects the split and the blob survives intact.
                best = []; bestChi = inf; bestN = NaN;
                for nSub = [3 2]
                    off = linspace(-wMax/2.5, wMax/2.5, nSub);
                    seeds = c; seeds(iw) = [];
                    seeds = sort([seeds(:); c(iw) + off(:)]);
                    fc = refit(f, v, e, seeds, opts);
                    [tk, wy] = acceptRefit(fit, fc, opts);
                    if opts.verboseRefine
                        fprintf(['  pass %d split %.1f MHz line at %.1f into %d ' ...
                                 '(seeds %s): %s\n'], pass, wMax, c(iw), nSub, ...
                                sprintf('%.1f ', c(iw) + off), ...
                                describeRefit(fit, fc, tk, wy));
                    end
                    if tk && fc.chi2Reduced < bestChi
                        best = fc; bestChi = fc.chi2Reduced; bestN = nSub;
                    end
                end
                f2 = best;
                if ~isempty(f2)
                    fprintf(['  pass %d: accepted the %d-way split of the %.1f MHz line at ' ...
                             '%.1f -> %d lines, chi2 %.1f -> %.1f, widest %.1f -> %.1f MHz\n'], ...
                            pass, bestN, wMax, c(iw), f2.nLines, ...
                            fit.chi2Reduced, f2.chi2Reduced, wMax, max([f2.lines.fwhm]));
                    fit = f2; improved = true;
                end
            end
            if ~improved; break; end
        end

        % ---- PHYSICS-SEEDED PASS: seed from the SOLVED FIELD, not from symmetry -----
        % Mirror seeding uses one fact -- a pair straddles D. The field vector implies far
        % more: three parameters (Bx, By, Bz) plus D fix ALL EIGHT line positions exactly.
        % So fit the field to whatever lines are already resolved, ask it where the rest
        % must be, and hand those positions to the fitter.
        %
        % This is what finally splits a blend. At 45 G the low side of this spectrum is one
        % 44 MHz blob; mirror seeding can only propose the reflection of a resolved line, so
        % it proposes at most as many components as the other side resolved -- and if the
        % other side is also blended it proposes nothing useful. The field solution has no
        % such limit: from the outer pair alone it predicts where all six inner lines sit,
        % including the ones no side has resolved. On the measured 45 G sweep, the five
        % resolved lines give |B| = 44.91 G and predict the missing components at 2830.3 and
        % 2845.3 -- inside the blob, which is precisely the seed the fitter could not
        % generate for itself.
        %
        % It is also self-checking: a seeded line that is not there fails its depth test and
        % is dropped, and if the field solution is wrong its predictions land on noise and
        % the refit is rejected by the same acceptance rule as everything else.
        if opts.physicsSeed
            for pass = 1:3
                c = [fit.lines.centre];
                sol = esrFieldSolve(c, [fit.lines.centreErr], ...
                        struct('D0', opts.D0, 'dTol', opts.dTol, 'verbose', false, ...
                               'sweepRange', [min(f) max(f)]));
                if ~sol.ok
                    if opts.verboseRefine
                        fprintf('  pass %d physics seed: field solve failed (%s)\n', pass, sol.why);
                    end
                    break
                end
                seeds = sol.predicted(sol.predicted >= min(f) & sol.predicted <= max(f));
                seeds = sort([seeds(:); sol.unmatched(:)]);
                keepS = [true; diff(seeds) > opts.pairTol];
                seeds = seeds(keepS);
                f2 = refit(f, v, e, seeds, opts);
                [take, why] = acceptRefit(fit, f2, opts);
                if opts.verboseRefine
                    fprintf(['  pass %d physics seed: |B| = %.2f G, tilt %.2f deg -> %d seeds ' ...
                             '(%s): %s\n'], pass, sol.B, sol.tiltDeg, numel(seeds), ...
                            sprintf('%.1f ', seeds), describeRefit(fit, f2, take, why));
                end
                if ~take; break; end
                fit = f2;
            end
        end
        out.fit = fit;

        % The spectrum is drawn HERE, before the field solve, so it appears even if the
        % solve fails -- but |B|, the tilt and D are not known yet, so the axes handle is
        % kept and annotateSurvey fills them in once they are. That is the only reason this
        % is two steps instead of one.
        axSurvey = [];
        if opts.plotIt
            esrMultiLineFit(f, v, e, struct('baselineOrder', opts.baselineOrder, ...
                'seedCentres', [fit.lines.centre], 'plot', true, 'name', 'field survey'));
            axSurvey = gca;
        end

        c  = [fit.lines.centre];  ce = [fit.lines.centreErr];
        w  = [fit.lines.fwhm];    we = [fit.lines.fwhmErr];
        dp = [fit.lines.depth];   n  = fit.nLines;
        fprintf('\n--- lines ---\n');
        fprintf('%3s %12s %10s %10s %10s\n', '#', 'centre/MHz', '+/-', 'FWHM', 'depth/%');
        for k = 1:n
            fprintf('%3d %12.3f %10.3f %10.2f %10.3f\n', k, c(k), ce(k), w(k), 100*dp(k));
        end
        out.lineCentres = c;

        % ---- THE FIELD, FROM THE EXACT HAMILTONIAN. This is the measurement. ----
        % Four parameters -- Bx, By, Bz, D -- predict all eight line positions exactly, so
        % the field is fitted to the lines rather than reconstructed from pair arithmetic.
        % Everything downstream (tilt, projections, the coil setting that puts a line on
        % target) comes from this, and the pairing below is kept only as a cross-check.
        %
        % IT REPLACES THE MIRROR-PARTNER CLASSIFIER, WHICH HAD A REAL FAILURE MODE. That
        % test dropped any line with no partner at 2D - f. But there is no single D: the
        % second-order transverse shift moves each orientation's pair mean up by an amount
        % that depends on its own B_perp, so at 45 G and 7 degrees the aligned pair straddles
        % 2870.6 while the three off-axis pairs straddle 2877.3 to 2878.4. A D estimated from
        % the data lands on 2878, because three of the four pairs agree on it -- and the
        % ALIGNED pair is then declared partnerless and thrown away. On a synthetic spectrum
        % with a known field that is exactly what happened: 2745.45 and 2995.94, both correct
        % to 0.1 MHz and the most reliable numbers in the sweep, were discarded as "not NV
        % lines" while the shifted pairs were kept and the reported D came out 7 MHz high.
        %
        % The Hamiltonian fit has the shift built in, so each orientation is allowed its own
        % pair mean and the geometry supplies the relation between them. A feature that
        % matches NO predicted transition -- the resonator dip is the candidate -- is
        % reported as unmatched, which is what the mirror test was really trying to do.
        sol = esrFieldSolve(c, ce, struct('D0', opts.D0, 'dTol', opts.dTol, ...
                    'huberMHz', opts.huberMHz, 'sweepRange', [min(f) max(f)], ...
                    'verbose', true));
        out.solve = sol;
        if sol.ok
            out.B = sol.B; out.tiltDeg = sol.tiltDeg; out.D = sol.D;
            out.projections = sol.projections;
            out.solvedSplittings = sol.splittings;
            out.residRms = sol.rms;
            if ~isempty(sol.unmatched)
                keepNV = ~ismember(c, sol.unmatched);
                c = c(keepNV); ce = ce(keepNV); w = w(keepNV);
                we = we(keepNV); dp = dp(keepNV); n = numel(c);
            end
            if sol.rms > 2
                fprintf(['\n  *** residual rms %.2f MHz is too large for a real NV pattern. The\n' ...
                         '      line list does not describe one field. Suspect a spurious line, a\n' ...
                         '      blend being reported as one line, or a fit that has not converged.\n'], ...
                        sol.rms);
            end
        else
            fprintf('\n  field solve failed: %s\n', sol.why);
            fprintf('  falling back to mirror pairing, which cannot handle the second-order shift.\n');
        end
        out.nvLines = c;

        % ---- pairing, as a CROSS-CHECK on the fit above -------------------------
        % Kept because it is independent of the Hamiltonian fit and reads directly off the
        % spectrum, and because the pair widths feed the blend-width alignment metric. It
        % is no longer allowed to discard anything.
        Dcand = opts.D0 + (-opts.dTol:0.05:opts.dTol);
        nMatch = zeros(size(Dcand)); resid2 = zeros(size(Dcand));
        for q = 1:numel(Dcand)
            for k = 1:n
                d = abs(c - (2*Dcand(q) - c(k)));
                d(k) = inf;                       % a line cannot pair with itself
                [dm, ~] = min(d);
                if dm <= opts.pairTol
                    nMatch(q) = nMatch(q) + 1; resid2(q) = resid2(q) + dm^2;
                end
            end
        end
        best = find(nMatch == max(nMatch));
        [~, ib] = min(resid2(best));
        Dfound = Dcand(best(ib));
        out.Dpairing = Dfound;
        % NOTE this Dfound is a PAIRING AID, not the physical D. The scan maximises how
        % many lines pair up, so it lands near the average of all the pair means -- and
        % those means are individually shifted upward by the second-order term. The
        % physical D is the mean of the LEAST-shifted pair, set below once the pairs exist.
        fprintf('\n--- symmetry ---\n');
        fprintf('  D = %.2f MHz pairs up %d of %d lines (a pairing aid; the physical D is\n', ...
                Dfound, max(nMatch), n);
        fprintf('  the least-shifted pair''s mean, reported with the pairs below)\n');

        partnered = false(1, n);
        for k = 1:n
            d = abs(c - (2*Dfound - c(k))); d(k) = inf;
            partnered(k) = min(d) <= opts.pairTol;
        end
        if any(~partnered)
            fprintf(['\n  NOTE %d line(s) have no mirror partner about the single D above:\n'], ...
                    sum(~partnered));
            for k = find(~partnered)
                fprintf('    %.2f MHz (mirror would be %.2f MHz)\n', c(k), 2*Dfound - c(k));
            end
            fprintf(['  That is EXPECTED when the tilt is large, because each orientation''s pair\n' ...
                     '  mean is shifted by its own transverse field. It is not grounds for\n' ...
                     '  dropping a line -- see the Hamiltonian fit above, which handles it.\n']);
        end

        % Pairing for the cross-check table uses only the lines that DO pair up, and by
        % rank: within that subset the kth from the bottom is the partner of the kth from
        % the top. An odd count means one line is pairing with a partner that is itself
        % paired elsewhere, so the last one is left out of the table rather than aborting
        % the whole survey -- the field measurement above does not depend on this.
        cp = c(partnered); cep = ce(partnered); wp = w(partnered); dpp = dp(partnered);
        np = numel(cp);
        if mod(np, 2) == 1
            fprintf(['\n  (%d pairable lines is odd, so the middle one is left out of the pair\n' ...
                     '   table; the field fit above used all of them.)\n'], np);
            np = np - 1;
            mid = ceil((np+2)/2);
            cp(mid) = []; cep(mid) = []; wp(mid) = []; dpp(mid) = [];
        end
        havePairs = np >= 4;
        if ~havePairs
            fprintf(['\n  fewer than two mirror pairs, so the cross-check table and the\n' ...
                     '  blend-width metric are skipped. The Hamiltonian fit above stands\n' ...
                     '  on its own -- it does not need pairs.\n']);
            out.pairs = [];
            out = finishNoPairs(out, opts, sol);
            annotateSurvey(axSurvey, out, sol);
            restore(a, saved); return
        end
        c = cp; ce = cep; w = wp; dp = dpp; n = np;
        nPair = n/2;
        P = struct('lo', {}, 'hi', {}, 'splitting', {}, 'splittingErr', {}, ...
                   'mean', {}, 'meanErr', {}, 'fwhmLo', {}, 'fwhmHi', {}, 'depth', {});
        for k = 1:nPair
            i = k; j = n + 1 - k;
            sp  = c(j) - c(i);
            spe = sqrt(ce(i)^2 + ce(j)^2);        % covariance is lost after the drop, so
                                                   % this is the conservative form
            P(k) = struct('lo', i, 'hi', j, 'splitting', sp, 'splittingErr', spe, ...
                'mean', (c(i)+c(j))/2, 'meanErr', spe/2, ...
                'fwhmLo', w(i), 'fwhmHi', w(j), 'depth', mean([dp(i) dp(j)]));
        end
        out.pairs = P;
        out.Dpairs = min([P.mean]);     % least-shifted pair mean: D only if that pair is
                                        % itself nearly aligned, which is why it does not
                                        % overwrite the fitted D below
        if ~(isstruct(sol) && isfield(sol, 'ok') && sol.ok)
            out.D = out.Dpairs;
        end
        fprintf('\n  least-shifted pair mean = %.3f MHz\n', out.Dpairs);
        if isstruct(sol) && isfield(sol, 'ok') && sol.ok
            fprintf(['  (the FITTED D is %.3f MHz. These differ by %.2f MHz because even the\n' ...
                     '   least-shifted pair in this table carries its own transverse shift --\n' ...
                     '   taking a pair mean as D is only right at perfect alignment.)\n'], ...
                    sol.D, abs(out.Dpairs - sol.D));
        end
        if numel(c) < numel(out.nvLines)
            fprintf(['\n  NOTE this table covers %d of the %d NV lines: the rest do not pair about\n' ...
                     '  a single D, so the cross-check below is PARTIAL and its |B| and tilt are\n' ...
                     '  not comparable to the fit. Read the fit for the field.\n'], ...
                    numel(c), numel(out.nvLines));
        end

        fprintf('\n--- pairs (each is one NV orientation) ---\n');
        fprintf('%4s %14s %12s %9s %12s %10s %12s\n', 'pair', 'lines/MHz', ...
                'splitting', '+/-', 'mean', 'B_par/G', 'FWHM lo/hi');
        for k = 1:nPair
            fprintf('%4d %6.1f/%-7.1f %12.3f %9.3f %12.3f %10.3f  %5.1f/%.1f\n', k, ...
                    c(P(k).lo), c(P(k).hi), P(k).splitting, P(k).splittingErr, ...
                    P(k).mean, P(k).splitting/(2*GAMMA), P(k).fwhmLo, P(k).fwhmHi);
        end
        wRatio = max([P.fwhmLo]./[P.fwhmHi], [P.fwhmHi]./[P.fwhmLo]);
        if max(wRatio) > 1.6
            fprintf(['\n  NOTE the two lines of a pair differ in width by up to %.1fx. They\n' ...
                     '  share an intrinsic linewidth, so that is the resonator envelope\n' ...
                     '  distorting whichever line sits on its steepest part. It biases WIDTHS\n' ...
                     '  much more than centres.\n'], max(wRatio));
        end

        % ---- THE ALIGNMENT KNOB, which works even when the inner lines are blended
        % The three off-axis orientations are degenerate ONLY at perfect alignment, and
        % their splitting spreads FIRST order in the tilt. When the tilt is small they
        % merge into one line at the intrinsic linewidth; when it is not, they broaden and
        % eventually resolve.
        %
        % That makes the INNER blend's width the alignment metric, and crucially it works
        % in the regime where the components CANNOT be separated. Tested on a synthetic
        % replica of this setup at a 6 degree tilt: three lines 5.5 MHz apart with 6.6 MHz
        % linewidths are not resolvable at all, but the blend they form is visibly wider
        % than the outer, single-orientation line. So minimise this number by adjusting
        % Bphi/Btheta -- do NOT try to maximise the outer splitting, which is second order
        % in the tilt and therefore flat.
        [~, iOut] = max([P.splitting]);
        [~, iIn]  = min([P.splitting]);
        wOut = min([P(iOut).fwhmLo P(iOut).fwhmHi]);   % single orientation, so this is
                                                        % the intrinsic linewidth
        wIn  = max([P(iIn).fwhmLo P(iIn).fwhmHi]);
        out.intrinsicFwhm = wOut; out.innerFwhm = wIn;
        out.excessWidth = wIn - wOut;
        fprintf('\n--- alignment knob ---\n');
        fprintf('  outer (single-orientation) FWHM %.2f MHz  <- the intrinsic linewidth\n', wOut);
        fprintf('  inner (blended) FWHM            %.2f MHz\n', wIn);
        fprintf('  excess                          %+.2f MHz\n', out.excessWidth);
        fprintf(['  MINIMISE the excess by adjusting Bphi/Btheta. It is FIRST order in the\n' ...
                 '  tilt, where the outer splitting is second order and therefore useless as\n' ...
                 '  a target. When the excess reaches zero the three off-axis orientations\n' ...
                 '  are degenerate, the spectrum becomes four clean lines, and everything\n' ...
                 '  below becomes exact rather than approximate.\n']);

        % ---- |B| and the tilt, from ALL the projections ----------------------
        % This replaces the aligned-versus-three-fold reading, which assumed the field was
        % already along a <111> axis. It is not: this spectrum shows three or four
        % DISTINCT splittings, which only happens when the three off-axis orientations
        % are non-degenerate, i.e. the field is tilted.
        %
        % The general result needs no alignment assumption. For a field B along n, the
        % projections onto the four <111> axes are p_i = B * (n . a_i), and those four
        % numbers obey two exact identities:
        %
        %       SUM p_i    = 0            (the four <111> vectors sum to zero)
        %       SUM p_i^2  = (4/3) B^2
        %
        % Each pair measures |p_i| = splitting / (2*gamma). So the second identity gives
        % |B| with no geometry fitted at all, and the first is a free consistency check
        % once the signs are chosen -- and there are only a handful of sign patterns.
        pAbs = [P.splitting] / (2*GAMMA);
        pErr = [P.splittingErr] / (2*GAMMA);
        out.pAbs = pAbs; out.pErr = pErr;

        if nPair >= 4
            pUse = pAbs(1:4); pUseE = pErr(1:4); inferred = false;
        elseif nPair == 3
            % The fourth projection follows from SUM p_i = 0 once the signs are set, so a
            % missing pair is recoverable rather than fatal.
            inferred = true;
        else
            inferred = true;
        end

        % This block is now a CROSS-CHECK on the Hamiltonian fit, not the measurement. It
        % uses only the two exact sum rules and the measured splittings, so it shares no
        % machinery with the fit -- which is what makes agreement between them meaningful.
        % Where they disagree, the fit wins: this route ignores the second-order shift, and
        % it needs three resolved pairs where the fit needs none.
        fprintf('\n--- cross-check: |B| from the sum rules, no alignment assumed ---\n');
        if nPair >= 3
            [Bfit, tiltDeg, signPattern, closure, pFull] = solveProjections(pAbs, GAMMA);
            out.Bsumrule = Bfit; out.tiltSumrule = tiltDeg; out.closure = closure;
            out.signPattern = signPattern;
            if ~(isstruct(sol) && isfield(sol, 'ok') && sol.ok)
                out.B = Bfit; out.tiltDeg = tiltDeg; out.projections = pFull;
            end
            fprintf('  |B|            = %.3f G      (from SUM p^2 = 4/3 B^2)\n', Bfit);
            fprintf('  projections    = %s G\n', mat2str(round(pFull*1000)/1000));
            fprintf('  closure SUM p  = %+.4f G   (must be 0; %.2f%% of |B|)\n', ...
                    closure, 100*abs(closure)/max(Bfit,eps));
            fprintf('  tilt from the nearest <111> axis = %.2f deg\n', tiltDeg);
            if out.excessWidth > 0.3*wOut
                fprintf(['  CAVEAT: the inner lines are blended (excess %.1f MHz), so the inner\n' ...
                         '  splittings are averages over unresolved components and |B| and the\n' ...
                         '  tilt are APPROXIMATE. The outer pair is clean, so the largest\n' ...
                         '  projection is reliable; the rest is not, until the alignment is\n' ...
                         '  improved enough to resolve or merge them.\n'], out.excessWidth);
            end
            if nPair == 3
                fprintf('  (only 3 pairs seen; the 4th projection was INFERRED from the closure)\n');
            end
            if abs(closure) > 0.05*Bfit
                fprintf(['  *** closure is %.1f%% of |B|, which is too large. Either a pair is\n' ...
                         '      misassigned, or a splitting is wrong, or there is a fifth\n' ...
                         '      feature in the sweep that is not an NV line.\n'], ...
                        100*abs(closure)/Bfit);
            end
            if isstruct(sol) && isfield(sol, 'ok') && sol.ok && numel(c) < numel(out.nvLines)
                % A partial pair table cannot give a meaningful |B|: the sum rules need all
                % four projections, and the missing pair is usually the ALIGNED one, whose
                % projection is the largest of the four. Saying nothing further is better
                % than editorialising on a number known to be incomplete.
                fprintf(['  ^ PARTIAL: %d of %d lines paired, so these two numbers are not a\n' ...
                         '    measurement of anything. They are shown only to make the omission\n' ...
                         '    visible.\n'], numel(c), numel(out.nvLines));
            elseif tiltDeg > 3
                fprintf(['\n  MISALIGNED by ~%.1f deg. That is why the off-axis orientations are\n' ...
                         '  resolved as separate lines instead of one degenerate line. For Rabi\n' ...
                         '  and Ramsey you want a single clean orientation, so bring the tilt\n' ...
                         '  down until the three off-axis lines merge -- that merging is a\n' ...
                         '  first-order indicator and much sharper than maximising the largest\n' ...
                         '  splitting, which is second order in the tilt.\n'], tiltDeg);
            else
                fprintf('\n  Within ~%.1f deg of a <111> axis.\n', tiltDeg);
            end
            if isstruct(sol) && isfield(sol, 'ok') && sol.ok
                fprintf(['\n  the two routes: |B| %.3f G (fit) vs %.3f G (sum rules), tilt %.2f\n' ...
                         '  vs %.2f deg. The sum-rule route neglects the second-order shift and\n' ...
                         '  reads blended splittings as single ones, so a gap of a few tenths of\n' ...
                         '  a degree is expected; a gap of degrees means a misassigned pair.\n'], ...
                        sol.B, Bfit, sol.tiltDeg, tiltDeg);
            end
        else
            fprintf(['  Only %d pair(s) resolved, so the sum-rule cross-check cannot run: a\n' ...
                     '  single splitting could belong to any of the four orientations, whose\n' ...
                     '  projections differ by up to 3x. The Hamiltonian fit above does not have\n' ...
                     '  this limitation and its answer stands.\n'], nPair);
        end

        % ---- second-order check, in its general form -------------------------
        % A transverse field pushes BOTH transitions of an orientation up by
        % 1.5*(gamma*B_perp)^2/D, and B_perp follows from |B| and that orientation's own
        % B_parallel: B_perp^2 = B^2 - B_par^2. So every pair's MEAN is predicted, with no
        % free parameter, once |B| is known. The most-aligned orientation has the smallest
        % B_perp and therefore the smallest shift, which is why its mean is the best
        % estimate of D.
        %
        % This is a genuine test and not decoration: the shifts here are ~2 MHz, ten times
        % the centre errors, so a wrong pairing shows up immediately.
        fprintf('\n--- second-order shift check (each pair''s mean is predicted) ---\n');
        fprintf('%4s %12s %12s %12s\n', 'pair', 'mean obs', 'mean pred', 'diff');
        Dbase = min([P.mean]);          % least-shifted pair
        worst = 0;
        for k = 1:nPair
            bPar = P(k).splitting/(2*GAMMA);
            bPerp2 = max(out.B^2 - bPar^2, 0);
            % Referenced to the LEAST-shifted pair, whose own shift is unknown but
            % common to all of them -- so only the DIFFERENCE of shifts is predicted,
            % which needs no absolute D.
            bParRef = max([P.splitting])/(2*GAMMA);
            bPerpRef2 = max(out.B^2 - bParRef^2, 0);
            pred = Dbase + 1.5*GAMMA^2*(bPerp2 - bPerpRef2)/Dbase;
            d = P(k).mean - pred;
            fprintf('%4d %12.3f %12.3f %+12.3f\n', k, P(k).mean, pred, d);
            worst = max(worst, abs(d));
        end
        out.secondOrderWorst = worst;
        if worst < 1.0
            fprintf('  => consistent to %.2f MHz. The pairing and |B| hold together.\n', worst);
        else
            fprintf(['  => off by up to %.2f MHz. Suspect the pairing, or blended inner lines\n' ...
                     '     making a splitting an average over unresolved components.\n'], worst);
        end

        % ---- what to change to hit the target frequency ---------------------
        out = reportTarget(out, opts, sol, P, GAMMA);

    catch err
        restore(a, saved);
        rethrow(err);
    end
    annotateSurvey(axSurvey, out, sol);
    restore(a, saved);
    fprintf('\n(no field setting was changed)\n');
end

% ======================================================================

function annotateSurvey(ax, out, sol)
%ANNOTATESURVEY Put the measurement on the figure: |B|, the tilt, D, and where each
% orientation's pair is predicted to sit.
%
% The spectrum plot on its own showed the fit but none of the result, so the numbers had to
% be read out of the printout and matched to the picture by eye. The whole point of the
% survey is |B| and the tilt, so those belong in the title.
%
% The four orientations are NOT drawn as four equal categories. One of them is the aligned
% axis -- the one every later calibration uses -- and the other three are the contamination
% you are trying to get away from. Drawing that hierarchy (aligned emphasised, the rest
% recessive grey) says more than four colours would, and it keeps identity off colour alone
% since each pair is labelled with its own projection.
%
% Silent no-op when there is no figure, and wrapped so a plotting problem can never lose a
% good measurement.
    if isempty(ax) || ~all(ishandle(ax)); return; end
    try
        COL.aligned = [0.00 0.45 0.35];      % green: the axis that matters
        COL.other   = [0.55 0.55 0.55];      % grey: the other three
        hold(ax, 'on');
        yl = ylim(ax);

        haveSol = ~isempty(sol) && isstruct(sol) && isfield(sol, 'ok') && sol.ok;
        alignedIdx = 1;
        if haveSol && isfield(sol, 'alignedAxisIdx') && ~isempty(sol.alignedAxisIdx)
            alignedIdx = sol.alignedAxisIdx;
        end

        hA = []; hO = [];
        if haveSol && isfield(sol, 'means') && isfield(sol, 'splittings')
            m = sol.means(:)'; sp = sol.splittings(:)';
            for i = 1:min(numel(m), numel(sp))
                isAligned = (i == alignedIdx);
                if isAligned; c = COL.aligned; lw = 1.6; ls = '-';
                else;         c = COL.other;   lw = 1.0; ls = ':';
                end
                for pos = m(i) + [-1 1]*sp(i)/2
                    hh = plot(ax, pos*[1 1], yl, ls, 'Color', c, 'LineWidth', lw);
                    if isAligned && isempty(hA); hA = hh;
                    elseif ~isAligned && isempty(hO); hO = hh;
                    else; set(hh, 'HandleVisibility', 'off');
                    end
                end
                % Label each pair with its own projection, so identity is never colour alone.
                if isfield(sol, 'projections') && numel(sol.projections) >= i
                    txt = sprintf(' %+.1f G', sol.projections(i));
                    if isAligned; txt = sprintf(' ALIGNED %+.1f G', sol.projections(i)); end
                    text(ax, m(i) + sp(i)/2, yl(2), txt, 'Color', c, 'FontSize', 8, ...
                         'VerticalAlignment', 'top', 'Rotation', 90);
                end
            end
            ylim(ax, yl);
        end

        % ---- the title carries the measurement ----
        if isfinite(out.B) && isfinite(out.tiltDeg)
            t = sprintf('field survey:   |B| = %.2f G     tilt = %.2f deg', out.B, out.tiltDeg);
            if isfield(out, 'D') && isfinite(out.D)
                t = sprintf('%s     D = %.2f MHz', t, out.D);
            end
        else
            t = 'field survey:   field NOT solved -- see the printout';
        end
        title(ax, t, 'FontWeight', 'bold');

        % Everything else that qualifies the number goes in the subtitle, which is where a
        % reader looks after the headline: how well it fits, and how far off aligned it is.
        bits = {};
        if haveSol && isfield(sol, 'rms') && isfinite(sol.rms)
            bits{end+1} = sprintf('line residual %.2f MHz', sol.rms);
        end
        if isfield(out, 'projections') && ~isempty(out.projections)
            bits{end+1} = sprintf('projections %s G', ...
                                  strtrim(sprintf('%+.1f ', out.projections)));
        end
        if isfinite(out.tiltDeg)
            if out.tiltDeg > 2
                bits{end+1} = 'ALIGN FIRST (tilt > 2 deg)';
            else
                bits{end+1} = 'aligned well enough to proceed';
            end
        end
        if ~isempty(bits)
            try
                subtitle(ax, strjoin(bits, '   |   '));   % R2020b+
            catch
                % No subtitle on this release: fold it into the xlabel rather than lose it.
                xlabel(ax, sprintf('%s        [%s]', get(get(ax,'XLabel'),'String'), ...
                                   strjoin(bits, ' | ')));
            end
        end

        if ~isempty(hA) || ~isempty(hO)
            hs = {hA, hO};
            ls = {'aligned pair (predicted)', 'other orientations (predicted)'};
            keep = cellfun(@(h) ~isempty(h) && all(ishandle(h)), hs);
            legend(ax, [hs{keep}], ls(keep), 'Location', 'southeast', 'Box', 'off');
        end
    catch err
        fprintf('  (could not annotate the survey plot: %s)\n', err.message);
    end
end

% ======================================================================

function out = finishNoPairs(out, opts, sol)
% The useful half of the report when fewer than two mirror pairs resolved. This used to
% be a dead end -- "at least two pairs are needed" and return -- which is why a sweep with
% one clean pair and one 44 MHz blend produced no field at all. The Hamiltonian fit needs
% three lines, not two pairs, so there is a measurement to report here.
    if ~(isstruct(sol) && isfield(sol, 'ok') && sol.ok)
        fprintf(['\n  and the field fit did not converge either, so this sweep yields nothing.\n' ...
                 '  Widen the sweep, or raise |B| until more lines separate.\n']);
        return
    end
    reportTarget(out, opts, sol, [], 2.8025);
end

function out = reportTarget(out, opts, sol, P, GAMMA)
% The coil change that puts a line on the target frequency, SOLVED rather than scaled.
%
% The old version scaled the measured splitting linearly: need/have. That is wrong at the
% fields in use, and wrong in the direction that matters. The splitting is linear in |B|
% only to first order; the transverse shift is quadratic, and at 45 G with a 7 degree tilt
% it moves the off-axis pair means by 7 MHz. Since this number is what the coil actually
% gets set from, it is solved by bisection on the exact Hamiltonian instead -- ask for the
% |B| scaling that puts a line ON the target, and let the model include every term.
    if ~(isstruct(sol) && isfield(sol, 'ok') && sol.ok)
        if isempty(P); return; end
        [~, iBig] = max([P.splitting]);
        need = 2*abs(out.D - opts.targetMHz);
        fprintf('\n--- to move the LARGEST-splitting line onto %.1f MHz (first order only) ---\n', ...
                opts.targetMHz);
        fprintf('  need splitting %.3f MHz, have %.3f -> scale the field by %.4f\n', ...
                need, P(iBig).splitting, need/max(P(iBig).splitting, eps));
        return
    end

    fprintf('\n--- to move a line onto %.1f MHz ---\n', opts.targetMHz);
    [scale, which, fAll, sols] = sol.fieldForLineAt(opts.targetMHz);
    % WHICH line goes on target is a physics choice. For Rabi and Ramsey you want the
    % ALIGNED orientation -- the one whose projection is largest -- because it is the only
    % one that is a single orientation rather than one of a near-degenerate trio. Taking
    % whichever line reaches the target at the smallest field change would often hand you an
    % off-axis line, which is the opposite of what the alignment work is for.
    if ~isempty(sols) && strcmpi(opts.whichLine, 'aligned') && ...
            isfield(sol, 'alignedAxisIdx')
        keep = find([sols.axis] == sol.alignedAxisIdx);
        if ~isempty(keep)
            [~, kk] = min(abs([sols(keep).scale] - 1));
            pick = sols(keep(kk));
            scale = pick.scale; which = pick.line; fAll = pick.lines;
            fprintf('  (aligned orientation requested, so the off-axis solutions are ignored)\n');
        else
            fprintf(['  (no |B| in range puts an ALIGNED line on target, so this is an\n' ...
                     '   off-axis line -- not what you want for Rabi)\n']);
        end
    end
    if ~isfinite(scale)
        fprintf(['  no |B| between 0.05x and 4x the present value puts any line on %.1f MHz\n' ...
                 '  with the field in its present DIRECTION. Change the target, or align\n' ...
                 '  first -- the tilt sets which lines exist at all.\n'], opts.targetMHz);
        return
    end
    out.targetScale = scale;
    out.targetLines = fAll;
    fprintf('  scale |B| by %.4f: %.3f G -> %.3f G  (direction held fixed)\n', ...
            scale, out.B, out.B*scale);
    fprintf('  that puts transition %d of orientation %d on target. Then the lines sit at\n', ...
            2 - mod(which, 2), ceil(which/2));
    fprintf('    %s MHz\n', sprintf('%.1f  ', sort(fAll)));
    fprintf(['  Solved on the exact Hamiltonian, not by scaling the splitting: the splitting\n' ...
             '  is linear in |B| but the transverse shift is quadratic, and at these fields\n' ...
             '  the difference is several MHz -- more than a linewidth.\n']);
    if isfinite(out.tiltDeg) && out.tiltDeg > 2
        fprintf(['\n  BUT ALIGN FIRST: the tilt is %.1f deg. Changing |B| moves every line, and\n' ...
                 '  there is no point putting a line on target while the orientation you want is\n' ...
                 '  still contaminated by its neighbours -- and the alignment itself changes\n' ...
                 '  which scaling is right.\n'], out.tiltDeg);
    end
    fprintf(['\n  And treat %.0f MHz as a HYPOTHESIS: it came from the circulators'' reverse\n' ...
             '  ports, which measure power failing to get INTO the antenna, not B1 at the NVs.\n' ...
             '  A Rabi at two or three line positions settles it.\n'], opts.targetMHz);
    if ~isempty(P) && numel(P) >= 1 && GAMMA > 0; end
end

function [take, why] = acceptRefit(fit, f2, opts)
% Decide whether a seeded refit replaces the current fit, and SAY WHY EITHER WAY.
%
% Three ways in, because "did it get better" means three different things here:
%
%   MORE LINES, not worse. A recovered component is the point of the exercise. But
%       "more lines" alone was the old mirror-pass rule, and that would accept an
%       invented line that made the fit worse, so chi2 must not degrade materially.
%
%   FEWER PARAMETERS FIT NO WORSE -- never accepted. Not our business here.
%
%   SAME LINES, materially better chi2. splitGain, as before.
%
%   THE WIDEST LINE GOT NARROWER. This is the criterion that was missing, and it is
%       the one that matters for a blend. Splitting a 44 MHz blob into two 12 MHz
%       lines barely moves reduced chi2 -- the blob already absorbs most of the
%       deficit, and the residual it leaves is a few tenths of a percent spread over
%       30 points against a chi2 already inflated to 2.1 by the envelope. So a
%       chi2-only test rejects a split that is obviously right when you look at the
%       widths: an NV line here is 7.5 to 15 MHz FWHM, and nothing physical is 44.
%       Requiring the widest FWHM to fall to widthGain of its old value, with chi2 not
%       degrading, accepts exactly that case and nothing else.
    take = false;
    if isempty(f2) || ~isstruct(f2) || ~isfield(f2, 'ok') || ~f2.ok
        if isstruct(f2) && isfield(f2, 'why') && ~isempty(f2.why)
            why = sprintf('fit failed (%s)', f2.why);
        else
            why = 'fit failed';
        end
        return
    end
    chi0 = fit.chi2Reduced;  chi1 = f2.chi2Reduced;
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
        why = sprintf('LOST %d line(s) -- the seeded component did not survive its %.1f-sigma depth test', ...
                      fit.nLines - f2.nLines, opts.minDepthSigma);
    else
        why = sprintf('no gain: %d lines, chi2 %.2f -> %.2f, widest %.1f -> %.1f MHz', ...
                      f2.nLines, chi0, chi1, w0, w1);
    end
end

function s = describeRefit(~, f2, take, why)
    if isstruct(f2) && isfield(f2, 'ok') && f2.ok
        if take; verdict = 'ACCEPTED'; else; verdict = 'rejected'; end
        s = sprintf('%d lines, chi2 %.2f, widest %.1f MHz -- %s (%s)', ...
                    f2.nLines, f2.chi2Reduced, max([f2.lines.fwhm]), verdict, why);
    else
        s = sprintf('REJECTED (%s)', why);
    end
end

function [f, v, e] = unpackReplay(r)
% Pull a swept spectrum out of whatever the caller handed over: an ESR EXPERIMENT, the
% output struct of a previous calibFieldSurvey, a bare struct with the three vectors, or a
% matrix.
%
% The experiment case used to be missing, so calibFieldSurvey(struct('replay', a)) failed
% with "replay must be a struct or a matrix, got a AWGExpESR" -- while the run script and
% the level 3 notes both said an experiment was accepted. Rather than copy the extraction
% in here, it lives in esrPullData and esrFitRobust calls the same function, so the two
% cannot drift apart again.
%
% verbose is off: a survey prints its own report and does not need a line about error bars
% in the middle of it.
    [f, v, e] = esrPullData(r, struct('verbose', false, 'name', 'calibFieldSurvey replay'));
end

function val = pickField(s, names)
    val = [];
    for k = 1:numel(names)
        if isfield(s, names{k}) && ~isempty(s.(names{k})) && isnumeric(s.(names{k}))
            val = s.(names{k}); return
        end
    end
end

function o = refit(f, v, e, seeds, opts)
    o = esrMultiLineFit(f, v, e, struct('baselineOrder', opts.baselineOrder, ...
            'seedCentres', seeds, 'minDepthSigma', opts.minDepthSigma, ...
            'maxLines', 12, 'minSepMHz', 3));
end

function [B, tiltDeg, bestSign, closure, pFull] = solveProjections(pAbs, GAMMA) %#ok<INUSD>
% |B| and the tilt from the measured projection magnitudes, assuming nothing about
% alignment.
%
% The four <111> unit vectors sum to zero and each has unit length, so for ANY field
% direction the four projections p_i = B*(n.a_i) satisfy exactly
%
%       SUM p_i = 0        and        SUM p_i^2 = (4/3) B^2
%
% Each measured pair gives |p_i| = splitting/(2*gamma) but not its sign. There are only
% 2^4 sign patterns and a global flip is irrelevant, so eight to try -- pick the one
% that makes the closure SUM p_i smallest. With three pairs the fourth projection is
% not measured at all, and the closure DEFINES it: p_4 = -(p_1+p_2+p_3).
%
% |B| then comes from the second identity with no geometry fitted, and the tilt from the
% largest projection: cos(tilt) = max|p_i| / B, because the axis with the largest
% projection is the one the field is closest to.
    pAbs = pAbs(:)';
    if numel(pAbs) >= 4
        p4 = pAbs(1:4);
        best = inf; bestSign = ones(1,4);
        for m = 0:15
            sg = 1 - 2*bitget(m, 1:4);
            if sg(1) < 0; continue; end          % global flip is not a distinct case
            cl = sum(sg .* p4);
            if abs(cl) < best; best = abs(cl); bestSign = sg; end
        end
        pFull   = bestSign .* p4;
        closure = sum(pFull);
    else
        % Three measured, fourth inferred from the closure. Try the sign patterns on the
        % three and keep the one whose inferred fourth is closest to a plausible
        % magnitude, i.e. the smallest inferred |p_4| consistent with the others.
        p3 = pAbs(1:3);
        best = inf; bestSign = ones(1,3); p4v = NaN;
        for m = 0:7
            sg = 1 - 2*bitget(m, 1:3);
            if sg(1) < 0; continue; end
            cand = -sum(sg .* p3);
            % A physical fourth projection cannot exceed the largest measured one by
            % much: all four are projections of the same vector onto unit axes.
            pen = abs(cand) - max(p3);
            if pen < best; best = pen; bestSign = sg; p4v = cand; end
        end
        pFull   = [bestSign .* p3, p4v];
        closure = sum(pFull);
        bestSign = [bestSign, sign(p4v)];
    end
    B = sqrt(3/4 * sum(pFull.^2));
    tiltDeg = acosd(min(max(max(abs(pFull))/max(B, eps), 0), 1));
end

function zeroFieldHint(c, dp, w, opts)
% When the pairing is refused, say whether the spectrum looks like ZERO FIELD rather
% than leaving the user with only a complaint. The signature is one line that both sits
% at D and dominates everything else in depth: at zero field all four orientations are
% degenerate, so the real line is several times deeper than any artefact.
    [~, iDeep] = max(dp);
    others = dp; others(iDeep) = [];
    dominates = isempty(others) || dp(iDeep) > 2*max(others);
    if abs(c(iDeep) - opts.D0) <= opts.dTol && dominates
        fprintf(['\n    This looks like ZERO FIELD: the dominant line is at %.2f MHz, which is\n' ...
                 '    D, and it is %.1fx deeper than anything else in the sweep. At zero field\n' ...
                 '    all four orientations are degenerate there, and the shallower broad\n' ...
                 '    features are the resonator envelope rather than split lines.\n' ...
                 '    Useful anyway: D = %.2f MHz, FWHM = %.2f MHz, contrast = %.2f%%.\n' ...
                 '    Turn the field ON and re-run.\n'], ...
                c(iDeep), dp(iDeep)/max([others eps]), c(iDeep), w(iDeep), 100*dp(iDeep));
    end
end

function a = makeExp(opts)
% A FRESH AWGExpESR for every run, and every mode flag set EXPLICITLY.
%
% The "parameters have changed" dialog is raised when an EXISTING experiment object has
% its parameters altered. A newly constructed object has no previous state to compare
% against, so the dialog never appears -- which is why this builds one rather than
% fetching the singleton with getExpByName.
%
% It also removes a whole class of bug we have already been bitten by: a fresh object
% cannot inherit useSequenceJump, useFgen or driveMWSwitches from whatever ran last.
% esrJumpValidate aborted on exactly that, and the FGEN work spent a round on a flag
% left behind by a previous cell.
%
% The AWG waveform cache is unaffected: loadedFingerprint lives on the AWG object, not
% on the experiment, so an unchanged frequency list still skips the upload.
    if opts.freshExperiment
        a = AWGExpESR;
    else
        a = getExpByName('AWG ESR');
    end
    if ~isempty(opts.jsonLocation)
        a.loadJsonParams('all', opts.jsonLocation);
    end
    a.frequency  = opts.freqs(:)';
    a.averages   = opts.averages;
    a.repeats    = opts.repeats;
    if ~isempty(opts.amplitude);  a.amplitude  = opts.amplitude;  end
    a.isTracking = opts.isTracking;
    % ---- the CW-ESR window: ALWAYS set, ALWAYS clamped ----
    % This used to be "if ~isempty(opts.detectionMHz)", which had two failure modes and both
    % were silent. With the option EMPTY nothing was set, so the window stayed at whatever
    % loadJsonParams had just copied in -- the few-microsecond PULSED value, which is far
    % below the transient floor and throws away most of the dip. And the clamp below sat
    % inside the same branch, so the empty case skipped it too and the json's pulsed
    % laserInitializationDuration went straight through.
    %
    % Empty now means "resolve it", not "leave it alone".
    det = opts.detectionMHz;
    if isempty(det)
        try
            det = esrWindowPick(struct('quiet', true, 'jsonLocation', opts.jsonLocation, ...
                                       'frequency', mean(opts.freqs))).detectionDuration;
        catch
            det = 150;      % a working CW window on this setup if the arithmetic is unavailable
        end
    end
    a.detectionDuration          = det;
    a.referenceDetectionDuration = det;

    % ---- the CW-ESR trigger constraint ----
    % loadJsonParams has just copied laserInitializationDuration in from
    % expProperties.json, where it is a real initialisation for the PULSED experiments and
    % runs to tens of microseconds. In CW ESR it is not an initialisation at all -- there is
    % no init event and the green laser is on throughout -- it only sets how much extra
    % microwave the AWG plays past the window. Leave the json value in place with a short
    % window and the burst becomes longer than one repeat, so triggers are dropped, those
    % repeats get no microwave, and the spectrum comes back as a COMB: alternate points at
    % exactly ratio 1 with the dip only on the others. It does not error, it does not warn,
    % and it will happily be fitted.
    %
    % This is the one place level 3 needs it: calibFieldAlign and calibB1Rabi both run their
    % ESR by calling this function, so neither needs its own copy.
    a.laserInitializationDuration = esrSafeInit(a.detectionDuration, ...
                                                a.delayBetweenDetections, ...
                                                a.laserInitializationDuration);
    a.useSequenceWalk = true;     % the fast, outlier-free path
    a.useSequenceJump = false;    % never: produces scrambled frequencies
    a.useFgen         = false;    % never yet: ~4 s of settling per point makes it slower
    a.driveMWSwitches = true;     % normal operation; false is a switch diagnostic
    % Optional, passed down from the caller: a new random order every average (walk mode).
    % Empty leaves the object's own value.
    if ~isempty(opts.newOrderEveryAverage)
        a.newOrderEveryAverage = logical(opts.newOrderEveryAverage);
    end
end

function restore(a, saved)
% Nothing to restore when the object was built fresh -- it is discarded either way.
    if isempty(a) || isempty(saved.frequency); return; end
    try
        a.frequency = saved.frequency; a.repeats = saved.repeats;
        a.averages = saved.averages;   a.amplitude = saved.amplitude;
    catch
    end
end

function o = withDefaults(o, d)
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(o, fn{i}) || isempty(o.(fn{i})); o.(fn{i}) = d.(fn{i}); end
    end
end
