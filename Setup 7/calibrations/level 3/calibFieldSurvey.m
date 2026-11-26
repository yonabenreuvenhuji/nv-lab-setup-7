function out = calibFieldSurvey(opts)
%CALIBFIELDSURVEY One wide ESR at the current field, and everything it can tell you.
%
%   out = calibFieldSurvey                     % target: the resonance this sweep finds
%   out = calibFieldSurvey(struct('freqs', 2740:2:2960, 'targetMHz', 2790))   % your own
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
% What replaces it: esrFitRobust with baselineOrder 2, so the envelope's curvature
% across a 200 MHz sweep is absorbed by the baseline rather than by the line centres,
% plus the pair-mean check above as an independent readout of whatever bias is left.
%
% ===================================================================
% WHAT esrFitRobust DOES FOR THIS, AND WHAT IS LEFT HERE
% ===================================================================
%
% The line fit, the resonator's own dip (fitted as its own component and kept out of the
% lines -- it moves, so it is detected, not looked up), the mirror / split / physics-seeded
% refits, the field solve and the figure all live in esrFitRobust; one ESR fitter for the
% whole setup. This file keeps what only a survey needs: running the sweep, the detailed
% field-solve report, the mirror-pair cross-check and alignment knob, the sum-rule and
% second-order checks, and the coil scaling that puts a line on targetMHz.
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
%   detectionMHz CW detection and reference window, us. Default [] = AWGExpESR's own
%                (esrWindowPick, set by its loadJsonParams), which is what you want.
%   jsonLocation expProperties.json, loaded with loadJsonParams('all', ...) before
%                anything else is set
%   isTracking   (default 0)
%   freshExperiment  construct a NEW AWGExpESR for every run (default true). See makeExp
%                -- this is what stops the "parameters have changed" dialog appearing.
%   newOrderEveryAverage  walk mode: draw a new random frequency order before every
%                average (true/false). Default [] = leave AWGExpESR's own setting (true
%                on a fresh object unless the json sets it). Each average then pays one
%                sequence upload; see AWGExpESR.newOrderEveryAverage.
%   targetMHz    the frequency you want a line moved to. Default [] = the RESONANCE this
%                sweep measured (out.resonance.centre, found by esrFitRobust), because that
%                is where the resonator puts its power and it moves: it was 2820 when this
%                default was a typed-in number, and has since been seen near 2763. With no
%                resonance in the sweep and no targetMHz given, the coil-target section says
%                so and is skipped. A number here always wins. out.target records which.
%   whichLine    'aligned' or 'triple': which pair to move onto targetMHz
%                (default 'aligned')
%   baselineOrder passed to esrFitRobust (default 2)
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
%   minDepthSigma passed to esrFitRobust as refineDepthSigma, the depth test for the
%                seeded refits (default 2.5)
%   physicsSeed, huberMHz   passed to esrFitRobust (defaults true, 6)
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
        'repeats', 50, 'averages', 3, 'detectionMHz', [], 'targetMHz', [], ...
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
        % uncertainty below inherits it. esrFitRobust's Birge scaling corrects the
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

        % ---- the lines, the resonance and the field: esrFitRobust ---------------
        % One call does what this section used to do by hand: the joint line fit, the
        % resonator's own dip fitted as its own component and kept OUT of the lines (it
        % used to be fitted as one, and once paired with the real line at D to report
        % 8.26 G with the field off), the mirror / split / physics-seeded refits that
        % recover missed and blended lines, and the Hamiltonian field solve. The refits are
        % forced on (refine = true) so a narrow window is treated as it always was here.
        %
        % seedCentres, when given, starts the fit from known line positions instead of
        % searching for them. During an angle scan the previous point's lines are only a
        % few MHz away, so every fit becomes a small refinement rather than a fresh search
        % from a peak-finder that can and does fail on a marginal line.
        fit = esrFitRobust(struct('f', f, 'y', v, 'e', e, ...
                    'baselineOrder', opts.baselineOrder, 'seedCentres', opts.seedCentres, ...
                    'maxLines', 10, 'D0', opts.D0, 'dTol', opts.dTol, ...
                    'huberMHz', opts.huberMHz, 'refine', true, 'pairTol', opts.pairTol, ...
                    'refineDepthSigma', opts.minDepthSigma, 'splitRatio', opts.splitRatio, ...
                    'splitGain', opts.splitGain, 'widthGain', opts.widthGain, ...
                    'physicsSeed', opts.physicsSeed, 'verboseRefine', opts.verboseRefine, ...
                    'plotIt', opts.plotIt, 'verbose', true, 'name', 'field survey'));
        out.fit = fit;
        out.resonance = fit.resonance;
        if ~fit.ok
            fprintf('\nFIT FAILED: %s\n', fit.why);
            restore(a, saved); return
        end

        % ---- the target: the resonance this sweep just measured, unless one was given --
        [opts, out] = resolveTarget(opts, out, fit);

        % ---- is the MODEL right, before trusting any number from it ----------
        % A reduced chi2 of 20 does not mean noisy data, it means the model is the wrong
        % shape -- and Birge scaling then inflates the error bars to cover the model
        % error, which makes them look honest while the CENTRES stay biased.
        if isfinite(fit.chi2Reduced) && fit.chi2Reduced > opts.maxChi2
            fprintf(['\n  reduced chi2 = %.1f, so the model is imperfect (errors scaled %.2fx\n' ...
                     '  to cover it). Usually the resonator envelope, or a blend the refits\n' ...
                     '  could not split. A narrow-window refit around the line in question is\n' ...
                     '  the fix for both.\n'], fit.chi2Reduced, fit.birge);
        end

        c  = [fit.lines.centre];  ce = [fit.lines.centreErr];
        w  = [fit.lines.fwhm];    we = [fit.lines.fwhmErr];
        dp = [fit.lines.depth];   n  = fit.nLines;
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
        % matches NO predicted transition is reported as unmatched, which is what the mirror
        % test was really trying to do. (The resonator's own dip no longer reaches this
        % list: esrFitRobust has already fitted it as the resonance, out.resonance.)
        %
        % esrFitRobust ran this same solve, quietly; it is repeated here with verbose on
        % because the per-line assignment table is the survey's main report.
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
    restore(a, saved);
    fprintf('\n(no field setting was changed)\n');
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
    if isempty(opts.targetMHz) || ~isfinite(opts.targetMHz)
        fprintf(['\n--- coil target: none ---\n' ...
                 '  No resonance was found in this sweep and no targetMHz was given, so there is\n' ...
                 '  no frequency to put a line on. Sweep across the resonator, or pass targetMHz,\n' ...
                 '  for the coil scaling.\n']);
        return
    end
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
    if isfield(out, 'target') && strcmp(out.target.source, 'resonance')
        fprintf(['\n  The target is the resonance THIS sweep measured (%.1f +- %.1f MHz): where the\n' ...
                 '  resonator takes its power, so the best guess for the B1 maximum -- but a dip in\n' ...
                 '  the ratio is not B1 itself. A Rabi at two or three line positions (calibB1Rabi)\n' ...
                 '  settles it.\n'], opts.targetMHz, out.target.err);
    else
        fprintf(['\n  And treat %.0f MHz as a HYPOTHESIS unless it came from a Rabi: a typed-in\n' ...
                 '  target goes stale (the resonator has moved from ~2820 to ~2763). Leave\n' ...
                 '  targetMHz empty to use the resonance the sweep measures.\n'], opts.targetMHz);
    end
    if ~isempty(P) && numel(P) >= 1 && GAMMA > 0; end
end

function [opts, out] = resolveTarget(opts, out, fit)
% targetMHz empty -> the resonance esrFitRobust found in this sweep. The number used to be
% typed in (2820, from the circulators' reverse ports), and the resonator does not stay
% put: it has since been measured near 2763, so a fixed default would have walked the
% field to a frequency the resonator had left.
    R = fit.resonance;
    out.target = struct('MHz', NaN, 'err', NaN, 'source', 'none', 'outside', false);
    if ~isempty(opts.targetMHz)
        out.target.MHz = opts.targetMHz; out.target.source = 'given';
    elseif R.found
        opts.targetMHz = R.centre;
        out.target = struct('MHz', R.centre, 'err', R.centreErr, 'source', 'resonance', ...
                            'outside', R.outside);
        fprintf('\n  coil target: the resonance at %.1f +- %.1f MHz (FWHM %.0f), measured in this sweep\n', ...
                R.centre, R.centreErr, R.fwhm);
        if R.outside
            fprintf(['  NOTE its centre lies OUTSIDE the sweep, extrapolated from the edge: extend\n' ...
                     '  the sweep past it to pin the target down.\n']);
        end
    end
    out.opts.targetMHz = opts.targetMHz;
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
    % CW window: AWGExpESR.loadJsonParams sets it (esrWindowPick) and prepare keeps
    % laserInitializationDuration short enough for it. detectionMHz overrides the window.
    if ~isempty(opts.detectionMHz)
        a.detectionDuration          = opts.detectionMHz;
        a.referenceDetectionDuration = opts.detectionMHz;
    elseif isempty(opts.jsonLocation)
        a.setCwWindow();
    end
    % Empty leaves the object's own value (true).
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
