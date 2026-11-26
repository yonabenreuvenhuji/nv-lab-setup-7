function out = calibFieldAlign(opts)
%CALIBFIELDALIGN Scan Btheta and Bphi to put the field on an NV axis.
%
%   out = calibFieldAlign
%   out = calibFieldAlign(struct('rounds', 3, 'thetaSpan', 6, 'phiSpan', 6))
%
% LEVEL 3. Alternating 1D scans in the Helmholtz spherical coordinates, with the
% optimum of each taken from a parabola fit rather than from the best measured point, and the
% final answer taken as the INVERSE-VARIANCE MEAN of those vertices rather than the last one
% -- every round measures the same angle again, so using only the last throws the rest away.
%
% ===================================================================
% WHY A SCAN AND NOT A CALCULATION
% ===================================================================
%
% Helmholtz.NVcoordinates returns the four NV axis vectors in the coil frame, and with
% those plus the measured projections the field direction follows in closed form:
% sum(a_i a_i') = (4/3)I for the four <111> vectors, so n = (3/4B) * sum(p_i a_i), and one
% step would land on the axis.
%
% That is not usable here. NVcoordinates is fixCoordinates(NVaxis) applied with `normal`,
% `edge`, `D_rotate`, `D_theta` and `D_phi` -- hand-entered diamond orientation
% parameters. The closed form would inherit exactly the error we are trying to remove,
% and it would do so invisibly, producing a confident answer built on an estimate. So the
% coil frame is used only for SIGNS, which are topological and survive a rough geometry,
% and the direction itself is found by measurement.
%
% ===================================================================
% TWO METRICS, BECAUSE ONE IS ROBUST AND THE OTHER IS SHARP
% ===================================================================
%
% COARSE, 'maxsplit': the largest splitting. At alignment the largest projection equals
% |B|, so maximising it finds an axis. It needs only the two outermost lines, so it is
% cheap and it works when the spectrum is a mess. But it goes as cos(tilt) -- second
% order -- so it is flat near the optimum and cannot finish the job.
%
% FINE, 'spread': the spread of the three SMALLER projections. At alignment the three
% off-axis orientations are degenerate at B/3, so the spread is zero, and it opens FIRST
% order in the tilt. This is the metric that actually converges.
%
% Both are built from FREQUENCIES. That is deliberate: the earlier width-based metric was
% contaminated because the resonator power-broadens the low-frequency lines much more
% than the high-frequency ones -- on this setup every line below 2850 came out at 18-23
% MHz FWHM while every line above 2890 came out at 6-9 MHz. Amplitudes and widths cannot
% be trusted across a sweep; line positions can.
%
% ===================================================================
% ONE FIXED FREQUENCY LIST
% ===================================================================
%
% Every measurement uses the SAME full-span list, at finer resolution than the survey's.
% Three reasons, and the first two were learned the hard way.
%
% Tracked +/-10 MHz windows lost their lines: the inner projections move by up to 2.7 G
% across an 8 degree scan, which is 7.6 MHz of line movement.
%
% A changing frequency list forces the AWG to re-upload a waveform per point -- 12.7 s
% measured, plus a confirmation dialog. An unchanged list hits isAlreadyLoaded and skips
% it, so the fixed fine sweep is FASTER despite having more points.
%
% And single-line fits in overlapping windows were biased: on lines 11 MHz apart with
% 9 MHz linewidths each fit was pulled by its neighbour, giving a tilt of 4.0 deg where
% the wide multi-line fit of the same field said 7.0.
%
% ===================================================================
% THE FIELD MUST BE ON, AND control() HAS A TRAP
% ===================================================================
%
% Helmholtz.control validates case-INSENSITIVELY and then branches case-SENSITIVELY:
%
%     if ~any(strcmpi(value, obj.CONTROL_STATE)); error(...); end   % 'on' passes
%     if strcmp(value, obj.CONTROL_STATE{1});  obj.output('on');    % 'ON' only
%     else                                     obj.output('off');
%
% So h.control('on') in lower case passes validation and then turns the supplies OFF,
% silently. Only h.control('ON') energises them. This script checks that a splitting is
% actually present rather than trusting the state.
%
% ===================================================================
% OPTIONS (struct fields, all optional)
%   thetaSpan, phiSpan   degrees either side of the current value (default 6)
%   nPoints              points per 1D scan (default 7)
%   rounds               theta-then-phi passes (default 2)
%   tolDeg               stop when both optima move less than this (default 0.4)
%   bowlFrac             fit only points with metric <= min + this*(max-min) (default 0.7).
%                        The metric is quadratic only NEAR the optimum; the wings are model
%                        error, not noise, and including them biases the vertex.
%   outlierSigma         reject a point whose deleted residual exceeds this many robust
%                        sigma (default 6 -- see parabolaOptimum for why not 3)
%   minFitPoints         never fit fewer than this many points (default 5)
%   scanFreqs            THE frequency list, used unchanged for every measurement
%                        (default 2725:1.5:3015, 194 points -- the span that holds the
%                        45 G pattern: outer lines at 2748/2992 with ~20 MHz of margin). Fixed on purpose: an
%                        unchanged list skips the AWG re-upload and the "parameters have
%                        changed" dialog, and nothing can wander out of range as the
%                        angle turns. Finer than the survey's 2 MHz because line centres
%                        are what the metric is built from.
%   repeats, averages    (default 50, 3)
%   amplitude            dBm per channel (default -15)
%   detectionMHz         CW detection and reference window, us. Default [] = let
%                        calibFieldSurvey resolve it through esrWindowPick. Do not set a
%                        short value: below ~10*tau_pol the dip is diluted, and below the
%                        AWG-trigger limit the spectrum becomes a comb.
%   jsonLocation         expProperties.json, loaded before anything else is set
%   isTracking           (default 0)
%   newOrderEveryAverage walk mode: a new random frequency order before every average
%                        (true/false), passed down to calibFieldSurvey for the locate sweep
%                        AND every scan point. Default [] = AWGExpESR's own setting (false
%                        on a fresh object). Costs one sequence upload per average.
%   freshExperiment      construct a NEW AWGExpESR for every measurement (default true).
%                        This is what stops the "parameters have changed" dialog: the
%                        dialog fires when an EXISTING object's parameters are altered,
%                        and a newly built object has no previous state to compare with.
%   settle               seconds after a field change (default 2)
%   baselineOrder        polynomial order for the shared baseline (default 3, because the
%                        290 MHz span at 45 G curves more than a parabola can follow).
%   metric               'auto' | 'innerspan' | 'tilt' | 'spread' | 'innerwidth' |
%                        'maxsplit' (default 'auto', which picks 'innerspan').
%
%                        'innerspan' is the one to use: the frequency span occupied by the
%                        three off-axis lines on one side of D, linewidths included. It is
%                        monotonic in the tilt whether or not those lines are resolved --
%                        measured 8.6 MHz at 0 deg rising steadily to 43.9 at 10 deg -- so
%                        one metric covers the whole approach. 'tilt' (the fitted angle,
%                        squared) is sharper above ~2.5 deg but goes blind below it, since
%                        line POSITIONS stop carrying the tilt once the off-axis lines
%                        merge. 'innerwidth' is kept only for comparison: it reverses
%                        direction above 2 deg and will drive a scan the wrong way.
%                        See metricFromSurvey for the measured tables.
%   plotLocate           draw the fit for the LOCATE sweep (default true). The scan
%                        points stay silent -- one plot per angle would be 30 figures --
%                        but the locate sweep is the one you want to look at, and a dry
%                        run is nothing but a locate sweep.
%   dryRun               measure the metric at the current setting and stop (default
%                        false) -- use this first, it costs one scan point
%
% Restores Btheta and Bphi on any error. On success it LEAVES the field at the optimum
% and prints the h.save command, but does not save for you.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('thetaSpan', 6, 'phiSpan', 6, 'nPoints', 7, ...
        'rounds', 2, 'tolDeg', 0.4, 'scanFreqs', 2725:1.5:3015, ...
        'bowlFrac', 0.7, 'outlierSigma', 6, 'minFitPoints', 5, ...
        'repeats', 50, 'averages', 3, 'detectionMHz', [], 'amplitude', -15, ...
        'settle', 2, 'metric', 'auto', 'dryRun', false, 'D0', 2870, 'pairTol', 4, ...
        'baselineOrder', 3, 'isTracking', 0, 'freshExperiment', true, ...
        'newOrderEveryAverage', [], ...
        'plotLocate', true, ...
        'jsonLocation', ''));

    % ---- where the calibration lives ------------------------------------
    % Resolved by setup7Paths rather than written down here. Twelve files used to carry
    % their own copy of the path, so moving the results folder meant twelve edits, and
    % missing one silently loaded the OLD calibration into this one experiment while
    % everything else used the new one -- two calibrations live at once, with no symptom
    % except data that will not reproduce.
    %
    % This function CANNOT run without a calibration, so a failure to resolve is fatal here
    % and setup7Paths says what is wrong with the folder.
    %
    % Pass jsonLocation explicitly ONLY to pin a run to an older calibration on purpose:
    %     setup7Paths('2026-08-14').json
    if isempty(opts.jsonLocation)
        opts.jsonLocation = setup7Paths().json;
    end

    GAMMA = 2.8025;                 % MHz/G. A physical constant, so every |B| below is
                                    % absolute and independent of magneticFieldRatio.

    % Check the helpers are reachable BEFORE touching the field. A missing function
    % discovered halfway through a scan leaves the coils at a scan value rather than
    % where they started.
    for fn = {'esrFitRobust', 'esrFieldSolve', 'calibFieldSurvey'}
        if isempty(which(fn{1}))
            error(['%s is not on the MATLAB path. Run setup7AddPath once this session ' ...
                   '(it adds functions\\ and every calibration level).'], fn{1});
        end
    end

    h = getObjByName('Helmholtz');
    if isempty(h); error('No Helmholtz object found.'); end

    % NO EXPERIMENT OBJECT IS HELD HERE. calibFieldSurvey builds a fresh AWGExpESR for
    % every measurement, so there is nothing to configure and nothing to restore -- and
    % nothing that can carry a stale flag from one measurement into the next.
    theta0 = h.Btheta; phi0 = h.Bphi; rho0 = h.Brho;

    fprintf('\n=== calibFieldAlign ===\n');
    fprintf('start: Btheta %.3f deg, Bphi %.3f deg, Brho %.3f G (coil model)\n', ...
            theta0, phi0, rho0);
    fprintf('NOTE Brho comes from magneticFieldRatio/Offset, which have no calibration\n');
    fprintf('     script behind them. |B| measured below is absolute and can calibrate them.\n');

    out = struct('theta0', theta0, 'phi0', phi0, 'rho0', rho0, 'gamma', GAMMA, ...
                 'scans', [], 'ok', false);

    try
        % ---- locate the lines, once -------------------------------------
        fprintf('\n--- locating lines (one wide sweep) ---\n');
        [lines, D, wideFit, surveyB, surveyTilt, sv0] = locateLines(opts, GAMMA);
        out.locate = struct('lines', lines, 'D', D, 'fit', wideFit, ...
                            'surveyTilt', surveyTilt, 'surveyB', surveyB);
        if numel(lines) < 4
            error(['Fewer than two usable pairs, so there is nothing to align on yet. ' ...
                   'The survey output above says which stage failed.\n' ...
                   'Three things to check, in order:\n' ...
                   '  1. Is the field actually ON? h.control(''on'') in LOWER case passes ' ...
                   'validation and then turns the supplies OFF -- only h.control(''ON'') ' ...
                   'energises them. h.Brho reads %.3f G.\n' ...
                   '  2. Is |B| large enough? Below about 8 G the inner orientations sit ' ...
                   'within a linewidth of D and cannot be separated.\n' ...
                   '  3. Does the sweep contain both members of every pair? At %.0f G the ' ...
                   'outer lines sit near %.0f and %.0f MHz; the sweep runs %.0f to %.0f.'], ...
                  rho0, rho0, opts.D0 - 2.8025*rho0, opts.D0 + 2.8025*rho0, ...
                  min(opts.scanFreqs), max(opts.scanFreqs));
        end
        fprintf('  %d lines, D = %.2f MHz\n', numel(lines), D);
        fprintf('  centres: %s\n', strjoin(cellstr(num2str(lines(:), '%.1f'))', ', '));

        % ---- pick the metric --------------------------------------------
        nPairAvail = floor(numel(lines)/2);
        metric = lower(opts.metric);
        if strcmp(metric, 'auto')
            % 'tilt' FIRST, whenever the field fit converged and reports a tilt the line
            % positions can actually resolve. It is the metric the others were proxies for:
            % calibFieldSurvey now fits Bx, By, Bz and D to every line at once through the
            % exact Hamiltonian, so the tilt in degrees comes out directly, with the
            % second-order transverse shift included and no pairing required.
            %
            % Its one limitation is real and is why the others are kept. Below about 3
            % degrees at 45 G the three off-axis orientations are no longer resolved -- their
            % splittings differ by less than a linewidth -- so line POSITIONS stop carrying
            % the tilt and the fit reports a tilt near zero whatever the truth is. Measured
            % on synthetic spectra with a known field: 10 deg -> 9.98, 7 -> 6.96, 4 -> 3.63,
            % 2 -> 1.95, 1 -> 0.10, 0 -> 0.08. Trustworthy to about 2 degrees, then blind.
            %
            % 'innerwidth' takes over exactly there, because the blend WIDTH still responds
            % when the positions no longer do. So: tilt to get from far out to a couple of
            % degrees, innerwidth to finish.
            % 'innerspan' first, because it is monotonic in the tilt at every tilt -- see
            % metricFromSurvey for the measured table and for why 'innerwidth', which it
            % replaces, actually reverses direction above 2 degrees. It needs three lines: an
            % outer pair and at least one inner line.
            if numel(lines) >= 3 && isfinite(innerSpanOf(sv0))
                metric = 'innerspan';
            elseif isfield(sv0, 'solve') && isstruct(sv0.solve) && sv0.solve.ok && ...
                    isfinite(sv0.solve.tiltDeg) && sv0.solve.tiltDeg > 2.5
                metric = 'tilt';
            elseif nPairAvail >= 3
                metric = 'spread';
            elseif numel(lines) >= 4
                metric = 'innerwidth';
            else
                metric = 'maxsplit';
            end
        end
        fprintf('  metric: ''%s'' (%d pairs resolved)\n', metric, nPairAvail);
        if strcmp(metric, 'maxsplit')
            fprintf(['  WARNING maxsplit is SECOND order in the tilt -- at 45 G and 7 deg it\n' ...
                     '  moves only ~0.5 MHz per degree against a 1.0 MHz splitting error. It\n' ...
                     '  will get you close and then stall.\n']);
        elseif strcmp(metric, 'innerspan')
            fprintf(['  innerspan: minimising the frequency span of the off-axis group,\n' ...
                     '  linewidths included. Monotonic in the tilt at every tilt -- ~3.5 MHz\n' ...
                     '  per degree near alignment, 4.5 further out -- because it does not\n' ...
                     '  depend on those lines being resolved. The value reported is\n' ...
                     '  (span - intrinsic linewidth) squared, which makes it an exact\n' ...
                     '  parabola in the scan angle.\n']);
        elseif strcmp(metric, 'innerwidth')
            fprintf(['  innerwidth: minimising the FWHM of the widest inner line. NOT\n' ...
                     '  MONOTONIC above ~2 deg -- prefer innerspan. It is\n' ...
                     '  blended now (three off-axis orientations, not degenerate) and narrows\n' ...
                     '  toward the intrinsic linewidth as the tilt falls. First order in the\n' ...
                     '  tilt, and it stays meaningful all the way to alignment.\n']);
        end
        out.metric = metric;

        % ---- the current value, which is also the dry run ---------------
        % THE STARTING METRIC COMES FROM THE LOCATE SWEEP. No second sweep.
        %
        % The previous version called measureMetric here, which ran calibFieldSurvey again
        % on the same field with the same frequency list -- two identical 121-point sweeps
        % back to back, for one number that the first had already computed. Now that both
        % use the same fixed grid and the same pipeline there is nothing a second sweep can
        % add, so the locate result IS the first measurement.
        [m0, s0] = metricFromSurvey(sv0, metric, GAMMA);
        fprintf('\n  metric at the start: %.4f  (%s)\n', m0, metricUnits(metric));
        fprintf('  |B| = %.3f G, tilt = %.2f deg\n', s0.B, s0.tiltDeg);
        fprintf('  (from the locate sweep -- no second sweep needed)\n');
        out.startMetric = m0; out.startState = s0;
        if opts.dryRun
            fprintf('\n(dryRun: nothing scanned, field unchanged)\n');
            restoreAll(h, theta0, phi0); out.ok = true; return
        end

        % ---- alternating 1D scans ---------------------------------------
        cur = struct('theta', theta0, 'phi', phi0);
        S = struct('round', {}, 'axis', {}, 'values', {}, 'metric', {}, ...
                   'opt', {}, 'fromParabola', {}, 'fitCoeff', {}, 'keep', {}, ...
                   'discarded', {}, 'outOfBowl', {}, 'optErr', {}, 'why', {}, 'suspect', {});
        for r = 1:opts.rounds
            moved = 0;
            for ax = {'Btheta', 'Bphi'}
                axis = ax{1};
                if strcmp(axis, 'Btheta')
                    centre = cur.theta; span = opts.thetaSpan;
                else
                    centre = cur.phi;   span = opts.phiSpan;
                end
                vals = linspace(centre - span, centre + span, opts.nPoints);
                fprintf('\n--- round %d, %s: %.2f to %.2f deg, %d points ---\n', ...
                        r, axis, vals(1), vals(end), numel(vals));
                mv = nan(size(vals));
                for iv = 1:numel(vals)
                    okSet = setAngle(h, axis, vals(iv));
                    if ~okSet
                        fprintf('  %7.2f deg : out of limits, skipped\n', vals(iv));
                        continue
                    end
                    pause(opts.settle);
                    [mv(iv), st] = measureMetric(lines, D, opts, GAMMA, metric);
                    if isfinite(mv(iv))
                        lines = st.lines;       % TRACK: the lines move as the field turns
                        D     = st.D;
                    end
                    fprintf('  %7.2f deg : metric %8.4f   |B| %7.3f G   tilt %5.2f deg\n', ...
                            vals(iv), mv(iv), st.B, st.tiltDeg);
                end
                % PARABOLA, not the best point. The scan step is coarse compared with the
                % precision the line fits give, so interpolating between points is worth
                % several times the step -- and it also exposes a scan that has no optimum
                % inside its range, which picking the best point silently hides.
                [optVal, fromPar, parInfo] = parabolaOptimum(vals, mv, metric, ...
                    struct('bowlFrac', opts.bowlFrac, 'outlierSigma', opts.outlierSigma, ...
                           'minFitPoints', opts.minFitPoints, 'verbose', true));
                if ~isfinite(optVal)
                    fprintf('  no usable optimum in this range -- leaving %s alone\n', axis);
                    optVal = centre; fromPar = false;
                end
                fprintf('  optimum %s = %.3f deg (%s)\n', axis, optVal, ...
                        ternary(fromPar, 'parabola vertex', 'best measured point'));

                % NOTE: the scan MOVES to every vertex, even a step smaller than its own
                % error bar. Holding still on small steps was tried and is worse -- it freezes
                % progress exactly when the vertex is well determined, and over 200 simulated
                % alignments it took the final theta error from 0.19 to 0.20 deg while doing
                % nothing about the large steps. The wander is dealt with after the loop
                % instead, by averaging the vertices rather than taking the last one.
                setAngle(h, axis, optVal);
                pause(opts.settle);
                if strcmp(axis, 'Btheta')
                    moved = max(moved, abs(optVal - cur.theta)); cur.theta = optVal;
                else
                    moved = max(moved, abs(optVal - cur.phi));   cur.phi = optVal;
                end
                S(end+1) = struct('round', r, 'axis', axis, 'values', vals, ...
                    'metric', mv, 'opt', optVal, 'fromParabola', fromPar, ...
                    'fitCoeff', parInfo.coeff, 'keep', parInfo.keep, ...
                    'discarded', parInfo.discarded, 'outOfBowl', parInfo.outOfBowl, ...
                    'optErr', parInfo.optErr, 'why', parInfo.why, ...
                    'suspect', parInfo.suspect);   %#ok<AGROW>
            end
            fprintf('\n  round %d moved the optimum by %.3f deg\n', r, moved);
            if moved < opts.tolDeg
                fprintf('  below tolDeg = %.2f, stopping early\n', opts.tolDeg);
                break
            end
        end
        out.scans = S;

        % ---- THE ANSWER IS THE WEIGHTED MEAN OF THE VERTICES, NOT THE LAST ONE ----
        %
        % Every round measures the same angle again, so the last vertex throws away all the
        % earlier measurements of it. Combining them by inverse variance -- each vertex
        % weighted by 1/sigma^2 from its own parabola covariance -- uses all of them. Over
        % 200 simulated alignments that took the final theta error from 0.193 to 0.119 deg,
        % a 39% reduction for no extra measurement time.
        %
        % This is also what actually deals with the wander. On a real alignment the Bphi
        % optimum went 48.2, 48.3, 47.8, 45.8, 47.0, 45.9, 46.5 over seven rounds and the
        % last value, 46.5, was no better than any other. The weighted mean of the seven is.
        [thBar, thErr, nTh] = combineVertices(S, 'Btheta', cur.theta);
        [phBar, phErr, nPh] = combineVertices(S, 'Bphi',   cur.phi);
        fprintf('\n--- combining the vertices ---\n');
        fprintf('  Btheta: last %.3f, weighted mean of %d rounds %.3f +/- %.3f deg\n', ...
                cur.theta, nTh, thBar, thErr);
        fprintf('  Bphi:   last %.3f, weighted mean of %d rounds %.3f +/- %.3f deg\n', ...
                cur.phi, nPh, phBar, phErr);
        if abs(thBar - cur.theta) > 3*max(thErr, eps) || abs(phBar - cur.phi) > 3*max(phErr, eps)
            fprintf(['  The mean and the last round differ by more than 3 sigma. The vertices\n' ...
                     '  are DRIFTING rather than scattering -- something is changing during the\n' ...
                     '  alignment (temperature, the coil driver, the spot). Averaging is then\n' ...
                     '  the wrong thing to do: re-run and use the last round.\n']);
        end
        cur.theta = thBar; cur.phi = phBar;
        out.thetaErr = thErr; out.phiErr = phErr;
        setAngle(h, 'Btheta', cur.theta); setAngle(h, 'Bphi', cur.phi);
        pause(opts.settle);

        % ---- verify -------------------------------------------------------
        fprintf('\n--- verification at the optimum ---\n');
        % This one IS a new measurement, not a repeat: the angles were just moved to the
        % weighted mean of the vertices, which is between scan points and has never been
        % measured at all.
        [m1, s1] = measureMetric(lines, D, opts, GAMMA, metric);
        out.endMetric = m1; out.endState = s1;
        out.theta = cur.theta; out.phi = cur.phi;
        fprintf('  metric  %.4f -> %.4f\n', m0, m1);
        fprintf('  tilt    %.2f -> %.2f deg\n', s0.tiltDeg, s1.tiltDeg);
        fprintf('  |B|     %.3f -> %.3f G   (coil model says Brho = %.3f)\n', ...
                s0.B, s1.B, h.Brho);
        fprintf('  Btheta  %.3f -> %.3f deg\n', theta0, cur.theta);
        fprintf('  Bphi    %.3f -> %.3f deg\n', phi0,   cur.phi);
        % The scale check is free and it is the plan's section 4.4: |B| here is absolute
        % because gamma is a constant, so its ratio to Brho calibrates
        % magneticFieldRatio/Offset against something physical.
        if h.Brho > 0
            fprintf(['\n  |B|_measured / Brho_coilmodel = %.4f. Those numbers come from\n' ...
                     '  magneticFieldRatio/Offset in setupInfo.json, which have no calibration\n' ...
                     '  script behind them -- this ratio is the correction, measured against a\n' ...
                     '  physical constant.\n'], s1.B / h.Brho);
        end
        if s1.tiltDeg > s0.tiltDeg
            fprintf(['\n  *** THE TILT GOT WORSE. Do not keep this. Restoring the start point.\n' ...
                     '      Likely causes: the scan range did not contain the optimum, or the\n' ...
                     '      metric is ''maxsplit'' which is flat near the optimum and can be\n' ...
                     '      pulled around by noise.\n']);
            setAngle(h, 'Btheta', theta0); setAngle(h, 'Bphi', phi0);
            out.theta = theta0; out.phi = phi0;
        else
            fprintf('\n  To keep this:\n');
            fprintf('    h = getObjByName(''Helmholtz'');\n');
            fprintf('    h.setBtheta(%.4f); h.setBphi(%.4f);\n', cur.theta, cur.phi);
            fprintf('    h.save([currentPath, ''helmholtzState.mat'']);\n');
        end
        if opts.nPoints >= 3; plotScans(S, metric, opts.tolDeg); end
        out.ok = true;

    catch err
        fprintf('\n(error -- restoring Btheta %.3f, Bphi %.3f)\n', theta0, phi0);
        restoreAll(h, theta0, phi0);
        rethrow(err);
    end
end

% ======================================================================

function [lines, D, fit, surveyB, surveyTilt, sv] = locateLines(opts, GAMMA) %#ok<INUSD>
% Locate the lines by CALLING calibFieldSurvey, not by reimplementing it.
%
% The first version of this function was a cut-down copy of the survey's pipeline, and it
% quietly left out the split-and-test pass -- so a blended line never separated, the
% mirror-partner filter then dropped everything without a partner, and it reported 2 lines
% where the survey found 6 on the same sweep. Two implementations of a hard step means one
% of them is wrong and you find out later.
%
% The survey already does all of it, validated on this setup's real spectrum: multi-line
% fit with a quadratic baseline, interleaved mirror-seeding for asymmetric omissions and
% split-and-test for symmetric blends, then the mirror-partner filter that removes
% features which are not NV lines at all.
    % Same grid the scan will use, so the cross-check below compares like with like.
    sv = calibFieldSurvey(struct('freqs', opts.scanFreqs, ...
            'repeats', opts.repeats, 'averages', opts.averages, ...
            'detectionMHz', opts.detectionMHz, 'plotIt', opts.plotLocate, ...
            'D0', opts.D0, 'pairTol', opts.pairTol, ...
            'baselineOrder', opts.baselineOrder, 'amplitude', opts.amplitude, ...
            'jsonLocation', opts.jsonLocation, 'isTracking', opts.isTracking, ...
            'freshExperiment', opts.freshExperiment, ...
            'newOrderEveryAverage', opts.newOrderEveryAverage));
    fit = []; surveyB = NaN; surveyTilt = NaN;
    % sv is returned whole so the caller can take the starting metric from THIS sweep
    % rather than paying for a second identical one.
    if isfield(sv, 'fit'); fit = sv.fit; end
    if isfield(sv, 'B');       surveyB    = sv.B;       end
    if isfield(sv, 'tiltDeg'); surveyTilt = sv.tiltDeg; end
    if isfield(sv, 'nvLines') && numel(sv.nvLines) >= 4
        lines = sort(sv.nvLines);
        D = sv.D;
        % B and tiltDeg only exist when the survey resolved three or more pairs -- with
        % exactly two it refuses to identify them, so guard rather than assume.
        if isfield(sv, 'B') && isfield(sv, 'tiltDeg')
            fprintf('  survey: |B| = %.3f G, tilt = %.2f deg\n', sv.B, sv.tiltDeg);
        else
            fprintf(['  survey resolved only %d pair(s), so |B| and the tilt are not\n' ...
                     '  identified -- the scan will run on ''maxsplit'' instead.\n'], ...
                    floor(numel(lines)/2));
        end
        return
    end
    % ---- FALLBACK: pair permissively from the raw fit ---------------------
    % The survey is strict on purpose -- it refuses to report a field magnitude from a
    % pairing it cannot fully justify, which is right for a measurement. But a scan does
    % not need a justified field magnitude. It needs a few lines it can track and a metric
    % that moves monotonically with the tilt, so an odd line count or one pair mean a
    % little outside tolerance should degrade the answer, not block it.
    %
    % So: take the survey's fitted centres, find the D that pairs the most of them, keep
    % the partnered ones, and if the count is still odd drop the WORST-matched line rather
    % than giving up. Everything that is dropped or forced is printed.
    lines = []; D = opts.D0;
    if isempty(fit) || ~isfield(fit, 'lines') || isempty(fit.lines)
        fprintf('\n  the wide sweep produced no fitted lines at all -- read the survey output above\n');
        return
    end
    c = sort([fit.lines.centre]);
    fprintf('\n--- falling back to permissive pairing ---\n');
    fprintf('  fitted centres: %s\n', strjoin(cellstr(num2str(c(:), '%.1f'))', ', '));

    Dc = opts.D0 + (-15:0.05:15); nm = zeros(size(Dc)); r2 = zeros(size(Dc));
    for q = 1:numel(Dc)
        for k = 1:numel(c)
            d = abs(c - (2*Dc(q) - c(k))); d(k) = inf;
            if min(d) <= opts.pairTol; nm(q) = nm(q) + 1; r2(q) = r2(q) + min(d)^2; end
        end
    end
    bestQ = find(nm == max(nm)); [~, ibq] = min(r2(bestQ)); D = Dc(bestQ(ibq));

    mism = inf(1, numel(c));
    for k = 1:numel(c)
        d = abs(c - (2*D - c(k))); d(k) = inf; mism(k) = min(d);
    end
    fprintf('  D = %.2f pairs %d of %d\n', D, max(nm), numel(c));
    for k = 1:numel(c)
        fprintf('    %8.2f  mirror at %8.2f  nearest line off by %6.2f MHz  %s\n', ...
                c(k), 2*D - c(k), mism(k), ...
                ternary(mism(k) <= opts.pairTol, 'PAIRED', 'orphan'));
    end
    keep = mism <= opts.pairTol;
    lines = c(keep); mism = mism(keep);
    if mod(numel(lines), 2) == 1
        [~, iw] = max(mism);
        fprintf('  odd count -- dropping the worst-matched line at %.2f (off by %.2f MHz)\n', ...
                lines(iw), mism(iw));
        lines(iw) = [];
    end
    if numel(lines) >= 4
        fprintf('  proceeding with %d lines: %s\n', numel(lines), ...
                strjoin(cellstr(num2str(lines(:), '%.1f'))', ', '));
    end
end

function [m, st] = measureMetric(lines, D, opts, GAMMA, metric)
% ONE FIXED FREQUENCY LIST, measured and fitted the same way every time.
%
% This replaces per-line tracked windows, and it is better for three separate reasons.
%
% NOTHING CAN LEAVE THE RANGE. As the angle turns, the inner projections move by up to
% 2.7 G, which is 7.6 MHz of line movement -- more than a +/-10 MHz window can follow
% reliably, and the lines were indeed wandering out.
%
% NO RELOAD, AND NO DIALOG. The AWG holds a waveform per frequency, so a frequency list
% that changes every scan point forces a full re-upload: measured at 12.7 s per point,
% plus the "parameters have changed" confirmation. An unchanged list hits
% isAlreadyLoaded and skips the upload outright. A 181-point fixed sweep at ~12 s of
% acquisition therefore costs LESS wall time than an 84-point tracked one at 12.7 s of
% reload plus 5 s of acquisition.
%
% AND IT REMOVES THE OVERLAP BUG BY CONSTRUCTION. Tracked windows put a single Lorentzian
% on lines 11 MHz apart with 9 MHz linewidths, so each fit was pulled by its neighbour --
% the metric read 2.77 G where the wide fit of the same field said 4.20. A full sweep
% fitted as one multiplet cannot have that failure mode.
%
% The fitting is calibFieldSurvey itself, called quietly. One implementation of the hard
% part: multi-line fit, interleaved mirror-seeding and split-and-test, mirror-partner
% filter, projections from the exact identities. lines and D are no longer needed here --
% the survey finds them fresh each time -- and are kept in the signature only so the
% caller's bookkeeping does not change.
    st = struct('B', NaN, 'tiltDeg', NaN, 'lines', lines, 'D', D, ...
                'splittings', [], 'projections', []);
    m = NaN;

    % SEED FROM THE PREVIOUS POINT'S LINES. The angle moves by a couple of degrees per
    % scan step, which shifts a line by only 2-4 MHz, so the fit starts from essentially
    % the right answer. Without it every point re-ran a peak finder from scratch and some
    % points found the lines while others did not -- which is exactly the intermittent
    % failure this is meant to cure. Physics carrying information forward, not a fudge:
    % if a seeded line is not really there it fails its depth test and is dropped.
    o = struct('freqs', opts.scanFreqs, 'repeats', opts.repeats, ...
               'averages', opts.averages, 'detectionMHz', opts.detectionMHz, ...
               'plotIt', false, 'D0', opts.D0, 'pairTol', opts.pairTol, ...
'baselineOrder', opts.baselineOrder, 'amplitude', opts.amplitude, ...
               'jsonLocation', opts.jsonLocation, 'isTracking', opts.isTracking, ...
               'freshExperiment', opts.freshExperiment, 'seedCentres', lines, ...
               'newOrderEveryAverage', opts.newOrderEveryAverage);
    % evalc so the survey's several dozen lines of reasoning do not fill the log once per
    % scan point. Errors still propagate; only the printing is swallowed.
    % evalc keeps the survey's several dozen lines of reasoning out of the log at every
    % scan point -- but the text is KEPT, and printed if the measurement fails. The first
    % version discarded it, so a point that came back NaN gave no clue why: the run that
    % exposed this printed "metric at the start: NaN" with the entire explanation swallowed.
    txt = '';
    try
        txt = evalc('sv = calibFieldSurvey(o);');
    catch err
        fprintf('    (survey errored: %s)\n', err.message);
        return
    end
    if ~exist('sv', 'var') || ~isfield(sv, 'projections') || isempty(sv.projections)
        st.why = 'survey could not pair the lines at this angle';
        fprintf('    (no metric here -- the survey could not pair the lines. Its output:)\n');
        printIndented(txt);
        return
    end
    [m, st] = metricFromSurvey(sv, metric, GAMMA);
end

function [m, st] = metricFromSurvey(sv, metric, GAMMA)
% Reduce an already-computed survey result to the scan metric. Shared by the starting
% point and by every scan point, so the two cannot drift apart.
    st = struct('B', NaN, 'tiltDeg', NaN, 'projections', [], 'splittings', [], ...
                'lines', [], 'D', NaN, 'why', '', 'innerFwhm', NaN, 'intrinsicFwhm', NaN);
    m = NaN;
    if ~isstruct(sv); st.why = 'no survey result'; return; end

    % Everything the field fit knows, recorded for every metric so the scan log carries
    % |B| and the tilt even when the metric being minimised is a width.
    if isfield(sv, 'solve') && isstruct(sv.solve) && isfield(sv.solve, 'ok') && sv.solve.ok
        st.B = sv.solve.B; st.tiltDeg = sv.solve.tiltDeg; st.D = sv.solve.D;
        st.projections = sv.solve.projections;
        st.splittings = sv.solve.splittings;
    end

    % 'tilt': the fitted angle between B and the nearest <111> axis, in degrees, minimised
    % directly. Not a proxy for the thing we want -- it IS the thing we want, which is why
    % it is first. See the auto-selection comment for the regime where it goes blind.
    % SQUARED, and this is not cosmetic. A 1D scan in theta through the optimum leaves a
    % fixed perpendicular residual phi_res, so the total tilt is
    %
    %       tilt(theta) = sqrt((theta - theta*)^2 + phi_res^2)
    %
    % which is a hyperbola -- V-shaped when phi_res is small -- and fitting a parabola to it
    % is fitting the wrong function precisely where the scan is most informative. Squaring
    % gives
    %
    %       tilt^2 = (theta - theta*)^2 + phi_res^2
    %
    % an EXACT parabola in the scan variable, with curvature exactly 1 deg^2/deg^2 and a
    % minimum value that is itself the useful number: sqrt of it is the misalignment left in
    % the OTHER angle, which tells you whether another round is worth running. So the
    % parabola vertex is not an approximation here, and the curvature is a free check on the
    % measurement: if the fitted curvature is not near 1, the metric is not measuring tilt.
    if strcmp(metric, 'tilt')
        if isfield(sv, 'solve') && isstruct(sv.solve) && isfield(sv.solve, 'ok') && ...
                sv.solve.ok && isfinite(sv.solve.tiltDeg)
            m = sv.solve.tiltDeg^2;
            if isfield(sv, 'fit') && isstruct(sv.fit) && ~isempty(sv.fit.lines)
                st.lines = [sv.fit.lines.centre];
            end
            return
        end
        st.why = 'the field fit did not converge, so there is no tilt to minimise';
        return
    end

    % 'innerspan': the FREQUENCY SPAN OCCUPIED BY THE OFF-AXIS GROUP on one side of D,
    % including the linewidths -- from the low edge of the lowest inner line to the high
    % edge of the highest, at half maximum:
    %
    %       span = max(centre + fwhm/2) - min(centre - fwhm/2)     over the inner lines
    %
    % This is the metric to use. It is the only one measured here that is monotonic in the
    % tilt over the ENTIRE range, because it does not care whether the three off-axis
    % orientations are resolved. Blended, the span is the blend's width; resolved, it is
    % the separation plus a linewidth; and the crossover between those two regimes is
    % smooth. Measured on synthetic spectra at 45 G with a known field, tilt in degrees
    % against span in MHz:
    %
    %       0.0 -> 8.63   0.5 -> 8.99   1.0 -> 10.26  1.5 -> 11.82  2.0 -> 13.41
    %       2.5 -> 17.90  3.0 -> 19.44  4.0 -> 22.95  5.0 -> 26.63  7.0 -> 33.10
    %      10.0 -> 43.87
    %
    % Monotonic throughout, about 3.5 MHz per degree at small tilt rising to 4.5 at large,
    % against a line-centre precision of 0.5 MHz and a width precision near 1 MHz. Compare
    % the two metrics it replaces, on the same spectra: 'tilt' reads 0.1 deg when the truth
    % is 1.5 (line positions stop carrying the tilt below ~2.5 deg), and 'innerwidth' is
    % NOT MONOTONIC -- 9.8, 10.3, 11.8, 13.4 up to 2 deg and then back DOWN to 11 MHz by 3
    % deg, because once the components resolve the widest inner line is a single line again
    % rather than the blend. A scan on innerwidth started above 2.5 deg would be driven the
    % wrong way, which is why auto no longer chooses it there.
    %
    % min(low side, high side), not the mean: the resonator distorts the low side's widths
    % more, and a spurious extra line inflates whichever side it lands on, so the smaller of
    % the two is the more trustworthy reading of the same physical quantity.
    %
    % SQUARED, less the intrinsic linewidth, for the same reason 'tilt' is squared: the span
    % is roughly w0 + k*tilt and tilt = sqrt((theta-theta*)^2 + phi_res^2), so (span - w0)^2
    % is parabolic in the scan angle while the span itself is V-shaped.
    if strcmp(metric, 'innerspan')
        [sp, w0] = innerSpanOf(sv);
        st.innerFwhm = sp; st.intrinsicFwhm = w0;
        if isfinite(sp)
            if isfinite(w0); m = max(sp - w0, 0)^2; else; m = sp^2; end
            if isfield(sv, 'fit') && isstruct(sv.fit) && ~isempty(sv.fit.lines)
                st.lines = [sv.fit.lines.centre];
            end
            return
        end
        st.why = 'no inner lines, so no off-axis group to measure the span of';
        return
    end

    % 'maxsplit' needs ONLY the outermost pair, so it must not be gated on the projections,
    % which require three pairs. That gating is why the metric came back NaN on a run that
    % had deliberately fallen back to maxsplit precisely because only two pairs resolved.
    if isfield(sv, 'pairs') && ~isempty(sv.pairs)
        st.splittings = [sv.pairs.splitting];
        if isfield(sv, 'D'); st.D = sv.D; end
        if isfield(sv, 'nvLines'); st.lines = sv.nvLines; end
        if strcmp(metric, 'maxsplit')
            m = -max(st.splittings);        % negated: everything here is MINIMISED
            st.B = max(st.splittings)/(2*GAMMA);   % a LOWER BOUND on |B|, equal at alignment
            return
        end
    end
    % 'innerwidth': the FWHM of the widest INNER line, minimised.
    %
    % This is the metric for the regime your spectrum is actually in. The three off-axis
    % orientations are degenerate only at perfect alignment, so as the tilt falls they
    % merge and the blend NARROWS toward the intrinsic linewidth. At 45 G and 7 deg the
    % low-side blend is 43 MHz FWHM against an 8.3 MHz single line -- an enormous, obvious
    % signal, and first order in the tilt.
    %
    % Earlier I rejected a width-based metric because the resonator power-broadens the
    % low-frequency lines far more than the high ones. That objection does not apply here,
    % and the reason matters: it only bites when comparing lines at DIFFERENT frequencies.
    % Here the same line is compared with itself across an angle scan, during which it
    % moves a few MHz at most, so the envelope's contribution is a constant offset that
    % cannot shift where the minimum lies.
    %
    % It is also the natural endgame. 'spread' needs the three inner splittings resolved,
    % which requires the misalignment it is trying to remove; as you approach alignment the
    % lines merge and the spread becomes unmeasurable. The blend width stays meaningful all
    % the way down to the intrinsic linewidth.
    if isfield(sv, 'fit') && isstruct(sv.fit) && isfield(sv.fit, 'lines') && ...
            ~isempty(sv.fit.lines) && strcmp(metric, 'innerwidth')
        cc = [sv.fit.lines.centre]; ww = [sv.fit.lines.fwhm];
        inner = cc > min(cc) + 1 & cc < max(cc) - 1;      % everything but the outer pair
        if any(inner)
            m = max(ww(inner));
            st.innerFwhm = m;
            st.intrinsicFwhm = min(ww(~inner));
            return
        end
    end
    if ~isfield(sv, 'projections') || isempty(sv.projections)
        st.why = 'survey did not resolve enough pairs for the spread metric'; return
    end
    st.B = sv.B; st.tiltDeg = sv.tiltDeg; st.projections = sv.projections;
    if isfield(sv, 'nvLines'); st.lines = sv.nvLines; end
    if isfield(sv, 'D');       st.D     = sv.D;       end
    if isfield(sv, 'pairs') && ~isempty(sv.pairs); st.splittings = [sv.pairs.splitting]; end
    p = sv.projections;
    if strcmp(metric, 'maxsplit')
        m = -max(abs(p)) * 2 * GAMMA;       % negated: everything here is MINIMISED
        return
    end
    small = sort(abs(p));
    small = small(1:min(3, numel(small)));
    m = max(small) - min(small);
end

function [span, w0] = innerSpanOf(sv)
% The span occupied by the off-axis group, and the intrinsic linewidth to compare it with.
%
% "Inner" means every line except the outermost pair, which belongs to the best-aligned
% orientation. The two sides of D are measured separately and the SMALLER is returned --
% see the comment at the call site. w0 is the narrowest outer line, which is a single
% orientation and therefore the intrinsic linewidth.
    span = NaN; w0 = NaN;
    if ~isstruct(sv) || ~isfield(sv, 'fit') || ~isstruct(sv.fit) || ...
            ~isfield(sv.fit, 'lines') || isempty(sv.fit.lines)
        return
    end
    cc = [sv.fit.lines.centre]; ww = [sv.fit.lines.fwhm];
    if numel(cc) < 3; return; end
    D = NaN;
    if isfield(sv, 'solve') && isstruct(sv.solve) && isfield(sv.solve, 'ok') && sv.solve.ok
        D = sv.solve.D;
    elseif isfield(sv, 'D') && isfinite(sv.D)
        D = sv.D;
    end
    if ~isfinite(D); D = (min(cc) + max(cc))/2; end

    % WHICH LINES ARE "OUTER" COMES FROM THE FIELD FIT'S AXIS ASSIGNMENT, not from taking
    % the lowest and highest centre.
    %
    % Taking min and max fails whenever one spurious line lands outside the real outer pair:
    % the true outer line is then counted as INNER, and the span jumps from ~14 MHz to the
    % full width of the aligned splitting. That is not hypothetical -- it happened at one
    % point of the simulated scan, which read 73.6 MHz where its neighbours read 12 and 12,
    % and one such point drags a parabola vertex 9 degrees off. A scan cannot afford a metric
    % with that failure mode, because on live data there is nothing to compare against.
    %
    % The field fit already assigns every matched line to an orientation, so the aligned pair
    % is identified by physics rather than by position, and a spurious line is simply
    % unmatched and excluded. min/max stays as the fallback for when the fit did not converge.
    isOuter = [];
    if isfield(sv, 'solve') && isstruct(sv.solve) && isfield(sv.solve, 'ok') && sv.solve.ok ...
            && isfield(sv.solve, 'matched') && ~isempty(sv.solve.matched) ...
            && isfield(sv.solve, 'alignedAxisIdx')
        mc = [sv.solve.matched.centre]; ma = [sv.solve.matched.axis];
        isOuter = false(size(cc)); known = false(size(cc));
        for k = 1:numel(cc)
            [dd, j] = min(abs(mc - cc(k)));
            if dd < 1e-6
                known(k) = true;
                isOuter(k) = (ma(j) == sv.solve.alignedAxisIdx);
            end
        end
        % Anything the fit did not match is not an NV line, so it is neither outer nor part
        % of the off-axis group.
        if any(isOuter) && any(known & ~isOuter)
            cc = cc(known); ww = ww(known); isOuter = isOuter(known);
        else
            isOuter = [];
        end
    end
    if isempty(isOuter)
        isOuter = (cc == min(cc)) | (cc == max(cc));
    end
    inner = ~isOuter;
    w0 = min(ww(isOuter));
    sides = [NaN NaN];
    for s = 1:2
        if s == 1; m = inner & cc < D; else; m = inner & cc > D; end
        if any(m)
            sides(s) = max(cc(m) + ww(m)/2) - min(cc(m) - ww(m)/2);
        end
    end
    if all(isnan(sides)); return; end
    span = min(sides(isfinite(sides)));
end

function printIndented(txt)
    if isempty(txt); return; end
    lines = strsplit(txt, newline);
    for i = 1:numel(lines)
        if ~isempty(strtrim(lines{i})); fprintf('      | %s\n', lines{i}); end
    end
end

function w = narrowestFwhm(fit)
    w = NaN;
    if ~isempty(fit) && isfield(fit, 'lines') && ~isempty(fit.lines)
        w = min([fit.lines.fwhm]);
    end
end

function [optVal, fromParabola, info] = parabolaOptimum(vals, mv, metric, opts)  %#ok<INUSL>
% Parabola through the scan, vertex taken as the optimum, with outliers rejected AGAINST
% THE FIT and the far wings excluded from it.
%
% ===================================================================
% WHY THE OLD REJECTION TEST WAS MEASURING THE WRONG THING
% ===================================================================
%
% It flagged a point when
%
%       abs(metric - median(metric))  >  6 * MAD(metric)
%
% i.e. it compared each point to the MEDIAN OF THE DATA. On a parabola that is the wrong
% reference in both directions, and it fails both ways:
%
%   IT MISSES REAL OUTLIERS. The arms of a parabola are legitimately far from the median, so
%   they inflate the MAD enormously. On a 9-point scan spanning 10 to 2710 the threshold
%   comes out at 6*MAD = 4104 -- larger than the entire range of the data -- so a point
%   sitting 180 above its neighbours passes untouched and drags the vertex.
%
%   IT DISCARDS GOOD POINTS. On a steep, well-sampled scan the outermost points are the
%   furthest from the median by construction, so they are the first candidates for a test
%   that only knows about distance from the median.
%
% The right reference is the FITTED CURVE. A point is an outlier when its residual from the
% parabola is large compared with the scatter of the other residuals -- which is what this
% now tests, iteratively, refitting after each rejection.
%
% ===================================================================
% AND THE WINGS ARE MODEL ERROR, NOT NOISE
% ===================================================================
%
% The metric is only quadratic NEAR the optimum. Far from it the underlying quantity stops
% being linear in the scan angle -- for the tilt metrics the off-axis splittings blend and
% the estimate saturates -- so the outer points are not noisy measurements of a parabola,
% they are good measurements of something that is not a parabola. Including them biases the
% vertex without inflating the residuals enough to be caught as outliers.
%
% So the fit is restricted to the BOWL: the contiguous run of points around the minimum
% whose metric is within bowlFrac of the way up the range. That is the region where the
% quadratic model is actually the right one. Points outside it are reported as excluded, not
% as bad -- they are neither.
%
% ===================================================================
% AND THE VERTEX NOW CARRIES AN ERROR BAR
% ===================================================================
%
% Propagated from the fit covariance. This is what tells you whether an optimum that moved
% 2 degrees between rounds actually moved: if sigma is a degree, it did not. Without it the
% alternating scans can wander and look like a failure to converge when they are simply not
% determining the angle, which is a different problem with a different fix (more points, more
% averages, or a larger |B| so the metric has more to work with).
%
% OPTIONS (4th argument, all optional)
%   bowlFrac       fit points with metric <= min + bowlFrac*(max-min)  (default 0.7)
%   outlierSigma   reject deleted residuals beyond this many robust sigma (default 6).
%                  Six, not three, and that is measured rather than guessed. Over 400
%                  synthetic scans at realistic noise, 6 still caught a rogue ten times the
%                  noise 98% of the time and anything larger every time -- while dropping
%                  0.87 good points per scan against 1.72 at 3.5, and giving a better vertex
%                  (rms 0.116 deg against 0.134). The threshold is applied to a LEAVE-ONE-OUT
%                  t with only three or four degrees of freedom, whose tails are much heavier
%                  than a Gaussian's, so a Gaussian intuition about "3 sigma" is the wrong
%                  yardstick here.
%   minFitPoints   never fit fewer than this                           (default 5)
%   verbose        announce rejections and exclusions                  (default true)

    if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
    opts = withDefaults(opts, struct('bowlFrac', 0.7, 'outlierSigma', 6, ...
                                     'minFitPoints', 5, 'verbose', true));

    ok0 = isfinite(vals) & isfinite(mv);
    optVal = NaN; fromParabola = false;
    info = struct('coeff', [], 'keep', ok0, 'discarded', false(size(ok0)), ...
                  'outOfBowl', false(size(ok0)), 'optErr', NaN, 'rms', NaN, ...
                  'why', '', 'suspect', false);
    if sum(ok0) < 4; info.why = 'fewer than 4 finite points'; return; end

    % ---- 1. the bowl -----------------------------------------------------
    keep = bowlMask(vals, mv, ok0, opts.bowlFrac, opts.minFitPoints);
    info.outOfBowl = ok0 & ~keep;
    if opts.verbose && any(info.outOfBowl)
        fprintf('    fitting the bowl only: %d of %d points, %s deg excluded as wings\n', ...
                sum(keep), sum(ok0), strtrim(sprintf('%.2f ', vals(info.outOfBowl))));
    end

    % ---- 2. rejection by LEAVE-ONE-OUT, worst point first ----------------
    %
    % Not "fit everything, then reject large residuals". That fails on exactly the case it is
    % meant for: a single wild point is a high-leverage point, so the fit chases it, the
    % residuals are large EVERYWHERE, the robust scale inflates to match, and nothing exceeds
    % the threshold. On a 9-point scan with one rogue at 660 among neighbours near 50, the
    % naive version rejected nothing and moved the vertex 1.7 degrees.
    %
    % So each candidate is judged against a fit that DOES NOT INCLUDE IT -- the deleted
    % residual, scaled by the scatter of the other points about that same fit. A rogue then
    % has nothing propping it up. One point is dropped per pass, the worst first, because
    % dropping several at once on a scale estimate that is still contaminated is how good
    % data gets deleted.
    b = []; resid = [];
    for iter = 1:4
        idx = find(keep);
        if numel(idx) <= opts.minFitPoints; break; end
        tWorst = 0; kWorst = 0;
        for kk = idx
            trial = keep; trial(kk) = false;
            [bt, rt, st] = fitPara(vals, mv, trial);
            if isempty(bt); continue; end
            st = max(st, 1e-4*max(range1(mv(keep)), eps));   % see the floor note below
            t = abs(mv(kk) - (bt(1) + bt(2)*vals(kk) + bt(3)*vals(kk)^2)) / st;
            if t > tWorst; tWorst = t; kWorst = kk; end
        end
        if kWorst == 0 || tWorst <= opts.outlierSigma; break; end
        if numel(idx) - 1 < opts.minFitPoints
            % Something still looks like an outlier but removing it would leave too few
            % points to fit. Say so and mark the result suspect rather than quietly
            % returning a vertex fitted through contaminated data.
            info.suspect = true;
            info.why = sprintf(['%d point(s) still look like outliers but dropping more ' ...
                                'would leave fewer than %d'], 1, opts.minFitPoints);
            if opts.verbose
                fprintf(['    STOPPING: the point at %.2f deg is %.1f sigma out but dropping ' ...
                         'it would leave fewer than %d points. The vertex below is fitted ' ...
                         'through contaminated data -- treat it as suspect.\n'], ...
                        vals(kWorst), tWorst, opts.minFitPoints);
            end
            break
        end
        if opts.verbose
            fprintf(['    discarding the point at %.2f deg: metric %.4g, %.1f sigma from a ' ...
                     'parabola fitted WITHOUT it -- that is a failed fit, not a measurement\n'], ...
                    vals(kWorst), mv(kWorst), tWorst);
        end
        info.discarded(kWorst) = true;
        keep(kWorst) = false;
    end
    [b, resid, ~] = fitPara(vals, mv, keep);
    % A scan that needed a third of its points thrown away is not a scan whose vertex you
    % should quote to three decimals, however well the survivors fit.
    if sum(info.discarded) >= max(2, ceil(0.3*sum(ok0)))
        info.suspect = true;
        if isempty(info.why)
            info.why = sprintf('%d of %d points rejected', sum(info.discarded), sum(ok0));
        end
        if opts.verbose
            fprintf(['    %d of %d points rejected. That is a lot: check the ESR fits at those ' ...
                     'angles before trusting this vertex.\n'], sum(info.discarded), sum(ok0));
        end
    end
    info.keep = keep;
    info.coeff = b;
    if isempty(b); info.why = 'parabola fit failed'; optVal = bestPoint(vals, mv, ok0); return; end
    info.rms = sqrt(mean(resid(keep).^2));

    % ---- 3. the vertex, and its error ------------------------------------
    x = vals(keep)'; y = mv(keep)';
    if b(3) <= 0
        info.why = 'curvature is not positive -- the range does not contain a minimum';
        optVal = bestPoint(vals, mv, ok0); return
    end
    xv = -b(2)/(2*b(3));
    if ~(xv > min(x) && xv < max(x))
        info.why = 'vertex falls outside the fitted points';
        optVal = bestPoint(vals, mv, ok0); return
    end

    % sigma(vertex) by propagation. xv = -b2/(2 b3), so
    %   d xv/d b2 = -1/(2 b3),   d xv/d b3 = b2/(2 b3^2)
    X = [ones(numel(x),1), x(:), x(:).^2];
    dof = max(numel(x) - 3, 1);
    s2 = sum(resid(keep).^2)/dof;
    R = X'*X;
    if rcond(R) > 1e-12
        Cb = s2 * inv(R);                                          %#ok<MINV>
        J = [0, -1/(2*b(3)), b(2)/(2*b(3)^2)];
        v = J*Cb*J';
        if v > 0; info.optErr = sqrt(v); end
    end

    optVal = xv; fromParabola = true;
    if opts.verbose
        span = max(vals(ok0)) - min(vals(ok0));
        if isfinite(info.optErr)
            fprintf('    vertex %.3f +/- %.3f deg from %d points (rms %.4g)\n', ...
                    xv, info.optErr, numel(x), info.rms);
            if info.optErr > 0.25*span/2
                fprintf(['    THE SCAN DOES NOT DETERMINE THIS ANGLE: the vertex error is %.0f%% ' ...
                         'of the half-span.\n      Wandering between rounds is then noise, not ' ...
                         'drift -- more points, more averages, or a larger |B| (the metric ' ...
                         'grows with it) rather than more rounds.\n'], ...
                        100*info.optErr/(span/2));
            end
        else
            fprintf('    vertex %.3f deg from %d points (no error estimate)\n', xv, numel(x));
        end
    end
end

function [xBar, xErr, n] = combineVertices(S, axisName, fallback)
% Inverse-variance mean of one axis's vertices, skipping the ones the fit flagged.
%
% A vertex is skipped when its fit was marked suspect, when it came from the best measured
% point rather than a vertex (no error bar, and it is quantised to the scan grid), or when it
% has no error estimate. If nothing survives, the caller's last value is returned unchanged
% -- averaging is an improvement, not a requirement.
    xBar = fallback; xErr = NaN; n = 0;
    if isempty(S); return; end
    sel = strcmp({S.axis}, axisName);
    if ~any(sel); return; end
    T = S(sel);
    x = []; w = [];
    for k = 1:numel(T)
        if ~T(k).fromParabola; continue; end
        if isfield(T(k), 'suspect') && T(k).suspect; continue; end
        e = NaN; if isfield(T(k), 'optErr'); e = T(k).optErr; end
        if ~isfinite(e) || e <= 0; continue; end
        x(end+1) = T(k).opt;  w(end+1) = 1/e^2;                        %#ok<AGROW>
    end
    if isempty(x); return; end
    n = numel(x);
    xBar = sum(x.*w)/sum(w);
    xErr = sqrt(1/sum(w));
    % If the vertices scatter by more than their error bars claim, inflate to match. Quoting
    % 1/sqrt(sum w) when the points disagree is the standard way to under-report an error.
    if n >= 2
        chi2r = sum(w.*(x - xBar).^2)/(n - 1);
        if chi2r > 1; xErr = xErr*sqrt(chi2r); end
    end
end

function keep = bowlMask(vals, mv, ok, frac, nMin)
% The contiguous run of points around the minimum whose metric is in the lowest frac of the
% range. Contiguity matters: a far-side point that happens to dip low must not join the
% bowl and stretch the fitted range back out to the wings.
    keep = ok;
    y = mv; y(~ok) = Inf;
    [~, iMin] = min(y);
    for f = [frac, 0.8, 0.9, 1.0]                 % widen if the bowl is too small
        lo = min(y); hi = max(mv(ok));
        if ~(hi > lo); keep = ok; return; end
        inBowl = ok & (mv <= lo + f*(hi - lo));
        m = false(size(ok));
        m(iMin) = true;
        for k = iMin-1:-1:1
            if inBowl(k); m(k) = true; else; break; end
        end
        for k = iMin+1:numel(ok)
            if inBowl(k); m(k) = true; else; break; end
        end
        if sum(m) >= nMin; keep = m; return; end
    end
    keep = ok;
end

function [b, resid, sRob] = fitPara(vals, mv, keep)
% Ordinary least squares parabola on the kept points; residuals returned for EVERY point so
% the caller can test the excluded ones too, and a robust scale from the kept residuals.
    b = []; resid = nan(size(mv)); sRob = 0;
    x = vals(keep)'; y = mv(keep)';
    if numel(x) < 3; return; end
    X = [ones(numel(x),1), x(:), x(:).^2];
    if rcond(X'*X) < 1e-14; return; end
    b = X \ y(:);
    resid = mv - (b(1) + b(2)*vals + b(3)*vals.^2);
    r = resid(keep);
    sRob = 1.4826*median(abs(r - median(r)));
    if ~(sRob > 0); sRob = std(r); end            % all residuals identical: fall back
    % A FLOOR ON THE SCALE, relative to the metric's own range. Without it a fit that happens
    % to be near-perfect gives sRob ~ 1e-9 and then EVERY point is thousands of sigma from
    % the curve -- on a noiseless synthetic parabola the first version declared all seven
    % points outliers. Anything within 1e-4 of the metric range is agreement, not deviation.
    sRob = max(sRob, 1e-4*max(range1(mv(keep)), eps));
end

function r = range1(v)
    v = v(isfinite(v));
    if isempty(v); r = 0; else; r = max(v) - min(v); end
end

function v = bestPoint(vals, mv, ok)
    y = mv; y(~ok) = Inf;
    [~, i] = min(y); v = vals(i);
end

function ok = setAngle(h, axis, value)
    ok = true;
    try
        if strcmp(axis, 'Btheta'); h.setBtheta(value); else; h.setBphi(value); end
    catch
        ok = false;
    end
end

function restoreAll(h, theta0, phi0)
% Only the FIELD needs restoring. The experiment object is built fresh per measurement
% and discarded, so there is no experiment state to put back.
    try; h.setBtheta(theta0); catch; end
    try; h.setBphi(phi0);     catch; end
end

function u = metricUnits(metric)
    switch metric
        case 'tilt';       u = 'SQUARED fitted tilt, deg^2 -- exactly parabolic in the scan angle';
        case 'innerspan';  u = 'SQUARED excess span of the off-axis group, MHz^2';
        case 'maxsplit';   u = 'negated largest splitting, MHz';
        case 'innerwidth'; u = 'FWHM of the widest inner line, MHz';
        otherwise;         u = 'spread of the three smaller projections, G';
    end
end

function plotScans(S, metric, tolDeg)
% One panel per scan, plus a convergence panel.
%
% The old version drew points joined by a line, a bare red vertical at the optimum, and
% "ylabel(metric)" -- which named the quantity but not its units, while metricUnits() had
% the description sitting unused two functions away. Three things were missing and all three
% are the reason you would look at this figure at all:
%
%   THE PARABOLA. The optimum is a fitted vertex, not the best measured point, and it lands
%   BETWEEN scan points. Without the curve drawn, the red line has no visible justification
%   and there is no way to see whether the fit was reasonable or whether one bad point
%   dragged it. The coefficients are now carried in S so this draws the fit that was
%   actually used.
%
%   THE DISCARDS. parabolaOptimum drops points more than 6 MADs from the median -- a failed
%   line fit is not a measurement of the metric. It announces them in text, but on the plot
%   they were indistinguishable from kept points, so a scan with a discard looked like a
%   scan whose fit ignored an obvious outlier.
%
%   THE CONVERGENCE. The whole procedure is alternating 1D scans that are supposed to
%   settle. Whether they did is the question the figure should answer first, and it was not
%   plotted anywhere.

    if isempty(S); return; end

    COL.pt   = [0.15 0.42 0.69];      % measured points
    COL.fit  = [0.85 0.33 0.10];      % the fitted parabola
    COL.opt  = [0.00 0.45 0.35];      % the chosen optimum
    COL.bad  = [0.65 0.65 0.65];      % discarded points
    COL.optBand = [0.86 0.94 0.91];   % +/- 1 sigma on the vertex
    COL.th   = [0.15 0.42 0.69];
    COL.ph   = [0.85 0.33 0.10];

    nS = numel(S);
    nCols = min(nS, 3);
    nRows = ceil(nS/nCols) + 1;       % last row: convergence

    figure('Name', 'calibFieldAlign', 'Color', 'w');

    for k = 1:nS
        ax = subplot(nRows, nCols, k); hold(ax, 'on'); tidyAlign(ax);
        v = S(k).values; y = S(k).metric;

        keep = trueOr(S(k), 'keep', isfinite(v) & isfinite(y));
        bad  = trueOr(S(k), 'discarded', false(size(v)));
        wing = trueOr(S(k), 'outOfBowl', false(size(v)));
        keep = keep & ~bad;
        % Wings are GOOD measurements outside the quadratic region, so they are drawn
        % hollow rather than crossed out. Conflating them with rejected points would say the
        % scan had five failures when it had none.

        % the fitted parabola, over the range it was fitted on
        hF = [];
        cf = [];
        if isfield(S(k), 'fitCoeff'); cf = S(k).fitCoeff; end
        if numel(cf) == 3 && any(keep)
            xf = linspace(min(v(keep)), max(v(keep)), 200);
            yf = cf(1) + cf(2)*xf + cf(3)*xf.^2;
            hF = plot(ax, xf, yf, '-', 'Color', COL.fit, 'LineWidth', 1.4);
        end

        hP = plot(ax, v(keep), y(keep), 'o', 'Color', COL.pt, 'MarkerFaceColor', COL.pt, ...
                  'MarkerSize', 5);
        hW = [];
        if any(wing & ~bad)
            hW = plot(ax, v(wing & ~bad), y(wing & ~bad), 'o', 'Color', COL.pt, ...
                      'MarkerFaceColor', 'none', 'MarkerSize', 5);
        end

        % SCALE TO THE KEPT POINTS, NOT THE DISCARDS. A discarded point is a failed line fit
        % and can be two orders of magnitude off -- one read 4527 against neighbours at 40 to
        % 180. Letting it set the axis compresses the real curve to a flat line and makes the
        % one panel that had a problem the one panel you cannot read. Off-scale discards are
        % marked at the axis edge with their value instead.
        yk = y(keep | (wing & ~bad));     % wings are real data and belong on the axis
        yl = [min(yk) max(yk)];
        if ~(diff(yl) > 0); yl = yl(1) + [-1 1]*max(abs(yl(1)), 1)*0.1; end
        yl = yl + 0.12*diff(yl)*[-1 1];
        ylim(ax, yl);

        hB = [];
        if any(bad)
            inRange = bad & y >= yl(1) & y <= yl(2);
            if any(inRange)
                hB = plot(ax, v(inRange), y(inRange), 'x', 'Color', COL.bad, ...
                          'MarkerSize', 9, 'LineWidth', 1.4);
            end
            for j = find(bad & ~(y >= yl(1) & y <= yl(2)))
                hh = plot(ax, v(j), yl(2), '^', 'Color', COL.bad, 'MarkerSize', 7, ...
                          'MarkerFaceColor', COL.bad);
                if isempty(hB); hB = hh; else; set(hh, 'HandleVisibility', 'off'); end
                text(ax, v(j), yl(2), sprintf(' %.3g ', y(j)), 'Color', COL.bad, ...
                     'FontSize', 7, 'VerticalAlignment', 'bottom', ...
                     'HorizontalAlignment', 'center');
            end
        end
        oe = NaN;
        if isfield(S(k), 'optErr'); oe = S(k).optErr; end
        if isfinite(oe) && oe > 0
            patch(ax, S(k).opt + oe*[-1 1 1 -1], [yl(1) yl(1) yl(2) yl(2)], COL.optBand, ...
                  'EdgeColor', 'none', 'HandleVisibility', 'off');
            set(ax, 'Children', flipud(get(ax, 'Children')));
        end
        hO = plot(ax, S(k).opt*[1 1], yl, '-', 'Color', COL.opt, 'LineWidth', 1.6);
        lbl = sprintf(' %.3f ', S(k).opt);
        if isfinite(oe); lbl = sprintf(' %.3f \\pm %.3f ', S(k).opt, oe); end
        text(ax, S(k).opt, yl(2), lbl, 'Color', COL.opt, ...
             'FontWeight', 'bold', 'FontSize', 8, 'VerticalAlignment', 'top');
        ylim(ax, yl);

        xlabel(ax, sprintf('%s / deg', S(k).axis));
        if mod(k-1, nCols) == 0; ylabel(ax, metric); end
        src = sprintf('vertex of %d pts', sum(keep));
        if ~S(k).fromParabola
            src = 'BEST POINT';
            if isfield(S(k), 'why') && ~isempty(S(k).why); src = ['BEST POINT: ' S(k).why]; end
        end
        extra = '';
        if any(bad); extra = sprintf(', %d dropped', sum(bad)); end
        if isfield(S(k), 'suspect') && S(k).suspect; extra = [extra ', SUSPECT']; end
        title(ax, sprintf('round %d  %s\n%.3f deg  (%s%s)', S(k).round, S(k).axis, ...
                          S(k).opt, src, extra), 'FontWeight', 'normal', 'FontSize', 9);

        if k == 1
            hs = {hP, hW, hF, hO, hB};
            ls = {'fitted (the bowl)', 'wings, not fitted', 'parabola fit', ...
                  'optimum \pm 1\sigma', 'discarded (failed fit)'};
            ok = cellfun(@(h) ~isempty(h) && all(ishandle(h)), hs);
            legend(ax, [hs{ok}], ls(ok), 'Location', 'best', 'Box', 'off', 'FontSize', 7);
        end
    end

    % ---- convergence: did the alternating scans settle? ----
    ax = subplot(nRows, nCols, (nRows-1)*nCols + (1:nCols)); hold(ax, 'on'); tidyAlign(ax);
    % PLOT THE STEPS, NOT THE ANGLES. Plotting the absolute optima puts theta near 6 deg and
    % phi near 4 on one axis, so the whole y range is the 2 deg gap BETWEEN them and the
    % 0.1 deg movement that actually matters is invisible. The loop's own stopping test is
    % the step size against tolDeg, so that is what the panel should show -- on a log axis,
    % because convergence is geometric and the interesting steps are the small ones.
    isTh = strcmp({S.axis}, 'Btheta');
    step = 1:nS;
    hT = []; hH = [];
    [xT, sT] = stepsOf(step(isTh), [S(isTh).opt]);
    [xP, sP] = stepsOf(step(~isTh), [S(~isTh).opt]);
    if ~isempty(sT)
        hT = plot(ax, xT, sT, 'o-', 'Color', COL.th, 'MarkerFaceColor', COL.th, ...
                  'LineWidth', 1.3, 'MarkerSize', 5);
    end
    if ~isempty(sP)
        hH = plot(ax, xP, sP, 's-', 'Color', COL.ph, 'MarkerFaceColor', COL.ph, ...
                  'LineWidth', 1.3, 'MarkerSize', 5);
    end
    hTol = [];
    if isfinite(tolDeg) && tolDeg > 0
        hTol = plot(ax, [0.5 nS+0.5], tolDeg*[1 1], '--', 'Color', [0.35 0.35 0.35], ...
                    'LineWidth', 1.2);
    end
    set(ax, 'YScale', 'log');
    xlabel(ax, 'scan number (theta and phi alternate)');
    ylabel(ax, 'step from previous scan / deg');
    xlim(ax, [0.5 nS+0.5]);
    hs = {hT, hH, hTol};
    ls = {'Btheta step', 'Bphi step', sprintf('tolDeg = %.3g', tolDeg)};
    ok = cellfun(@(h) ~isempty(h) && all(ishandle(h)), hs);
    if any(ok)
        legend(ax, [hs{ok}], ls(ok), 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
    end

    mvT = lastMove([S(isTh).opt]);
    mvP = lastMove([S(~isTh).opt]);
    done = 'CONVERGED';
    if ~(isfinite(tolDeg) && max([mvT mvP]) < tolDeg); done = 'NOT converged -- rounds ran out'; end
    if ~isfinite(mvT) && ~isfinite(mvP); done = 'only one scan per axis -- no step to judge'; end
    title(ax, sprintf(['convergence: %s  --  last theta step %s, last phi step %s' ...
                       '     [%s = %s]'], done, degStr(mvT), degStr(mvP), ...
                      metric, metricUnits(metric)), 'FontWeight', 'bold', 'FontSize', 9);
end

function [x, d] = stepsOf(idx, v)
% Consecutive absolute changes, on a log axis, so a zero step has to become something
% plottable rather than vanishing. A step below a milli-degree is converged by any standard
% here, so flooring at 1e-3 loses nothing and keeps the marker visible.
    x = []; d = [];
    if numel(v) < 2; return; end
    d = max(abs(diff(v)), 1e-3);
    x = idx(2:end);
end

function tidyAlign(ax)
    set(ax, 'Box', 'off', 'TickDir', 'out', 'Layer', 'top', ...
            'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.12, 'FontSize', 8);
end

function v = trueOr(S, field, dflt)
% A logical mask from S if it is there and the right size, otherwise the default. Keeps this
% plotter working on scan structs saved before fitCoeff/keep/discarded existed.
    v = dflt;
    if isfield(S, field) && ~isempty(S.(field)) && numel(S.(field)) == numel(dflt)
        v = logical(S.(field));
    end
end

function d = lastMove(v)
    if numel(v) < 2; d = NaN; else; d = abs(v(end) - v(end-1)); end
end

function s = degStr(d)
    if isfinite(d); s = sprintf('%.3f deg', d); else; s = 'n/a'; end
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
