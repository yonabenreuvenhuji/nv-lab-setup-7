function out = calibB1Rabi(opts)
%CALIBB1RABI Find the B1 maximum by MEASURING B1, using the Rabi frequency as the probe.
%
%   out = calibB1Rabi
%   out = calibB1Rabi(struct('rhoCoarse', 13:1:23))
%   out = calibB1Rabi(struct('twoPass', false, 'rhoCoarse', 16:0.5:20))
%
% LEVEL 3, last step -- the handover to the microwave calibrations. Replaces calibB1Peak.
%
% At each field setting: run an ESR to find where the lowest NV line actually is, then run
% a Rabi AT that frequency and read off the Rabi frequency. The Rabi frequency is
% proportional to B1, so the scan maps out the antenna response directly. Two passes, wide
% then narrow, and the answer is the field setting that puts the line at the response peak.
%
% ===================================================================
% WHY THIS REPLACES THE LINEWIDTH METHOD
% ===================================================================
%
% calibB1Peak inferred B1 from a CW ESR dip. That needs a saturation model to get from a
% width to a field:
%
%       depth = Cmax * s/(1+s)          FWHM = FWHM0 * sqrt(1 + s)
%
% so B1 only enters through s, the two observables saturate at different rates, and the
% intrinsic width FWHM0 has to be known or fitted alongside. Worse, the measured width
% also carries inhomogeneous broadening, strain and hyperfine structure, none of which
% has anything to do with how hard the NVs are being driven. The best that method can do
% is a monotone proxy.
%
% The Rabi frequency is not a proxy. It IS the drive:
%
%       Omega = gamma * B1 * <matrix element>
%
% and on this scan the matrix element does not move. The state mixing that would change
% it goes as (gamma*Bperp/D)^2, which at 20 G is (56/2870)^2 = 4e-4 even for a fully
% off-axis orientation and far less for the aligned one this measures. So Omega is B1 to
% about one part in 5000, with no model in between, and it comes with a real error bar
% from the fit rather than from a saturation curve nobody measured.
%
% The cost is time: a Rabi per field point on top of the ESR, and the Rabi cannot reuse
% its uploaded waveforms because the carrier frequency changes at every point and the
% carrier is baked into them. Budget about two minutes per field setting.
%
% ===================================================================
% WHAT IS FITTED, AND WHY NOT A PARABOLA IN FIELD
% ===================================================================
%
% The ESR gives the line frequency at every point, so the natural x axis is FREQUENCY, not
% the coil setting. The resonator responds to frequency; the field is only the knob that
% slides the NV line along it. And because a resonator passes POWER as a Lorentzian while
% B1 goes as the square root of power,
%
%       Omega^2(f) = A / (1 + 4 (f - f0)^2 / W^2)  +  c
%
% is the model, fitted by resonatorLorentzianFit. That gives the centre f0, the width W
% and therefore Q, and the on-resonance Rabi frequency -- four physical numbers from the
% same points a parabola would have reduced to one. The coil calibration then only has to
% convert f0 back into a field setting, and it is fitted from this scan's own (Brho, line)
% pairs rather than trusted from an earlier one.
%
% A weighted parabola is computed too, as an independent check on f0. It is the honest
% fallback when the scan turns out too short to determine a width: on synthetic data with
% a known centre, a 12 MHz span was too short for the Lorentzian (its width pinned at a
% bound and the fit was rejected) while the parabola still found the centre to 0.2 MHz.
%
% ===================================================================
% TWO PASSES
% ===================================================================
%
% COARSE, wide field range, wide ESR sweep. The wide sweep is what lets calibFieldSurvey
% fit the field vector, which is what identifies WHICH orientation the lowest line belongs
% to. That identification is not optional: on the first live run of calibB1Peak, three of
% nine field points had the aligned pair outside the sweep, the "lowest line" was a
% different orientation's line, and the accepted junk dragged the reported peak by 18 MHz
% and turned the coil calibration into 0.13*Brho + 9.7. Those gates are kept here.
%
% FINE, a few gauss around the coarse peak. By then the (Brho -> line frequency) relation
% is measured, so the ESR sweep can be a narrow window that follows the prediction --
% about a quarter of the points, and the same protection comes from CONTINUITY with the
% coarse relation instead of from the field fit. A line that jumped to a different
% orientation would sit far off the predicted position and be rejected on that. Set
% fineNarrowEsr = false to keep the wide sweep on the fine pass too, at roughly double the
% time.
%
% The response curve is then fitted to BOTH passes together, not to the fine pass alone.
% That matters: a 4 G fine scan spans about 11 MHz, which is far too short to determine a
% width on its own, so fitted by itself it would give a good centre and a meaningless Q.
% Combined, the coarse points supply the wings and the fine points supply resolution where
% the curvature is, and one fit gets both.
%
% ===================================================================
% THE ONE THING THAT MUST NOT CHANGE DURING THE SCAN
% ===================================================================
%
% The drive amplitude. Omega is proportional to B1, so any change in the level being sent
% is indistinguishable from a change in what the antenna delivers -- it would appear as a
% real feature of the response curve. rabiAmplitude is applied identically at every point
% and nothing is allowed to re-derive it mid-scan. If reduced chi2 on the Lorentzian comes
% back large, this is the first thing to suspect.
%
% ===================================================================
% REPLAY -- finishing a scan that already ran
% ===================================================================
%
% The scan is the expensive part: half an hour of ESR and Rabi runs. The fit afterwards costs
% milliseconds. So anything that kills this function AFTER the scan -- a bug in the fit, a bad
% option, a plotting failure -- destroys thirty minutes of measurement to save nothing.
%
%   out = calibB1Rabi(struct('replayPoints', tab, 'plotIt', true))
%
% skips the hardware entirely and goes straight to the fit. replayPoints takes either a
% previous run's out.points (a struct array -- exact, nothing assumed) or a NUMERIC MATRIX
% laid out exactly like the scan table this function prints:
%
%   Brho/G   line/MHz   FWHM   fRabi/MHz   pi/ns   contr%   R^2   R2max
%
% which is the recovery path when the run died and the Command Window is all that is left.
% Select the table, paste it between brackets, pass it in. 8 columns is the full table;
% 4 is [rho fLine fRabi fRabiErr]; 3 is [rho fLine fRabi]. Rejected rows -- the ones whose
% fRabi printed as NaN -- can be pasted in as they are and stay rejected, so the usable count
% matches what the run reported.
%
% The table has no error bars in it, so a matrix replay assumes them, and then MEASURES
% whether that mattered: it refits across an eightfold range of assumed error and prints how
% far the answer moved. On a well-sampled curve it moves by nothing, because equal relative
% weights on points of equal quality is the right fit anyway. When it does move, that spread
% is the honest uncertainty and the printout says to quote it instead of the fit's.
%
% A replay needs no instruments at all: with no readable Brho it says so, turns leaveAtPeak
% off, and tells you the field to set by hand. verify is forced off -- it would put a Rabi
% back on the hardware, which is the one thing a replay exists to avoid.
%
% ===================================================================
% OPTIONS
% ===================================================================
%
%   replayPoints   [] = measure (default). A struct array or a numeric table skips the
%                  hardware and fits what is already measured -- see REPLAY above.
%   replayFErrFrac assumed fRabiErr/fRabi for a matrix replay (default 0.01)
%   replayFLineErr assumed fLineErr, MHz, for a matrix replay (default 0.2)
%
%   rhoCoarse      coarse field list, G (default: current Brho -7 to +7 in 11 steps --
%                  +/- 7 G is the narrowest scan that determines Q every time, see the
%                  table where the default is set)
%   twoPass        run the fine pass (default true)
%   rhoFine        fine field list, G (default: rhoFineSpan wide, nFine points, centred
%                  on the coarse peak)
%   rhoFineSpan, nFine        (defaults 4 G, 9)
%
%   RABI -- these default to the walk block that works on this setup, and they are the
%   standard for the calibrations now. Changing them changes what B1 means.
%   rabiAmplitude  dBm per channel before amplitudeDiff (default -16)
%   repeats, averages         (defaults 300, 2)
%   tau            us (default 0.02:0.01:1.0, 99 points)
%   useSequenceWalk, passesPerRun, walkWarmupTriggers, resyncEveryPass
%                  (defaults true, [] = auto, 1, true)
%   newOrderEveryAverage  Rabi runs (walk mode only): a new random tau order before every
%                  average (true/false). Default [] = AWGExpRabi's own setting (true).
%                  Costs one sequence upload + arm per average. See
%                  AWGExpRabi.newOrderEveryAverage.
%   rabiDetectionMHz          detection window, us (default 0.3, calibDetectionWindow's)
%   nChannels      (default 2)
%   nFreq          frequencies in the Rabi model (default 1; set 2 if the residuals are
%                  autocorrelated, which usually means a second orientation is being
%                  driven off resonance)
%   minContrast    reject a Rabi below this contrast (default 0.005 = 0.5%). Catches the
%                  case where the microwave never reached the sample, which otherwise
%                  returns a ratio near 1 at every tau with no error at all.
%   maxFErrFrac    reject when fErr/f exceeds this (default 0.05). THIS IS THE REAL GATE.
%                  It asks the only question that matters -- is the Rabi frequency
%                  determined -- and it asks it in the units the answer is used in.
%   minR2          absolute floor on R^2. DEFAULT 0, i.e. OFF, and that is deliberate.
%
%                  R^2 is the fraction of the variance the model explains, so it is bounded
%                  above by how much of the variance is signal rather than noise. Write A
%                  for the fitted half-contrast and sigma for the error on a point; the
%                  oscillation contributes about A^2/2 to the variance, so
%
%                        R^2  <=  (A^2/2) / (A^2/2 + sigma^2)
%
%                  no matter how perfect the fit. At the 2.5% contrast this setup gives at
%                  300 repeats, that ceiling is around 0.8 -- so the 0.90 default this
%                  function shipped with rejected sound fits for having a low signal-to-
%                  noise ratio, which is not a fitting failure and is not fixed by fitting
%                  differently. A LOW R^2 HERE MEANS FEWER PHOTONS, NOT A WORSE FIT.
%
%                  Every point now prints its R^2 next to that ceiling (R2max), so you can
%                  see whether a fit is as good as the noise allows. Set minR2 to something
%                  like 0.3 if you want a floor that catches genuine garbage; anything near
%                  0.9 will throw away good data.
%   maxResidAutocorr   flag (not reject) when residuals are autocorrelated (default 0.5)
%
%   ESR
%   esrFreqs       THE WIDE SWEEP, as an explicit MHz list, e.g. 2800:2:2942. Default []
%                  means work it out from the field range: D0 +/- gamma*coilRatio*max(Brho)
%                  plus esrMargin, in esrStepCoarse steps. Whichever way it is set, the
%                  sweep is checked against every field in the scan BEFORE the scan starts,
%                  and if it is too narrow at the top the printout names the largest Brho it
%                  supports and prints the esrFreqs line that would fix it.
%
%                  This is the sweep used on the coarse pass, and on the fine pass too when
%                  fineNarrowEsr is false. The fine pass's narrow windows cannot be given as
%                  a fixed list because they have to FOLLOW the line as the field changes;
%                  esrNarrowHalf and esrStepFine set their width and resolution instead.
%   esrStepCoarse, esrStepFine   MHz per point (defaults 2, 0.5). esrStepCoarse is ignored
%                  when esrFreqs is given, since the list carries its own spacing.
%   esrMargin      extra MHz beyond the outermost expected line, so the fitter has baseline
%                  on both sides (default 8). Ignored when esrFreqs is given.
%   esrNarrowHalfMax   HARD CAP on the fine-pass half-width, MHz ([] = none). Use this when
%                  you want to force a narrow window; esrNarrowHalf cannot do it because it
%                  is a floor, so setting it small changes nothing.
%   esrNeighbourGuard  the fine-pass window may reach this fraction of the way to the
%                  neighbouring orientation (default 0.45). The off-axis group sits
%                  (2/3)*|D0 - fLine| away, so this is what stops the window swallowing it.
%                  See narrowWindow.
%   esrNarrowHalf  FLOOR on the half-width of a fine-pass window, MHz (default 15). The
%                  window actually used is max(esrNarrowHalf, narrowHalfFwhmFactor * the
%                  median linewidth measured on the coarse pass), so it adapts to how broad
%                  the lines are at your drive level instead of being a constant that
%                  happens to be too small.
%   narrowHalfFwhmFactor   how many linewidths of half-width a fine window gets (default
%                  2.5). Below about 1 the line is wider than the window and gets rejected
%                  as baseline no matter how good the fit is; 2.5 also leaves the fitter's
%                  baseline some clean wings to sit on.
%   esrAmplitude   dBm (default -35, the level the CW calibrations ran at)
%   esrRepeats, esrAverages    (defaults 50, 3)
%   esrDetectionMHz  CW window for the ESR step, us. Default [] = resolved through
%                  esrWindowPick by calibFieldSurvey. Leave it empty.
%   fineNarrowEsr  narrow window on the fine pass (default true)
%   esrNewOrderEveryAverage  the ESR steps: a new random frequency order before every
%                  average (true/false), passed down to calibFieldSurvey as its
%                  newOrderEveryAverage. Default [] = AWGExpESR's own setting (true).
%                  Separate from newOrderEveryAverage, which is the Rabi runs'.
%   baselineOrder, freshExperiment, minFwhm, edgeGuard      as calibB1Peak
%
%   D0, coilRatio  only used to SIZE the first sweep, never to compute a result
%   settle         seconds after setting the field (default 3)
%   baselineInFit  let the Lorentzian have a floor (default true)
%   leaveAtPeak    set the field to the answer at the end (default true)
%   verify         one more ESR + Rabi there to confirm (default true)
%   plotIt, controlAuto, isTracking, jsonLocation
%
% Restores Brho on any error. Btheta and Bphi are never touched -- align first with
% calibFieldAlign; only the magnitude moves here.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct( ...
        'rhoCoarse', [], 'rhoFine', [], 'rhoFineSpan', 4, 'nFine', 9, 'twoPass', true, ...
        'rabiAmplitude', -16, 'repeats', 300, 'averages', 2, 'tau', [], ...
        'useSequenceWalk', true, 'passesPerRun', [], 'walkWarmupTriggers', 1, ...
        'resyncEveryPass', true, 'rabiDetectionMHz', 0.3, 'nChannels', 2, 'nFreq', 1, ...
        'minContrast', 0.005, 'minR2', 0, 'maxFErrFrac', 0.05, 'maxResidAutocorr', 0.5, ...
        'esrFreqs', [], 'esrAmplitude', -35, 'esrRepeats', 50, 'esrAverages', 3, ...
        'esrDetectionMHz', [], ...
        'esrStepCoarse', 2, 'esrStepFine', 0.5, 'esrMargin', 8, 'esrNarrowHalf', 15, ...
        'esrNarrowHalfMax', [], 'esrNeighbourGuard', 0.45, ...
        'fineNarrowEsr', true, 'narrowHalfFwhmFactor', 2.5, ...
        'baselineOrder', 3, 'freshExperiment', true, ...
        'esrNewOrderEveryAverage', [], 'newOrderEveryAverage', [], ...
        'minFwhm', 3, 'edgeGuard', 4, ...
        'D0', 2870.5, 'coilRatio', 0.97, 'settle', 3, 'baselineInFit', true, ...
        'leaveAtPeak', true, 'verify', true, 'plotIt', true, 'controlAuto', true, ...
        'isTracking', 0, ...
        'replayPoints', [], 'replayFErrFrac', 0.01, 'replayFLineErr', 0.2, ...
        'jsonLocation', ''), ...
        'calibB1Rabi');

    % ---- replay: finish a scan whose data is already measured -------------
    % See the REPLAY section of the header. Nothing is measured, so the point count, the
    % window checks and the runtime estimate below are all meaningless -- and verify would
    % put a Rabi back on the hardware, which is the one thing a replay is meant to avoid.
    replaying = ~isempty(opts.replayPoints);
    if replaying
        opts.verify = false;
    end

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

    GAMMA = 2.8025;
    if isempty(opts.tau); opts.tau = round(0.02 : 0.01 : 1.0, 3); end
    out = struct('ok', false, 'opts', opts, 'gamma', GAMMA);

    h = getObjByName('Helmholtz');
    if isempty(h)
        try; h = Helmholtz; catch; end
    end
    if isempty(h) && ~replaying; error('calibB1Rabi: no Helmholtz object found.'); end

    rho0 = NaN; theta0 = NaN; phi0 = NaN;
    if ~isempty(h); [rho0, theta0, phi0] = readSpherical(h); end
    fprintf('\n=== calibB1Rabi ===\n');
    fprintf('current field: Brho = %.3f G, Btheta = %.3f deg, Bphi = %.3f deg\n', ...
            rho0, theta0, phi0);
    if ~isfinite(rho0) || rho0 <= 0
        % NOT fatal when replaying. A replay measures nothing, so it does not need the coils
        % at all -- it should be runnable from a laptop with no instruments in the room. The
        % only thing rho0 is for is restoring the field afterwards, and there is nothing to
        % restore. leaveAtPeak has to go off with it, since it would try to set the coils.
        if replaying
            fprintf(['No readable Brho, which is fine for a replay: nothing is measured and\n' ...
                     'there is no field to restore. leaveAtPeak is off, so SET THE FIELD BY\n' ...
                     'HAND to the Brho this prints at the end.\n']);
            opts.leaveAtPeak = false;
        else
            error(['calibB1Rabi: could not read the current Brho (%g). Pass rhoCoarse ' ...
                   'explicitly, and check the field is ON -- h.control(''on'') in LOWER case ' ...
                   'passes validation and then turns the supplies OFF. Only h.control(''ON'') ' ...
                   'energises them.'], rho0);
        end
    end
    out.rho0 = rho0; out.theta0 = theta0; out.phi0 = phi0;

    % THE COARSE HALF-WIDTH IS SET BY WHAT THE FIT NEEDS, NOT BY TASTE.
    %
    % A Lorentzian's width is only determined once the scan covers a decent fraction of it.
    % Below that its shoulder is indistinguishable from a parabola: the CENTRE survives but
    % (f0, W, A) trade off against each other and Q comes out biased. Measured on synthetic
    % scans of a known Q = 94 resonator, 60 noise realisations at each width:
    %
    %       +/- 5 G    27 MHz span    fit usable 77 % of the time,   Q came out 101  (7 % high)
    %       +/- 6 G    33 MHz         usable 88 %,                   Q 94.5
    %       +/- 7 G    38 MHz         usable 100 %,                  Q 93.9
    %       +/- 10 G   55 MHz         usable 100 %,                  Q 93.9
    %
    % So +/- 7 G is the smallest scan that pins the width every time, and that is the
    % default. A NARROWER scan still finds the peak -- the parabola fallback gets the centre
    % to better than 1 MHz down to an 8 MHz span -- so pass rhoCoarse yourself if you only
    % want the field setting and not Q. It just should be a decision rather than an accident.
    rhosC = opts.rhoCoarse;
    if isempty(rhosC); rhosC = linspace(max(2, rho0 - 7), rho0 + 7, 11); end
    rhosC = sort(unique(round(rhosC(:)'*1000)/1000));

    nRabiPoints = numel(rhosC) + opts.twoPass*opts.nFine + double(opts.verify);
    if ~replaying
    fprintf(['%d field points in the coarse pass%s, so %d Rabi runs in all.\n' ...
             'Each one reloads the AWG because the carrier frequency changes, so budget\n' ...
             'about two minutes per point -- roughly %.0f minutes total. Ctrl-C is safe:\n' ...
             'the field is restored on any error.\n'], ...
            numel(rhosC), BooleanHelper.ifTrueElse(opts.twoPass, ...
                sprintf(' and %d in the fine pass', opts.nFine), ''), ...
            nRabiPoints, 2*nRabiPoints);
    fprintf('Rabi standard: %d taus, repeats %d, averages %d, amplitude %g dBm, walk %d.\n', ...
            numel(opts.tau), opts.repeats, opts.averages, opts.rabiAmplitude, opts.useSequenceWalk);
    end

    if opts.controlAuto && ~isempty(h) && ~replaying
        try; h.control('auto'); catch; end
    end

    try
        if replaying
        % ================= REPLAY, NO HARDWARE =================
        [P, PC] = pointsFromReplay(opts.replayPoints, opts);
        out.coarse = PC;
        rel = lineVsField(PC, opts, 'replayed points', GAMMA);
        freqsC = [];
        else
        % ================= COARSE PASS =================
        freqsC = esrWindowFor(max(rhosC), opts, GAMMA);
        checkWindow(freqsC, rhosC, opts, GAMMA);

        printHeader('COARSE');
        PC = [];
        for i = 1:numel(rhosC)
            p = measureAt(h, rhosC(i), freqsC, [], 'coarse', opts);
            printRow(p);
            PC = [PC p];                                                    %#ok<AGROW>
        end
        out.coarse = PC;

        rel = lineVsField(PC, opts, 'coarse pass', GAMMA);

        P = PC;

        % ================= FINE PASS =================
        if opts.twoPass
            [rhoMid, howMid] = coarsePeakRho(PC, rel);
            if ~isfinite(rhoMid)
                fprintf(['\nSKIPPING THE FINE PASS: the coarse pass gave nothing to centre it on\n' ...
                         '(%s). Fix the coarse pass first -- a fine scan around a peak that was\n' ...
                         'never found is nine more field points spent on nothing.\n'], howMid);
            else
                rhosF = opts.rhoFine;
                if isempty(rhosF)
                    rhosF = linspace(rhoMid - opts.rhoFineSpan/2, ...
                                     rhoMid + opts.rhoFineSpan/2, opts.nFine);
                end
                rhosF = sort(unique(round(rhosF(:)'*1000)/1000));
                fprintf('\nFine pass: %.2f to %.2f G in %d steps, centred on %.2f G (%s).\n', ...
                        min(rhosF), max(rhosF), numel(rhosF), rhoMid, howMid);

                freqsF = freqsC; relF = [];
                if opts.fineNarrowEsr && rel.ok
                    % THE WINDOW IS SIZED FROM THE LINEWIDTH WE JUST MEASURED, not from a
                    % constant. A window only a little wider than the line is wrong twice
                    % over: the fitter's baseline has no clean wings to sit on, so it eats
                    % part of the line; and esrFitRobust rejects a component whose FWHM
                    % exceeds the span, which is the correct test but bites hard when the
                    % span is small. A fixed +/- 15 MHz did exactly this on a 12 MHz line --
                    % the fit found the dip at the predicted frequency with a sensible depth
                    % and width, and it was thrown out anyway.
                    fwhmMed = medianFwhmOf(PC);
                    relF = rel;  freqsF = [];
                    relF.narrowFwhm   = fwhmMed;
                    relF.narrowFloor  = opts.esrNarrowHalf;
                    relF.narrowFactor = opts.narrowHalfFwhmFactor;
                    relF.narrowMax    = opts.esrNarrowHalfMax;
                    relF.guard        = opts.esrNeighbourGuard;
                    % The window is sized PER FIELD POINT, because the distance to the
                    % neighbouring orientation is proportional to the field and the fine pass
                    % spans a range of fields. A single number computed here would be right at
                    % one end of the scan and wrong at the other -- which is exactly the bug:
                    % a +/- 30 MHz window is clear of the neighbour at 24 G and 8 MHz past it
                    % at 12 G. narrowWindow does the arithmetic; this just reports it at the
                    % middle of the scan so the printout means something.
                    fMid = rel.a*mean([min(rhosF) max(rhosF)]) + rel.b;
                    [loMid, hiMid, gapMid] = narrowWindow(fMid, relF, opts);
                    nPts = 1 + round((loMid + hiMid)/opts.esrStepFine);
                    tolMid = min(0.5*gapMid, max(loMid, hiMid));
                    fprintf(['Narrow ESR windows, ASYMMETRIC: -%.1f / +%.1f MHz around the\n' ...
                             'predicted line, %g MHz steps (%d points instead of %d).\n'], ...
                            loMid, hiMid, opts.esrStepFine, nPts, numel(freqsC));
                    fprintf(['Asymmetric because the two sides are limited by different things. The\n' ...
                             'aligned lower line is the OUTERMOST line in the spectrum, so below it\n' ...
                             'there is nothing and the window can be as generous as the fit wants --\n' ...
                             '%.1f x the %.1f MHz linewidth measured on the coarse pass, which leaves\n' ...
                             'the baseline clean wings to sit on. Above it sits the off-axis group at\n' ...
                             '(2/3)*gamma*|B| = %.1f MHz away, and a window that reaches it fits the\n' ...
                             'wrong dip. So that side is capped at %.0f%% of the gap.\n'], ...
                            loMid/max(fwhmMed, eps), fwhmMed, gapMid, 100*opts.esrNeighbourGuard);
                    if hiMid < 1.2*fwhmMed
                        fprintf(['NOTE the capped side is only %.1f x the linewidth. The neighbour is\n' ...
                                 'close at this field, so there is no window that is both wide enough\n' ...
                                 'for a clean baseline and clear of it. Expect some points to be\n' ...
                                 'rejected; raising the field would separate them.\n'], ...
                                hiMid/max(fwhmMed, eps));
                    end
                    fprintf(['The field fit is not available in a window this narrow, so the\n' ...
                             'orientation check becomes CONTINUITY with the coarse relation: a line\n' ...
                             'more than %.1f MHz from the prediction is rejected, which is what a\n' ...
                             'jump to a different orientation looks like.\n'], tolMid);
                elseif opts.fineNarrowEsr
                    fprintf(['Wide ESR windows on the fine pass after all: the coarse pass did not\n' ...
                             'give a usable line-vs-field relation, so there is nothing to point a\n' ...
                             'narrow window at.\n']);
                end

                printHeader('FINE');
                PF = [];
                for i = 1:numel(rhosF)
                    p = measureAt(h, rhosF(i), freqsF, relF, 'fine', opts);
                    printRow(p);
                    PF = [PF p];                                            %#ok<AGROW>
                end
                out.fine = PF;
                P = [PC PF];
            end
        end
        end

        out.points = P;
        good = logical([P.good]);
        out.good = good;
        fprintf('\n%d of %d field points usable.\n', sum(good), numel(P));
        listRejects(P);
        if sum(good) < 5
            error(['calibB1Rabi: only %d usable points -- not enough to fit a response ' ...
                   'curve. The rejection reasons above say what to fix; a bad sweep window ' ...
                   'or a Rabi with no contrast are both cheap to correct.'], sum(good));
        end

        % ================= THE FIT =================
        G = P(good);
        fL   = [G.fLine];      fLerr = [G.fLineErr];
        wR   = [G.fRabi];      wRerr = [G.fRabiErr];
        rhoG = [G.rho];

        fit = resonatorLorentzianFit(fL, wR, wRerr, ...
                  struct('baseline', opts.baselineInFit, 'verbose', true));
        out.fit = fit;

        relAll = lineVsField(G, opts, 'all usable points', GAMMA);
        out.relation = relAll;

        % Which frequency to aim at: the Lorentzian centre when it is usable, otherwise the
        % parabola vertex. Both are reported either way so the choice is visible rather
        % than buried.
        [fTarget, fTargetErr, howTarget] = pickTarget(fit);
        out.resonantFrequency = fTarget;
        out.resonantFrequencyErr = fTargetErr;
        out.targetFrom = howTarget;

        [rhoStar, rhoStarErr] = relAll.invert(fTarget, fTargetErr);
        out.rhoOptimum = rhoStar; out.rhoOptimumErr = rhoStarErr;

        fprintf('\n=== the answer ===\n');
        fprintf('  drive the NVs at        %.2f +/- %.2f MHz   (%s)\n', ...
                fTarget, fTargetErr, howTarget);
        fprintf('  set the field to Brho = %.3f +/- %.3f G\n', rhoStar, rhoStarErr);
        if fit.ok
            fprintf('  expect fRabi           %.3f +/- %.3f MHz  -> pi pulse %.1f +/- %.1f ns\n', ...
                    fit.wMax, fit.wMaxErr, fit.piMinNs, fit.piMinErrNs);
            out.fRabiAtOptimum = fit.wMax; out.piNsAtOptimum = fit.piMinNs;
        end
        [bestW, iBest] = max(wR);
        fprintf('  best point MEASURED    %.3f MHz at %.2f MHz, Brho %.3f\n', ...
                bestW, fL(iBest), rhoG(iBest));

        % A fitted optimum outside the fields actually visited is an extrapolation through
        % the coil calibration as well as through the response curve, and both extrapolate
        % badly. Say so.
        if rhoStar < min(rhoG) || rhoStar > max(rhoG)
            fprintf(['\n  *** THE OPTIMUM FIELD (%.2f G) IS OUTSIDE THE FIELDS SCANNED (%.2f to\n' ...
                     '  %.2f G). Both the response curve and the coil calibration are being\n' ...
                     '  extrapolated. Re-run with rhoCoarse extended in that direction.\n'], ...
                    rhoStar, min(rhoG), max(rhoG));
        end
        checkAgainstOldMethod(out, fTarget);
        if replaying && ~isstruct(opts.replayPoints)
            replaySensitivity(P, good, opts, fTarget, rhoStar, relAll);
        end

        if opts.plotIt; drawIt(P, good, fit, relAll, fTarget, rhoStar, opts); end

        % ================= FINISH =================
        out.verifyPoint = [];
        if opts.leaveAtPeak && isfinite(rhoStar) && rhoStar > 0
            fprintf('\nsetting Brho = %.3f G\n', rhoStar);
            relV = []; freqsV = freqsC;
            if opts.fineNarrowEsr && relAll.ok; relV = relAll; freqsV = []; end
            if replaying; relV = []; freqsV = []; end   % verify is off anyway when replaying
            if opts.verify
                fprintf('verifying with one ESR + Rabi there...\n');
                pv = measureAt(h, rhoStar, freqsV, relV, 'verify', opts);
                printHeader('VERIFY'); printRow(pv);
                out.verifyPoint = pv;
                reportVerify(pv, fit, fTarget);
            else
                h.setBrho(rhoStar); pause(opts.settle);
            end
            fprintf('\n  field LEFT at Brho = %.3f G, direction untouched.\n', rhoStar);
            out.rhoLeftAt = rhoStar;
        elseif isfinite(rho0)
            restoreRho(h, rho0);
            fprintf('\n  field restored to Brho = %.3f G (leaveAtPeak is off).\n', rho0);
            out.rhoLeftAt = rho0;
        else
            fprintf('\n  no field to restore (replay).  SET Brho = %.3f G BY HAND.\n', rhoStar);
            out.rhoLeftAt = NaN;
        end

        out.ok = isfinite(fTarget) && isfinite(rhoStar) && ...
                 (fit.usable || (fit.parabola.ok && fit.parabola.inRange));
        printHandover(out);

    catch err
        if isfinite(rho0)
            restoreRho(h, rho0);
            fprintf('\n(field restored to Brho = %.3f G after the error)\n', rho0);
        end
        rethrow(err);
    end
end

% ======================================================================
% REPLAY
% ======================================================================

function [P, PC] = pointsFromReplay(R, opts)
%POINTSFROMREPLAY Rebuild a point struct array from data that has already been measured.
%
% Two accepted forms.
%
%   A STRUCT ARRAY, i.e. a previous run's out.points. Everything is present, nothing is
%   assumed, and the answer is bit-identical to what that run would have produced.
%
%   A NUMERIC MATRIX, one row per field point, laid out as the scan table printed during the
%   run -- which is the point of accepting it. If the run died after the scan, the table in
%   the Command Window is the only surviving record, and it can be selected and pasted
%   straight into a matrix literal:
%
%       Brho/G   line/MHz   FWHM   fRabi/MHz   pi/ns   contr%   R^2   R2max
%       [ 1 ]      [ 2 ]    [ 3 ]    [ 4 ]     [ 5 ]    [ 6 ]   [ 7 ]  [ 8 ]
%
%   8 columns is the full table. 4 columns is [rho fLine fRabi fRabiErr]. 3 columns is
%   [rho fLine fRabi], and then the error bars are assumed -- see below.
%
% WHAT THE TABLE DOES NOT CONTAIN, and what that costs.
%
% The printed table has no error bars: fLineErr and fRabiErr were computed and used but never
% printed. So a replay from a pasted table has to assume them, from replayFErrFrac (relative,
% on fRabi) and replayFLineErr (absolute, MHz).
%
% This matters LESS than it looks, and it is worth being precise about why rather than just
% asserting it. The weights enter the Lorentzian fit only through their RATIOS. Assuming one
% relative error for every point makes the fit equally weighted, which is the correct fit
% when the points really are of similar quality -- and they are, because every point that
% survived rejection passed the same maxFErrFrac gate. What the assumption does corrupt is
% the ABSOLUTE scale of the parameter errors, i.e. the error bar on the peak position and on
% Q. resonatorLorentzianFit Birge-scales its errors by the reduced chi2, which recovers much
% of that, but not the part that depends on knowing the noise independently.
%
% So: the FREQUENCY and the FIELD are sound. Their error bars, and Q, are indicative. The
% sensitivity block printed below measures exactly this rather than leaving you to trust it.
%
% A rejected row -- one whose fRabi is NaN because the point failed -- can be pasted in as is
% and is marked bad, so the count of usable points still matches what the run reported.

    if isstruct(R)
        P = R(:)';
        need = {'rho', 'fLine', 'fLineErr', 'fRabi', 'fRabiErr', 'good', 'pass'};
        miss = need(~ismember(need, fieldnames(P)));
        if ~isempty(miss)
            error(['calibB1Rabi: replayPoints is a struct array but has no "%s" field. ' ...
                   'Pass a previous run''s out.points, or a numeric table -- see ' ...
                   'pointsFromReplay.'], miss{1});
        end
        PC = P(strcmp({P.pass}, 'coarse'));
        if isempty(PC); PC = P; end
        fprintf('\nREPLAY from a struct array: %d points, no hardware, nothing assumed.\n', ...
                numel(P));
        return
    end

    if ~isnumeric(R) || isempty(R) || ndims(R) > 2
        error('calibB1Rabi: replayPoints must be a numeric matrix or a struct array.');
    end
    if size(R, 2) == 1 && size(R, 1) > 1
        error(['calibB1Rabi: replayPoints has one column. Rows are field points and columns ' ...
               'are quantities -- you have probably pasted a single column of Brho values.']);
    end
    nc = size(R, 2);
    if ~ismember(nc, [3 4 8])
        error(['calibB1Rabi: replayPoints has %d columns. Accepted: 8 (the printed scan ' ...
               'table), 4 ([rho fLine fRabi fRabiErr]), or 3 ([rho fLine fRabi]).'], nc);
    end

    n = size(R, 1);
    fprintf('\n=== REPLAY: no hardware, no measurement ===\n');
    fprintf('%d rows, %d columns.\n', n, nc);

    for i = n:-1:1
        q = struct('rho', R(i,1), 'pass', 'replay', 'good', false, 'why', '', ...
                   'fLine', R(i,2), 'fLineErr', opts.replayFLineErr, ...
                   'fwhm', NaN, 'depth', NaN, ...
                   'B', NaN, 'tiltDeg', NaN, 'aligned', true, 'fPred', NaN, ...
                   'windowLo', NaN, 'windowHi', NaN, 'neighbourGap', NaN, ...
                   'fRabi', NaN, 'fRabiErr', NaN, 'contrast', NaN, 'r2', NaN, ...
                   'r2Ceiling', NaN, 'r2Deficit', NaN, ...
                   'chi2', NaN, 'residAutocorr', NaN, 'piNs', NaN, 'piErrNs', NaN, ...
                   'rabiSeconds', NaN, 'shapeWarn', false, 'skippedMHz', NaN);
        switch nc
            case 8
                q.fwhm      = R(i,3);
                q.fRabi     = R(i,4);
                q.piNs      = R(i,5);
                q.contrast  = R(i,6)/100;
                q.r2        = R(i,7);
                q.r2Ceiling = R(i,8);
                q.r2Deficit = R(i,8) - R(i,7);
                q.fRabiErr  = opts.replayFErrFrac * abs(R(i,4));
            case 4
                q.fRabi    = R(i,3);
                q.fRabiErr = R(i,4);
            case 3
                q.fRabi    = R(i,3);
                q.fRabiErr = opts.replayFErrFrac * abs(R(i,3));
        end
        % A row is usable if it carries the three numbers the fit needs. A NaN fRabi is how a
        % rejected point appears in the printed table, so pasting the whole table including
        % its failures gives the same usable count the run reported.
        q.good = isfinite(q.rho) && isfinite(q.fLine) && isfinite(q.fRabi) && ...
                 isfinite(q.fRabiErr) && q.fRabi > 0;
        if ~q.good; q.why = 'replayed row has no usable fRabi (it was rejected in the original run)'; end
        P(i) = q;
    end
    PC = P;

    if nc < 4
        fprintf(['\nfRabiErr ASSUMED at %.2g%% of fRabi, and fLineErr at %.2g MHz -- the printed\n' ...
                 'table does not carry them. The peak position is insensitive to this; its\n' ...
                 'error bar and Q are not. The sensitivity block below measures which.\n'], ...
                100*opts.replayFErrFrac, opts.replayFLineErr);
    elseif nc == 8
        fprintf(['\nfRabiErr ASSUMED at %.2g%% of fRabi (the table has R^2 but not the error\n' ...
                 'bar). fLineErr assumed at %.2g MHz. See the sensitivity block below.\n'], ...
                100*opts.replayFErrFrac, opts.replayFLineErr);
    end
    printHeader('REPLAYED');
    for i = 1:numel(P); printRow(P(i)); end
end

function replaySensitivity(P, good, opts, fTarget, rhoStar, relAll)
%REPLAYSENSITIVITY Refit at three error assumptions and report the spread.
%
% The honest way to present an assumed error bar: change it by a factor of four and show what
% moves. If the answer does not move, the assumption did not matter and saying so is worth
% more than a caveat. If it does move, that is the real uncertainty and it should be quoted
% instead of the fit's.
    fracs = [0.005 0.01 0.02 0.04];
    fT = nan(size(fracs)); rS = nan(size(fracs));
    G = P(good);
    fL = [G.fLine]; wR = [G.fRabi];
    for k = 1:numel(fracs)
        try
            fitK = resonatorLorentzianFit(fL, wR, fracs(k)*abs(wR), ...
                       struct('baseline', opts.baselineInFit, 'verbose', false));
            ft = pickTarget(fitK);
            fT(k) = ft;
            if relAll.ok
                rS(k) = relAll.invert(ft, 0);
            end
        catch
        end
    end
    ok = isfinite(fT);
    fprintf('\n--- how much did the ASSUMED error bar matter? ---\n');
    fprintf('  %-14s %12s %12s\n', 'fRabiErr/fRabi', 'freq/MHz', 'Brho/G');
    for k = 1:numel(fracs)
        mark = '';
        if abs(fracs(k) - opts.replayFErrFrac) < 1e-12; mark = '   <- used'; end
        fprintf('  %12.1f%% %12.2f %12.3f%s\n', 100*fracs(k), fT(k), rS(k), mark);
    end
    if sum(ok) >= 2
        sf = max(fT(ok)) - min(fT(ok));
        sr = max(rS(isfinite(rS))) - min(rS(isfinite(rS)));
        if sf < 0.5
            fprintf(['\n  The frequency moves %.2f MHz across an eightfold change in the assumed\n' ...
                     '  error bar, and the field %.3f G. The assumption did not matter: an equally\n' ...
                     '  weighted fit to points of similar quality is the right fit, and %.2f MHz is\n' ...
                     '  far inside the linewidth. Use these numbers.\n'], sf, sr, sf);
        else
            fprintf(['\n  The frequency moves %.2f MHz and the field %.3f G across an eightfold\n' ...
                     '  change in the assumed error bar. That is NOT negligible -- it means one or\n' ...
                     '  two points are carrying the fit, so the weighting decides the answer. Quote\n' ...
                     '  %.2f MHz as the uncertainty rather than the fit''s own, and consider\n' ...
                     '  re-measuring a few field points near the peak.\n'], sf, sr, sf);
        end
    end
    fprintf('  (the answer above was computed at %.2g%%)\n', 100*opts.replayFErrFrac);
end

% ======================================================================
% ONE FIELD POINT
% ======================================================================

function p = measureAt(h, rho, freqs, rel, tag, opts)
% Set the field, find the lowest line, drive it, read the Rabi frequency.
    p = struct('rho', rho, 'pass', tag, 'good', false, 'why', '', ...
               'fLine', NaN, 'fLineErr', NaN, 'fwhm', NaN, 'depth', NaN, ...
               'B', NaN, 'tiltDeg', NaN, 'aligned', false, 'fPred', NaN, ...
               'windowLo', NaN, 'windowHi', NaN, 'neighbourGap', NaN, ...
               'fRabi', NaN, 'fRabiErr', NaN, 'contrast', NaN, 'r2', NaN, ...
               'r2Ceiling', NaN, 'r2Deficit', NaN, ...
               'chi2', NaN, 'residAutocorr', NaN, 'piNs', NaN, 'piErrNs', NaN, ...
               'rabiSeconds', NaN, 'shapeWarn', false, 'skippedMHz', NaN);

    % windowLo/windowHi/neighbourGap ARE SEEDED HERE, and that is the whole point of these
    % three lines. They are only assigned in the narrow branch below, i.e. only on the fine
    % pass -- so before this, coarse points had 25 fields and fine points had 28, and
    % P = [PC PF] died with
    %     "Number of fields in structure arrays being concatenated do not match"
    % after the entire scan had finished. Half an hour of measurement, thrown away at the
    % concatenation. Every field a struct array might carry has to exist in the seed, even
    % when the value is never filled in.

    h.setBrho(rho);
    pause(opts.settle);

    % ---- the ESR window --------------------------------------------------
    narrow = isempty(freqs);
    halfLo = opts.esrNarrowHalf; halfHi = halfLo; nbGap = Inf;
    if narrow
        p.fPred = rel.a*rho + rel.b;
        [halfLo, halfHi, nbGap] = narrowWindow(p.fPred, rel, opts);
        freqs = round((p.fPred - halfLo)*2)/2 : opts.esrStepFine : ...
                round((p.fPred + halfHi)*2)/2;
        p.windowLo = halfLo; p.windowHi = halfHi; p.neighbourGap = nbGap;
    end

    % ---- the ESR ---------------------------------------------------------
    % evalc so a 200-line survey report per field point does not bury the scan table, but
    % KEPT so it can be printed when something goes wrong. A swallowed explanation has
    % cost a debugging session before.
    bo = opts.baselineOrder;
    if narrow; bo = 1; end      % a cubic baseline across a 30 MHz window with one line in
                                % it will happily absorb the line itself
    sv = []; txt = '';
    try
        txt = evalc(['sv = calibFieldSurvey(struct(''freqs'', freqs, ''amplitude'', ' ...
                     'opts.esrAmplitude, ''repeats'', opts.esrRepeats, ''averages'', ' ...
                     'opts.esrAverages, ''detectionMHz'', opts.esrDetectionMHz, ' ...
                     '''isTracking'', opts.isTracking, ''freshExperiment'', ' ...
                     'opts.freshExperiment, ''baselineOrder'', bo, ''jsonLocation'', ' ...
                     'opts.jsonLocation, ''plotIt'', false, ''verboseRefine'', false, ' ...
                     '''D0'', opts.D0, ''newOrderEveryAverage'', opts.esrNewOrderEveryAverage));']);
    catch err
        p.why = ['ESR failed: ' err.message];
        printIndented(txt);
        return
    end
    if ~isstruct(sv) || ~isfield(sv, 'fit') || ~isstruct(sv.fit) || isempty(sv.fit.lines)
        p.why = 'ESR fitted no lines';
        printIndented(txt);
        return
    end

    L = sv.fit.lines;
    if narrow
        % In a narrow window there is a PREDICTION, so use it: take the line closest to
        % where the coarse relation says the line should be.
        %
        % Not the deepest, which is what this did before. The off-axis group is three
        % overlapping orientations, so where the window reaches it that blob is DEEPER than
        % the single aligned line -- and "take the deepest" then confidently returns the
        % neighbour. Not the lowest either: the lowest in a wide-ish window can be a
        % hyperfine satellite or a baseline wiggle the fitter promoted to a line. Closest to
        % the prediction is the only rule that says "the line we came for".
        [~, iPick] = min(abs([L.centre] - p.fPred));
    else
        % Coarse: the line the FIELD FIT assigned to the lower aligned transition, not simply
        % the lowest fitted centre. The resonator's own response at the low edge of the sweep
        % used to be fitted as a broad "line" (2767 MHz, FWHM 33, at Brho 31.5 on setup 7)
        % below the real aligned line at 2786, and min([L.centre]) picked it. esrFitRobust now
        % recognises that dip as the resonance and keeps it out of L (sv.resonance), so this
        % is the SECOND guard rather than the only one: identify the line by physics, fall
        % back to position only when the field fit did not converge (and the gate below then
        % rejects it, as before). Same rule as calibFieldAlign.
        [~, iPick] = min([L.centre]);
        iAl = alignedLowerLine(L, sv);
        if ~isempty(iAl) && iAl ~= iPick
            p.skippedMHz = L(iPick).centre;
            iPick = iAl;
        end
    end
    p.fLine = L(iPick).centre;  p.fLineErr = L(iPick).centreErr;
    p.fwhm  = L(iPick).fwhm;    p.depth    = L(iPick).depth;

    % ---- ESR gates -------------------------------------------------------
    if p.fwhm < opts.minFwhm
        p.why = sprintf('line FWHM %.1f MHz is below %.1f -- that is a noise spike, not a line', ...
                        p.fwhm, opts.minFwhm);
        return
    end
    if min(abs(p.fLine - [min(freqs) max(freqs)])) < opts.edgeGuard
        p.why = sprintf('line at %.1f MHz is within %.1f MHz of a sweep edge, so its centre and width are unreliable', ...
                        p.fLine, opts.edgeGuard);
        return
    end

    if narrow
        % CONTINUITY, in place of the field fit. A line that has jumped to a different
        % orientation lands tens of MHz from where the coarse relation predicts.
        % Half the distance to the neighbouring orientation, which is what the gate is
        % actually for -- and never wider than the window, since a line outside it could not
        % have been fitted anyway. The old narrowHalf/2 was unrelated to where the neighbour
        % is, so it was loose at low field and tight at high field, both by accident.
        tol = min(0.5*nbGap, max(halfLo, halfHi));
        if abs(p.fLine - p.fPred) > tol
            p.why = sprintf(['line at %.1f MHz is %.1f MHz from the predicted %.1f (tolerance ' ...
                             '%.1f) -- probably a different orientation'], ...
                            p.fLine, abs(p.fLine - p.fPred), p.fPred, tol);
            return
        end
        p.aligned = true;       % by continuity with a point that WAS checked properly
    else
        % The full check, as calibB1Peak: is the lowest line really the aligned
        % orientation, and was the aligned pair fully inside the sweep?
        if ~(isfield(sv, 'solve') && isstruct(sv.solve) && isfield(sv.solve, 'ok') && sv.solve.ok)
            p.why = 'no field fit, so the lowest line cannot be attributed to an orientation';
            return
        end
        p.B = sv.solve.B; p.tiltDeg = sv.solve.tiltDeg;
        mc = [sv.solve.matched.centre]; ma = [sv.solve.matched.axis];
        [dd, j] = min(abs(mc - p.fLine));
        p.aligned = dd < 1e-6 && ma(j) == sv.solve.alignedAxisIdx;
        if ~p.aligned
            p.why = 'the lowest line is not the aligned orientation';
            return
        end
        if ~isfield(sv.solve, 'outerVisible') || ~sv.solve.outerVisible
            p.why = 'the aligned pair is not fully inside the sweep';
            return
        end
    end

    % ---- the Rabi, at the frequency just measured ------------------------
    r = runRabiAt(p.fLine, opts);
    p.fRabi = r.f; p.fRabiErr = r.fErr; p.contrast = r.contrast;
    p.r2 = r.r2; p.chi2 = r.chi2; p.residAutocorr = r.residAutocorr;
    p.r2Ceiling = r.r2Ceiling; p.r2Deficit = r.r2Deficit;
    p.piNs = r.piNs; p.piErrNs = r.piErrNs; p.rabiSeconds = r.seconds;
    p.shapeWarn = r.shapeWarn;
    if ~r.ok
        p.why = r.why;
        return
    end

    p.good = true;
end

function iPick = alignedLowerLine(L, sv)
% Index into L of the line the field fit matched to the LOWER line of the aligned pair, or []
% when the fit did not converge or did not match that transition. Lower specifically: the
% scan follows the ms=0 -> ms=-1 aligned line, and if only the upper aligned line were matched,
% falling back to "lowest" (and letting the gate reject it) is right, while taking the upper
% one would silently follow the wrong transition.
    iPick = [];
    if ~(isfield(sv, 'solve') && isstruct(sv.solve) && isfield(sv.solve, 'ok') && sv.solve.ok ...
            && isfield(sv.solve, 'matched') && ~isempty(sv.solve.matched) ...
            && isfield(sv.solve, 'alignedAxisIdx') && isfield(sv.solve, 'alignedLines'))
        return
    end
    M = sv.solve.matched;
    lowPred = min(sv.solve.alignedLines);
    k = find([M.axis] == sv.solve.alignedAxisIdx & abs([M.predicted] - lowPred) < 1e-9, 1);
    if isempty(k); return; end
    [dd, j] = min(abs([L.centre] - M(k).centre));
    if dd < 1e-6; iPick = j; end
end

function w = medianFwhmOf(P)
% The typical fitted linewidth over the points whose line position was trusted. Median, not
% mean: one point that latched onto a blend of two orientations reports a very wide line and
% would drag a mean upward, which would then widen every fine window for no reason.
    v = [];
    for i = 1:numel(P)
        if isfinite(P(i).fwhm) && P(i).aligned; v = [v P(i).fwhm]; end             %#ok<AGROW>
    end
    if isempty(v); w = NaN; else; w = median(v); end
end

function r = runRabiAt(freq, opts)
% One Rabi at one frequency, with the walk standard, fitted and gated.
    r = struct('ok', false, 'why', '', 'f', NaN, 'fErr', NaN, 'contrast', NaN, ...
               'r2', NaN, 'r2Ceiling', NaN, 'r2Deficit', NaN, ...
               'chi2', NaN, 'residAutocorr', NaN, 'piNs', NaN, ...
               'piErrNs', NaN, 'seconds', NaN, 'shapeWarn', false);

    b = AWGExpRabi;
    b.loadJsonParams('all', opts.jsonLocation);
    % nChannels BEFORE amplitude, always: set.nChannels re-enters set.amplitude, which
    % would otherwise add amplitudeDiff twice.
    b.nChannels = opts.nChannels;
    if opts.nChannels == 2
        b.amplitude = opts.rabiAmplitude + [0 0];
    else
        b.amplitude = opts.rabiAmplitude;
    end
    b.frequency  = freq;
    b.repeats    = opts.repeats;
    b.averages   = opts.averages;
    b.isTracking = opts.isTracking;
    b.balancedSequence = 1;
    b.detectionDuration = opts.rabiDetectionMHz;
    b.referenceDetectionDuration = b.detectionDuration;
    b.shouldAutosave = false;               % ~7 s per run writing a .mat nothing reads
    b.tau = opts.tau;
    b.useSequenceWalk    = opts.useSequenceWalk;
    b.passesPerRun       = opts.passesPerRun;
    b.walkWarmupTriggers = opts.walkWarmupTriggers;
    b.resyncEveryPass    = opts.resyncEveryPass;
    if ~isempty(opts.newOrderEveryAverage)
        b.newOrderEveryAverage = logical(opts.newOrderEveryAverage);
    end

    txt = '';
    tId = tic;      % ID FORM. A bare tic is one global stopwatch that any tic inside
                    % b.run resets, and that reported 1.9 s for a 68 s run once already.
    try
        txt = evalc('b.run');
    catch err
        r.seconds = toc(tId);
        r.why = ['Rabi run failed: ' err.message];
        printIndented(txt);
        return
    end
    r.seconds = toc(tId);

    try
        o = rabiFitRobust(struct('exp', b, 'plotIt', false, 'verbose', false, ...
                                 'nFreq', opts.nFreq));
    catch err
        r.why = ['Rabi fit threw: ' err.message];
        return
    end
    if ~isstruct(o) || ~isfield(o, 'ok') || ~o.ok
        r.why = 'Rabi fit failed';
        if isstruct(o) && isfield(o, 'why'); r.why = ['Rabi fit failed: ' o.why]; end
        return
    end

    r.f = o.f; r.fErr = o.fErr; r.contrast = o.contrast; r.r2 = o.r2;
    r.chi2 = o.chi2Reduced; r.residAutocorr = o.residAutocorr;
    r.piNs = 1e3*o.pulses(2); r.piErrNs = 1e3*o.pulsesErr(2);

    % ---- how high could R^2 possibly have been? --------------------------
    % If the model is right, the residuals are the error bars, so the smallest achievable
    % residual sum is sum(e^2) and
    %
    %       R^2_ceiling = 1 - sum(e^2) / sum((y - mean(y))^2)
    %
    % which is set entirely by the signal-to-noise ratio of the sweep and not at all by the
    % fit. Comparing the achieved R^2 with this is the difference between "the model is
    % wrong" and "there were not many photons": a fit at R^2 = 0.82 against a ceiling of
    % 0.84 is as good as the data permit, while the same 0.82 against a ceiling of 0.99 is a
    % real failure. Two numbers say what one cannot, and the absolute R^2 gate that shipped
    % here could not tell those cases apart -- it rejected both.
    try
        ssTot = sum((o.y - mean(o.y)).^2);
        r.r2Ceiling = 1 - sum(o.e.^2)/max(ssTot, eps);
    catch
        r.r2Ceiling = NaN;
    end
    r.r2Deficit = r.r2Ceiling - r.r2;

    % ---- gates -----------------------------------------------------------
    if r.contrast < opts.minContrast
        r.why = sprintf(['Rabi contrast %.2f %% is below %.2f %% -- check the MW actually ' ...
                         'reached the sample'], 100*r.contrast, 100*opts.minContrast);
        return
    end
    if r.r2 < opts.minR2
        % SAY WHETHER THE GATE OR THE SIGNAL IS THE PROBLEM. "R^2 0.82 is below 0.90" leaves
        % you no way to tell a bad fit from a quiet one, and the answer changes what you do:
        % lower the gate, or raise the repeats.
        if isfinite(r.r2Ceiling) && r.r2Ceiling < opts.minR2
            r.why = sprintf(['Rabi R^2 %.3f is below the minR2 gate of %.2f -- but the NOISE ' ...
                             'CEILING here is only %.3f, so no fit could have passed. This is ' ...
                             'a signal-to-noise limit at %.2f%% contrast, not a fitting ' ...
                             'failure. Lower minR2 (or set it to 0 and rely on maxFErrFrac), ' ...
                             'or raise repeats'], ...
                            r.r2, opts.minR2, r.r2Ceiling, 100*r.contrast);
        else
            r.why = sprintf(['Rabi R^2 %.3f is below %.2f, and the noise allowed %.3f -- so ' ...
                             'this one really is a poor fit, not just a quiet one. A wrong ' ...
                             'shape is the usual cause; try nFreq = 2'], ...
                            r.r2, opts.minR2, r.r2Ceiling);
        end
        return
    end
    if r.fErr/max(r.f, eps) > opts.maxFErrFrac
        r.why = sprintf('Rabi frequency %.2f +/- %.2f MHz is only %.0f%% determined', ...
                        r.f, r.fErr, 100*(1 - r.fErr/max(r.f, eps)));
        return
    end

    % Autocorrelated residuals mean the SHAPE is wrong -- most often a second orientation
    % being driven off resonance, which adds a component at sqrt(Omega^2 + Delta^2). That
    % biases a single-frequency fit, but it does NOT make the point worthless, and throwing
    % away every point with hyperfine structure in it would empty the scan. Flagged, kept,
    % and the fix (nFreq = 2) is named at the end.
    r.shapeWarn = abs(r.residAutocorr) > opts.maxResidAutocorr;
    r.ok = true;
end

% ======================================================================
% THE LINE-VS-FIELD RELATION AND ITS INVERSE
% ======================================================================

function rel = lineVsField(P, opts, tag, GAMMA)
% Which points feed the coil relation, and then coilLineRelation does the fitting.
%
% An ESR that passed its gates counts here even if its Rabi later failed: the line POSITION
% is measured by the ESR alone, and the Rabi has no bearing on it. Using only fully-good
% points would throw away information for no reason -- and the relation is what points the
% narrow fine-pass windows, so it wants every position it can get.
    use = arrayfun(@(q) isfinite(q.fLine) && q.aligned, P);
    if ~any(use)
        rel = coilLineRelation([], [], [], struct('tag', tag, 'gamma', GAMMA, ...
                                                  'coilRatio', opts.coilRatio));
        return
    end
    rel = coilLineRelation([P(use).rho], [P(use).fLine], [P(use).fLineErr], ...
                           struct('tag', tag, 'gamma', GAMMA, 'coilRatio', opts.coilRatio));
end

function [rho, how] = coarsePeakRho(P, rel)
% Where to centre the fine pass. The Lorentzian is overkill on a coarse pass and can be
% rejected outright there, so this uses the cheapest thing that works and says which.
    rho = NaN; how = '';
    good = logical([P.good]);
    if ~any(good); how = 'no usable coarse points'; return; end
    G = P(good);
    w = [G.fRabi]; r = [G.rho]; f = [G.fLine];

    if numel(G) >= 5
        fit = resonatorLorentzianFit(f, w, [G.fRabiErr], ...
                  struct('verbose', false, 'retryNoBaseline', false));
        [fT, ~, howT] = pickTarget(fit);
        if isfinite(fT) && rel.ok
            rho = rel.invert(fT, 0);
            how = ['coarse ' howT];
            % Only trust it as a CENTRE if it is inside the coarse fields. Outside, the
            % best measured point is a better place to put nine more measurements than an
            % extrapolation is.
            if rho >= min(r) && rho <= max(r); return; end
            how = [how ', but outside the coarse range'];
        end
    end
    [~, i] = max(w);
    rho = r(i);
    if isempty(how); how = 'best measured point'; else; how = [how ' -> best measured point']; end
end

function [f, fErr, how] = pickTarget(fit)
% The frequency to aim at. The Lorentzian when it is usable; the parabola vertex when it
% is not but the vertex is inside the scan. Never a fit that has been flagged unusable
% without saying so.
    f = NaN; fErr = NaN; how = 'nothing usable';
    if ~isstruct(fit); return; end
    if isfield(fit, 'usable') && fit.usable
        f = fit.f0; fErr = fit.f0Err; how = 'Lorentzian centre'; return
    end
    if isfield(fit, 'parabola') && fit.parabola.ok && fit.parabola.inRange
        f = fit.parabola.f0; fErr = fit.parabola.f0Err;
        how = 'parabola vertex -- the Lorentzian was not usable';
        return
    end
    if isfield(fit, 'f0') && isfinite(fit.f0)
        f = fit.f0; fErr = fit.f0Err;
        how = 'Lorentzian centre, FLAGGED UNUSABLE -- treat as a direction, not a value';
    end
end

% ======================================================================
% REPORTING
% ======================================================================

function printHeader(tag)
    fprintf('\n%s\n', tag);
    % R2max is the ceiling the noise imposes. Printing R^2 without it invites reading a
    % low value as a bad fit when it is a quiet one.
    fprintf('%8s %10s %7s %9s %9s %7s %7s %7s  %s\n', ...
            'Brho/G', 'line/MHz', 'FWHM', 'fRabi/MHz', 'pi/ns', 'contr%', 'R^2', 'R2max', 'note');
end

function printRow(p)
    note = '';
    if p.good
        note = 'ok';
        if p.shapeWarn; note = 'ok, residuals autocorrelated'; end
    else
        note = ['REJECTED: ' p.why];
    end
    if isfield(p, 'skippedMHz') && isfinite(p.skippedMHz)
        note = sprintf('%s (skipped unmatched line at %.1f)', note, p.skippedMHz);
    end
    fprintf('%8.3f %10.2f %7.2f %9.3f %9.1f %7.2f %7.3f %7.3f  %s\n', ...
            p.rho, p.fLine, p.fwhm, p.fRabi, p.piNs, 100*p.contrast, p.r2, p.r2Ceiling, note);
end

function listRejects(P)
    bad = P(~logical([P.good]));
    if isempty(bad); return; end
    fprintf('\n--- rejected points, and why ---\n');
    for i = 1:numel(bad)
        fprintf('  Brho %7.3f (%s): %s\n', bad(i).rho, bad(i).pass, bad(i).why);
    end
end

function checkAgainstOldMethod(out, fTarget)
% calibB1Peak's FWHM^2 metric put the peak at 2823.1 MHz on this setup. That was a proxy
% and this is a measurement, so they are allowed to differ -- but a LARGE difference is
% information, and quietly replacing one number with another that disagrees by 15 MHz
% without remarking on it is how a setup ends up with two contradictory calibrations.
    fOld = 2823.1;
    d = fTarget - fOld;
    fprintf('\n--- against the earlier linewidth scan ---\n');
    fprintf('  calibB1Peak (FWHM^2 proxy) said %.1f MHz; this says %.2f MHz. Difference %+.2f MHz.\n', ...
            fOld, fTarget, d);
    if abs(d) < 3
        fprintf(['  They agree. The proxy was monotone in B1 after all, and the two methods\n' ...
                 '  now support each other.\n']);
    else
        fprintf(['  They disagree by more than 3 MHz. Believe THIS one: the Rabi frequency is\n' ...
                 '  B1 with no saturation model in between, while a CW width also carries\n' ...
                 '  inhomogeneous broadening and hyperfine structure. Worth noting in the log\n' ...
                 '  that the linewidth number is superseded rather than leaving both on file.\n']);
    end
    if isfield(out, 'fit') && isstruct(out.fit) && isfield(out.fit, 'usable') && out.fit.usable
        fprintf('  This scan also gives Q = %.0f +/- %.0f, which the linewidth scan could not.\n', ...
                out.fit.Q, out.fit.QErr);
    end
end

function reportVerify(pv, fit, fTarget)
    fprintf('\n--- verification at the chosen field ---\n');
    if ~pv.good
        fprintf(['  THE VERIFICATION POINT WAS REJECTED: %s\n' ...
                 '  The field has still been set. Re-run one Rabi by hand before trusting it.\n'], ...
                pv.why);
        return
    end
    fprintf('  line landed at %.2f MHz, target was %.2f (%+.2f MHz off)\n', ...
            pv.fLine, fTarget, pv.fLine - fTarget);
    if isfield(fit, 'wMax') && isfinite(fit.wMax)
        d = pv.fRabi - fit.wMax;
        s = abs(d)/max(sqrt(pv.fRabiErr^2 + fit.wMaxErr^2), eps);
        fprintf('  fRabi %.3f +/- %.3f MHz against the predicted %.3f +/- %.3f (%.1f sigma)\n', ...
                pv.fRabi, pv.fRabiErr, fit.wMax, fit.wMaxErr, s);
        if s > 3
            fprintf(['    *** the confirmation disagrees with the fit by %.1f sigma. The most\n' ...
                     '    likely causes, in order: the drive level moved between the scan and\n' ...
                     '    now; the response is not a single Lorentzian so the peak height is\n' ...
                     '    wrong even though the position is right; or the field did not settle.\n' ...
                     '    The POSITION is confirmed by the line landing where predicted, so the\n' ...
                     '    chosen field is probably still right.\n'], s);
        else
            fprintf('    consistent -- the response curve predicts the point it was not fitted to.\n');
        end
    end
    fprintf('  pi pulse here: %.1f +/- %.1f ns. That is the number the pulsed calibrations want.\n', ...
            pv.piNs, pv.piErrNs);
end

function printHandover(out)
    fprintf('\n=== hand over to the microwave calibrations ===\n');
    if ~out.ok
        fprintf(['  NOT SETTLED. The response curve could not be pinned down; read the flags\n' ...
                 '  above and extend the field scan before running anything that depends on\n' ...
                 '  B1. Nothing downstream should use these numbers yet.\n']);
        return
    end
    fprintf(['  frequency  %.2f MHz   drive here. level3.m cell 3.5 writes it into\n' ...
             '                          expProperties.json (key "frequency") with expPropsWrite.\n'], ...
            out.resonantFrequency);
    fprintf(['  Brho       %.3f G       cell 3.4c sets it and saves it to expProperties.json\n' ...
             '                          (a coil setting, not a json field).\n'], out.rhoOptimum);
    if isfield(out, 'piNsAtOptimum') && isfinite(out.piNsAtOptimum)
        fprintf('  pi pulse about %.1f ns\n', out.piNsAtOptimum);
    end
    anyShape = any([out.points.shapeWarn]);
    if anyShape
        fprintf(['\n  Some points had autocorrelated Rabi residuals, which usually means a\n' ...
                 '  second orientation is being driven off resonance and adding a component at\n' ...
                 '  sqrt(Omega^2 + Delta^2). Re-run with nFreq = 2 to fit both and see whether\n' ...
                 '  the peak position moves; if it does not, ignore it.\n']);
    end
end

% ======================================================================
% PLOT
% ======================================================================

function drawIt(P, good, fit, rel, fTarget, rhoStar, opts)
% Three axes, because there are three different quantities: the antenna response against
% frequency (the measurement), the line position against field (the conversion), and the
% residuals of that conversion (the only place a wrong-orientation point is visible).
%
% They are NOT stacked on shared y scales -- that would imply a relationship between the
% axes that does not exist.
%
% What was added over the two-panel version, and why each one earns its space:
%
%   THE FWHM, DRAWN. Q was in the title as a number. The resonator's width is the thing you
%   are actually measuring, and a half-maximum line with the FWHM shaded says in one glance
%   what "Q = 94" does not: whether the scan spanned the width or only guessed at it. The
%   grader's TEST B exists because a scan narrower than ~35 MHz cannot determine the width
%   at all, and this makes that visible instead of implied.
%
%   THE ANSWERS, LABELLED. f0 and the optimum B_rho were a dashed line and an unlabelled
%   star. They are the two outputs of the whole calibration and they now carry their values
%   on the axes.
%
%   THE RESIDUAL PANEL. TEST D in testB1RabiChain makes the point directly: one point that
%   tracked a different NV orientation bends the coil line without breaking the slope test,
%   and reduced chi2 is what catches it (0.6 clean, 2745 with the rogue point). chi2 was
%   printed but the residuals were not plotted, so there was no way to see WHICH point. Now
%   there is.
    C = struct('surface', [252 252 251]/255, 'ink', [11 11 11]/255, ...
               'ink2', [82 81 78]/255, 'data', [42 120 214]/255, ...
               'fit', [235 104 52]/255, 'drop', [140 138 133]/255, ...
               'band', [246 224 210]/255, 'ok', [0 115 89]/255);

    G = P(good); B = P(~good);
    f = [G.fLine]; w = [G.fRabi]; we = [G.fRabiErr];
    y = w.^2; ye = 2*w.*we;

    fig = figure('Name', 'calibB1Rabi', 'NumberTitle', 'off', 'Color', C.surface);

    % =================================================================
    % 1. the antenna response
    % =================================================================
    ax1 = subplot(3, 1, 1); hold(ax1, 'on');

    usable = isstruct(fit) && isfield(fit, 'usable') && fit.usable;
    fwhm = NaN;
    if usable && isfinite(fit.Q) && fit.Q > 0; fwhm = fit.f0/fit.Q; end

    % The FWHM band goes down first so the data draws over it.
    hBand = [];
    if usable && isfinite(fwhm)
        yl0 = [0, max(y + ye)*1.15];
        hBand = patch(ax1, fit.f0 + fwhm/2*[-1 1 1 -1], [yl0(1) yl0(1) yl0(2) yl0(2)], ...
                      C.band, 'EdgeColor', 'none');
        ylim(ax1, yl0);
    end

    hE = errorbar(ax1, f, y, ye);
    set(hE, 'LineStyle', 'none', 'Color', C.data, 'LineWidth', 0.75);
    try; set(hE, 'CapSize', 4); catch; end
    hD = plot(ax1, f, y, 'o', 'MarkerSize', 5, 'LineWidth', 1, 'Color', C.data, ...
              'MarkerFaceColor', C.surface);
    hF = []; hH = []; hT = [];
    if usable && isfield(fit, 'model') && isa(fit.model, 'function_handle')
        pad = 0.08*(max(f) - min(f));
        ff = linspace(min(f) - pad, max(f) + pad, 800);
        hF = plot(ax1, ff, fit.model(ff), '-', 'LineWidth', 2, 'Color', C.fit);
    end
    yl = ylim(ax1);
    if usable && isfinite(fwhm)
        % The half-maximum level is what DEFINES the width, so draw it rather than leaving
        % the reader to trust Q.
        hH = plot(ax1, fit.f0 + fwhm/2*[-1 1], (fit.wMax^2/2)*[1 1], '-', ...
                  'Color', C.fit, 'LineWidth', 1.2);
        text(ax1, fit.f0, fit.wMax^2/2, sprintf(' FWHM %.1f MHz ', fwhm), ...
             'Color', C.fit, 'FontSize', 8, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'BackgroundColor', C.surface);
    end
    if isfinite(fTarget)
        hT = plot(ax1, [fTarget fTarget], yl, '--', 'LineWidth', 1.2, 'Color', C.ok);
        text(ax1, fTarget, yl(1), sprintf(' %.1f MHz ', fTarget), 'Color', C.ok, ...
             'FontWeight', 'bold', 'FontSize', 8, 'VerticalAlignment', 'bottom');
        ylim(ax1, yl);
    end
    % Omega_max and t_pi are NOT annotated at the peak. They are already in the title, and a
    % text box at (f0, wMax^2) lands in the same corner as the legend and under the title --
    % three things fighting for the top of the axes. One place each.
    styleAxes(ax1, C);
    xlabel(ax1, 'NV line frequency  (MHz)', 'Color', C.ink2);
    ylabel(ax1, '\Omega^2  (MHz^2)', 'Color', C.ink2);
    % 'best' rather than a fixed corner: the Lorentzian peak can sit anywhere in the scan, so
    % which corner is free depends on the data.
    legendOf(ax1, {hD, hF, hH, hBand, hT}, ...
             {'measured \Omega^2', 'Lorentzian fit', 'half maximum', 'FWHM', ...
              'operating frequency'}, 'best');
    if usable
        title(ax1, sprintf(['f_0 = %.2f \\pm %.2f MHz     Q = %.0f \\pm %.0f     ' ...
                            'FWHM = %.1f MHz     \\Omega_{max} = %.2f MHz     ' ...
                            't_\\pi = %.1f \\pm %.1f ns'], ...
                     fit.f0, fit.f0Err, fit.Q, fit.QErr, fwhm, fit.wMax, ...
                     fit.piMinNs, fit.piMinErrNs), 'Color', C.ink, 'FontWeight', 'bold');
    else
        title(ax1, sprintf(['peak near %.2f MHz -- LORENTZIAN NOT USABLE. The scan probably ' ...
                            'does not span the width; see the report.'], fTarget), ...
              'Color', C.ink, 'FontWeight', 'bold');
    end

    % =================================================================
    % 2. the coil conversion
    % =================================================================
    ax2 = subplot(3, 1, 2); hold(ax2, 'on');
    hE2 = errorbar(ax2, [G.rho], [G.fLine], [G.fLineErr]);
    set(hE2, 'LineStyle', 'none', 'Color', C.data, 'LineWidth', 0.75);
    try; set(hE2, 'CapSize', 4); catch; end
    hD2 = plot(ax2, [G.rho], [G.fLine], 'o', 'MarkerSize', 5, 'LineWidth', 1, ...
               'Color', C.data, 'MarkerFaceColor', C.surface);
    % Rejected points shown, in grey, with a cross. A scan with three silent holes in it
    % looks like a scan that was planned that way.
    hB2 = []; nDrop = 0;
    if ~isempty(B)
        bf = [B.fLine]; br = [B.rho]; k = isfinite(bf); nDrop = sum(k);
        if any(k)
            hB2 = plot(ax2, br(k), bf(k), 'x', 'MarkerSize', 9, 'LineWidth', 1.2, ...
                       'Color', C.drop);
        end
    end
    hL2 = [];
    if rel.ok
        rr = linspace(min([P.rho]), max([P.rho]), 100);
        hL2 = plot(ax2, rr, rel.a*rr + rel.b, '-', 'LineWidth', 2, 'Color', C.fit);
    end
    hS2 = [];
    if isfinite(rhoStar) && isfinite(fTarget)
        yl2 = ylim(ax2);
        plot(ax2, [min([P.rho]) max([P.rho])], fTarget*[1 1], ':', 'Color', C.ok, ...
             'LineWidth', 1, 'HandleVisibility', 'off');
        hS2 = plot(ax2, rhoStar, fTarget, 'p', 'MarkerSize', 14, 'LineWidth', 1.2, ...
                   'Color', C.ok, 'MarkerFaceColor', C.ok);
        % THE answer of the whole calibration. It was an unlabelled star.
        text(ax2, rhoStar, fTarget, sprintf('  B_\\rho = %.3f G  ', rhoStar), ...
             'Color', C.ok, 'FontWeight', 'bold', 'FontSize', 9, ...
             'VerticalAlignment', 'top', 'BackgroundColor', C.surface);
        ylim(ax2, yl2);
    end
    styleAxes(ax2, C);
    xlabel(ax2, 'B_\rho requested  (G)', 'Color', C.ink2);
    ylabel(ax2, 'NV line frequency  (MHz)', 'Color', C.ink2);
    legendOf(ax2, {hD2, hL2, hS2, hB2}, ...
             {sprintf('used (%d)', numel(G)), 'straight-line fit', 'optimum', ...
              sprintf('rejected (%d)', nDrop)}, 'northeast');
    if rel.ok
        title(ax2, sprintf(['f = %.3f B_\\rho %+.2f MHz     optimum B_\\rho = %.3f G     ' ...
                            'slope %s'], rel.a, rel.b, rhoStar, ...
                     ternary(rel.slopeSane, 'sane', 'SUSPECT')), ...
              'Color', C.ink, 'FontWeight', 'bold');
    else
        title(ax2, 'coil relation NOT fitted -- see the report', 'Color', C.ink, ...
              'FontWeight', 'bold');
    end

    % =================================================================
    % 3. residuals of the coil line -- where a rogue orientation shows
    % =================================================================
    ax3 = subplot(3, 1, 3); hold(ax3, 'on');
    haveRes = rel.ok && isfield(rel, 'resid') && isfield(rel, 'rho') && ...
              numel(rel.resid) == numel(rel.rho) && ~isempty(rel.resid);
    if haveRes
        rx = rel.rho(:)'; rr2 = rel.resid(:)';
        yl3 = max(abs(rr2))*1.35; if ~(yl3 > 0); yl3 = 1; end
        patch(ax3, [min(rx) max(rx) max(rx) min(rx)], ...
              rel.rms*[-1 -1 1 1], C.band, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax3, [min(rx) max(rx)], [0 0], '-', 'Color', C.ink2, 'LineWidth', 0.75, ...
             'HandleVisibility', 'off');
        hR = plot(ax3, rx, rr2, 'o-', 'MarkerSize', 5, 'LineWidth', 1, 'Color', C.data, ...
                  'MarkerFaceColor', C.surface);
        % Name the worst point: with a rogue orientation it is the one to re-measure.
        [mx, im] = max(abs(rr2));
        if mx > 3*rel.rms
            plot(ax3, rx(im), rr2(im), 'o', 'MarkerSize', 11, 'LineWidth', 1.4, ...
                 'Color', C.fit, 'HandleVisibility', 'off');
            text(ax3, rx(im), rr2(im), sprintf('  %.2f G: %+.1f MHz  ', rx(im), rr2(im)), ...
                 'Color', C.fit, 'FontWeight', 'bold', 'FontSize', 8, ...
                 'VerticalAlignment', 'bottom', 'BackgroundColor', C.surface);
        end
        ylim(ax3, yl3*[-1 1]);
        legendOf(ax3, {hR}, {'residual'}, 'northeast');
        verdict = 'consistent with one orientation';
        if rel.chi2Reduced > 10
            verdict = 'chi2 >> 1: a point is on a DIFFERENT orientation, or the errors are too small';
        elseif rel.chi2Reduced > 3
            verdict = 'chi2 high: check the marked point';
        end
        title(ax3, sprintf('residuals: rms %.2f MHz, reduced \\chi^2 %.1f  --  %s', ...
                           rel.rms, rel.chi2Reduced, verdict), ...
              'Color', C.ink, 'FontWeight', 'bold');
    else
        title(ax3, 'no residuals to show (the coil line was not fitted)', ...
              'Color', C.ink, 'FontWeight', 'bold');
    end
    styleAxes(ax3, C);
    xlabel(ax3, 'B_\rho requested  (G)', 'Color', C.ink2);
    ylabel(ax3, 'f - fit  (MHz)', 'Color', C.ink2);

    try; set(fig, 'Position', [200 80 820 860]); catch; end
end

function v = ternary(c, a, b)
    if c; v = a; else; v = b; end
end

function legendOf(ax, handles, labels, loc)
% Handles and labels kept in CELLS and filtered together, so an absent series cannot shift
% the labels onto the wrong lines -- which is what [h1 h2 h3] does when one of them is [].
    keep = cellfun(@(h) ~isempty(h) && all(ishandle(h)), handles);
    if ~any(keep); return; end
    try
        legend(ax, [handles{keep}], labels(keep), 'Location', loc, 'Box', 'off', ...
               'FontSize', 8);
    catch
    end
end

function styleAxes(ax, C)
    grid(ax, 'on');
    set(ax, 'Color', C.surface, 'XColor', C.ink2, 'YColor', C.ink2, 'Box', 'off', ...
            'GridLineStyle', '-', 'Layer', 'top', 'TickDir', 'out', 'FontSize', 10);
    try; set(ax, 'GridAlpha', 0.18); catch; end
end

% ======================================================================
% SMALL HELPERS
% ======================================================================

function [halfLo, halfHi, gap] = narrowWindow(fPred, rel, opts)
%NARROWWINDOW The fine-pass ESR window around a predicted line: asymmetric, and capped on
% the side where the next orientation lives.
%
% ===================================================================
% WHY IT IS ASYMMETRIC
% ===================================================================
%
% The line being tracked is the ALIGNED pair's outer line -- the leftmost line in the whole
% spectrum when it sits below D. So the two sides are limited by completely different things:
%
%   AWAY FROM D   nothing is there. The window can be as generous as the fit wants, and it
%                 wants a couple of linewidths so the baseline has clean wings.
%   TOWARDS D     the off-axis group sits at
%
%                     gap  =  gamma*|B| - gamma*|B|/3  =  (2/3)*|D0 - fPred|
%
%                 because the aligned projection is |B| and the three off-axis projections
%                 are |B|/3. A window that reaches it fits the wrong dip -- and the off-axis
%                 group is three overlapping orientations, so it is DEEPER than the line we
%                 want.
%
% Note the gap is proportional to the field, so a single symmetric half-width is wrong across
% a field scan by construction. At 24 G the gap is 44.8 MHz and +/-30 is fine; at 12 G the gap
% is 22.4 MHz and the same window overshoots by 8 MHz. That is the bug this replaces: the
% fine pass caught the second dip from the left at the low end of the scan and not at the top.
%
% ===================================================================
% WHY esrNarrowHalf ALONE COULD NOT FIX IT
% ===================================================================
%
% esrNarrowHalf is a FLOOR:  want = max(esrNarrowHalf, narrowHalfFwhmFactor*fwhm).  Setting
% it to 7 changes nothing when 2.5 linewidths is 30, which is the whole point of the floor --
% a window narrower than the line gets the line rejected as baseline by esrFitRobust's own
% width test. So the fix is not a smaller floor, it is a CAP, and the cap has to come from
% the neighbour rather than from a number typed in. esrNarrowHalfMax is there as a hard
% override for when you want one anyway.
    fwhm = 12;
    if isfield(rel, 'narrowFwhm') && isfinite(rel.narrowFwhm) && rel.narrowFwhm > 0
        fwhm = rel.narrowFwhm;
    end
    floorHalf = opts.esrNarrowHalf;
    if isfield(rel, 'narrowFloor') && isfinite(rel.narrowFloor); floorHalf = rel.narrowFloor; end
    factor = opts.narrowHalfFwhmFactor;
    if isfield(rel, 'narrowFactor') && isfinite(rel.narrowFactor); factor = rel.narrowFactor; end
    guard = opts.esrNeighbourGuard;
    if isfield(rel, 'guard') && isfinite(rel.guard); guard = rel.guard; end
    hardMax = opts.esrNarrowHalfMax;
    if isfield(rel, 'narrowMax'); hardMax = rel.narrowMax; end

    want = max(floorHalf, factor*fwhm);
    explicit = ~isempty(hardMax) && isfinite(hardMax);
    if explicit; want = min(want, hardMax); end

    gap = (2/3)*abs(fPred - opts.D0);
    if ~isfinite(gap) || gap <= 0; gap = Inf; end
    capped = min(want, guard*gap);
    % Never let the NEIGHBOUR-DERIVED cap collapse to nothing: below about one linewidth
    % there is no line left to fit, and a window that cannot contain the line is worse than
    % one that risks the neighbour.
    %
    % But an EXPLICIT esrNarrowHalfMax is not overridden. If you ask for 7 MHz you get 7 MHz,
    % because silently returning 12 when you typed 7 is how the old floor behaved and is the
    % reason this function exists. It warns instead.
    if ~explicit
        capped = max(capped, 1.0*fwhm);
    elseif hardMax < fwhm
        warning('calibB1Rabi:narrowWindow', ...
                ['esrNarrowHalfMax = %.1f MHz is narrower than the %.1f MHz linewidth. ' ...
                 'esrFitRobust rejects a component whose FWHM exceeds the span, so ' ...
                 'expect the line to be thrown out as baseline. This is what you asked for; ' ...
                 'leave esrNarrowHalfMax empty and let esrNeighbourGuard size it instead.'], ...
                hardMax, fwhm);
    end

    if fPred <= opts.D0
        halfLo = want;   halfHi = capped;      % neighbours are above
    else
        halfLo = capped; halfHi = want;        % neighbours are below
    end
end

function freqs = esrWindowFor(rhoMax, opts, GAMMA)
% The wide ESR sweep: esrFreqs verbatim if you gave one, otherwise sized from the field
% range.
%
% The sweep has to hold BOTH lines of the aligned pair at the LARGEST field, because the
% field fit needs them and the orientation check needs the field fit. So the half-width is
% gamma*|B| plus a margin for the fitter to have some baseline on each side.
%
% D0 and coilRatio are used HERE and nowhere else that produces a number. They only decide
% how wide to look; every reported quantity comes from what was measured inside that
% window. So a coilRatio that is off by a few percent costs a slightly wider sweep and
% nothing else -- which is the only safe way to use a constant you do not trust.
    if ~isempty(opts.esrFreqs)
        freqs = sort(opts.esrFreqs(:)');
        return
    end
    half = GAMMA*rhoMax*opts.coilRatio + opts.esrMargin;
    freqs = round(opts.D0 - half) : opts.esrStepCoarse : round(opts.D0 + half);
end

function checkWindow(freqs, rhos, opts, GAMMA)
% Does this sweep actually reach the lines it needs to, at every field in the scan?
%
% Checked BEFORE the scan starts, not discovered forty minutes in. calibB1Peak's first live
% run lost the aligned pair at three of nine fields because the sweep stopped short at the
% top, and the points were accepted anyway -- they moved the reported peak 18 MHz and turned
% the coil fit into 0.13*Brho + 9.7. The gates inside measureAt would now reject those
% points, which is better, but a scan that spends its time collecting rejections is still a
% wasted scan.
    lo = min(freqs); hi = max(freqs);
    step = median(diff(freqs));
    fprintf('\nWide ESR sweep: %.1f to %.1f MHz, %g MHz steps, %d points%s.\n', ...
            lo, hi, step, numel(freqs), ...
            BooleanHelper.ifTrueElse(isempty(opts.esrFreqs), '', ' (from esrFreqs)'));

    rhoMax = max(rhos);
    needHi = opts.D0 + GAMMA*rhoMax*opts.coilRatio;
    needLo = opts.D0 - GAMMA*rhoMax*opts.coilRatio;
    fprintf(['It has to reach about %.0f MHz (upper line) and %.0f MHz (lower line) at the\n' ...
             'top field of %.2f G. Both lines of the aligned pair have to be inside the\n' ...
             'sweep or calibFieldSurvey cannot fit the field vector, and without that there\n' ...
             'is no check on which orientation the lowest line belongs to.\n'], ...
            needHi, needLo, rhoMax);

    shortAtTop = hi < needHi + opts.edgeGuard;
    shortAtBot = lo > needLo - opts.edgeGuard;
    if shortAtTop || shortAtBot
        % The largest field this sweep can actually support, so the choice is concrete:
        % widen the sweep or shorten the scan.
        rhoOkHi = (hi - opts.edgeGuard - opts.D0)/(GAMMA*opts.coilRatio);
        rhoOkLo = (opts.D0 - lo - opts.edgeGuard)/(GAMMA*opts.coilRatio);
        rhoOk = min(rhoOkHi, rhoOkLo);
        nBad = sum(rhos > rhoOk);
        fprintf(['\n*** THIS SWEEP IS TOO NARROW FOR THE TOP OF THE SCAN. It supports fields up\n' ...
                 '    to about %.2f G, and %d of your %d coarse points are above that -- they\n' ...
                 '    will be rejected for "the aligned pair is not fully inside the sweep",\n' ...
                 '    which costs about %.0f minutes for nothing. Either\n' ...
                 '      esrFreqs = %.0f : %g : %.0f       (widen the sweep)\n' ...
                 '    or stop the field scan at Brho = %.2f.\n'], ...
                rhoOk, nBad, numel(rhos), 2*nBad, ...
                round(opts.D0 - GAMMA*rhoMax*opts.coilRatio - opts.esrMargin), ...
                opts.esrStepCoarse, ...
                round(opts.D0 + GAMMA*rhoMax*opts.coilRatio + opts.esrMargin), rhoOk);
    end

    % A step comparable with a linewidth leaves too few points on each line to fit it. The
    % lines here are 8 to 25 MHz wide depending on drive level, so 2 MHz gives four to
    % twelve points across one, which is enough; 6 MHz would not be.
    if step > 4
        fprintf(['\n*** %g MHz STEPS ARE TOO COARSE. An ESR line here is 8 to 25 MHz wide, so\n' ...
                 '    this puts only two or three points on each one and the fitted centre and\n' ...
                 '    width will both be poor. 2 MHz is the default for a reason.\n'], step);
    end
    if any(abs(diff(diff(freqs))) > 1e-6)
        fprintf(['\nNote: the sweep is not evenly spaced. calibFieldSurvey handles that, but its\n' ...
                 'baseline polynomial is fitted in frequency, so a big gap in the middle can let\n' ...
                 'the baseline wander through it. Even spacing is safer.\n']);
    end
end

function [rho, theta, phi] = readSpherical(h)
% Read Brho/Btheta/Bphi without assuming which of the several accessors this class version
% has. Bsphere is the documented one; norm(h.B) is the fallback and cannot be wrong.
    rho = NaN; theta = NaN; phi = NaN;
    try
        s = h.Bsphere;
        if numel(s) >= 3; rho = s(1); theta = s(2); phi = s(3); end
    catch
    end
    if ~isfinite(rho)
        try; rho = norm(h.B); catch; end
    end
    if ~isfinite(theta)
        try; theta = h.Btheta; phi = h.Bphi; catch; end
    end
end

function restoreRho(h, rho0)
    try; h.setBrho(rho0); catch; end
end

function printIndented(txt)
    if isempty(txt); return; end
    L = strsplit(txt, newline);
    for i = 1:numel(L)
        if ~isempty(strtrim(L{i})); fprintf('      | %s\n', L{i}); end
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
