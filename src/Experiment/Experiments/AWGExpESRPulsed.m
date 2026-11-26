classdef AWGExpESRPulsed < ExperimentAWG
    %AWGEXPESRPULSED Pulsed ESR (pulsed ODMR): polarise, wait in the dark, ONE microwave
    %pulse at frequency f, read out. Swept over f.
    %
    %   p = AWGExpESRPulsed;
    %   p.loadJsonParams('all', jsonLocation);      % piTime, amplitude, windows, delays
    %   p.frequency = 2870 + (-20:0.5:20);
    %   p.useSequenceWalk = true;                   % the fast path (see below)
    %   p.newOrderEveryAverage = true;              % walk: new random order every average
    %   p.run;
    %
    % To SEE the sequence (no hardware touched):  p.plotPulsedSequence
    %
    % To check it against CW ESR, and to find its best settings, run
    %   out = esrPulsedVsCwCompare;     (Setup 7\experiments\_diagnostics)
    %
    % ===================================================================
    % WHY PULSED, AND WHAT IT BUYS
    % ===================================================================
    %
    % In CW ESR the laser and the microwave are on together, so the laser keeps
    % re-polarising while the MW mixes: the line is power-broadened by BOTH, and you
    % trade contrast against width with no way to have both. Pulsed ESR separates them
    % in time. The laser polarises into ms = 0 and is switched off; one MW pulse acts on
    % the spin in the dark; the laser comes back on and the first part of the readout
    % shows how much population the pulse moved to ms = +-1.
    %
    % With the pulse a pi pulse (duration tau, amplitude matched so the rotation is 180
    % degrees on resonance), the line shape is the pulse's own spectral response:
    %
    %       FWHM ~ 0.8 / tau    (MHz, tau in us)      -- plus 1/(pi T2*) inhomogeneous width
    %
    % So a LONGER, WEAKER pi pulse gives a NARROWER line at full contrast. 55 ns gives
    % ~15 MHz; to resolve the 14N hyperfine triplet (2.16 MHz) you want tau ~ 0.5-1 us,
    % i.e. about 20-25 dB less MW power than a 55 ns pi pulse. Beyond tau ~ T2* the line
    % stops narrowing and the contrast falls. amplitudeForPulse() gives the amplitude.
    %
    % ===================================================================
    % ONE REPEAT  (all times in us, as seen at the sample -- the pulse generator applies
    % the level-1 channel delays, so fixDelays must stay true)
    % ===================================================================
    %
    %   SIGNAL HALF
    %   |-- darkBeforeMW --|-- pulse --|-- lastDelay --|------- laserInitializationDuration -------|
    %    laser off: the     the AWG      dark gap        readout laser, which is also the
    %    metastable singlet  pulse,      before the      polarisation for the next half:
    %    empties (~0.2-0.3   mwPulse-    readout         |detectionStartDelay|detectionDuration|...
    %    us lifetime)        Duration                                        signal window
    %
    %   REFERENCE HALF, referenceMode = 'noMW' (default): the SAME half with no trigger
    %   and no MW switch. Same optical history, same durations, differs only by the MW.
    %   This is the clean reference: anything that is not the microwave cancels.
    %
    %   referenceMode = 'laserTail': no second half. The reference window sits near the
    %   END of the readout laser (referenceStartDelay after it comes on), where the spin
    %   is re-polarised. Half the repeat time, but the two windows see different parts of
    %   the laser pulse, so slow optical transients enter the ratio.
    %
    %   MW switch: opens mwSwitchLead before the pulse and closes mwSwitchTail after it.
    %   AWG trigger: rises at the pulse start (TRIG_DURATION long; its length does not
    %   matter, the AWG waveform sets the pulse length).
    %
    % Every half has the same duration at every frequency, so the duty cycle and any
    % heating are constant across the sweep, and nothing is changed per point but the
    % AWG sequence.
    %
    % ===================================================================
    % WHAT IS NOT INHERITED FROM AWGExpESR, AND WHY
    % ===================================================================
    %
    % This class derives from ExperimentAWG, as AWGExpESR does, not from AWGExpESR. Two
    % reasons. (1) MATLAB cannot redefine a Constant in a subclass, so a child of
    % AWGExpESR would be forced to be called 'AWG ESR' as well: constructing it would
    % delete any CW ESR object of that name, and its files would be saved under the CW
    % name. (2) Much of AWGExpESR is jump/FGEN/burst machinery for CW, none of which
    % applies to a pulse -- inherited, those switches would exist here and silently do
    % the wrong thing.
    %
    % What IS copied from AWGExpESR (2026-11-08), with the same names and meanings:
    %   useSequenceWalk        one sequence holds every frequency, armed ONCE per pg.run
    %   newOrderEveryAverage   walk: a new random visit order before every average
    %   passesPerRun, walkWarmupTriggers, resyncEveryPass, randomizeFrequencyOrder
    %   the load cache that also says WHY it missed, and the AWG/FGEN mode check.
    % passesPerRun = [] picks itself from the digitizer depth and a run-length cap, as
    % AWGExpRabi/Ramsey/Echo do (CW ESR still takes it by hand).
    %
    % NOT copied, and why: useSequenceJump (scrambled frequencies on this AWG), useFgen
    % (no sequencer, no pulses), mirrorSweepAround and the short looped burst (CW-only:
    % here the waveform IS one pulse, so the frequency is exact and nothing loops).
    %
    % Classic path (useSequenceWalk = false): every point re-arms unconditionally after
    % the sequence switch (Run waits for *OPC? first). AWGExpESR asks the AWG whether it
    % stopped and re-arms only then; that answer can come back before the switch has
    % taken effect. The switch stops the AWG anyway, so Run is needed either way and
    % costs the same. It does not depend on the magnetic field.
    %
    % Must live in prod\Experiment\Experiments, like every other experiment: SaveLoad only
    % listens to experiments found there, so anywhere else it would not autosave.

    properties (Constant)
        NAME = 'AWG ESR Pulsed'
    end

    properties
        % All of the parameters have a set function that raises changeFlag.
        mwPulseDuration = [];       % us. The MW pulse. [] = piTime (from the json).
        darkBeforeMW = 1;           % us. Laser off -> MW pulse start. Lets the metastable
                                    % singlet (lifetime ~0.2-0.3 us) empty back into the
                                    % triplet; shorter costs contrast.
        detectionStartDelay = 0;    % us. Readout laser on -> signal window opens.
        referenceMode = 'noMW';     % 'noMW' (default) | 'laserTail'. See the header.
        referenceStartDelay = [];   % us. laserTail only: readout laser on -> reference
                                    % window opens. [] = as late as the laser allows.
        mwSwitchLead = 0.05;        % us. MW switch opens this long before the pulse.
        mwSwitchTail = 0.05;        % us. ... and closes this long after it.
        randomizeFrequencyOrder = true; % Random visit order. Keep it on for real data: in
                                    % ascending order any slow drift of the light level
                                    % becomes a slope in frequency.
                                    %   classic: a NEW order every average.
                                    %   walk: ONE order per upload, reused by every pass and
                                    %     average, unless newOrderEveryAverage is on.
        useSequenceWalk = true;     % logical. The fast path (default), as in AWGExpESR. ONE sequence
                                    % per channel holds every frequency as a block
                                    % [zero pad (waits for trigger), pulse], the AWG is armed
                                    % once, and each trigger plays the next block. One pg.run
                                    % is then whole passes over the sweep, so the per-point
                                    % re-arm (~0.6 s) is gone. 'repeats' counts PASSES here.
                                    % With repeats = 50 on 201 points the classic path spends
                                    % ~2 min re-arming per average for ~1 s of data; walk
                                    % does the same average in about one pg.run.
                                    % Check it once against classic: esrPulsedVsCwCompare
                                    % does that (validateWalk).
        newOrderEveryAverage = true;  % logical. WALK ONLY. Draw a new random order before
                                    % every average and re-upload ONLY the sequence (the
                                    % waveforms stay). Costs one sequence upload + arm per
                                    % average, printed. Effects tied to position in the pass
                                    % then average down instead of sitting on the same
                                    % frequencies all run. Needs randomizeFrequencyOrder.
        passesPerRun = [];          % walk: passes over the sweep per pg.run. [] = choose:
                                    %   min(repeats, digitizer ceiling, run-length ceiling).
                                    % Each pg.run costs one re-arm, so bigger is faster.
        maxDigitizerSegments = 32000; % walk: ceiling on 2*(nFreq+warmup)*passesPerRun for
                                    % one pg.run (the value AWGExpRamsey/Echo measured, with
                                    % 1.5 us windows). A digitizer segment is as long as the
                                    % window, so for longer windows the ceiling used is
                                    % scaled down by 1.5/detectionDuration.
        maxRunSeconds = 10;         % walk: ceiling on one pg.run's length, so the plot
                                    % updates and the emergency stop responds.
        walkWarmupTriggers = 1;     % walk: throwaway blocks at the START of each pass, so
                                    % the first kept block is not the trigger after the
                                    % arm/idle gap (the classic path drops repeat 1 at every
                                    % point for the same reason).
        resyncEveryPass = true;     % walk: stop and re-arm before every pg.run, so a pass
                                    % always starts on block 1 and a retried pass cannot
                                    % shift the frequency mapping of the rest.
        channelPhase = [0 0];       % DEGREES, per AWG channel, written into the waveform.
                                    % Set timeDelay = 0 when using it (they add).
        driveMWSwitches = true;     % false = never open the MW switches (diagnostic).
        minDetectionGap = 5;        % us. Warn when two detection windows are closer than
                                    % this: the digitizer needs time to re-arm between
                                    % segments (the lab's experiments leave 5-10 us).
        piTimeAmplitude = [];       % dBm (channel 1) at which piTime is a pi pulse.
                                    % loadJsonParams sets it to the json's amplitude, so
                                    % amplitudeForPulse and the rotation note stay right
                                    % after you change amplitude. [] = current amplitude.
        % Inherited and used:
        %   piTime                       us, the calibrated pi pulse at 'amplitude'
        %   amplitude                    dBm (amplitudeDiff applied per channel)
        %   lastDelay                    us, MW pulse END -> readout laser on
        %   detectionDuration            us, both windows (the digitizer needs them equal)
        %   laserInitializationDuration  us, the whole readout laser per half
        %   frequency                    MHz, the sweep
    end

    properties (SetAccess = protected)
        pulseDurationActual = [];   % us. mwPulseDuration rounded to whole AWG samples.
        lastVisitOrder = [];        % frequency indices in the order the last average ran.
        timing = struct();          % the repeat as built, for the record and the printout.
        passesPerRunUsed = [];      % walk: what the automatic choice settled on.
        seqAssetIndex = [];         % walk: the ONE sequence per channel (cached).
        walkOrder = [];             % walk: walkOrder(k) = frequency index of block k.
        walkParts = [];             % walk: what a re-order needs (waveform indices).
        walkOrderFresh = false;     % walk: loaded order not yet used by an average.
    end

    %% Construction and setters
    methods
        function obj = AWGExpESRPulsed()
            obj@ExperimentAWG(AWGExpESRPulsed.NAME);
            obj.parameterName = 'frequencies';

            obj.repeats = 10000;
            obj.averages = 10;
            obj.frequency = 2870 + (-20:0.5:20);     % MHz
            obj.amplitude = -16 + [0 0];             % dBm
            obj.piTime = 0.055;                      % us -- the json overrides both
            obj.lastDelay = 0.5;                     % us
            obj.detectionDuration = 2;               % us
            obj.referenceDetectionDuration = 2;
            obj.laserInitializationDuration = 30;    % us

            obj.isPlotAlternateAvailable = true;
            obj.displayType2 = 'Raw signal (kcps)';   % the window shows 'Raw signal (V)' with the photodiode (Spcm.localiseUnits)
            obj.mCurrentXAxisParam = ExpParamDoubleVector('Frequency', [], [], 'MHz', obj.NAME);
            % Labels as in the robust fit's figure (ExpPlotStyle / fitStyle)
            obj.plotXLabel = ExpPlotStyle.X_FREQUENCY;
            obj.plotYLabel = ExpPlotStyle.Y_NORMALISED;
            obj.signalParam = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME);
        end

        function set.mwPulseDuration(obj, v)
            if ~isempty(v); obj.checkTimeScalar(v); end
            obj.mwPulseDuration = v;  obj.changeFlag = true;
        end
        function set.darkBeforeMW(obj, v)
            obj.checkNonNegScalar(v, 'darkBeforeMW');
            obj.darkBeforeMW = v;  obj.changeFlag = true;
        end
        function set.detectionStartDelay(obj, v)
            obj.checkNonNegScalar(v, 'detectionStartDelay');
            obj.detectionStartDelay = v;  obj.changeFlag = true;
        end
        function set.referenceMode(obj, v)
            if ~(ischar(v) || isstring(v)) || ~any(strcmpi(v, {'noMW', 'laserTail'}))
                obj.sendError('referenceMode must be ''noMW'' or ''laserTail''');
            end
            if strcmpi(v, 'noMW'); v = 'noMW'; else; v = 'laserTail'; end
            obj.referenceMode = v;  obj.changeFlag = true;
        end
        function set.referenceStartDelay(obj, v)
            if ~isempty(v); obj.checkNonNegScalar(v, 'referenceStartDelay'); end
            obj.referenceStartDelay = v;  obj.changeFlag = true;
        end
        function set.mwSwitchLead(obj, v)
            obj.checkNonNegScalar(v, 'mwSwitchLead');
            obj.mwSwitchLead = v;  obj.changeFlag = true;
        end
        function set.mwSwitchTail(obj, v)
            obj.checkNonNegScalar(v, 'mwSwitchTail');
            obj.mwSwitchTail = v;  obj.changeFlag = true;
        end
        function set.randomizeFrequencyOrder(obj, v)
            if ~ValidationHelper.isTrueOrFalse(v)
                obj.sendError('randomizeFrequencyOrder must be true or false');
            end
            obj.randomizeFrequencyOrder = logical(v);
        end
        function set.channelPhase(obj, v)
            if ~isnumeric(v) || ~isvector(v)
                obj.sendError('channelPhase must be a numeric vector, in degrees');
            end
            obj.channelPhase = mod(v(:)', 360);  obj.changeFlag = true;
        end
        function set.driveMWSwitches(obj, v)
            if ~ValidationHelper.isTrueOrFalse(v)
                obj.sendError('driveMWSwitches must be true or false');
            end
            obj.driveMWSwitches = logical(v);  obj.changeFlag = true;
        end
        function set.minDetectionGap(obj, v)
            obj.checkNonNegScalar(v, 'minDetectionGap');
            obj.minDetectionGap = v;
        end
        function set.useSequenceWalk(obj, v)
            if ~ValidationHelper.isTrueOrFalse(v)
                obj.sendError('useSequenceWalk must be true or false');
            end
            obj.useSequenceWalk = logical(v);  obj.changeFlag = true;
        end
        function set.newOrderEveryAverage(obj, v)
            if ~ValidationHelper.isTrueOrFalse(v)
                obj.sendError('newOrderEveryAverage must be true or false');
            end
            obj.newOrderEveryAverage = logical(v);
        end
        function set.passesPerRun(obj, v)
            if ~isempty(v) && ~ValidationHelper.isValuePositiveInteger(v)
                obj.sendError('passesPerRun must be a positive integer, or [] to choose it automatically');
            end
            obj.passesPerRun = v;
        end
        function set.maxDigitizerSegments(obj, v)
            if ~ValidationHelper.isValuePositiveInteger(v)
                obj.sendError('maxDigitizerSegments must be a positive integer');
            end
            obj.maxDigitizerSegments = v;
        end
        function set.maxRunSeconds(obj, v)
            if ~(isnumeric(v) && isscalar(v) && isfinite(v) && v > 0)
                obj.sendError('maxRunSeconds must be a positive number of seconds');
            end
            obj.maxRunSeconds = v;
        end
        function set.walkWarmupTriggers(obj, v)
            if ~(isnumeric(v) && isscalar(v) && v >= 0 && v == round(v))
                obj.sendError('walkWarmupTriggers must be a non-negative integer');
            end
            obj.walkWarmupTriggers = v;  obj.changeFlag = true;
        end
        function set.resyncEveryPass(obj, v)
            if ~ValidationHelper.isTrueOrFalse(v)
                obj.sendError('resyncEveryPass must be true or false');
            end
            obj.resyncEveryPass = logical(v);
        end
    end

    methods
        function n = getTotalNumberOfParams(obj)
            n = length(obj.frequency);
        end

        function loadJsonParams(obj, params, jsonLocation)
            % As Experiment.loadJsonParams, and then remember which amplitude the json's
            % piTime belongs to (piTimeAmplitude), so changing amplitude later does not
            % silently change what "pi" means.
            if ~exist('params', 'var'); params = 'all'; end
            if ~exist('jsonLocation', 'var'); jsonLocation = ''; end
            loadJsonParams@ExperimentAWG(obj, params, jsonLocation);
            p = params; if ~iscell(p); p = {p}; end
            if any(strcmpi(p, 'all')) || any(strcmp(p, 'piTime'))
                obj.piTimeAmplitude = obj.amplitude(1);
            end
        end

        function dBm = amplitudeForPulse(obj, tau)
            % Amplitude (dBm) that makes a pulse of length tau (us) a pi pulse, scaled
            % from the calibrated (piTime, piTimeAmplitude): the Rabi frequency goes as the
            % field amplitude, so tau_pi ~ 10^(-dBm/20). Assumes the MW chain is linear
            % over the range -- check with a Rabi at the new amplitude before trusting it.
            if isempty(obj.piTime)
                obj.sendError('piTime is empty: load the json or set it first');
            end
            dBm = obj.refAmplitude() + 20*log10(obj.piTime / tau);
        end

        function deg = rotationDegrees(obj, tau)
            % On-resonance rotation (degrees) of a pulse of tau us at the current amplitude.
            if ~exist('tau', 'var') || isempty(tau); tau = obj.pulseLength(); end
            deg = 180 * (tau / obj.piTime) * 10^((obj.amplitude(1) - obj.refAmplitude()) / 20);
        end

        function w = expectedLinewidth(obj)
            % FWHM (MHz) of a pi-pulse line from the pulse alone, 0.8/tau. The real line
            % is at least this wide; T2* adds to it.
            w = 0.8 / obj.pulseLength();
        end

        function fig = plotPulsedSequence(obj, fig)
            % One repeat as the NVs see it, drawn from the SAME timeline the experiment runs
            % (buildTimeline), without touching hardware: nothing is prepared or uploaded,
            % so it is safe while another experiment holds the AWG.
            %   p.plotPulsedSequence          new figure
            %   p.plotPulsedSequence(fig)     redraw into this figure
            % Panels
            %   top      the whole repeat, every pulse-generator channel, both halves and
            %            both detection windows, durations written on the bars
            %   middle   zoom on the MW pulse: dark, switch gate, trigger, the RF waveform
            %            the AWG plays (per channel, with its phase), lastDelay, laser on,
            %            start of the readout window -- in ns from the pulse start
            %   bottom   the pulse's spectral response with the sweep points on it (does the
            %            step resolve the line?), and the order the frequencies are visited
            %   right    the numbers: rotation angle, expected width, duty cycle, window gap,
            %            time per average walk vs classic, warnings
            % Times are AT THE SAMPLE (fixDelays shifts each channel by its own delay).
            % The inherited plotSequence still exists; it runs prepare (and so uploads to
            % the AWG) and draws the pulse generator's channels only.
            if ~exist('fig', 'var') || isempty(fig)
                fig = figure('Color', 'w', 'Position', [70 50 1300 820], ...
                             'Name', 'AWG ESR Pulsed - sequence');
            else
                clf(fig); set(fig, 'Color', 'w');
            end

            % ---- the repeat, exactly as prepare builds it ----
            S_rate = obj.AWG.SampleRate;                     % samples per us
            tReq = obj.pulseLength();
            tp = round(tReq * S_rate) / S_rate;
            P = obj.timelineParams();
            P.pulse = tp;
            [iv, T] = AWGExpESRPulsed.buildTimeline(P);
            rot = 180;
            if ~isempty(obj.piTime); rot = obj.rotationDegrees(tp); end

            % rows, top to bottom: triggers, MW switches, detector, laser
            present = {};
            for k = 1:numel(iv); present = [present, iv(k).ch]; end %#ok<AGROW>
            want = [{'greenLaser', 'detector'}, P.gate, P.trig];
            rows = want(ismember(want, present));
            rows = [rows, setdiff(unique(present), rows, 'stable')];
            nR = numel(rows);
            C = AWGExpESRPulsed.seqColors();
            rowName = cell(1, nR); rowCol = zeros(nR, 3);
            for r = 1:nR
                [rowName{r}, rowCol(r, :)] = AWGExpESRPulsed.seqRowStyle(rows{r}, C);
            end

            % ================= top: the whole repeat =================
            ax1 = axes('Parent', fig, 'Position', [0.07 0.67 0.90 0.26]);
            hold(ax1, 'on');
            span = T.repeatDuration;
            yTop = nR + 0.75;
            if strcmp(obj.referenceMode, 'noMW')
                patch(ax1, [T.halfDuration span span T.halfDuration], [0.4 0.4 yTop yTop], ...
                      C.refShade, 'EdgeColor', 'none');
                text(ax1, T.halfDuration/2, yTop + 0.3, 'SIGNAL half: dark, MW pulse, readout', ...
                     'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9);
                text(ax1, 1.5*T.halfDuration, yTop + 0.3, 'REFERENCE half: identical, no MW', ...
                     'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9, ...
                     'Color', C.dim);
            else
                text(ax1, span/2, yTop + 0.3, 'signal and reference windows in the same laser pulse (laserTail)', ...
                     'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9);
            end
            minW = 0.004 * span;
            for k = 1:numel(iv)
                for c = 1:numel(iv(k).ch)
                    r = find(strcmp(rows, iv(k).ch{c}));
                    t = iv(k).t; d = iv(k).d;
                    if d >= minW
                        patch(ax1, [t t+d t+d t], r + [-0.32 -0.32 0.32 0.32], rowCol(r, :), ...
                              'EdgeColor', 'none');
                        if d >= 0.07 * span
                            text(ax1, t + d/2, r, AWGExpESRPulsed.seqDur(d), 'Color', 'w', ...
                                 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
                        end
                    else
                        % too short to see at this scale: a tick and a marker (see the zoom)
                        plot(ax1, [t t], r + [-0.32 0.32], '-', 'Color', rowCol(r, :), 'LineWidth', 2);
                        plot(ax1, t, r + 0.47, 'v', 'Color', rowCol(r, :), ...
                             'MarkerFaceColor', rowCol(r, :), 'MarkerSize', 5);
                    end
                end
            end
            % window labels on the detector row
            rD = find(strcmp(rows, 'detector'));
            if ~isempty(rD)
                text(ax1, T.sigWindow, rD + 0.55, 'signal', 'FontSize', 8, 'Color', C.det);
                text(ax1, T.refWindow, rD + 0.55, 'reference', 'FontSize', 8, 'Color', C.det);
            end
            set(ax1, 'YTick', 1:nR, 'YTickLabel', rowName, 'YLim', [0.4 yTop + 0.6], ...
                'XLim', [0 span], 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);
            xlabel(ax1, 'time in the repeat (\mus, at the sample)');
            title(ax1, sprintf(['AWG ESR Pulsed: one repeat = %.3f \\mus   (%s path, %d frequencies, ' ...
                '%d %s, %d averages)'], span, BooleanHelper.ifTrueElse(obj.useSequenceWalk, 'walk', 'classic'), ...
                numel(obj.frequency), obj.repeats, ...
                BooleanHelper.ifTrueElse(obj.useSequenceWalk, 'passes', 'repeats per point'), obj.averages), ...
                'FontSize', 11);

            % ================= middle: zoom on the pulse =================
            ax2 = axes('Parent', fig, 'Position', [0.07 0.35 0.56 0.24]);
            hold(ax2, 'on');
            before = min(obj.darkBeforeMW, max(0.12, 3*obj.mwSwitchLead));
            tA = T.pulseStart - before;
            tB = T.laserOn + obj.detectionStartDelay + min(obj.detectionDuration, ...
                 max(0.15, 1.5*(tp + obj.lastDelay)));
            ns = @(t) (t - T.pulseStart) * 1e3;
            nCh = obj.nChannels;
            rfRows = nR + (1:nCh);
            for k = 1:numel(iv)
                a = max(iv(k).t, tA); b = min(iv(k).t + iv(k).d, tB);
                if b <= a; continue; end
                for c = 1:numel(iv(k).ch)
                    r = find(strcmp(rows, iv(k).ch{c}));
                    patch(ax2, ns([a b b a]), r + [-0.32 -0.32 0.32 0.32], rowCol(r, :), 'EdgeColor', 'none');
                    if iv(k).t >= tA && iv(k).t + iv(k).d <= tB && iv(k).d*1e3 >= 0.06*(tB - tA)*1e3
                        text(ax2, ns(iv(k).t + iv(k).d/2), r, sprintf('%.0f ns', iv(k).d*1e3), ...
                             'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
                    end
                end
            end
            % the RF the AWG plays, per channel (phase and inter-channel delay included)
            fC = mean([min(obj.frequency), max(obj.frequency)]);
            cp = zeros(1, nCh);
            cp(1:min(numel(obj.channelPhase), nCh)) = obj.channelPhase(1:min(numel(obj.channelPhase), nCh));
            later = BooleanHelper.ifTrueElse(obj.timeDelay > 0, 1, 2);
            nS = min(60000, max(400, ceil(tp * fC * 14)));
            tt = linspace(0, tp, nS);
            rfName = cell(1, nCh);
            for c = 1:nCh
                sh = 0; if c == later && nCh == 2; sh = abs(obj.timeDelay); end
                x = (tt + sh) * 1e3;
                y = rfRows(c) + 0.36 * cos(2*pi*fC*tt - deg2rad(cp(c)));
                patch(ax2, [x(1) x(end) x(end) x(1)], rfRows(c) + [-0.36 -0.36 0.36 0.36], ...
                      C.rfShade(c, :), 'EdgeColor', 'none');
                plot(ax2, x, y, '-', 'Color', C.rf(c, :), 'LineWidth', 0.5);
                if nCh == 2
                    rfName{c} = sprintf('AWG ch%d  (%.0f deg)', obj.activeChannels(c), cp(c));
                else
                    rfName{c} = sprintf('AWG ch%d', obj.activeChannels(c));
                end
            end
            yTop2 = rfRows(end) + 0.75;
            text(ax2, mean(ns([tA tB])), yTop2, sprintf(['MW pulse %.1f ns = %d samples, %.1f dBm -> %.0f deg ' ...
                 'on resonance  (drawn at %.1f MHz)'], tp*1e3, round(tp*S_rate), obj.amplitude(1), rot, fC), ...
                 'HorizontalAlignment', 'center', 'FontSize', 8.5, 'FontWeight', 'bold', 'Color', C.rf(1, :));
            % the gaps that are physics: dark (singlet empties), lastDelay (pulse -> laser).
            % Written between the laser and detector rows, left of laser-on, where it is empty.
            yG = 1.5;
            text(ax2, ns(tA), yG, sprintf(' \\leftarrow dark %.2f \\mus before the pulse (laser off)', ...
                 obj.darkBeforeMW), 'FontSize', 8, 'Color', C.dim);
            text(ax2, ns(T.laserOn), yG, sprintf('lastDelay %.0f ns \\rightarrow ', obj.lastDelay*1e3), ...
                 'FontSize', 8, 'Color', C.dim, 'HorizontalAlignment', 'right');
            plot(ax2, ns(T.laserOn)*[1 1], [0.4 yTop2 - 0.25], ':', 'Color', C.laser);
            plot(ax2, ns(T.sigWindow)*[1 1], [0.4 yTop2 - 0.25], ':', 'Color', C.det);
            set(ax2, 'YTick', 1:(nR + nCh), 'YTickLabel', [rowName, rfName], 'YLim', [0.3 yTop2 + 0.45], ...
                'XLim', ns([tA tB]), 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);
            xlabel(ax2, 'time from the MW pulse start (ns)');
            title(ax2, 'zoom on the pulse', 'FontSize', 10);

            % ================= bottom left: the pulse's spectral response =================
            ax3 = axes('Parent', fig, 'Position', [0.07 0.07 0.25 0.19]);
            hold(ax3, 'on');
            fr = obj.frequency(:)';
            Om = deg2rad(rot) / tp;                          % rad/us on resonance
            resp = @(dMHz) (Om^2 ./ (Om^2 + (2*pi*dMHz).^2)) .* sin(sqrt(Om^2 + (2*pi*dMHz).^2) * tp/2).^2;
            halfW = max([2.5*0.8/tp, 1.05*(max(fr) - min(fr))/2, 2]);
            dd = linspace(-halfW, halfW, 2000);
            pr = resp(dd);
            hm = find(pr >= max(pr)/2);
            fwhm = dd(hm(end)) - dd(hm(1));
            for hf = [-2.16 0 2.16]
                plot(ax3, [hf hf], [0 1.05], ':', 'Color', C.dim);
            end
            plot(ax3, dd, pr, '-', 'Color', C.rf(1, :), 'LineWidth', 1.5);
            plot(ax3, fr - fC, resp(fr - fC), '.', 'Color', [0.1 0.1 0.1], 'MarkerSize', 9);
            step = NaN; if numel(fr) > 1; step = median(diff(sort(fr))); end
            set(ax3, 'XLim', [-halfW halfW], 'YLim', [0 1.08], 'TickDir', 'out', 'Box', 'off', 'FontSize', 8);
            xlabel(ax3, 'detuning from sweep centre (MHz)');
            ylabel(ax3, 'flip probability');
            title(ax3, sprintf('pulse response: FWHM %.2f MHz, %.1f points per FWHM', fwhm, fwhm/step), ...
                  'FontSize', 9);
            text(ax3, 2.16, 1.04, ' 14N', 'FontSize', 7, 'Color', C.dim);

            % ================= bottom middle: visit order =================
            ax4 = axes('Parent', fig, 'Position', [0.39 0.07 0.24 0.19]);
            hold(ax4, 'on');
            N = numel(fr);
            if obj.useSequenceWalk && numel(obj.walkOrder) == N
                order = obj.walkOrder; src = 'as loaded';
            elseif numel(obj.lastVisitOrder) == N
                order = obj.lastVisitOrder; src = 'last average';
            elseif obj.randomizeFrequencyOrder
                order = randperm(N); src = 'example';
            else
                order = 1:N; src = 'ascending';
            end
            if obj.useSequenceWalk
                nW = obj.walkWarmupTriggers;
                blk = [repmat(order(1), 1, nW), order];
                plot(ax4, 1:numel(blk), fr(blk), '-', 'Color', [0.8 0.8 0.8]);
                plot(ax4, nW + (1:N), fr(order), '.', 'Color', C.rf(1, :), 'MarkerSize', 8);
                if nW > 0
                    plot(ax4, 1:nW, fr(blk(1:nW)), 'x', 'Color', C.dim, 'MarkerSize', 7);
                end
                xlabel(ax4, 'block in the walk sequence (x = warm-up, discarded)');
                ttl = sprintf('walk order (%s): %s', src, BooleanHelper.ifTrueElse( ...
                    obj.newOrderEveryAverage && obj.randomizeFrequencyOrder, 'new every average', ...
                    'same all run'));
            else
                plot(ax4, 1:N, fr(order), '.-', 'Color', C.rf(1, :), 'MarkerSize', 8);
                xlabel(ax4, 'point in the average');
                ttl = sprintf('classic order (%s): %s', src, BooleanHelper.ifTrueElse( ...
                    obj.randomizeFrequencyOrder, 'new every average', 'ascending'));
            end
            ylabel(ax4, 'MHz');
            set(ax4, 'TickDir', 'out', 'Box', 'off', 'FontSize', 8);
            title(ax4, ttl, 'FontSize', 9);

            % ================= right: the numbers =================
            ax5 = axes('Parent', fig, 'Position', [0.665 0.04 0.33 0.56], 'Visible', 'off');
            L = obj.sequenceSummary(tReq, tp, rot, fwhm, step, T, S_rate);
            text(ax5, 0, 1, sprintf('%s\n', L{:}), 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                 'FontName', 'Courier New', 'FontSize', 8.5, 'Interpreter', 'none');
        end

        function L = sequenceSummary(obj, tReq, tp, rot, fwhm, step, T, S_rate)
            % The text panel of plotPulsedSequence, as a cell of lines (also handy to print).
            fr = obj.frequency; N = numel(fr);
            L = {};
            L{end+1} = 'PULSE';
            L{end+1} = sprintf('  %.2f ns asked -> %.2f ns played (%d samples, %.0f GS/s)', ...
                               tReq*1e3, tp*1e3, round(tp*S_rate), S_rate/1e3);
            if ~isempty(obj.piTime)
                L{end+1} = sprintf('  %.1f dBm; pi = %.2f ns at %.1f dBm (json)', obj.amplitude(1), ...
                                   obj.piTime*1e3, obj.refAmplitude());
            end
            L{end+1} = sprintf('  rotation on resonance: %.0f deg', rot);
            L{end+1} = sprintf('  line from the pulse alone: FWHM %.2f MHz', fwhm);
            L{end+1} = '';
            L{end+1} = 'SWEEP';
            if N > 1
                L{end+1} = sprintf('  %d points, %.2f - %.2f MHz, step %.3f MHz', N, min(fr), max(fr), step);
                L{end+1} = sprintf('  %.1f points per pulse FWHM', fwhm/step);
            else
                L{end+1} = sprintf('  1 point at %.3f MHz', fr);
            end
            L{end+1} = '';
            L{end+1} = 'REPEAT';
            nHalves = BooleanHelper.ifTrueElse(strcmp(obj.referenceMode, 'noMW'), 2, 1);
            laserFrac = nHalves * obj.laserInitializationDuration / T.repeatDuration;
            winFrac = 2 * obj.detectionDuration / T.repeatDuration;
            L{end+1} = sprintf('  %.3f us = %s', T.repeatDuration, ...
                BooleanHelper.ifTrueElse(nHalves == 2, '2 halves (signal, reference)', '1 half (laserTail)'));
            L{end+1} = sprintf('  laser on %.0f%% of the time, windows %.1f%%', 100*laserFrac, 100*winFrac);
            L{end+1} = sprintf('  windows 2 x %.2f us, closest gap %.2f us', obj.detectionDuration, T.minWindowGap);
            L{end+1} = sprintf('  reference: %s', obj.referenceMode);
            L{end+1} = '';
            L{end+1} = 'TIME PER AVERAGE (estimate)';
            tData = N * obj.repeats * T.repeatDuration * 1e-6;
            L{end+1} = sprintf('  sequence time      %.1f s', tData);
            L{end+1} = sprintf('  classic: + %d re-arms (~0.6 s) = ~%.1f s', N, tData + 0.6*N);
            Pw = obj.choosePassesPerRun(T.repeatDuration, false);
            nRun = ceil(obj.repeats / Pw);
            L{end+1} = sprintf('  walk:    + %d pg.runs (P = %d) = ~%.1f s', nRun, Pw, ...
                               tData * (N + obj.walkWarmupTriggers) / max(N, 1) + 0.3*nRun);
            L{end+1} = sprintf('  -> walk is %s', BooleanHelper.ifTrueElse(obj.useSequenceWalk, 'ON', ...
                               'OFF (p.useSequenceWalk = true)'));
            L{end+1} = '';
            L{end+1} = 'CHANNELS';
            L{end+1} = sprintf('  nChannels %d, AWG outputs %s, phase %s deg', obj.nChannels, ...
                               mat2str(obj.activeChannels), mat2str(round(obj.channelPhase)));
            L{end+1} = sprintf('  timeDelay %.4g us, MW switches %s', obj.timeDelay, ...
                               BooleanHelper.ifTrueElse(obj.driveMWSwitches, 'driven', 'NOT driven'));
            % warnings
            W = {};
            if abs(rot/180 - 1) > 0.05
                W{end+1} = sprintf('  not a pi pulse: amplitude for pi = %.1f dBm', obj.amplitudeForPulse(tp));
            end
            if N > 1 && step > fwhm/3
                W{end+1} = sprintf('  step %.2f MHz > FWHM/3: line under-sampled', step);
            end
            if obj.darkBeforeMW < 0.5
                W{end+1} = sprintf('  darkBeforeMW %.2f us < 0.5: singlet may not be empty', obj.darkBeforeMW);
            end
            if T.minWindowGap < obj.minDetectionGap
                W{end+1} = sprintf('  window gap %.2f us < %g us (digitizer re-arm)', T.minWindowGap, obj.minDetectionGap);
            end
            if round(tp*S_rate) < 8
                W{end+1} = '  pulse under 8 samples: its edges are a big part of it';
            end
            if ~isempty(W)
                L{end+1} = '';
                L{end+1} = 'WARNINGS';
                L = [L, W];
            end
        end
    end

    %% Overridden from Experiment
    methods (Access = protected)
        function prepare(obj)
            obj.detectionPeriodsPerRepeat = 2;      % signal, reference
            obj.runsPerPerform = 1;
            % The digitizer path needs the two windows equal; they are, by construction.
            obj.referenceDetectionDuration = obj.detectionDuration;

            % The pulse, in whole AWG samples -- that is the pulse the NVs get.
            S_rate = obj.AWG.SampleRate;                            % samples per us
            obj.pulseDurationActual = round(obj.pulseLength() * S_rate) / S_rate;

            [iv, T] = AWGExpESRPulsed.buildTimeline(obj.timelineParams());
            obj.timing = T;

            S = Sequence(Pulse(T.repeatDuration, '', ''));
            for k = 1:numel(iv)
                S.addEventAtGivenTime(iv(k).t, iv(k).d, iv(k).ch);
            end
            S.repairSequence();

            % Re-enable the outputs on every run (wrapUp turns them off). activeChannels,
            % not 1:nChannels, so nChannels = 1 with MWChannel = 'MW2' drives output 2.
            obj.AWG.workChannels = obj.activeChannels;
            obj.prepareInternal(S);
            obj.LoadAWG();          % cheap when the AWG already holds this set

            if obj.useSequenceWalk
                % prepareInternal sized the pulse generator and digitizer for 'repeats'
                % triggers at ONE frequency. Walk wants one trigger per BLOCK, whole passes
                % per pg.run. Override those two numbers (as AWGExpESR/AWGExpRabi do).
                pg   = getObjByName(PulseGenerator.NAME);
                spcm = getObjByName(Spcm.NAME);
                nF = length(obj.frequency);
                blocksPerPass = nF + obj.walkWarmupTriggers;
                obj.passesPerRunUsed = obj.choosePassesPerRun(T.repeatDuration, true);
                triggersPerRun = blocksPerPass * obj.passesPerRunUsed;
                pg.repeats = triggersPerRun;
                numScans = obj.detectionPeriodsPerRepeat * triggersPerRun;
                seqTime = pg.sequenceDuration() * 1e-6;
                spcm.prepareExperimentCount(numScans, 2 * seqTime * triggersPerRun, 1);
                fprintf(['Walk: %d triggers per pg.run (%d frequencies + %d warm-up, x %d ' ...
                         'passes), %d digitizer segments, %.2f s of data per pg.run.\n'], ...
                        triggersPerRun, nF, obj.walkWarmupTriggers, obj.passesPerRunUsed, ...
                        numScans, triggersPerRun * T.repeatDuration * 1e-6);
            end

            obj.reportTiming();
            obj.mCurrentXAxisParam.value = obj.frequency;
        end

        function P = choosePassesPerRun(obj, repeatUs, announce)
            % Passes per pg.run, derived (same rule as AWGExpEcho):
            %   digitizer:  2 * (nFreq + warmup) * P  <=  maxDigitizerSegments (scaled, above)
            %   clock:      (nFreq + warmup) * P * repeat  <=  maxRunSeconds
            % capped at repeats (there is nothing more to fit) and floored at 1. Walk costs
            % one re-arm per pg.run, classic one per frequency, so P as big as allowed wins.
            if ~exist('announce', 'var'); announce = false; end
            blocksPerPass = length(obj.frequency) + obj.walkWarmupTriggers;
            segs = obj.maxDigitizerSegments * min(1, 1.5 / obj.detectionDuration);
            pDig = floor(segs / (2 * blocksPerPass));
            pClk = floor(obj.maxRunSeconds / max(blocksPerPass * repeatUs * 1e-6, eps));
            if ~isempty(obj.passesPerRun)
                P = obj.passesPerRun;
            else
                P = min(pDig, pClk);
            end
            P = max(1, min(P, obj.repeats));
            if announce
                fprintf(['Walk: passesPerRun = %d (%s; ceilings %d digitizer, %d clock). ' ...
                         're-arms per average: %d walk vs %d classic.\n'], P, ...
                        BooleanHelper.ifTrueElse(isempty(obj.passesPerRun), 'auto', 'set by hand'), ...
                        pDig, pClk, ceil(obj.repeats / P), length(obj.frequency));
            end
        end

        function LoadAWG(obj)
            % Classic: one sequence per (frequency, channel) = [zero pad that waits for the
            % trigger, pulse waveform (pulse, then zeros to the minimum length)] -- the
            % structure the level-1 trigger calibration and AWGExpRabi use, so the
            % calibrated trigger delay applies as it is.
            % Walk: the SAME waveforms, but ONE sequence per channel in which every
            % frequency is a block of those same two steps (see uploadWalkSequence).
            obj.ensureInstrumentMode('AWG');
            [fp, parts] = obj.awgFingerprint();
            if obj.tryLoadFromCache(fp, parts); return; end

            N   = length(obj.frequency);
            nCh = obj.nChannels;
            tp  = obj.pulseDurationActual;

            obj.AWG.ResetStoredData;
            obj.AWG.workChannels = obj.activeChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            A     = obj.AWG.MaxNormPowerAllowed;
            minWF = obj.AWG.minWaveformDuration;
            S_rate = obj.AWG.SampleRate;

            % Inter-channel delay: whole samples go into the leading zero pad of the
            % LATER channel, the sub-sample remainder into its phase (as
            % ExperimentAWG.channelsDifference does, per frequency here because the
            % remainder's phase depends on f). channelPhase is added on top.
            dAbs   = abs(obj.timeDelay);
            dWhole = floor(dAbs * S_rate) / S_rate;
            later  = BooleanHelper.ifTrueElse(obj.timeDelay > 0, 1, 2);
            prePad = minWF * ones(1, nCh);
            if dAbs > 0 && later <= nCh; prePad(later) = minWF + dWhole; end
            cp = zeros(1, nCh);
            cp(1:min(numel(obj.channelPhase), nCh)) = obj.channelPhase(1:min(numel(obj.channelPhase), nCh));
            phi0 = deg2rad(cp);

            obj.AWG.progressBegin(sprintf(['Pulsed ESR waveforms: %d frequencies x %d ' ...
                'channel(s), one %.1f ns pulse each (+%.3f us zero tail), plus %d zero pads.'], ...
                N, nCh, 1e3*tp, minWF, nCh), N*nCh + nCh, 'waveform');
            zeroIdx = zeros(1, nCh);
            for chan = 1:nCh
                zeroIdx(chan) = obj.AWG.AddWaveformByFunction(@(t) 0*t, prePad(chan), 0, 1);
                obj.AWG.progressStep();
            end
            wfIdx = zeros(nCh, N);
            for i = 1:N
                w = 2*pi*obj.frequency(i);
                for chan = 1:nCh
                    phi = phi0(chan);
                    if chan == later && dAbs > 0; phi = phi + w*(dAbs - dWhole); end
                    f = @(t) A(chan) * (0 <= t & t < tp) .* cos(w*t - phi);
                    wfIdx(chan, i) = obj.AWG.AddWaveformByFunction(f, tp + minWF, 0, 1);
                    obj.AWG.progressStep();
                end
            end
            obj.AWG.progressEnd();
            obj.AWG.LoadWaveform();

            if obj.useSequenceWalk
                % FIELD BY FIELD: struct(...) would expand a cell value into an array.
                p = struct();
                p.zeroIdx = zeroIdx;
                p.wfIdx   = wfIdx;
                p.nWarm   = obj.walkWarmupTriggers;
                obj.walkParts = p;
                tSeq = tic;
                nSteps = obj.uploadWalkSequence(obj.drawWalkOrder());
                obj.stampLoaded(fp, parts);
                fprintf(['Walk load: %d frequencies + %d warm-up block(s), %d steps per channel, ' ...
                         '%s order; sequence + first arm %.1f s.\n'], N, p.nWarm, nSteps, ...
                        BooleanHelper.ifTrueElse(isequal(obj.walkOrder, 1:N), 'ascending', 'random'), ...
                        toc(tSeq));
                return
            end

            obj.AWG.progressBegin('Pulsed ESR sequences: zero pad (waits for trigger) + pulse.', ...
                N*nCh, 'sequence');
            obj.indexSeq = zeros(nCh, N);
            for i = 1:N
                for chan = 1:nCh
                    obj.indexSeq(chan, i) = obj.AWG.AddSequence([zeroIdx(chan), wfIdx(chan, i)], ...
                                                                [1, 1], [1, 0]);
                    obj.AWG.progressStep();
                end
            end
            obj.AWG.progressEnd();
            obj.AWG.LoadSequence();
            obj.walkParts = [];  obj.walkOrder = [];  obj.seqAssetIndex = [];
            obj.stampLoaded(fp, parts);
        end

        function tf = tryLoadFromCache(obj, fp, parts)
            % True when the AWG already holds exactly this set (and it was re-armed if walk).
            tf = false;
            if ~obj.AWG.isAlreadyLoaded(fp)
                % Say WHY it missed: a hash alone cannot name the field that moved.
                cached = obj.AWG.getLoadedCache();
                if isstruct(cached) && isfield(cached, 'parts') && isstruct(cached.parts) && ...
                        isfield(cached.parts, 'cls') && strcmp(cached.parts.cls, parts.cls)
                    obj.reportCacheMiss(cached.parts, parts);
                elseif isempty(cached)
                    fprintf('AWG load cache: empty -- uploading.\n');
                else
                    fprintf('AWG holds another experiment''s set -- uploading.\n');
                end
                return
            end
            % Verify against the instrument rather than trusting the hash alone.
            if obj.useSequenceWalk
                nWanted = numel(obj.activeChannels);
            else
                nWanted = length(obj.frequency) * obj.nChannels;
            end
            nLoaded = str2double(obj.AWG.queryRaw('SLISt:SIZE?', 10));
            if nLoaded ~= nWanted
                obj.sendWarning(sprintf(['AWG holds %d sequences but this sweep needs %d -- ' ...
                    'the load cache was stale. Reloading.'], nLoaded, nWanted));
                obj.AWG.invalidateLoaded();
                return
            end
            cached = obj.AWG.getLoadedCache();
            obj.indexSeq = cached.indexSeq;
            % Amplitude and the channel enables live in the instrument front end, not in the
            % waveform data, so re-apply them rather than assume they survived.
            obj.AWG.workChannels = obj.activeChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            if obj.useSequenceWalk
                obj.seqAssetIndex = cached.seqAssetIndex;
                obj.walkOrder     = cached.walkOrder;
                obj.walkParts     = cached.walkParts;
                obj.walkOrderFresh = false;         % used by an earlier run already
                % Walk arms once and never touches the AWG during the sweep, so a cache hit
                % still has to assign and arm.
                ch = obj.activeChannels;
                for c = 1:numel(ch)
                    obj.AWG.assignSequence(obj.seqAssetIndex(c), ch(c));
                end
                obj.AWG.Run;
            end
            fprintf('AWG already holds this pulse set (%d sequences) - skipping upload.\n', nLoaded);
            tf = true;
        end

        function reportCacheMiss(obj, was, now)                             %#ok<INUSL>
            % Name the fields that differ, so a miss is one line instead of a mystery.
            f = fieldnames(now);
            changed = {};
            for i = 1:numel(f)
                a = []; b = now.(f{i});
                if isfield(was, f{i}); a = was.(f{i}); end
                if ~isequal(a, b)
                    if isnumeric(b) && isnumeric(a) && numel(a) == numel(b) && numel(b) > 3
                        changed{end+1} = sprintf('%s (%d values, max diff %.4g)', ...
                            f{i}, numel(b), max(abs(double(b(:)) - double(a(:)))));  %#ok<AGROW>
                    elseif isnumeric(b) && isscalar(b) && isnumeric(a) && isscalar(a)
                        changed{end+1} = sprintf('%s (%g -> %g)', f{i}, a, b);       %#ok<AGROW>
                    else
                        changed{end+1} = f{i};                                       %#ok<AGROW>
                    end
                end
            end
            if isempty(changed)
                fprintf('AWG load cache missed although every field matches -- uploading.\n');
            else
                fprintf('AWG load cache miss, because: %s\n', strjoin(changed, ', '));
            end
        end

        function stampLoaded(obj, fp, parts)
            % Call only after a load (or a walk re-order) has completed.
            cache = struct();
            cache.indexSeq      = obj.indexSeq;
            cache.seqAssetIndex = obj.seqAssetIndex;
            cache.walkOrder     = obj.walkOrder;
            cache.walkParts     = obj.walkParts;
            cache.parts         = parts;
            obj.AWG.markLoaded(fp, cache);
        end

        function ensureInstrumentMode(obj, wanted)
            % As in AWGExpESR: an ESR run with useFgen leaves the instrument in FGEN mode,
            % where an uploaded sequence cannot play. Switch only when needed (a switch
            % stops playback) and drop the load cache when it happens.
            try
                current = strtrim(obj.AWG.queryRaw('INST:MODE?'));
            catch
                return          % cannot ask: carry on as before
            end
            if strcmpi(current, wanted); return; end
            fprintf('AWG instrument mode: %s -> %s\n', current, wanted);
            obj.AWG.setValue('mode', wanted);
            obj.AWG.invalidateLoaded;
        end

        function order = drawWalkOrder(obj)
            N = length(obj.frequency);
            if obj.randomizeFrequencyOrder
                order = randperm(N);
            else
                order = 1:N;
            end
        end

        function nSteps = uploadWalkSequence(obj, order)
            % ONE sequence per channel playing the frequencies in 'order' (order(k) = the
            % frequency index of block k), uploaded, assigned and armed. Only the sequence
            % list is replaced: the waveforms in obj.walkParts must already be on the AWG.
            %
            % Block = the classic sequence's two steps:
            %   [zero pad, waits for trigger]  [pulse waveform]      trigger flags [1 0]
            % plus walkWarmupTriggers copies of the first block in front. No GOTO: the last
            % step falls back to the first by default, so a pass ends where the next begins.
            p = obj.walkParts;
            N = numel(order);
            blockFreq = [repmat(order(1), 1, p.nWarm), order];
            nBlocks = numel(blockFreq);
            nSteps  = 2 * nBlocks;
            if nSteps > obj.AWG.MAX_SEQUENCE_STEPS
                obj.sendError(sprintf(['%d frequencies need %d sequence steps, over the %d ' ...
                    'step limit. Split the sweep or use the classic path.'], N, nSteps, ...
                    obj.AWG.MAX_SEQUENCE_STEPS));
            end

            % Nothing may still reference the sequence that is about to be replaced.
            obj.AWG.sendCommand('awgcontrol:stop:immediate');
            obj.AWG.sequences = struct('waveforms', {{}}, 'repeats', {{}}, ...
                                       'trigger', {{}}, 'goto', {{}});   % sequences only
            ch = obj.activeChannels;
            obj.indexSeq = zeros(numel(ch), N);     % STEP numbers here, not sequences
            obj.seqAssetIndex = zeros(1, numel(ch));
            for c = 1:numel(ch)
                wfList = zeros(1, nSteps);
                for blk = 1:nBlocks
                    i  = blockFreq(blk);
                    s0 = 2 * (blk - 1);
                    wfList(s0 + (1:2)) = [p.zeroIdx(c), p.wfIdx(c, i)];
                    if blk > p.nWarm
                        obj.indexSeq(c, i) = s0 + 1;
                    end
                end
                obj.seqAssetIndex(c) = obj.AWG.AddSequence(wfList, ones(1, nSteps), ...
                    repmat([1 0], 1, nBlocks), 0, zeros(1, nSteps));
            end
            obj.AWG.LoadSequence();

            % Assign and arm ONCE for the whole sweep.
            for c = 1:numel(ch)
                obj.AWG.assignSequence(obj.seqAssetIndex(c), ch(c));
            end
            obj.AWG.Run;
            obj.walkOrder = order;
            obj.walkOrderFresh = true;      % not used by an average yet
        end

        function newWalkOrderIfAsked(obj)
            % newOrderEveryAverage: before every average except the first one after a fresh
            % load, draw a new order and re-upload ONLY the sequence. The load cache is
            % re-stamped with the new order, so a later cache hit maps its data correctly.
            if obj.newOrderEveryAverage && obj.randomizeFrequencyOrder && ~obj.walkOrderFresh
                if isempty(obj.walkParts)
                    obj.sendWarning(['newOrderEveryAverage: the loaded sequence cannot be ' ...
                        'rebuilt on its own. Keeping its order for this average; change any ' ...
                        'setting to force a reload.']);
                else
                    tNew = tic;
                    obj.uploadWalkSequence(obj.drawWalkOrder());
                    [fp, parts] = obj.awgFingerprint();
                    obj.stampLoaded(fp, parts);
                    fprintf('Average %d: new visit order -- sequence re-upload + arm %.1f s.\n', ...
                            obj.currIter, toc(tNew));
                end
            end
            obj.walkOrderFresh = false;
        end

        function perform(obj)
            if obj.useSequenceWalk
                obj.performWalk();
                return
            end
            pg = getObjByName(PulseGenerator.NAME);
            spcm = getObjByName(Spcm.NAME);
            tracker = getObjByName(Tracker.NAME);
            if isempty(tracker); throwBaseObjException(Tracker.NAME); end

            len = length(obj.frequency);
            if obj.randomizeFrequencyOrder
                order = randperm(len);
            else
                order = 1:len;
            end
            obj.lastVisitOrder = order;
            ch = obj.activeChannels;

            for k = 1:len
                i = order(k);
                if obj.currIter <= size(obj.signal, 3) && sum(obj.signal(:, i, obj.currIter)) ~= 0
                    continue        % already measured (resumed average)
                end
                success = false;
                for trial = 1:5
                    if obj.checkEmergencyStop(); return; end
                    try
                        for c = 1:numel(ch)
                            obj.AWG.assignSequence(obj.indexSeq(c, i), ch(c));
                        end
                        % Assigning a sequence STOPS the AWG and completes in the
                        % background. Run waits for it (*OPC?), arms, and waits until the
                        % instrument reports armed. Never skip this on the strength of a
                        % state query: the query can be answered before the switch has
                        % happened, and the point is then measured at the previous
                        % frequency.
                        obj.AWG.Run;

                        data = obj.getRawData(pg, spcm);
                        [sig, sterr] = obj.processData(data);
                        obj.signal(:, i, obj.currIter) = sig;
                        obj.sterr(:, i, obj.currIter) = sterr;

                        success = true;
                        obj.currParamIter = obj.currParamIter + 1;
                        if obj.currParamIter == 1 || mod(obj.currParamIter, 5) == 0 || ...
                                obj.currParamIter == obj.totalNumberOfParams
                            sendEventParamIterationDone(obj);
                        end

                        if obj.isTracking
                            if obj.checkEmergencyStop(); return; end
                            if tracker.compareReference(sig(2), sterr(2), ...
                                    Tracker.REFERENCE_TYPE_KCPS, obj.trackThreshhold)
                                tracker.trackUsing(TrackablePosition.NAME, obj);
                                obj.prepare;
                            end
                        end
                        break
                    catch err
                        err2warning(err);
                        fprintf('Experiment failed at trial %d, attempting again.\n', trial);
                        try; spcm.stopExperimentCount; catch; end
                    end
                end
                if ~success
                    obj.emergencyStop;
                    return
                end
            end
            obj.updateSignalParam();
        end

        function performWalk(obj)
            % One average, acquired as whole passes over the sweep. A pg.run of
            % (nFreq + warmup) * P triggers returns [s1 r1 s2 r2 ... ] in BLOCK order
            % (obj.walkOrder), P times over; it is filed back under frequency here and
            % reduced by walkReduce, a transcription of processData. Nothing is sent to
            % the AWG inside the loop except the per-pg.run re-arm (resyncEveryPass).
            nF = length(obj.frequency);
            pg   = getObjByName(PulseGenerator.NAME);
            spcm = getObjByName(Spcm.NAME);
            tracker = getObjByName(Tracker.NAME);
            if isempty(tracker); throwBaseObjException(Tracker.NAME); end

            obj.newWalkOrderIfAsked();              % newOrderEveryAverage
            if numel(obj.walkOrder) ~= nF; obj.walkOrder = 1:nF; end
            obj.lastVisitOrder = obj.walkOrder;

            passesWanted = obj.repeats;             % passes play the role of repeats
            P = obj.passesPerRunUsed;
            if isempty(P)
                P = obj.choosePassesPerRun(obj.timing.repeatDuration, false);
                obj.passesPerRunUsed = P;
            end
            nRuns  = ceil(passesWanted / P);
            nWarm  = obj.walkWarmupTriggers;
            blocksPerPass = nF + nWarm;
            nDet = obj.detectionPeriodsPerRepeat;
            readsPerRun = nDet * blocksPerPass * P;

            raw = nan(nDet, nF, nRuns * P);
            passesDone = 0;  failedPasses = 0;
            opts = obj.walkReduceOptions(spcm);

            for r = 1:nRuns
                if obj.checkEmergencyStop(); return; end
                success = false;
                for trial = 1:5
                    try
                        if obj.resyncEveryPass
                            % Back to step 1, so this pg.run is known to start on block 1.
                            obj.AWG.sendCommand('awgcontrol:stop:immediate');
                            obj.AWG.Run;
                        else
                            obj.AWG.EnsureRunning;
                        end
                        data = obj.getRawData(pg, spcm);
                        if numel(data) < readsPerRun
                            error('Walk pass returned %d reads, expected %d', numel(data), readsPerRun);
                        end
                        passData = reshape(data(1:readsPerRun), nDet, blocksPerPass, P);
                        % Drop the warm-up blocks; block k played frequency walkOrder(k).
                        raw(:, obj.walkOrder, passesDone + (1:P)) = passData(:, nWarm+1:end, :);
                        passesDone = passesDone + P;
                        success = true;
                        break
                    catch err
                        err2warning(err);
                        failedPasses = failedPasses + 1;
                        fprintf('Walk pass %d failed at trial %d, attempting again.\n', r, trial);
                        try; spcm.stopExperimentCount; catch; end
                    end
                end
                if ~success
                    obj.emergencyStop;
                    return
                end

                % Reduce what has arrived so far, so the plot fills in during the average.
                usable = min(passesDone, passesWanted);
                [sig, err] = walkReduce(raw(:, :, 1:usable), nDet, nF, opts);
                obj.signal(:, :, obj.currIter) = sig;
                obj.sterr(:, :, obj.currIter)  = err;
                obj.updateSignalParam();

                obj.currParamIter = min(nF, round(nF * usable / passesWanted));
                if r == 1 || mod(r, 5) == 0 || r == nRuns
                    sendEventParamIterationDone(obj);
                end

                % Tracking once per pg.run: the whole-sweep mean reference is the brightness.
                if obj.isTracking
                    if obj.checkEmergencyStop(); return; end
                    if tracker.compareReference(mean(sig(2, :)), mean(err(2, :)), ...
                            Tracker.REFERENCE_TYPE_KCPS, obj.trackThreshhold)
                        tracker.trackUsing(TrackablePosition.NAME, obj);
                        obj.prepare;        % re-arms the walk sequence (cache hit)
                    end
                end
            end

            if failedPasses == 0
                fprintf('Walk: %d pg.run(s), none retried.\n', nRuns);
            elseif obj.resyncEveryPass
                fprintf(['Walk: %d of %d pg.runs retried. resyncEveryPass re-armed each from ' ...
                         'block 1, so the frequency mapping is intact.\n'], failedPasses, nRuns);
            else
                obj.sendWarning(sprintf(['Walk: %d of %d pg.runs retried with resyncEveryPass OFF. ' ...
                    'An aborted pass leaves the sequencer on an unknown block, so later passes ' ...
                    'can be shifted. Turn resyncEveryPass on.'], failedPasses, nRuns));
            end
        end

        function opts = walkReduceOptions(obj, spcm)
            % What walkReduce needs to reproduce processData exactly (as in AWGExpESR).
            opts.hasPhotodiode = spcm.hasPhotodiode;
            opts.detectionWithGI = false;
            opts.GIgain = 1;
            if opts.hasPhotodiode && isprop(spcm, 'detectionWithGI')
                opts.detectionWithGI = logical(spcm.detectionWithGI);
                if opts.detectionWithGI && isprop(spcm, 'GIgain')
                    opts.GIgain = spcm.GIgain;
                end
            end
            opts.detectionDuration = obj.detectionDuration;
            opts.referenceDetectionDuration = obj.referenceDetectionDuration;
        end

        function wrapUp(obj)
            obj.wrapUpInternal();
        end

        function params = alternateSignal(obj)
            % Alternate view: the signal and reference windows themselves.
            if obj.currIter == 0
                params = ExpResultDoubleVector('fluorescence', [], [], Spcm.signalUnits(), obj.NAME);
                return
            end
            S1 = squeeze(obj.signal(1, :, 1:obj.currIter));
            S1e = squeeze(obj.sterr(1, :, 1:obj.currIter));
            S2 = squeeze(obj.signal(2, :, 1:obj.currIter));
            S2e = squeeze(obj.sterr(2, :, 1:obj.currIter));
            if obj.currIter ~= 1
                S1e = obj.getCombinedSterr(S1, S1e);
                S2e = obj.getCombinedSterr(S2, S2e);
                S1 = mean(S1, 2);
                S2 = mean(S2, 2);
            end
            refName = BooleanHelper.ifTrueElse(strcmp(obj.referenceMode, 'noMW'), ...
                                               'Reference (no MW)', 'Reference (laser tail)');
            p1 = ExpResultDoubleVector('fluorescence', S1, S1e, Spcm.signalUnits(), obj.NAME, 'After MW pulse');
            p2 = ExpResultDoubleVector('fluorescence', S2, S2e, Spcm.signalUnits(), obj.NAME, refName);
            params = {p1, p2, {'yParam', 'yParam2'}};
        end
    end

    %% Helpers
    methods (Access = protected)
        function a = refAmplitude(obj)
            a = obj.piTimeAmplitude;
            if isempty(a); a = obj.amplitude(1); end
        end

        function tp = pulseLength(obj)
            tp = obj.mwPulseDuration;
            if isempty(tp); tp = obj.piTime; end
            if isempty(tp)
                obj.sendError('No MW pulse length: set mwPulseDuration, or piTime (loadJsonParams sets it)');
            end
        end

        function p = timelineParams(obj)
            if obj.seperatedTrigger
                trig = {'trigger', 'trigger2'};
            else
                trig = {'trigger'};
            end
            if obj.nChannels == 2
                gate = {'MW', 'MW2'};
            else
                gate = {obj.MWChannel};
            end
            if obj.AWGorSRSswitch; gate = [gate, {'AWG'}]; end
            if ~obj.driveMWSwitches; gate = {}; end
            p = struct('darkBeforeMW', obj.darkBeforeMW, ...
                       'pulse', obj.pulseDurationActual, ...
                       'lastDelay', obj.lastDelay, ...
                       'detectionStartDelay', obj.detectionStartDelay, ...
                       'detectionDuration', obj.detectionDuration, ...
                       'laser', obj.laserInitializationDuration, ...
                       'lead', obj.mwSwitchLead, 'tail', obj.mwSwitchTail, ...
                       'trigDuration', obj.AWG.TRIG_DURATION, ...
                       'referenceMode', obj.referenceMode, ...
                       'referenceStartDelay', obj.referenceStartDelay, ...
                       'gate', {gate}, 'trig', {trig});
        end

        function [fp, parts] = awgFingerprint(obj)
            % Everything that decides what is uploaded, as named parts (so a miss can say
            % which one moved) and their hash. The class name keeps it from ever matching
            % a CW ESR or Rabi set on the same instrument. Walk-only settings enter only in
            % walk mode, so toggling them never forces a classic reload.
            % FIELD BY FIELD: struct(...) would expand a cell value into an array.
            walk = obj.useSequenceWalk;
            parts = struct();
            parts.cls            = 'AWGExpESRPulsed';
            parts.frequency      = obj.frequency;
            parts.amplitude      = obj.amplitude;
            parts.pulse          = obj.pulseDurationActual;
            parts.nChannels      = obj.nChannels;
            parts.activeChannels = obj.activeChannels;
            parts.timeDelay      = obj.timeDelay;
            parts.channelPhase   = obj.channelPhase;
            parts.sampleRate     = obj.AWG.SampleRate;
            parts.minWaveform    = obj.AWG.minWaveformDuration;
            parts.walk           = walk;
            parts.walkWarmup     = double(walk) * obj.walkWarmupTriggers;
            parts.walkRandom     = walk && obj.randomizeFrequencyOrder;
            fp = DataHash(parts);
        end

        function updateSignalParam(obj)
            N  = obj.getTotalNumberOfParams();
            S1 = squeeze(obj.signal(1, :, 1:obj.currIter));
            S1e = squeeze(obj.sterr(1, :, 1:obj.currIter));
            S2 = squeeze(obj.signal(2, :, 1:obj.currIter));
            S2e = squeeze(obj.sterr(2, :, 1:obj.currIter));
            if exist('ratioAcrossAverages', 'file')
                % Ratio per average first, where level drift cancels, then combined.
                [value, sterr] = ratioAcrossAverages(S1, S2, S1e, S2e, N);
            else
                [value, sterr] = obj.getRatioDistributionValues(S1, S2, S1e, S2e);
            end
            obj.signalParam.value = value;
            obj.signalParam.sterr = sterr;
        end

        function reportTiming(obj)
            T  = obj.timing;
            tp = obj.pulseDurationActual;
            fprintf('\n--- %s ---\n', obj.NAME);
            fprintf(['repeat %.3f us: [dark %.3f][MW pulse %.4f][wait %.3f][laser %.3f: ' ...
                     'signal window %.3f-%.3f'], T.repeatDuration, obj.darkBeforeMW, tp, ...
                    obj.lastDelay, obj.laserInitializationDuration, T.sigWindow - T.laserOn, ...
                    T.sigWindow - T.laserOn + obj.detectionDuration);
            if strcmp(obj.referenceMode, 'noMW')
                fprintf('] x2, the second half without MW (reference)\n');
            else
                fprintf(', reference window %.3f-%.3f]\n', T.refWindow - T.laserOn, ...
                        T.refWindow - T.laserOn + obj.detectionDuration);
            end

            % Is it a pi pulse? The rotation scales with the pulse length AND with the field
            % amplitude, relative to the (piTime, piTimeAmplitude) pair from the json.
            if ~isempty(obj.piTime)
                rot = obj.rotationDegrees(tp);
                if abs(rot / 180 - 1) > 0.05
                    fprintf(['NOTE: %.4f us at %.1f dBm is a ~%.0f-degree rotation on resonance, ' ...
                             'not pi.\n     For a pi pulse of %.4f us use amplitude ~%.1f dBm ' ...
                             '(amplitudeForPulse; check with a Rabi).\n'], tp, obj.amplitude(1), ...
                            rot, tp, obj.amplitudeForPulse(tp));
                end
            end
            fw = 0.8 / tp;
            df = median(abs(diff(sort(obj.frequency))));
            fprintf('line FWHM from the pulse alone ~%.2f MHz (0.8/tau); frequency step %.3g MHz', fw, df);
            if ~isempty(df) && df > fw / 3
                fprintf('  <-- coarse: fewer than 3 points across the line');
            end
            fprintf('\n');
            if obj.darkBeforeMW < 0.5
                fprintf(['NOTE darkBeforeMW = %.2f us: the metastable singlet (~0.2-0.3 us) may ' ...
                         'not have emptied, which costs contrast.\n'], obj.darkBeforeMW);
            end
            if T.minWindowGap < obj.minDetectionGap
                obj.sendWarning(sprintf(['Two detection windows are only %.2f us apart (minDetectionGap ' ...
                    '= %g us). The digitizer may miss a trigger re-arming: lengthen ' ...
                    'laserInitializationDuration or use referenceMode ''noMW''.'], ...
                    T.minWindowGap, obj.minDetectionGap));
            end
        end

        function checkNonNegScalar(obj, v, name)
            if ~(isnumeric(v) && isscalar(v) && isfinite(v) && v >= 0)
                obj.sendError(sprintf('%s must be a non-negative number of us', name));
            end
        end
    end

    methods (Static, Hidden)
        % plotPulsedSequence's styling
        function C = seqColors()
            C = struct();
            C.laser    = [0.18 0.62 0.29];
            C.det      = [0.42 0.27 0.72];
            C.gate     = [0.93 0.55 0.13];
            C.trig     = [0.35 0.35 0.38];
            C.other    = [0.55 0.55 0.55];
            C.dim      = [0.45 0.45 0.45];
            C.refShade = [0.94 0.94 0.95];
            C.rf       = [0.80 0.15 0.15; 0.12 0.38 0.75];
            C.rfShade  = [1.00 0.90 0.90; 0.88 0.93 1.00];
        end

        function [name, col] = seqRowStyle(ch, C)
            switch ch
                case 'greenLaser'; name = 'laser';        col = C.laser;
                case 'detector';   name = 'detector';     col = C.det;
                case 'MW';         name = 'MW switch';    col = C.gate;
                case 'MW2';        name = 'MW2 switch';   col = C.gate*0.85;
                case 'AWG';        name = 'AWG switch';   col = C.gate*0.7;
                case 'trigger';    name = 'AWG trigger';  col = C.trig;
                case 'trigger2';   name = 'AWG trigger2'; col = C.trig;
                otherwise;         name = ch;             col = C.other;
            end
        end

        function s = seqDur(d)
            if d >= 1
                s = sprintf('%.2f \\mus', d);
            else
                s = sprintf('%.0f ns', d*1e3);
            end
        end
    end

    methods (Static)
        function [iv, T] = buildTimeline(p)
            % The repeat as absolute-time intervals [t, d, channels], laid onto an empty
            % sequence with addEventAtGivenTime (which unions the channels where they
            % overlap -- the same way PulseGenerator.fixSequenceDelays builds sequences).
            % Static and free of hardware, so it can be checked on its own.
            %
            % Computed in WHOLE NANOSECONDS (the pulse generator's resolution) and only
            % converted to us at the end, so edges that should coincide coincide exactly
            % instead of leaving 1e-16 us slivers for the pulse generator to choke on.
            ns = @(x) round(x * 1e3);
            dark = ns(p.darkBeforeMW);  tp = ns(p.pulse);      lastD = ns(p.lastDelay);
            las  = ns(p.laser);         dStart = ns(p.detectionStartDelay);
            det  = ns(p.detectionDuration);
            lead = ns(p.lead);          tail = ns(p.tail);     trigD = ns(p.trigDuration);
            if tp < 1
                error('AWGExpESRPulsed:timing', 'The MW pulse rounds to 0 ns.');
            end
            if lead > dark
                error('AWGExpESRPulsed:timing', ['mwSwitchLead (%.3f us) is longer than ' ...
                      'darkBeforeMW (%.3f us): the switch would open while the laser is on.'], ...
                      lead/1e3, dark/1e3);
            end
            if dStart + det > las
                error('AWGExpESRPulsed:timing', ['The signal window ends %.3f us after the ' ...
                      'readout laser comes on, but the laser is only %.3f us long.'], ...
                      (dStart + det)/1e3, las/1e3);
            end

            half     = dark + tp + lastD + las;
            tPulse   = dark;
            laserOn  = dark + tp + lastD;                      % within a half
            sigStart = laserOn + dStart;

            iv = struct('t', {}, 'd', {}, 'ch', {});
            add = @(iv, t, d, ch) [iv, struct('t', t, 'd', d, 'ch', {ch})];

            % --- signal half: MW gate, trigger, readout ---
            if ~isempty(p.gate)
                iv = add(iv, tPulse - lead, lead + tp + tail, p.gate);
            end
            iv = add(iv, tPulse, trigD, p.trig);
            iv = add(iv, laserOn, las, {'greenLaser'});
            iv = add(iv, sigStart, det, {'detector'});

            switch p.referenceMode
                case 'noMW'
                    % identical half, no MW
                    iv = add(iv, half + laserOn, las, {'greenLaser'});
                    refStart = half + sigStart;
                    iv = add(iv, refStart, det, {'detector'});
                    repeatDuration = 2 * half;
                case 'laserTail'
                    if isempty(p.referenceStartDelay)
                        rs = las - det - 50;                    % 50 ns before the laser goes off
                    else
                        rs = ns(p.referenceStartDelay);
                    end
                    if rs < dStart + det
                        error('AWGExpESRPulsed:timing', ['The reference window (from %.3f us) ' ...
                              'overlaps the signal window (until %.3f us).'], rs/1e3, ...
                              (dStart + det)/1e3);
                    end
                    if rs + det > las
                        error('AWGExpESRPulsed:timing', ['The reference window ends after the ' ...
                              'readout laser (%.3f > %.3f us).'], (rs + det)/1e3, las/1e3);
                    end
                    refStart = laserOn + rs;
                    iv = add(iv, refStart, det, {'detector'});
                    repeatDuration = half;
                otherwise
                    error('AWGExpESRPulsed:timing', 'Unknown referenceMode %s', p.referenceMode);
            end

            % Nothing may run past the end of the repeat (the trigger or the switch tail
            % could, with a very short repeat): clip it.
            for k = 1:numel(iv)
                iv(k).d = min(iv(k).d, repeatDuration - iv(k).t);
            end
            iv = iv([iv.d] > 0);

            % Gaps between consecutive windows, including signal-after-reference across
            % the repeat boundary.
            g1 = refStart - (sigStart + det);
            g2 = repeatDuration - (refStart + det) + sigStart;

            % to us
            for k = 1:numel(iv)
                iv(k).t = iv(k).t / 1e3;  iv(k).d = iv(k).d / 1e3;
            end
            T = struct('repeatDuration', repeatDuration/1e3, 'halfDuration', half/1e3, ...
                       'pulseStart', tPulse/1e3, 'laserOn', laserOn/1e3, ...
                       'sigWindow', sigStart/1e3, 'refWindow', refStart/1e3, ...
                       'minWindowGap', min(g1, g2)/1e3);
        end
    end
end
