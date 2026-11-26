classdef AWGExpESR_gapAom < ExperimentAWG
    %AWGEXPESR_GAPAOM ESR (electron spin resonance) Experiment with a separate laser
    % power (AOM value) during delayBetweenDetections.
    %
    % Identical to AWGExpESR except for the gap power. With gapAomValue = [] (the
    % default) it IS AWGExpESR: same sequence, same AWG upload, same data.
    %
    % ===================================================================
    % HOW IT WORKS, AND THE ONE CABLE IT NEEDS
    % ===================================================================
    % On Setup 7 the green AOM level is a static voltage from the NI-DAQ (ao0, set by
    % the laser panel). Software cannot change it within a 10 us gap, so this class
    % makes ao0 HARDWARE-TIMED for the duration of the run:
    %
    %   * the pulse sequence emits a short marker pulse on PulseStreamer output
    %     gapAomPgChannel at the START of each gap and again gapAomReturnLead before
    %     its END (4 marker pulses per repeat);
    %   * that output is cabled to the DAQ input gapAomDaqTerminal, and ao0 runs as a
    %     regenerating 2-level waveform clocked by it: every rising edge steps ao0 to
    %     the next value, alternating gap value / window value.
    %
    %   CABLE (BNC): PulseStreamer output gapAomPgChannel (default 1)  ->  DAQ
    %                gapAomDaqTerminal (default PFI0).
    %
    % Edges come in pairs, so every pg.run ends on the window value. The task is
    % re-armed at the start of every average (which also re-aligns it if a pair was
    % ever broken), and on wrapUp ao0 is released and set back to the panel value.
    %
    % While a run is armed, the laser panel cannot write the AOM (the DAQ channel is
    % in use); change the power between runs. After a crash or Ctrl-C, call
    % a.releaseGapAom to free ao0 and restore the window value.
    %
    % USAGE
    %   a = AWGExpESR_gapAom;  ... same settings as AWGExpESR ...
    %   a.delayBetweenDetections = 20;   % us, both gaps
    %   a.gapAomValue = 0;               % panel units (Setup 7: volts, -0.6 .. 0)
    %   a.run
    %
    % CHECK IT ONCE ON A SCOPE: put the AOM control line (ao0) and the marker on a
    % scope, or a photodiode in the beam. Adjust gapAomReturnLead so the power is back
    % at the window value before each window starts. The class also counts the edges
    % the DAQ received and warns if it saw none (cable missing or wrong terminal).

    properties (Constant)
        NAME = 'AWG ESR_gapAom'

        GAP_AOM_PG_NAME = 'gapAomClock'   % PulseStreamer channel name this class registers

        ZERO_FIELD_SPLITTING = 2.87e3     % in Mhz
    end

    properties
        % All of the parameters must have a set function, and when they are
        % changed, the parameter obj.changeFlag needs to be set to true
        % (unless the change is minor and the experiment can continue
        % running);
        mirrorSweepAround   % double. in MHz

        mode                % string. Either 'CW' or 'pulsed'
                            % ('pulsed is to be implemented in the future, if needed)
        freqMirrored

        isSingleMeasurement
        delayBetweenDetections

        temp

        useSequenceJump = false;    % logical. Selects how the sweep changes frequency.
                    %
                    % false (default) - CLASSIC path. One AWG sequence per
                    %   (channel, frequency); each point calls assignSequence.
                    %   Measured on the AWG70002A: assigning an asset STOPS
                    %   playback and re-arming costs ~500 ms whatever is loaded,
                    %   so a 191-point x 2-average sweep spent 176 s of its 178 s
                    %   inside waitForRunning.
                    %
                    % true - JUMP path. ONE sequence per channel holds every
                    %   frequency as a consecutive 3-step block, each block's
                    %   last step pointing back at its own first step so the
                    %   sequencer parks there for all 'repeats' triggers. Each
                    %   point then costs a single SOURce:JUMP:FORCe (~3 ms) and
                    %   the AWG is armed once for the whole sweep.
                    %
                    % Start with false, run your reference spectrum, then flip to
                    % true and check you get the SAME spectrum before trusting it.
                    %
                    % Requires repeats > 10. A forced jump can land mid-burst, so
                    % the first repeat after a jump may be corrupt -- which is
                    % exactly what processData already discards (start = 2 when
                    % repeats > 10). Below that threshold nothing is discarded and
                    % a bad repeat would enter the average.
                    %
                    % Not yet supported together with mirrorSweepAround.

        useSequenceWalk = false;    % logical. The fast acquisition mode.
                    %
                    % Measured cause of the slow sweep: assignSequence stops the
                    % AWG, and re-arming after an ASSET CHANGE costs ~499 ms,
                    % against ~97 ms after a plain stop. So the per-point cost was
                    % never arming latency -- it was re-assigning the sequence.
                    %
                    % Walk mode assigns ONE sequence per channel holding every
                    % frequency as a consecutive 3-step block, arms once, and lets
                    % the sequencer fall through to the next block on each trigger.
                    % A single pg.run therefore covers a whole pass over the sweep,
                    % and the per-point AWG cost is zero.
                    %
                    % It uses only mechanisms the classic path already proves work
                    % on this instrument: a step that waits for a trigger, and
                    % fall-through to the next step with one 'goto first' at the
                    % end. Deliberately NOT useSequenceJump, which relied on
                    % SOURce:JUMP:FORCe plus a self-referencing GOTO and produced
                    % frequency-scrambled data (correlation +0.04 with the classic
                    % spectrum at every settle time tried).
                    %
                    % The data comes back interleaved as [s1 r1 s2 r2 ... sN rN]
                    % per pass and is reduced by walkReduce, which is a faithful
                    % transcription of Experiment.processData (see testWalkReduce).
                    %
                    % Validate with esrWalkValidate before trusting it.
                    % Not supported together with mirrorSweepAround or
                    % useSequenceJump.

        passesPerRun = 1;           % How many full passes over the sweep per pg.run.
                    % 1 needs 2*nFreq digitizer segments -- 192 for a 96-point
                    % sweep, essentially the same as the 200 the classic path
                    % already commits, so it is safe without knowing the card's
                    % memory limit. Raising it cuts pulse-generator round trips
                    % but multiplies the segment count; CsMl_Commit will error if
                    % you exceed the board's memory.

        resyncEveryPass = true;     % logical. Stop and re-arm the AWG at the start of
                    % each pass, so the sequencer is known to be on frequency 1.
                    % Costs the ~97 ms plain-stop arm measured on the bench, and
                    % bounds any trigger-count error to a single pass instead of
                    % shifting every pass after it. Turn off only once you trust
                    % the sync.

        driveMWSwitches = true;     % logical. Drive the MW switch channels at all.
                    %
                    % DIAGNOSTIC. Leave it true for real measurements.
                    %
                    % The setup json declares two switches on pulseStreamer channels
                    % 5 and 6 ("MW", "MW2"), and the sequence asserts them during the
                    % signal window only. So on paper there are TWO gates in series:
                    % those switches, and the AWG playing zeros between bursts.
                    %
                    % FGEN mode removes the AWG's zeros and relies on the switches
                    % alone -- and produces no ESR signal at all, at verifiably the
                    % same output amplitude, with the same switches and amplifiers
                    % downstream. That is only consistent with the switches not
                    % actually gating.
                    %
                    % Setting this false tests that directly, on the path that works.
                    % The MW switch channels are left unasserted for the whole
                    % sequence, so the switches are never told to open:
                    %
                    %   spectrum UNCHANGED -> the switches are not gating anything.
                    %     The AWG's zeros have been the only gate all along, and FGEN
                    %     cannot work on this setup until a working switch sits in
                    %     the AWG's output path.
                    %
                    %   spectrum goes FLAT -> the switches are gating, and correctly:
                    %     with them never opened, no microwave reaches the NVs. Then
                    %     FGEN's failure is upstream of them -- the unexplained
                    %     'run state = 0' being the obvious suspect.
                    %
                    % Note this does NOT change what the AWG uploads, so it costs no
                    % reload: the pulse sequence is rebuilt on every prepare anyway.

        useFgen = false;            % logical. Use the AWG's built-in function
                    % generator instead of uploaded waveforms. AWG mode is the
                    % default and nothing else changes unless you set this true.
                    %
                    % In FGEN mode the instrument synthesises the sine itself, so
                    % there are NO waveforms, NO sequence, NO build and NO first
                    % arm: the ~130 s cold load disappears completely. A frequency
                    % change is one SCPI write, so the per-point cost is a few ms.
                    %
                    % Measured on this AWG70002A by awgFgenProbe:
                    %   sine accepted to at least 6 GHz, well past the sweep
                    %   per-channel amplitude independent, in dBm (FGEN:CHAN<n>:AMPL:POW)
                    %   per-channel phase independent, 0.1 deg step
                    %   frequency COUPLED across channels -- irrelevant here, both
                    %     channels always play the same frequency anyway
                    %
                    % timeDelay is carried by PHASE rather than by padding a zero
                    % waveform, since there is no waveform. For a sine those are the
                    % same thing: delay = phase / (360 f). At 2870 MHz the 0.1 deg
                    % step is 0.097 ps, against 62.5 ps (one sample) in AWG mode --
                    % so the delay is finer here, not coarser. The catch is RANGE:
                    % phase wraps at +/-180 deg, which is +/-174 ps at 2870 MHz. A
                    % larger timeDelay cannot be expressed and is rejected rather
                    % than silently wrapped.
                    %
                    % WHY THE CONTRAST SHOULD SURVIVE. In FGEN mode the MW is on
                    % continuously and the signal/reference contrast rests entirely
                    % on the switch the pulse generator drives. That is already the
                    % case: the switch is on for detectionDuration while the AWG
                    % burst runs for detectionDuration + laserInitializationDuration,
                    % so the burst OVERHANGS the switch window and both edges of the
                    % effective MW window are cut by the switch, not by the AWG.
                    % Continuous MW is just an unlimited overhang.
                    %
                    % The residual risk is leakage: MW now sits at the switch input
                    % during the reference window too, so whatever gets through a
                    % closed switch adds a little MW where there should be none,
                    % which would shrink contrast. That is a number -- compare
                    % contrast against an AWG-mode run before trusting it.
                    %
                    % Mutually exclusive with useSequenceWalk and useSequenceJump.

        fgenRunTimeout = 10;        % s. How long to wait for FGEN to reach run state.
                    % Measured by awgFgenBringUp: after AWGCONTROL:RUN the state is
                    % still 0 at 0.2 s and only reads 2 at about 1.2 s. Reading it
                    % immediately -- which is what the first version did -- always
                    % says "not running", whether it is or not.

        fgenSettleTime = 0.001;     % s. Pause after changing the FGEN frequency.
                    % Jump mode is the cautionary tale here: it looked fast and
                    % returned frequency-scrambled data. If an FGEN spectrum looks
                    % shifted or smeared, raise this before concluding anything
                    % else -- and if raising it changes nothing, settling is not the
                    % problem, which is exactly how jump mode was ruled out.
                    %
                    % MEASURED, at the scope, by fgenStepTime: 1 ms is far too short.
                    % A bare frequency write leaves no usable tone at the new frequency
                    % for at least 5 s; a full restart needs ~4.1 s to be stable. So
                    % this default is why FGEN ESR read 0.039% -- not a code fault, a
                    % settling fault. It is left at 1 ms because FGEN is off by default
                    % and the value only matters to a deliberate diagnostic.

        randomizeFrequencyOrder = true;   % logical. Visit frequencies in random order.
                    % TRUE for all real data: random order decorrelates slow drift from
                    % frequency, so a warming amplifier cannot look like a spectral
                    % feature. Set FALSE only for diagnostics, where it converts a
                    % one-step visit lag into a plain frequency shift that the
                    % best-shift metric can see.
                    %
                    % What "random" means depends on the path:
                    %   classic (useSequenceWalk = false): a NEW order every average.
                    %   walk: ONE random order per upload, used by every pass and every
                    %     average -- free, because it is the same sequence, just in a
                    %     different order. The data is put back in frequency order before
                    %     anything sees it. Same scheme as AWGExpRabi/Ramsey/Echo walk.
                    %     For a new order every average in walk mode, see
                    %     newOrderEveryAverage.

        newOrderEveryAverage = true;    % logical. WALK MODE ONLY. Draw a new random order
                    % before every average and re-upload ONLY the sequence (the waveforms
                    % stay on the instrument). Each average then pays one sequence upload
                    % + arm, printed as it happens -- that is the same cost as the first
                    % "sequence" + "first arm" numbers in the Cold load breakdown.
                    % What it buys over one order per upload: effects tied to position in
                    % the pass (start-of-pass transient, the previous block's frequency)
                    % average down across averages instead of sitting on the same
                    % frequencies for the whole run. Needs randomizeFrequencyOrder = true.

        lastVisitOrder = [];        % read-only record of the last run's visit order.
                    % lastVisitOrder(k) is the frequency INDEX visited kth. Needed to
                    % test the jump-mode lag hypothesis without fitting anything: if the
                    % value filed under frequency lastVisitOrder(k) actually belongs to
                    % lastVisitOrder(k-1), that remap is fully determined by this vector
                    % and has no free parameters at all.

        fgenFullRestart = false;    % logical. Repeat the whole generator bring-up --
                    % RUN:IMMEDIATE, poll the run state, re-enable the outputs -- at
                    % EVERY frequency, instead of writing the frequency alone.
                    %
                    % Needed because a bare FGEN:CHANn:FREQuency write does not produce
                    % a stable tone at the new frequency on this instrument (fgenStepTime
                    % watched it at the scope: noise floor for ~2.5 s, a brief
                    % appearance, then gone again, at every step from 2 to 800 MHz).
                    % With this true the tone is stable ~4.1 s after each change.
                    %
                    % This is a DIAGNOSTIC, not a mode to run in. ~4.5 s per point makes
                    % a 96-point sweep over seven minutes of pure settling against ~30 s
                    % for the whole AWG-mode run, so FGEN cannot win on speed even when
                    % it works. Its only purpose is to answer whether settling was the
                    % entire explanation for the missing contrast.

        channelPhase = [0 0];       % DEGREES, per AWG channel. The relative phase of
                    % the two microwave drives, put into the WAVEFORM instead of
                    % into the length of a padding waveform.
                    %
                    % This replaces timeDelay as the phase knob, and the reason is
                    % arithmetic. timeDelay is realised by lengthening one channel's
                    % leading zero waveform, and AddWaveformByFunction rounds that
                    % duration to whole samples. At 16 GS/s one sample is 62.5 ps,
                    % and one microwave period at 2870 MHz is 348.4 ps -- so one
                    % sample is 64.6 degrees of phase. The achievable settings are
                    % roughly 0, 65, 129, 194, 258, 323 degrees: about 5.6 points
                    % per period, so the best you can do is +/-32 degrees. That was
                    % the accuracy ceiling on channel balance, and it was invisible
                    % because the scan was in picoseconds while the rounding
                    % happened silently inside the AWG class.
                    %
                    % Written into the cosine, the resolution is limited only by
                    % floating point. No extra samples, no extra upload, nothing
                    % else changes.
                    %
                    % Note also that the stored timeDelay of 327 ps is 337.9
                    % degrees, i.e. within 6% of a full period -- physically almost
                    % the same as no phase shift at all. Calibrating in picoseconds
                    % hides that the quantity is periodic; in degrees it is obvious.
                    %
                    % Leave timeDelay at 0 when using this. Both mechanisms are
                    % additive and having both non-zero means the actual phase is
                    % the sum, which nothing records.

        walkWarmupTriggers = 1;     % Throwaway blocks at the START of each pass.
                    %
                    % This is the fix for the "first frequency disagrees by 5.6
                    % sigma" result. The two paths were NOT discarding the same
                    % thing. The classic path arms, then fires 'repeats' triggers
                    % AT ONE FREQUENCY, and processData drops repeat 1 -- so at
                    % every frequency the trigger that follows the arm/idle gap is
                    % thrown away. Walk mode fires one trigger per frequency and
                    % walkReduce drops the first PASS, which is a different thing
                    % entirely: within every pass the trigger that follows the gap
                    % is the one belonging to frequency 1, so frequency 1 carried
                    % an un-discarded warm-up on all 100 passes while the other 95
                    % frequencies were clean. Dropping a pass cannot repair that.
                    %
                    % So the sequence now begins with N extra blocks that are
                    % copies of frequency 1, each pass fires nFreq+N triggers, and
                    % the first N read pairs of every pass are discarded. Every
                    % frequency that is kept, including frequency 1, then sits
                    % behind a played-and-thrown-away block exactly like the
                    % classic path.
                    %
                    % 1 is the right number unless the instrument needs more than
                    % one block to settle. 0 restores the old, biased behaviour.

        walkCompactSteps = true;    % logical. Two sequence steps per frequency
                    % instead of three. This is where the walk-mode startup time
                    % goes.
                    %
                    % Measured: the first arm after a load is linear in step count
                    % at ~0.16 s/step, and building the sequence inside the
                    % instrument is superlinear. A 96-point sweep at 3 steps each
                    % is 288 steps -> 47 s to build plus 54 s to arm, which is the
                    % ~100 s of dead time before acquisition starts.
                    %
                    % The three steps are [leading zero (waits for trigger), MW
                    % burst, trailing zero]. In walk mode the trailing zero is
                    % redundant: it plays immediately after the burst and is then
                    % followed by a step that waits for the next trigger anyway,
                    % so the AWG outputs zero there either way. Dropping it leaves
                    %   [leading zero (waits for trigger), MW burst]
                    % and 2*(nFreq+walkWarmupTriggers) steps -- 194 instead of 288
                    % for 96 points, so roughly a third off both numbers.
                    %
                    % Jump mode always keeps three steps: there the trailing zero
                    % carries the GOTO that makes a block loop on itself.
                    %
                    % If a walk spectrum ever looks wrong with this on, set it
                    % false and re-run esrWalkValidate -- that isolates the layout
                    % from everything else in one step.

        burstTargetDuration = [];   % us. Length of ONE looped MW burst waveform.
                    %
                    % This is the lever on the upload time, and it is a much better
                    % one than the sample rate.
                    %
                    % The MW burst is not uploaded as a whole 110 us tone. One
                    % short waveform holding a whole number of MW cycles is
                    % uploaded, and the sequencer loops it numSegments times. Only
                    % the SHORT waveform costs transfer time; the looping is free.
                    %
                    % The original code fixed that short waveform at ~1 us
                    % (round(f)/f, with f in MHz), which at 16 GS/s is 16000
                    % samples -- and there are 2 channels x nFreq of them, so a
                    % 96-point sweep transfers about 3.1 million samples. That is
                    % where the ~59 s of waveform time in a cold load goes.
                    %
                    % The instrument's floor is MIN_WAVEFORM_LENGTH = 2400 samples.
                    % [] (auto) picks the shortest whole number of MW cycles that
                    % clears that floor with 5% margin -- about 0.158 us and 2520
                    % samples at 16 GS/s, so roughly 6x less data to move.
                    %
                    % Nothing about the emitted microwave changes: same frequency,
                    % a whole number of cycles so the loop is seamless, and the
                    % same total burst length to within one segment. Sample rate,
                    % per-channel amplitude and the timeDelay resolution are all
                    % untouched -- which matters, because those are what you are
                    % about to calibrate.
                    %
                    % Set to 1 for the original ~1 us behaviour.

        armCheckEvery = [];         % How often the sweep verifies the AWG is still armed.
                    %
                    % Each check is one AWGControl:RSTate? query, and on the
                    % AWG70002A over VISA TCPIP that query costs ~1.2 s. Asking
                    % on every one of 191 points therefore cost 235 s -- which
                    % after jump mode was 80% of the entire run, dwarfing both
                    % the acquisition and the sweep itself.
                    %
                    % []  auto: 1 on the classic path (where the AWG really does
                    %     stop at every point, so the check is load-bearing) and
                    %     25 on the jump path (where it stopped 0 times out of
                    %     193, so the check is nearly pure cost).
                    % N   check on the first point of each average and then every
                    %     Nth point. A retry always forces a check.
                    % 0   never check. Fastest, but if the AWG ever does stop
                    %     mid-sweep you get silent garbage rather than a re-arm.
                    %     Only once you trust the jump path.

        % ---- Laser power during the gaps (this class only) --------------------------
        gapAomValue = [];           % AOM value during delayBetweenDetections, in the laser
                    % panel's units (Setup 7: volts, -0.6 .. 0 -- the Green Laser aomControl
                    % range in the json). [] = off: identical to AWGExpESR.
                    % The detection windows keep the panel value at the start of the run;
                    % it is restored when the run ends.
        gapAomReturnLead = 2;       % us. The AOM is switched back to the window value this
                    % long BEFORE each detection window, so the window starts at the right
                    % power. Covers the DAQ update (~1 us) and the AOM response. Verify once
                    % on a scope and lengthen it if the window start is still moving.
        gapAomEdgeWidth = 0.2;      % us. Width of each marker pulse.
        gapAomPgChannel = 1;        % PulseStreamer output (0-7) carrying the marker. Must
                    % be free: the Setup 7 json uses 0, 4, 5, 6, 7 (and lists 2 as 'dummy';
                    % the beam-sample attenuator json also uses 2).
        gapAomDaqTerminal = 'PFI0'; % DAQ terminal that receives the marker (ao0's sample
                    % clock). Any free PFI line on the DAQ.
    end

    properties (SetAccess = protected)
        seqAssetIndex       % 1 x nChannels. Index of the single big sequence per
                            % channel, used only on the useSequenceJump path.
                            % NOTE on that path obj.indexSeq holds STEP numbers,
                            % not sequence indices.
        walkOrder = []      % 1 x nFreq. Walk/jump: the frequency INDEX each block plays,
                            % in sequence order. performWalk uses it to put the data back
                            % in frequency order. Stored in the load cache with the upload.
        walkParts = []      % struct. What uploadWalkSequence needs to rebuild the sequence
                            % in a new order without re-uploading the waveforms.
        walkOrderFresh = false  % true while the loaded order has not been used for an
                            % average yet, so newOrderEveryAverage does not redraw it.

        gapAomTask = []         % DAQmx task handle while the gap power is armed
        gapAomWindowPanel = []  % panel value held by the windows; restored after the run
        gapAomVolts = []        % [gap, window] in volts, as written to the DAQ
        gapAomEdges = NaN       % edges the DAQ received in the last armed period
    end

    methods

        function obj = AWGExpESR_gapAom(passedExpName)

            if nargin == 1
               expName = passedExpName;
            else
               expName = AWGExpESR_gapAom.NAME;
            end

            obj@ExperimentAWG(expName);

            obj.parameterName = 'frequencies';

            obj.repeats = 100;
            obj.averages = 1000;

            obj.frequency = obj.ZERO_FIELD_SPLITTING + (-2 : 2 : 2);     %in MHz
            obj.amplitude = -25 + [6 0];        % dBm
            obj.mode = 'CW';            % Can only be 'CW' for now.

            obj.detectionDuration = 500;
            obj.referenceDetectionDuration = 500;
            obj.laserInitializationDuration = 10; % laser initialization
            spcm = getObjByName(Spcm.NAME);
            if isa(spcm, 'PhotoDiodeDigitizerNiDaqControlled')
                obj.delayBetweenDetections = 10;
            else
                obj.delayBetweenDetections = 0;
            end

            obj.mirrorSweepAround = []; % Use this for a single frequency range

            obj.displayType2 = 'Raw signal (kcps)';   % the window shows 'Raw signal (V)' with the photodiode (Spcm.localiseUnits)
            obj.isPlotAlternateAvailable = true;

            obj.mCurrentXAxisParam = ExpParamDoubleVector('Frequency', [], [], 'MHz', obj.NAME);
            % Labels as in the robust fit's figure (ExpPlotStyle / fitStyle)
            obj.plotXLabel = ExpPlotStyle.X_FREQUENCY;
            obj.plotYLabel = ExpPlotStyle.Y_NORMALISED;
            obj.signalParam = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME, 'Normal');
            obj.signalParam2 = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME, 'Mirrored');
        end
    end

    %% Setters
    methods
        function set.mode(obj, newVal)
            checkMode(obj, newVal)
            % If we got here, then the mode was changed
            obj.mode = newVal;
            obj.changeFlag = true;
        end

        function set.useSequenceWalk(obj, newVal)
            if ~(islogical(newVal) || (isnumeric(newVal) && isscalar(newVal) && ismember(newVal, [0 1])))
                obj.sendError('useSequenceWalk must be true or false')
            end
            obj.useSequenceWalk = logical(newVal);
            obj.changeFlag = true;      % the AWG has to be reprogrammed either way
        end

        function set.passesPerRun(obj, newVal)
            if ~ValidationHelper.isValuePositiveInteger(newVal)
                obj.sendError('passesPerRun must be a positive integer')
            end
            obj.passesPerRun = newVal;
            obj.changeFlag = true;      % changes pg.repeats and the digitizer segment count
        end

        function set.resyncEveryPass(obj, newVal)
            if ~(islogical(newVal) || (isnumeric(newVal) && isscalar(newVal) && ismember(newVal, [0 1])))
                obj.sendError('resyncEveryPass must be true or false')
            end
            obj.resyncEveryPass = logical(newVal);
        end

        function set.channelPhase(obj, newVal)
            if ~isnumeric(newVal) || ~isvector(newVal)
                obj.sendError('channelPhase must be a numeric vector, in degrees')
            end
            % Wrapped into [0,360) on the way in. The phase is periodic, and letting
            % it run outside one period is how a scan ends up with several equivalent
            % optima and reports whichever one it happened to reach first.
            obj.channelPhase = mod(newVal(:)', 360);
            obj.changeFlag = true;      % the waveform samples change
        end

        function set.walkWarmupTriggers(obj, newVal)
            if ~(isnumeric(newVal) && isscalar(newVal) && newVal >= 0 && mod(newVal,1) == 0)
                obj.sendError('walkWarmupTriggers must be a non-negative integer')
            end
            obj.walkWarmupTriggers = newVal;
            obj.changeFlag = true;      % changes the sequence AND pg.repeats
        end

        function set.driveMWSwitches(obj, newVal)
            if ~(islogical(newVal) || (isnumeric(newVal) && isscalar(newVal) && ismember(newVal, [0 1])))
                obj.sendError('driveMWSwitches must be true or false')
            end
            obj.driveMWSwitches = logical(newVal);
            % Deliberately does NOT set changeFlag: this only affects the pulse
            % sequence, which prepare rebuilds every time, so forcing an AWG
            % reload would cost ~40 s for nothing.
        end

        function set.useFgen(obj, newVal)
            if ~(islogical(newVal) || (isnumeric(newVal) && isscalar(newVal) && ismember(newVal, [0 1])))
                obj.sendError('useFgen must be true or false')
            end
            obj.useFgen = logical(newVal);
            obj.changeFlag = true;      % the instrument has to change mode
        end

        function set.fgenRunTimeout(obj, newVal)
            if ~(isnumeric(newVal) && isscalar(newVal) && newVal > 0)
                obj.sendError('fgenRunTimeout must be a positive number of seconds')
            end
            obj.fgenRunTimeout = newVal;
        end

        function set.fgenSettleTime(obj, newVal)
            if ~(isnumeric(newVal) && isscalar(newVal) && newVal >= 0)
                obj.sendError('fgenSettleTime must be a non-negative number of seconds')
            end
            obj.fgenSettleTime = newVal;
        end

        function set.randomizeFrequencyOrder(obj, newVal)
            if ~(islogical(newVal) || (isnumeric(newVal) && isscalar(newVal) && ismember(newVal, [0 1])))
                obj.sendError('randomizeFrequencyOrder must be true or false')
            end
            obj.randomizeFrequencyOrder = logical(newVal);
            if ~obj.randomizeFrequencyOrder
                obj.sendWarning(['randomizeFrequencyOrder is OFF. Sequential visit order ' ...
                    'lets slow drift look like spectral structure -- use it for diagnostics ' ...
                    'only, and turn it back on before taking data you intend to keep.']);
            end
            % In walk mode the visit order IS the uploaded sequence, so this now
            % changes what is loaded. (In the classic path it does not, and the
            % one extra reload it causes there is harmless.)
            obj.changeFlag = true;
        end

        function set.newOrderEveryAverage(obj, newVal)
            if ~(islogical(newVal) || (isnumeric(newVal) && isscalar(newVal) && ismember(newVal, [0 1])))
                obj.sendError('newOrderEveryAverage must be true or false')
            end
            obj.newOrderEveryAverage = logical(newVal);
            % No changeFlag: the first upload is the same either way; the redraws
            % happen inside performWalk, one per average.
        end

        function set.fgenFullRestart(obj, newVal)
            if ~(islogical(newVal) || (isnumeric(newVal) && isscalar(newVal) && ismember(newVal, [0 1])))
                obj.sendError('fgenFullRestart must be true or false')
            end
            obj.fgenFullRestart = logical(newVal);
            % No changeFlag: this changes only what selectFgenFrequency sends per
            % point, not what is loaded on the instrument.
        end

        function set.burstTargetDuration(obj, newVal)
            if ~isempty(newVal) && ~(isnumeric(newVal) && isscalar(newVal) && newVal > 0)
                obj.sendError('burstTargetDuration must be a positive number of us, or [] for auto')
            end
            obj.burstTargetDuration = newVal;
            obj.changeFlag = true;      % changes every uploaded waveform
        end

        function set.walkCompactSteps(obj, newVal)
            if ~(islogical(newVal) || (isnumeric(newVal) && isscalar(newVal) && ismember(newVal, [0 1])))
                obj.sendError('walkCompactSteps must be true or false')
            end
            obj.walkCompactSteps = logical(newVal);
            obj.changeFlag = true;      % changes the sequence layout
        end

        function set.useSequenceJump(obj, newVal)
            if ~(islogical(newVal) || (isnumeric(newVal) && isscalar(newVal) && ismember(newVal, [0 1])))
                obj.sendError('useSequenceJump must be true or false')
            end
            obj.useSequenceJump = logical(newVal);
            obj.changeFlag = true;      % the AWG has to be reprogrammed either way
        end

        function set.mirrorSweepAround(obj, newVal)
            if ~isempty(newVal)
                % If it is empty, no further checking is needed
                checkFrequencyScalar(obj, newVal)
            end
            % If we got here, then newVal is OK.
            obj.mirrorSweepAround = newVal;
            obj.changeFlag = true;
        end

end

    methods
        function totalParamNum = getTotalNumberOfParams(obj)
            totalParamNum = length(obj.frequency);
        end
    end

    %% Helper functions
    methods
        function f = mirrorFrequency(obj)
            if isempty(obj.mirrorSweepAround)
                f = [];
            else
                % newFrequencies = |mirrorFreq - (oldFrequencies - mirrorFreq)|
                %                = |2 * mirrorFreq - oldFrequencies|
                % or, in proper code:
                f = abs(2*obj.mirrorSweepAround - obj.frequency);
                if min(f) < obj.mirrorSweepAround && max(f) > obj.mirrorSweepAround
                    f = flip(f);
                end
            end
        end

        function phase = calculatePhase(obj, f)
            phase = obj.timeDelay * f;
        end
    end

    %% Overridden from Experiment
    methods (Access = protected)
        function prepare(obj)
            % Initialize devices (SPCM, PulseGenerator, etc.)

            % Sequence
            %%% Useful parameters for what follows
            obj.isSingleMeasurement = isempty(obj.mirrorSweepAround);
            obj.detectionPeriodsPerRepeat = 2;
            obj.runsPerPerform = (1 + ~obj.isSingleMeasurement);
            obj.freqMirrored = obj.mirrorFrequency;
            MWChannel = obj.MWChannel;

            %%% Create
            if obj.seperatedTrigger                         % AWG trigger
                triggers = {'trigger', 'trigger2'};
            else
                triggers = {'trigger'};
            end
            switch obj.nChannels
                case 1
                    MW = {MWChannel};
                case 2
                    MW = {'MW', 'MW2'};
                otherwise
                    obj.sendError('What should we do here?')
            end
            if obj.AWGorSRSswitch
                MW = [MW{:}, {'AWG'}];
            end
            if ~obj.driveMWSwitches
                % Diagnostic: never tell the switches to open. See driveMWSwitches.
                MW = {};
                fprintf(['MW SWITCHES NOT DRIVEN (driveMWSwitches = false). If the spectrum ' ...
                         'is unchanged,\n  the switches are not gating and the AWG''s zeros ' ...
                         'are the only gate.\n']);
            end

            S = Sequence;
%             S.addEvent(obj.AWG.TRIG_DURATION,               [{'greenLaser'}, triggers(:)', MW(:)']);
%             S.addEvent(obj.laserInitializationDuration,     [MW{:}, {'greenLaser'}]);
%             S.addEvent(obj.detectionDuration,               [MW{:}, {'greenLaser'}, {'detector'}]);
%             S.addEvent(obj.laserInitializationDuration,     {'greenLaser'});
%             S.addEvent(obj.referenceDetectionDuration,      {'greenLaser','detector'})
%
            % The only difference from AWGExpESR: each delayBetweenDetections gap is
            % split into four events carrying two marker pulses (see addGap). The total
            % length of the repeat is unchanged, so the AWG upload is unchanged too.
            obj.registerGapAomChannel();
            S.addEvent(obj.AWG.TRIG_DURATION,                           [{'greenLaser'}, triggers(:)', MW(:)']);
            S.addEvent(obj.detectionDuration - obj.AWG.TRIG_DURATION,   [{'greenLaser'}, MW{:}, {'detector'}]);
            obj.addGap(S);
            S.addEvent(obj.AWG.TRIG_DURATION,                           {'greenLaser'});
            S.addEvent(obj.detectionDuration - obj.AWG.TRIG_DURATION,   {'greenLaser','detector'})
            obj.addGap(S);

            if obj.changeFlag % || obj.restartFlag
                obj.LoadAWG;
            elseif obj.useFgen
                % Cheap and idempotent, and it has to hold on EVERY run: a run that
                % skips LoadAWG would otherwise inherit whatever mode the last
                % experiment left behind.
                obj.SetupFgen;
            end

            % Re-enable the AWG outputs on EVERY run, not just the ones that
            % reload. Experiment.wrapUpInternal calls fg.disconnect at the end of
            % every run, and disconnect sends enableOutput 0. That command used to
            % be malformed ('AWGC:STOP; OUTP1 OFF' -- a bare ';' keeps the parser
            % in the AWGControl subtree, so the instrument rejected it with -113)
            % and therefore never did anything. Now that it works, something has
            % to turn the outputs back on, and prepareInternal does not: its
            % fg.output = true branch is skipped for a TektronixAWG by
            %   (~contains(name,'TektronixAWG') + ~contains(name,'synthhd...') == 2)
            % which is 0 + 1 = 1, never 2.
            %
            % Without this line, any run that skips LoadAWG -- a.restart, or any
            % run with changeFlag false -- emits no MW at all and every point
            % comes back at ratio ~= 1.
            %
            % Two SCPI commands, idempotent.
            obj.AWG.workChannels = 1:obj.nChannels;

            obj.prepareInternal(S)

            if obj.useSequenceWalk && obj.useFgen
                obj.sendError('useFgen and useSequenceWalk are mutually exclusive -- pick one')
            end
            if obj.useSequenceWalk
                % prepareInternal sizes the pulse generator and digitizer for
                % "repeats triggers at ONE frequency". Walk mode instead wants
                % "one trigger per frequency, passesPerRun times over", so the
                % sequencer advances a frequency per trigger and a single pg.run
                % returns a whole pass. Override those two numbers after the fact
                % rather than duplicating prepareInternal.
                pg = getObjByName(PulseGenerator.NAME);
                spcm = getObjByName(Spcm.NAME);
                nF = length(obj.frequency);
                % Each pass fires walkWarmupTriggers throwaway blocks before the
                % sweep proper, so that frequency 1 is not the trigger that
                % follows the arm/idle gap. performWalk discards those reads.
                blocksPerPass = nF + obj.walkWarmupTriggers;
                triggersPerRun = blocksPerPass * obj.passesPerRun;
                pg.repeats = triggersPerRun;
                numScans = obj.detectionPeriodsPerRepeat * triggersPerRun;
                seqTime = pg.sequenceDuration() * 1e-6;
                spcm.prepareExperimentCount(numScans, 2 * seqTime * triggersPerRun, 1);
                fprintf(['Walk mode: %d triggers per pg.run (%d frequencies + %d warm-up, ' ...
                         'x %d passes), %d digitizer segments.\n'], ...
                         triggersPerRun, nF, obj.walkWarmupTriggers, obj.passesPerRun, numScans);
            end

            % Set parameter, for saving
            obj.mCurrentXAxisParam.value = obj.frequency;

            % Gap power: free any task a crashed run left behind (restoring its window
            % value), then take the window value from the panel as it is NOW.
            obj.releaseGapAom();
            obj.captureGapAomValues();
        end

        function ensureInstrumentMode(obj, wanted)
            % Put the AWG in 'AWG' or 'FGEN' mode, but only if it is not already
            % there -- the query is ~2 ms on visadev and a needless mode switch
            % stops playback. Switching invalidates the upload cache, because the
            % assets are no longer being played and the cheap thing to do is
            % reload rather than reason about what survived.
            current = strtrim(obj.AWG.queryRaw('INST:MODE?'));
            if strcmpi(current, wanted)
                return
            end
            fprintf('AWG instrument mode: %s -> %s\n', current, wanted);
            obj.AWG.setValue('mode', wanted);
            obj.AWG.invalidateLoaded;
        end

        function SetupFgen(obj)
            % Configure the built-in function generator. No uploads at all.
            if ~obj.isSingleMeasurement
                obj.sendError('useFgen does not support mirrorSweepAround yet')
            end

            obj.ensureInstrumentMode('FGEN');

            % Reject a timeDelay that phase cannot express, rather than wrapping it
            % silently. Worst case over the sweep is the highest frequency.
            if obj.timeDelay ~= 0
                worstPhase = 360 * max(obj.frequency) * abs(obj.timeDelay);   % MHz * us = cycles
                if worstPhase > 180
                    obj.sendError(sprintf(['timeDelay = %g us is %.1f deg at %g MHz, past the ' ...
                        '+/-180 deg the FGEN phase control allows (+/-%.1f ps at that frequency). ' ...
                        'Use AWG mode for a delay this large.'], obj.timeDelay, worstPhase, ...
                        max(obj.frequency), 0.5/max(obj.frequency)*1e6))
                end
            end

            cmds = {'FGEN:CHAN1:TYPE SINE'};
            if obj.nChannels == 2
                cmds{end+1} = 'FGEN:CHAN2:TYPE SINE';
            end
            % Amplitude in FGEN is dBm (AMPL:POW), the same units as obj.amplitude,
            % so the calibrated per-channel difference carries over directly and
            % does NOT go through MaxNormPowerAllowed.
            ampl = obj.amplitude;
            if isscalar(ampl); ampl = ampl * ones(1, obj.nChannels); end
            for chan = 1:obj.nChannels
                cmds{end+1} = sprintf('FGEN:CHAN%d:AMPL:POW %g', chan, ampl(chan));  %#ok<AGROW>
            end
            % Set a starting frequency before starting the generator, so this
            % matches the sequence awgFgenBringUp actually validated.
            cmds{end+1} = sprintf('FGEN:CHAN1:FREQ %.6f', obj.frequency(1)*1e6);
            obj.AWG.sendCommandBatch(cmds);

            % START THE GENERATOR *BEFORE* ENABLING THE OUTPUTS, then WAIT.
            %
            % Both halves of that matter, and getting either wrong is why FGEN
            % measured nothing for four rounds. awgFgenBringUp tried four orderings:
            %
            %   outputs on -> RUN:IMMEDIATE   run state 0 at 0.0, 0.2 and 1.2 s
            %   RUN:IMMEDIATE -> outputs on   0, 0, then 2
            %   OUTP1 ON;:AWGCONTROL:RUN:IMM  0, 0, then 2
            %   outputs on -> RUN             0, 0, then 2
            %
            % The first is exactly what this method used to do, and it is the only
            % one that never runs. On top of that, even the working orderings read 0
            % until about 1.2 s, so checking the state immediately after the command
            % reports "not running" whether it is or not. No SCPI error is raised in
            % either case -- the error queue was clean throughout.
            obj.AWG.sendCommand('awgcontrol:run:immediate');
            runState = obj.waitForFgenRunning();
            obj.AWG.workChannels = 1:obj.nChannels;     % outputs on, after the run

            outState = strtrim(obj.AWG.queryRaw('OUTPUT1:STATE?'));
            freqBack = str2double(obj.AWG.queryRaw('FGEN:CHAN1:FREQ?'))/1e6;
            amplBack = str2double(obj.AWG.queryRaw('FGEN:CHAN1:AMPL:POW?'));
            fprintf(['FGEN mode: sine, amplitude %s dBm, timeDelay %g us as phase.\n' ...
                     '  readback: output1 = %s, run state = %s, %.3f MHz, %.2f dBm\n'], ...
                    mat2str(ampl), obj.timeDelay, outState, runState, freqBack, amplBack);
            if ~any(strcmp(outState, {'1', 'ON'}))
                obj.sendWarning(['FGEN: output 1 reads back as "' outState '" rather than on. ' ...
                    'With no microwave there is nothing to measure.']);
            end
        end

        function runState = waitForFgenRunning(obj)
            % Poll until the generator reports a non-zero run state, or give up.
            % Non-zero rather than a specific value: 2 is what this instrument
            % returns in FGEN mode, 1 means waiting for trigger, and either counts
            % as started. Only 0 -- stopped -- is a failure.
            t0 = tic;
            runState = '0';
            while toc(t0) < obj.fgenRunTimeout
                runState = strtrim(obj.AWG.queryRaw('AWGCONTROL:RSTATE?'));
                if ~strcmp(runState, '0')
                    return
                end
                pause(0.05);
            end
            obj.sendError(sprintf(['FGEN did not start: run state still 0 after %g s. ' ...
                'Nothing will be measured in this state -- every point would come back ' ...
                'at ratio 1. Run awgFgenBringUp to see which start ordering this ' ...
                'firmware accepts.'], obj.fgenRunTimeout));
        end

        function selectFgenFrequency(obj, freqMHz)
            % Move the sweep to one frequency: one write for the tone, plus the
            % phase that carries timeDelay. Phase for a FIXED time delay depends on
            % frequency (delay = phase / 360 f), so it is recomputed at every point
            % rather than set once.
            %
            % A BARE FREQUENCY WRITE DOES NOT GIVE A STABLE TONE ON THIS INSTRUMENT.
            % fgenStepTime measured it directly at the scope: after writing only
            % FGEN:CHANn:FREQuency, the amplitude at the NEW frequency sits at the noise
            % floor for about 2.5 s, appears briefly, and does not persist -- out to the
            % full 5 s watched, at every step size from 2 to 800 MHz. Repeating the whole
            % bring-up (RUN:IMMEDIATE, poll the run state, re-enable the outputs) gives a
            % tone in ~2.8 s that is stable by ~4.1 s.
            %
            % That is almost certainly why FGEN ESR measured 0.039% against AWG mode's
            % 3.9%: with fgenSettleTime at its 1 ms default, every point was read with no
            % microwave present at the frequency being swept.
            %
            % obj.fgenFullRestart switches this method between the two. It defaults to
            % FALSE, which preserves the original behaviour exactly, because the honest
            % summary is that neither setting makes FGEN worth using for a sweep: at
            % ~4.5 s per point, 96 points is over seven minutes of settling against ~30 s
            % for the entire working AWG-mode run. It is here so the contrast question can
            % be closed by measurement rather than left as a suspicion.
            cmds = {};
            for chan = 1:obj.nChannels
                % Frequency is coupled across channels on this instrument, so
                % writing both is redundant but makes the result independent of
                % which channel was written last.
                cmds{end+1} = sprintf('FGEN:CHAN%d:FREQ %.6f', chan, freqMHz*1e6);  %#ok<AGROW>
            end
            for chan = 1:obj.nChannels
                cmds{end+1} = sprintf('FGEN:CHAN%d:PHAS %.4f', chan, ...
                                      obj.fgenPhaseFor(chan, freqMHz));             %#ok<AGROW>
            end
            obj.AWG.sendCommandBatch(cmds);

            if obj.fgenFullRestart
                % Same ordering SetupFgen validated: run first, poll, THEN outputs.
                % Enabling outputs before the run fails silently on this firmware.
                obj.AWG.sendCommand('awgcontrol:run:immediate');
                obj.waitForFgenRunning();
                obj.AWG.workChannels = 1:obj.nChannels;
            end

            if obj.fgenSettleTime > 0
                pause(obj.fgenSettleTime);
            end
        end

        function phaseDeg = fgenPhaseFor(obj, chan, freqMHz)
            % Which channel gets the delay follows the AWG-mode convention exactly:
            % a positive timeDelay delays channel 1, a negative one delays channel 2
            % (see LoadAWG, where the same test picks which leading zero is padded).
            phaseDeg = 0;
            if obj.timeDelay > 0 && chan == 1
                phaseDeg = 360 * freqMHz * obj.timeDelay;        % MHz * us = cycles
            elseif obj.timeDelay < 0 && chan == 2
                phaseDeg = 360 * freqMHz * abs(obj.timeDelay);
            end
        end

        function [segmentTime, nCycles, nSamples] = burstSegmentTime(obj, freqMHz)
            % Length of the short looped waveform for one frequency, in us.
            %
            % Must hold a WHOLE number of MW cycles, or the loop seam is a phase
            % discontinuity repeated numSegments times. So the length is n/f for
            % integer n, and the only question is how small n can be.
            %
            % The binding constraint is the instrument's minimum waveform length
            % (MIN_WAVEFORM_LENGTH = 2400 samples), not the physics. Auto takes
            % the smallest n that clears it with 5% margin.
            S = obj.AWG.SampleRate;                 % MHz, so samples = us * S
            if isempty(obj.burstTargetDuration)
                minSamples = 1.05 * obj.AWG.MIN_WAVEFORM_LENGTH;
                nCycles = obj.seamlessCycles(freqMHz, minSamples, 4000);
            else
                nCycles = round(freqMHz * obj.burstTargetDuration);
            end
            nCycles = max(nCycles, 1);
            segmentTime = nCycles / freqMHz;
            nSamples = segmentTime * S;
            if nSamples < obj.AWG.MIN_WAVEFORM_LENGTH
                obj.sendError(sprintf(['burstTargetDuration = %g us gives %.0f samples at ' ...
                    '%g MS/s, under the AWG minimum of %d. Raise it to at least %.4f us ' ...
                    'or leave it empty for auto.'], obj.burstTargetDuration, nSamples, S, ...
                    obj.AWG.MIN_WAVEFORM_LENGTH, obj.AWG.MIN_WAVEFORM_LENGTH/S))
            end
        end

        function fp = awgLoadFingerprint(obj)
            % Hash of everything that determines the waveform/sequence content.
            % If none of these changed, the AWG already holds the right data and
            % re-uploading it is pure waste (17.8 s on a 191-point sweep).
            % Deliberately does NOT include repeats or averages -- those live in
            % the pulse generator and the software loop, not in the AWG.
            % EVERY setting that changes what gets uploaded must be in here.
            % useSequenceWalk was missing once, and the result was a false cache
            % hit: a classic run stamped the fingerprint, the walk run computed
            % the same hash, skipped the upload, and spent the whole sweep playing
            % one classic single-frequency sequence. That produced a flat spectrum
            % which looked like a sequencer failure and was not one.
            fp = DataHash({obj.frequency, obj.amplitude, obj.detectionDuration, ...
                           obj.laserInitializationDuration, obj.nChannels, ...
                           obj.timeDelay, obj.channelPhase, obj.AWG.SampleRate, obj.mirrorSweepAround, ...
                           obj.useSequenceJump, obj.useSequenceWalk, ...
                           obj.walkWarmupTriggers, obj.walkCompactSteps, ...
                           obj.burstTargetDuration, ...
                           obj.MWChannel, obj.AWGorSRSswitch, ...
                           obj.randomizeFrequencyOrder});   % walk: shuffled vs ascending sequence
        end

        function n = expectedSequenceCount(obj)
            % How many sequences the AWG should be holding for the current mode.
            % Cheap independent check that a cache hit is telling the truth.
            if obj.useSequenceJump || obj.useSequenceWalk
                n = obj.nChannels;                              % one big concatenated sequence each
            else
                n = length(obj.frequency) * obj.nChannels;      % one per (frequency, channel)
            end
        end

        function LoadAWG(obj)
            if obj.useFgen
                if obj.useSequenceJump || obj.useSequenceWalk
                    obj.sendError(['useFgen cannot be combined with useSequenceJump or ' ...
                                   'useSequenceWalk -- FGEN mode has no sequencer at all'])
                end
                obj.SetupFgen();
                return
            end

            % Not FGEN: make sure the instrument is actually in AWG mode. Nothing
            % used to assert this, because nothing used to leave it in FGEN. Now
            % that an ESR run can, a following AWG-mode run would otherwise upload
            % a sequence to an instrument that cannot play one.
            obj.ensureInstrumentMode('AWG');

            if obj.useSequenceJump && obj.useSequenceWalk
                obj.sendError('useSequenceJump and useSequenceWalk are mutually exclusive -- pick one')
            end
            concatenated = obj.useSequenceJump || obj.useSequenceWalk;

            % --- skip the upload if the AWG already holds this exact set ------
            fp = obj.awgLoadFingerprint();
            if obj.AWG.isAlreadyLoaded(fp)
                % Verify against the instrument rather than trusting the hash on
                % its own. One query, ~2.3 ms on visadev, and it catches a stale
                % or wrong-mode cache before it can silently produce garbage data.
                nLoaded = str2double(obj.AWG.queryRaw('SLISt:SIZE?', 10));
                nWanted = obj.expectedSequenceCount();
                if nLoaded ~= nWanted
                    obj.sendWarning(sprintf(['AWG holds %d sequences but this mode needs %d -- ' ...
                        'the load cache was stale. Reloading.'], nLoaded, nWanted));
                    obj.AWG.invalidateLoaded();
                else
                    cached = obj.AWG.getLoadedCache();
                    obj.indexSeq = cached.indexSeq;
                    if concatenated
                        obj.seqAssetIndex = cached.seqAssetIndex;
                        % The order the loaded sequence plays in, and what is needed to
                        % rebuild it. An older cache has neither: its order was ascending.
                        obj.walkOrder = 1:length(obj.frequency);
                        obj.walkParts = [];
                        if isfield(cached, 'walkOrder') && ~isempty(cached.walkOrder)
                            obj.walkOrder = cached.walkOrder;
                        end
                        if isfield(cached, 'walkParts'); obj.walkParts = cached.walkParts; end
                        obj.walkOrderFresh = false;     % it was used by an earlier run
                    end
                    % Amplitude and channel enables live in the instrument's front
                    % end, not in the waveform data, so re-apply them cheaply rather
                    % than assume they survived.
                    obj.AWG.workChannels = 1:obj.nChannels;
                    obj.AWG.setAmplitude(obj.amplitude);
                    if concatenated
                        % Both concatenated modes need their one big sequence
                        % assigned and the instrument armed -- the sweep never
                        % assigns anything again.
                        for chan = 1:obj.nChannels
                            obj.AWG.assignSequence(obj.seqAssetIndex(chan), chan);
                        end
                        obj.AWG.Run;
                    end
                    fprintf('AWG already holds this waveform set (%d sequences) - skipping upload.\n', nLoaded);
                    return
                end
            end

            if concatenated
                % perBlockGoto TRUE  -> each block loops on itself, addressed by a
                %                       forced jump (the jump path, unreliable here)
                % perBlockGoto FALSE -> blocks fall through to the next one and the
                %                       last step returns to the first, so the
                %                       sequencer WALKS the sweep on triggers
                obj.LoadAWGConcatenated(obj.useSequenceJump);
                obj.markWalkLoaded(fp);
                return
            end
            obj.AWG.ResetStoredData;
            obj.AWG.workChannels = 1:obj.nChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            minWaveformTime = obj.AWG.minWaveformDuration;
            N = length(obj.frequency);
            WaveformsIndex = zeros(obj.nChannels, N, obj.runsPerPerform);
            obj.indexSeq = zeros(obj.nChannels, N);

            % creating waveforms:
            A = obj.AWG.MaxNormPowerAllowed;
            % Phase in radians, per channel. See the channelPhase property.
            phi = zeros(1, obj.nChannels);
            cp  = obj.channelPhase(:)';
            phi(1:min(numel(cp), obj.nChannels)) = ...
                deg2rad(cp(1:min(numel(cp), obj.nChannels)));
            w = 2*pi*obj.frequency;
            numSegments = zeros(obj.runsPerPerform,N);
            % THE OLD COUNTS WERE BOTH WRONG. The first said N*runsPerPerform*nChannels, but
            % the loop below adds one waveform per frequency per channel -- runsPerPerform has
            % nothing to do with it. The second said 2*N*nChannels for a loop that adds
            % N*nChannels, and adds NOTHING AT ALL when isSingleMeasurement is true, so it
            % announced a load that never happened. Neither counted the zero pads.
            %
            % One session for the whole load, with the true total: the mirrored set is only
            % built when there is a reference measurement, and there are nChannels + 1 zero
            % pads at the end.
            nMirror = 0;
            if ~obj.isSingleMeasurement; nMirror = N*obj.nChannels; end
            obj.AWG.progressBegin(sprintf(['ESR waveforms (classic, one sequence per point): ' ...
                '%d frequency burst(s) x %d channel(s)%s, plus %d zero pads. Each burst is a ' ...
                'whole number of RF cycles so it can be repeated seamlessly to fill the ' ...
                'detection window.'], N, obj.nChannels, ...
                BooleanHelper.ifTrueElse(obj.isSingleMeasurement, '', ...
                    sprintf(' + %d mirrored reference', nMirror)), obj.nChannels + 1), ...
                N*obj.nChannels + nMirror + obj.nChannels + 1, 'waveform');
            for i = 1:N
                segmentTime = obj.burstSegmentTime(obj.frequency(i));
                numSegments(1,i) = round((obj.detectionDuration+obj.laserInitializationDuration) / segmentTime);
                for chan = 1:obj.nChannels
                    f = @(t) A(chan) * (0<=t & t<segmentTime) .* cos(w(i)*t - phi(chan));
                    WaveformsIndex(chan,i,1) = obj.AWG.AddWaveformByFunction(f, segmentTime, 0, 1);
                end
            end
            if ~obj.isSingleMeasurement
                w = 2*pi*obj.freqMirrored;
                for i = 1:N
                    segmentTime = obj.burstSegmentTime(obj.freqMirrored(i));
                    numSegments(2,i) = round(obj.detectionDuration/segmentTime);
                    for chan = 1:obj.nChannels
                        f = @(t) A(chan) * (0<=t & t<segmentTime) .* cos(w(i)*t - phi(chan));
                        WaveformsIndex(chan,i,2) = obj.AWG.AddWaveformByFunction(f, segmentTime, 0, 1);
                    end
                end
            end
            f = @(t) 0*t;
            firstZeroIndex = zeros(obj.nChannels,1);
            for chan = 1:obj.nChannels
                if obj.timeDelay > 0
                    prePulseTime = minWaveformTime + obj.timeDelay*(chan==1);
                else
                    prePulseTime = minWaveformTime + obj.timeDelay*(chan==2);
                end
                firstZeroIndex(chan) = obj.AWG.AddWaveformByFunction(f, prePulseTime, 0, 1);
            end
            lastZeroIndex  = obj.AWG.AddWaveformByFunction(f, minWaveformTime, 0, 1);
            obj.AWG.progressEnd();
            obj.AWG.LoadWaveform();

            % creating sequences:
            obj.AWG.progressBegin(sprintf(['ESR sequences: one per frequency per channel, ' ...
                '%d steps each.'], BooleanHelper.ifTrueElse(obj.isSingleMeasurement, 3, 6)), ...
                N*obj.nChannels, 'sequence');
            for i = 1:N
                for chan = 1:obj.nChannels
                    if obj.isSingleMeasurement
                        obj.indexSeq(chan,i) = obj.AWG.AddSequence([firstZeroIndex(chan) WaveformsIndex(chan,i,1) lastZeroIndex], ...
                                                                   [1 numSegments(1,i) 1], [1 0 0]);
                    else
                        obj.indexSeq(chan,i) = obj.AWG.AddSequence([firstZeroIndex(chan) WaveformsIndex(chan,i,1) lastZeroIndex, ...
                                                                    firstZeroIndex(chan) WaveformsIndex(chan,i,2) lastZeroIndex], ...
                                                                   [1 numSegments(1,i) 1, 1 numSegments(2,i) 1], [1 0 0, 1 0 0], 0);
                    end
                    obj.AWG.progressStep();
                end
            end
            obj.AWG.progressEnd();
            obj.AWG.LoadSequence()
            obj.AWG.markLoaded(fp, struct('indexSeq', obj.indexSeq, ...
                                          'seqAssetIndex', obj.seqAssetIndex));
        end

        function LoadAWGConcatenated(obj, perBlockGoto)
            % ONE sequence per channel holding every frequency, so the sweep can
            % move between frequencies with a forced jump instead of an asset
            % reassign. See the useSequenceJump property for why.
            %
            % Layout, per channel: frequency i occupies steps 3i-2 .. 3i
            %   3i-2  leading zero,  rcount 1,              waits for trigger
            %   3i-1  MW burst,      rcount numSegments(i), no trigger
            %   3i    trailing zero, rcount 1,              GOTO -> 3i-2
            % so once parked on frequency i the sequencer replays that block on
            % every trigger and never wanders to another frequency by itself.

            if ~obj.isSingleMeasurement
                obj.sendError('Neither useSequenceJump nor useSequenceWalk supports mirrorSweepAround yet')
            end
            if obj.repeats <= 10
                obj.sendWarning(sprintf(['repeats = %d: a mid-burst trigger boundary can corrupt the first ' ...
                    'repeat, and only repeats > 10 discards it. ' ...
                    'Raise repeats above 10.'], obj.repeats))
            end

            obj.AWG.ResetStoredData;
            obj.AWG.workChannels = 1:obj.nChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            minWaveformTime = obj.AWG.minWaveformDuration;
            N = length(obj.frequency);
            A = obj.AWG.MaxNormPowerAllowed;
            % Phase in radians, per channel. See the channelPhase property.
            phi = zeros(1, obj.nChannels);
            cp  = obj.channelPhase(:)';
            phi(1:min(numel(cp), obj.nChannels)) = ...
                deg2rad(cp(1:min(numel(cp), obj.nChannels)));
            w = 2*pi*obj.frequency;

            % --- block layout ------------------------------------------------
            % Jump mode is fixed at three steps per block: its trailing zero is
            % what carries the self-referencing GOTO. Walk mode can drop that
            % trailing zero (see walkCompactSteps) and prepends warm-up blocks
            % (see walkWarmupTriggers).
            if perBlockGoto
                stepsPerBlock = 3;
                nWarmup = 0;
            else
                stepsPerBlock = BooleanHelper.ifTrueElse(obj.walkCompactSteps, 2, 3);
                nWarmup = obj.walkWarmupTriggers;
            end
            nSteps  = stepsPerBlock * (N + nWarmup);
            if nSteps > obj.AWG.MAX_SEQUENCE_STEPS
                obj.sendError(sprintf(['%d frequencies need %d sequence steps, over the %d step ' ...
                    'limit. Split the sweep or turn the concatenated modes off.'], ...
                    N, nSteps, obj.AWG.MAX_SEQUENCE_STEPS))
            end

            % --- waveforms: identical to the classic path -------------------
            % Timed in three separate phases, because the two settle times the
            % loader already prints do NOT include the waveform generation or the
            % binary transfer. Inferring that share by subtracting from the total
            % is how I mis-attributed a cold load as 73 s when it was ~130 s.
            tWave = tic;
            burstIndex  = zeros(obj.nChannels, N);
            numSegments = zeros(1, N);
            [~, ~, nSamp1] = obj.burstSegmentTime(obj.frequency(1));
            fprintf('Burst waveform: %.0f samples each (%.4f us at %g MS/s), %d of them.\n', ...
                    nSamp1, nSamp1/obj.AWG.SampleRate, obj.AWG.SampleRate, N*obj.nChannels);
            % The zero pads were never counted here either: there are nChannels + 1 of them.
            obj.AWG.progressBegin(sprintf(['ESR waveforms (concatenated, ONE sequence per ' ...
                'channel): %d frequency burst(s) x %d channel(s), plus %d zero pads. The ' ...
                'sweep lives in the SEQUENCE, not in separate sequences, which is what lets ' ...
                'the AWG be armed once for the whole scan.'], N, obj.nChannels, ...
                obj.nChannels + 1), N*obj.nChannels + obj.nChannels + 1, 'waveform');
            for i = 1:N
                segmentTime = obj.burstSegmentTime(obj.frequency(i));
                numSegments(i) = round((obj.detectionDuration + obj.laserInitializationDuration) / segmentTime);
                for chan = 1:obj.nChannels
                    f = @(t) A(chan) * (0<=t & t<segmentTime) .* cos(w(i)*t - phi(chan));
                    burstIndex(chan,i) = obj.AWG.AddWaveformByFunction(f, segmentTime, 0, 1);
                end
            end
            fZero = @(t) 0*t;
            firstZeroIndex = zeros(obj.nChannels,1);
            for chan = 1:obj.nChannels
                if obj.timeDelay > 0
                    prePulseTime = minWaveformTime + obj.timeDelay*(chan==1);
                else
                    prePulseTime = minWaveformTime + obj.timeDelay*(chan==2);
                end
                firstZeroIndex(chan) = obj.AWG.AddWaveformByFunction(fZero, prePulseTime, 0, 1);
            end
            lastZeroIndex = obj.AWG.AddWaveformByFunction(fZero, minWaveformTime, 0, 1);
            obj.AWG.progressEnd();
            obj.AWG.LoadWaveform();
            tWaveform = toc(tWave);

            % --- the sequence, in the visit order ------------------------------
            % Kept apart from the waveforms so newOrderEveryAverage can rebuild it
            % in a new order without re-uploading a single waveform.
            obj.walkParts = struct('burstIndex', burstIndex, 'numSegments', numSegments, ...
                'firstZeroIndex', firstZeroIndex, 'lastZeroIndex', lastZeroIndex, ...
                'stepsPerBlock', stepsPerBlock, 'nWarmup', nWarmup, 'perBlockGoto', perBlockGoto);
            [tSequence, tArm] = obj.uploadWalkSequence(obj.drawWalkOrder());
            fprintf(['AWG armed once for %d frequencies + %d warm-up block(s), ' ...
                     '%d steps/channel (%d per block), %s mode, %s order.\n'], ...
                    N, nWarmup, nSteps, stepsPerBlock, ...
                    BooleanHelper.ifTrueElse(perBlockGoto, 'jump', 'walk'), ...
                    BooleanHelper.ifTrueElse(isequal(obj.walkOrder, 1:N), 'ascending', 'random'));
            fprintf(['Cold load breakdown: waveforms %.1f s, sequence %.1f s, ' ...
                     'first arm %.1f s  ->  %.1f s total.\n'], ...
                    tWaveform, tSequence, tArm, tWaveform + tSequence + tArm);
        end

        function order = drawWalkOrder(obj)
            % Visit order for the concatenated sequence: random when
            % randomizeFrequencyOrder is on, ascending otherwise. Jump mode
            % addresses blocks by step number and keeps the ascending order.
            N = length(obj.frequency);
            if obj.randomizeFrequencyOrder && ~obj.useSequenceJump
                order = randperm(N);
            else
                order = 1:N;
            end
        end

        function [tSequence, tArm] = uploadWalkSequence(obj, order)
            % Build ONE sequence per channel that plays the frequencies in 'order'
            % (order(k) = the frequency index of the k-th block), upload it, assign
            % it and arm. Only the sequence list is replaced: the waveforms named in
            % obj.walkParts must already be on the instrument.
            p = obj.walkParts;
            N = length(obj.frequency);
            blockFreq = [repmat(order(1), 1, p.nWarmup), order];  % warm-up copies the first block
            nBlocks = numel(blockFreq);
            nSteps  = p.stepsPerBlock * nBlocks;

            tSeq = tic;
            % Nothing may still reference the sequence that is about to be deleted.
            obj.AWG.sendCommand('awgcontrol:stop:immediate');
            obj.AWG.sequences = struct('waveforms', {{}}, 'repeats', {{}}, ...
                                       'trigger', {{}}, 'goto', {{}});   % sequences only
            obj.indexSeq      = zeros(obj.nChannels, N);    % STEP numbers here
            obj.seqAssetIndex = zeros(1, obj.nChannels);
            for chan = 1:obj.nChannels
                wfList   = zeros(1, nSteps);
                rptList  = zeros(1, nSteps);
                trgList  = zeros(1, nSteps);
                gotoList = zeros(1, nSteps);
                for blk = 1:nBlocks
                    i  = blockFreq(blk);
                    s0 = p.stepsPerBlock*(blk-1);
                    if p.stepsPerBlock == 3
                        wfList(s0+1:s0+3)  = [p.firstZeroIndex(chan), p.burstIndex(chan,i), p.lastZeroIndex];
                        rptList(s0+1:s0+3) = [1, p.numSegments(i), 1];
                        trgList(s0+1:s0+3) = [1, 0, 0];
                        if p.perBlockGoto
                            gotoList(s0+3) = s0+1;      % block loops on itself (jump mode)
                        end
                    else
                        wfList(s0+1:s0+2)  = [p.firstZeroIndex(chan), p.burstIndex(chan,i)];
                        rptList(s0+1:s0+2) = [1, p.numSegments(i)];
                        trgList(s0+1:s0+2) = [1, 0];
                    end
                    if blk > p.nWarmup
                        obj.indexSeq(chan,i) = s0+1;    % jump here to select frequency i
                    end
                end
                obj.seqAssetIndex(chan) = obj.AWG.AddSequence(wfList, rptList, trgList, 0, gotoList);
            end
            obj.AWG.LoadSequence();
            tSequence = toc(tSeq);

            % --- assign and arm ONCE for the whole sweep ---------------------
            tArmStart = tic;
            for chan = 1:obj.nChannels
                obj.AWG.assignSequence(obj.seqAssetIndex(chan), chan);
            end
            obj.AWG.Run;
            tArm = toc(tArmStart);

            obj.walkOrder = order;
            obj.walkOrderFresh = true;      % not used by an average yet
        end

        function markWalkLoaded(obj, fp)
            % Stamp the load cache with everything a later cache hit needs,
            % including the order the loaded sequence plays in.
            obj.AWG.markLoaded(fp, struct('indexSeq', obj.indexSeq, ...
                'seqAssetIndex', obj.seqAssetIndex, 'walkOrder', obj.walkOrder, ...
                'walkParts', obj.walkParts));
        end

        function updateSignalParam(obj)
            % Collapse obj.signal/obj.sterr across completed averages into the
            % plotted normalized signal. Shared by the classic and walk paths so
            % they cannot drift apart.
            S1 = squeeze(obj.signal(1, :, 1:obj.currIter));
            S1sterr = squeeze(obj.sterr(1, :, 1:obj.currIter));
            S2 = squeeze(obj.signal(2, :, 1:obj.currIter));
            S2sterr = squeeze(obj.sterr(2, :, 1:obj.currIter));

            [value, sterr] = getRatioDistributionValues(obj, S1, S2, S1sterr, S2sterr);
            obj.signalParam.value = value;
            obj.signalParam.sterr = sterr;

            if obj.isSingleMeasurement
                obj.signalParam2.value = [];
            else
                S3 = squeeze(obj.signal(3, :, 1:obj.currIter));
                S3sterr = squeeze(obj.sterr(3, :, 1:obj.currIter));
                S4 = squeeze(obj.signal(4, :, 1:obj.currIter));
                S4sterr = squeeze(obj.sterr(4, :, 1:obj.currIter));

                [value, sterr] = getRatioDistributionValues(obj, S3, S4, S3sterr, S4sterr);
                obj.signalParam2.value = value;
                obj.signalParam2.sterr = sterr;
            end
        end

        function performWalk(obj)
            % One average, acquired as whole passes over the sweep.
            %
            % The AWG holds one sequence with every frequency as a consecutive
            % block, and advances one block per trigger. So a pg.run of
            % nFreq*passesPerRun triggers returns
            %   [s1 r1 s2 r2 ... sN rN]  repeated passesPerRun times
            % and there is no per-point AWG traffic and no per-point re-arm at all.
            nF = length(obj.frequency);
            pg = getObjByName(PulseGenerator.NAME);
            spcm = getObjByName(Spcm.NAME);
            tracker = getObjByName(Tracker.NAME);
            if isempty(tracker); throwBaseObjException(Tracker.Name); end

            passesWanted = obj.repeats;      % passes here play the role of repeats
            nRuns = ceil(passesWanted / obj.passesPerRun);
            % Each pass plays walkWarmupTriggers throwaway blocks first, so a pass
            % is nF+nWarmup blocks long and the leading nWarmup read pairs are cut.
            % This is what makes frequency 1 comparable to the other 95: the
            % classic path throws away the first repeat AT EVERY frequency, so
            % nothing it keeps is ever the trigger straight after the idle gap.
            nWarmup = obj.walkWarmupTriggers;
            blocksPerPass = nF + nWarmup;
            readsPerRun = obj.detectionPeriodsPerRepeat * blocksPerPass * obj.passesPerRun;

            raw = nan(obj.detectionPeriodsPerRepeat, nF, nRuns * obj.passesPerRun);
            passesDone = 0;
            failedPasses = 0;       % counted and reported: a pass that had to be
                                    % retried aborts mid-acquisition, and that is
                                    % the leading suspect for walk data that
                                    % correlates ~0.8 with the classic spectrum
                                    % instead of matching it.

            opts = obj.walkReduceOptions(spcm);

            % ---- visit order for this average -------------------------------
            % Block k of the loaded sequence plays frequency obj.walkOrder(k); the
            % data is filed back under its frequency below.
            if numel(obj.walkOrder) ~= nF
                obj.walkOrder = 1:nF;      % an upload made before walkOrder existed was ascending
            end
            if obj.newOrderEveryAverage && obj.randomizeFrequencyOrder && ~obj.walkOrderFresh
                if isempty(obj.walkParts)
                    obj.sendWarning(['newOrderEveryAverage: the loaded sequence predates this ' ...
                        'feature, so it cannot be rebuilt on its own. Keeping its order for this ' ...
                        'average; change any setting (or invalidate the AWG cache) to reload.']);
                else
                    [tSeqAvg, tArmAvg] = obj.uploadWalkSequence(obj.drawWalkOrder());
                    obj.markWalkLoaded(obj.awgLoadFingerprint());
                    fprintf('Average %d: new visit order -- sequence %.1f s + arm %.1f s.\n', ...
                            obj.currIter, tSeqAvg, tArmAvg);
                end
            end
            obj.walkOrderFresh = false;
            obj.lastVisitOrder = obj.walkOrder;

            for r = 1:nRuns
                if obj.checkEmergencyStop()
                    return
                end
                success = false;
                for trial = 1:5
                    try
                        if obj.resyncEveryPass
                            % Restart the sequencer at its first step so this pass
                            % is known to begin on block 1. A plain stop-then-run
                            % costs ~97 ms, against ~499 ms after an asset change,
                            % and it bounds any trigger miscount to one pass.
                            obj.AWG.sendCommand('awgcontrol:stop:immediate');
                            obj.AWG.Run;
                        else
                            obj.AWG.EnsureRunning;
                        end

                        data = obj.getRawData(pg, spcm);
                        if numel(data) < readsPerRun
                            error('Walk pass returned %d reads, expected %d', numel(data), readsPerRun);
                        end
                        passData = reshape(data(1:readsPerRun), ...
                                           obj.detectionPeriodsPerRepeat, blocksPerPass, obj.passesPerRun);
                        % Drop the warm-up blocks and file each block under the
                        % frequency it played: block k played obj.walkOrder(k).
                        raw(:, obj.walkOrder, passesDone + (1:obj.passesPerRun)) = ...
                            passData(:, nWarmup+1:end, :);
                        passesDone = passesDone + obj.passesPerRun;
                        success = true;
                        break
                    catch err
                        err2warning(err);
                        failedPasses = failedPasses + 1;
                        fprintf('Walk pass %d failed at trial %d, attempting again.\n', r, trial);
                        try
                            spcm.stopExperimentCount;
                        catch
                            % Maybe the resources did not need clearing. That is fine.
                        end
                    end
                end
                if ~success
                    obj.emergencyStop;
                    return
                end

                % Reduce what we have so far, so the plot fills in as passes
                % accumulate rather than only at the end of the average.
                usable = min(passesDone, passesWanted);
                % Pass the dimensions explicitly: raw(:,:,1:1) arrives as
                % [2 x nFreq] because MATLAB drops trailing singletons.
                [sig, err] = walkReduce(raw(:, :, 1:usable), ...
                                        obj.detectionPeriodsPerRepeat, nF, opts);
                obj.signal(:, :, obj.currIter) = sig;
                obj.sterr(:, :, obj.currIter) = err;
                obj.updateSignalParam();

                obj.currParamIter = min(nF, round(nF * usable / passesWanted));
                if r == 1 || mod(r, 5) == 0 || r == nRuns
                    sendEventParamIterationDone(obj);
                end

                % Tracking runs once per pass instead of once per point -- there is
                % no per-point pause any more. The whole-sweep mean reference is
                % the brightness metric.
                if obj.isTracking
                    if obj.checkEmergencyStop()
                        return
                    end
                    refLevel = mean(sig(2, :));
                    refError = mean(err(2, :));
                    if tracker.compareReference(refLevel, refError, ...
                            Tracker.REFERENCE_TYPE_KCPS, obj.trackThreshhold)
                        tracker.trackUsing(TrackablePosition.NAME, obj)
                        obj.prepare;    % Before the next pass
                    end
                end
            end

            if failedPasses > 0
                % A retry aborts a pass part-way through, so the sequencer is left
                % on some arbitrary frequency. That is exactly what
                % resyncEveryPass repairs: the next pass stops the AWG and re-runs
                % it, which puts the sequencer back on step 1. So with
                % resyncEveryPass on, a retry costs one pass of data, not the
                % frequency mapping. Only worry about it when it is off.
                if obj.resyncEveryPass
                    fprintf(['Walk mode: %d of %d passes retried. resyncEveryPass is on, so each ' ...
                             'retry re-armed from step 1 and the mapping is intact.\n'], ...
                             failedPasses, nRuns);
                else
                    obj.sendWarning(sprintf(['Walk mode: %d of %d passes had to be retried with ' ...
                        'resyncEveryPass OFF. An aborted pass leaves the sequencer on an unknown ' ...
                        'frequency and nothing puts it back, so every later pass can be shifted. ' ...
                        'Turn resyncEveryPass on.'], failedPasses, nRuns));
                end
            else
                fprintf('Walk mode: %d passes, none retried.\n', nRuns);
            end
        end

        function opts = walkReduceOptions(obj, spcm)
            % Gather everything walkReduce needs to reproduce processData exactly.
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

        function perform(obj)
            % Gap power: (re)arm ao0 at the start of every average. Re-arming resets the
            % 2-level waveform to "gap value next", so a broken edge pair can never
            % survive past one average.
            obj.armGapAom();

            if obj.useSequenceWalk
                obj.performWalk();
                return
            end

            % Initialization
            len = length(obj.frequency);
            % VISIT ORDER. Random by default and it must stay that way for real data:
            % sequential order lets slow drift masquerade as spectral structure.
            %
            % But it is recorded now, and it can be made sequential, because the jump-mode
            % failure cannot be diagnosed without it. Jump mode collects the RIGHT amount
            % of contrast (2.16% against the classic path's 2.10%, min 0.9784, median
            % 0.9921 -- nearly the same distribution) and files it under the WRONG
            % frequencies, giving a correlation of +0.044. The leading explanation is a
            % one-step lag: the value read at visit k belongs to the frequency visited at
            % k-1. In RANDOM visit order such a lag scatters the spectrum into something
            % uncorrelated, which is exactly what +0.044 looks like -- and it is invisible
            % to any frequency-space shift search, which is why "best shift" kept
            % reporting 0.
            %
            % With the order recorded, the lag hypothesis can be tested by UNDOING it
            % exactly, with no fitted parameters. With the order sequential, the same lag
            % becomes a plain shift that the existing best-shift metric finds on its own.
            if obj.randomizeFrequencyOrder
                f1 = randperm(len);
            else
                f1 = 1:len;
            end
            obj.lastVisitOrder = f1;    % f1(k) is the frequency INDEX visited kth

            %%% Devices (+ Tracker)
            pg = getObjByName(PulseGenerator.NAME);
            spcm = getObjByName(Spcm.NAME);
            tracker = getObjByName(Tracker.NAME);
            if isempty(tracker); throwBaseObjException(Tracker.Name); end

            % How often to verify the AWG is still armed. See armCheckEvery.
            if isempty(obj.armCheckEvery)
                checkEvery = BooleanHelper.ifTrueElse(obj.useSequenceJump, 25, 1);
            else
                checkEvery = obj.armCheckEvery;
            end

            % Some magic numbers
            n = obj.detectionPeriodsPerRepeat * obj.runsPerPerform;

            sig = zeros(1, n);
            sterr = zeros(1, n);

            obj.temp = [];

            % Run - Go over all frequencies, in random order
            for k = 1:len
                if obj.currIter <= size(obj.signal, 3) && ...
                        sum(obj.signal(:, f1(k), obj.currIter)) ~= 0
                    continue
                end
                success = false;

                for trial = 1 : 5
                    if obj.checkEmergencyStop()
                        return;
                    end
                    try
                        i = f1(k);
                        if obj.useFgen
                            % No sequencer: just retune the built-in generator.
                            obj.selectFgenFrequency(obj.frequency(i));
                        elseif obj.useSequenceJump
                            % One sequence already assigned and armed; just move
                            % the sequencer to this frequency's block.
                            obj.AWG.jumpToStep(obj.indexSeq(1,i), 1);
                            if obj.nChannels == 2
                                obj.AWG.jumpToStep(obj.indexSeq(2,i), 2);
                            end
                        else
                            obj.AWG.assignSequence(obj.indexSeq(1,i),1)
                            if obj.nChannels == 2
                                obj.AWG.assignSequence(obj.indexSeq(2,i),2);
                            end
                        end
                        % Only poll the AWG's run state when it might tell us
                        % something. Always on the first point of the average,
                        % always on a retry, otherwise every checkEvery points.
                        % Nothing to arm in FGEN mode -- the generator is simply on.
                        if ~obj.useFgen && checkEvery > 0 && ...
                                (k == 1 || trial > 1 || mod(k-1, checkEvery) == 0)
                            obj.AWG.EnsureRunning;
                        end
                        data = obj.getRawData(pg, spcm);

                        [sig(1:2), sterr(1:2)] = obj.processData(data);

                        if ~obj.isSingleMeasurement % run another sweep with the same source
                            data = obj.getRawData(pg, spcm);
                            [sig(3:4), sterr(3:4)] = obj.processData(data);
                        end

                        obj.signal(:, i, obj.currIter) = sig;
                        obj.sterr(:, i, obj.currIter) = sterr;

                        success = true;
                        obj.currParamIter = obj.currParamIter + 1;
                        if obj.currParamIter == 1 || mod(obj.currParamIter, 5) == 0 || obj.currParamIter == obj.totalNumberOfParams
                            sendEventParamIterationDone(obj);
                        end

                        if obj.isTracking
                            if obj.checkEmergencyStop()
                                return;
                            end
                            isTrackingNeeded = tracker.compareReference(...
                                sig(2), sterr(2), ... % 2 for reference
                                Tracker.REFERENCE_TYPE_KCPS, obj.trackThreshhold);
                            if isTrackingNeeded
                                tracker.trackUsing(TrackablePosition.NAME, obj)
                                obj.prepare;    % Before next measurement
                            end
                        end

                        break;
                    catch err
                        err2warning(err);
                        fprintf('Experiment failed at trial %d, attempting again.\n', trial);
                        try
                            % Maybe we need to manually clear the resources
                            spcm.stopExperimentCount;
                        catch
                            % But maybe we don't, and that's perfectly ok.
                        end
                    end
                end
                if ~success
                    obj.emergencyStop;
                    return
                end
            end

            % Shared with performWalk so the two paths cannot drift apart.
            obj.updateSignalParam();
        end

        function [value, sterr] = getRatioDistributionValues(obj, S1, S2, S1sterr, S2sterr)
            % Overrides Experiment.getRatioDistributionValues.
            %
            % The inherited version combined S1 and S2 across averages
            % SEPARATELY (via getCombinedSterr) and only formed the ratio at the
            % end. getCombinedSterr is a pooled-variance formula, so its
            % cross-average term charges the difference in mean LEVEL between
            % averages as if it were statistical scatter -- which is exactly the
            % thing dividing by the reference is meant to be immune to. Measured
            % up to ~10x error-bar inflation at 20% drift while the value stayed
            % correct. Worse on the photodiode path, where processData sets the
            % reference column to a constant so its sterr is identically 0, and
            % getCombinedSterr then returned a number driven purely by drift.
            %
            % ratioAcrossAverages forms the ratio per average first -- where the
            % level drift cancels -- and only then combines. See
            % prod/Util/ratioAcrossAverages.m and testRatioAcrossAverages.m.
            %
            % nParams is passed explicitly because the incoming shape is
            % ambiguous: squeeze() hands over a 1 x nParams ROW at one average
            % and an nParams x nAverages MATRIX at two or more.
            [value, sterr] = ratioAcrossAverages(S1, S2, S1sterr, S2sterr, ...
                                                obj.getTotalNumberOfParams());
        end

        function wrapUp(obj)
            % Things that need to happen when the experiment is done; a
            % counterpart for obj.prepare.
            % In the future, it will also analyze results and fit from it
            % the resonance frequency/ies

            obj.wrapUpInternal()

            % Gap power: free ao0 and put the window value back on it.
            obj.releaseGapAom();

            % Never leave the instrument in FGEN mode. Rabi, Echo, T1 and Ramsey
            % all upload sequences, which an instrument in FGEN mode cannot play,
            % and they would only find out by producing meaningless data -- worse,
            % a run that skips its reload would not even try to fix the mode. So
            % ESR puts it back itself.
            if obj.useFgen
                try
                    obj.ensureInstrumentMode('AWG');
                catch err
                    obj.sendWarning(sprintf(['Could not return the AWG to AWG mode: %s. ' ...
                        'Set it back before running Rabi, Echo, T1 or Ramsey.'], err.message));
                end
            end
        end

        function params = alternateSignal(obj)
            % Returns alternate view of the data.
            % Return a cell array ExpParams. The first cell is dataParam.
            % The second cell is either another dataParam or a cell of
            % three dataParams. The the last cell is a cell of strings
            % which indiciates what are the ExpParams. It is either
            % {'yParam', 'yParam2'} or {'yParam', 'yParam2+'}

            if obj.currIter == 0
                params = ExpResultDoubleVector('fluorescence', [], [], Spcm.signalUnits(), obj.NAME);
                return
            end
            S1 = squeeze(obj.signal(1, :, 1:obj.currIter));
            S1sterr = squeeze(obj.sterr(1, :, 1:obj.currIter));
            S2 = squeeze(obj.signal(2, :, 1:obj.currIter));
            S2sterr = squeeze(obj.sterr(2, :, 1:obj.currIter));
            if obj.currIter ~= 1
                % Calculate the mean
                S1sterr = obj.getCombinedSterr(S1, S1sterr);
                S2sterr = obj.getCombinedSterr(S2, S2sterr);
                S1 = mean(S1, 2);
                S2 = mean(S2, 2);
            end
            dataParam = ExpResultDoubleVector('fluorescence', S1, S1sterr, Spcm.signalUnits(), obj.NAME, 'With MW');
            dataParam2 = ExpResultDoubleVector('fluorescence', S2, S2sterr, Spcm.signalUnits(), obj.NAME, 'Without MW');

            if ~obj.isSingleMeasurement
                S3 = squeeze(obj.signal(3, :, 1:obj.currIter));
                S3sterr = squeeze(obj.sterr(3, :, 1:obj.currIter));
                S4 = squeeze(obj.signal(4, :, 1:obj.currIter));
                S4sterr = squeeze(obj.sterr(4, :, 1:obj.currIter));

                if obj.currIter ~= 1
                    % Calculate the mean
                    S3sterr = obj.getCombinedSterr(S3, S3sterr);
                    S4sterr = obj.getCombinedSterr(S4, S4sterr);
                    S3 = mean(S3, 2);
                    S4 = mean(S4, 2);
                end
                dataParam3 = ExpResultDoubleVector('fluorescence', S3, S3sterr, Spcm.signalUnits(), obj.NAME, 'With MW - Mirrored');
                dataParam4 = ExpResultDoubleVector('fluorescence', S4, S4sterr, Spcm.signalUnits(), obj.NAME, 'Without MW - Mirrored');
                params = {dataParam, {dataParam2, dataParam3, dataParam4}, {'yParam', 'yParam2+'}};
            else
                params = {dataParam, dataParam2, {'yParam', 'yParam2'}};
            end
        end
    end

    %% Laser power during the gaps (this class only)
    methods
        function releaseGapAom(obj)
            % Stop the hardware-timed ao0 task (if armed), report how many edges the
            % DAQ received, and write the window value back through the laser panel.
            % Safe to call at any time, e.g. after a crash or Ctrl-C.
            if isempty(obj.gapAomTask)
                return
            end
            task = obj.gapAomTask;
            obj.gapAomTask = [];
            n = NaN;
            try
                [st, n] = daq.ni.NIDAQmx.DAQmxGetWriteTotalSampPerChanGenerated(task, uint64(0));
                if st ~= 0; n = NaN; end
                n = double(n);
            catch
            end
            obj.gapAomEdges = n;
            try; DAQmxStopTask(task); catch; end
            try; DAQmxClearTask(task); catch; end
            % Put the window value back with an ordinary static write (also covers
            % devices that drop to 0 V when a task is cleared).
            if ~isempty(obj.gapAomWindowPanel)
                try
                    part = AWGExpESR_gapAom.gapAomPart();
                    part.value = obj.gapAomWindowPanel;
                catch err
                    obj.sendWarning(sprintf(['Gap AOM: could not restore the window value ' ...
                        '(%g): %s. Set the laser panel by hand.'], obj.gapAomWindowPanel, err.message));
                end
            end
            if n == 0
                obj.sendWarning(sprintf(['Gap AOM: the DAQ received NO marker edges, so the gap ' ...
                    'power was never applied. Check the cable from PulseStreamer output %d to DAQ %s.'], ...
                    obj.gapAomPgChannel, obj.gapAomDaqTerminal));
            elseif isfinite(n) && mod(n, 2) ~= 0
                obj.sendWarning(sprintf(['Gap AOM: the DAQ received an odd number of edges (%d) -- ' ...
                    'an edge pair was broken. Re-arming each average re-aligns it, but check the ' ...
                    'marker signal on a scope.'], n));
            end
        end
    end

    methods (Access = protected)
        function registerGapAomChannel(obj)
            % Register the marker output with the pulse generator, once. The json's
            % pulseGenerator.channelNames list is not what registers channels -- the
            % devices do -- so a free output has to be registered here before a
            % sequence may use it.
            if isempty(obj.gapAomValue)
                return
            end
            pg = getObjByName(PulseGenerator.NAME);
            if isempty(pg); throwBaseObjException(PulseGenerator.NAME); end
            names = pg.channelNames;
            k = find(strcmp(names, obj.GAP_AOM_PG_NAME), 1);
            if ~isempty(k)
                if pg.channelName2Address(obj.GAP_AOM_PG_NAME) ~= obj.gapAomPgChannel
                    obj.sendError(sprintf(['Gap AOM: "%s" is already registered on another ' ...
                        'PulseStreamer output. Restart MATLAB''s setup, or set gapAomPgChannel = %d.'], ...
                        obj.GAP_AOM_PG_NAME, pg.channelName2Address(obj.GAP_AOM_PG_NAME)));
                end
                return
            end
            if ismember(obj.gapAomPgChannel, pg.channelAddresses)
                obj.sendError(sprintf(['Gap AOM: PulseStreamer output %d is already used by ' ...
                    'another channel. Pick a free output for gapAomPgChannel.'], obj.gapAomPgChannel));
            end
            pg.registerChannel(Channel.Digital(obj.GAP_AOM_PG_NAME, obj.gapAomPgChannel, 0, 0));
        end

        function addGap(obj, S)
            % One delayBetweenDetections gap. Laser on throughout, exactly as in
            % AWGExpESR. With the gap power on, two marker pulses: the first steps ao0
            % to the gap value, the second steps it back gapAomReturnLead before the
            % next window.
            G = obj.delayBetweenDetections;
            if isempty(obj.gapAomValue)
                S.addEvent(G, {'greenLaser'});
                return
            end
            [dur, on] = AWGExpESR_gapAom.gapAomEvents(G, obj.gapAomReturnLead, ...
                obj.gapAomEdgeWidth, obj.GAP_AOM_PG_NAME);
            for k = 1:numel(dur)
                S.addEvent(dur(k), on{k});
            end
        end

        function captureGapAomValues(obj)
            % Window value = the panel value now; gap value = gapAomValue. Both in
            % panel units (volts on Setup 7).
            obj.gapAomVolts = [];
            if isempty(obj.gapAomValue)
                obj.gapAomWindowPanel = [];
                return
            end
            part = AWGExpESR_gapAom.gapAomPart();
            if ~ValidationHelper.isInBorders(obj.gapAomValue, part.minValue, part.maxValue)
                obj.sendError(sprintf('Gap AOM: gapAomValue = %g is outside the AOM range [%g, %g].', ...
                    obj.gapAomValue, part.minValue, part.maxValue));
            end
            % The window value: what the panel last wrote, or a readback if it never did.
            win = [];
            if isprop(part, 'valueInternal'); win = part.valueInternal; end
            if isempty(win); win = part.value; end
            win = round(win, 4);
            obj.gapAomWindowPanel = win;
            % Panel units ARE volts here: LaserPartAbstract's percentage calibration
            % (calibFunc) is never set, so the panel writes its number as a voltage.
            obj.gapAomVolts = [obj.gapAomValue, win];
            fprintf(['Gap AOM: windows at %g, gaps at %g (%s); marker on PulseStreamer output %d ' ...
                     '-> DAQ %s; back to the window value %g us before each window.\n'], ...
                     win, obj.gapAomValue, strtrim(part.units), obj.gapAomPgChannel, ...
                     obj.gapAomDaqTerminal, obj.gapAomReturnLead);
        end

        function armGapAom(obj)
            % (Re)start ao0 as a 2-level waveform clocked by the marker edges.
            obj.releaseGapAom();
            if isempty(obj.gapAomValue)
                return
            end
            if isempty(obj.gapAomVolts)
                obj.captureGapAomValues();
            end
            part = AWGExpESR_gapAom.gapAomPart();
            niDaq = getObjByName(NiDaq.NAME);
            if isempty(niDaq); throwBaseObjException(NiDaq.NAME); end
            if niDaq.dummyMode
                obj.sendWarning('Gap AOM: the NiDaq is in dummy mode -- gap power not armed.');
                return
            end
            dev = niDaq.deviceName;
            phys = sprintf('/%s/%s', dev, part.niDaqChannel);
            clk = sprintf('/%s/%s', dev, obj.gapAomDaqTerminal);
            v = obj.gapAomVolts;                    % [gap, window]
            lo = min([part.minValue, part.maxValue, v]);
            hi = max([part.minValue, part.maxValue, v]);

            % The edges come gapAomReturnLead and (gap - lead) apart; the DAC must keep up.
            maxRate = 1e5;                          % conservative default
            try
                [st, r] = daq.ni.NIDAQmx.DAQmxGetDevAOMaxRate(dev, double(0));
                if st == 0 && r > 0; maxRate = double(r); end
            catch
            end
            minGap = min(obj.gapAomReturnLead, obj.delayBetweenDetections - obj.gapAomReturnLead);
            if minGap * 1e-6 < 1 / maxRate
                obj.sendError(sprintf(['Gap AOM: marker edges only %.2f us apart, but the DAQ ' ...
                    'updates at most every %.2f us. Lengthen delayBetweenDetections or ' ...
                    'gapAomReturnLead.'], minGap, 1e6 / maxRate));
            end

            [st, ~, task] = DAQmxCreateTask([]);
            obj.daqCheck(st, [], 'create the task');
            try
                st = DAQmxCreateAOVoltageChan(task, phys, lo, hi, []);
                obj.daqCheck(st, task, sprintf('use %s (is the laser panel writing it?)', phys));
                st = DAQmxCfgSampClkTiming(task, clk, maxRate, ...
                    daq.ni.NIDAQmx.DAQmx_Val_Rising, daq.ni.NIDAQmx.DAQmx_Val_ContSamps, 8);
                obj.daqCheck(st, task, sprintf('clock %s from %s', phys, clk));
                buf = repmat(v, 1, 4);              % gap, window, gap, window, ... (even length)
                st = DAQmxWriteAnalogF64(task, numel(buf), 0, 10, ...
                    daq.ni.NIDAQmx.DAQmx_Val_GroupByChannel, buf, 0);
                obj.daqCheck(st, task, 'write the 2-level buffer');
                st = DAQmxStartTask(task);
                obj.daqCheck(st, task, 'start the task');
            catch err
                try; DAQmxClearTask(task); catch; end
                rethrow(err)
            end
            obj.gapAomTask = task;
        end

        function daqCheck(obj, status, task, what)
            % Negative DAQmx status = error (clear the task, stop); positive = warning.
            if status == 0
                return
            end
            msg = '';
            try
                buf = char(ones(1, 500));
                [st2, msg] = daq.ni.NIDAQmx.DAQmxGetErrorString(int32(status), buf, uint32(500));
                if st2 ~= 0; msg = ''; end
            catch
            end
            msg = strtrim(char(msg));
            if status < 0
                if ~isempty(task)
                    try; DAQmxClearTask(task); catch; end
                end
                obj.sendError(sprintf('Gap AOM: could not %s. DAQmx error %d: %s', what, status, msg));
            else
                obj.sendWarning(sprintf('Gap AOM: while trying to %s, DAQmx warning %d: %s', what, status, msg));
            end
        end
    end

    methods (Static)
        function [dur, on] = gapAomEvents(G, lead, w, name)
            % The four events of one gap: marker pulse, gap power, marker pulse, lead.
            % Laser on in all four. Their durations add up to G exactly.
            if lead <= w + 1e-9
                error('Gap AOM: gapAomReturnLead (%g us) must be longer than gapAomEdgeWidth (%g us).', lead, w);
            end
            if G - lead - w < w - 1e-9
                error(['Gap AOM: delayBetweenDetections (%g us) is too short: it needs at least ' ...
                       'gapAomReturnLead + 2*gapAomEdgeWidth = %g us.'], G, lead + 2*w);
            end
            dur = [w, G - lead - w, w, lead - w];
            on  = {{'greenLaser', name}, {'greenLaser'}, {'greenLaser', name}, {'greenLaser'}};
        end

        function part = gapAomPart()
            % The Green Laser's AOM part -- the same one checkGreenLaserPower uses.
            laser = getObjByName('Green Laser');
            if isempty(laser); throwBaseObjException('Green Laser'); end
            [~, parts] = laser.getContollableParts();
            if isempty(parts)
                error('Gap AOM: the Green Laser has no settable part.');
            end
            part = parts{1};
            if ~isprop(part, 'niDaqChannel')
                error(['Gap AOM: the Green Laser''s first part (%s) is not a NiDaq-controlled AOM, ' ...
                       'so its level cannot be hardware-timed.'], class(part));
            end
        end
    end
end
