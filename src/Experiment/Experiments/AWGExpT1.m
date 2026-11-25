classdef AWGExpT1 < ExperimentAWG
    %EXPECHO Echo experiment

    properties (Constant)
        NAME = 'AWG T1';
    end
    
    properties
        % All of the parameters must have a set function, and when they are
        % changed, the parameter obj.changeFlag needs to be set to true
        % (unless the change is minor and the experiment can continue
        % running);
        tau                 % in us 
        doubleMeasurement   % logical
        channelPhase = [0 0];       % DEGREES, per AWG channel. The relative phase of the two
                                    % microwave drives, from calibMwPair (level 4.2).
                                    %
                                    % WHY THIS EXISTS AND timeDelay DOES NOT DO THE JOB.
                                    % timeDelay expresses a phase as a LENGTH of silent
                                    % padding in front of one channel's waveform, so its
                                    % resolution is the sample grid and it can only ever
                                    % delay -- there is no negative padding. channelPhase is
                                    % written into the cosine argument itself, so it is exact
                                    % to double precision and signed.
                                    %
                                    % The two are ADDITIVE where they are applied, so with
                                    % both non-zero the real phase is their sum and nothing
                                    % records that. Level 4.3 sets timeDelay = 0 for exactly
                                    % this reason; leave it there.

        mwGateMode = 'tau';         % 'tau' | 'pulse' | 'all'. WHEN THE MICROWAVE SWITCH IS OPEN.
                                    %
                                    % ==========================================================
                                    % WHY THIS IS A KNOB AND NOT JUST A FIX
                                    % ==========================================================
                                    %
                                    % T1 has ONE microwave pulse and it is at the START of tau.
                                    % The AWG plays for 0.335 us after each trigger -- a 0.15 us
                                    % zero pad, then piTime = 0.035 plus another 0.15 of pad --
                                    % and then it is idle until the next trigger. Nothing needs
                                    % the switch open after that.
                                    %
                                    % It was open for the whole of tau anyway. The sequence read
                                    %
                                    %     S.addEvent(TRIG_DURATION, [triggers, MW]);
                                    %     S.addEvent(tau,           MW,          'tau');
                                    %
                                    % which at tau = 4900 us leaves the switch open for 4900 us
                                    % where 0.3 is needed. For scale: echoGateTest bounded the
                                    % gate as harmless out to 240 us of excess, and Ramsey's
                                    % measured-harmless case was 4.05 us. This is 20x the former
                                    % and 1200x the latter.
                                    %
                                    % Two things make it worse here than anywhere else, and both
                                    % are about T1 specifically rather than about size.
                                    %
                                    % THE ARTEFACT HAS THE SHAPE OF THE SIGNAL. The excess is
                                    % PROPORTIONAL to tau, so anything leaking through removes
                                    % contrast as exp(-r*tau) and simply adds to the relaxation
                                    % rate: 1/T1_measured = 1/T1_true + r. There is no feature to
                                    % notice. Echo's fixed gate at least put its worst excess at
                                    % SHORT tau, against the direction of the decay.
                                    %
                                    % AN OPEN SWITCH PASSES NOISE, NOT ONLY CARRIER. Broadband
                                    % noise out of the amplifier chain drives spin flips off
                                    % resonance, and driven spin flips are exactly what T1
                                    % measures. So the mechanism is not a subtle one.
                                    %
                                    % AND THE LAB'S OWN NON-AWG ExpT1 CLOSES IT: that file has
                                    % S.addEvent(obj.tau(end), '', 'tau') -- MW off -- and opens
                                    % the gate only for piTime. This class inherited Ramsey's
                                    % gate window, and Ramsey genuinely needs it, because its
                                    % second pi/2 lands at the END of tau.
                                    %
                                    % ==========================================================
                                    % THE THREE MODES, AND WHY THE DEFAULT IS STILL THE OLD ONE
                                    % ==========================================================
                                    %
                                    %   'tau'    MW open for the whole tau. What this class has
                                    %            always done. Gate-open time = tau.
                                    %   'pulse'  MW open for mwHoldUs only, then closed for the
                                    %            rest of the wait. What ExpT1 does, and what the
                                    %            physics asks for. Gate-open time = mwHoldUs.
                                    %   'all'    MW open through tau AND lastDelay. Only sensible
                                    %            with constantTime = true, where it gives a
                                    %            CONSTANT gate-open time of lastDelay + max(tau)
                                    %            on every point. This is the deliberately
                                    %            over-dosed arm for t1GateTest -- see below.
                                    %
                                    % The default is 'tau', i.e. unchanged. That is deliberate:
                                    % switching it would silently change every T1 this setup
                                    % measures, and there is no measurement yet saying which way
                                    % is better. Run t1GateTest first. If it finds no leak,
                                    % 'pulse' is still the right default on the noise argument
                                    % alone and costs nothing; if it finds one, every T1 taken
                                    % before now was short and by a known factor.
                                    %
                                    % HOW t1GateTest USES THE THREE MODES. With a leak rate r per
                                    % microsecond of open gate, the extra decay factor is
                                    % exp(-r*G) where G is the gate-open time, so the three modes
                                    % carry three DIFFERENT and independently measurable
                                    % signatures against the same physics:
                                    %
                                    %       'pulse'  G ~ 0            -> the baseline
                                    %       'tau'    G = tau          -> log-ratio SLOPE of -r
                                    %       'all'    G = const        -> constant OFFSET, no slope
                                    %
                                    % A real leak has to show up as the same r in all three. A
                                    % drift, a tracking step or a laser change cannot fake that,
                                    % which is the point of having a third arm rather than two.

        mwHoldUs = [];              % how long the switch stays open in 'pulse' mode, in us.
                                    % [] = derived in prepare as 2*minWaveformDuration + piTime
                                    % plus a margin, i.e. long enough to cover everything the AWG
                                    % actually plays after a trigger. Set it by hand only to test
                                    % the margin; making it SHORTER than the playback clips the pi
                                    % pulse, which costs contrast at every tau equally and so
                                    % looks like a calibration error rather than a timing one.
                                    % prepare prints the value it used.
    end

    properties (SetAccess = protected)
        mwGateOpenUs                % gate-open time per repeat-half at the LONGEST tau, filled in
                                    % by prepare. Recorded rather than recomputed so that saved
                                    % data carries the dose it was taken at.
        mwHoldUsed                  % the hold prepare actually used, in us. perform subtracts it
                                    % from every tau it writes, and reading it from here rather
                                    % than recomputing it is what stops the two drifting apart.
        usePulseGate = false;        % whether the sequence has a separate 'mwHold' event
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = AWGExpT1()
            % AWGExpT1.NAME, not AWGExpRamsey.NAME.
            %
            % Experiment's constructor does removeObjIfExists(name) and then addBaseObject(obj)
            % under that name, so this class was registering itself as 'AWG Ramsey'. Every
            % consequence of that is real and none of them look like a bug at the time:
            %
            %   * constructing an AWGExpT1 DELETED any live AWGExpRamsey, printing "Creating a
            %     new AWG Ramsey Experiment. All of its parameters were reset to default." The
            %     run script makes a Ramsey and then a T1 in one session, so it fired on every
            %     run, and the message names the wrong experiment so it reads as noise.
            %   * getObjByName('AWG Ramsey') returned the T1 afterwards.
            %   * EventSender(name) routed T1's progress events to the Ramsey view.
            %   * Savable(name) filed autosaved T1 data under 'AWG Ramsey'. The T1 cell in the
            %     run script never sets shouldAutosave = false, so that actually happened.
            %
            % Five other classes in this folder do the same thing -- AWGExpRabiCont and
            % AWGExpRabiHyperfineCont both take AWGExpRabi.NAME, the two PiSpinLockDualChannel
            % variants take AWGExpPiSpinLock.NAME, and AWGExpIcoDynamicalBase takes
            % AWGExpIcoBase.NAME. Some of those may be deliberate aliasing so a pair of classes
            % shares one GUI slot; this one was not, because AWGExpT1 has its own NAME constant
            % and uses it for every ExpParam below.
            obj@ExperimentAWG(AWGExpT1.NAME);
            obj.parameterName = 'taus';
            
            % Set properties inherited from Experiment
            obj.repeats = 100;
            obj.averages = 10;
            
            obj.frequency = 2687;           % in MHz
            obj.amplitude = -15;            % in dBm
            obj.tau = 10:100:3000;                   
            obj.piTime = 0.035;         % in us
            obj.lastDelay = Experiment.DEFAULT_LAST_DELAY;
            obj.AWGorSRSswitch = 0;

            obj.constantTime = false;       % logical
            obj.doubleMeasurement = true;   % logical
            
            obj.detectionDuration = 0.5;           % detection window, in us
            obj.referenceDetectionDuration = 0.5;     % in us. Detection duration of the reference read
            obj.laserInitializationDuration = 10;   % laser initialization in pulsed experiments in us
            
            obj.mCurrentXAxisParam = ExpParamDoubleVector('Time', [], [], StringHelper.MICROSEC, obj.NAME);
            obj.topParam = ExpParamDoubleVector(StringHelper.TAU, [], [], StringHelper.MICROSEC, obj.NAME);
            % THE LABELS ARE THE OTHER WAY ROUND FROM ECHO AND RAMSEY, because the sequence is.
            %
            % LoadAWG builds [zero, noPi, zero, pi] with triggers [1 0 1 0], so the FIRST
            % pulse-generator half plays the ZEROED pi -- the spin sits in |0> through tau --
            % and the second plays the real one. Detection rows 1-2 are therefore the |0>
            % relaxation and rows 3-4 the |1> relaxation. Echo and Ramsey put the signal in the
            % first half, and these labels had been copied from them, so signalParam was
            % announcing rows 1-2 as |1>.
            %
            % Nothing about the measurement changes here -- only which curve gets which name --
            % but the name is the whole content of a T1 plot: the |1> trace is the one that
            % decays and the |0> trace is the flat control. Read backwards, a flat line looks
            % like an infinite T1 and a real decay looks like a broken reference.
            %
            % This is a reading of the sequence rather than a measurement, so checkProjectionLabels
            % below tests it against the DATA on every run and warns if it disagrees.
            obj.signalParam = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME, '|0>');
            obj.signalParam2 = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME, '|1>');
        end
    end
    
    %% Setters
    methods

        function set.channelPhase(obj, newVal)
            if ~isnumeric(newVal) || ~isvector(newVal)
                obj.sendError('channelPhase must be a numeric vector, in degrees')
            end
            % Wrapped into [0,360) on the way in. The phase is periodic, and letting it run
            % outside one period is how a scan ends up with several equivalent optima and
            % reports whichever one it reached first.
            obj.channelPhase = mod(newVal(:)', 360);
            obj.changeFlag = true;      % the waveform SAMPLES change, so it must reload
        end
        function set.doubleMeasurement(obj, newVal)
            checkBoolean(obj, newVal)
            % If we got here, then newVal is OK.
            obj.doubleMeasurement = newVal;
            obj.changeFlag = true;
        end

        function set.mwGateMode(obj, newVal)
            if ~ischar(newVal) && ~isstring(newVal)
                obj.sendError('mwGateMode must be one of ''tau'', ''pulse'', ''all''.')
            end
            v = lower(char(newVal));
            if ~ismember(v, {'tau', 'pulse', 'all'})
                obj.sendError(sprintf(['mwGateMode was ''%s''. It must be ''tau'' (open for the ' ...
                    'whole wait, the historical behaviour), ''pulse'' (open only for the pi, what ' ...
                    'ExpT1 does) or ''all'' (open through tau AND lastDelay, the over-dosed arm ' ...
                    'for t1GateTest).'], v))
            end
            obj.mwGateMode = v;
            % The pulse-generator sequence changes SHAPE, not just a duration, so prepare has to
            % rebuild it. changeFlag also forces LoadAWG, which is not needed -- the AWG side is
            % tau-independent and identical in all three modes -- but a spurious waveform reload
            % costs one upload and getting this wrong costs a whole run.
            obj.changeFlag = true;
        end

        function set.mwHoldUs(obj, newVal)
            if ~isempty(newVal) && (~isnumeric(newVal) || ~isscalar(newVal) || newVal <= 0)
                obj.sendError('mwHoldUs must be a positive scalar in microseconds, or [] to derive it.')
            end
            obj.mwHoldUs = newVal;
            obj.changeFlag = true;
        end
    end
    
    methods
        function totalParamNum = getTotalNumberOfParams(obj)
            totalParamNum = length(obj.tau);
        end
    end
    
    %% Overridden from Experiment
    %% channelPhase
    methods
        function phi = channelPhaseRadians(obj)
            % channelPhase padded or trimmed to nChannels, in radians. Padded rather than
            % erroring: a two-element phase left in place while nChannels drops to 1 is a
            % normal thing to do mid-calibration, and it should not throw.
            %
            % Same implementation as AWGExpRabi's, deliberately duplicated rather than put in
            % ExperimentAWG: that file is loaded by every setup in the lab, and this change
            % only needs to reach Setup 7.
            phi = zeros(1, obj.nChannels);
            cp = obj.channelPhase(:)';
            k = min(numel(cp), obj.nChannels);
            phi(1:k) = deg2rad(cp(1:k));
        end
    end

    %% Channel names
    methods
        function names = mwChannelNames(obj)
            %MWCHANNELNAMES The microwave switch's pulse-generator channels, as a flat cellstr.
            %
            % ==============================================================================
            % WHY THIS EXISTS: THE ONE LINE THAT MADE EVERY T1 RUN FAIL
            % ==============================================================================
            %
            % prepare used to build the channel list as
            %
            %       case 2;  MW = {MWChannel, 'MW2'};       % MWChannel = obj.MWChannel
            %
            % and on this setup obj.MWChannel is NOT the plain char 'MW' that Experiment.m
            % defaults it to. Wrapping whatever it holds inside another cell gives a first element
            % that is itself a cell -- {{'MW','MW2'}, 'MW2'} -- and Pulse/changeChannels then
            % reaches
            %
            %       oc.(channel) = true;
            %       -> Argument to dynamic structure reference must evaluate to a valid field name
            %
            % which is where every T1 run was ending. Nothing in that message points at MWChannel,
            % and nothing points at T1 either: the traceback names Pulse, Sequence and prepare, so
            % it reads as a broken sequence rather than a mistyped channel list.
            %
            % There is more than one way for the property to stop being a char. jsondecode turns a
            % json array of strings into a cell, so loadJsonParams can do it; ExpDoubleRabi sets
            % obj.MWChannel = {'MW','MW2'} outright; and a numeric value would fail the same way.
            % This function does not care which -- it accepts a char or a cell of chars and refuses
            % anything else by name.
            %
            % AWGExpT1 WAS THE ONLY ONE OF THE FIVE THAT COULD HIT IT. Echo, Ramsey, Rabi and ESR
            % all write the two-channel case out as {'MW', 'MW2'} and never read the property
            % there. That is not an accident of style, it is the correct reading of what the
            % property means: at nChannels = 2 BOTH microwave channels are driven, so there is
            % nothing for MWChannel to select. It carries information only at nChannels = 1, where
            % it says WHICH of the two is live -- exactly how ExperimentAWG.activeChannels uses it,
            % via strcmp(obj.MWChannel, 'MW'). Reading it at two channels was a copy-paste that
            % happened to work for as long as the property held its default.
            %
            % So: two channels gives the literal pair, one channel flattens whatever MWChannel
            % holds, and every name is checked with isvarname before it leaves here. A bad name
            % caught HERE names the property and its value; a bad name that escapes fails four
            % frames deeper with a message about dynamic structure references.
            if obj.nChannels >= 2
                names = {'MW', 'MW2'};
                return
            end

            raw = obj.MWChannel;
            names = {};
            if ischar(raw)
                names = {raw};
            elseif iscell(raw)
                % One level of flattening, which is all jsondecode or a hand-set cell produces.
                for i = 1:numel(raw)
                    v = raw{i};
                    if ischar(v)
                        names{end+1} = v;                                       %#ok<AGROW>
                    elseif iscell(v)
                        for j = 1:numel(v)
                            if ischar(v{j}); names{end+1} = v{j}; end           %#ok<AGROW>
                        end
                    end
                end
            end
            if isempty(names)
                obj.sendError(sprintf(['MWChannel does not name a pulse-generator channel: it is ' ...
                    'a %s. It has to be a channel name such as ''MW'' or ''MW2'', or a cell of ' ...
                    'them. Experiment.m defaults it to the char ''MW''; a json load or another ' ...
                    'experiment class can leave it as something else.'], class(raw)))
            end
            % At one channel there is one live channel. More than one name here means MWChannel is
            % describing a two-channel setup while nChannels says otherwise, and driving both would
            % contradict activeChannels, which resolves to a SINGLE output.
            if numel(names) > 1
                obj.sendWarning(sprintf(['MWChannel lists %d channels (%s) but nChannels = 1. ' ...
                    'Using ''%s''. Set nChannels = 2 if both are meant to be driven.'], ...
                    numel(names), strjoin(names, ', '), names{1}))
                names = names(1);
            end
            for i = 1:numel(names)
                if ~isvarname(names{i})
                    obj.sendError(sprintf(['MWChannel gives ''%s'', which is not a usable ' ...
                        'pulse-generator channel name -- Pulse stores channels as struct fields, ' ...
                        'so a name has to be a valid identifier.'], names{i}))
                end
            end
        end
    end

    methods (Access = protected)
        function prepare(obj)
            % Initialize devices (SPCM, PulseGenerator, etc.)
            
            % Sequence
            %%% Useful parameters for what follows
            initDuration = obj.laserInitializationDuration-obj.detectionDuration;
            detectionSequnceDuration = obj.detectionDuration + obj.stabilizationDuration + obj.acquisitionDuration + obj.GIresetDuration;
            obj.detectionPeriodsPerRepeat = 2 * (1 + double(obj.doubleMeasurement));
            obj.runsPerPerform = 1;

            %%% Creating the sequence
            if obj.seperatedTrigger                         % AWG trigger
                triggers = {'trigger', 'trigger2'};
            else
                triggers = {'trigger'};
            end
            MW = obj.mwChannelNames();
            if obj.AWGorSRSswitch
                MW = [MW, {'AWG'}];
            end
            
            % ---- how long the microwave switch stays open -------------------------------
            % See the mwGateMode property comment for why this is a knob. The hold is what the AWG
            % actually needs: the zero pad, the pi pulse, and the pad the pi waveform carries,
            % plus a margin. minWaveformDuration is 0.15 us on this instrument, so the hold comes
            % out to 0.435 us against a tau that reaches 4900.
            minWaveformTime = obj.AWG.minWaveformDuration;
            holdUs = obj.mwHoldUs;
            if isempty(holdUs)
                % CEILED TO A WHOLE NANOSECOND, and both halves of that matter.
                %
                % Nanosecond, because the pulse generator streams durations as integer
                % nanoseconds -- PulseStreamerNewClass.sendToHardware does P(p.duration*1e3, ...)
                % -- and piTime is on the AWG's 62.5 ps grid, so 0.3 + piTime + 0.1 lands on a
                % fraction of a nanosecond (0.44475 us for a 44.75 ns pi). Every other duration
                % this class hands the pulse generator is already a whole nanosecond; this was
                % the one that would not have been, and only in 'pulse' mode, which is the mode
                % nobody has run yet.
                %
                % CEIL rather than round, because the failure is one-sided: a hold SHORTER than
                % the AWG's playback clips the pi pulse, which costs contrast at every tau equally
                % and so reads as a calibration error rather than a timing one. Rounding up costs
                % at most a nanosecond of extra gate.
                holdUs = ceil((2*minWaveformTime + obj.piTime + 0.1) * 1000) / 1000;
            end
            tauMax = max(obj.tau);              % max(), NOT tau(end)
            tauMin = min(obj.tau);
            % max(obj.tau) rather than obj.tau(end), the same latent bug that was fixed in
            % AWGExpEcho. Nothing sorts the tau list, and (1:3:70).^2 only happens to be
            % ascending. perform overwrites this duration at every point, so an unsorted list
            % would not break the run -- but repairSequence merges by nickname and duration
            % arithmetic below depends on the initial value, and 'the longest tau' is what this
            % number is meant to be in every line that reads it.

            usePulseGate = strcmp(obj.mwGateMode, 'pulse');
            if usePulseGate && holdUs >= tauMin
                obj.sendError(sprintf(['mwGateMode = ''pulse'' needs the shortest tau (%.4f us) to ' ...
                    'exceed the gate hold (%.4f us), because the wait is split into a held part ' ...
                    'and a closed part and the closed part cannot be negative. Either drop the ' ...
                    'taus below %.4f us or shorten mwHoldUs -- but read its comment first, a hold ' ...
                    'shorter than the AWG playback clips the pi pulse.'], tauMin, holdUs, holdUs))
            end

            S = Sequence;
            for k = 1:1+obj.doubleMeasurement
                switch obj.mwGateMode
                    case 'pulse'
                        % Open for the pulse, then CLOSED for the rest of the wait. Two events
                        % rather than one, and the nickname 'tau' stays on the long one so
                        % perform's changeSequence call still finds it -- it now sets
                        % tau - hold, which is why the hold has to be known here too.
                        S.addEvent(obj.AWG.TRIG_DURATION,   [triggers(:)', MW(:)']);
                        S.addEvent(holdUs,                    MW,                 'mwHold');
                        S.addEvent(tauMax - holdUs,           '',                 'tau');
                        S.addEvent(obj.lastDelay,           '',                 'lastDelay');
                    case 'all'
                        % Open through the wait AND the last delay. With constantTime = true the
                        % two sum to lastDelay + max(tau) at every point, so the gate-open time
                        % is CONSTANT across the sweep -- which is the whole reason this mode
                        % exists. Its signature is an offset with no slope.
                        S.addEvent(obj.AWG.TRIG_DURATION,   [triggers(:)', MW(:)']);
                        S.addEvent(tauMax,                  MW,                 'tau');
                        S.addEvent(obj.lastDelay,           MW,                 'lastDelay');
                    otherwise   % 'tau'
                        % The historical sequence, unchanged.
                        S.addEvent(obj.AWG.TRIG_DURATION,   [triggers(:)', MW(:)']);
                        S.addEvent(tauMax,                  MW,                 'tau');
                        S.addEvent(obj.lastDelay,           '',                 'lastDelay');
                end

                S.addEvent(obj.detectionDuration,       {'greenLaser', 'detector'});                    % Detection
                S.addEvent(initDuration,                'greenLaser');                                  % Initialization
                S.addEvent(obj.referenceDetectionDuration,...
                    {'greenLaser', 'detector'});                                                        % Reference detection
                S.addEvent(1,                           'greenLaser')
            end

            % ---- record and announce the dose -------------------------------------------
            switch obj.mwGateMode
                case 'pulse'; obj.mwGateOpenUs = obj.AWG.TRIG_DURATION + holdUs;
                case 'all';   obj.mwGateOpenUs = obj.AWG.TRIG_DURATION + tauMax + obj.lastDelay;
                otherwise;    obj.mwGateOpenUs = obj.AWG.TRIG_DURATION + tauMax;
            end
            needUs = obj.AWG.TRIG_DURATION + 2*minWaveformTime + obj.piTime;
            fprintf('T1 gate: mwGateMode = ''%s''.\n', obj.mwGateMode);
            fprintf('         The switch is open for %.3f us per repeat-half at the longest tau,\n', ...
                    obj.mwGateOpenUs);
            fprintf('         against %.3f us that the AWG actually plays. Excess %.1f us.\n', ...
                    needUs, obj.mwGateOpenUs - needUs);
            if strcmp(obj.mwGateMode, 'tau')
                fprintf('         That excess is PROPORTIONAL TO TAU, so anything leaking adds a\n');
                fprintf('         constant rate to 1/T1 and is invisible in the shape.\n');
                fprintf('         echoGateTest bounded the gate out to 240 us; this is %.0fx that.\n', ...
                        (obj.mwGateOpenUs - needUs)/240);
                fprintf('         Run t1GateTest before quoting a T1.\n');
            elseif strcmp(obj.mwGateMode, 'all')
                if ~obj.constantTime
                    fprintf('         *** ''all'' with constantTime = false does NOT give a constant dose:\n');
                    fprintf('         lastDelay is fixed, so the gate-open time still tracks tau and this\n');
                    fprintf('         arm loses the flat signature that makes it a separate measurement.\n');
                    fprintf('         Set constantTime = true.\n');
                else
                    fprintf('         constantTime is on, so this dose is the SAME at every tau -- an\n');
                    fprintf('         offset with no slope, which is what makes it an independent\n');
                    fprintf('         measurement of the leak rate.\n');
                end
            end

            if ~obj.doubleMeasurement
                fprintf(['         doubleMeasurement is OFF, so there is one half and it now carries ' ...
                         'the PI pulse\n         (it used to carry the no-pi one and measure ' ...
                         'nothing). signalParam therefore holds\n         the |1> relaxation ' ...
                         'despite being labelled |0>, and signalParam2 is empty.\n']);
            end

            if obj.changeFlag
                obj.LoadAWG;
            end

            obj.prepareInternal(S)

            % The hold has to reach perform, which subtracts it from every tau it writes.
            % Stored rather than recomputed there so the two cannot disagree.
            obj.mwHoldUsed = holdUs;
            obj.usePulseGate = usePulseGate;

            % Set parameter, for saving
            obj.mCurrentXAxisParam.value = obj.tau;
            obj.topParam.value = obj.tau;
            obj.isPlotAlternateAvailable = obj.doubleMeasurement;
        end
    
        
        function LoadAWG(obj)
            obj.AWG.ResetStoredData;
            obj.AWG.SampleRate = obj.AWG.MAX_SAMPLE_RATE;
            % activeChannels, NOT 1:nChannels. With nChannels = 1 the lab convention is
            % that MWChannel selects WHICH physical channel is live, and activeChannels
            % resolves that to [1] or [2]; 1:nChannels is always [1]. So a single-channel T1
            % on channel 2 armed channel 1 and drove nothing, while the run looked healthy
            % and simply returned no contrast. Echo and Ramsey already do this correctly.
            obj.AWG.workChannels = obj.activeChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            minWaveformTime = obj.AWG.minWaveformDuration;

            % ---- piTime onto the AWG's sample grid -------------------------------------
            % AddWaveformByFunction refuses a duration that is not a whole number of samples:
            %
            %     The parameter 'duration' not divided by the AWG sample time. Fix it or change
            %     the parameter 'toRound' to 'True'
            %
            % and every T1 run was ending there, because a piTime that came out of calibPulseTimes
            % is a fitted number and lands wherever the fit put it. At 16 GS/s the grid is 62.5 ps,
            % so an arbitrary value misses it essentially always.
            %
            % AWGExpT1 was the only one of the five that neither rounded nor asked to be rounded.
            % AWGExpEcho does exactly the line below to halfPiTime, and AWGExpRamsey passes
            % toRound = 1 into the waveform call. This does BOTH: rounding the property is the
            % substantive fix -- it means the value stored, saved and reported is the value the
            % instrument actually played -- and toRound below is then a no-op guard against a
            % floating-point residue in a sum.
            %
            % The physical cost is nil. Rounding a ~45 ns pulse to the nearest 62.5 ps moves it by
            % at most 0.07%, far inside what calibPulseTimes can resolve. What would NOT be nil is
            % leaving obj.piTime disagreeing with the waveform, so the round is written back.
            %
            % Written back only when it actually moves, because set.piTime raises changeFlag and
            % an unconditional assignment here would mark the experiment dirty on every load --
            % which forces a full waveform rebuild on the next prepare for no reason.
            piTimeOnGrid = obj.RoundDurationByFrequencies(obj.piTime, obj.AWG.SampleRate);
            if abs(piTimeOnGrid - obj.piTime) > 1e-12
                fprintf(['T1: piTime %.6f -> %.6f us, onto the AWG sample grid (%.4f ns). ' ...
                         'A fitted\n    pulse time does not land on it by itself, and ' ...
                         'AddWaveformByFunction refuses\n    a duration that is not a whole ' ...
                         'number of samples.\n'], ...
                        obj.piTime, piTimeOnGrid, 1e3/obj.AWG.SampleRate);
                obj.piTime = piTimeOnGrid;
            end
            
            initZeroIndex = zeros(1, obj.nChannels);
            piWaveformsIndex = zeros(1, obj.nChannels);
            noPiWaveformsIndex = zeros(1, obj.nChannels);   % was grown inside the loop
            
            prePulseTime = minWaveformTime * ones(1, obj.nChannels);
            phase = zeros(1, obj.nChannels);

            % creating waveforms:
            A = obj.AWG.MaxNormPowerAllowed;
            w = 2*pi*obj.frequency;
            % '~2 per channel' was a guess and it was low: the loop adds THREE per channel --
            % the leading zero pad, the pi pulse, and a zeroed copy of it for the reference.
            % T1's whole waveform set is tau-INDEPENDENT, which is why there are only a handful
            % and why perform assigns the sequence once outside the tau loop.
            obj.AWG.progressBegin(sprintf(['T1 waveforms: 3 per channel (leading zero pad, ' ...
                'pi pulse, zeroed pi for the reference) x %d channel(s). None depend on tau -- ' ...
                'the wait is a pulse-generator duration, not an AWG waveform.'], obj.nChannels), ...
                3*obj.nChannels, 'waveform');
            
            % the resolution of the time delay is the sample rate, so we add also aditional phase delay
            timeDelay = obj.RoundDurationByFrequencies(abs(obj.timeDelay), obj.AWG.SampleRate, -1);
            phaseDelay = w * (abs(obj.timeDelay) - timeDelay);
            if obj.timeDelay > 0
                prePulseTime(1) = minWaveformTime + timeDelay;
                phase(1) = phaseDelay;
            elseif obj.nChannels >= 2
                % GUARDED on nChannels. Unguarded, this branch runs whenever timeDelay <= 0
                % -- which is always now, since level 4.3 sets it to 0 -- and assigning
                % index 2 of a one-element array GROWS it to two elements. Harmless today
                % because the loop below only reads 1:nChannels, but it is the kind of
                % silent reshape that makes a later indexing bug impossible to find.
                prePulseTime(2) = minWaveformTime + timeDelay;
                phase(2) = phaseDelay;
            end

            % THE CALIBRATED RELATIVE PHASE. See AWGExpEcho for the full reasoning; briefly,
            % channelsDifference-style phase comes from timeDelay alone, level 4.3 sets
            % timeDelay = 0, so without this line the two channels are driven with zero
            % relative phase while piTime in the json was measured with the calibrated phase
            % applied -- and the pi pulse is then not a pi pulse.
            phase = phase + obj.channelPhaseRadians();
            for chan = 1:obj.nChannels
                % The 4th argument is toRound. With piTime already on the grid these are no-ops;
                % they are here so that a residue from adding two on-grid numbers, or a future
                % duration that is not pre-rounded, degrades to a 62.5 ps adjustment instead of
                % an error four frames deep. AWGExpRamsey passes it for the same reason.
                % progressStep() takes a SAMPLE COUNT, not a step count -- it increments the item
                % counter by one on every call regardless -- so it is called once per waveform.
                % progressBegin promised 3*nChannels items and this is what delivers them.
                f = @(t) 0 * t;
                initZeroIndex(chan) = obj.AWG.AddWaveformByFunction(f, prePulseTime(chan), 0, 1);
                obj.AWG.progressStep();
                f = @(t) A(chan) * cos(w*t - phase(chan)) .* (t < obj.piTime);
                piWaveformsIndex(chan) = obj.AWG.AddWaveformByFunction(f, obj.piTime + minWaveformTime, 0, 1);
                obj.AWG.progressStep();
                f = @(t) A(chan) * cos(w*t - phase(chan)) .* (t < obj.piTime) * 0;
                noPiWaveformsIndex(chan) = obj.AWG.AddWaveformByFunction(f, obj.piTime + minWaveformTime, 0, 1);
                obj.AWG.progressStep();
            end
            
            obj.AWG.progressEnd();
            obj.AWG.LoadWaveform();
            
            % creating sequences:
            obj.indexSeq = zeros(1, obj.nChannels);
            if obj.doubleMeasurement
                seqTxt = 'T1 sequences: one per channel, 4 steps each (zero, no-pi, zero, pi).';
            else
                seqTxt = 'T1 sequences: one per channel, 2 steps each (zero, pi) -- single measurement.';
            end
            obj.AWG.progressBegin(seqTxt, obj.nChannels, 'sequence');
            for chan = 1:obj.nChannels
                if obj.doubleMeasurement
                    % Two trigger-gated blocks. Trigger 1 plays [zero, noPi] -- NO pulse, the
                    % spin stays in |0> -- and trigger 2 plays [zero, pi]. That order is why
                    % detection rows 1-2 are the |0> relaxation and rows 3-4 the |1> one, which
                    % is the opposite of Echo and Ramsey; see the constructor.
                    obj.indexSeq(chan) = obj.AWG.AddSequence( ...
                        [initZeroIndex(chan), noPiWaveformsIndex(chan), initZeroIndex(chan), piWaveformsIndex(chan)], ...
                        [1 1 1 1], [1 0 1 0]);
                else
                    % SINGLE MEASUREMENT: the pi block, not the no-pi one.
                    %
                    % This branch is new, and what it replaces did not work. The four-step
                    % sequence above was built regardless, while prepare only emits ONE
                    % trigger-bearing half when doubleMeasurement is false -- so the run fired
                    % trigger 1, played [zero, noPi], and never reached the pi pulse at all. The
                    % result was a T1 measurement of the |0> state: a flat line at full
                    % fluorescence, with error bars, that looks exactly like an excellent sample.
                    % Nothing errored and nothing warned.
                    %
                    % The single measurement that is worth taking is the |1> relaxation, so that
                    % is what this builds. Note the consequence for naming: with doubleMeasurement
                    % off, signalParam carries |1> even though it is labelled |0>, because the
                    % labels describe the two-half ordering. prepare prints this.
                    obj.indexSeq(chan) = obj.AWG.AddSequence( ...
                        [initZeroIndex(chan), piWaveformsIndex(chan)], [1 1], [1 0]);
                end
                obj.AWG.progressStep();
            end
            obj.AWG.progressEnd();
            obj.AWG.LoadSequence();
        end
        
        
        function perform(obj)
            %%% Initialization
            
            % Devices (+ Tracker)
            pg = getObjByName(PulseGenerator.NAME);
            spcm = getObjByName(Spcm.NAME);
            tracker = getObjByName(Tracker.NAME);
                if isempty(tracker); throwBaseObjException(Tracker.Name); end
            
            % Some magic numbers
            maxLastDelay = obj.lastDelay + max(obj.tau);


            % obj.activeChannels, NOT the literals 1 and 2. The sequence list is built 1-based
            % (indexSeq(1) is the first ACTIVE channel), and activeChannels is what maps that onto
            % the physical outputs: [1 2] at nChannels = 2, but a single 2 when nChannels = 1 with
            % MWChannel = 'MW2'.
            %
            % The old form hardcoded output 1 for the first sequence. At two channels that is right
            % by coincidence. At one channel on output 2 it armed the output that had no waveform
            % loaded -- LoadAWG sets workChannels from activeChannels, so the waveforms went to 2
            % and the sequence to 1 -- and the run then looked healthy and returned no contrast.
            % The comment at the top of LoadAWG describes that exact trap for workChannels; this
            % is the other half of it. Echo, Ramsey and Rabi all index activeChannels here.
            for i = 1:obj.nChannels
                obj.AWG.assignSequence(obj.indexSeq(i), obj.activeChannels(i));
            end
            obj.AWG.Run;

            %%% Run - Go over all tau's, in random order
            for t = randperm(length(obj.tau))
                if obj.currIter <= size(obj.signal, 3) && ...
                        sum(obj.signal(:, t, obj.currIter)) ~= 0
                    continue
                end
                success = false;
                
                for trial = 1 : 5
                    if obj.checkEmergencyStop()
                        return;
                    end
                    try
                        % In 'pulse' mode the wait is split: an 'mwHold' event with the switch
                        % open, then the 'tau' event with it closed. The nickname 'tau' still
                        % names the long one, so this call does not change -- except that what
                        % it must write is the REMAINDER. Writing the full tau here would
                        % lengthen every point by the hold, which is 0.435 us against taus up to
                        % 4900 and so would never have been noticed at the long end and would
                        % have been a 44% error at tau = 1.
                        if obj.usePulseGate
                            pg.changeSequence('tau', 'duration', obj.tau(t) - obj.mwHoldUsed);
                        else
                            pg.changeSequence('tau', 'duration', obj.tau(t));
                        end
                        if obj.constantTime
                            pg.changeSequence('lastDelay', 'duration', maxLastDelay - obj.tau(t));
                        end
                        
                        data = obj.getRawData(pg, spcm);
                        [sig, sterr] = obj.processData(data);
                        
                        obj.signal(:, t, obj.currIter) = sig;
                        obj.sterr(:, t, obj.currIter) = sterr;
                        
                        success = true;
                        obj.currParamIter = obj.currParamIter + 1;
                        sendEventParamIterationDone(obj);
                        
                        if obj.isTracking
                            if obj.checkEmergencyStop()
                                return;
                            end
                            isTrackingNeeded = tracker.compareReference(...
                                sig(2), sterr(2), ... % 2 for reference
                                Tracker.REFERENCE_TYPE_KCPS, obj.trackThreshhold);
                            if isTrackingNeeded
                                tracker.trackUsing(TrackablePosition.NAME)
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
            
            % Saving results in the Experiment parameters
            S1 = squeeze(obj.signal(1, :, 1:obj.currIter));
            S1sterr = squeeze(obj.sterr(1, :, 1:obj.currIter));
            S2 = squeeze(obj.signal(2, :, 1:obj.currIter));
            S2sterr = squeeze(obj.sterr(2, :, 1:obj.currIter));
            
            [value, sterr] = getRatioDistributionValues(obj, S1, S2, S1sterr, S2sterr);
            obj.signalParam.value = value;       % rows 1-2: the NO-pi half, i.e. |0>
            obj.signalParam.sterr = sterr;
            
            if obj.doubleMeasurement
                S3 = squeeze(obj.signal(3, :, 1:obj.currIter));
                S3sterr = squeeze(obj.sterr(3, :, 1:obj.currIter));
                S4 = squeeze(obj.signal(4, :, 1:obj.currIter));
                S4sterr = squeeze(obj.sterr(4, :, 1:obj.currIter));
                
                [value, sterr] = getRatioDistributionValues(obj, S3, S4, S3sterr, S4sterr);
                obj.signalParam2.value = value;  % rows 3-4: the pi half, i.e. |1>
                obj.signalParam2.sterr = sterr;
                obj.checkProjectionLabels();
            else    
                obj.signalParam2.value = [];
            end
        end

        function checkProjectionLabels(obj)
            %CHECKPROJECTIONLABELS Test the |0>/|1> assignment against the DATA, every run.
            %
            % The assignment in the constructor comes from reading LoadAWG's waveform order, and
            % a reading is not a measurement. This turns it into one, cheaply: in a T1 the |1>
            % trace is the one that MOVES with tau and the |0> trace is the flat control, so
            % whichever trace varies more is |1>. If that disagrees with the labels, something
            % is wrong -- either this comment or the sequence -- and it should say so at the
            % moment the data exists rather than being discovered in a plot months later.
            %
            % Deliberately quiet unless the two differ by a clear factor. Early in a run, or on a
            % sample whose T1 is far longer than the tau range, BOTH traces are flat and the
            % comparison means nothing; demanding a factor of two keeps it from crying wolf.
            try
                y0 = obj.signalParam.value(:)';
                y1 = obj.signalParam2.value(:)';
                if numel(y0) < 5 || numel(y1) ~= numel(y0); return; end
                v0 = max(y0) - min(y0);
                v1 = max(y1) - min(y1);
                if ~isfinite(v0) || ~isfinite(v1) || v1 <= 0; return; end
                if v0 > 2*v1
                    fprintf(['\n*** T1 PROJECTION LABELS: the trace labelled |0> varies %.1fx more ' ...
                             'with tau than the\n    one labelled |1>, and it should be the other ' ...
                             'way round -- |0> is the flat control.\n    Either LoadAWG''s waveform ' ...
                             'order changed, or the labels set in the constructor are\n    wrong. ' ...
                             'Do not read which curve is which off the legend until this is ' ...
                             'settled.\n\n'], v0/max(v1, eps));
                end
            catch
                % A label check must never be the thing that kills a finished run.
            end
        end
        
        function wrapUp(obj)
            % Things that need to happen when the experiment is done; a
            % counterpart for obj.prepare.
            % In the future, it will also analyze results and fit from it
            % the coherence time.
            
            obj.wrapUpInternal()
        end
        
        function dataParam = alternateSignal(obj)
            % Returns alternate view ("referenced") of the data, as an
            % ExpParam.
            
            % signalParam is |0> and signalParam2 is |1> -- see the constructor. The variable
            % names here used to say the opposite, which made the subtraction below read as
            % |0> - |1> while computing |1> - |0>.
            %
            % The orientation is kept as it was INTENDED: |0> minus |1>, which is positive after
            % a fresh pi pulse and decays to zero. That is the curve worth looking at, and it is
            % also the one whose zero crossing is meaningful. Flipping it would have made every
            % saved alternate view change sign relative to the ones taken before today.
            N0 = obj.signalParam.value;
            N1 = obj.signalParam2.value;

            value = N0 - N1;
            sterr = sqrt(obj.signalParam.sterr.^2 + obj.signalParam2.sterr.^2);
            dataParam = ExpResultDoubleVector('FL', value, sterr, 'Normalized', obj.NAME);
        end
    end
end

