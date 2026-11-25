classdef AWGExpLastDelayDuration < ExperimentAWG
    
    properties (Constant)
        NAME = 'AWG last delay';
    end
    
    properties
        % All of the parameters must have a set function, and when they are changed, the parameter obj.changeFlag needs to be set to true
        % (unless the change is minor and the experiment can continue running);
        tau             % in us
        totalZero
        smallOnDelay    % in us
        smallOffDelay   % in us
        dontLoadAWG

        % ---- isolation switches for hunting detector artefacts (added 2026-11-04) ----
        % Defaults reproduce the normal measurement.
        mwOff = false       % true: the AWG plays ZEROS. Switches and triggers unchanged, so
                            %   whatever is left is not caused by microwave power.
        gapSwitch = true    % true: MW switches stay open through the gap and close when the
                            %   window opens (the normal sequence). false: they close at the
                            %   END OF THE BURST, so the switch edge is no longer at the window.
                            %   That can clip the tail of the AWG burst: DIAGNOSTIC ONLY (off
                            %   resonance), never for measuring the step.
        noTrigger = false   % true: the AWG is never triggered (no MW at all). Together with
                            %   refSwitch the two halves are identical in every channel, so any
                            %   remaining "contrast" comes from the readout, not the sequence.
        refSwitch = false   % (balanced sequence only) true: the REFERENCE half opens and closes the MW switches exactly
                            %   like the signal half, with the AWG not triggered (silent). Any
                            %   switch-correlated detector artefact is then common to both
                            %   halves and cancels in the ratio.
    end

    methods
        function obj = AWGExpLastDelayDuration()
            obj@ExperimentAWG(AWGExpLastDelayDuration.NAME);
            obj.parameterName = 'taus';
            
            obj.repeats = 10000;
            obj.averages = 100;
            obj.frequency = 2760;           % in MHz
            obj.amplitude = -10 + [0 0];    % in dBm
            obj.tau = 0.005:0.02:1;       % us
            obj.constantTime = true;
            obj.smallOnDelay = 0;         % us
            obj.smallOffDelay = 0; %1.75;        % us
            obj.lastDelay = Experiment.DEFAULT_LAST_DELAY;
            obj.dontLoadAWG = 0;
            % BALANCED by default (2026-11-04). Experiment's default is 0, which gives a reference
            % window straight after continuous laser while the signal window follows TRIG +
            % piTime + gap of dark. Fluorescence after dark is lower than under steady light
            % (0.14% after ~1 us, ~0.7% after 100 us on Setup 7), so the unbalanced scan showed a
            % slow 'contrast' rise with the gap that had nothing to do with the microwave.
            obj.balancedSequence = 1;
            
            obj.detectionDuration = 0.25;               % detection window, in us
            obj.referenceDetectionDuration = 5;         % in us. Detection duration of the reference read
            obj.laserInitializationDuration = 10;       % laser initialization in pulsed experiments
            
            obj.mCurrentXAxisParam = ExpParamDoubleVector('Time', [], [], StringHelper.MICROSEC, obj.NAME);
            obj.signalParam = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME);
        end
    end

%% Setters
    methods
        function set.mwOff(obj, newVal)
            obj.mwOff = logical(newVal);
            obj.changeFlag = true;          % the AWG waveform changes
        end

        function set.gapSwitch(obj, newVal)
            obj.gapSwitch = logical(newVal);
        end

        function set.noTrigger(obj, newVal)
            obj.noTrigger = logical(newVal);
        end

        function set.refSwitch(obj, newVal)
            obj.refSwitch = logical(newVal);
        end

        function set.tau(obj, newVal) % newVal is in us
            checkTimeVector(obj, newVal)
            % If we got here, then newVal is OK.
            obj.tau = newVal;
            obj.changeFlag = true;
        end
        
    end
    
    methods
        function totalParamNum = getTotalNumberOfParams(obj)
            totalParamNum = length(obj.tau);
        end
    end
    
%% Overridden from Experiment
    methods (Access = protected)
        
        function prepare(obj)
            % Initialize devices (SPCM, PulseGenerator, etc.)
            
            % Sequence
            %%% Useful parameters for what follows
            initDuration = obj.laserInitializationDuration-obj.detectionDuration;
            obj.detectionPeriodsPerRepeat = 2;
            obj.runsPerPerform = 1;
            MWChannel = obj.MWChannel;
            checkDetectionDuration(obj, obj.detectionDuration); % The mode might have changed; before running, we check
            
            %%% Creating the sequence
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
            
            % Switch states for the gap and for the reference half; see the isolation switches.
            gapMW = '';
            if obj.gapSwitch; gapMW = MW; end
            refMW = ''; refGapMW = '';
            if obj.refSwitch
                refMW = MW;
                if obj.gapSwitch; refGapMW = MW; end
            end

            trigCh = [triggers(:)', MW(:)'];
            if obj.noTrigger; trigCh = MW; end          % isolation switch: AWG never plays

            S = Sequence;
            S.addEvent(obj.AWG.TRIG_DURATION,           trigCh);                                        % AWG trigger
            S.addEvent(obj.piTime,                      MW)               % MW
            S.addEvent(obj.lastDelay,                   gapMW,                      'lastDelay');       % Last delay

            S.addEvent(obj.detectionDuration,       {'greenLaser', 'detector'});                    % Detection
            S.addEvent(initDuration,                'greenLaser');                                  % Initialization
            if ~obj.balancedSequence
                S.addEvent(obj.referenceDetectionDuration,...
                    {'greenLaser', 'detector'});                                                        % Reference detection
                S.addEvent(1,                           'greenLaser');                                  % Initialization
            else
                % No AWG trigger here, ever: the reference must get no microwave. With
                % refSwitch the MW switches still open and close exactly as in the signal half.
                S.addEvent(obj.AWG.TRIG_DURATION,           refMW);                     % AWG trigger
                % Same dark time as the signal half (TRIG + piTime + lastDelay). This used to be
                % obj.tau(end) -- the LONGEST GAP -- so the reference window followed
                % tau(end) - piTime more dark than the signal window, and any dark-time effect
                % (laser/AOM turn-on, charge state) read as contrast. Fixed 2026-11-04.
                S.addEvent(obj.piTime,                      refMW,                      'MW');               % MW
                S.addEvent(obj.lastDelay,                   refGapMW,                   'lastDelay');       % Last delay
                S.addEvent(obj.referenceDetectionDuration,...
                    {'greenLaser', 'detector'});                                                        % Reference detection
                S.addEvent(initDuration,                'greenLaser');                                  % Initialization
            end
            
            obj.prepareInternal(S)
            if (obj.changeFlag || obj.restartFlag) && (~obj.dontLoadAWG)
                obj.LoadAWG;
            end
            
            % Set parameter, for saving
            obj.mCurrentXAxisParam.value = obj.tau;
        end
        
        function LoadAWG(obj)
            obj.AWG.ResetStoredData;
            obj.AWG.SampleRate = obj.AWG.MAX_SAMPLE_RATE;
            obj.AWG.workChannels = 1:obj.nChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            minWaveformTime = obj.AWG.minWaveformDuration;
            
            initZeroIndex = zeros(1, obj.nChannels);
            piWaveformsIndex = zeros(1, obj.nChannels);
            
            prePulseTime = minWaveformTime * ones(1, obj.nChannels);
            phase = zeros(1, obj.nChannels);

            % creating waveforms:
            A = obj.AWG.MaxNormPowerAllowed;
            w = 2*pi*obj.frequency;
            fprintf('Loading ~%d waveforms: ', 2 * obj.nChannels);
            
            % the resolution of the time delay is the sample rate, so we add also aditional phase delay
            timeDelay = obj.RoundDurationByFrequencies(abs(obj.timeDelay), obj.AWG.SampleRate, -1);
            phaseDelay = w * (abs(obj.timeDelay) - timeDelay);
            if obj.timeDelay > 0
                prePulseTime(1) = minWaveformTime + timeDelay;
                phase(1) = phaseDelay;
            else
                prePulseTime(2) = minWaveformTime + timeDelay;
                phase(2) = phaseDelay;
            end
            amp = A * double(~obj.mwOff);          % mwOff: same waveform, all zeros
            for chan = 1:obj.nChannels
                f = @(t) 0 * t;
                initZeroIndex(chan) = obj.AWG.AddWaveformByFunction(f, prePulseTime(chan), 0);
                f = @(t) amp(chan) * cos(w*t - phase(chan)) .* (t < obj.piTime);
                piWaveformsIndex(chan) = obj.AWG.AddWaveformByFunction(f, obj.piTime + minWaveformTime, 0);
            end
            
            fprintf('\n');
            obj.AWG.LoadWaveform();
            
            % creating sequences:
            obj.indexSeq = zeros(1, obj.nChannels);
            % ONE triggered block: wait for trigger, play the burst. The sequence used to be
            % [zero, noPi(zeros), zero, pi] with triggers [1 0 1 0] -- two triggered blocks -- while
            % the pulse sequence triggers the AWG only once per repeat (signal half). So the
            % sequencer alternated: every OTHER repeat played zeros and got no microwave, which
            % halved the measured contrast. Fixed 2026-11-04.
            for chan = 1:obj.nChannels
                obj.indexSeq(chan) = obj.AWG.AddSequence([initZeroIndex(chan), piWaveformsIndex(chan)], [1 1], [1 0]);
            end
            obj.AWG.LoadSequence();
        end
       
        function perform(obj)
            %%% Initialization
            
            % Devices (+ Tracker)
            pg = getObjByName(PulseGenerator.NAME);
            spcm = getObjByName(Spcm.NAME);
            tracker = getObjByName(Tracker.NAME);
                if isempty(tracker); throwBaseObjException(Tracker.Name); end

            obj.AWG.assignSequence(obj.indexSeq(1),1)
            if obj.nChannels == 2
                obj.AWG.assignSequence(obj.indexSeq(2),2);
            end
            obj.AWG.Run;
            
            k = 0; %added by rotem 18.4.21
            f1 = []; %added by rotem 18.4.21
            
            for t = randperm(length(obj.tau))
                if obj.currIter <= size(obj.signal, 3) && ...
                        sum(obj.signal(:, t, obj.currIter)) ~= 0
                    continue
                end
                success = false;
                
                f1 = [f1, t]; %added by rotem 18.4.21
                k = k + 1; %added by rotem 18.4.21
                
                for trial = 1 : 10
                    if obj.checkEmergencyStop()
                        return;
                    end
                    try
                        pg.changeSequence('lastDelay', 'duration', obj.tau(t));
                        
                        data = obj.getRawData(pg, spcm);
                        [sig, sterr] = obj.processData(data);
                        obj.signal(:, t, obj.currIter) = sig;
                        obj.sterr(:, t, obj.currIter) = sterr;
                        
                        success = true;
                        obj.currParamIter = obj.currParamIter + 1;
                        if obj.currParamIter == 1 || mod(obj.currParamIter, 5) == 0 || obj.currParamIter == obj.totalNumberOfParams
                            if obj.currParamIter > obj.totalNumberOfParams
                                '';
                            end
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
            
            % for time tagger setup
            if length(data) == obj.detectionPeriodsPerRepeat*obj.repeats*obj.getTotalNumberOfParams %if true we're proccessing data at the end of the average
                for k = 1:length(obj.tau)
                    [obj.signal(:, f1(k), obj.currIter), obj.sterr(:, f1(k), obj.currIter)] = obj.processData(data(obj.detectionPeriodsPerRepeat*obj.repeats*(k-1)+(1:obj.detectionPeriodsPerRepeat*obj.repeats)));
                end
            end
            
            % Saving results in the Experiment parameters
            S1 = squeeze(obj.signal(1, :, 1:obj.currIter));
            S1sterr = squeeze(obj.sterr(1, :, 1:obj.currIter));
            S2 = squeeze(obj.signal(2, :, 1:obj.currIter));
            S2sterr = squeeze(obj.sterr(2, :, 1:obj.currIter));
            
            [value, sterr] = getRatioDistributionValues(obj, S1, S2, S1sterr, S2sterr);
            obj.signalParam.value = value;
            obj.signalParam.sterr = sterr;         
        end
        
        function wrapUp(obj)
            % Things that need to happen when the experiment is done; a
            % counterpart for obj.prepare.
            % In the future, it will also analyze results.
            
            obj.wrapUpInternal()
        end
        
        function dataParam = alternateSignal(obj) %#ok<MANU>
            % Returns alternate view ("referenced") of the data, as an
            % ExpParam, if possible. If not, it returns an empty variable.
            dataParam = [];
        end
    end

end
