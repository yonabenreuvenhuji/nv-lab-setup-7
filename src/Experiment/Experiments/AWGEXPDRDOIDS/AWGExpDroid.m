classdef AWGExpDroid < ExperimentAWG
    %EXPECHO Echo experiment

    properties (Constant)
        NAME = 'AWG Droid';
    end
    
    properties
        % All of the parameters must have a set function, and when they are
        % changed, the parameter obj.changeFlag needs to be set to true
        % (unless the change is minor and the experiment can continue
        % running);
        amplitude           % in dBm
        
        tau                 % in us 
        pulsePhase          % Phase of the pi pulses. Can be vector (like 'XYYX' for XY4) and the phases will repeats the order.
        cycles              % Number of pi pulses in each sequence. 1 for regular echo.
        
        constantTime        % logical
        doubleMeasurement   % logical
        pulseSequence
        pulseAngles
        directions
        mormSegmentDuration
        segmentDurs
        initPulseDirection
        initPulseAngle
        nRepetitions
        EchoLike
        piPulseLocs
    end

    properties (Dependent)
        pulseTrainTime
    end


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = AWGExpDroid()
            obj@ExperimentAWG(AWGExpDroid.NAME);
            obj.parameterName = 'taus';
            
            % Set properties inherited from Experiment
            obj.repeats = 5000;
            obj.averages = 1000;
            obj.frequency = 3000;           % in MHz
            obj.amplitude = -3;            % in dBm
            obj.nRepetitions = 1:1:3;                   
            obj.halfPiTime = 0.0170;         % in us
            obj.piTime = 0.0346;
            obj.lastDelay = Experiment.DEFAULT_LAST_DELAY;
            obj.AWGorSRSswitch = 0;
            obj.isTracking = 0;

            obj.pulsePhase = 'Y';           % Phase of the pi pulses. Can be vector (like 'XYYX' for XY4) and the phases will repeats the order.
            obj.cycles = 1;                 % Number of pi pulses in each sequence. 1 for regular echo.
            obj.nChannels = 2;
            
            phi = (sqrt(5) + 1)/ 2;
            obj.pulseAngles = [1]*pi;

            obj.directions = [0,1,0]; 
            obj.mormSegmentDuration = 5;        % us
            obj.segmentDurs = [1,1]*obj.mormSegmentDuration;
            obj.initPulseDirection = [1, 0, 0];
            obj.initPulseAngle = pi / 2;
            obj.RabiFreq = 14.12;
            obj.RabiX0 = 0;
            obj.EchoLike = 0;
            obj.piPulseLocs = length(obj.pulseAngles) / 2 + 0.5;

            obj.constantTime = false;       % logical
            obj.doubleMeasurement = true;   % logical
            
            obj.detectionDuration = 0.5;           % detection window, in us
            obj.referenceDetectionDuration = 0.5;     % in us. Detection duration of the reference read
            obj.laserInitializationDuration = 10;   % laser initialization in pulsed experiments in us
            
            obj.tau = sum(obj.segmentDurs)*obj.nRepetitions;
            
            obj.mCurrentXAxisParam = ExpParamDoubleVector('Time', [], [], StringHelper.MICROSEC, obj.NAME);
            obj.topParam = ExpParamDoubleVector(StringHelper.TAU, [], [], StringHelper.MICROSEC, obj.NAME);
            obj.signalParam = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME, '|1>');
            obj.signalParam2 = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME, '|0>');
        end
    end
    
    %% Setters
    methods
        function set.amplitude(obj, newVal) % newVal is in dBm
            checkAmplitude(obj, newVal)
            % If we got here, then newVal is OK.
            obj.amplitude = newVal;
            obj.changeFlag = true;
        end
        
        function set.pulsePhase(obj, newVal)	% newVal in microsec
            checkPhases(obj, newVal)
            % If we got here, then newVal is OK.
            obj.pulsePhase = newVal;
            obj.changeFlag = true;
                end
        
        function set.cycles(obj, newVal)	% newVal in microsec
            checkPositive(obj, newVal)
            % If we got here, then newVal is OK.
            obj.cycles = newVal;
            obj.changeFlag = true;
        end

        function set.constantTime(obj, newVal)
            checkBoolean(obj, newVal)
            % If we got here, then newVal is OK.
            obj.constantTime = newVal;
            obj.changeFlag = true;
        end
        
        function set.doubleMeasurement(obj, newVal)
            checkBoolean(obj, newVal)
            % If we got here, then newVal is OK.
            obj.doubleMeasurement = newVal;
            obj.changeFlag = true;
        end

        function set.pulseSequence(obj, pulseSequence)
            obj.pulseAngles = pulseSequence.angles;
            obj.directions = pulseSequence.axes;
            obj.segmentDurs = pulseSequence.segmentDurs;
            obj.segmentDurs = [obj.segmentDurs, 0];
            obj.changeFlag = true;
        end
        
        function pulseTrainTime = get.pulseTrainTime(obj)
            pulsesDuration = 0;
            for i = 1:length(obj.pulseAngles)
                [durs, ~] = obj.pulseToDurations(obj.directions(i,:), obj.pulseAngles(i));
                pulsesDuration = pulsesDuration + sum(durs);
            end
            pulseTrainTime = sum(obj.segmentDurs) + pulsesDuration;
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
            detectionSequnceDuration = obj.detectionDuration + obj.stabilizationDuration + obj.acquisitionDuration + obj.GIresetDuration;
            obj.detectionPeriodsPerRepeat = 2 * (1 + double(obj.doubleMeasurement));
            obj.tau = sum(obj.segmentDurs)*obj.nRepetitions;
            obj.runsPerPerform = 1;
            MWChannel = obj.MWChannel;
            [durs, ~] = obj.pulseToDurations(obj.initPulseDirection, obj.initPulseAngle);
            totalTime = 2*sum(durs) + obj.pulseTrainTime;         

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
                    MW = {MWChannel, 'MW2'};
                otherwise
                    obj.sendError('What should we do here?')
            end
            if obj.AWGorSRSswitch
                MW = [MW{:}, {'AWG'}];
            end
            
            S = Sequence;
            for k = 1:1+obj.doubleMeasurement
                S.addEvent(obj.AWG.TRIG_DURATION,           [triggers(:)', MW(:)']);   
                S.addEvent(totalTime,                       MW,                         'tau');
                S.addEvent(obj.lastDelay,                   '',                         'lastDelay');   % Last delay

                S.addEvent(obj.detectionDuration,       {'greenLaser', 'detector'});                    % Detection
                S.addEvent(initDuration,                'greenLaser');                                  % Initialization
                S.addEvent(obj.referenceDetectionDuration,...
                                                        {'greenLaser', 'detector'});                    % Reference detection
                S.addEvent(1,                           'greenLaser');                                  % Initialization
            end
            
            if obj.changeFlag
                obj.LoadAWG;
            end

            obj.prepareInternal(S)
            
            % Set parameter, for saving
            obj.mCurrentXAxisParam.value = obj.nRepetitions * obj.pulseTrainTime;
            obj.isPlotAlternateAvailable = obj.doubleMeasurement;
        end
    
        
        function LoadAWG(obj)
            obj.AWG.SampleRate = 15158;
            totalMainPulseTime = obj.RoundDurationByFrequencies(obj.pulseTrainTime, [obj.frequency,obj.AWG.SampleRate]);
            temp_time = obj.pulseTrainTime;
            numNotZeros = sum(obj.segmentDurs ~= 0);
            for i=1:length(obj.segmentDurs)
                if obj.segmentDurs(i)
                    obj.segmentDurs(i) = obj.segmentDurs(i) + (totalMainPulseTime - temp_time) / numNotZeros;
                end
            end
            obj.tau = sum(obj.segmentDurs)*obj.nRepetitions;
            obj.AWG.ResetStoredData;
            obj.AWG.workChannels = 1:obj.nChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            obj.AWG.ChangeSampleRate(totalMainPulseTime);
            
            N = length(obj.nRepetitions);
            initWaveformIndex = zeros(obj.nChannels, 1);
            mainWaveformsIndex = zeros(obj.nChannels, 1);
            finalWaveformIndex = zeros(obj.nChannels, 1);
            finalWaveformRefIndex = zeros(obj.nChannels, 1);
            phaseZero = zeros(obj.nChannels, 1);
            
            A = obj.AWG.MaxNormPowerAllowed;
            w = 2*pi*obj.frequency;
            fprintf('Loading 5 waveforms: ');
            
            minWaveformTime = obj.AWG.minWaveformDuration;
                    
        

            [durs, signs] = obj.pulseToDurations(obj.initPulseDirection, obj.initPulseAngle);
            initTime = minWaveformTime + sum(durs);
            minWaveformTime = minWaveformTime - (initTime - obj.RoundDurationByFrequencies(initTime, obj.AWG.SampleRate, -1));
            
            % the resolution of the time delay is the sample rate, so we add also aditional phase delay
            prePulseTime = minWaveformTime * ones(1,obj.nChannels);
            phase = zeros(1,obj.nChannels);
            timeDelay = obj.RoundDurationByFrequencies(abs(obj.timeDelay), obj.AWG.SampleRate, -1);
            phaseDelay = w * (abs(obj.timeDelay) - timeDelay);
            if obj.timeDelay > 0
                prePulseTime(1) = minWaveformTime + timeDelay;
                phase(1) = phaseDelay;
            else
                prePulseTime(2) = minWaveformTime + timeDelay;
                phase(2) = phaseDelay;
                
            end
            
            % init pulse:
            pulseTime = sum(durs);
            for chan = 1:obj.nChannels
                phaseZero(chan) = w * (prePulseTime(chan) + pulseTime) + phase(chan);
                f = @(t) A(chan) * ((prePulseTime(chan)                <=t & t<prePulseTime(chan)+durs(1))                 * signs(1) .* cos(w*t - phaseZero(chan)) + ...
                                    (prePulseTime(chan)+durs(1)        <=t & t<prePulseTime(chan)+durs(1)+durs(2))         * signs(2) .* cos(w*t + pi/2 - phaseZero(chan)) + ...
                                    (prePulseTime(chan)+durs(1)+durs(2)<=t & t<prePulseTime(chan)+durs(1)+durs(2)+durs(3)) * signs(3) .* cos(w*t - phaseZero(chan)));
                initWaveformIndex(chan) = obj.AWG.AddWaveformByFunction(f, prePulseTime(chan)+pulseTime, 0, 0);
            end
            
            % last pulse:
            [durs, signs] = obj.pulseToDurations(-obj.initPulseDirection, obj.initPulseAngle);
            pulseTime = sum(durs);
            for chan = 1:obj.nChannels
                f = @(t) A(chan) * ((0              <=t & t<durs(1))                 * signs(1) .* cos(w*t - phase(chan)) + ...
                                    (durs(1)        <=t & t<durs(1)+durs(2))         * signs(2) .* cos(w*t + pi/2 - phase(chan)) + ...
                                    (durs(1)+durs(2)<=t & t<durs(1)+durs(2)+durs(3)) * signs(3) .* cos(w*t - phase(chan)));
                finalWaveformIndex(chan) = obj.AWG.AddWaveformByFunction(f, 3*minWaveformTime+pulseTime,  0, 1);
                f = @(t) A(chan) * ((0                       <=t & t<durs(1))                 * signs(1) .* cos(w*t - phase(chan)) + ...
                                    (durs(1)                 <=t & t<durs(1)+durs(2))         * signs(2) .* cos(w*t + pi/2 - phase(chan)) + ...
                                    (durs(1)+durs(2)         <=t & t<durs(1)+durs(2)+durs(3)) * signs(3) .* cos(w*t - phase(chan))+...
                                    (durs(1)+durs(2)+durs(3) <=t & t<durs(1)+durs(2)+durs(3)+obj.piTime) .* cos(w*t - phase(chan)));
                finalWaveformRefIndex(chan) = obj.AWG.AddWaveformByFunction(f, 3*minWaveformTime+obj.piTime+pulseTime,  0, 1);
            end
            
            % main pulse train:
            for chan = 1:obj.nChannels
                total_duration = obj.segmentDurs(1);
                sumMainFunction = @(t) 0*t;
                % Calculate the total duration of the waveform
                for i = 1:length(obj.pulseAngles)
                    % if ~obj.EchoLike
                    %     powerOn = 1;
                    % else
                    %     powerOn = 0;
                    % end
                    [durs, signs] = obj.pulseToDurations(obj.directions(i,:), obj.pulseAngles(i));
                    if sum(durs)
                        f = @(t) A(chan) * ((total_duration                <=t & t<total_duration+durs(1))                 * signs(1) .* cos(w*t - phase(chan)) + ...
                                            (total_duration+durs(1)        <=t & t<total_duration+durs(1)+durs(2))         * signs(2) .* cos(w*t + pi/2 - phase(chan)) + ...
                                            (total_duration+durs(1)+durs(2)<=t & t<total_duration+durs(1)+durs(2)+durs(3)) * signs(3) .* cos(w*t - phase(chan)));
                        sumMainFunction = @(t) sumMainFunction(t) + f(t);
                    end
                    total_duration = total_duration + obj.segmentDurs(i+1)+sum(durs);
                end
                %if total_duration~=obj.pulseTrainTime
                %    error('problem in calculation')
                %end
                if ~obj.EchoLike
                    mainWaveformsIndex(chan) = obj.AWG.AddWaveformByFunction(sumMainFunction, total_duration, 0, 0);
                else
                    [durs, signs] = obj.pulseToDurations([1, 0, 0], pi);
                    f = @(t) A(chan) * (((total_duration-sum(durs))/2                <=t & t<(total_duration-sum(durs))/2+durs(1))                 * signs(1) .* cos(w*t - phase(chan)) + ...
                                        ((total_duration-sum(durs))/2+durs(1)        <=t & t<(total_duration-sum(durs))/2+durs(1)+durs(2))         * signs(2) .* cos(w*t + pi/2 - phase(chan)) + ...
                                        ((total_duration-sum(durs))/2+durs(1)+durs(2)<=t & t<(total_duration-sum(durs))/2+durs(1)+durs(2)+durs(3)) * signs(3) .* cos(w*t - phase(chan)));
                    mainWaveformsIndex(chan) = obj.AWG.AddWaveformByFunction(f, total_duration, 0, 0);
                end
             end
             
        
            obj.AWG.LoadWaveform();
                
            % creating sequences:
            obj.indexSeq = zeros(obj.nChannels, N);
            fprintf('Loading %d sequences: ', obj.nChannels*N);
            for chan = 1:obj.nChannels
                for i = 1:N
                    obj.indexSeq(chan,i) = obj.AWG.AddSequence([initWaveformIndex(chan), mainWaveformsIndex(chan), finalWaveformIndex(chan), ...
                                                                initWaveformIndex(chan), mainWaveformsIndex(chan), finalWaveformRefIndex(chan)], ...
                                                                [1 obj.nRepetitions(i) 1 1 obj.nRepetitions(i) 1], [1 0 0 1 0 0]);%, [w, A(chan)/2, phaseZero(chan)/w]);
                    fprintf('%d ', N*(chan-1)+i);
                end
            end
            fprintf('\n');
            obj.AWG.LoadSequence()
        end
        
        
        function perform(obj)
            %%% Initialization
            
            % Devices (+ Tracker)
            pg = getObjByName(PulseGenerator.NAME);
            spcm = getObjByName(Spcm.NAME);
            tracker = getObjByName(Tracker.NAME);
                if isempty(tracker); throwBaseObjException(Tracker.Name); end
            
            % Some magic numbers
            [durs, ~] = obj.pulseToDurations(obj.initPulseDirection, obj.initPulseAngle);
            maxLastDelay = obj.lastDelay + max(obj.nRepetitions) * obj.pulseTrainTime + 2*sum(durs);

            %%% Run - Go over all tau's, in random order
            for t = randperm(length(obj.nRepetitions))
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
                        expDuration = obj.nRepetitions(t) * obj.pulseTrainTime + 2*sum(durs);
                        pg.changeSequence('tau', 'duration', expDuration);
                        if obj.constantTime
                            pg.changeSequence('lastDelay', 'duration', maxLastDelay - expDuration);
                        end
                        
                        obj.AWG.assignSequence(obj.indexSeq(1,t),1)
                        if obj.nChannels == 2
                            obj.AWG.assignSequence(obj.indexSeq(2,t),2);
                        end
                        obj.AWG.Run;
                        
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
            obj.signalParam.value = value;
            obj.signalParam.sterr = sterr;
            
            if obj.doubleMeasurement
                S3 = squeeze(obj.signal(3, :, 1:obj.currIter));
                S3sterr = squeeze(obj.sterr(3, :, 1:obj.currIter));
                S4 = squeeze(obj.signal(4, :, 1:obj.currIter));
                S4sterr = squeeze(obj.sterr(4, :, 1:obj.currIter));
                
                [value, sterr] = getRatioDistributionValues(obj, S3, S4, S3sterr, S4sterr);
                obj.signalParam2.value = value;
                obj.signalParam2.sterr = sterr;
            else    
                obj.signalParam2.value = [];
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
            
            N1 = obj.signalParam.value;
            N0 = obj.signalParam2.value;

            value = N0 - N1;
            sterr = sqrt(obj.signalParam.sterr.^2 + obj.signalParam2.sterr.^2);
            dataParam = ExpResultDoubleVector('FL', value, sterr, 'Normalized', obj.NAME);
        end
    end
end

