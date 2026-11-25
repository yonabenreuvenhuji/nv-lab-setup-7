classdef AWGPulsedDynamicalDecoupling < ExperimentAWG
    %EXPECHO Echo experiment

    properties (Constant)
        NAME = 'AWG Pulsed Dynamical Decoupling';
    end
    
    properties
        % All of the parameters must have a set function, and when they are
        % changed, the parameter obj.changeFlag needs to be set to true
        % (unless the change is minor and the experiment can continue
        % running);
        tau                 % in us 
        pulsePhase          % Phase of the pi pulses. Can be vector (like 'XYYX' for XY4) and the phases will repeats the order.
        cycles              % Number of pi pulses in each sequence. 1 for regular echo.
        
        doubleMeasurement   % logical
        
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = AWGPulsedDynamicalDecoupling()
            obj@ExperimentAWG(AWGPulsedDynamicalDecoupling.NAME);
            obj.parameterName = 'pulses';
            
            % Set properties inherited from Experiment
            obj.repeats = 5000;
            obj.averages = 10;
            
            obj.frequency = 2766;           % in MHz
            obj.amplitude = -15;            % in dBm
            obj.tau = 80;                   
            obj.halfPiTime = 0.0180;        % in us
            obj.piTime = 0.0367;            % in us
            obj.lastDelay = Experiment.DEFAULT_LAST_DELAY;
            
            obj.pulsePhase = 'XYYX';        % Phase of the pi pulses. Can be vector (like 'XYYX' for XY4) and the phases will repeats the order.
            obj.cycles = 4:4:80;            % Number of pi pulses in each sequence. 1 for regular echo.
            
            obj.constantTime = false;       % logical
            obj.doubleMeasurement = true;   % logical
            
            obj.detectionDuration = 0.5;           % detection window, in us
            obj.referenceDetectionDuration = 0.5;     % in us. Detection duration of the reference read
            obj.laserInitializationDuration = 10;   % laser initialization in pulsed experiments in us
            
            obj.mCurrentXAxisParam = ExpParamDoubleVector('Cycles', [], [], '', obj.NAME);
            obj.topParam = ExpParamDoubleVector(StringHelper.TAU, [], [], StringHelper.MICROSEC, obj.NAME);
            obj.signalParam = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME, '|1>');
            obj.signalParam2 = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME, '|0>');
        end
    end
    
    %% Setters
    methods
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

        function set.doubleMeasurement(obj, newVal)
            checkBoolean(obj, newVal)
            % If we got here, then newVal is OK.
            obj.doubleMeasurement = newVal;
            obj.changeFlag = true;
        end
    end
    
    methods
        function totalParamNum = getTotalNumberOfParams(obj)
            totalParamNum = length(obj.cycles);
        end
    end
    
    %% Overridden from Experiment
    methods (Access = protected)
        function prepare(obj)
            % Initialize devices (SPCM, PulseGenerator, etc.)
            
            % Sequence
            %%% Useful parameters for what follows
            initDuration = obj.laserInitializationDuration-obj.detectionDuration;
            obj.detectionPeriodsPerRepeat = 2 * (1 + double(obj.doubleMeasurement));
            obj.runsPerPerform = 1;
            MWChannel = obj.MWChannel;

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
            
            S = Sequence;
            for k = 1:1+obj.doubleMeasurement
                S.addEvent(obj.smallDelay,                  '');
                S.addEvent(obj.AWG.TRIG_DURATION,           [triggers(:)', MW(:)']);
                S.addEvent(2*obj.tau*obj.cycles(end) + 2*obj.halfPiTime, MW,                    '2tau');
                S.addEvent(obj.smallDelay,                  '');
                S.addEvent(obj.lastDelay,                   '',             'lastDelay');       % Last delay
                
                S.addEvent(obj.detectionDuration,       {'greenLaser', 'detector'});                    % Detection
                S.addEvent(initDuration,                'greenLaser');                                  % Initialization
                if ~obj.balancedSequence
                    S.addEvent(obj.referenceDetectionDuration,...
                        {'greenLaser', 'detector'});                                                        % Reference detection
                    S.addEvent(1,                           'greenLaser');                                  % Initialization
                else
                    S.addEvent(obj.smallDelay,                  '');
                    S.addEvent(obj.AWG.TRIG_DURATION,           '');                        % AWG trigger
                S.addEvent(2*obj.tau*obj.cycles(end) + 2*obj.halfPiTime, '',                    '2tau');
                    S.addEvent(obj.smallDelay,                  '');
                    S.addEvent(obj.lastDelay,                   '',                         'lastDelay');       % Last delay
                    S.addEvent(obj.referenceDetectionDuration,...
                        {'greenLaser', 'detector'});                                                        % Reference detection
                    S.addEvent(initDuration,                'greenLaser');                                  % Initialization
                end
            end
            
            if obj.changeFlag
                obj.LoadAWG;
            end

            obj.prepareInternal(S)
            
            % Set parameter, for saving
            obj.mCurrentXAxisParam.value = obj.cycles;
            obj.topParam.value = 2*obj.tau*obj.cycles;
            
            obj.isPlotAlternateAvailable = obj.doubleMeasurement;
        end
    
        
        function LoadAWG(obj)
            baseWFtime = obj.tau/round(obj.tau/0.5);
            baseWFtime = obj.RoundDurationByFrequencies(baseWFtime, obj.frequency, 1);
            obj.AWG.ResetStoredData;
            obj.AWG.workChannels = obj.activeChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            obj.AWG.ChangeSampleRate(baseWFtime);
            minWaveformTime = obj.AWG.minWaveformDuration;
            
            obj.halfPiTime = obj.RoundDurationByFrequencies(obj.halfPiTime , obj.AWG.SampleRate);
            obj.tau = ceil(obj.tau/baseWFtime)*baseWFtime;
            fprintf('The ''tau'' value rounded to be a complete period of the baseWFtime - %d.', baseWFtime);
            
            N = length(obj.cycles);
            A = obj.AWG.MaxNormPowerAllowed;
            w = 2*pi*obj.frequency;
            
            % creating waveforms:
            fprintf('Loading %d waveforms: ', 1 + obj.nChannels * 5);
            firstPiHalfIndex = zeros(1, obj.nChannels);
            lastPiHalfIndex = zeros(1, obj.nChannels);
            lastPiHalfRefIndex = zeros(1, obj.nChannels);
            centerPulseXIndex = zeros(1, obj.nChannels);
            centerPulseYIndex = zeros(1, obj.nChannels);
            
            [prePulseTime, phase] = obj.channelsDifference;
            
            f = @(t) 0 * t;
            zeroStepIndex = obj.AWG.AddWaveformByFunction(f, baseWFtime, 0);
            for chan = 1:obj.nChannels
                phaseZero = w * (prePulseTime(chan)+obj.halfPiTime) + phase(chan);
                f = @(t) A(chan) * (prePulseTime(chan)<=t & t<prePulseTime(chan)+obj.halfPiTime) .* cos(w*t - phaseZero);
                firstPiHalfIndex(chan)   = obj.AWG.AddWaveformByFunction(f, prePulseTime(chan)+obj.halfPiTime, 0, 0);
                f = @(t) A(chan) * (0<=t & t<obj.halfPiTime) .* cos(w*t - phase(chan));
                lastPiHalfIndex(chan)    = obj.AWG.AddWaveformByFunction(f, minWaveformTime+obj.halfPiTime,  0, 1);
                f = @(t) A(chan) * (0<=t & t<obj.halfPiTime) .* cos(w*t  - phase(chan)+ pi);
                lastPiHalfRefIndex(chan) = obj.AWG.AddWaveformByFunction(f, minWaveformTime+obj.halfPiTime,  0, 1);
                centerPulseTime = 2*baseWFtime;
                f = @(t) A(chan) * ((centerPulseTime/2-obj.piTime/2)<=t & t<(centerPulseTime/2+obj.piTime/2)) .* cos(w*t - phase(chan));
                centerPulseXIndex(chan)  = obj.AWG.AddWaveformByFunction(f, centerPulseTime, 0, 0);
                f = @(t) A(chan) * ((centerPulseTime/2-obj.piTime/2)<=t & t<(centerPulseTime/2+obj.piTime/2)) .* cos(w*t - phase(chan) + pi/2);
                centerPulseYIndex(chan)  = obj.AWG.AddWaveformByFunction(f, centerPulseTime, 0, 0);
            end
            fprintf('\n');
            obj.AWG.LoadWaveform();
            
            % creating sequences:
            obj.indexSeq = zeros(obj.nChannels, N);
            X_pulses = obj.pulsePhase=='X';
            X = repmat(X_pulses, 1, ceil(obj.cycles(end)/length(obj.pulsePhase)));
            fprintf('Loading %d sequences: ', obj.nChannels*N*2);
            for chan = 1:obj.nChannels
                for i = 1:N
                    totalZeroTime  = obj.tau-centerPulseTime/2;
                    zeroSegNum = totalZeroTime/baseWFtime;
                    piPulses = centerPulseXIndex(chan)*X + centerPulseYIndex(chan)*~X;
                    piPulses = piPulses(1:obj.cycles(i));
                    seqBody = reshape([repmat(zeroStepIndex,1,obj.cycles(i)); piPulses; repmat(zeroStepIndex,1,obj.cycles(i))], 1, 3*obj.cycles(i));
                    waveFormsIndexes = [firstPiHalfIndex(chan) seqBody lastPiHalfIndex(chan) ...
                                        firstPiHalfIndex(chan) seqBody lastPiHalfRefIndex(chan)];
                    repeats =  [1 repmat([zeroSegNum 1 zeroSegNum],1,obj.cycles(i)) 1, ...
                                1 repmat([zeroSegNum 1 zeroSegNum],1,obj.cycles(i)) 1];
                    triggers = [1 repmat([0          0 0         ],1,obj.cycles(i)) 0, ...
                                1 repmat([0          0 0         ],1,obj.cycles(i)) 0];
                    
                    %obj.indexSeq(chan, i) = obj.AWG.AddSequence(waveFormsIndexes, repeats, triggers);
                    obj.indexSeq(chan, i) = obj.AWG.AddSequence(waveFormsIndexes, fix(repeats), triggers);
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
            maxLastDelay = obj.lastDelay + 2 * obj.tau*max(obj.cycles)+ 2*obj.halfPiTime;
            
            %%% Run - Go over all tau's, in random order
            for t = randperm(length(obj.cycles))
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
                        pg.changeSequence('2tau', 'duration', 2*obj.tau*obj.cycles(t) + 2*obj.halfPiTime);
                        if obj.constantTime
                            pg.changeSequence('lastDelay', 'duration', maxLastDelay - 2*obj.tau*obj.cycles(t) - 2*obj.halfPiTime);
                        end
                        
                        for i = 1:length(obj.activeChannels)
                            obj.AWG.assignSequence(obj.indexSeq(i,t), obj.activeChannels(i));
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

