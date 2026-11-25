classdef AWGExpWeakMeasCorr < ExperimentAWG
    
    properties (Constant)
        NAME = 'AWG Weak Measurement Correlation';
    end
    
    properties
        % All of the parameters must have a set function, and when they are changed, the parameter obj.changeFlag needs to be set to true
        % (unless the change is minor and the experiment can continue running);
        tau             % in us
        totalZero
        smallOnDelay    % in us
        smallOffDelay   % in us
        dontLoadAWG
                projection      % projection orientation - 1 for z or  0 for x (Sz or Sz operator)
        theta           % Initial state of experiment (angle of Bloch vector [rad] around y axis)
        tauProbe
        param_idx
        weakValue
        weakSterr
        strongValue
        strongSterr
        correlation
          
    end

    methods
        function obj = AWGExpWeakMeasCorr()
            obj@ExperimentAWG(AWGExpWeakMeasCorr.NAME);
            obj.parameterName = 'taus';
            

            obj.smallOnDelay = 0;         % us
            obj.smallOffDelay = 0; %1.75;        % us
            obj.lastDelay = Experiment.DEFAULT_LAST_DELAY;
            obj.dontLoadAWG = 0;
            

            obj.repeats = 10000;
            obj.averages = 1000;
            obj.nChannels=1;
            obj.frequency = 2870; %in MHz
            obj.amplitude = -5; % in dBm
            obj.tauProbe = 0.02; % weak measurement probing duration
            obj.constantTime = true;
            obj.tau = 0.01:0.02:0.3;  %0.005:0.001:0.05; % in us (the running parameter)
            obj.piTime = 0.035;
            obj.theta = 0;
            obj.projection = 0;  % Sz default, for Sx insert 1
%             obj.weakValue = zeros(1, size(obj.tau));
            
            obj.detectionDuration = 5; % detection window, in us
            obj.referenceDetectionDuration = 5; % in us. Detection duration of the reference read
            obj.laserInitializationDuration = 15; % laser initialization in pulsed experiments

            obj.displayType1 = 'Weak Signal';
            obj.displayType2 = 'Both Signals';
            obj.isPlotAlternateAvailable = true;
            
            obj.mCurrentXAxisParam = ExpParamDoubleVector('Time', [], [], StringHelper.MICROSEC, obj.NAME);
            obj.signalParam = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME);
% Correlation array 05.05.2025
            obj.correlation = zeros(length(obj.tau), obj.averages);
        end
    end

%% Setters
    methods
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
               S = Sequence;
           
            S.addEvent(obj.AWG.TRIG_DURATION,           [triggers(:)', MW(:)']);                        % AWG trigger
            % state MW initialization (theta pulse)
            S.addEvent(obj.piTime * obj.theta / pi+1,        MW,                         'initPulse')               % MW state initialization
               
            % weak detection
%             S.addEvent(0.06, {'detector'}) % adding detection before laser
            S.addEvent(obj.tauProbe+1,        {'greenLaser','detector'}, 'weakSignal');  % add greelaser
%             S.addEvent(0.2, {'detector'}) % adding detection after laser


            % System Hamiltonian evolution (correlation time)
%             S.addEvent(0.06, {'detector'}) % adding detection before laser
            S.addEvent(max(obj.tau)+1,        MW, 'MW');  % add greelaser
%             S.addEvent(0.2, {'detector'}) % adding detection after laser

            % small delay
            S.addEvent(Experiment.DEFAULT_LAST_DELAY,           '',                         'lastDelay');       % Last delay

            % strong detection
%             S.addEvent(0.02, {'detector'}) % adding detection before laser
            S.addEvent(obj.detectionDuration, {'greenLaser', 'detector'}, 'StrongMeasurements');                    % Detection
%             S.addEvent(0.045, {'detector'}) % adding detection after laser

            % small delay
%             S.addEvent((3/2)*obj.piTime,                             '',                         '');

            % reference preperation and detection
            S.addEvent(obj.laserInitializationDuration,        'greenLaser');                                  % Initialization
            S.addEvent(obj.referenceDetectionDuration,...
                                            {'greenLaser', 'detector'});                    % Reference detection
            S.addEvent(1,                   'greenLaser');                                  % 1usec delay to ensure full overlap of ref detection and laser, and to init first repaet (less important), Yachel 27.02.22

            obj.prepareInternal(S)
            
            if (obj.changeFlag || obj.restartFlag) && (~obj.dontLoadAWG)
                obj.LoadAWG;
            end

            % Set parameter, for saving
            obj.mCurrentXAxisParam.value = obj.tau;
            pairName = 'Weak Signal';
            obj.signalParam.desc = pairName;
            
            obj.isPlotAlternateAvailable = true;
        end
        
        function LoadAWG(obj)
            baseWFtime = 0.2;
            baseWFtime = obj.RoundDurationByFrequencies(baseWFtime, obj.frequency);
            obj.AWG.ResetStoredData;
            obj.AWG.workChannels = obj.activeChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            obj.AWG.ChangeSampleRate(baseWFtime);
            minWaveformTime = obj.AWG.minWaveformDuration;
            
            N = length(obj.tau);
            shortWaveformsIndex = zeros(obj.nChannels,0);
            shortWaveformsTime = [];
            waveformsIndex = zeros(obj.nChannels,N);
            numZeroStepsInPoint = zeros(1,N);
            initZeroIndex = zeros(1,obj.nChannels);
            constStepIndex = zeros(1,obj.nChannels);
            
            % creating waveforms:
            A = obj.AWG.MaxNormPowerAllowed;
            w = 2*pi*obj.frequency;
            fprintf('Loading ~%d waveforms: ', (sum((obj.tau<baseWFtime)) + 2)*obj.nChannels);
            
            % the resolution of the time delay is the sample rate, so we add also aditional phase delay
            [prePulseTime, phase] = obj.channelsDifference;
            for chan = 1:obj.nChannels
                f = @(t) 0 * t;
                initZeroIndex(chan) = obj.AWG.AddWaveformByFunction(f, prePulseTime(chan), 0);
                f = @(t) A(chan) * cos(w*t - phase(chan));
                constStepIndex(chan) = obj.AWG.AddWaveformByFunction(f, baseWFtime, 0);
            end
            
            for i = 1:N
                numZeroStepsInPoint(i) = floor(obj.tau(i) / baseWFtime);
                residue = mod(obj.tau(i), baseWFtime);
                [closeWFdiff, closeWFindex] = min(abs((residue - shortWaveformsTime)));
                ind = length(shortWaveformsTime) + 1;
                for chan = 1:obj.nChannels
                    if closeWFdiff < 0.01*obj.tau(i)
                        waveformsIndex(chan, i) = shortWaveformsIndex(chan, closeWFindex);
                        obj.tau(i) = numZeroStepsInPoint(i) * baseWFtime + shortWaveformsTime(closeWFindex);
                    else
                        shortWaveformsTime(ind) = residue;
                        f = @(t) A(chan) * (0<=t & t<residue) .* cos(w*t - phase(chan));
                        shortWaveformsIndex(chan,ind) = obj.AWG.AddWaveformByFunction(f, residue + minWaveformTime, 0, 1);
                        waveformsIndex(chan, i) = shortWaveformsIndex(chan, ind);
                    end
                end
            end
            fprintf('\n');
            obj.AWG.LoadWaveform();
           % delete double points in tau (because of the round):
            while any(diff(obj.tau) == 0)
                i = find(diff(obj.tau) == 0, 1);
                obj.tau(i) = [];
                waveformsIndex(:, i) = [];
                numZeroStepsInPoint(i) = [];
            end
            
            % creating sequences:
            N = length(obj.tau);
            obj.indexSeq = zeros(obj.nChannels,N);
            for i = 1:N
                for chan = 1:obj.nChannels
                    
                    draw = 0; %(chan-1)*mod(i,2);
%                     draw = chan==2;
%                     if i == N && chan==2
%                         draw = [w, A(chan)*0.8, minWaveformTime-obj.timeDelay*(chan==1)-5e-5+4e-4];
% %                         figure;
%                     end
                    

                    if numZeroStepsInPoint(i) == 0
                        obj.indexSeq(chan,i) = obj.AWG.AddSequence([initZeroIndex(chan), waveformsIndex(chan,i)], ...
                                                                   [1                  , 1], ...
                                                                   [1                  , 0], draw);
                    else
                        obj.indexSeq(chan,i) = obj.AWG.AddSequence([initZeroIndex(chan), constStepIndex(chan),    waveformsIndex(chan,i)], ...
                                                                   [1                  , numZeroStepsInPoint(i), 1], ...
                                                                   [1                  , 0                       0], draw);
                    end
                end    
            end
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
            maxLastDelay = obj.lastDelay + max(obj.tau);
            
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
                        pg.changeSequence('MW', 'duration', obj.tau(t));
                        if obj.constantTime
                            pg.changeSequence('lastDelay', 'duration', maxLastDelay - obj.tau(t));
                        end
                        for i = 1:length(obj.activeChannels)
                            obj.AWG.assignSequence(obj.indexSeq(i,t), obj.activeChannels(i));
                        end
                        obj.AWG.Run;
                        pause(1);
                        
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
            
             % autocorrelation calculation
            rawDataReshaped = (reshape(data, 3, round(length(data)/3)))';
            for k=1:length(obj.weakDetectionDuration)
                cur_index = f1(k);
                cur_data = rawDataReshaped(obj.repeats*(k-1)+1:obj.repeats*k,:);
                cur_multipication = cur_data(:,2).*cur_data(:,1);
                cur_corr = mean(cur_multipication);
                obj.correlation(cur_index,obj.currIter) = cur_corr;

                
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
