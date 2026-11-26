classdef AWGExpEcho < ExperimentAWG
    %EXPECHO Echo experiment

    properties (Constant)
        NAME = 'AWG Echo';
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

        useSequenceWalk = true;     % logical. (default) Arm ONCE and let the sequencer advance one tau
                                    % per pulse-generator repeat, instead of MATLAB assigning
                                    % a new sequence at every point.
                                    %
                                    % ECHO IS THE EASIEST OF THE THREE TO WALK, and for a
                                    % reason worth knowing: tau already lives in the AWG
                                    % SEQUENCE, as the repeat count of a zero-amplitude step
                                    % (see LoadAWG). The waveform set is fixed and tiny -- one
                                    % zero step plus five per channel -- and it does not grow
                                    % with the tau list at all. So walking Echo costs no extra
                                    % upload whatsoever; it only concatenates sequences that
                                    % were being built anyway.
                                    %
                                    % Each tau's existing sequence is ALREADY a walk block:
                                    % 3*cycles+2 steps for the signal, the same again for the
                                    % reference, and exactly two trigger-gated steps among
                                    % them -- one per pulse-generator repeat. Concatenating
                                    % them changes who advances between taus (the sequencer
                                    % rather than MATLAB) and nothing about what the NV sees.
                                    %
                                    % WHAT IT DOES CHANGE, and this is the one thing to test
                                    % rather than assume: the microwave gate can no longer be
                                    % rewritten per point, so it stays at the LONGEST tau. At
                                    % tau = 1 us in a list reaching 50 us that leaves ~98 us
                                    % of gate open with the AWG idling at zero. The Ramsey
                                    % gate test cleared 4 us of excess with no measurable
                                    % contrast loss; 98 us is twenty-five times beyond
                                    % anything that has been measured on this setup. Run
                                    % echoGateTest before trusting a walked Echo, and read the
                                    % warning prepare prints.
        passesPerRun = [];          % full passes over the tau list per pg.run. [] = derived
                                    % from maxDigitizerSegments and maxRunSeconds.
        passesPerRunUsed = [];      % what the auto choice settled on.
        walkWarmupTriggers = 1;     % throwaway blocks at the START of each pass, so the first
                                    % real block is not the one that follows the idle gap.
        resyncEveryPass = true;     % stop and re-arm at the start of every pass, so a pass
                                    % that aborts cannot shift the tau mapping of the rest.
        maxDigitizerSegments = 32000;
                                    % ceiling on nDet * blocksPerPass * passesPerRun for ONE
                                    % pg.run.
        maxRunSeconds = 30;         % ceiling on how long ONE pg.run may take, which Ramsey
                                    % does not need and Echo very much does.
                                    %
                                    % Echo's sequence duration is 2*max(tau)*cycles, and with
                                    % the gate fixed that is what EVERY block costs. At
                                    % tau -> 200 us with cycles = 8 one block is 3.2 ms, so a
                                    % 300-pass run over 25 taus would hold the instruments for
                                    % four minutes in a single uninterruptible call -- no
                                    % plot updates, no emergency stop, and a digitizer timeout
                                    % if the estimate is off. Capping the run length trades a
                                    % few extra re-arms for a responsive experiment.
        walkOrder = [];             % the order the blocks were built in, so the data can be
                                    % put back in tau order. Set by buildWalkSequence.
        newOrderEveryAverage = true;    % logical. WALK MODE ONLY. Draw a new random tau
                                    % order before every average and re-upload ONLY the
                                    % sequence -- the waveforms stay on the AWG. Costs one
                                    % sequence upload + arm per average, printed as it happens.
                                    % Without it the order is drawn once per waveform load and
                                    % every pass and every average reuse it, so anything tied
                                    % to position in the pass (the blocks right after each
                                    % re-arm, the previous block's tau) sits on the SAME taus
                                    % for the whole run and cannot average down. With it, that
                                    % becomes noise that does. Same option as AWGExpESR. The
                                    % classic path already draws a new order every average.

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
    end

    
    properties (SetAccess = protected)
        walkParts = []          % struct. What uploadWalkSequence needs to rebuild the walk
                                % sequence in a new order without touching the waveforms.
        walkOrderFresh = false  % true while the loaded order has not been used by an average
                                % yet, so newOrderEveryAverage does not redraw it.
        seqAssetIndex               % the ONE sequence asset per channel, in walk mode. Kept
                                    % apart from indexSeq because in walk mode indexSeq holds
                                    % STEP numbers inside that single asset, whereas on the
                                    % classic path it holds one asset index per tau. Two
                                    % different things; two different properties.
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = AWGExpEcho()
            obj@ExperimentAWG(AWGExpEcho.NAME);
            obj.parameterName = 'taus';
            
            % Set properties inherited from Experiment
            obj.repeats = 5000;
            obj.averages = 10;
            
            obj.frequency = 2766;           % in MHz
            obj.amplitude = -15;            % in dBm
            obj.tau = 0.2:1:300;                   
            obj.halfPiTime = 0.0180;         % in us
            obj.piTime = 0.0367;              % in us
            obj.lastDelay = Experiment.DEFAULT_LAST_DELAY;
            
            obj.pulsePhase = 'Y';           % Phase of the pi pulses. Can be vector (like 'XYYX' for XY4) and the phases will repeats the order.
            obj.cycles = 1;                 % Number of pi pulses in each sequence. 1 for regular echo.
            
            obj.constantTime = false;       % logical
            obj.doubleMeasurement = true;   % logical
            
            obj.detectionDuration = 0.5;           % detection window, in us
            obj.referenceDetectionDuration = 0.5;     % in us. Detection duration of the reference read
            obj.laserInitializationDuration = 10;   % laser initialization in pulsed experiments in us
            
            obj.mCurrentXAxisParam = ExpParamDoubleVector('Time', [], [], StringHelper.MICROSEC, obj.NAME);
            % Labels as in the robust fit's figure (ExpPlotStyle / fitStyle)
            obj.plotXLabel = ExpPlotStyle.X_ECHO;     % the x axis is 2*tau*cycles
            obj.plotYLabel = ExpPlotStyle.Y_NORMALISED;
            obj.plotYLabelAlternate = ExpPlotStyle.Y_CONTRAST;   % alternateSignal = |0> - |1>
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
        function set.useSequenceWalk(obj, newVal)
            if ~ValidationHelper.isTrueOrFalse(newVal)
                obj.sendError('useSequenceWalk must be true or false')
            end
            obj.useSequenceWalk = logical(newVal);
            obj.changeFlag = true;
        end

        function set.passesPerRun(obj, newVal)
            if ~isempty(newVal) && ~ValidationHelper.isValuePositiveInteger(newVal)
                obj.sendError('passesPerRun must be a positive integer, or [] to derive it')
            end
            obj.passesPerRun = newVal;
            obj.changeFlag = true;
        end

        function set.walkWarmupTriggers(obj, newVal)
            if ~ValidationHelper.isValueInteger(newVal) || newVal < 0
                obj.sendError('walkWarmupTriggers must be a non-negative integer')
            end
            obj.walkWarmupTriggers = newVal;
            obj.changeFlag = true;
        end

        function set.resyncEveryPass(obj, newVal)
            if ~ValidationHelper.isTrueOrFalse(newVal)
                obj.sendError('resyncEveryPass must be true or false')
            end
            obj.resyncEveryPass = logical(newVal);
        end

        function set.newOrderEveryAverage(obj, newVal)
            if ~ValidationHelper.isTrueOrFalse(newVal)
                obj.sendError('newOrderEveryAverage must be true or false')
            end
            obj.newOrderEveryAverage = logical(newVal);
            % No changeFlag: the first upload is the same either way; the redraws happen
            % inside performWalk, one per average.
        end

        function set.maxDigitizerSegments(obj, newVal)
            if ~ValidationHelper.isValuePositiveInteger(newVal)
                obj.sendError('maxDigitizerSegments must be a positive integer')
            end
            obj.maxDigitizerSegments = newVal;
            obj.changeFlag = true;
        end

        function set.maxRunSeconds(obj, newVal)
            if ~isnumeric(newVal) || ~isscalar(newVal) || ~(newVal > 0)
                obj.sendError('maxRunSeconds must be a positive number')
            end
            obj.maxRunSeconds = newVal;
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
                S.addEvent(2*max(obj.tau)*obj.cycles + 2*obj.halfPiTime, MW,               '2tau');   % MAX, not tau(end)
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
                    S.addEvent(2*max(obj.tau)*obj.cycles + 2*obj.halfPiTime, '',               '2tau');   % MAX, not tau(end)
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
            % max(obj.tau), NOT obj.tau(end), in the two gate events above. The classic path
            % rewrites the gate at every point so it never noticed, but nothing sorts the tau
            % list -- LoadAWG rounds it and removes duplicates and leaves the order alone -- so
            % with an unsorted list tau(end) can be SHORTER than the longest tau. The gate
            % would then close while the AWG was still playing, cutting the tail off the pulse
            % train, and no error would say so. In walk mode the gate is never rewritten, so
            % this is not a tidy-up: it is the difference between working and silently wrong.

            if obj.useSequenceWalk
                % THE GATE IS NOW FIXED AT THE LONGEST TAU, and unlike Ramsey the excess is
                % not small. tau is realised inside the AWG sequence as repeats of a zero
                % step, so the NV sees the right Echo whatever the gate does -- what is being
                % trusted is that the AWG really idles at zero while the switch is open, for
                % as long as the difference between this tau and the longest one.
                pg = getObjByName(PulseGenerator.NAME);
                spcm = getObjByName(Spcm.NAME);
                N = length(obj.tau);
                blocksPerPass = N + obj.walkWarmupTriggers;
                obj.passesPerRunUsed = obj.choosePassesPerRun(true);
                repeatsPerRun = blocksPerPass * obj.passesPerRunUsed;
                pg.repeats = repeatsPerRun;
                nDet = obj.detectionPeriodsPerRepeat;
                numScans = nDet * repeatsPerRun;
                seqTime = pg.sequenceDuration() * 1e-6;
                spcm.prepareExperimentCount(numScans, 2 * seqTime * repeatsPerRun, 1);

                gateUs = 2*max(obj.tau)*obj.cycles + 2*obj.halfPiTime;
                shortUs = 2*min(obj.tau)*obj.cycles + 2*obj.halfPiTime;
                fprintf(['Walk mode: %d repeats per pg.run (%d taus + %d warm-up, x %d ' ...
                         'passes), %d digitizer scans. MW gate fixed at %.2f us.\n'], ...
                        repeatsPerRun, N, obj.walkWarmupTriggers, obj.passesPerRunUsed, ...
                        numScans, gateUs);
                fprintf(['      Excess gate: 0 us at the longest tau, %.2f us at the shortest. ' ...
                         '\n'], gateUs - shortUs);
                if gateUs - shortUs > 10
                    fprintf(['      *** THAT IS UNTESTED TERRITORY. ramseyGateTest cleared 4.05 us\n' ...
                             '      of excess gate with no measurable contrast loss (106%% +/- 5%%).\n' ...
                             '      %.0f us is %.0fx beyond anything measured on this setup, and the\n' ...
                             '      assumption being made is that the AWG idles at zero while the\n' ...
                             '      switch is open for that long. Run echoGateTest before quoting a\n' ...
                             '      T2 from a walked Echo. If it fails, band the tau list: walk a\n' ...
                             '      few narrow ranges, each with its own gate.\n'], ...
                            gateUs - shortUs, (gateUs - shortUs)/4.05);
                end
            end

            obj.mCurrentXAxisParam.value = 2*obj.tau*obj.cycles;
            obj.topParam.value = obj.tau;
            
            obj.isPlotAlternateAvailable = obj.doubleMeasurement;
        end
    
        
        function P = choosePassesPerRun(obj, announce)
            % How many passes to put in one pg.run, DERIVED rather than guessed, from TWO
            % ceilings. Ramsey has only the first; Echo needs both.
            %
            %   the digitizer:  nDet * blocksPerPass * P  <=  maxDigitizerSegments
            %   the clock:      blocksPerPass * P * sequenceDuration  <=  maxRunSeconds
            %
            % The second is what Echo adds. Its sequence duration is 2*max(tau)*cycles and the
            % fixed gate means every block pays it, so the run length grows with the tau range
            % in a way Ramsey's never does. Without the cap, a long-tau XY8 run disappears into
            % one uninterruptible call for minutes: no plot updates, no emergency stop, and a
            % digitizer timeout if the estimate is off.
            %
            % Walk pays only when it re-arms LESS than classic, i.e. when
            %       ceil(repeats/P)  <  nTau
            % because classic costs one arm per tau and walk costs one per pg.run. On a short
            % tau list with a long sequence both ceilings bite and walk can end up slower --
            % which is worth being told rather than discovering from the clock.
            N = length(obj.tau);
            nDet = obj.detectionPeriodsPerRepeat;
            if isempty(nDet) || nDet < 1; nDet = 2 * (1 + double(obj.doubleMeasurement)); end
            blocksPerPass = N + obj.walkWarmupTriggers;

            seqUs = NaN;
            try
                pg = getObjByName(PulseGenerator.NAME);
                seqUs = pg.sequenceDuration();
            catch
            end
            if ~isfinite(seqUs) || seqUs <= 0
                % Fall back to the arithmetic rather than refusing: gate plus the four
                % detection-and-initialise pairs is the whole repeat.
                seqUs = 2*max(obj.tau)*obj.cycles + 2*obj.halfPiTime ...
                        + nDet*(obj.laserInitializationDuration + obj.detectionDuration);
            end

            pDig = floor(obj.maxDigitizerSegments / (nDet * blocksPerPass));
            pClk = floor(obj.maxRunSeconds / max(blocksPerPass * seqUs * 1e-6, eps));
            if ~isempty(obj.passesPerRun)
                P = obj.passesPerRun;                       % explicit wins, but still reported
            else
                P = min(pDig, pClk);
            end
            P = max(1, min(P, obj.repeats));

            if announce
                fprintf(['Walk: passesPerRun = %d (%s). Ceilings: %d from the digitizer ' ...
                         '(%d scans of %d), %d from maxRunSeconds (%.0f s at %.0f us per ' ...
                         'block).\n'], ...
                        P, BooleanHelper.ifTrueElse(isempty(obj.passesPerRun), 'auto', 'set by hand'), ...
                        pDig, nDet*blocksPerPass*P, obj.maxDigitizerSegments, ...
                        pClk, obj.maxRunSeconds, seqUs);
                armsWalk = ceil(obj.repeats / max(P,1));
                fprintf('      re-arms per average: %d walk vs %d classic.\n', armsWalk, N);
                if armsWalk >= N
                    fprintf(['      *** WALK WILL NOT BE FASTER HERE: %d re-arms against\n' ...
                             '      classic''s %d. Walk pays only when ceil(repeats/P) < nTau.\n' ...
                             '      At this tau range one block costs %.0f us, so maxRunSeconds\n' ...
                             '      (%g s) caps P at %d. Lengthen the tau list, raise\n' ...
                             '      maxRunSeconds if you can accept a less responsive run, or\n' ...
                             '      set useSequenceWalk = false. Running anyway.\n'], ...
                            armsWalk, N, seqUs, obj.maxRunSeconds, pClk);
                end
            end
        end

        function LoadAWG(obj)
            baseWFtime = min([0.5, obj.tau(1), diff(obj.tau)]);      % the minimum length of waveform. All waveforms time must be a complete period of it!
            minWaveformTime = obj.AWG.MIN_WAVEFORM_LENGTH / obj.AWG.MAX_SAMPLE_RATE;
            if baseWFtime < minWaveformTime
                baseWFtime = minWaveformTime;
            end
            
            baseWFtime = obj.RoundDurationByFrequencies(baseWFtime, obj.frequency, 1);
            obj.AWG.ResetStoredData;
            obj.walkParts = [];     % set again below on the walk path
            obj.AWG.workChannels = obj.activeChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            obj.AWG.ChangeSampleRate(baseWFtime);
            minWaveformTime = obj.AWG.minWaveformDuration;
            
            obj.halfPiTime = obj.RoundDurationByFrequencies(obj.halfPiTime , obj.AWG.SampleRate);
            obj.tau = ceil(obj.tau/baseWFtime)*baseWFtime;
            fprintf('The ''tau'' values rounded to be a complete period of the baseWFtime - %d.', baseWFtime);
            % delete double points in tau:
            while any(diff(obj.tau) == 0)
                i = find(diff(obj.tau) == 0, 1);
                obj.tau(i) = [];
            end
            
            N = length(obj.tau);
            A = obj.AWG.MaxNormPowerAllowed;
            w = 2*pi*obj.frequency;
            
            % creating waveforms:
            % A FIXED, SMALL SET -- unlike Ramsey, Echo does not build one waveform per tau.
            % tau is realised by REPEATING the zero step inside the sequence, which is why the
            % sequence count grows with tau and the waveform count does not.
            obj.AWG.progressBegin(sprintf(['Echo waveforms: 1 zero step + 5 per channel ' ...
                '(first pi/2, last pi/2, last pi/2 reference, pi about X, pi about Y) ' ...
                'x %d channel(s). tau is built from repeats of the zero step, not from ' ...
                'separate waveforms.'], obj.nChannels), 1 + obj.nChannels*5, 'waveform');
            firstPiHalfIndex = zeros(1, obj.nChannels);
            lastPiHalfIndex = zeros(1, obj.nChannels);
            lastPiHalfRefIndex = zeros(1, obj.nChannels);
            centerPulseXIndex = zeros(1, obj.nChannels);
            centerPulseYIndex = zeros(1, obj.nChannels);
            
            [prePulseTime, phase] = obj.channelsDifference;
            % THE CALIBRATED RELATIVE PHASE, added to whatever channelsDifference produced.
            % Without this line channelPhase is a property nothing reads: channelsDifference
            % builds its phase from timeDelay alone, and level 4.3 sets timeDelay = 0 -- so
            % the two channels would be driven with ZERO relative phase while piTime and
            % halfPiTime in the json were measured with the calibrated phase applied. The
            % pulses would then not be pi and pi/2 at all, and nothing would say so.
            phase = phase + obj.channelPhaseRadians();
            
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
            obj.AWG.progressEnd();
            obj.AWG.LoadWaveform();
            
            % creating sequences:
            obj.indexSeq = zeros(obj.nChannels, N);
            X_pulses = obj.pulsePhase == 'X';
            X = repmat(X_pulses, 1, ceil(obj.cycles/length(obj.pulsePhase)));
            % THE COUNT WAS WRONG HERE TOO, the other way round: it said nChannels*N*2, but the
            % loop adds ONE sequence per (channel, tau) which happens to contain both the
            % signal and the reference block. nChannels*N is the number of sequences.
            obj.AWG.progressBegin(sprintf(['Echo sequences: one per tau per channel, ' ...
                '%d steps each (signal block + reference block at cycles = %d, phases %s).'], ...
                2*(3*obj.cycles + 2), obj.cycles, obj.pulsePhase), ...
                obj.nChannels*N, 'sequence');
            % The per-tau block is now built by blockFor, because BOTH paths need exactly the
            % same steps and the only difference is where they are put: the classic path gives
            % each block its own sequence asset, walk concatenates them all into one. Building
            % it twice would be two chances to let the paths drift apart, and a drift here
            % would not error -- it would measure a slightly different Echo on each path.
            idx = struct('zero', zeroStepIndex, 'firstHalf', firstPiHalfIndex, ...
                         'lastHalf', lastPiHalfIndex, 'lastHalfRef', lastPiHalfRefIndex, ...
                         'piX', centerPulseXIndex, 'piY', centerPulseYIndex);

            if obj.useSequenceWalk
                obj.AWG.progressEnd();
                % Kept apart from the sequence so newOrderEveryAverage can rebuild it in a
                % new order without re-uploading a single waveform. Field by field: struct()
                % would expand a cell value into a struct array.
                p = struct();
                p.idx             = idx;
                p.X               = X;
                p.baseWFtime      = baseWFtime;
                p.centerPulseTime = centerPulseTime;
                p.nWarm           = obj.walkWarmupTriggers;
                obj.walkParts = p;
                nSteps = obj.uploadWalkSequence(randperm(N));
                fprintf(['Walk load: %d taus + %d warm-up, %d steps per channel, %d AddSequence ' ...
                         'calls total (the classic path would make %d). The WAVEFORM count is ' ...
                         'unchanged -- Echo never had one per tau.\n'], ...
                        N, p.nWarm, nSteps, obj.nChannels, obj.nChannels*N);
                return
            end

            for chan = 1:obj.nChannels
                for i = 1:N
                    [waveFormsIndexes, repeats, triggers] = ...
                        obj.blockFor(chan, i, idx, X, baseWFtime, centerPulseTime);
                    obj.indexSeq(chan, i) = obj.AWG.AddSequence(waveFormsIndexes, fix(repeats), triggers);
                    obj.AWG.progressStep();
                end
            end
            obj.AWG.progressEnd();
            obj.AWG.LoadSequence()
        end

        function [wf, rpt, trg] = blockFor(obj, chan, i, idx, X, baseWFtime, centerPulseTime)
            % ONE tau's block: the signal half then the reference half, 3*cycles+2 steps each.
            %
            % tau is NOT a waveform here -- it is the repeat count of the zero step, which is
            % why Echo's waveform set is fixed at one zero plus five per channel however long
            % the tau list is, and why walking Echo costs no extra upload at all.
            %
            % Exactly TWO steps carry a trigger, the first of each half, so one block consumes
            % two AWG triggers -- which is precisely one pulse-generator repeat. That is what
            % makes this block a walk block without modification.
            totalZeroTime = obj.tau(i) - centerPulseTime/2;
            zeroSegNum = totalZeroTime/baseWFtime;
            piPulses = idx.piX(chan)*X + idx.piY(chan)*~X;
            piPulses = piPulses(1:obj.cycles);
            seqBody = reshape([repmat(idx.zero,1,obj.cycles); piPulses; ...
                               repmat(idx.zero,1,obj.cycles)], 1, 3*obj.cycles);
            wf  = [idx.firstHalf(chan) seqBody idx.lastHalf(chan) ...
                   idx.firstHalf(chan) seqBody idx.lastHalfRef(chan)];
            rpt = [1 repmat([zeroSegNum 1 zeroSegNum],1,obj.cycles) 1, ...
                   1 repmat([zeroSegNum 1 zeroSegNum],1,obj.cycles) 1];
            trg = [1 repmat([0          0 0         ],1,obj.cycles) 0, ...
                   1 repmat([0          0 0         ],1,obj.cycles) 0];
        end

        function nSteps = uploadWalkSequence(obj, order)
            % Build ONE sequence per channel that plays the taus in 'order' (order(k) = the
            % tau index of the k-th block), upload it, assign it and arm. Only the sequence
            % is replaced: the waveforms in obj.walkParts must already be on the AWG.
            % LoadAWG calls this for the first order and newWalkOrderIfAsked for every later
            % one, so both build identical blocks.
            %
            % ONE sequence per channel holding every tau, instead of N sequences.
            %
            % Each tau contributes the block blockFor already builds, unchanged. The last step
            % points back at step 1, so a pass ends where the next begins and pass two needs no
            % re-arm.
            %
            % THE BLOCK ORDER IS RANDOMISED ONCE, HERE. The classic path draws a fresh randperm
            % every average; walk cannot, because the visit order IS the sequence order. So the
            % randomisation moves into the sequence and the order is recorded, so performWalk
            % can put the data back in tau order. Without it, slow drift across a pass would
            % appear as a systematic trend in tau -- and on an Echo that is indistinguishable
            % from the T2 decay being measured.
            % (The order now comes in as 'order': drawn per load by LoadAWG and, with
            % newOrderEveryAverage, redrawn before every later average.)
            p = obj.walkParts;
            idx = p.idx; X = p.X; baseWFtime = p.baseWFtime; centerPulseTime = p.centerPulseTime;
            N = numel(order);
            nWarm = p.nWarm;
            blockOrder = [repmat(order(1), 1, nWarm), order];

            stepsPerBlock = 2*(3*obj.cycles + 2);
            nSteps = stepsPerBlock * numel(blockOrder);
            if nSteps > obj.AWG.MAX_SEQUENCE_STEPS
                obj.sendError(sprintf(['Walk needs %d sequence steps (%d per tau at cycles = ' ...
                    '%d, %d taus plus %d warm-up) and the AWG allows %d. Shorten the tau ' ...
                    'list, lower cycles, or use classic mode.'], nSteps, stepsPerBlock, ...
                    obj.cycles, N, nWarm, obj.AWG.MAX_SEQUENCE_STEPS));
            end

            % Nothing may still reference the sequence that is about to be replaced, and
            % LoadSequence deletes every sequence on the instrument before it uploads. The
            % list is cleared here so a re-upload REPLACES the walk sequence instead of
            % appending a second copy after it.
            obj.AWG.sendCommand('awgcontrol:stop:immediate');
            obj.AWG.sequences = struct('waveforms', {{}}, 'repeats', {{}}, ...
                                       'trigger', {{}}, 'goto', {{}});   % sequences only
            obj.indexSeq      = zeros(obj.nChannels, N);
            obj.seqAssetIndex = zeros(1, obj.nChannels);
            obj.AWG.progressBegin(sprintf(['Echo WALK sequence: ONE asset per channel, %d ' ...
                'steps (%d per tau x %d taus + %d warm-up block(s)) at cycles = %d, phases ' ...
                '%s. The sweep lives in the sequence, so the AWG is armed once per pass.'], ...
                nSteps, stepsPerBlock, N, nWarm, obj.cycles, obj.pulsePhase), ...
                obj.nChannels, 'sequence');
            for chan = 1:obj.nChannels
                wfList = zeros(1, nSteps); rptList = ones(1, nSteps);
                trgList = zeros(1, nSteps); gotoList = zeros(1, nSteps);
                s0 = 0;
                for blk = 1:numel(blockOrder)
                    i = blockOrder(blk);
                    [wf, rpt, trg] = obj.blockFor(chan, i, idx, X, baseWFtime, centerPulseTime);
                    wfList(s0+1:s0+stepsPerBlock)  = wf;
                    rptList(s0+1:s0+stepsPerBlock) = fix(rpt);
                    trgList(s0+1:s0+stepsPerBlock) = trg;
                    s0 = s0 + stepsPerBlock;
                    if blk > nWarm
                        obj.indexSeq(chan, i) = s0;     % last step of tau i's block
                    end
                end
                gotoList(nSteps) = 1;                   % wrap, so passes chain without a re-arm
                obj.seqAssetIndex(chan) = obj.AWG.AddSequence(wfList, rptList, trgList, 0, gotoList);
                obj.AWG.progressStep();
            end
            obj.AWG.progressEnd();
            obj.AWG.LoadSequence();

            for chan = 1:obj.nChannels
                obj.AWG.assignSequence(obj.seqAssetIndex(chan), obj.activeChannels(chan));
            end
            obj.AWG.Run;
            obj.walkOrder = order;
            obj.walkOrderFresh = true;      % not used by an average yet
        end
        
        
        function perform(obj)
            if obj.useSequenceWalk
                obj.performWalk;
                return
            end

            %%% Initialization
            
            % Devices (+ Tracker)
            pg = getObjByName(PulseGenerator.NAME);
            spcm = getObjByName(Spcm.NAME);
            tracker = getObjByName(Tracker.NAME);
                if isempty(tracker); throwBaseObjException(Tracker.Name); end
            
            % Some magic numbers
            maxLastDelay = obj.lastDelay + 2 * max(obj.tau)*obj.cycles+ 2*obj.halfPiTime;
            
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
                        pg.changeSequence('2tau', 'duration', 2*obj.tau(t)*obj.cycles + 2*obj.halfPiTime);
                        if obj.constantTime
                            pg.changeSequence('lastDelay', 'duration', maxLastDelay - 2*obj.tau(t)*obj.cycles - 2*obj.halfPiTime);
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
            
            obj.publishEcho();
        end

        function publishEcho(obj)
            % The two ratios, exactly as the classic path used to form them inline at the end
            % of perform. Factored out so walk and classic CANNOT drift apart: rows 1-2 are
            % the |1> projection and rows 3-4 the |0> one, and there is now one place that
            % decides that rather than two.
            S1 = squeeze(obj.signal(1, :, 1:obj.currIter));
            S1sterr = squeeze(obj.sterr(1, :, 1:obj.currIter));
            S2 = squeeze(obj.signal(2, :, 1:obj.currIter));
            S2sterr = squeeze(obj.sterr(2, :, 1:obj.currIter));

            [value, sterr] = getRatioDistributionValues(obj, S1, S2, S1sterr, S2sterr);
            obj.signalParam.value = value;
            obj.signalParam.sterr = sterr;

            if obj.doubleMeasurement && size(obj.signal, 1) >= 4
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

        function opts = walkReduceOptionsEcho(obj, spcm)
            % Everything walkReduce needs to reproduce processData exactly. Same contents as
            % the Rabi and Ramsey versions; duplicated rather than shared because these are
            % separate experiment classes and none should reach into another.
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

        function newWalkOrderIfAsked(obj)
            % newOrderEveryAverage: before every average except the first one after a fresh
            % load, draw a new tau order and re-upload ONLY the sequence.
            if obj.newOrderEveryAverage && ~obj.walkOrderFresh
                if isempty(obj.walkParts)
                    obj.sendWarning(['newOrderEveryAverage: the loaded sequence predates this ' ...
                        'option, so it cannot be rebuilt on its own. Keeping its order for this ' ...
                        'average; change any setting to force a reload.']);
                else
                    tNew = tic;
                    obj.uploadWalkSequence(randperm(length(obj.tau)));
                    fprintf('Average %d: new tau order -- sequence re-upload + arm %.1f s.\n', ...
                            obj.currIter, toc(tNew));
                end
            end
            obj.walkOrderFresh = false;
        end

        function performWalk(obj)
            % One average, acquired as whole PASSES over the tau list.
            %
            % The sequencer advances one BLOCK per pulse-generator repeat, so a pg.run of
            % (N + warmup) * passesPerRun repeats returns the whole sweep repeated
            % passesPerRun times, in BLOCK order -- which is obj.walkOrder, not tau order.
            % Nothing is sent to the AWG inside this loop and nothing to the pulse generator
            % either: the gate is fixed, so there is no changeSequence and no assignSequence
            % per point. That is the whole saving.
            %
            % Note what `repeats` means here: PASSES. Each tau's repeats arrive spread through
            % the run rather than back to back, which is strictly better against drift and is
            % why walkReduce exists to do the reduction processData would have done.
            N = length(obj.tau);
            pg   = getObjByName(PulseGenerator.NAME);
            spcm = getObjByName(Spcm.NAME);
            tracker = getObjByName(Tracker.NAME);
            if isempty(tracker); throwBaseObjException(Tracker.Name); end

            % A new visit order for this average, if asked (newOrderEveryAverage).
            obj.newWalkOrderIfAsked();

            nDet = obj.detectionPeriodsPerRepeat;
            passesWanted = obj.repeats;
            P = obj.passesPerRunUsed;
            if isempty(P); P = obj.choosePassesPerRun(false); obj.passesPerRunUsed = P; end
            nRuns = ceil(passesWanted / P);
            nWarm = obj.walkWarmupTriggers;

            raw = nan(nDet, N, nRuns * P);
            passesDone = 0; failedPasses = 0;
            opts = obj.walkReduceOptionsEcho(spcm);

            for r = 1:nRuns
                if obj.checkEmergencyStop(); return; end
                success = false;
                for trial = 1:5
                    try
                        if obj.resyncEveryPass
                            % Put the sequencer back on step 1 so this pass is known to start
                            % on the first block. A stop-then-run is far cheaper than an asset
                            % change, and it bounds any trigger miscount to a single pass
                            % instead of corrupting the tau mapping for everything after it.
                            obj.AWG.sendCommand('awgcontrol:stop:immediate');
                            obj.AWG.Run;
                        else
                            obj.AWG.EnsureRunning;
                        end

                        data = obj.getRawData(pg, spcm);
                        % walkUnscramble, not four lines here: this is the one place a walk can
                        % corrupt data SILENTLY -- a wrong reordering gives a shuffled Echo
                        % that still fits and still decays, and returns a confident wrong T2.
                        % So it lives in a function with an offline test. See testEchoWalk.
                        raw(:, :, passesDone + (1:P)) = ...
                            walkUnscramble(data, nDet, nWarm, N, P, obj.walkOrder);
                        passesDone = passesDone + P;
                        success = true;
                        break
                    catch err
                        err2warning(err);
                        failedPasses = failedPasses + 1;
                        fprintf('Walk pass %d failed at trial %d, attempting again.\n', r, trial);
                        try
                            spcm.stopExperimentCount;
                        catch
                            % The resources may not have needed clearing, which is fine.
                        end
                    end
                end
                if ~success
                    obj.emergencyStop; return
                end

                % Reduce what has arrived so far, so the plot fills in as passes accumulate
                % rather than only at the end of the average.
                usable = min(passesDone, passesWanted);
                [sig, err] = walkReduce(raw(:, :, 1:usable), nDet, N, opts);
                obj.signal(:, :, obj.currIter) = sig;
                obj.sterr(:, :, obj.currIter)  = err;

                obj.publishEcho();

                obj.currParamIter = min(N, round(N * usable / passesWanted));
                if r == 1 || mod(r, 5) == 0 || r == nRuns
                    sendEventParamIterationDone(obj);
                end

                % Tracking once per pass rather than once per point: there is no per-point
                % pause left to hang it on, and the whole-sweep mean reference is a better
                % brightness metric than any single tau's.
                if obj.isTracking
                    if obj.checkEmergencyStop(); return; end
                    if tracker.compareReference(mean(sig(2, :)), mean(err(2, :)), ...
                            Tracker.REFERENCE_TYPE_KCPS, obj.trackThreshhold)
                        tracker.trackUsing(TrackablePosition.NAME, obj)
                        obj.prepare;    % before the next pass
                    end
                end
            end

            if failedPasses > 0
                if obj.resyncEveryPass
                    fprintf(['Walk mode: %d of %d passes retried. resyncEveryPass is on, so each ' ...
                             'retry re-armed from block 1 and the tau mapping is intact.\n'], ...
                            failedPasses, nRuns);
                else
                    fprintf(['Walk mode WARNING: %d of %d passes retried with resyncEveryPass OFF. ' ...
                             'An aborted pass leaves the sequencer on an unknown tau and nothing ' ...
                             'puts it back, so every later pass can be shifted. Turn it on.\n'], ...
                            failedPasses, nRuns);
                end
            else
                fprintf('Walk mode: %d passes, none retried.\n', nRuns);
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

