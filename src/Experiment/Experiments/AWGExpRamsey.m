classdef AWGExpRamsey < ExperimentAWG
    %EXPECHO Echo experiment

    properties (Constant)
        NAME = 'AWG Ramsey';
    end
    
    properties
        % All of the parameters must have a set function, and when they are
        % changed, the parameter obj.changeFlag needs to be set to true
        % (unless the change is minor and the experiment can continue
        % running);
        tau                 % in us 
        
        doubleMeasurement   % logical

        % ---- sequence walk: the fast path -------------------------------------
        useSequenceWalk = true;     % logical. (default) Arm ONCE and let the sequencer advance one tau
                                    % per pg repeat, instead of arm-run-read per tau.
                                    %
                                    % WHAT IT REMOVES, measured by ramseyStageTimer on this
                                    % setup: 0.588 s per tau, of which awg.Run is only 38%.
                                    % The rest is pg.changeSequence twice, assignSequence per
                                    % channel and the digitizer arm -- walk deletes ALL of it,
                                    % not just the re-arm. One tau's overhead is worth about
                                    % 4000 shots at 148 us each.
                                    %
                                    % REQUIRES A FIXED MW GATE -- see prepare.
        passesPerRun = [];          % full passes over the tau list per pg.run. [] = derived
                                    % from maxDigitizerSegments.
        passesPerRunUsed = [];      % what the auto choice settled on.
        walkWarmupTriggers = 1;     % throwaway blocks at the START of each pass, so the first
                                    % measured block is not the one that followed an idle
                                    % sequencer.
        resyncEveryPass = true;     % stop and re-arm at the start of every pass, so a pass
                                    % always begins on block 1. Cheap, and it bounds a trigger
                                    % miscount to one pass instead of shifting the tau mapping
                                    % for the rest of the run.
        maxDigitizerSegments = 32000;
                                    % ceiling on nDet * blocksPerPass * passesPerRun for ONE
                                    % pg.run. Measured by ramseySpeedProbe: the digitizer
                                    % accepted 32000 and the probe's try-list ran out before
                                    % the instrument refused.
                                    %
                                    % WHY NOT THE SAFER 12000. Because at 12000 walk is worse
                                    % than useless. Walk wins only when
                                    %
                                    %       passesPerRun  >  repeats / nTau
                                    %
                                    % since classic pays nTau arms and walk pays
                                    % ceil(repeats/P). At 34 taus and 3000 repeats that
                                    % threshold is P > 88, and 12000 gives P = 85: walk would
                                    % use 36 arms against classic's 34. Shipping a fast path
                                    % whose default setting is slower than the thing it
                                    % replaces is worse than not shipping it.
                                    %
                                    % ACCEPTANCE IS NOT CORRECTNESS, though, and that is the
                                    % honest caveat: the probe showed the digitizer allocating
                                    % 32000 scans, not that it returns good data for all of
                                    % them. THE FIRST WALK RUN MUST BE COMPARED AGAINST A
                                    % CLASSIC ONE at the same settings. If they disagree,
                                    % halve this and compare again -- that bisection finds the
                                    % real ceiling in two or three runs.
        walkOrder = [];             % the order the blocks were built in, so the data can be
                                    % put back into tau order.
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

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function spec = fitSpec(obj)
            % The experiment window's fit: ramseyFitRobust (see ExpFit)
            spec = ExpFit.ramsey(obj);
        end

        function names = rawWindowNames(~, ~)
            % Rows of obj.signal, for the window's raw view (see publishRamsey)
            names = {'|1> signal', '|1> reference', '|0> signal', '|0> reference'};
        end

        function obj = AWGExpRamsey()
            obj@ExperimentAWG(AWGExpRamsey.NAME);
            obj.parameterName = 'taus';
            
            % Set properties inherited from Experiment
            obj.repeats = 5000;
            obj.averages = 10;
            
            obj.frequency = 2687;           % in MHz
            obj.amplitude = -15;            % in dBm
            obj.tau = 0.001:0.01:2;                   
            obj.halfPiTime = 0.0198;         % in us
            obj.lastDelay = Experiment.DEFAULT_LAST_DELAY;
            obj.AWGorSRSswitch = 0;

            obj.constantTime = false;       % logical
            obj.doubleMeasurement = true;   % logical
            
            obj.detectionDuration = 0.5;           % detection window, in us
            obj.referenceDetectionDuration = 0.5;     % in us. Detection duration of the reference read
            obj.laserInitializationDuration = 10;   % laser initialization in pulsed experiments in us
            
            obj.mCurrentXAxisParam = ExpParamDoubleVector('Time', [], [], StringHelper.MICROSEC, obj.NAME);
            % Labels as in the robust fit's figure (ExpPlotStyle / fitStyle)
            obj.plotXLabel = ExpPlotStyle.X_RAMSEY;
            obj.plotYLabel = ExpPlotStyle.Y_NORMALISED;
            obj.plotYLabelAlternate = ExpPlotStyle.Y_CONTRAST;   % alternateSignal = |0> - |1>
            obj.topParam = ExpParamDoubleVector(StringHelper.TAU, [], [], StringHelper.MICROSEC, obj.NAME);
            obj.signalParam = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME, '|1>');
            obj.signalParam2 = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME, '|0>');
        end
    end
    
    properties (SetAccess = protected)
        walkParts = []          % struct. What uploadWalkSequence needs to rebuild the walk
                                % sequence in a new order without touching the waveforms.
        walkOrderFresh = false  % true while the loaded order has not been used by an average
                                % yet, so newOrderEveryAverage does not redraw it.
        seqAssetIndex               % the ONE sequence asset per channel, in walk mode
    end

    %% Setters
    methods

        function set.useSequenceWalk(obj, newVal)
            if ~ValidationHelper.isTrueOrFalse(newVal)
                obj.sendError('useSequenceWalk must be true or false')
            end
            obj.useSequenceWalk = logical(newVal);
            obj.changeFlag = true;      % the SEQUENCE layout changes, so it must reload
        end

        function set.passesPerRun(obj, newVal)
            if ~isempty(newVal) && ~ValidationHelper.isValuePositiveInteger(newVal)
                obj.sendError('passesPerRun must be a positive integer, or [] to choose it automatically')
            end
            obj.passesPerRun = newVal;
            obj.changeFlag = true;      % pg.repeats is derived from it
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
            detectionSequnceDuration = obj.detectionDuration + obj.stabilizationDuration + obj.acquisitionDuration + obj.GIresetDuration;
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
                S.addEvent(obj.tau(end) + 2*obj.halfPiTime, MW,                    'tau');
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
                    S.addEvent(obj.tau(end),                    '',                         'MW');               % MW
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

            if obj.useSequenceWalk
                % THE MW GATE IS LEFT FIXED AT max(tau), AND THAT IS THE ONE PHYSICS
                % DECISION IN WALK MODE.
                %
                % The sequence above was already built with
                %     S.addEvent(obj.tau(end) + 2*obj.halfPiTime, MW, 'tau')
                % i.e. the longest tau. The classic path rewrites that gate per point with
                % pg.changeSequence and pads lastDelay to keep the total constant. Walk
                % cannot: the pulse generator fires the same repeat over and over and never
                % changes between blocks, so one gate has to suit every tau.
                %
                % It can, because the AWG WAVEFORM is what defines the pulses here: each
                % waveform holds both pi/2 pulses at their tau spacing and then zeros, so the
                % NV sees the right sequence regardless of how long the gate stays open. The
                % waveform ends with minWaveformTime of zeros, so the AWG idles at zero for
                % the rest of the gate rather than holding a carrier.
                %
                % WHAT YOU ARE TRADING. With the gate fixed, a SHORT tau finishes its second
                % pi/2 early and then waits out the remainder of the gate before the readout
                % window opens. That wait is tau-dependent, up to max(tau). But after the
                % second pi/2 the state is a POPULATION, not a coherence, so it decays at T1
                % (milliseconds) and not T2* (microseconds): 2.5 us against a 3 ms T1 is under
                % 0.1%, far below the Ramsey contrast. It is a real assumption though, and the
                % check is one run each way -- walk against classic, same settings.
                %
                % A side benefit: the sequence duration is now identical for every tau with no
                % lastDelay juggling, so constantTime has nothing left to do and is ignored.
                pg   = getObjByName(PulseGenerator.NAME);
                spcm = getObjByName(Spcm.NAME);
                N = length(obj.tau);
                nDet = obj.detectionPeriodsPerRepeat;
                blocksPerPass = N + obj.walkWarmupTriggers;
                obj.passesPerRunUsed = obj.choosePassesPerRun(true);
                repeatsPerRun = blocksPerPass * obj.passesPerRunUsed;

                % ONE pg repeat = ONE tau block, because the block is two trigger-gated
                % steps and a repeat fires exactly two AWG triggers (the measurement half and
                % the doubleMeasurement half; the balanced halves carry no trigger).
                pg.repeats = repeatsPerRun;
                numScans = nDet * repeatsPerRun;
                seqTime  = pg.sequenceDuration() * 1e-6;
                spcm.prepareExperimentCount(numScans, 2 * seqTime * repeatsPerRun, 1);
                fprintf(['Walk mode: %d pg repeats per run (%d taus + %d warm-up, x %d ' ...
                         'passes), %d digitizer scans. MW gate fixed at %.3f us.\n'], ...
                        repeatsPerRun, N, obj.walkWarmupTriggers, obj.passesPerRunUsed, ...
                        numScans, obj.tau(end) + 2*obj.halfPiTime);
            end

            % Set parameter, for saving
            obj.mCurrentXAxisParam.value = obj.tau;
            obj.isPlotAlternateAvailable = obj.doubleMeasurement;
        end

        function P = choosePassesPerRun(obj, announce)
            % How many passes to put in one pg.run, DERIVED rather than guessed.
            %
            % Walk re-arms once per pg.run rather than once per tau, so the whole average
            % costs ceil(repeats/P) arms instead of nTau -- P as large as possible is always
            % better. The limit is that one pg.run must fit in the digitizer: it returns
            % nDet * (nTau + warmup) * P scans, and nDet is FOUR here where Rabi has two, so
            % the same ceiling buys Ramsey half the passes.
            if ~exist('announce', 'var'); announce = false; end
            N = length(obj.tau);
            nDet = obj.detectionPeriodsPerRepeat;
            if isempty(nDet) || nDet < 1; nDet = 2 * (1 + double(obj.doubleMeasurement)); end
            blocksPerPass = N + obj.walkWarmupTriggers;
            if ~isempty(obj.passesPerRun)
                P = obj.passesPerRun;                       % explicit wins
            else
                P = floor(obj.maxDigitizerSegments / (nDet * blocksPerPass));
                P = max(1, min(P, obj.repeats));
            end
            if announce
                fprintf(['Walk: passesPerRun = %d (%s), %d scans per pg.run of %d allowed.\n' ...
                         '      arms this average: %d walk vs %d classic.\n'], ...
                        P, BooleanHelper.ifTrueElse(isempty(obj.passesPerRun), 'auto', 'set by hand'), ...
                        nDet*blocksPerPass*P, obj.maxDigitizerSegments, ...
                        ceil(obj.repeats/P), N);
                % The break-even condition, stated rather than implied: walk pays only when
                % P > repeats/nTau, because classic costs nTau arms and walk costs
                % ceil(repeats/P). It is easy to be the wrong side of this without noticing.
                if ceil(obj.repeats/P) >= N
                    fprintf(['\n      WALK IS NOT HELPING AT THESE SETTINGS: %d arms against\n' ...
                             '      classic''s %d. Walk pays only when passesPerRun exceeds\n' ...
                             '      repeats/nTau = %.0f, and P is %d.\n' ...
                             '      Raise maxDigitizerSegments (now %d), lengthen the tau list,\n' ...
                             '      or set useSequenceWalk = false. Running anyway.\n\n'], ...
                            ceil(obj.repeats/P), N, obj.repeats/max(N,1), P, ...
                            obj.maxDigitizerSegments);
                end
            end
        end
    
        
        function LoadAWG(obj)
            obj.AWG.ResetStoredData;
            obj.walkParts = [];     % set again below on the walk path
            obj.AWG.workChannels = obj.activeChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            minWaveformTime = obj.AWG.minWaveformDuration;
                        
            obj.halfPiTime = obj.RoundDurationByFrequencies(obj.halfPiTime , obj.AWG.SampleRate);

            N = length(obj.tau);
            A = obj.AWG.MaxNormPowerAllowed;
            w = 2*pi*obj.frequency;
            
            % creating waveforms:
            % THE COUNT WAS WRONG HERE. It said nChannels*N, but the loop below adds TWO
            % waveforms per (channel, tau) -- signal and reference -- so the true count is
            % 2*nChannels*N. That is why the old output said "Loading 334 waveforms:" and
            % then "Loading 668 waveforms to the AWG" two lines later.
            obj.AWG.progressBegin(sprintf(['Ramsey waveforms: %d tau x 2 (signal, reference) ' ...
                'x %d channel(s).  Each one holds both pi/2 pulses, so tau is encoded in the ' ...
                'SAMPLES and the waveform lengthens with it.'], N, obj.nChannels), ...
                2*obj.nChannels*N, 'waveform');
            waveformsIndex = zeros(obj.nChannels, N);
            waveformsRefIndex = zeros(obj.nChannels, N);

            [prePulseTime, phase] = obj.channelsDifference;
            % THE CALIBRATED RELATIVE PHASE, added to whatever channelsDifference produced.
            % Without this line channelPhase is a property nothing reads: channelsDifference
            % builds its phase from timeDelay alone, and level 4.3 sets timeDelay = 0 -- so
            % the two channels would be driven with ZERO relative phase while piTime and
            % halfPiTime in the json were measured with the calibrated phase applied. The
            % pulses would then not be pi and pi/2 at all, and nothing would say so.
            phase = phase + obj.channelPhaseRadians();
            for chan = 1:obj.nChannels
                for i = 1:N
                    f = @(t) A(chan) * ((prePulseTime(chan)                           <= t & t < prePulseTime(chan)+obj.halfPiTime)              .* cos(w*t - phase(chan))       + ...
                                        (prePulseTime(chan)+obj.halfPiTime+obj.tau(i) <= t & t < prePulseTime(chan)+2*obj.halfPiTime+obj.tau(i)) .* cos(w*t - phase(chan)));
                    waveformsIndex(chan,i)= obj.AWG.AddWaveformByFunction(   f, prePulseTime(chan)+2*obj.halfPiTime+obj.tau(i)+minWaveformTime, 0, 1);
                    f = @(t) A(chan) * ((prePulseTime(chan)                           <= t & t < prePulseTime(chan)+obj.halfPiTime)              .* cos(w*t - phase(chan))       + ...
                                        (prePulseTime(chan)+obj.halfPiTime+obj.tau(i) <= t & t < prePulseTime(chan)+2*obj.halfPiTime+obj.tau(i)) .* cos(w*t - phase(chan) + pi));
                    waveformsRefIndex(chan,i)= obj.AWG.AddWaveformByFunction(f, prePulseTime(chan)+2*obj.halfPiTime+obj.tau(i)+minWaveformTime, 0, 1);
                end
            end
            obj.AWG.progressEnd();
            obj.AWG.LoadWaveform();

            if obj.useSequenceWalk
                % Kept apart from the sequence so newOrderEveryAverage can rebuild it in a
                % new order without re-uploading a single waveform.
                p = struct();
                p.waveformsIndex    = waveformsIndex;
                p.waveformsRefIndex = waveformsRefIndex;
                p.nWarm             = obj.walkWarmupTriggers;
                obj.walkParts = p;
                nSteps = obj.uploadWalkSequence(randperm(N));
                fprintf(['Walk load: %d taus + %d warm-up, %d steps per channel, %d AddSequence ' ...
                         'calls total (the classic path would make %d).\n'], ...
                        N, p.nWarm, nSteps, obj.nChannels, obj.nChannels*N);
                return
            end

            % creating sequences:
            obj.indexSeq = zeros(obj.nChannels, N);
            obj.AWG.progressBegin(sprintf(['Ramsey sequences: one per tau per channel, two ' ...
                'steps each (signal then reference), both trigger-gated.']), ...
                obj.nChannels*N, 'sequence');
            for chan = 1:obj.nChannels
                for i = 1:N
                    obj.indexSeq(chan,i) = obj.AWG.AddSequence([waveformsIndex(chan,i), waveformsRefIndex(chan,i)], [1 1], [1 1]);
                    obj.AWG.progressStep();
                end
            end
            obj.AWG.progressEnd();
            obj.AWG.LoadSequence()
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
            % Each tau is a two-step block -- signal waveform then reference waveform, both
            % trigger-gated -- which is byte for byte the block the classic path puts in its
            % own little sequence. Concatenating them changes nothing about what the NV sees;
            % it changes who advances between taus. Classic: MATLAB, via assignSequence, once
            % per point. Walk: the sequencer itself, on the trigger that was going to be sent
            % anyway.
            %
            % The last step points back at step 1, so a pass ends where the next one begins
            % and pass two needs no re-arm.

            % ---- the block order, randomised once ----------------------------
            % The classic path draws a fresh randperm every average. Walk cannot: the visit
            % order IS the sequence order. So the randomisation moves here, and the order is
            % recorded so performWalk can put the data back in tau order. Without it, slow
            % drift over a pass would appear as a systematic trend in tau -- on a Ramsey that
            % is indistinguishable from the T2* envelope we are trying to measure.
            % (The order now comes in as 'order': drawn per load by LoadAWG and, with
            % newOrderEveryAverage, redrawn before every later average.)
            p = obj.walkParts;
            waveformsIndex = p.waveformsIndex; waveformsRefIndex = p.waveformsRefIndex;
            N = numel(order);
            nWarm = p.nWarm;
            blockOrder = [repmat(order(1), 1, nWarm), order];
            nSteps = 2 * numel(blockOrder);

            if nSteps > obj.AWG.MAX_SEQUENCE_STEPS
                obj.sendError(sprintf(['Walk needs %d sequence steps (2 per tau, %d taus plus ' ...
                    '%d warm-up) and the AWG allows %d. Shorten the tau list or use classic ' ...
                    'mode.'], nSteps, N, nWarm, obj.AWG.MAX_SEQUENCE_STEPS));
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
            obj.AWG.progressBegin(sprintf(['Ramsey WALK sequence: ONE asset per channel, %d ' ...
                'steps (2 per tau x %d taus + %d warm-up block(s)). The sweep lives in the ' ...
                'sequence, so the AWG is armed once for the whole pass.'], ...
                nSteps, N, nWarm), obj.nChannels, 'sequence');
            for chan = 1:obj.nChannels
                wfList = zeros(1, nSteps); rptList = ones(1, nSteps);
                trgList = ones(1, nSteps); gotoList = zeros(1, nSteps);
                s0 = 0;
                for blk = 1:numel(blockOrder)
                    i = blockOrder(blk);
                    wfList(s0+1:s0+2) = [waveformsIndex(chan,i), waveformsRefIndex(chan,i)];
                    s0 = s0 + 2;
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

            % ---- assign and arm ONCE for the whole sweep ---------------------
            for chan = 1:obj.nChannels
                obj.AWG.assignSequence(obj.seqAssetIndex(chan), obj.activeChannels(chan));
            end
            obj.AWG.Run;
            obj.walkOrder = order;
            obj.walkOrderFresh = true;      % not used by an average yet
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
            % One average, acquired as whole passes over the tau list.
            %
            % The sequencer advances one BLOCK per pg repeat, so a pg.run of
            % (N + warmup) * passesPerRun repeats returns the whole sweep, repeated
            % passesPerRun times, in BLOCK order -- which is obj.walkOrder, not tau order.
            % Nothing is sent to the AWG inside the loop and nothing to the pulse generator
            % either: the gate is fixed, so there is no changeSequence per point. That is the
            % 0.588 s per tau ramseyStageTimer measured, gone.
            %
            % Note what `repeats` means here: PASSES. Each tau's repeats arrive spread out in
            % time rather than back to back, which is strictly better for drift and is why
            % walkReduce exists to do the same reduction processData would.
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
            blocksPerPass = N + nWarm;
            readsPerRun = nDet * blocksPerPass * P;

            raw = nan(nDet, N, nRuns * P);
            passesDone = 0; failedPasses = 0;
            opts = obj.walkReduceOptionsRamsey(spcm);

            for r = 1:nRuns
                if obj.checkEmergencyStop(); return; end
                success = false;
                for trial = 1:5
                    try
                        if obj.resyncEveryPass
                            % Put the sequencer back on step 1 so this pass is known to start
                            % on the first block. A stop-then-run is far cheaper than an asset
                            % change, and it bounds any trigger miscount to a single pass
                            % instead of corrupting the tau mapping for the rest of the run.
                            obj.AWG.sendCommand('awgcontrol:stop:immediate');
                            obj.AWG.Run;
                        else
                            obj.AWG.EnsureRunning;
                        end

                        data = obj.getRawData(pg, spcm);
                        % walkUnscramble, not four lines here: this is the one place a walk
                        % can corrupt data silently -- a wrong reordering gives a shuffled
                        % Ramsey that still fits and still decays -- so it lives in a
                        % function with an offline test. See testRamseyWalk.
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

                obj.publishRamsey();

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

        function publishRamsey(obj)
            % The two ratios, exactly as the classic path forms them at the end of perform.
            % Factored out so walk and classic cannot drift apart: rows 1-2 are the |1>
            % projection and rows 3-4 the |0> one.
            S1 = squeeze(obj.signal(1, :, 1:obj.currIter));
            S1e = squeeze(obj.sterr(1, :, 1:obj.currIter));
            S2 = squeeze(obj.signal(2, :, 1:obj.currIter));
            S2e = squeeze(obj.sterr(2, :, 1:obj.currIter));
            [value, sterrOut] = getRatioDistributionValues(obj, S1, S2, S1e, S2e);
            obj.signalParam.value = value;
            obj.signalParam.sterr = sterrOut;

            if obj.doubleMeasurement && size(obj.signal, 1) >= 4
                S3 = squeeze(obj.signal(3, :, 1:obj.currIter));
                S3e = squeeze(obj.sterr(3, :, 1:obj.currIter));
                S4 = squeeze(obj.signal(4, :, 1:obj.currIter));
                S4e = squeeze(obj.sterr(4, :, 1:obj.currIter));
                [value, sterrOut] = getRatioDistributionValues(obj, S3, S4, S3e, S4e);
                obj.signalParam2.value = value;
                obj.signalParam2.sterr = sterrOut;
            else
                obj.signalParam2.value = [];
            end
        end

        function opts = walkReduceOptionsRamsey(obj, spcm)
            % Everything walkReduce needs to reproduce processData exactly. Same contents as
            % AWGExpRabi's version; duplicated rather than shared because these are two
            % experiment classes and neither should reach into the other.
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
            maxLastDelay = obj.lastDelay + max(obj.tau) + 2*obj.halfPiTime;
            
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
                        pg.changeSequence('tau', 'duration', obj.tau(t) + 2*obj.halfPiTime);
                        if obj.constantTime
                            pg.changeSequence('lastDelay', 'duration', maxLastDelay - obj.tau(t));
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
            
            % Saving results in the Experiment parameters. Same call the walk path makes,
            % so the two cannot form the ratios differently.
            obj.publishRamsey();
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

