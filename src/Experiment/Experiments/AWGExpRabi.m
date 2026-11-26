classdef AWGExpRabi < ExperimentAWG
    
    properties (Constant)
        NAME = 'AWG Rabi';
    end
    
    properties
        % All of the parameters must have a set function, and when they are changed, the parameter obj.changeFlag needs to be set to true
        % (unless the change is minor and the experiment can continue running);
        tau             % in us
        totalZero
        smallOnDelay    % in us
        smallOffDelay   % in us
        dontLoadAWG

        % ---- WALK MODE, ported from AWGExpESR -----------------------------
        % One AWG sequence holds EVERY tau as a consecutive block, the first step of
        % each block waits for a trigger, and the sequencer advances one block per
        % trigger. So a single pg.run walks the whole sweep and there is no per-point
        % AWG traffic and no per-point re-arm at all.
        %
        % Two costs disappear together, and the load is the bigger of them. The classic
        % path calls AddSequence once per tau per channel -- 194 calls for 97 taus, which
        % measured 249 s of a 265 s run, at about 173 ms per SCPI command. Walk mode calls
        % it TWICE, once per channel, with the whole step list. That is the same change
        % that took the ESR load from 130 s to 6 s.
        useSequenceWalk = true;     % logical. The fast path (default). Requires a FIXED MW gate --
                                    % see prepare, and read the note there before using it
                                    % on data you care about.
        passesPerRun = [];          % full passes over the tau list per pg.run.
                                    % [] = CHOOSE IT AUTOMATICALLY, which is the sane
                                    % default: the right value follows from repeats, the
                                    % tau count and the digitizer depth, and 25 was a
                                    % number I picked by hand with no reasoning behind it.
                                    % See choosePassesPerRun.
                                    % THIS IS THE KNOB THAT DECIDES WHETHER WALK WINS,
                                    % and shipping it at 1 is what made walk lose its
                                    % first comparison. resyncEveryPass re-arms once per
                                    % pg.run, so nRuns = ceil(repeats/passesPerRun) is the
                                    % number of re-arms for the whole average. At
                                    % passesPerRun = 1 and repeats = 300 that is 300
                                    % re-arms, against 97 for the classic path's one per
                                    % tau -- three times worse, measured as 88 s against
                                    % 69 s. At 25 it is 12.
                                    %
                                    % The ceiling is the digitizer: one pg.run has to hold
                                    % 2*(nTau+warmup)*passesPerRun segments, so 25 passes
                                    % over 97 taus is 4900. Raise it until
                                    % prepareExperimentCount complains, then back off one
                                    % step. 10 is a deliberately safe starting point, not
                                    % an optimum.
        channelPhase = [0 0];       % DEGREES, per AWG channel. The relative phase of the
                                    % two microwave drives, written into the waveform rather
                                    % than into the length of a padding waveform. See where
                                    % it is applied for why timeDelay cannot do this job.
                                    % Set timeDelay = 0 when using it.

        maxDigitizerSegments = 5000; % the ceiling on 2*(nTau+warmup)*passesPerRun for one
                                    % pg.run. 5000 is not a datasheet number -- it is the
                                    % largest value MEASURED to work here (99 taus at 25
                                    % passes ran clean). Raise it until
                                    % prepareExperimentCount complains, then back off, and
                                    % put the answer here so the auto choice can use it.
        passesPerRunUsed = [];      % what the auto choice settled on, for performWalk and
                                    % for the record.
        walkWarmupTriggers = 1;     % throwaway blocks at the START of each pass, so the
                                    % first real tau is not the trigger that follows the
                                    % arm/idle gap. Their reads are discarded.
        resyncEveryPass = true;     % stop and re-arm at the start of every pass, so the
                                    % sequencer is known to begin on block 1. A retry
                                    % aborts a pass part-way and otherwise leaves the
                                    % sequencer on an arbitrary tau with nothing to put
                                    % it back -- which silently shifts every later pass.
        seqAssetIndex = [];         % walk mode: the ONE sequence asset per channel. Cached,
                                    % because a cache hit has to re-assign and re-arm it.
        walkOrder = [];             % the order the blocks were built in, so the data can
                                    % be mapped back to tau. Randomised per load: with
                                    % walk the visit order is the SEQUENCE order, so the
                                    % randperm the classic path does per average has to
                                    % move into the sequence itself or slow drift becomes
                                    % a trend in tau.
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
    end

    properties (SetAccess = protected)
        walkParts = []          % struct. What uploadWalkSequence needs to rebuild the walk
                                % sequence in a new order without touching the waveforms.
        walkOrderFresh = false  % true while the loaded order has not been used by an average
                                % yet, so newOrderEveryAverage does not redraw it.
        walkCacheKey = ''       % the AWG load-cache key the walk sequence is stamped under,
                                % so a per-average re-order can re-stamp the SAME entry.
    end

    methods
        function obj = AWGExpRabi()
            obj@ExperimentAWG(AWGExpRabi.NAME);
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
            
            obj.detectionDuration = 0.25;               % detection window, in us
            obj.referenceDetectionDuration = 5;         % in us. Detection duration of the reference read
            obj.laserInitializationDuration = 10;       % laser initialization in pulsed experiments
            
            obj.mCurrentXAxisParam = ExpParamDoubleVector('Time', [], [], StringHelper.MICROSEC, obj.NAME);
            % Labels as in the robust fit's figure (ExpPlotStyle / fitStyle)
            obj.plotXLabel = ExpPlotStyle.X_RABI;
            obj.plotYLabel = ExpPlotStyle.Y_NORMALISED;
            obj.signalParam = ExpResultDoubleVector('FL', [], [], 'Normalized', obj.NAME);
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

        function set.useSequenceWalk(obj, newVal)
            if ~ValidationHelper.isTrueOrFalse(newVal)
                obj.sendError('useSequenceWalk must be true or false')
            end
            obj.useSequenceWalk = logical(newVal);
            obj.changeFlag = true;      % the sequence layout itself changes
        end

        function set.passesPerRun(obj, newVal)
            if ~isempty(newVal) && ~ValidationHelper.isValuePositiveInteger(newVal)
                obj.sendError('passesPerRun must be a positive integer, or [] to choose it automatically')
            end
            obj.passesPerRun = newVal;
            obj.changeFlag = true;      % pg.repeats is derived from it
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
            S.addEvent(obj.tau(end),                    MW,                         'MW')               % MW
            S.addEvent(obj.lastDelay,                   MW,                         'lastDelay');       % Last delay

            S.addEvent(obj.detectionDuration,       {'greenLaser', 'detector'});                    % Detection
            S.addEvent(initDuration,                'greenLaser');                                  % Initialization
            if ~obj.balancedSequence
                S.addEvent(obj.referenceDetectionDuration,...
                    {'greenLaser', 'detector'});                                                        % Reference detection
                S.addEvent(1,                           'greenLaser');                                  % Initialization
            else
                S.addEvent(obj.AWG.TRIG_DURATION,           '');                        % AWG trigger
                S.addEvent(obj.tau(end),                    '',                         'MW');               % MW
                S.addEvent(obj.lastDelay,                   '',                         'lastDelay');       % Last delay
                S.addEvent(obj.referenceDetectionDuration,...
                    {'greenLaser', 'detector'});                                                        % Reference detection
                S.addEvent(initDuration,                'greenLaser');                                  % Initialization
            end
            
            % ---------------------------------------------------------------------
            % TURN THE MW OUTPUTS BACK ON, on every run.
            %
            % This is the counterpart of a real bug fixed in
            % FrequencyGeneratorTektronixAWG. The stop command used to be
            % 'AWGCONTROL:STOP:IMMEDIATE; OUTP1 OFF' -- a bare ';' leaves the Tektronix
            % parser inside the AWGControl subtree, so 'OUTP1 OFF' was parsed as
            % 'AWGCONTROL:OUTP1 OFF' and rejected with error -113. Nothing reads the error
            % queue, so the failure was invisible: the outputs never actually went off at
            % the end of a run, and every experiment silently relied on that.
            %
            % The separator is now ';:' and the outputs really do go off. So something has
            % to switch them back on, and prepareInternal does NOT: its fg.output = true
            % branch is skipped for a TektronixAWG by
            %   (~contains(name,'TektronixAWG') + ~contains(name,'synthhd...') == 2)
            % which is 0 + 1 = 1, never 2.
            %
            % LoadAWG below does set workChannels -- but it only runs when
            % (changeFlag || restartFlag) && ~dontLoadAWG. So a restart, a rerun with
            % nothing changed, or dontLoadAWG = 1 emits NO MICROWAVES AT ALL, and every
            % point comes back at ratio ~= 1. That is a data-integrity failure that looks
            % like a physics result.
            %
            % activeChannels, not 1:nChannels: with nChannels = 1 and MWChannel = 'MW2',
            % activeChannels is 2 and 1:nChannels would enable the wrong output.
            %
            % Two SCPI commands, idempotent, safe to send every run.
            obj.AWG.workChannels = obj.activeChannels;

            obj.prepareInternal(S)
            if (obj.changeFlag || obj.restartFlag) && (~obj.dontLoadAWG)
                obj.LoadAWG;
            end

            if obj.useSequenceWalk
                % THE MW GATE IS LEFT FIXED AT max(tau), AND THAT IS THE ONE PHYSICS
                % DECISION IN WALK MODE.
                %
                % The sequence above was already built with S.addEvent(obj.tau(end), MW,
                % 'MW'), i.e. the longest tau. The classic path then rewrites that gate
                % per point with pg.changeSequence, and pads lastDelay to keep the total
                % duration constant. Walk mode cannot: the pulse generator fires one
                % trigger per block and never changes between them, so the gate has to
                % suit every tau at once.
                %
                % It can, because in THIS class the AWG waveform is what defines the
                % pulse: each block plays initZero, then numZeroSteps whole cycles, then
                % a residue waveform, and after that the AWG output is zero for the rest
                % of the gate. So the pulse the NV sees is tau regardless of the gate.
                %
                % What you are trusting instead is that the AWG really does idle at zero
                % while the switch is still open. Leakage during that tail is AWG noise
                % floor rather than full carrier, and it is the same tail the reference
                % window sees, but it is a real assumption and it is worth one check:
                % run a known pi time in both modes and compare the contrast.
                %
                % A side benefit: the sequence duration is now identical for every tau
                % without lastDelay juggling, so the duty cycle is constant across the
                % sweep. constantTime is doing that job in the classic path.
                pg = getObjByName(PulseGenerator.NAME);
                spcm = getObjByName(Spcm.NAME);
                nT = length(obj.tau);
                blocksPerPass = nT + obj.walkWarmupTriggers;
                obj.passesPerRunUsed = obj.choosePassesPerRun(true);
                triggersPerRun = blocksPerPass * obj.passesPerRunUsed;
                pg.repeats = triggersPerRun;
                numScans = obj.detectionPeriodsPerRepeat * triggersPerRun;
                seqTime = pg.sequenceDuration() * 1e-6;
                spcm.prepareExperimentCount(numScans, 2 * seqTime * triggersPerRun, 1);
                fprintf(['Walk mode: %d triggers per pg.run (%d taus + %d warm-up, x %d ' ...
                         'passes), %d digitizer segments. MW gate fixed at %.3f us.\n'], ...
                        triggersPerRun, nT, obj.walkWarmupTriggers, obj.passesPerRunUsed, ...
                        numScans, max(obj.tau));
            end

            % Set parameter, for saving
            obj.mCurrentXAxisParam.value = obj.tau;
        end
        
        function P = choosePassesPerRun(obj, announce)
            % How many passes to put in one pg.run, DERIVED rather than guessed.
            %
            % Two facts fix it. First, walk re-arms once per pg.run, not once per tau, so
            % the whole average costs ceil(repeats/P) re-arms at ~0.6 s each -- P as large
            % as possible is always better. Second, one pg.run has to fit in the digitizer:
            % it returns 2*(nTau+warmup)*P segments. So
            %
            %       P = min( repeats, floor(maxDigitizerSegments / (2*(nTau+warmup))) )
            %
            % capped at repeats because past that there is nothing left to fit in, and
            % floored at 1. On 99 taus with maxDigitizerSegments = 5000 that gives 25 --
            % which is where the 25 I handed you came from, except that I had picked it by
            % hand and it only happened to be right for this tau count. Change nTau and the
            % hand-picked number silently becomes wrong; this does not.
            if ~exist('announce', 'var'); announce = false; end
            nT = length(obj.tau);
            blocksPerPass = nT + obj.walkWarmupTriggers;
            if ~isempty(obj.passesPerRun)
                P = obj.passesPerRun;                       % explicit wins
            else
                P = floor(obj.maxDigitizerSegments / (2*blocksPerPass));
                P = max(1, min(P, obj.repeats));
            end

            reArmsWalk    = ceil(obj.repeats / P);
            reArmsClassic = nT;
            if announce
                fprintf(['Walk: passesPerRun = %d (%s), %d segments per pg.run.\n' ...
                         '      re-arms this average: %d walk vs %d classic.\n'], ...
                        P, BooleanHelper.ifTrueElse(isempty(obj.passesPerRun), 'auto', 'set by hand'), ...
                        2*blocksPerPass*P, reArmsWalk, reArmsClassic);
                if reArmsWalk >= reArmsClassic
                    fprintf(['      *** WALK WILL BE SLOWER THAN CLASSIC HERE. Re-arms are what\n' ...
                             '      costs (~0.6 s each), and walk needs ceil(repeats/P) of them\n' ...
                             '      against classic''s one per tau. Either drop repeats to about\n' ...
                             '      %d, or raise maxDigitizerSegments if the digitizer allows it,\n' ...
                             '      or just use the classic path for this measurement.\n'], ...
                            max(P*(nT-1), 1));
                end
            end
        end

        function phi = channelPhaseRadians(obj)
            % channelPhase padded or trimmed to nChannels, in radians. Padded rather than
            % erroring: a two-element phase left in place while nChannels drops to 1 is a
            % normal thing to do mid-calibration, and it should not throw.
            phi = zeros(1, obj.nChannels);
            cp = obj.channelPhase(:)';
            k = min(numel(cp), obj.nChannels);
            phi(1:k) = deg2rad(cp(1:k));
        end

        function [fp, parts] = awgLoadFingerprint(obj)
            % Hash of everything that decides what gets uploaded, plus the PARTS it was
            % built from so a miss can be explained instead of guessed at.
            %
            % THE KEY IS THE REQUESTED STATE, NOT THE LOADED STATE, and getting that
            % backwards is why the cache never hit. LoadAWG REWRITES obj.tau: values get
            % snapped onto existing residue waveforms and duplicates are deleted. The
            % first version stamped the cache by recomputing the hash after the load, so
            % the key described the SNAPPED tau -- while the next run starts from the
            % list in the run script, hashes that, and never matches. Every run missed,
            % every run re-uploaded. Now the hash is taken once on entry and that same
            % value is what gets stamped; the snapped tau travels in the payload instead.
            %
            % AWG.SampleRate is deliberately NOT in here, for the same reason: LoadAWG
            % calls ChangeSampleRate, so a hash taken on entry and a hash taken on exit
            % would disagree about it. It is redundant anyway -- it follows from
            % baseWFtime and frequency, and frequency is in the hash.
            %
            % detectionDuration and laserInitializationDuration ARE included even though
            % this class's waveforms do not depend on them. They cost an unnecessary
            % reload if changed alone; a false hit costs a run of wrong data.
            % FIELD BY FIELD, not struct('a', x, 'b', y, ...). struct() EXPANDS a cell
            % value into a struct ARRAY -- struct('c', {1,2}) is 1x2, not a scalar with a
            % cell in it -- and one of these properties is a cell on this setup. The
            % result was a 1xN parts, and the next line's parts.tau = ... then failed with
            % "Scalar structure required for this assignment". Assigning fields one at a
            % time cannot expand, whatever the value turns out to be.
            parts = struct();
            parts.tau                = obj.tau;
            parts.frequency          = obj.frequency;
            parts.amplitude          = obj.amplitude;
            parts.nChannels          = obj.nChannels;
            parts.timeDelay          = obj.timeDelay;
            parts.amplitudeDiff      = obj.amplitudeDiff;
            parts.channelPhase       = obj.channelPhase;
            parts.detectionDuration  = obj.detectionDuration;
            parts.refDetection       = obj.referenceDetectionDuration;
            parts.laserInit          = obj.laserInitializationDuration;
            parts.useSequenceWalk    = obj.useSequenceWalk;
            parts.walkWarmupTriggers = obj.walkWarmupTriggers;
            parts.MWChannel          = obj.MWChannel;
            parts.AWGorSRSswitch     = obj.AWGorSRSswitch;
            parts.smallOnDelay       = obj.smallOnDelay;
            parts.smallOffDelay      = obj.smallOffDelay;
            fp = DataHash(parts);
        end

        function n = expectedSequenceCount(obj)
            % How many sequences the AWG should be holding. A cheap independent check
            % that a cache hit is telling the truth, rather than trusting the hash alone.
            if obj.useSequenceWalk
                n = obj.nChannels;                          % one concatenated sequence each
            else
                n = length(obj.tau) * obj.nChannels;        % one per (tau, channel)
            end
        end

        function [tf, fp] = tryLoadFromCache(obj)
            % True when the upload was skipped. Also returns the ENTRY fingerprint, which
            % is what stampLoaded must record -- see awgLoadFingerprint.
            tf = false;
            [fp, parts] = obj.awgLoadFingerprint();
            if ~obj.AWG.isAlreadyLoaded(fp)
                % SAY WHY IT MISSED. A hash tells you nothing about which field moved, and
                % "it still loads every time" cost a round trip to diagnose. Comparing the
                % stored parts against the current ones names the culprit immediately.
                cached = obj.AWG.getLoadedCache();
                if isstruct(cached) && isfield(cached, 'parts')
                    obj.reportCacheMiss(cached.parts, parts);
                else
                    fprintf('AWG load cache: empty (first load of this session).\n');
                end
                return
            end

            % Verify against the instrument rather than trusting the hash alone. One
            % query, a couple of ms, and it catches a stale or wrong-mode cache before it
            % can silently produce garbage.
            cached = obj.AWG.getLoadedCache();
            nLoaded = str2double(obj.AWG.queryRaw('SLISt:SIZE?', 10));
            if obj.useSequenceWalk
                nWanted = obj.nChannels;
            else
                % Against the CACHED tau count, not the requested one: the loaded
                % sequences correspond to the snapped list, which can be shorter.
                nWanted = numel(cached.tau) * obj.nChannels;
            end
            if nLoaded ~= nWanted
                fprintf(['AWG holds %d sequences but this mode needs %d -- the load cache ' ...
                         'was stale. Reloading.\n'], nLoaded, nWanted);
                obj.AWG.invalidateLoaded();
                return
            end

            % Restore the DERIVED state. obj.tau has to become the snapped list the
            % loaded waveforms actually correspond to, or the x axis will not match the
            % data and every later fit is against slightly wrong times.
            obj.tau       = cached.tau;
            obj.indexSeq  = cached.indexSeq;
            if isfield(cached, 'seqAssetIndex'); obj.seqAssetIndex = cached.seqAssetIndex; end
            if isfield(cached, 'walkOrder');     obj.walkOrder     = cached.walkOrder;     end
            % What a per-average re-order needs. An older cache has no walkParts, and then
            % newWalkOrderIfAsked keeps the loaded order and says so.
            obj.walkParts = [];
            if isfield(cached, 'walkParts');     obj.walkParts     = cached.walkParts;     end
            obj.walkOrderFresh = false;          % this order was used by an earlier run
            obj.walkCacheKey   = fp;

            % Amplitude and the channel enables live in the instrument front end, not in
            % the waveform data, so re-apply them rather than assume they survived.
            obj.AWG.workChannels = obj.activeChannels;
            obj.AWG.setAmplitude(obj.amplitude);

            if obj.useSequenceWalk
                % Walk arms once and never touches the AWG again during the sweep, so a
                % cache hit still has to assign and arm.
                for chan = 1:obj.nChannels
                    obj.AWG.assignSequence(obj.seqAssetIndex(chan), obj.activeChannels(chan));
                end
                obj.AWG.Run;
            end
            fprintf('AWG already holds this waveform set (%d sequences) - skipping upload.\n', nLoaded);
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
                fprintf(['AWG load cache MISSED but every hashed field matches -- the hash ' ...
                         'itself is the problem, not the settings.\n']);
            else
                fprintf('AWG load cache miss, because: %s\n', strjoin(changed, ', '));
            end
        end

        function stampLoaded(obj, fpAtEntry)
            % Stamp the cache with the fingerprint taken ON ENTRY to LoadAWG, and carry
            % the derived state as the payload.
            %
            % Recomputing the hash here instead is what broke it: by this point obj.tau
            % has been snapped and possibly shortened, so the key would describe a tau
            % list no run script ever asks for.
            [~, parts] = obj.awgLoadFingerprint();
            if ~exist('fpAtEntry', 'var') || isempty(fpAtEntry)
                fpAtEntry = DataHash(parts);        % only reachable if a caller forgot
            end
            % parts here already describes the POST-load state, including the snapped tau,
            % because awgLoadFingerprint reads obj.* as they are now. It is stored only so
            % a future miss can be explained; the KEY is fpAtEntry.
            cache = struct();
            cache.indexSeq      = obj.indexSeq;
            cache.tau           = obj.tau;
            cache.seqAssetIndex = obj.seqAssetIndex;
            cache.walkOrder     = obj.walkOrder;
            cache.walkParts     = obj.walkParts;
            cache.parts         = parts;
            obj.AWG.markLoaded(fpAtEntry, cache);
            obj.walkCacheKey    = fpAtEntry;
        end

        function LoadAWGWalk(obj, fpAtEntry)
            % ONE sequence per channel, every tau a consecutive block.
            %
            % The classic LoadAWG below calls AddSequence once per tau per channel. On 97
            % taus and 2 channels that is 194 calls, and it measured 249 s of a 265 s run
            % -- about 173 ms per SCPI command. This calls AddSequence twice, once per
            % channel, with the whole step list. The waveforms are identical to the
            % classic path; only the sequence assembly differs.
            %
            % Block layout, matching AWGExpESR's walk blocks:
            %   [initZero, constStep x numZeroSteps, residue]   trigger flags [1 0 0]
            %   [initZero, residue]                             trigger flags [1 0]
            % The first step of each block WAITS for a trigger, the rest play straight
            % through, so one trigger plays exactly one tau and the sequencer is left
            % sitting on the next block.
            %
            % Unlike the ESR version the blocks have DIFFERENT step counts, because
            % numZeroSteps varies with tau. That is fine for walking, which only needs the
            % blocks in order; it would not be fine for jump mode, which needs a known
            % step index per block to jump to.
            baseWFtime = 0.2;
            baseWFtime = obj.RoundDurationByFrequencies(baseWFtime, obj.frequency);
            obj.AWG.ResetStoredData;
            obj.AWG.workChannels = obj.activeChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            obj.AWG.ChangeSampleRate(baseWFtime);
            minWaveformTime = obj.AWG.minWaveformDuration;

            N = length(obj.tau);
            A = obj.AWG.MaxNormPowerAllowed;
            w = 2*pi*obj.frequency;
            [prePulseTime, phase] = obj.channelsDifference;
            % channelPhase, in RADIANS, added to whatever channelsDifference produced.
            %
            % channelsDifference lives in ExperimentAWG (shared, not ours to touch) and the
            % `phase` it returns is only the SUB-SAMPLE REMAINDER of timeDelay after the AWG
            % class rounds that delay to whole samples. It is not a phase knob and cannot be
            % used as one: at 16 GS/s one sample is 62.5 ps and one period at 2870 MHz is
            % 348.4 ps, so the settable delays land about 64.6 degrees apart. Roughly 0, 65,
            % 129, 194, 258, 323 -- about 5.6 points per period, best case +/-32 degrees, and
            % the rounding happens silently so a scan in picoseconds cannot see it.
            %
            % Writing the phase into the cosine instead costs nothing: no extra samples, no
            % extra upload, resolution limited only by floating point. Same patch as
            % AWGExpESR already carries.
            %
            % SET timeDelay = 0 WHEN USING THIS. The two mechanisms are additive -- this line
            % adds to what channelsDifference returned -- so having both non-zero means the
            % real phase is their sum and nothing records that.
            phase = phase + obj.channelPhaseRadians();

            % ---- waveforms: exactly as the classic path builds them ----------
            initZeroIndex  = zeros(1, obj.nChannels);
            constStepIndex = zeros(1, obj.nChannels);
            for chan = 1:obj.nChannels
                fZero = @(t) 0 * t;
                initZeroIndex(chan)  = obj.AWG.AddWaveformByFunction(fZero, prePulseTime(chan), 0);
                fCos = @(t) A(chan) * cos(w*t - phase(chan));
                constStepIndex(chan) = obj.AWG.AddWaveformByFunction(fCos, baseWFtime, 0);
            end

            shortWaveformsIndex = zeros(obj.nChannels, 0);
            shortWaveformsTime  = [];
            waveformsIndex      = zeros(obj.nChannels, N);
            numZeroStepsInPoint = zeros(1, N);
            fprintf('Walk load: building waveforms for %d taus...\n', N);
            for i = 1:N
                numZeroStepsInPoint(i) = floor(obj.tau(i) / baseWFtime);
                residue = mod(obj.tau(i), baseWFtime);
                % min() OF AN EMPTY VECTOR RETURNS EMPTY, and && will not accept that:
                % "Operands to || and && must be convertible to logical scalar values".
                % The classic LoadAWG has the same min() and gets away with it only
                % because it feeds the result straight to `if closeWFdiff < ...`, and in
                % MATLAB `if []` is quietly false. That is an accident, not a decision --
                % it happens to be the behaviour wanted on the first tau, when no residue
                % waveform exists yet. Made explicit here so the guard is a real scalar.
                if isempty(shortWaveformsTime)
                    closeWFdiff = inf; closeWFindex = 1;
                else
                    [closeWFdiff, closeWFindex] = min(abs(residue - shortWaveformsTime));
                end
                ind = length(shortWaveformsTime) + 1;
                for chan = 1:obj.nChannels
                    if closeWFdiff < 0.01*obj.tau(i)
                        waveformsIndex(chan, i) = shortWaveformsIndex(chan, closeWFindex);
                        obj.tau(i) = numZeroStepsInPoint(i)*baseWFtime + shortWaveformsTime(closeWFindex);
                    else
                        shortWaveformsTime(ind) = residue;
                        fRes = @(t) A(chan) * (0<=t & t<residue) .* cos(w*t - phase(chan));
                        shortWaveformsIndex(chan, ind) = ...
                            obj.AWG.AddWaveformByFunction(fRes, residue + minWaveformTime, 0, 1);
                        waveformsIndex(chan, i) = shortWaveformsIndex(chan, ind);
                    end
                end
            end
            obj.AWG.LoadWaveform();

            % The classic path snaps taus onto existing residue waveforms and then deletes
            % duplicates. Same here, and the deletion has to take the parallel arrays with
            % it or the blocks stop matching the taus they are supposed to play.
            while any(diff(obj.tau) == 0)
                i = find(diff(obj.tau) == 0, 1);
                obj.tau(i) = [];
                waveformsIndex(:, i) = [];
                numZeroStepsInPoint(i) = [];
            end
            N = length(obj.tau);

            % ---- the block order, randomised per load ------------------------
            % The classic path draws a fresh randperm every average. Walk cannot: the
            % visit order IS the sequence order. So the randomisation moves into the
            % sequence, and the order is recorded so performWalk can put the data back in
            % tau order. Without this, slow drift over a pass would appear as a systematic
            % trend in tau -- which on a Rabi is indistinguishable from a decay envelope.
            % With newOrderEveryAverage, performWalk redraws it before every later average
            % by rebuilding only the sequence from walkParts.
            %
            % FIELD BY FIELD, not struct(...): struct() expands a cell value into a struct
            % array (see awgLoadFingerprint).
            p = struct();
            p.initZeroIndex       = initZeroIndex;
            p.constStepIndex      = constStepIndex;
            p.waveformsIndex      = waveformsIndex;
            p.numZeroStepsInPoint = numZeroStepsInPoint;
            p.nWarm               = obj.walkWarmupTriggers;
            obj.walkParts = p;
            nSteps = obj.uploadWalkSequence(randperm(N));
            obj.stampLoaded(fpAtEntry);
            fprintf(['Walk load: %d taus + %d warm-up block(s), %d steps per channel, ' ...
                     '%d AddSequence calls total (the classic path would make %d).\n'], ...
                    N, p.nWarm, nSteps, obj.nChannels, obj.nChannels*N);
        end

        function nSteps = uploadWalkSequence(obj, order)
            % Build ONE sequence per channel that plays the taus in 'order' (order(k) = the
            % tau index of the k-th block), upload it, assign it and arm. Only the sequence
            % list is replaced: the waveforms named in obj.walkParts must already be on the
            % instrument. LoadAWGWalk calls this for the first order and
            % newWalkOrderIfAsked for every later one, so both build identical blocks.
            %
            % Block layout, matching AWGExpESR's walk blocks:
            %   [initZero, constStep x numZeroSteps, residue]   trigger flags [1 0 0]
            %   [initZero, residue]                             trigger flags [1 0]
            p = obj.walkParts;
            N = numel(order);
            nWarm = p.nWarm;
            blockOrder = [repmat(order(1), 1, nWarm), order];
            stepsFor = @(i) 2 + (p.numZeroStepsInPoint(i) > 0);
            nSteps = sum(arrayfun(stepsFor, blockOrder));

            % Nothing may still reference the sequence that is about to be replaced, and
            % LoadSequence deletes every sequence on the instrument before it uploads.
            obj.AWG.sendCommand('awgcontrol:stop:immediate');
            obj.AWG.sequences = struct('waveforms', {{}}, 'repeats', {{}}, ...
                                       'trigger', {{}}, 'goto', {{}});   % sequences only
            obj.indexSeq = zeros(obj.nChannels, N);
            obj.seqAssetIndex = zeros(1, obj.nChannels);
            for chan = 1:obj.nChannels
                wfList = zeros(1, nSteps); rptList = zeros(1, nSteps);
                trgList = zeros(1, nSteps); gotoList = zeros(1, nSteps);
                s0 = 0;
                for blk = 1:numel(blockOrder)
                    i = blockOrder(blk);
                    if p.numZeroStepsInPoint(i) == 0
                        wfList(s0+1:s0+2)  = [p.initZeroIndex(chan), p.waveformsIndex(chan,i)];
                        rptList(s0+1:s0+2) = [1, 1];
                        trgList(s0+1:s0+2) = [1, 0];
                        s0 = s0 + 2;
                    else
                        wfList(s0+1:s0+3)  = [p.initZeroIndex(chan), p.constStepIndex(chan), ...
                                              p.waveformsIndex(chan,i)];
                        rptList(s0+1:s0+3) = [1, p.numZeroStepsInPoint(i), 1];
                        trgList(s0+1:s0+3) = [1, 0, 0];
                        s0 = s0 + 3;
                    end
                    if blk > nWarm
                        obj.indexSeq(chan, i) = s0;   % last step of tau i's block
                    end
                end
                % The final step sends the sequencer back to the first step, so a pass
                % ends where the next one begins. Without it the sequencer stops at the
                % end and pass two would need a re-arm.
                gotoList(nSteps) = 1;
                obj.seqAssetIndex(chan) = obj.AWG.AddSequence(wfList, rptList, trgList, 0, gotoList);
            end
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
            % The load cache is re-stamped under the SAME key with the new order, so a later
            % run that hits the cache maps its data back with the order the AWG really plays.
            if obj.newOrderEveryAverage && ~obj.walkOrderFresh
                if isempty(obj.walkParts)
                    obj.sendWarning(['newOrderEveryAverage: the loaded sequence predates this ' ...
                        'option, so it cannot be rebuilt on its own. Keeping its order for this ' ...
                        'average; change any setting to force a reload.']);
                else
                    tNew = tic;
                    obj.uploadWalkSequence(randperm(length(obj.tau)));
                    obj.stampLoaded(obj.walkCacheKey);
                    fprintf('Average %d: new tau order -- sequence re-upload + arm %.1f s.\n', ...
                            obj.currIter, toc(tNew));
                end
            end
            obj.walkOrderFresh = false;
        end

        function performWalk(obj)
            % One average, acquired as whole passes over the tau list.
            %
            % The sequencer advances one tau per trigger, so a pg.run of
            % (N + warmup) * passesPerRun triggers returns
            %     [s1 r1 s2 r2 ... sN rN]  repeated passesPerRun times
            % in BLOCK order, which is obj.walkOrder, not tau order. Nothing is sent to
            % the AWG inside the loop and nothing is sent to the pulse generator either --
            % the gate is fixed, so there is no changeSequence per point.
            %
            % Note what `repeats` means here: passes, exactly as in AWGExpESR's walk. Each
            % tau's repeats arrive spread out in time rather than back to back, which is
            % strictly better for drift and is why walkReduce exists to do the same
            % reduction processData would.
            N = length(obj.tau);
            pg   = getObjByName(PulseGenerator.NAME);
            spcm = getObjByName(Spcm.NAME);
            tracker = getObjByName(Tracker.NAME);
            if isempty(tracker); throwBaseObjException(Tracker.Name); end

            % A new visit order for this average, if asked (newOrderEveryAverage).
            obj.newWalkOrderIfAsked();

            passesWanted = obj.repeats;
            P = obj.passesPerRunUsed;
            if isempty(P); P = obj.choosePassesPerRun(false); obj.passesPerRunUsed = P; end
            nRuns  = ceil(passesWanted / P);
            nWarm  = obj.walkWarmupTriggers;
            blocksPerPass = N + nWarm;
            readsPerRun = obj.detectionPeriodsPerRepeat * blocksPerPass * P;

            raw = nan(obj.detectionPeriodsPerRepeat, N, nRuns * P);
            passesDone = 0; failedPasses = 0;
            opts = obj.walkReduceOptionsRabi(spcm);

            for r = 1:nRuns
                if obj.checkEmergencyStop(); return; end
                success = false;
                for trial = 1:5
                    try
                        if obj.resyncEveryPass
                            % Put the sequencer back on step 1 so this pass is known to
                            % start on the first block. A stop-then-run is much cheaper
                            % than an asset change, and it bounds any trigger miscount to
                            % a single pass instead of corrupting the tau mapping for the
                            % rest of the run.
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
                                           obj.detectionPeriodsPerRepeat, blocksPerPass, P);
                        passData = passData(:, nWarm+1:end, :);     % drop the warm-up blocks
                        % Back into tau order. passData column k holds the tau that block
                        % k played, which is obj.walkOrder(k).
                        raw(:, obj.walkOrder, passesDone + (1:P)) = passData;
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

                % Reduce what has arrived so far, so the plot fills in as passes
                % accumulate rather than only at the end of the average.
                usable = min(passesDone, passesWanted);
                [sig, err] = walkReduce(raw(:, :, 1:usable), ...
                                        obj.detectionPeriodsPerRepeat, N, opts);
                obj.signal(:, :, obj.currIter) = sig;
                obj.sterr(:, :, obj.currIter)  = err;

                S1 = squeeze(obj.signal(1, :, 1:obj.currIter));
                S1e = squeeze(obj.sterr(1, :, 1:obj.currIter));
                S2 = squeeze(obj.signal(2, :, 1:obj.currIter));
                S2e = squeeze(obj.sterr(2, :, 1:obj.currIter));
                [value, sterrOut] = getRatioDistributionValues(obj, S1, S2, S1e, S2e);
                obj.signalParam.value = value;
                obj.signalParam.sterr = sterrOut;

                obj.currParamIter = min(N, round(N * usable / passesWanted));
                if r == 1 || mod(r, 5) == 0 || r == nRuns
                    sendEventParamIterationDone(obj);
                end

                % Tracking once per pass rather than once per point: there is no
                % per-point pause left to hang it on, and the whole-sweep mean
                % reference is a better brightness metric than one tau's.
                if obj.isTracking
                    if obj.checkEmergencyStop(); return; end
                    refLevel = mean(sig(2, :));
                    refError = mean(err(2, :));
                    if tracker.compareReference(refLevel, refError, ...
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

        function opts = walkReduceOptionsRabi(obj, spcm)
            % Everything walkReduce needs to reproduce processData exactly. Same contents
            % as AWGExpESR's version; duplicated rather than shared because these are two
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

        function LoadAWG(obj)
            % --- skip the whole upload if the AWG already holds this exact set ---
            % fpAtEntry is captured BEFORE anything mutates obj.tau, and it is the key
            % both branches stamp with.
            [hit, fpAtEntry] = obj.tryLoadFromCache();
            if hit; return; end
            obj.walkParts = [];     % set again by LoadAWGWalk; a classic load has none

            if obj.useSequenceWalk
                obj.LoadAWGWalk(fpAtEntry);
                return
            end
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
            % NO TOTAL, and the '~' in the old header was an admission of the same thing:
            % Rabi's waveform count CANNOT be known before the loop runs. Each tau is split
            % into whole base steps plus a residue, and a residue within 1% of one already
            % built is SNAPPED onto it instead of adding a waveform. How many get snapped
            % depends on the tau list, so the total is discovered, not predicted.
            %
            % progressBegin with total = 0 reports a live count and rate without pretending to
            % know the end, and progressEnd then says how many were actually built -- which,
            % next to the number of taus, is the reuse figure and the thing worth knowing here.
            obj.AWG.progressBegin(sprintf(['Rabi waveforms: %d channel(s) x (1 zero pad + 1 ' ...
                'full base step) + one SHORT waveform per distinct tau residue. Residues ' ...
                'within 1%% of an existing one are reused, so the count below is discovered ' ...
                'rather than predicted -- %d taus is the upper bound.'], ...
                obj.nChannels, numel(obj.tau)), 0, 'waveform');
            
            % the resolution of the time delay is the sample rate, so we add also aditional phase delay
            [prePulseTime, phase] = obj.channelsDifference;
            % channelPhase, in RADIANS, added to whatever channelsDifference produced.
            %
            % channelsDifference lives in ExperimentAWG (shared, not ours to touch) and the
            % `phase` it returns is only the SUB-SAMPLE REMAINDER of timeDelay after the AWG
            % class rounds that delay to whole samples. It is not a phase knob and cannot be
            % used as one: at 16 GS/s one sample is 62.5 ps and one period at 2870 MHz is
            % 348.4 ps, so the settable delays land about 64.6 degrees apart. Roughly 0, 65,
            % 129, 194, 258, 323 -- about 5.6 points per period, best case +/-32 degrees, and
            % the rounding happens silently so a scan in picoseconds cannot see it.
            %
            % Writing the phase into the cosine instead costs nothing: no extra samples, no
            % extra upload, resolution limited only by floating point. Same patch as
            % AWGExpESR already carries.
            %
            % SET timeDelay = 0 WHEN USING THIS. The two mechanisms are additive -- this line
            % adds to what channelsDifference returned -- so having both non-zero means the
            % real phase is their sum and nothing records that.
            phase = phase + obj.channelPhaseRadians();
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
            nTauBefore = numel(obj.tau);
            obj.AWG.progressEnd();
            % The reuse figure, per channel, because that is the number that explains the
            % waveform count. A dense tau list on a coarse base step reuses almost everything;
            % a sparse one reuses nothing and the upload is nTau waveforms long.
            % round(): the division is exact in practice -- every residue is added once per
            % channel -- but a fractional count printed with %d would silently truncate.
            nDistinct = round(max((obj.AWG.prog.buildCount - 2*obj.nChannels), 0) / max(obj.nChannels,1));
            fprintf(['  %d distinct tau residues for %d taus, so %d tau(s) reused a waveform ' ...
                     'already built (%.0f%% reuse).\n'], ...
                    nDistinct, nTauBefore, max(nTauBefore - nDistinct, 0), ...
                    100*max(nTauBefore - nDistinct, 0)/max(nTauBefore, 1));
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
            % This loop printed NOTHING before, so a long Rabi looked hung between the
            % waveform upload and the sequence upload.
            obj.AWG.progressBegin(sprintf(['Rabi sequences: one per tau per channel, 2 or 3 ' ...
                'steps each (zero pad, base step repeated floor(tau/base) times, residue).']), ...
                obj.nChannels*N, 'sequence');
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
                    obj.AWG.progressStep();
                end    
            end
            obj.AWG.progressEnd();
            obj.AWG.LoadSequence()
            obj.stampLoaded(fpAtEntry);
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
                        % AWG.Run ALREADY BLOCKS UNTIL THE INSTRUMENT IS ARMED, so the
                        % pause(1) that used to sit here was pure waiting.
                        %
                        % Run() calls deviceReady() (an *OPC? barrier), then
                        % 'awgcontrol:run:immediate', then waitForRunning(), which polls
                        % AWGControl:RSTate? until the state is no longer Stopped. When it
                        % returns, the AWG reports armed or running. There is nothing left
                        % for a fixed delay to wait for.
                        %
                        % What it cost: one second per point per average. At 97 taus and 3
                        % averages that is 291 s -- about five minutes of a run in which the
                        % photons themselves take roughly 75 ms per point. Removing it does
                        % not change what is measured; the trigger from the pulse generator
                        % is what starts the burst, and that cannot arrive before Run has
                        % returned.
                        %
                        % NOTE the remaining per-point cost is NOT removable this way.
                        % assignSequence STOPS playback on this instrument (measured: ~500 ms
                        % to re-arm, regardless of how much is loaded), which is exactly why
                        % Run has to be called again every point. The fix for that is the
                        % single-sequence scheme used by AWGExpESR's jump mode: build ONE
                        % sequence whose steps are the taus and use jumpToStep, which does
                        % not stop playback and costs one SCPI write (~3 ms). That is a
                        % larger change and is not attempted here.
                        obj.AWG.Run;

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
