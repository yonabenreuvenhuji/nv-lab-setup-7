classdef FrequencyGeneratorTektronixAWG < FrequencyGenerator
        %FREQUENCYGENERATOR Tektronix AWG class
        %
        % ---------------------------------------------------------------------
        % PERFORMANCE CHANGES (see also obj.commDelayTime / obj.skipRedundantDeletes)
        %
        % 1. Redundant deletes removed. LoadWaveform/LoadSequence already wipe
        %    the whole waveform/sequence list with a single 'delete ALL', so the
        %    per-waveform and per-sequence delete inside WriteWaveformToAWG /
        %    WriteSequenceToAWG were deleting something that had just been
        %    wiped. Skipping them saves one SCPI command per waveform AND per
        %    sequence (~240 commands for a 120-point ESR sweep).
        %    Controlled by obj.skipRedundantDeletes (default true).
        %
        % 2. The unconditional inter-command delay is now tunable. It used to be
        %    a hard-coded 5 ms busy-wait after EVERY SCPI write, which for a
        %    120-point ESR sweep is ~2200 commands => ~11 s of pure spinning.
        %    Set obj.commDelayTime = 0 to remove it entirely.
        %    Status-polling loops now use their own POLL_DELAY so that setting
        %    commDelayTime = 0 does not turn them into tight spin loops.
        %
        %    obj.commandCount / obj.resetCommandCount are provided so you can
        %    measure the effect instead of guessing.
        % ---------------------------------------------------------------------


    properties (Constant, Hidden)
        TYPE = 'TektronixAWG';
        MODE_FG = 'FGEN';               % function generator mode
        MODE_AWG = 'AWG';               % AWG mode
        
        NEEDED_FIELDS = {'address', 'serialNumber', 'minFrequency', 'maxFrequency', 'minAmplitude', 'maxAmplitude', ...
                         'trigger1Channel', 'trigger1ChannelName',  'onDelay1', 'offDelay1'};
        OPTIONAL_FIELDS = {'trigger2Channel', 'trigger2ChannelName',  'onDelay2', 'offDelay2'};
        MAX_SAMPLE_RATE = 16e3;         % MHz
        MIN_SAMPLE_RATE = 8e3;          % MHz. recommended, not physical
        MAX_VPP = 0.50;                 % The constant Vpp we work with
        TRIG_SOURCE = 'EXT';            % 'INT' for internal trigger. 'EXT' for external trig. Use internal for testing.
        TIME_CONVERT = 1e6;
        MAX_REAPETS = 10e5;
        MIN_WAVEFORM_LENGTH = 2400;     % minimum samples for waveform
        MAX_SEQUENCE_STEPS = 16384;     % max steps in a single sequence (AWG70000 family)
        MAX_BATCH_COMMANDS = 40;        % max SCPI commands joined into one compound message.
                                        %   Tektronix compound messages need ';:' between commands --
                                        %   a bare ';' keeps the parser in the current subtree, which
                                        %   is exactly what made 'AWGC:STOP; OUTP1 OFF' fail. Kept
                                        %   modest so no single message gets unreasonably long.
        VISA_BRAND = 'ni';              % VISA brand
        BUFFER = 5e7;                   % buffer length
        TRIG_DURATION = 0.1;
        COMM_DELAY = 0.005;             % double. (Default value) time (in seconds) between consecutive commands.
                                        %   Kept for backwards compatibility with any external reference;
                                        %   the value actually used at runtime is obj.commDelayTime.
        POLL_DELAY = 0.001;             % double. Interval (s) between consecutive polls of a status query.
                                        %   Independent of commDelayTime, so zeroing that does not create a tight spin.
        TRANSFER_TIMEOUT = 5;           % double. Max time (s) to wait for a VISA write to report 'idle'
        NOMINAL_TIMEOUT = 100;          % double. The ambient VISA read timeout, asserted on
                                        %   connect. visadev defaults to 10 s and nothing
                                        %   raised it, which is why a legitimate wait for the
                                        %   instrument to drain showed up as "timed out after
                                        %   10 s". Also the fallback for ioGetTimeout: a
                                        %   save/restore pair that cannot read the current
                                        %   value must NOT restore a 2 s probe timeout, or
                                        %   every such pair silently downgrades the ambient.

        OPC_TIMEOUT = 2;                % double. Max time (s) to wait for *OPC? to confirm previous command finished
        LOAD_SETTLE_TIMEOUT = 600;      % double. Max time (s) to wait for the instrument to finish
                                        %   digesting a waveform/sequence upload. Needs to be far
                                        %   larger than OPC_TIMEOUT: measured 15.7 s of *OPC? wait
                                        %   after a 24.6 MB upload on the AWG70002A, and building a
                                        %   single 288-step concatenated sequence takes longer still.
                                        %   deviceReady raises the VISA read timeout to match, so this
                                        %   is the real ceiling.
        RUN_TIMEOUT = 120;              % double. Max time (s) to wait for AWGControl:RSTate? to report "Running".
                                        %   MEASURED on the AWG70002A: the first arm after loading a
                                        %   sequence scales LINEARLY with its step count, about 0.16 s
                                        %   per step --
                                        %        36 steps  5.8 s
                                        %        72 steps 11.1 s
                                        %       144 steps 23.4 s
                                        %       288 steps ~46 s (extrapolated)
                                        %   Subsequent arms of an already-assigned sequence are fast
                                        %   (<2 s), so this only bites on the first arm after a load.
                                        %   Earlier values of 2, 5 and 30 s each produced a spurious
                                        %   "device does not enter into RUN mode" that looked like the
                                        %   instrument refusing a long sequence when it was simply
                                        %   still working. Do not lower this without re-measuring.
        RSTATE_STOPPED = 0;             % AWGControl:RSTate? return value meaning "Stopped" (0=Stopped, 1=Waiting for trigger, 2=Running)
        RUN_MODE_CONT = 'continuous';
        RUM_MODE_TRIG = 'triggered';
        RUN_MODE_TCON = 'tcontinuous';
    end
       
    properties (Access = private)
        address
        v                       % visa object
        mode                    % AWG / FGEN (function generator)
    end
    
    properties
        MaxNormPowerAllowed
        Vpp
        workChannels

        SampleRate              % in MHz
        waveforms               % segmentLibrary{channel}, where channel is 1 or 2. Each segment must by of length 1(DC), of length >=16 and devides by 4.
        sequences
        prog = struct('active',false, 'label','', 'unit','item', 'total',0, 'n',0, ...
                      'samples',0, 'maxLen',0, 'lastLen',0, 'lastAt',0, 'mile',0, ...
                      'buildSeconds',NaN, 'buildSamples',0, 'buildCount',0, 'tic',[]);
                                % struct. State of the current progress report. buildSeconds /
                                %   buildSamples / buildCount survive progressEnd so LoadWaveform
                                %   can separate the time spent BUILDING samples in MATLAB from
                                %   the time spent SENDING them, which are the two things you
                                %   would act on differently.

        commDelayTime = 0.005;  % double. Delay (s) inserted after every SCPI write.
                                %   Historically hard-coded at 0.005 (see COMM_DELAY). This is a
                                %   BUSY-WAIT, not pause(), and it is paid once per command --
                                %   a 120-point ESR load is ~2200 commands, i.e. ~11 s of spinning.
                                %   Set to 0 to remove it. If the instrument starts dropping or
                                %   mis-ordering commands, raise it back towards 0.005 rather than
                                %   leaving it somewhere arbitrary.
                                %   NOTE: deliberately not pause(). On Windows, pause() granularity
                                %   is far coarser than a few ms, so pause(0.005) can cost ~15 ms
                                %   and would make things worse.

        useVisadev = false;     % logical. Which MATLAB instrument interface to talk through.
                                %
                                % false - the legacy visa()/fopen/fwrite/fscanf/query interface.
                                %         Deprecated since ~R2021a. Every property access on the
                                %         object goes through instrument.subsref ->
                                %         ICTRemoveFunctionalityHelper.showWarning, which builds a
                                %         TABLE just to decide whether to warn. Profiled on a
                                %         191-point ESR: 11,532 calls, 10.9 s; plus 5764 connect()
                                %         status checks for 5.7 s and 5.8 s polling TransferStatus.
                                %
                                % true  - visadev()/write/writeline/writeread (R2021a onward, so
                                %         fine on R2023b). Writes are synchronous, so
                                %         verifyOperationComplete becomes a no-op, and the binary
                                %         waveform transfer becomes one message with no EOIMode
                                %         juggling.
                                %
                                % Expect ~10-15 s back on a 96-point sweep. Will NOT change the
                                % ~450 ms arm latency -- that is instrument firmware, and our poll
                                % granularity (~30 ms) is well inside it.
                                %
                                % Switching this at runtime closes the connection, reopens it on the
                                % other interface, and resets the instrument -- so any loaded
                                % waveforms are gone and the load cache is invalidated. Set it once,
                                % in the speed-knobs cell, right after main.

        armWithoutChecking = false;  % logical. In EnsureRunning, skip the RSTate? query and
                                %   just re-arm unconditionally.
                                %
                                %   Measured on the AWG70002A: the RSTate? query costs ~1.2 s while
                                %   the sequencer is armed (and was once seen at 21.6 s), whereas
                                %   Run() -- which polls the SAME query internally -- costs ~0.5 s.
                                %   On the classic path assignSequence has already stopped the AWG
                                %   at 88-100% of points, so asking first is strictly more expensive
                                %   than just arming:
                                %
                                %       check-then-arm : 1.2 + 0.88*0.5 = 1.64 s per point
                                %       arm always     :             0.5 = 0.50 s per point
                                %
                                %   Safe because EnsureRunning is only ever called before pg.run,
                                %   so nothing is mid-playback when the re-arm happens.
                                %
                                %   Default false so no other experiment's behaviour changes by
                                %   surprise. Turn it on in your speed-knobs cell.

        jumpSettleTime = 0.010; % double. Delay (s) after a forced jump, BEFORE the pulse
                                %   generator is allowed to trigger. This is load-bearing for
                                %   data correctness, not a performance knob.
                                %
                                %   sendCommand uses fwrite, which returns as soon as the bytes
                                %   leave the PC -- not when the instrument has acted on them. So
                                %   without a wait here, pg.run can fire the trigger while the
                                %   sequencer is still parked on the PREVIOUS frequency's block,
                                %   and the point gets measured at the wrong frequency.
                                %
                                %   Measured signature when this was 0: jump mode collected the
                                %   right amount of contrast (2.16% vs 2.10% on the classic path)
                                %   but its correlation with the classic spectrum was +0.04, and
                                %   point-to-point scatter was 5.5x the statistical noise --
                                %   real data, filed under the wrong frequencies. Because the
                                %   sweep order is randperm, a one-point lag scrambles the
                                %   spectrum rather than shifting it.
                                %
                                %   Earlier jump runs looked fine only because EnsureRunning's
                                %   ~1.2 s RSTate? query happened to sit here and act as the
                                %   delay. Calibrate the real minimum with esrJumpValidate.

        skipRedundantDeletes = true;    % logical. Skip the per-waveform / per-sequence 'delete'
                                        %   inside WriteWaveformToAWG / WriteSequenceToAWG, which is
                                        %   redundant after the bulk 'delete ALL' that LoadWaveform /
                                        %   LoadSequence already issue. Saves ~1 command per waveform
                                        %   and per sequence. Set false to restore the old behaviour.
    end

    properties (SetAccess = private)
        loadedFingerprint = '';     % char. Hash of the parameters the currently loaded
                                    %   waveforms/sequences were built from, so an identical
                                    %   sweep can skip the whole upload. Cleared by Reset and
                                    %   ResetStoredData, so any other experiment loading the AWG
                                    %   invalidates it automatically.
        loadedCache = [];           % Whatever the experiment needs to reconstruct its own
                                    %   bookkeeping on a cache hit (indexSeq etc.). Opaque here.
    end

    properties (Hidden)
        commandCount = 0;       % double. Number of SCPI commands sent since the last resetCommandCount.
                                %   Instrumentation only -- use it to verify that a change actually
                                %   reduced the traffic:
                                %       awg = getObjByName('TektronixAWG'); awg.resetCommandCount;
                                %       tic; exp.run; toc; awg.commandCount

        progressStyle = 'inplace';   % char. How the waveform/sequence build reports itself.
                                %   'inplace'   one self-overwriting line with count, rate and
                                %               ETA. The default.
                                %   'milestone' a line at 25/50/75/100%, for a diary or when
                                %               something else is printing at the same time and
                                %               the backspaces of 'inplace' would corrupt it.
                                %   'count'     the legacy '1 2 3 4 ...' per index.
                                %   'quiet'     nothing until the summary.
                                %
                                %   Why this exists: a 167-point Ramsey builds 668 waveforms, and
                                %   the legacy style prints 668 integers -- four screens of digits
                                %   that say nothing about WHAT is being built or how long it will
                                %   take, and that push the actual result out of the scrollback.

        timing = struct('assign',0, 'assignCalls',0, ...
                        'run',0,    'runCalls',0, ...
                        'opcWait',0,'armWait',0, ...
                        'ensureCalls',0, 'ensureStopped',0, ...
                        'stateQuery',0);
                                % struct. Wall-clock spent in the per-point hot path, accumulated
                                %   during a real sweep so there is no probing artefact. Read it with
                                %   printInstrumentation. Note 'run' CONTAINS 'opcWait' + 'armWait'.
                                %
                                %   ensureStopped/ensureCalls is the key ratio: it is the fraction of
                                %   sweep points at which the AWG had dropped to Stopped and needed a
                                %   full re-arm. If that is near 1, assigning a sequence stops the
                                %   instrument and the per-point re-assign has to go. If it is near
                                %   0, the AWG stays armed during a real sweep and the overhead is
                                %   somewhere else.
    end

    properties (Dependent)
        minWaveformDuration     % in mus
    end

    %% Performance knobs and instrumentation
    methods
        function set.commDelayTime(obj, newVal)
            if ~isnumeric(newVal) || ~isscalar(newVal) || ~isreal(newVal) || newVal < 0 || newVal > 1
                error('commDelayTime must be a real scalar between 0 and 1 (seconds), the current is %s', num2str(newVal));
            end
            obj.commDelayTime = newVal;
        end

        function set.useVisadev(obj, newVal)
            if ~(islogical(newVal) || (isnumeric(newVal) && isscalar(newVal) && ismember(newVal, [0 1])))
                error('useVisadev must be true or false');
            end
            newVal = logical(newVal);
            if newVal == obj.useVisadev
                return
            end
            wasOpen = obj.ioIsOpen();

            % DISCARD the handle, do not merely close it. The two interfaces use
            % different object classes, so ioOpen has to be free to construct the
            % right one -- and ioOpen only constructs when obj.v is empty or
            % invalid. fclose() leaves a legacy object non-empty AND valid, so
            % ioOpen would skip construction and we would end up driving a legacy
            % object with visadev calls. That was the original hang.
            obj.ioDiscard();
            obj.useVisadev = newVal;
            obj.invalidateLoaded();

            if wasOpen
                obj.connect();
                obj.ioSetTimeoutRaw(100);
                % Re-apply the state this class caches. Deliberately NOT *RST:
                % reopening a VISA session does not reset the instrument, and
                % *RST on an AWG70002A is slow enough to look like a hang.
                if ~isempty(obj.SampleRate)
                    obj.SampleRate = obj.SampleRate;      %#ok<ASGSL> re-sends clock:srate
                end
                if ~isempty(obj.workChannels)
                    obj.workChannels = obj.workChannels;  %#ok<ASGSL> re-sends output state
                end
            end
            fprintf('AWG transport switched to %s.\n', ...
                    BooleanHelper.ifTrueElse(newVal, 'visadev', 'legacy visa'));
        end

        function set.armWithoutChecking(obj, newVal)
            if ~(islogical(newVal) || (isnumeric(newVal) && isscalar(newVal) && ismember(newVal, [0 1])))
                error('armWithoutChecking must be true or false');
            end
            obj.armWithoutChecking = logical(newVal);
        end

        function set.jumpSettleTime(obj, newVal)
            if ~isnumeric(newVal) || ~isscalar(newVal) || ~isreal(newVal) || newVal < 0 || newVal > 1
                error('jumpSettleTime must be a real scalar between 0 and 1 (seconds), the current is %s', num2str(newVal));
            end
            obj.jumpSettleTime = newVal;
        end

        function set.skipRedundantDeletes(obj, newVal)
            if ~(islogical(newVal) || (isnumeric(newVal) && isscalar(newVal) && ismember(newVal, [0 1])))
                error('skipRedundantDeletes must be true or false');
            end
            obj.skipRedundantDeletes = logical(newVal);
        end

        function resource = getResource(obj)
            % The VISA resource string this object is talking through.
            resource = obj.address;
        end

        function setResource(obj, newAddress)
            % Point this object at a DIFFERENT VISA resource for the SAME
            % instrument -- e.g. swapping VXI-11 (…inst0::INSTR) for USBTMC or a
            % raw socket. visadevlist shows what is available.
            %
            % Per-transaction latency differs a lot between these, and this class
            % polls AWGControl:RSTate? repeatedly while arming, so the transport
            % can dominate what looks like instrument latency. Use
            % awgTransportBench to measure before committing.
            %
            % Loaded waveforms live in the instrument, not the session, so they
            % survive this. The load cache is invalidated anyway, out of caution.
            if ~(ischar(newAddress) || isstring(newAddress))
                error('setResource expects a VISA resource string');
            end
            newAddress = char(newAddress);
            wasOpen = obj.ioIsOpen();
            obj.ioDiscard();
            obj.address = newAddress;
            obj.invalidateLoaded();
            if wasOpen
                obj.connect();
                obj.ioSetTimeoutRaw(100);
                if ~isempty(obj.SampleRate)
                    obj.SampleRate = obj.SampleRate;      %#ok<ASGSL>
                end
                if ~isempty(obj.workChannels)
                    obj.workChannels = obj.workChannels;  %#ok<ASGSL>
                end
            end
            fprintf('AWG resource is now %s\n', newAddress);
        end

        function response = queryRaw(obj, command, timeoutSeconds)
            % Send a raw SCPI query and return the response, for diagnostics.
            %
            %   awg.queryRaw('SOURce1:RMODe?')
            %   awg.queryRaw('SYSTem:ERRor?')
            %
            % Use a short timeoutSeconds when probing a mnemonic you are not
            % sure the instrument supports -- the class default is 100 s, and an
            % unrecognised query will otherwise hang MATLAB for that long.
            % Restores the previous timeout afterwards, including on error.
            obj.connect();
            oldTimeout = obj.ioGetTimeout();
            if exist('timeoutSeconds', 'var') && ~isempty(timeoutSeconds)
                obj.ioSetTimeoutRaw(timeoutSeconds);
            end
            try
                response = strtrim(obj.ioQuery(command));
            catch err
                obj.ioSetTimeoutRaw(oldTimeout);
                rethrow(err);
            end
            obj.ioSetTimeoutRaw(oldTimeout);
        end

        function tf = isAlreadyLoaded(obj, fingerprint)
            % True when the AWG already holds exactly this waveform/sequence set.
            % Lets an unchanged sweep skip the upload entirely -- measured at
            % 17.8 s on a 191-point ESR, paid again on every single re-run.
            tf = ~isempty(obj.loadedFingerprint) && ...
                 ~isempty(fingerprint) && ...
                 strcmp(obj.loadedFingerprint, fingerprint);
        end

        function cache = getLoadedCache(obj)
            cache = obj.loadedCache;
        end

        function markLoaded(obj, fingerprint, cache)
            % Call ONLY after a load has completed successfully.
            obj.loadedFingerprint = fingerprint;
            if exist('cache', 'var')
                obj.loadedCache = cache;
            else
                obj.loadedCache = [];
            end
        end

        function invalidateLoaded(obj)
            obj.loadedFingerprint = '';
            obj.loadedCache = [];
        end

        function resetCommandCount(obj)
            % Zero the SCPI command counter and the hot-path timers.
            % Call immediately before a timed operation.
            obj.commandCount = 0;
            obj.timing = struct('assign',0, 'assignCalls',0, ...
                                'run',0,    'runCalls',0, ...
                                'opcWait',0,'armWait',0, ...
                                'ensureCalls',0, 'ensureStopped',0, ...
                                'stateQuery',0, 'wallStart',tic);
            % wallStart lets printInstrumentation report the true elapsed time.
            % Do NOT measure a run with a bare "tic; a.run; toc" -- this class
            % calls tic internally many times, and a bare toc measures from the
            % most recent one, so it reports a small meaningless number. Use
            % t0 = tic; a.run; toc(t0), or just read the line below.
        end

        function printInstrumentation(obj)
            % Report where the per-point time actually went during the last run.
            %
            %   awg.resetCommandCount; a.run; awg.printInstrumentation
            %
            % Unlike awgSpeedCheck's probe loop, these numbers are accumulated
            % DURING the sweep, while the pulse generator is actually firing, so
            % they reflect the real behaviour rather than an idle instrument.
            t = obj.timing;
            fprintf('\n=== AWG hot-path timing ===\n');
            if isfield(t, 'wallStart') && ~isempty(t.wallStart)
                fprintf('  elapsed since reset     : %7.2f s   <-- the real wall clock\n', ...
                        toc(t.wallStart));
            end
            fprintf('  SCPI commands sent      : %d\n', obj.commandCount);
            fprintf('  assignSequence          : %7.2f s over %d calls (%.2f ms each)\n', ...
                    t.assign, t.assignCalls, 1e3*t.assign/max(t.assignCalls,1));
            fprintf('  EnsureRunning           : %d calls, %d found the AWG STOPPED (%.0f%%)\n', ...
                    t.ensureCalls, t.ensureStopped, 100*t.ensureStopped/max(t.ensureCalls,1));
            fprintf('    RSTate? state query   : %7.2f s over %d calls (%.0f ms each)\n', ...
                    t.stateQuery, t.ensureCalls, 1e3*t.stateQuery/max(t.ensureCalls,1));
            fprintf('  Run (re-arm)            : %7.2f s over %d calls (%.0f ms each)\n', ...
                    t.run, t.runCalls, 1e3*t.run/max(t.runCalls,1));
            fprintf('    of which *OPC? wait   : %7.2f s\n', t.opcWait);
            fprintf('    of which RSTate? wait : %7.2f s\n', t.armWait);
            fprintf('  ---------------------------------------------\n');
            fprintf('  total in AWG hot path   : %7.2f s\n', t.assign + t.run + t.stateQuery);
            if t.stateQuery > 0.5*(t.assign + t.run + t.stateQuery)
                fprintf('  >> Most of the cost is the RSTate? poll, not the sweep. If the AWG\n');
                fprintf('     never actually stops, raise AWGExpESR.armCheckEvery so the sweep\n');
                fprintf('     stops asking on every point.\n');
            end
            if t.ensureCalls > 0 && t.ensureStopped/t.ensureCalls > 0.5
                fprintf('  >> The AWG is dropping to Stopped on most points, so every point\n');
                fprintf('     pays a full re-arm. Raising repeats cannot fix this -- the\n');
                fprintf('     per-point sequence re-assign itself has to go.\n');
            end
            fprintf('\n');
        end
    end


    methods (Access = private)
        function obj = FrequencyGeneratorTektronixAWG(name, address, frequencyLimits, amplitudeLimits, triggersNames, channels, delays)
            % All models are the same in regards to controlling them, but
            % the limitations on the amplitude and on the allowed frequencies may vary.
            obj@FrequencyGenerator(name, frequencyLimits, amplitudeLimits);
            obj.address = address;
            obj.connect;
            if strcmp(obj.queryValue('mode'), obj.MODE_AWG)
                obj.mode = obj.MODE_AWG;
            else
                obj.mode = obj.MODE_FG;
            end
            obj.RegisterPgChannels(triggersNames, channels, delays);
            obj.Reset;
            obj.initialize;
        end
        
        
    end
    
%%  General functions     
    methods
        function connect(obj)
            obj.ioOpen();
            % ASSERT THE AMBIENT READ TIMEOUT. visadev objects default to Timeout = 10 s,
            % and nothing on this code path ever raised it -- setTimeOut's 30..1000
            % validation only ever applied to values a caller passed in, never to the
            % value the object was born with. So the real ambient timeout was MATLAB's
            % 10 s, not the 100 s queryRaw's own help text claims.
            %
            % That is what turned a legitimate wait into a failure: "INST:MODE? timed out
            % after 10 s" while the AWG was still draining a jump-mode backlog. The wait
            % was correct; the deadline was a default nobody chose.
            if obj.ioGetTimeout() < obj.NOMINAL_TIMEOUT
                obj.ioSetTimeoutRaw(obj.NOMINAL_TIMEOUT);
            end
        end

        function disconnect(obj)
            % Respect keepOn, which the rest of the framework already treats as
            % "leave this generator's output alone" -- Experiment.wrapUpInternal
            % checks it before setting fg.output = false, but this method used to
            % switch the outputs off unconditionally regardless.
            %
            % That inconsistency was invisible while the enableOutput command was
            % malformed and silently rejected. It is not invisible any more, so
            % honour the flag.
            if obj.ioIsOpen()
                if ~obj.keepOn
                    obj.setValue('enableOutput', '0', 1);
                    obj.setValue('enableOutput', '0', 2);
                end
                obj.ioClose();
            end
        end

        function delete(obj)
            obj.disconnect;
            if ~isempty(obj.v) && isvalid(obj.v)
                delete(obj.v)
            end
        end
        
        function sendCommand(obj, command, precision)
%             Actually sends command to hardware
            obj.connect() % Only happens if it is closed
            if nargin < 3
                try
                    obj.ioWriteText(command);
                catch err
                    obj.disconnect();
                    obj.connect();
                    obj.sendError(err);
                end
            else
                if strcmp(precision, 'single')
                    obj.ioWriteSingles(command);
                else
                    error('can not send command: "precision" does not right')
                end
            end
            obj.commandCount = obj.commandCount + 1;
            obj.commDelay();
            obj.verifyOperationComplete();
        end

        function sendCommandBatch(obj, commands)
            % Send a cell array of SCPI commands as compound messages, chunked
            % to MAX_BATCH_COMMANDS each.
            %
            % Every command in the list must be a SET command with no response.
            % Never mix a query in here -- the replies would queue up and
            % desynchronise every later read.
            %
            % The separator is ';:' not ';'. The leading colon resets the command
            % tree to root; a bare ';' resumes inside the previous command's
            % subtree, which silently produced -113 "Undefined header" errors
            % elsewhere in this class for years.
            if isempty(commands)
                return
            end
            if ~iscell(commands)
                commands = {commands};
            end
            for first = 1:obj.MAX_BATCH_COMMANDS:length(commands)
                last = min(first + obj.MAX_BATCH_COMMANDS - 1, length(commands));
                obj.sendCommand(strjoin(commands(first:last), ';:'));
            end
        end

        function value = readOutput(obj, what)
            % Get value returned from object
            value = obj.ioReadLine();
            switch what
                case {'frequency', 'freq', 'f'}
                    value = num2str(str2double(value)/1e6);  % convert Hz to MHz
                case {'enableoutput', 'output', 'enabled', 'enable'}
                    value = num2str(strcmp(value,'1;FGEN'));
            end
        end

        function setTimeOut (obj, TimeOut)
            if (TimeOut < 30) || (TimeOut > 1000)
                warning('timeout is not changed, it must be more than 30 and less than 1000')
            else
                obj.ioSetTimeout(TimeOut);
            end
        end

        function setAmplitude (obj, Power)
            % The actual amplitude is const, and we just change the 'MaxNormPowerAllowed'
            % just for back support:
            if nargin == 1
                obj.MaxNormPowerAllowed = 1;
            else
                if max(Power) > obj.maxAmpl
                    error('max value for power is %d dBm, the current is %d.', obj.maxAmpl, max(Power))
                end
                Power_Volt  = obj.dBmToVpp(Power);
                obj.MaxNormPowerAllowed = Power_Volt/obj.MAX_VPP;
            end
            for i = 1:length(Power)
                AmpCmd = sprintf('SOURCE%d:VOLTAGE:AMPLITUDE %d', i, obj.MAX_VPP);
                obj.sendCommand(AmpCmd);
            end
        end
        
        function set.SampleRate(obj, newSampleRate)
            if (newSampleRate < obj.MIN_SAMPLE_RATE) || (newSampleRate > obj.MAX_SAMPLE_RATE)
                error('Sample rate must be between %d and %d Ms, the current is %d', obj.MIN_SAMPLE_RATE, obj.MAX_SAMPLE_RATE, newSampleRate);
            end
            obj.SampleRate = newSampleRate;
            RateCommand = sprintf('clock:srate %f', obj.SampleRate*obj.TIME_CONVERT);
            obj.sendCommand(RateCommand);
        end
        
        function set.workChannels(obj, n)
            obj.workChannels = n;
            if ismember(1, n)
                obj.sendCommand('output1:state on');
            else
                obj.sendCommand('output1:state off');
            end
            if ismember(2, n)
                obj.sendCommand('output2:state on');
            else
                obj.sendCommand('output2:state off');
            end
        end
        
        function t = get.minWaveformDuration(obj)
            minTime = 1/obj.SampleRate * obj.MIN_WAVEFORM_LENGTH;
            t = ceil(minTime * obj.SampleRate) / obj.SampleRate;
        end
    end

%% Progress reporting
%
% Three calls: progressBegin(label, total, unit), progressStep(nSamples), progressEnd.
% Nothing here talks to the instrument, and the step call is throttled to ~8 updates a second,
% so a 668-waveform build spends microseconds on reporting rather than 668 fprintf round trips
% to the Command Window -- which on a slow console is itself a measurable part of the build.
    methods
        function progressBegin(obj, label, total, unit)
            % label: what is being built, in words. This is the part the old output was missing
            %        entirely -- '334' does not tell you it is one waveform per tau per channel
            %        for signal and reference.
            if ~exist('unit', 'var'); unit = 'waveform'; end
            obj.prog = struct('active',true, 'label',label, 'unit',unit, ...
                              'total',max(total,0), 'n',0, 'samples',0, 'maxLen',0, ...
                              'lastLen',0, 'lastAt',0, 'mile',0, ...
                              'buildSeconds',NaN, 'buildSamples',0, 'buildCount',0, 'tic',tic);
            if strcmp(obj.progressStyle, 'quiet'); return; end
            fprintf('%s\n', label);
            if total > 0 && ~strcmp(obj.progressStyle, 'count')
                fprintf('  %d %ss to build\n', total, unit);
            end
        end

        function progressStep(obj, nSamples)
            if ~obj.prog.active; return; end
            if ~exist('nSamples', 'var'); nSamples = 0; end
            p = obj.prog;
            p.n = p.n + 1;
            p.samples = p.samples + nSamples;
            p.maxLen  = max(p.maxLen, nSamples);
            obj.prog = p;

            switch obj.progressStyle
                case 'quiet'
                    return
                case 'count'
                    fprintf('%d ', p.n);
                    return
                case 'milestone'
                    if p.total <= 0; return; end
                    frac = p.n / p.total;
                    if frac >= (p.mile + 1)/4 - 1e-9 || p.n == p.total
                        obj.prog.mile = floor(frac*4);
                        el = toc(p.tic);
                        fprintf('  %3.0f%%  %d/%d  %s elapsed', 100*frac, p.n, p.total, ...
                                obj.humanTime(el));
                        if frac > 0 && frac < 1
                            fprintf(', ~%s left', obj.humanTime(el*(1-frac)/frac));
                        end
                        fprintf('\n');
                    end
                    return
            end

            % ---- 'inplace' ----
            % Throttled on WALL CLOCK, not on a count, because the whole point is to be cheap
            % and a modulo on the index still prints hundreds of times on a long build.
            el = toc(p.tic);
            if p.n < p.total && el - p.lastAt < 0.12; return; end
            obj.prog.lastAt = el;

            if p.total > 0
                frac = p.n / p.total;
                txt = sprintf('  %5.1f%%  %d/%d %ss', 100*frac, p.n, p.total, p.unit);
            else
                frac = 0;
                txt = sprintf('  %d %ss', p.n, p.unit);
            end
            if el > 0.4
                txt = [txt sprintf('  |  %.0f/s', p.n/el)];
                if frac > 0.02 && p.n < p.total
                    txt = [txt sprintf('  |  %s left', obj.humanTime(el*(1-frac)/frac))];
                end
            end
            % Erase exactly what was printed last time. \b is the only cursor motion that
            % behaves in the MATLAB Command Window -- \r does not.
            fprintf('%s%s', repmat(sprintf('\b'), 1, obj.prog.lastLen), txt);
            obj.prog.lastLen = numel(txt);
            if p.n >= p.total && p.total > 0; fprintf('\n'); obj.prog.lastLen = 0; end
        end

        function progressEnd(obj)
            if ~obj.prog.active; return; end
            p = obj.prog;
            el = toc(p.tic);
            if strcmp(obj.progressStyle, 'inplace') && p.lastLen > 0
                fprintf('%s', repmat(sprintf('\b'), 1, p.lastLen));
            elseif strcmp(obj.progressStyle, 'count')
                fprintf('\n');
            end
            if ~strcmp(obj.progressStyle, 'quiet')
                if p.samples > 0
                    fprintf('  built %d %ss, %s samples, %.1f MB, in %s  (%s samples/s)\n', ...
                            p.n, p.unit, obj.humanCount(p.samples), 4*p.samples/2^20, ...
                            obj.humanTime(el), obj.humanCount(p.samples/max(el,eps)));
                    fprintf('  longest %s %s samples = %.3f us\n', p.unit, ...
                            obj.humanCount(p.maxLen), p.maxLen/obj.SampleRate);
                else
                    fprintf('  built %d %ss in %s\n', p.n, p.unit, obj.humanTime(el));
                end
            end
            % Kept for LoadWaveform, which reports build against transfer.
            obj.prog.active       = false;
            obj.prog.buildSeconds = el;
            obj.prog.buildSamples = p.samples;
            obj.prog.buildCount   = p.n;
            obj.prog.lastLen      = 0;
        end
    end

    methods (Static)
        function s = humanTime(x)
            if ~isfinite(x);        s = '?';                            return; end
            if x < 1;               s = sprintf('%.0f ms', 1e3*x);      return; end
            if x < 90;              s = sprintf('%.1f s', x);           return; end
            if x < 5400;            s = sprintf('%.0f min %02.0f s', floor(x/60), mod(x,60)); return; end
            s = sprintf('%.1f h', x/3600);
        end

        function s = humanCount(x)
            if ~isfinite(x);        s = '?';                            return; end
            if x < 1e3;             s = sprintf('%.0f', x);             return; end
            if x < 1e6;             s = sprintf('%.1f k', x/1e3);       return; end
            if x < 1e9;             s = sprintf('%.2f M', x/1e6);       return; end
            s = sprintf('%.2f G', x/1e9);
        end
    end

%% Adding waveforms and sequences:
    methods

        function index = AddWaveformByPoints(obj, newWaveform, initPadding, finalPadding, indexPlot)
            % Stores a new waveform at the end of obj.waveforms. Use LoadWaveforms to load to AWG.
            % The waveform is given by an array of amplitudes in the range of -1 to 1.
            % The optional padding parameters if for zero points around the waveform. The value is in number op points.
            % this function returns the index of the new waveform.
            if ~exist('initPadding', 'var')
                initPadding = 0;
            end
            if ~exist('finalPadding', 'var')
                finalPadding = 0;
            end
            if ~exist('indexPlot', 'var')
                indexPlot = 1;
            end
            if length(newWaveform) + initPadding + finalPadding < obj.MIN_WAVEFORM_LENGTH
                error('The minimum length of waveform is %d, the current is %d', obj.MIN_WAVEFORM_LENGTH, length(newWaveform) + initPadding + finalPadding)
            end
            if max(abs(newWaveform)) > 1
                error('Amplitude waveform must be between -1 to 1');
            end
            s = obj.waveforms;
            s{end+1} = [zeros(1,initPadding), newWaveform, zeros(1,finalPadding)];
            obj.waveforms = s;
            index = length(obj.waveforms);
            if indexPlot
                if obj.prog.active
                    % The sample count goes with the step, so progressEnd can report MB and
                    % samples/s rather than just a count of things.
                    obj.progressStep(numel(s{end}));
                elseif ~strcmp(obj.progressStyle, 'quiet')
                    fprintf('%d ', index);      % legacy path, for callers that never begin a
                end                              % progress session. Unchanged on purpose.
            end
        end
        
        function [index,wf,t] = AddWaveformByFunction(obj,f, duration, shift, toRound, toReverse)
            % Adding waveform of f(t-shift), for t=[0,duration].
            % t in mus!!
            % the parameter 'toRound' is 1 if you want round the 'duration' to divided by the sample time
            % this function return the waveform and time vectors for debuging. 
            if ~exist('toRound', 'var')
                toRound = 0;
            end
            if ~exist('toReverse', 'var')
                toReverse = 0;
            end
            t_sample = 1/obj.SampleRate; %sign(duration) is for t left to zero, without using shift
            if toRound
                duration = round(duration/t_sample)*t_sample;
            end
            if sign(duration) == 1
                t = 0 : t_sample : duration-t_sample;
            else
                t = 0-t_sample : -t_sample : duration;
            end
            % tests:
            MAX_ERROR = 1e-9;
            if mod(duration,t_sample) > MAX_ERROR && t_sample - mod(duration,t_sample) > MAX_ERROR 
                error('The parameter ''duration'' not divided by the AWG sample time. Fix it or change the parameter ''toRound'' to ''True''')
            end
            if length(t) < obj.MIN_WAVEFORM_LENGTH
                error('Tha minimum waveform length the AWG can get is %d samples, The current is %d.', obj.MIN_WAVEFORM_LENGTH, length(t))
            end
            
            % creating the waveform:
            wf = f(t-shift);
            if toReverse
                wf = flip(wf);
            end
            index = AddWaveformByPoints(obj,wf);
        end

        function [index,WF,t] = AddWaveformByDuration(obj,newWaveformI,newWaveformQ,f,duration,t0)
            % Adds a new waveform with each point caried out for a duration given by "duration" - in \mus
            % newWaveformI/Q must be between -1:1. "f"- in \MHz

            if isempty(newWaveformI) && ~isempty(newWaveformQ)
                newWaveformI = zeros(size(newWaveformQ));
            elseif ~isempty(newWaveformI) && isempty(newWaveformQ)
                newWaveformQ = zeros(size(newWaveformI));
            end
            if sum(size(newWaveformI) ~= [length(f),length(duration)]) || sum(size(newWaveformQ) ~= [length(f),length(duration)])
                error('Inputs must have the right size')
            end
            
            % change from duration to points
            points = round(duration*obj.SampleRate);
            timeChange = abs(sum((points)/obj.SampleRate - duration));
            if timeChange > 1e-6
                warning('Change in time due to round ups detected when converting waveform duration to # of points. total change of %f',timeChange)
            end
            
            % rewrite waveform as points for the AWG, (with each point of period 1/obj.SampleRate)
            t = (0 : 1/(sum(points)) : 1-1/(sum(points))).*sum(duration) + t0;
            if length(t) < obj.MIN_WAVEFORM_LENGTH
                error('Tha minimum point at waveform the AWG can get is %d , The current is %d. (At maximum sumple rate is consider toal waveform duration of 0.15us)',obj.MIN_WAVEFORM_LENGTH,length(t))
            end
            waveformI = zeros(length(f),length(points));
            waveformQ = waveformI;
            startIndex = 1;
            for k = 1:length(points)
                endIndex = startIndex + points(k) - 1;
                for j=1:length(f)
                    waveformI(j,startIndex:endIndex) = newWaveformI(j,k);
                    waveformQ(j,startIndex:endIndex) = newWaveformQ(j,k);
                end
                startIndex = endIndex+1;
            end
            for j=1:length(f)
                waveformI(j,:) = waveformI(j,:) .* sin(2*pi*f(j)*(t));
                waveformQ(j,:) = waveformQ(j,:) .* cos(2*pi*f(j)*(t));
            end
            waveformI = sum(waveformI,1);
            waveformQ = sum(waveformQ,1);
            WF = waveformI + waveformQ;
            index = AddWaveformByPoints(obj,WF);
        end

        function [index,WF,t] = AddWaveformByDurationPhaseModulationSpinLock(obj,newWaveformI,newWaveformQ,FuncI,FuncQ,f,duration,t0,alpha,RF,spinLock)
            % Adds a new waveform, with each point caried out for a duration given by "duration" - in \mus
            % newWaveformI/Q must be between -1:1. "f"- in \MHz
            % 'spinLock' gets 1 for spin lock and 0 for spin unLock
            
            if isempty(newWaveformI) && ~isempty(newWaveformQ)
                newWaveformI = zeros(size(newWaveformQ));
            elseif ~isempty(newWaveformI) && isempty(newWaveformQ)
                newWaveformQ = zeros(size(newWaveformI));
            end
            if sum(size(newWaveformI) ~= [length(f),length(duration)]) || sum(size(newWaveformQ) ~= [length(f),length(duration)])
                error('Inputs must have the right size')
            end
            
            % change from duration to points
            points = round(duration*obj.SampleRate);
            timeChange = abs(sum((points)/obj.SampleRate - duration));
            if timeChange > 1e-6
                warning('Change in time due to round ups detected when converting waveform duration to # of points. total change of %f',timeChange)
            end
            
            % rewrite waveform as points for the AWG, (with each point of period 1/obj.SampleRate)
            t = (0 : 1/(sum(points)) : 1-1/(sum(points))).*sum(duration) +t0;
            if length(t) < obj.MIN_WAVEFORM_LENGTH
                error('Tha minimum point at waveform the AWG can get is %d , The current is %d. (At maximum sumple rate is consider toal waveform duration of 0.15us)',obj.MIN_WAVEFORM_LENGTH,length(t))
            end
            waveformI = zeros(length(f),sum(points));
            waveformQ = waveformI;
            waveformIFunc = waveformI;
            waveformQFunc = waveformI;
            startIndex = 1;
            for k = 1:length(points)
                endIndex = startIndex + points(k) - 1;
                for j=1:length(f)
                    waveformI(j,startIndex:endIndex) = newWaveformI(j,k);
                    waveformQ(j,startIndex:endIndex) = newWaveformQ(j,k);
                    waveformIFunc(j,startIndex:endIndex) = FuncI(j,k);
                    waveformQFunc(j,startIndex:endIndex) = FuncQ(j,k);
                end
                startIndex = endIndex+1;
            end
            for j=1:length(f)
                if spinLock
                    waveformI(j,:) = (waveformI(j,:) + waveformIFunc(j,:) .* sin(alpha*sin(2*pi*RF*t))) .* sin(2*pi*f(j)*(t));
                    waveformQ(j,:) = (waveformQ(j,:) + waveformQFunc(j,:) .* cos(alpha*sin(2*pi*RF*t))) .* cos(2*pi*f(j)*(t));
                else % spin nuLock
                    waveformI(j,:) = (waveformI(j,:) + waveformIFunc(j,:) .* cos(alpha*sin(2*pi*RF*t))) .* sin(2*pi*f(j)*(t));
                    waveformQ(j,:) = (waveformQ(j,:) + waveformQFunc(j,:) .* sin(alpha*sin(2*pi*RF*t))) .* cos(2*pi*f(j)*(t));
                end
            end
            waveformI = sum(waveformI,1);
            waveformQ = sum(waveformQ,1);
            WF = waveformI + waveformQ;
            index = AddWaveformByPoints(obj,WF);
        end

        function [index, time, total_WF] = AddSequence(obj, WaveformsIndeces, Repeats, Triggers, WaveformDrawing, GotoTargets)
            % Add a new sequence to obj.sequences:
            % WaveformsIndeces - vector of the waveforms indexes.
            % Repeats - vector of number of reapets for any waveform.
            % Triggers - vector of boolean values for all waveform if waiting to trigger
            % WaveformDrawing - boolean. If to draw the full sequence points.
            % GotoTargets - optional vector, one entry per step. 0 means "leave
            %      the default flow" (fall through to the next step, with the
            %      LAST step of the sequence looping back to the first -- the
            %      original behaviour). A positive value N makes that step jump
            %      to step N instead.
            %
            %      This is what makes a single sequence able to hold many
            %      independently-addressable blocks: point each block's last
            %      step back at its own first step and the sequencer parks on
            %      that block until something forces a jump elsewhere. See
            %      jumpToStep and AWGExpESR's useSequenceJump path.

            if nargin < 5
                WaveformDrawing = 0;
            end
            if ~exist('GotoTargets', 'var') || isempty(GotoTargets)
                GotoTargets = zeros(1, length(WaveformsIndeces));
            end
            if length(GotoTargets) ~= length(WaveformsIndeces)
                error('GotoTargets must have one entry per step (%d), got %d', ...
                      length(WaveformsIndeces), length(GotoTargets));
            end
            if length(WaveformsIndeces) > obj.MAX_SEQUENCE_STEPS
                error('A sequence may have at most %d steps, this one has %d', ...
                      obj.MAX_SEQUENCE_STEPS, length(WaveformsIndeces));
            end

            s = obj.sequences;
            s.waveforms{end+1} = WaveformsIndeces;
            s.repeats{end+1}   = Repeats;
            s.trigger{end+1}   = Triggers;
            s.goto{end+1}      = GotoTargets;
            obj.sequences = s;
            index = length(obj.sequences.waveforms);
            
            % drawing the wavefotm
            time = 0;
            total_WF = 0;
            if WaveformDrawing(1) ~= 0
                % ploting waveform:
                total_WF = [];
                for q = 1:length(WaveformsIndeces)
                    total_WF = [total_WF, repmat(obj.waveforms{WaveformsIndeces(q)},1, Repeats(q))];
                end
                time = 0 : 1/obj.SampleRate : (length(total_WF)-1)/obj.SampleRate;
                figure;
                hold on; plot(time,total_WF);
                
                % optional ploting of constant cosine, to control the phase:
                if length(WaveformDrawing) > 1
                    w = WaveformDrawing(1);
                    A = WaveformDrawing(2);
                    if length(WaveformDrawing) == 3
                        t0 = WaveformDrawing(3);
                    else
                        t0 = 0;
                    end
                    time = 0 : 1/(10*obj.SampleRate) : (length(total_WF)-1)/obj.SampleRate;
                    f = @(t) A*cos(w*(t-t0));
                    plot(time, f(time), 'red');
                end
            end
        end
        
        function waveform = getWaveform(obj, index)
            if nargin == 1
                waveform = obj.waveforms;
            else
                waveform = obj.waveforms{index};
            end
        end

        function sequence = getSequence(obj, index)
            if nargin == 1
                sequence = obj.sequences;
            else
                S.waveforms = obj.sequences.waveforms{index};
                S.repeats   = obj.sequences.repeats{index};
                S.trigger   = obj.sequences.trigger{index};
                sequence = S;
            end
        end
        
        
        function LoadWaveform(obj)
            % Load stored waveforms to memory.
            %
            % THE TIME BREAKDOWN IS THE POINT OF THE REPORTING HERE. A slow load is one of
            % three different problems and they have nothing to do with each other:
            %
            %   BUILD    MATLAB evaluating the waveform functions into sample vectors. Scales
            %            with total samples. Fix by shortening waveforms or thinning tau.
            %   TRANSFER pushing those samples over the wire. Scales with MB. Fix by changing
            %            transport (see setTransport) or sending less.
            %   SETTLE   the instrument digesting what it was sent, after the last byte. Fixed
            %            cost per load, largely independent of size.
            %
            % Printing one total tells you none of that, which is why all three are separated.
            nWf = length(obj.waveforms);
            nSamp = 0;
            for k = 1:nWf; nSamp = nSamp + numel(obj.waveforms{k}); end
            mb = 4*nSamp/2^20;              % single precision on the wire
            if obj.prog.buildCount == nWf && isfinite(obj.prog.buildSeconds)
                fprintf('Uploading to the AWG...');      % the build line just gave the counts
            else
                fprintf('Uploading %d waveforms, %s samples, %.1f MB...', ...
                        nWf, obj.humanCount(nSamp), mb);
            end
            tXfer = tic;
            if obj.skipRedundantDeletes
                % We are about to rely on 'wlist:waveform:delete ALL' being
                % authoritative, so that the per-waveform delete inside
                % WriteWaveformToAWG can be skipped. A waveform that is still
                % referenced by an assigned sequence on a running channel can
                % survive the bulk delete -- which is very likely why the
                % per-waveform delete was there in the first place. So drop the
                % references first: stop the instrument and clear the sequence
                % list. Costs 2 commands, saves one per waveform.
                %
                % Stopping here is safe: the caller (Experiment.prepare ->
                % LoadAWG) is followed by perform(), whose EnsureRunning sees
                % RSTate == Stopped and re-arms the AWG before the first point.
                obj.sendCommand('awgcontrol:stop:immediate');
                obj.deleteSequences;
            end
            obj.deleteWaveforms;
            for k = 1:length(obj.waveforms)
                currentWaveform = obj.waveforms{k};
                if ~(max(abs(currentWaveform)) <= 1)
                     error('The max value allowed to insert to the AWG is +-1, the current is : %g' ,max(abs(currentWaveform)))
                end
                try
                    obj.WriteWaveformToAWG(currentWaveform, k, obj.skipRedundantDeletes, true)
                catch err
                    fprintf('Caught error trying to write waveform %d\n', k);
                    rethrow(err)
                end
            end
            xferSeconds = toc(tXfer);
            tSettle = tic;
            obj.sendCommand('*WAI');    % one barrier instead of one per waveform
            % *WAI only orders the INSTRUMENT's queue -- it does not block MATLAB.
            % Without an explicit wait here the sweep starts while the AWG is
            % still digesting the upload, and that cost reappears smeared across
            % the per-point re-arms: measured 15.7 s of *OPC? wait spread over 97
            % points, which also slowed each arm (691 ms vs 459 ms when idle).
            % Wait once, in one place, where it is visible.
            obj.deviceReady(obj.LOAD_SETTLE_TIMEOUT);
            setSeconds = toc(tSettle);
            fprintf(' done\n');
            obj.reportLoad('waveform', nWf, nSamp, mb, xferSeconds, setSeconds);
            obj.warnOnInstrumentErrors('after LoadWaveform');
        end


        function reportLoad(obj, what, nItems, nUnits, mb, xferSeconds, setSeconds)
            % One block, always the same shape, so two runs can be compared by eye.
            build = obj.prog.buildSeconds;
            sameSet = isfinite(build) && obj.prog.buildCount == nItems && strcmp(what, 'waveform');
            total = xferSeconds + setSeconds;
            if sameSet; total = total + build; end

            fprintf('  %-9s %8s', 'total', obj.humanTime(total));
            if nItems > 0
                fprintf('   = %s per %s\n', obj.humanTime(total/nItems), what);
            else
                fprintf('\n');
            end
            if sameSet
                fprintf('    build   %8s   MATLAB computing samples (%s samples/s)\n', ...
                        obj.humanTime(build), obj.humanCount(obj.prog.buildSamples/max(build,eps)));
            end
            if isfinite(mb) && mb > 0
                fprintf('    send    %8s   %.1f MB over the link (%.1f MB/s)\n', ...
                        obj.humanTime(xferSeconds), mb, mb/max(xferSeconds,eps));
            else
                fprintf('    send    %8s   %d steps (%s per step)\n', ...
                        obj.humanTime(xferSeconds), nUnits, obj.humanTime(xferSeconds/max(nUnits,1)));
            end
            fprintf('    settle  %8s   instrument digesting it, after the last byte\n', ...
                    obj.humanTime(setSeconds));

            % Say which term dominates, because the fix differs for each and the numbers alone
            % do not make that obvious at a glance.
            terms = [xferSeconds setSeconds];
            names = {'the link', 'the instrument settling'};
            if sameSet; terms = [build terms]; names = [{'building samples in MATLAB'}, names]; end
            [~, iw] = max(terms);
            if total > 2 && terms(iw) > 0.55*total
                fprintf('    -> %.0f%% of it is %s.\n', 100*terms(iw)/total, names{iw});
            end
        end

        function LoadSequence(obj)
            % Load stored sequences to memory.
            nSeq = length(obj.sequences.waveforms);
            nStep = 0;
            for k = 1:nSeq; nStep = nStep + numel(obj.sequences.waveforms{k}); end
            % Steps, not sequences, is what the instrument actually has to swallow: an Echo
            % sequence at cycles = 8 is 52 steps where a Ramsey sequence is 2, so a sequence
            % count on its own does not predict the time even to an order of magnitude.
            fprintf('Uploading %d sequences, %d steps total...', nSeq, nStep);
            tXfer = tic;
            obj.deleteSequences;
            for k = 1:length(obj.sequences.waveforms)
                obj.WriteSequenceToAWG(obj.sequences.waveforms{k}, obj.sequences.repeats{k}, ...
                                       obj.sequences.trigger{k}, k, obj.skipRedundantDeletes, ...
                                       obj.sequences.goto{k}, true);
            end
            xferSeconds = toc(tXfer);
            tSettle = tic;
            obj.sendCommand('*WAI');    % one barrier for the whole batch instead of one per sequence
            obj.deviceReady(obj.LOAD_SETTLE_TIMEOUT);   % see the note in LoadWaveform
            setSeconds = toc(tSettle);
            fprintf(' done\n');
            obj.reportLoad('sequence', nSeq, nStep, NaN, xferSeconds, setSeconds);
            % Say straight away if anything was rejected. SCPI execution errors
            % are invisible to write(), so a malformed sequence otherwise only
            % shows up later as a mysterious refusal to arm.
            obj.warnOnInstrumentErrors('after LoadSequence');
        end
        
        function set.waveforms(obj, newWaveform)
            % waveforms will be added to obj.wavefoems{end+1}
            % waveform names are stored in AWG memory as arbseg# where # is as given in the waveforms{#}
            if ~isa(newWaveform,'cell')
                error('Input waveforms must be a cell array');
            end
            for k = 1:length(newWaveform)
                WF = newWaveform{k};
                if ~isnumeric(WF)
                    error('input must be numeric')
                end
                if length(WF) < obj.MIN_WAVEFORM_LENGTH
                    error('Tha minimum waveform length the AWG can get is %d samples, The current is %d.', obj.MIN_WAVEFORM_LENGTH, length(t))
                end
                if max(abs(WF)) > 1
                    error('Waveform''s maximum allowed total value is 1. the current maximum total value : %g', max(abs(WF)))
                end
            end
            obj.waveforms = newWaveform;
        end
        
        function set.sequences(obj,newSequence)
            % inserts values to sequence property. Composed of a cell
            % (stores different sequences), with each one having 3 vector fields:
            % waveform numbers - "waveforms", repeats and trigger
            % sequence names are stored in AWG memory as arbseq# where # is
            % as given in the sequence{#}
            if ~isfield(newSequence,'waveforms') ||...
                    ~isfield(newSequence,'repeats')...
                    || ~isfield(newSequence,'trigger')
                error('sequence must include waveforms, repeats and triggers')
            end
            % 'goto' is optional, so anything written before it existed still
            % works: synthesise all-zero targets, which means "default flow".
            if ~isfield(newSequence,'goto') || length(newSequence.goto) ~= length(newSequence.waveforms)
                newSequence.goto = cellfun(@(x) zeros(1,length(x)), newSequence.waveforms, ...
                                           'UniformOutput', false);
            end
            if range([length(newSequence.waveforms), length(newSequence.repeats), length(newSequence.trigger)]) > 0
                error('Input fields must have the same length')
            end
            if ~isa(newSequence.waveforms,'cell') || ...
                    ~isa(newSequence.repeats,'cell') || ~isa(newSequence.trigger,'cell')
                error('input fields must be cell arrays')
            end
            for k = 1:length(newSequence.waveforms) % test that each element in the cell has the currect structure
                % test if input is empty or numeric
                if (~isempty(newSequence.waveforms{k}) && ~isnumeric(newSequence.waveforms{k})) || ...
                        (~isempty(newSequence.repeats{k}) && ~isnumeric(newSequence.repeats{k}))
                    error('input must be numeric')
                end
                % test if input is empty or numeric or logical
                if ~isempty(newSequence.trigger{k}) && ~isnumeric(newSequence.trigger{k}) && ~islogical(newSequence.trigger{k})
                    error('input must be numeric or logical')
                end
                % test equal length
                if range([length(newSequence.waveforms{k}), length(newSequence.repeats{k}), length(newSequence.trigger{k})]) > 0
                    error('Input fields must have the same length')
                end
                % test number of repeats
                if sum(newSequence.repeats{k} < 0) || sum(newSequence.repeats{k} > obj.MAX_REAPETS)
                    error('repeats must be between 0 and %d', obj.MAX_REAPETS)
                end
            end
            obj.sequences.waveforms = newSequence.waveforms;
            obj.sequences.repeats   = newSequence.repeats;
            obj.sequences.trigger   = newSequence.trigger;
            obj.sequences.goto      = newSequence.goto;
        end
        
    end
    
%% User control functions
    methods
        
        function Run(obj)
            tRun = tic;
            obj.deviceReady();                          % waits for previous overlapping command to actually finish
            obj.sendCommand('awgcontrol:run:immediate');
            obj.waitForRunning();                       % waits for the AWG to actually report "Running"
            obj.timing.run      = obj.timing.run + toc(tRun);
            obj.timing.runCalls = obj.timing.runCalls + 1;
        end

        function EnsureRunning(obj)
            % Cheaper alternative to Run() for use between points of a
            % sweep where the AWG is expected to already be armed: only
            % pays the full arm-and-wait cost if the instrument has
            % actually dropped to Stopped (RSTate = 0). Reassigning the
            % active sequence/waveform on an already-armed channel takes
            % effect on the next trigger without a fresh RUN:IMMEDIATE.
            obj.timing.ensureCalls = obj.timing.ensureCalls + 1;

            if obj.armWithoutChecking
                % Don't ask, just arm. See armWithoutChecking -- the question
                % costs more than the answer is worth on this instrument.
                obj.timing.ensureStopped = obj.timing.ensureStopped + 1;
                obj.Run();
                return
            end

            % Time the state query itself. Measured on the AWG70002A over VISA
            % TCPIP: this single query costs ~1.2 s, which made it by far the
            % most expensive thing in the whole sweep once jump mode had removed
            % the re-arms -- 193 calls, 235 s. Earlier versions of this
            % instrumentation timed only Run() and so reported 2.63 s for a run
            % that actually took 296 s. Do not remove this timer.
            tStateQuery = tic;
            state = str2double(obj.ioQuery('AWGControl:RSTate?'));
            obj.timing.stateQuery = obj.timing.stateQuery + toc(tStateQuery);
            if state == obj.RSTATE_STOPPED
                obj.timing.ensureStopped = obj.timing.ensureStopped + 1;
                obj.Run();
            end
        end
        
        function setRunMode(obj, mode, channels)
            % set the run mode of the channels
            % channels: 1, 2, or [1 2]
            if ~exist('channels', 'var')
                channels = 1;
            end
            % NOTE was 'srtcmpi' (typo) in every branch below, so this method
            % threw an undefined-function error if it was ever called -- which
            % is presumably why the run mode has never been configured and the
            % AWG has been left in whatever *RST selected. Fixed to strcmpi.
            if ~strcmpi(mode, obj.RUN_MODE_CONT) && ~strcmpi(mode, obj.RUM_MODE_TRIG) && ~strcmpi(mode, obj.RUN_MODE_TCON)
                error('AWG run mode can be just "%s", "%s", or "%s"', obj.RUN_MODE_CONT, obj.RUM_MODE_TRIG, obj.RUN_MODE_TCON);
            end
            for i = 1:length(channels)
                runModeCommand = sprintf('SOURce%d:RMODe %s', channels(i), mode);
                obj.sendCommand(runModeCommand);
            end
        end
        
        function assignSequence(obj, sequenceNum, channel)
            % assign sequenct to a channel
            if nargin < 3
                channel = 1;
            end
            tAssign = tic;
            assign_command = sprintf('source%d:casset:sequence "%s", 1',channel, num2str(sequenceNum)); % assign sequence to track 1
            obj.sendCommand(assign_command);
            obj.timing.assign      = obj.timing.assign + toc(tAssign);
            obj.timing.assignCalls = obj.timing.assignCalls + 1;
        end
        
        function jumpToStep(obj, stepNum, channel)
            % Force the sequencer to a given step WITHOUT reassigning the asset.
            %
            % This is the whole point of the single-sequence scheme. Measured on
            % the AWG70002A: assigning an asset (assignSequence) STOPS playback,
            % and re-arming costs ~500 ms regardless of how much is loaded, so a
            % per-point reassign cost 509 ms x 382 points = 176 s. A forced jump
            % does not stop the instrument, so it costs one SCPI write (~3 ms).
            %
            % Timed against the same counters as assignSequence so that
            % printInstrumentation keeps working.
            if nargin < 3
                channel = 1;
            end
            tJump = tic;
            obj.sendCommand(sprintf('source%d:jump:force %d', channel, stepNum));
            % Let the instrument actually act on the jump before anyone triggers
            % it. See jumpSettleTime -- omitting this silently measures the
            % previous frequency.
            if obj.jumpSettleTime > 0
                tSettle = tic;
                while toc(tSettle) < obj.jumpSettleTime
                end
            end
            obj.timing.assign      = obj.timing.assign + toc(tJump);
            obj.timing.assignCalls = obj.timing.assignCalls + 1;
        end

        function assignWaveform(obj, waveformNum, channel)
            % assign waveform to a channel
            if nargin < 3
                channel = 1;
            end
            assign_command = sprintf('source%d:casset:waveform "%s"',channel, num2str(waveformNum)); % assign waveform to track 1
            obj.sendCommand(assign_command);
        end
        
        function ChangeSampleRate(obj, Durations)
            % change the sample rate of the AWG to the highest number that divide (1/f) by all the times in the durations vector
            MAX_ERROR = 1e-10;
            f = obj.MAX_SAMPLE_RATE;
            while sum(mod(Durations,1/f)) > MAX_ERROR
                f = f - 1e-6;
                if f < obj.MIN_SAMPLE_RATE
                    error('AWG sample rate not found');
                end
            end
            obj.SampleRate = round(f,6);
        end
        
        function maxfreq = CheckMinTimeFotReapets(obj,f)
            if length(f) == 6
                syms m1 m2 m3 m4 m5 m6
                eqns = [f(1)/m1-f(2)/m2,f(2)/m2-f(3)/m3,f(3)/m3-f(4)/m4,f(4)/m4-f(5)/m5,f(5)/m5-f(6)/m6,f(6)/m6-f(1)/m1];
                vars = [m1,m2,m3,m4,m5,m6];
                sol = solve(eqns,vars);
                y = [(sol.m1) (sol.m2) (sol.m3) (sol.m4) (sol.m5) (sol.m6)] ;
                [~,D] = numden(y);
                multiple = single(max(D));
                y = multiple*single(y);
                if y(1)/f(1) ~= y(2)/f(2) &&  y(3)/f(3) ~= y(4)/f(4) && y(1)/f(1) ~= y(4)/f(4) && y(5)/f(5) ~= y(6)/f(6) && y(1)/f(1) ~= y(6)/f(6)
                    error('can not find corporate pereiod time')
                end
                commonFrequnecy = f(1)/y(1);
            elseif length(f) == 2
                syms m1 m2
                eqns = f(1)/m1-f(2)/m2;
                vars=[m1,m2];
                sol = solve(eqns,vars);
                y = [(sol.m1) (sol.m2)] ;
                [~,D] = numden(y);
                multiple = single(max(D));
                y = multiple*single(y);
                if y(1)/f(1) ~= y(2)/f(2)
                    error('can not find corporate pereiod time')
                end
                commonFrequnecy = f(1)/y(1);
            elseif length(f) == 1
                commonFrequnecy = f;
            else
                error('there is no option to check period for %g frequencies. just for : 1 , 2 or 6 per waveform', length(f))
            end
            maxfreq = gcd(commonFrequnecy*obj.TIME_CONVERT, obj.SampleRate*obj.TIME_CONVERT) / obj.TIME_CONVERT;
        end
        
        function Reset(obj)
            obj.sendCommand('*RST');
            obj.setTimeOut(100);
            obj.SampleRate = obj.MAX_SAMPLE_RATE;
            obj.workChannels = 1;
            obj.ResetStoredData;
        end
        
        function ResetStoredData(obj)
            % Anything that wipes the stored data must also drop the cache, or a
            % later experiment could get a false hit.
            obj.invalidateLoaded();
            obj.waveforms = {};
            s.waveforms = {};
            s.repeats = {};
            s.trigger = {};
            s.goto = {};
            obj.sequences = s;
        end
        
    end

%% Function generator mode functions
    methods
        function command = createCommand(obj, what, value, channel)
            space = repmat(' ',1,~strcmp(value,'?'));   % '?' shold be without space and value with space
            chan =  num2str(channel);
            switch lower(what)
                case {'enableoutput', 'output', 'enabled', 'enable'}
                    switch value
                        % NOTE both compound commands below used a bare ';' as
                        % the separator. In SCPI a bare ';' keeps the parser in
                        % the CURRENT subtree, so 'AWGC:STOP; OUTP1 OFF' was
                        % parsed as AWGControl:OUTP1 and rejected:
                        %
                        %   -113,"Undefined header;Command not found - AWGC:STOP; OUTP1 "
                        %
                        % sendCommand uses fwrite, which cannot see an
                        % instrument-side rejection, so this failed silently on
                        % every call -- meaning disconnect() has never actually
                        % switched the outputs off. Found 20 of these queued on
                        % the AWG70002A (10 disconnects x 2 channels).
                        %
                        % Fixed by using ';:' (the leading colon resets to root)
                        % and by spelling the run/stop headers out in the form
                        % this instrument is known to accept -- awgcontrol:
                        % stop:immediate and :run:immediate are already used
                        % elsewhere in this class and produce no error.
                        case {'1', 'ON'}
                            obj.setValue('mode', 'FGEN');
                            command = ['OUTP', chan, ' ON;:AWGCONTROL:RUN:IMMEDIATE'];
                        case {'0', 'OFF'}
                            command = ['AWGCONTROL:STOP:IMMEDIATE;:OUTP', chan, ' OFF'];
                        case '?'
                            command = ['OUTP', chan, '?; INST:MODE?'];
                        otherwise
                            error('no type of OUTPUT syntax "%s"',value)
                    end
                case {'frequency', 'freq'}
                    if ~strcmp(value, '?')
                        value = value*1e6; % convert MHz to Hz
                    end
                    command = ['FGEN:CHAN', chan, ':FREQ' ,    space , num2str(value)];
                case {'amplitude', 'ampl', 'amp'}
                    command = ['FGEN:CHAN', chan, ':AMPL:POW', space , num2str(value)];
                case {'phase'}
                    command = ['FGEN:CHAN', chan, ':PHAS',     space , num2str(value)];
                case {'mode'}
                    if sum(strncmp(value, {obj.MODE_FG, obj.MODE_AWG, '?'},10))
                        command = ['INST:MODE', space, value];
                        if ~strcmp(value, '?')
                            obj.mode = value;
                        end
                    else
                        error('no type of mode %s', value)
                    end
                otherwise
                    error('Unknown command type %s', what)
            end
        end
    end
    
%% AWG programming:
    methods (Access = private)
        
        function deleteWaveforms(obj)
            command = 'wlist:waveform:delete ALL';
            obj.sendCommand(command); 
        end
        
        function deleteSequences(obj)
            command = 'slist:sequence:delete ALL';
            obj.sendCommand(command); 
        end
       
        function WriteWaveformToAWG(obj, waveform, waveformName, alreadyCleared, deferWait)
            % loads a single waveform (a vector of numbers) to trace number 'traceNumber' of the AWG.
            % Segment is in bytesarray, as given by the 'arbitraryWaveform'.
            %
            % alreadyCleared - logical, optional (default false). True when the
            %   caller has already wiped the whole waveform list, in which case
            %   the per-waveform delete below is redundant and is skipped.

            if ~exist('alreadyCleared', 'var')
                alreadyCleared = false;
            end
            if ~exist('deferWait', 'var')
                deferWait = false;      % true = caller will send one *WAI at the end
            end

            if max(abs(waveform)) > 1
                error('The maximal total value for waveform must be less than 1')
            end
            if length(waveform) < obj.MIN_WAVEFORM_LENGTH
                error('Tha minimum point at waveform the AWG can get is %d , The current is %d. (At maximal sumple rate is consider toal waveform duration of 0.15us)',obj.MIN_WAVEFORM_LENGTH,length(waveform))
            end                
            
            Waveform_Name_1 = sprintf('%g',waveformName);                   % Waveform name
            Waveform_Data_1 = single(waveform);                             % Waveform data

            % Generate commands:
            delete_waveform = sprintf('wlist:waveform:delete "%s"', Waveform_Name_1);   %delete any previous instances of the waveform we're sending
            data_prefix     = sprintf('wlist:waveform:data "%s",', Waveform_Name_1);    %the binblock header and samples get appended to this
            new_waveform    = sprintf('wlist:waveform:new "%s", %d', Waveform_Name_1, length(Waveform_Data_1));

            % Send commands:
            if ~alreadyCleared
                obj.sendCommand(delete_waveform);       % send delete waveform command
            end
            obj.sendCommand(new_waveform);              % send new waveform command
            obj.ioWriteBinaryBlock(data_prefix, Waveform_Data_1);   % header + samples as ONE message
            if ~deferWait
                obj.sendCommand('*WAI');
            end
        end

        function WriteSequenceToAWG(obj, WaveformsIndeces, Repeats, Triggers, traceNum, alreadyCleared, GotoTargets, deferWait)
            % Uploads a sequence of waveforms (vectors), with a given number of repeats for each one and a trigger mode (0 or 1).
            % This will be stored as "arbseq#' where # is the traceNum
            %
            % alreadyCleared - logical, optional (default false). True when the
            %   caller has already wiped the whole sequence list, in which case
            %   the per-sequence delete below is redundant and is skipped.

            if ~exist('alreadyCleared', 'var')
                alreadyCleared = false;
            end
            if ~exist('GotoTargets', 'var') || isempty(GotoTargets)
                GotoTargets = zeros(1, length(WaveformsIndeces));
            end
            if ~exist('deferWait', 'var')
                deferWait = false;      % true = caller will send one *WAI at the end
            end

            if range([length(WaveformsIndeces), length(Repeats), length(Triggers), length(GotoTargets)]) > 0
                error('All of the vectors must be in the same length')
            end

            seq_name = sprintf('%g',traceNum);
            steps = length(WaveformsIndeces);

            if ~alreadyCleared
                delete = sprintf('slist:seq:delete "%s"', seq_name);
                obj.sendCommand(delete);
            end
            create_seq = sprintf('slist:seq:new "%s", %d, 1', seq_name, steps); %2 steps, 1 track
            obj.sendCommand(create_seq);
            
            % Build every step command first, then send them as a handful of
            % compound messages instead of 3-4 separate writes per step. For a
             % 3-step sequence that is 2 messages instead of 12; for the
            % 573-step single sequence it is ~44 instead of 1912.
            cmds = cell(1, 0);
            for r = 1:steps
                if (Triggers(r) == 1) && (Repeats(r) ~= 1)
                    error('last sequence is a loop. This is not supported (trig problem)')
                end

                currentWavform = sprintf('%d',WaveformsIndeces(r));
                cmds{end+1} = sprintf('slist:seq:step%d:tasset1:wav "%s", "%s"',r, seq_name, currentWavform); %#ok<AGROW>
                cmds{end+1} = sprintf('slist:seq:step%d:rcount "%s", %d',r ,seq_name, Repeats(r)); %#ok<AGROW>
                if Triggers(r) == 0
                    cmds{end+1} = sprintf('slist:seq:step%d:winput "%s", off',r, seq_name); %#ok<AGROW>
                else
                    cmds{end+1} = sprintf('slist:seq:step%d:winput "%s", atrigger',r, seq_name); %#ok<AGROW>
                end

                if GotoTargets(r) > 0
                    % Explicit per-step jump target. Verified accepted on the
                    % AWG70002A (FV 8.0.0262.0): a numeric target returns
                    % 0,"No error", so a sequence can be carved into blocks that
                    % each loop on themselves.
                    cmds{end+1} = sprintf('slist:seq:step%d:goto "%s", %d', ...
                                          r, seq_name, GotoTargets(r)); %#ok<AGROW>
                end
            end
            if all(GotoTargets == 0)
                % Original behaviour: no explicit targets, so make the whole
                % sequence loop by sending the last step back to the first.
                cmds{end+1} = sprintf('slist:seq:step%d:goto "%s", first', steps, seq_name);
            end
            obj.sendCommandBatch(cmds);

            if ~deferWait
                obj.sendCommand('*WAI');    %Make sure no other commands are exectued until arb is done downloadin
            end
        end
        
        function waitForRunning(obj)
            % In triggered run mode the AWG arms into "Waiting for Trigger"
            % (RSTate = 1) rather than "Running" (RSTate = 2) until the
            % external trigger actually arrives, so both count as success.
            % Only "Stopped" (RSTate = 0) means it failed to arm/run.
            timer = tic;
            state = str2double(obj.ioQuery('AWGControl:RSTate?'));
            while state == obj.RSTATE_STOPPED
                if toc(timer) > obj.RUN_TIMEOUT
                    obj.timing.armWait = obj.timing.armWait + toc(timer);
                    % Say WHY. A bare "did not enter RUN mode" cannot distinguish
                    % "this sequence is slow to arm" from "the instrument rejected
                    % it", and SCPI execution errors are invisible to fwrite/write
                    % -- they only sit in the error queue. So read the queue and
                    % put it in the exception.
                    error(['AWG time out - device does not enter into RUN mode after %.1f s ' ...
                           '(RSTate = %g, 0 = stopped). %s'], ...
                          toc(timer), state, obj.describeInstrumentErrors());
                end
                obj.pollDelay();
                state = str2double(obj.ioQuery('AWGControl:RSTate?'));
            end
            armSeconds = toc(timer);
            obj.timing.armWait = obj.timing.armWait + armSeconds;
            if armSeconds > 2
                % Worth knowing about: a long arm points at sequence size rather
                % than transport, and it is the thing that decides whether the
                % concatenated modes are practical.
                fprintf('AWG took %.1f s to arm (RSTate = %g).\n', armSeconds, state);
            end
        end

        function warnOnInstrumentErrors(obj, context)
            % Non-fatal check: drain the error queue and warn if it is not empty.
            % One query when clean, which is ~2 ms on visadev.
            %
            % Plain warning(), NOT obj.sendWarning. FrequencyGenerator supplies
            % sendError, but there is no sendWarning anywhere in this hierarchy -- the
            % call raised "Unrecognized method, property, or field 'sendWarning'". And
            % because it only fires when the queue is NON-empty, it sat here latent until
            % the AWG actually had something to report, at which point it turned a
            % non-fatal notice into an abort in the middle of LoadWaveform. The rest of
            % this class uses plain warning() for exactly this; match it.
            %
            % The '%s' form matters: SCPI error strings can contain characters like %
            % that warning() would otherwise read as format specifiers.
            txt = obj.describeInstrumentErrors();
            if ~contains(txt, 'queue is empty')
                warning('NVLab:awgInstrumentError', '%s', ...
                        sprintf('AWG %s: %s', context, txt));
            end
        end

        function txt = describeInstrumentErrors(obj, maxEntries)
            % Drain up to maxEntries from the SCPI error queue and format them for
            % an exception message. Wrapped in try/catch because this runs on the
            % failure path, where the instrument may already be unresponsive.
            if ~exist('maxEntries', 'var') || isempty(maxEntries)
                maxEntries = 5;
            end
            lines = cell(1, 0);
            try
                for k = 1:maxEntries
                    resp = char(strtrim(obj.ioQuery('SYSTem:ERRor?')));
                    if isempty(resp) || strncmp(resp, '0,', 2) || contains(resp, 'No error')
                        break
                    end
                    lines{end+1} = resp;     %#ok<AGROW>
                end
            catch
                lines{end+1} = '(could not read the error queue)';
            end
            if isempty(lines)
                txt = 'The instrument error queue is empty, so nothing was rejected.';
            else
                txt = ['The instrument reports: ' strjoin(lines, ' | ')];
            end
        end

        function deviceReady(obj, timeoutSeconds)
            % tests if the the previous overlapping command was executed
            if ~exist('timeoutSeconds', 'var') || isempty(timeoutSeconds)
                timeoutSeconds = obj.OPC_TIMEOUT;
            end

            % *OPC? BLOCKS the instrument until all pending operations complete,
            % so the VISA READ has to be allowed to wait at least as long as our
            % own loop does. Bounding only the loop is not enough: the read times
            % out first and throws before the loop ever gets to check. That is
            % what happened building a 288-step concatenated sequence, where the
            % instrument stayed busy longer than the ambient timeout.
            oldTimeout = obj.ioGetTimeout();
            raised = false;
            if oldTimeout < timeoutSeconds
                obj.ioSetTimeoutRaw(timeoutSeconds + 5);
                raised = true;
            end

            timer = tic;
            try
                while ~contains(obj.ioQuery('*OPC?'), '1')
                    if toc(timer) > timeoutSeconds
                        error('AWG time out - overlapping commands has not completed within %g s', timeoutSeconds);
                    end
                    obj.pollDelay();
                end
            catch err
                obj.timing.opcWait = obj.timing.opcWait + toc(timer);
                if raised; obj.ioSetTimeoutRaw(oldTimeout); end
                rethrow(err);
            end
            obj.timing.opcWait = obj.timing.opcWait + toc(timer);
            if raised; obj.ioSetTimeoutRaw(oldTimeout); end
        end

        function verifyOperationComplete(obj)
            % Wait for the VISA write to drain. For a synchronous fwrite this is
            % already 'idle' on entry, so the loop body normally never executes.
            % Bounded now, so a wedged transfer raises an error instead of
            % hanging MATLAB forever.
            if obj.useVisadev
                % write()/writeline() on a visadev object are synchronous -- there
                % is no asynchronous transfer to drain, and no TransferStatus
                % property. This whole method was ~5.8 s of legacy-interface
                % property-access overhead per sweep.
                return
            end
            if strcmp(obj.v.TransferStatus, 'idle')
                return
            end
            timer = tic;
            while ~strcmp(obj.v.TransferStatus, 'idle')
                if toc(timer) > obj.TRANSFER_TIMEOUT
                    error('AWG time out - VISA transfer did not complete within %g s', obj.TRANSFER_TIMEOUT);
                end
                obj.pollDelay();
            end
        end

        %% ---- transport layer -------------------------------------------------
        % Everything that touches the instrument handle goes through here, so the
        % legacy visa() interface and the modern visadev() interface differ in one
        % place instead of twenty. Selected by obj.useVisadev.

        function ioOpen(obj)
            if obj.useVisadev
                if isempty(obj.v) || ~isvalid(obj.v)
                    try
                        obj.v = visadev(obj.address);
                    catch err
                        error(['visadev could not open "%s": %s\n' ...
                               'Run visadevlist to see the resources MATLAB can see, and set ' ...
                               'awg.useVisadev = false to go back to the legacy interface.'], ...
                              obj.address, err.message);
                    end
                    % Binary waveform samples are written as raw little-endian
                    % float32, matching what the legacy fwrite(...,'single') sent
                    % from a PCWIN64 host. Set it explicitly rather than inherit.
                    obj.v.ByteOrder = "little-endian";
                end
            else
                if isempty(obj.v) || ~isvalid(obj.v)
                    obj.v = visa(obj.VISA_BRAND, obj.address, ...
                                 'InputBufferSize', obj.BUFFER, 'OutputBufferSize', obj.BUFFER);
                end
                if strcmp(obj.v.Status, 'closed')
                    fopen(obj.v);
                end
            end
        end

        function tf = ioIsOpen(obj)
            tf = ~isempty(obj.v) && isvalid(obj.v);
            if tf && ~obj.useVisadev
                tf = strcmp(obj.v.Status, 'open');
            end
        end

        function ioClose(obj)
            if isempty(obj.v) || ~isvalid(obj.v)
                return
            end
            if obj.useVisadev
                delete(obj.v);      % a visadev object closes when deleted
                obj.v = [];
            else
                if strcmp(obj.v.Status, 'open')
                    fclose(obj.v);
                end
            end
        end

        function ioDiscard(obj)
            % Close AND destroy the handle, leaving obj.v empty so that ioOpen
            % will construct a fresh object of whichever class the current
            % useVisadev setting calls for. Use this when switching transports;
            % ioClose alone is not enough (see set.useVisadev).
            if isempty(obj.v)
                return
            end
            if isvalid(obj.v)
                try
                    if ~obj.useVisadev && strcmp(obj.v.Status, 'open')
                        fclose(obj.v);
                    end
                catch
                    % Already unusable -- nothing to salvage, just destroy it.
                end
                delete(obj.v);
            end
            obj.v = [];
        end

        function ioWriteText(obj, command)
            % One SCPI message, no trailing terminator -- the VISA END indicator
            % terminates it. This matches exactly what fwrite() did before.
            if obj.useVisadev
                write(obj.v, command, 'char');
            else
                fwrite(obj.v, command);
            end
        end

        function ioWriteSingles(obj, data)
            % Raw float32 payload with no command prefix. Kept for the legacy
            % sendCommand(...,'single') path; the visadev path builds the whole
            % message in ioWriteBinaryBlock instead.
            if obj.useVisadev
                write(obj.v, single(data), 'single');
            else
                fwrite(obj.v, data, 'single');
            end
        end

        function ioWriteBinaryBlock(obj, commandPrefix, data)
            % Send  <prefix>#<n><nbytes><raw float32>  as ONE instrument message.
            %
            % The legacy path had to drop EOIMode to 'off' so the header and the
            % samples formed a single message, then put it back. visadev has no
            % EOIMode for TCPIP, so instead the whole message is assembled as one
            % uint8 array and written in a single call -- byte-for-byte the same
            % thing on the wire, without the mode juggling.
            data = single(data);
            if obj.useVisadev
                bytes  = typecast(data(:).', 'uint8');      % native LE on PCWIN64
                lenStr = num2str(numel(bytes));
                header = ['#' num2str(length(lenStr)) lenStr];
                write(obj.v, [uint8(commandPrefix), uint8(header), bytes], 'uint8');
                obj.commandCount = obj.commandCount + 1;
                obj.commDelay();
            else
                wfm_bytes = numel(data) * 4;
                header = sprintf('#%d%d', length(num2str(wfm_bytes)), wfm_bytes);
                obj.v.EOIMode = 'off';
                obj.sendCommand([commandPrefix header]);
                obj.v.EOIMode = 'on';
                obj.sendCommand(data, 'single');
            end
        end

        function response = ioQuery(obj, command)
            % A TIMED-OUT QUERY MUST FLUSH, or it corrupts every read after it.
            %
            % A SCPI query cannot answer until the instrument has finished processing
            % everything queued ahead of it, and sendCommand's fwrite returns as soon as
            % the bytes leave the PC. So MATLAB can run far ahead -- during a jump-mode
            % sweep the arm-state check at point 26 was measured waiting ~30 s for a
            % 25-point backlog of forced jumps to drain. When that wait exceeds the VISA
            % timeout the read aborts, which is the visible failure.
            %
            % The invisible failure is worse. The instrument eventually FINISHES and
            % sends its reply, which then sits unread in the input buffer. The next query
            % collects that stale answer, and every reading after it is shifted by one
            % question -- run state, output state, frequency readback, all of them
            % answering the previous query while looking perfectly plausible. This is the
            % same mechanism that once made the instrument report FGEN mode while it was
            % demonstrably playing sequences in AWG mode.
            %
            % So on timeout: flush, drain the error queue, and rethrow with a message
            % that names the real cause. Failing loudly is the whole point -- the
            % alternative is silently wrong data.
            try
                if obj.useVisadev
                    % writeread returns a STRING, whereas the legacy query returns a
                    % char array. Normalise to char so every caller behaves the same
                    % on both transports -- a cell of strings is NOT valid input to
                    % strjoin, which crashed describeInstrumentErrors on the one code
                    % path that most needed to work.
                    response = char(writeread(obj.v, command));
                else
                    response = query(obj.v, command);
                end
            catch err
                if contains(lower(err.message), 'timeout') || ...
                   contains(lower(err.identifier), 'timeout')
                    obj.ioFlushInput();
                    error('NVLab:awgQueryTimeout', ...
                        ['AWG query "%s" timed out after %g s.\n' ...
                         'The input buffer has been flushed, so later reads are not ' ...
                         'silently shifted by one.\n' ...
                         'A query cannot answer until the instrument has drained every ' ...
                         'command queued before it, so this usually means the PC ran ' ...
                         'ahead of the AWG -- most often a jump-mode sweep accumulating ' ...
                         'forced jumps faster than the instrument executes them. Either ' ...
                         'synchronise per point (armCheckEvery = 1, or a large enough ' ...
                         'jumpSettleTime) or raise the timeout with ioSetTimeout.'], ...
                        command, obj.ioGetTimeout());
                else
                    rethrow(err);
                end
            end
        end

        function ioFlushInput(obj)
            % Discard anything waiting to be read, then empty the SCPI error queue.
            % Runs on the failure path, so every step is individually guarded: the
            % instrument may already be unresponsive and a flush that throws would
            % replace a diagnosable timeout with a confusing secondary error.
            try
                if obj.useVisadev
                    flush(obj.v);
                else
                    flushinput(obj.v);
                end
            catch
            end
            % Short timeout for the drain: if the instrument is still busy there is no
            % point waiting the full ambient timeout again, once per queued entry.
            oldTimeout = [];
            try
                oldTimeout = obj.ioGetTimeout();
                obj.ioSetTimeoutRaw(30);
            catch
            end
            for k = 1:10
                try
                    r = strtrim(obj.ioQueryRawNoFlush('SYSTEM:ERROR?'));
                catch
                    break
                end
                if isempty(r) || startsWith(r, '0') || startsWith(r, '+0')
                    break
                end
            end
            if ~isempty(oldTimeout)
                try; obj.ioSetTimeoutRaw(oldTimeout); catch; end
            end
        end

        function response = ioQueryRawNoFlush(obj, command)
            % The plain query, with no timeout handling. Used by ioFlushInput so that
            % draining cannot recurse back into the flush path.
            if obj.useVisadev
                response = char(writeread(obj.v, command));
            else
                response = query(obj.v, command);
            end
        end

        function value = ioReadLine(obj)
            if obj.useVisadev
                value = readline(obj.v);
                value = char(value);
            else
                value = fscanf(obj.v, '%s');
            end
        end

        function ioSetTimeout(obj, seconds)
            % Public-facing setter path, keeps the legacy 30..1000 s validation
            % that setTimeOut already applied.
            obj.ioSetTimeoutRaw(seconds);
        end

        function ioSetTimeoutRaw(obj, seconds)
            if isempty(obj.v) || ~isvalid(obj.v)
                return
            end
            if obj.useVisadev
                obj.v.Timeout = seconds;
            else
                obj.v.timeout = seconds;
            end
        end

        function seconds = ioGetTimeout(obj)
            if isempty(obj.v) || ~isvalid(obj.v)
                % NOMINAL_TIMEOUT, not OPC_TIMEOUT. Callers use this to save a value
                % they will later restore, so returning the 2 s *OPC? probe timeout on a
                % failed read would make every save/restore pair a permanent downgrade of
                % the ambient timeout -- and the symptom would be exactly what we saw: a
                % correct wait killed by a deadline nobody set.
                seconds = obj.NOMINAL_TIMEOUT;
                return
            end
            if obj.useVisadev
                seconds = obj.v.Timeout;
            else
                seconds = obj.v.timeout;
            end
        end

        function commDelay(obj)
            % Communication Delay -- inserted after every SCPI write.
            %
            % Busy-wait on purpose: MATLAB's pause() granularity on Windows is
            % much coarser than the few ms wanted here, so pause(0.005) can
            % easily cost ~15 ms and would be slower than what it replaces.
            % The way to make this cheap is to set obj.commDelayTime = 0, which
            % skips it outright.
            if obj.commDelayTime <= 0
                return
            end
            delay = tic;
            while toc(delay) < obj.commDelayTime
            end
        end

        function pollDelay(obj)
            % Delay between consecutive polls of an instrument status query.
            % Kept separate from commDelay so that commDelayTime = 0 does not
            % turn the status-polling loops into tight spins that hammer VISA.
            delay = tic;
            while toc(delay) < obj.POLL_DELAY
            end
        end
    end
    
%%
    methods (Static)
        
        function RegisterPgChannels(names, channels, delays)
            PG = getObjByName(PulseGenerator.NAME);
                if isempty(PG); throwBaseObjException(PulseGenerator.NAME); end
            for i = 1:length(names)
                channel = Channel.Digital(names{i}, channels(i), delays(2*i-1), delays(2*i));
                PG.registerChannel(channel);
            end
        end
        
        function obj = getInstance(struct)
            missingField = FactoryHelper.usualChecks(struct, ...
                FrequencyGeneratorTektronixAWG.NEEDED_FIELDS);
            if ~isnan(missingField)
                EventStation.anonymousError(...
                    'Trying to create a Tektronix AWG, encountered missing field - "%s". Aborting',...
                    missingField);
            end
            struct = FactoryHelper.supplementStruct(struct, FrequencyGeneratorTektronixAWG.OPTIONAL_FIELDS);
            
            name = [FrequencyGeneratorTektronixAWG.TYPE, '-', struct.serialNumber];
            frequencyLimits = [struct.minFrequency, struct.maxFrequency];
            amplitudeLimits = [struct.minAmplitude, struct.maxAmplitude];
            if ~isempty(struct.trigger2ChannelName)
                triggersNames = {struct.trigger1ChannelName, struct.trigger2ChannelName};
                channels = [struct.trigger1Channel, struct.trigger2Channel];
                delays = [struct.onDelay1, struct.offDelay1, struct.onDelay2, struct.offDelay2];
            else
                triggersNames = {struct.trigger1ChannelName};
                channels = struct.trigger1Channel;
                delays = [struct.onDelay1, struct.offDelay1];
            end
            
            obj = FrequencyGeneratorTektronixAWG(name, struct.address, frequencyLimits, amplitudeLimits, triggersNames, channels, delays);
            addBaseObject(obj);
        end
        
        function Vpp = dBmToVpp(P_dBm)
            P_Watt = 1e-3*(10.^(P_dBm/10));
            Vpp = 2*sqrt(P_Watt*50*2);
        end

    end
end