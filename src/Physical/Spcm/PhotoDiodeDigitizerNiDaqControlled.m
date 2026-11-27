classdef PhotoDiodeDigitizerNiDaqControlled < Spcm & NiDaqControlled
    % Scan fix, 2026-11-23: the detector owns the scan data.
    % The stage moves to a pixel, waits until it is there, and calls
    % acquireScanPixel(tPixel). readFromScan returns the line. The stage
    % no longer prepares the digitizer or writes into this object.
    %   scanAcquisitionMode = 'software' (default): one forced acquisition
    %     per pixel. Checked on the hardware (scanDebug, 2026-11-23 14:19).
    %   scanAcquisitionMode = 'trigger': the digitizer is armed once per
    %     line for one segment per pixel, each started by a rising edge on
    %     the pulse generator's 'detector' channel, and the line is read at
    %     once. On Setup 7 the digitizer also triggers when that line goes
    %     low after a segment has ended, so the pulse is now ended while
    %     the segment is still being recorded (as in the experiments), and
    %     a line where the digitizer finishes early is refused. Not yet
    %     checked on the hardware; on 2026-11-23 it was also slower than
    %     software mode.
    
    properties (Access = protected)
        % Backup, for NiDaq reset
        isEnabled   % logical
        gateOnTic = []      % tic of the last off -> on switch of the photodiode gate
        
        % For time measure
        voltageIntegrationTime
        nTimeIntegration
        voltageTimeTask
        
        % For scanning
        nScanIntegration            % pixels per line
        scanningStageName
        pixelTime                   % s, from prepareCountByStage (1D) or the stage
        scanPixelIndex = 0          % pixels measured so far in this line
        scanArmed = false           % digitizer configured/armed for this line
        scanArmedPixelTime = []     % s, pixel time the digitizer is armed for
        scanSegmentSeconds = 0      % s, length of one acquired segment (trigger mode)
        scanArmTic                  % tic at arming (trigger mode)
        scanTriggerTimes = []       % s since arming, one per pixel trigger (trigger mode)
        
        % For Experiment
        nExpIntegration
        expTimeoutTime
        measureExpTask
        
        % Channel Names
        niDaqGateChannelName            % torn on the photodiode power supply
        niDaqOffsetChannelName
    end
    
    properties (SetAccess = protected)
        current1DScanVoltageValues = [];    % last scan line, one value per pixel (filled by this object only)
        current1DScanSterrValues = [];
    end

    properties
        % ---- scan acquisition, can be changed from the command line ----
        scanAcquisitionMode = 'software'    % 'software' or 'trigger', see the class header
        scanTriggerChannel = 'detector'     % pulse generator channel wired to the digitizer trigger input
        scanTriggerMargin = 0.003           % s. Extra wait after a pixel's segment before the stage may move on
        scanMaxSamplesPerPixel = 1e5        % caps the digitizer sample rate (<= 1 MHz) for long pixel times
        gateSettleTime = 0.1                % s. Minimum time between switching the photodiode gate on and the first scan pixel
        
        currentChannelGUI = 1;     % The channnel to read by time and stage. we have 3 channels: #1: fluorescence, #2: laser sampling, #3: laser block (optional)
        % Options:
        % 1: normal input           % #1
        % 2: reference              % #2
        % 3: both                   % #1 and #2
        % 4: diffrential input      % #1 - #2
        % 5: referenced input       % #1 / <#2>
        % 6: laser block            % #3
        
        defultInputConfiguration
        inputConfiguration
        digitizerChannelNames
        digitizerChannelValues
        segmentDuration
        detectionWithGI = 0;

        % For cases when required the full acquired data:
        fullDataAcquisition = 0;
        recordVoltageSpan = 0       % boolean. When true the object saves the measured voltage histogram
        voltageHistogram            % struct with the bins and values of the measures voltage
        daqVolatgeOffset = [];      % Optional way to generate DC voltage from DAQ as an offset voltage to differential input.
    end
    
    properties (Dependent)
        activeChannelsGUI
        activeChannelsMeas
    end
    
    properties (Constant, Hidden)
        NEEDED_FIELDS = [Spcm.SPCM_NEEDED_FIELDS, {'nidaq_channel_gate', 'digitizer_channel_names', 'digitizer_channel_values', ...
                         'defultd_input_configuration', 'voltage_offset', 'max_voltage', 'input_impedance'}];
        OPTIONAL_FIELDS = {'nidaq_channel_min_val', 'nidaq_channel_max_val', 'channel_delay', 'nidaq_channel_gate', 'nidaq_channel_offset'};
        INPUT_CONFIGURATIONS = {'nomral', 'differential', 'referenced'};             % Three optind for the voltage input: {IN1, IN1-IN2, IN1/<IN2>}
        PHOTODIODE_OPTIONS = {'normal input', 'reference input', 'both inputs', 'differential', 'referenced', 'laser block'};
        TYPE_OPTIONS = {'diffrential', 'normal', 'diff - both'};
    end

    methods
        function obj = PhotoDiodeDigitizerNiDaqControlled(name, niDaqGateChannel, digitizerChannelNames, digitizerChannelValues, defultInputConfiguration, ...
                                                          channelDelay, maxVoltage, inputImpedance, voltageOffset, nidaqChannelOffset, channelMinValue, channelMaxValue)
            % Contructor, creates the object and registers the channels in the DAQ.
            obj@Spcm(name);
            niDaqGateChannelName = sprintf('%s_gate', name);
            niDaqOffsetChannelName = sprintf('%s_offset', name);
            if ~isempty(nidaqChannelOffset)
                channelsNames = {niDaqGateChannelName, niDaqOffsetChannelName};
                channels = {niDaqGateChannel, nidaqChannelOffset};
            else
                channelsNames = {niDaqGateChannelName};
                channels = {niDaqGateChannel};
            end
            obj@NiDaqControlled(channelsNames, channels, num2cell(channelMinValue), num2cell(channelMaxValue));
            digi = Digitizer();
            
            obj.availableProperties.(obj.HAS_PHOTODIODE) = true;
            obj.availableProperties.(obj.HAS_DUAL_INPUT) = any(strcmp(digitizerChannelNames, 'input_2'));
            if ~any(strcmp(defultInputConfiguration, obj.INPUT_CONFIGURATIONS))
                obj.sendError('Digitizer input configuration can be only ''normal'', ''differential'', or ''referenced''!');
            end
            obj.defultInputConfiguration = defultInputConfiguration;
            obj.inputConfiguration = defultInputConfiguration;


            obj.niDaqGateChannelName = niDaqGateChannelName;
            obj.niDaqOffsetChannelName = niDaqOffsetChannelName;
            obj.digitizerChannelNames = digitizerChannelNames;
            obj.digitizerChannelValues = digitizerChannelValues;
            
            daq = getObjByName(NiDaq.NAME);
            obj.isEnabled = daq.readDigital(obj.niDaqGateChannelName);
            digi.maxVoltage = maxVoltage;
            digi.voltageOffset = voltageOffset;
            digi.inputImpedance = inputImpedance;
            if isempty(channelDelay)
                channelDelay = digitizerChannelValues * 0;
            end
            if any(mod(channelDelay * digi.MAX_SAMPLE_RATE, 1))
                obj.sendError('Digitizer channel delay durration must be a complete period of the maximum digitizer sample rate - %d MHz', digi.MAX_SAMPLE_RATE);
            end
            digi.channelDelay = channelDelay;
        end
    end
    
    methods % Integrator functions
        function setSPCMEnable(obj, newBooleanValue)
            % Enables/Disables the Integrator.
            daq = getObjByName(NiDaq.NAME);
            daq.writeDigital(obj.niDaqGateChannelName, newBooleanValue)
            if newBooleanValue && ~(~isempty(obj.isEnabled) && obj.isEnabled)
                obj.gateOnTic = tic;    % photodiode supply just switched on: scan pixels wait gateSettleTime
            end
            obj.isEnabled = newBooleanValue;
        end
        
    %%% Read by time %%%
        function prepareReadByTime(obj, integrationTimeInSec)
            % Prepare the Integrator to a scan by timer, with integration time of
            % integrationTime in seconds.
            digi = getObjByName(Digitizer.NAME);
            obj.voltageIntegrationTime = integrationTimeInSec * digi.UNIT_CONV;
            sampleRate = 1;     % 1 MHz should be enough
            digi.prepareAcquisition(obj.activeChannelsGUI, 1, obj.voltageIntegrationTime, 0, sampleRate);
        end
        
        function [voltage, sterr] = readFromTime(obj)
            % Reads the voltage for the integration time and returns a
            % single point which is the voltage and also the standard error.
            digi = getObjByName(Digitizer.NAME);
            digi.startExperiment;
            digi.forceAcquisition;
            digi.waitForDeviceReady;
            voltageFull = digi.readExperimentData(obj.activeChannelsGUI, true, true);
            [voltage, sterr] = obj.reduceSegments(voltageFull);     % one segment: scalar, or 2x1 for 'both inputs'
        end

        function [value, sterr] = reduceSegments(obj, voltageFull)
            % One value per acquired segment, for the photodiode option in
            % currentChannelGUI. Same arithmetic as readFromTime used before.
            % voltageFull: [active channels x segments x samples], with
            % the rows in the order of activeChannelsGUI.
            % value, sterr: [1 x segments], or [2 x segments] for option 3.
            % (Rows are indexed by position: the old code indexed them by
            % channel number, which failed for option 2, 'reference'.)
            nSeg = size(voltageFull, 2);
            row = @(r) reshape(voltageFull(r, :, :), nSeg, []);     % [segments x samples]
            switch obj.currentChannelGUI
                case {1, 2, 6}  % normal input, reference, laser block: one channel acquired
                    v = {row(1)};
                case 3          % both channels
                    v = {row(1), row(2)};
                case 4          % differential input
                    v = {row(1) - row(2)};
                case 5          % referenced input: signal / (reference normalised to its mean in the segment)
                    refSignal = row(2);
                    refSignal = refSignal ./ mean(refSignal, 2, "omitnan");
                    v = {row(1) ./ refSignal};
                otherwise
                    error('Unknown photodiode option %d', obj.currentChannelGUI);
            end
            value = zeros(numel(v), nSeg);
            sterr = zeros(numel(v), nSeg);
            for k = 1:numel(v)
                value(k, :) = mean(v{k}, 2, "omitnan")';
                sterr(k, :) = ste(v{k}, 0, 2)';
            end
        end
        
        function clearTimeRead(obj)
            % Clears the task for reading voltage by time.
            digi = getObjByName(Digitizer.NAME);
            digi.stopAcquisition;
        end
    %%% End (by time) %%%
    
    
    %%% By Stage %%%
    % Order of calls (StageScanner, unchanged):
    %   prepareCountByStage -> [startScanCount -> stage.ScanNextLine /
    %   stage.ScanX -> readFromScan] per line -> clearScanRead.
    % During ScanNextLine the stage calls acquireScanPixel(tPixel) once per
    % pixel, after it has arrived and settled.
        function prepareCountByStage(obj, stageName, nPixels, timeout, fastScan, pixelTime) %#ok<INUSD>
            % Prepare the photodiode for a scan by a stage: nPixels per
            % line. Called once per scan (or chunk). pixelTime (s) is only
            % passed by 1D scans; otherwise the stage gives it per pixel.
            if ~ValidationHelper.isValuePositiveInteger(nPixels)
                obj.sendError(sprintf('Can''t prepare for reading %d points, only positive integers allowed!', nPixels));
            end
            if ~any(strcmpi(obj.scanAcquisitionMode, {'trigger', 'software'}))
                obj.sendError(sprintf('scanAcquisitionMode must be ''trigger'' or ''software'', not ''%s''.', obj.scanAcquisitionMode));
            end
            if obj.currentChannelGUI == 3
                obj.sendError('The ''both inputs'' photodiode option gives two values per pixel and cannot be scanned. Choose another option.');
            end
            obj.scanningStageName = stageName;
            obj.nScanIntegration = nPixels;
            if exist('pixelTime', 'var') && ~isempty(pixelTime)
                obj.pixelTime = pixelTime;
            else
                obj.pixelTime = [];
            end
            obj.scanArmedPixelTime = [];
            obj.resetScanLine();
        end
        
        function startScanCount(obj)
            % Called before every line: forget the previous line's pixels
            % and disarm. The digitizer is armed at the line's first pixel,
            % after the stage has reached it, so no trigger timeout can
            % run out while the stage travels to the line start.
            obj.resetScanLine();
        end
        
        function acquireScanPixel(obj, tPixel)
            % Measure one pixel of the current line (tPixel in s). The
            % stage calls this when it is at the pixel and settled. Returns
            % when the measurement is complete, so the stage may move on.
            if isempty(obj.nScanIntegration) || obj.nScanIntegration < 1
                error('acquireScanPixel called before prepareCountByStage.');
            end
            if obj.scanPixelIndex >= obj.nScanIntegration
                error('Scan line has only %d pixels, but pixel %d was requested.', obj.nScanIntegration, obj.scanPixelIndex + 1);
            end
            if ~exist('tPixel', 'var') || isempty(tPixel)
                tPixel = obj.pixelTime;
            end
            if obj.scanPixelIndex > 0 && ~isequal(tPixel, obj.scanArmedPixelTime)
                error('Pixel time changed in the middle of a scan line.');
            end
            obj.waitForGate();
            
            i = obj.scanPixelIndex + 1;
            if strcmpi(obj.scanAcquisitionMode, 'trigger')
                if i == 1
                    obj.armScanLine(tPixel);
                end
                obj.fireScanTrigger(i);         % data stays on the digitizer until readFromScan
            else
                if ~obj.scanArmed || ~isequal(tPixel, obj.scanArmedPixelTime)
                    obj.prepareReadByTime(tPixel);
                    obj.scanArmed = true;
                    obj.scanArmedPixelTime = tPixel;
                end
                [v, e] = obj.readFromTime();
                obj.current1DScanVoltageValues(i) = v(1);
                obj.current1DScanSterrValues(i) = e(1);
            end
            obj.scanPixelIndex = i;
        end
        
        function [voltage, sterr] = readFromScan(obj)
            % Returns the current line, one value per pixel (row vectors).
            % Errors if the line is incomplete or data are missing, so the
            % scanner rescans the line.
            n = obj.nScanIntegration;
            if obj.scanPixelIndex ~= n
                measured = obj.scanPixelIndex;
                obj.resetScanLine();
                error('Scan line incomplete: %d of %d pixels were measured.', measured, n);
            end
            if strcmpi(obj.scanAcquisitionMode, 'trigger')
                digi = getObjByName(Digitizer.NAME);
                obj.checkTriggerGaps(digi);
                % Every pixel waited for its segment, so the acquisition
                % should be complete now. If the digitizer still waits for
                % a trigger, one was missed: discard the line now, rather
                % than let the digitizer trigger itself after its timeout
                % and shift the pixels.
                status = CsMl_QueryStatus(digi.handle);
                t = tic;
                while status == 2 && toc(t) < 1      % 2 = last segment still being acquired
                    status = CsMl_QueryStatus(digi.handle);
                end
                if status ~= 0
                    digi.stopAcquisition;
                    obj.scanArmed = false;
                    error('Digitizer still waits for a trigger after the last pixel (status %d): a pixel trigger was missed, line discarded.', status);
                end
                digi.waitForDeviceReady;        % returns at once; restarts the timer readExperimentData uses
                voltageFull = digi.readExperimentData(obj.activeChannelsGUI, true, true);
                obj.scanArmed = false;
                if size(voltageFull, 2) ~= n
                    error('Digitizer returned %d segments for %d pixels.', size(voltageFull, 2), n);
                end
                [v, e] = obj.reduceSegments(voltageFull);
                if any(isnan(v(1, :)))
                    error('No digitizer data for pixel(s) %s of this line (missed trigger?).', mat2str(find(isnan(v(1, :)))));
                end
                obj.current1DScanVoltageValues = v(1, :);
                obj.current1DScanSterrValues = e(1, :);
            end
            voltage = obj.current1DScanVoltageValues;
            sterr = obj.current1DScanSterrValues;
        end
        
        function clearScanRead(obj)
            % Clear the task that scans from stage.
            digi = getObjByName(Digitizer.NAME);
            digi.stopAcquisition;
            obj.scanArmed = false;
            obj.scanArmedPixelTime = [];
            if strcmpi(obj.scanAcquisitionMode, 'trigger')
                obj.setTrigger(false);          % in case a pixel was interrupted
            end
        end
    %%% End (By stage) %%%%
        
        
    %%% By PulseGenerator (Experiment) %%%
        function prepareExperimentCount(obj, nReads, timeout, nParemeters)
            % Prepare to read voltage from opening the detector window
            if ~ValidationHelper.isValuePositiveInteger(nReads)
                obj.sendError(sprintf('Can''t prepare for reading %d times, only positive integers allowed! Igonring.', nReads));
            end
            digi = getObjByName(Digitizer.NAME);
            obj.nExpIntegration = nReads;
%             if obj.fullDataAcquisition
%                 obj.nExpIntegration = nReads / (obj.segmentDuration*digi.MAX_SAMPLE_RATE);
%             end
            obj.expTimeoutTime = timeout;
            digi.prepareAcquisition(obj.activeChannelsMeas, obj.nExpIntegration, obj.segmentDuration)
            if obj.recordVoltageSpan
                obj.voltageHistogram = struct;
                obj.voltageHistogram.bins = -digi.maxVoltage(1):0.001:digi.maxVoltage(1);
                obj.voltageHistogram.values = zeros(1+obj.hasDualInput, length(obj.voltageHistogram.bins)-1);
            end
            if ~isempty(obj.daqVolatgeOffset)
                daq = getObjByName(NiDaq.NAME);
                daq.writeVoltage(obj.niDaqOffsetChannelName, obj.daqVolatgeOffset);
            end

        end
        
        function startExperimentCount(obj)
            % Actually start the process
            digi = getObjByName(Digitizer.NAME);
            digi.startExperiment;
        end
        
        function voltage = readFromExperiment(obj)
            digi = getObjByName(Digitizer.NAME);
            voltageFull = digi.readExperimentData(obj.activeChannelsMeas, true, true);
            
            if size(voltageFull, 2) ~= obj.nExpIntegration
                error('wrong number of points acquired')
            end
            
            if obj.recordVoltageSpan
                vFull = voltageFull(1,:,:);
                voltageHist = histcounts(vFull(:), obj.voltageHistogram.bins);
                obj.voltageHistogram.values(1,:) = obj.voltageHistogram.values(1,:) + voltageHist;
                if obj.hasDualInput
                    vFull = voltageFull(2,:,:);
                    voltageHist = histcounts(vFull(:), obj.voltageHistogram.bins);
                    obj.voltageHistogram.values(2,:) = obj.voltageHistogram.values(2,:) + voltageHist;
                end
            end

            input_1 = obj.digitizerChannelValues(strcmp(obj.digitizerChannelNames, 'input_1'));
            input_2 = obj.digitizerChannelValues(strcmp(obj.digitizerChannelNames, 'input_2'));
            switch obj.inputConfiguration
                case obj.INPUT_CONFIGURATIONS{1}                    % nomral
                    v = voltageFull(input_1, :, :);
                case obj.INPUT_CONFIGURATIONS{2}                    % differential
                    v = voltageFull(input_1, :, :) - voltageFull(input_2, :, :);
                case obj.INPUT_CONFIGURATIONS{3}                    % referenced
                    ref = voltageFull(input_2, :, :);
                    ref = mean(ref, 3, "omitnan");                  % mean over fast noise - we avarage it anyway
                    ref = ref ./ mean(ref, 2, "omitnan");           % normalization to the mean value
                    v = voltageFull(input_1, :, :) ./ ref;
            end

            v = squeeze(v);
            if ~obj.fullDataAcquisition
                voltage = mean(v, 2, "omitnan");
            else
                voltage = v;
            end
        end
        
        function stopExperimentCount(obj)
            digi = getObjByName(Digitizer.NAME);
            digi.stopAcquisition;
        end
        
        function clearExperimentRead(obj)
            if ~isempty(obj.daqVolatgeOffset)
                daq = getObjByName(NiDaq.NAME);
                daq.writeVoltage(obj.niDaqOffsetChannelName, 0);
            end
        end
    %%% End (by PulseGenerator) %%%
    end
    
    methods (Access = protected) % scan helpers
        function values = partialScanLine(obj)
            % The pixels of the current line measured so far, in the order they were
            % measured (the main window draws them live). Empty in 'trigger' mode, where
            % the line is only read from the digitizer at its end.
            values = [];
            if strcmpi(obj.scanAcquisitionMode, 'software') && obj.scanPixelIndex > 0 ...
                    && numel(obj.current1DScanVoltageValues) >= obj.scanPixelIndex
                values = obj.current1DScanVoltageValues(1:obj.scanPixelIndex);
            end
        end

        function resetScanLine(obj)
            % Empty line buffer, disarm the digitizer if it was armed.
            if obj.scanArmed && strcmpi(obj.scanAcquisitionMode, 'trigger')
                digi = getObjByName(Digitizer.NAME);
                digi.stopAcquisition;
                obj.scanArmed = false;
            end
            obj.scanPixelIndex = 0;
            obj.scanTriggerTimes = [];
            n = obj.nScanIntegration;
            if isempty(n); n = 0; end
            obj.current1DScanVoltageValues = nan(1, n);
            obj.current1DScanSterrValues = nan(1, n);
        end

        function armScanLine(obj, tPixel)
            % Trigger mode: configure the digitizer for one segment of
            % tPixel per pixel, external trigger, and arm it.
            digi = getObjByName(Digitizer.NAME);
            obj.setTrigger(false);                  % so the first pixel gives a rising edge
            tUs = tPixel * digi.UNIT_CONV;
            rate = obj.scanSampleRate(tUs, digi);
            digi.prepareAcquisition(obj.activeChannelsGUI, obj.nScanIntegration, tUs, 1, rate);
            % actual segment: rounded up to whole 32-sample blocks (+ channel delay)
            obj.scanSegmentSeconds = digi.segmentPostTrigSize / digi.sampleRate / digi.UNIT_CONV;
            digi.startExperiment;                   % armed: waits for nScanIntegration triggers
            obj.scanArmed = true;
            obj.scanArmedPixelTime = tPixel;
            obj.scanArmTic = tic;
            obj.scanTriggerTimes = [];
        end

        function fireScanTrigger(obj, i)
            % Trigger pixel i of the line with a short pulse: the rising
            % edge starts the pixel's segment, and the line goes low again
            % while that segment is still being recorded, so the digitizer
            % ignores that edge, as in the experiments. (On Setup 7 the
            % low-going edge starts an extra segment if it comes after the
            % segment has ended: scanDebug, 2026-11-23.) Returns once the
            % segment is complete. Errors, so the line is rescanned, if
            % the digitizer finished early or the pulse may have outlasted
            % the segment.
            tCall = toc(obj.scanArmTic);
            obj.setTrigger(true);
            tOn = toc(obj.scanArmTic);          % the rising edge happened before this
            obj.setTrigger(false);
            tOff = toc(obj.scanArmTic);         % ... and the falling edge before this
            obj.scanTriggerTimes(end+1) = tOn;
            obj.waitSeconds(obj.scanSegmentSeconds + obj.scanTriggerMargin - (tOff - tOn));

            digi = getObjByName(Digitizer.NAME);
            n = obj.nScanIntegration;
            if i < n && CsMl_QueryStatus(digi.handle) == 0
                digi.stopAcquisition;
                obj.scanArmed = false;
                error('Digitizer finished after pixel %d of %d (it took extra segments): line discarded.', i, n);
            end
            if tOff - tCall >= obj.scanSegmentSeconds
                digi.stopAcquisition;
                obj.scanArmed = false;
                error(['Trigger pulse for pixel %d took up to %.1f ms, longer than the %.1f ms segment, so the digitizer ', ...
                    'may have taken an extra segment: line discarded. Use scanAcquisitionMode = ''software'' or a longer pixel time.'], ...
                    i, 1e3 * (tOff - tCall), 1e3 * obj.scanSegmentSeconds);
            end
        end

        function setTrigger(obj, value)
            pg = getObjByName(PulseGenerator.NAME);
            if isempty(pg)
                error('Scan trigger mode needs the pulse generator (channel ''%s''). Set scanAcquisitionMode = ''software''.', obj.scanTriggerChannel);
            end
            if value
                pg.On(obj.scanTriggerChannel);
            else
                pg.Off(obj.scanTriggerChannel);
            end
        end

        function checkTriggerGaps(obj, digi)
            % The digitizer triggers itself after acqInfo.TriggerTimeout
            % (us) without a trigger, which would shift the pixels. Refuse
            % the line if any wait came close to that.
            limit = 0.8 * double(digi.acqInfo.TriggerTimeout) / 1e6;
            gaps = diff([0, obj.scanTriggerTimes]);
            if limit > 0 && any(gaps > limit)
                error('%.1f s between two pixel triggers (digitizer timeout %.0f s): line discarded.', max(gaps), limit / 0.8);
            end
        end

        function waitForGate(obj)
            % After the photodiode supply is switched on, give it
            % gateSettleTime before the first pixel.
            if ~isempty(obj.gateOnTic)
                remaining = obj.gateSettleTime - toc(obj.gateOnTic);
                if remaining > 0
                    pause(remaining);
                end
                if toc(obj.gateOnTic) >= obj.gateSettleTime
                    obj.gateOnTic = [];
                end
            end
        end

        function rate = scanSampleRate(obj, tUs, digi)
            % Highest rate <= 1 MHz (the read-by-time rate) that keeps a
            % pixel below scanMaxSamplesPerPixel samples, preferring rates
            % that give a whole number of samples.
            rates = digi.AVALIABLE_SAMPLE_RATES(digi.AVALIABLE_SAMPLE_RATES <= 1);    % MHz, descending
            ok = rates(tUs .* rates <= obj.scanMaxSamplesPerPixel);
            if isempty(ok)
                ok = rates(end);
            end
            whole = abs(tUs .* ok - round(tUs .* ok)) < 1e-6 & round(tUs .* ok) >= 1;
            if any(whole)
                rate = ok(find(whole, 1));
            else
                rate = ok(1);
            end
        end
    end

    methods (Static, Access = protected)
        function waitSeconds(seconds)
            % pause() for most of it, then spin: never returns early.
            t = tic;
            if seconds > 0.005
                pause(seconds - 0.003);
            end
            while toc(t) < seconds
            end
        end
    end

    methods
        function onNiDaqReset(obj, niDaq) %#ok<INUSD> NiDaqControlled.onEvent passes the NiDaq
            % This function jumps when the NiDaq resets
%             if obj.nScanIntegration > 0
%                 prepareMeasureByStageInternal(obj, niDaq);
%             end
            if obj.isEnabled
                % When reset, the NiDaq no longer remembers whether the
                % channel was off or on. We need to set the value, but only
                % if needed (since writeDigital is costly)
                obj.setSPCMEnable(obj.isEnabled)
            end
        end
        
        function channels = get.activeChannelsGUI(obj)
            signalChan = find(strcmp({'input_1'}, obj.digitizerChannelNames));
            referenceChan = find(strcmp({'input_2'}, obj.digitizerChannelNames));
            laserBlockChan = find(strcmp({'laser_block'}, obj.digitizerChannelNames));
            
            switch obj.currentChannelGUI
                case 1
                    channels = signalChan;
                case 2
                    channels = referenceChan;
                case {3, 4, 5}
                    channels = [signalChan, referenceChan];
                case 6
                    channels = laserBlockChan;
            end
        end

        function channels = get.activeChannelsMeas(obj)
            channels = find(strcmp({'input_1'}, obj.digitizerChannelNames));
            if obj.hasDualInput
                channels = [channels, find(strcmp({'input_2'}, obj.digitizerChannelNames))];
            end
        end

        
    end
    
    methods (Static)
        function finalValue = roundOrError(value, digits, epsilon, errorTxt)
            % check wheather the value is not rounded. If the error larger
            % than epsilon its return error, otherwise it rounded
            if ~prod((value - round(value, digits)) == 0)
                if max(abs((value - round(value, digits)))) > epsilon
                    error(errorTxt)
                else
                    finalValue = round(value, digits);
                end
            else
                finalValue = value;
            end
        end
        
        function PhotoDiodeObj = create(PhotoDiodeName, PhotoDiodeStruct)
            missingField = FactoryHelper.usualChecks(PhotoDiodeStruct, PhotoDiodeDigitizerNiDaqControlled.NEEDED_FIELDS);
            if ~isnan(missingField)
                EventStation.anonymousError(...
                    'Can''t initialize NiDaq-controlled PhotoDiode - required field "%s" was not found in initialization struct!', ...
                    missingField);
            end
            
            % We want to get either values set in json, or empty variables
            % (which will be handled by NiDaqControlled constructor):
            PhotoDiodeStruct = FactoryHelper.supplementStruct(PhotoDiodeStruct, PhotoDiodeDigitizerNiDaqControlled.OPTIONAL_FIELDS);
            
            gate = PhotoDiodeStruct.nidaq_channel_gate;
            digitizerChannelNames = PhotoDiodeStruct.digitizer_channel_names;
            digitizerChannelValues = PhotoDiodeStruct.digitizer_channel_values;
            maxVoltage = PhotoDiodeStruct.max_voltage;
            inputImpedance = PhotoDiodeStruct.input_impedance;
            minVal = PhotoDiodeStruct.nidaq_channel_min_val;
            maxVal = PhotoDiodeStruct.nidaq_channel_max_val;
            voltageOffset = PhotoDiodeStruct.voltage_offset;
            defultInputConfiguration = PhotoDiodeStruct.defultd_input_configuration;
            nidaqChannelOffset = PhotoDiodeStruct.nidaq_channel_offset;
            channelDelay = PhotoDiodeStruct.channel_delay;
            PhotoDiodeObj = PhotoDiodeDigitizerNiDaqControlled(PhotoDiodeName, gate, digitizerChannelNames, digitizerChannelValues, ...
                        defultInputConfiguration, channelDelay, maxVoltage, inputImpedance, voltageOffset, nidaqChannelOffset, minVal, maxVal);
        end
    end
    
    %% Overridden from spcm
    methods (Static)
        % Auxilary function, for parallel reading: we need to fetch objects
        % before sending task to workers
        function measuringObj = variablesForTimeRead
            measuringObj = getObjByName(NiDaq.NAME);
        end
    end

end