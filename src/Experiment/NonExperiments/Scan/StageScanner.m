classdef StageScanner < EventSender & EventListener & Savable
    %STAGESCANNER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mScan
        mScanExtraInfo % Hist, if available.
        mStageName
        mStageScanParams
        mCurrentlyScanning = false
        mProgress = struct('done', 0, 'inLine', 0, 'total', 0, 'tic', [], 'lastSent', [], ...
            'line', 0, 'nLines', 0, 'slice', 0, 'nSlices', 0, 'pass', 0, 'dims', 0)
    end

    properties (Access = private)
        mVolume = []    % 3D scan in progress: struct(vol, slice) -- see scan3D / publishScan
        mHalted = false % the stages were halted during this scan (then no return to start)
        mLive = []      % the line being scanned, for drawing it pixel by pixel -- see livePreview
    end

    properties (Constant)
        NAME = 'stageScanner'

        EVENT_SCAN_UPDATED = 'scanUpdated'
        EVENT_SCAN_PROGRESS = 'scanProgress'    % extraInfo.scanProgress = struct (see progressInfo)
        PROGRESS_MIN_INTERVAL = 0.2             % s between progress events
        LIVE_MIN_INTERVAL = 0.1                 % s between pixel-by-pixel redraws of the image
        EVENT_SCAN_STARTED = 'scanStarted'
        EVENT_SCAN_FINISHED = 'scanFinished'
        EVENT_SCAN_STOPPED_MANUALLY = 'scanStoppedManually'
        % events can be sent with both EVENT_SCAN_FINISHED and EVENT_SCAN_STOPPED_MANUALLY
        % such events will always have EVENT_SCAN_FINISHED = true,
        %                and will have EVENT_SCAN_STOPPED_MANUALLY = (true or false)
        
        PROPERTY_SCAN_PARAMS = 'scanParameters'
        
        TRIALS_AMOUNT_ON_ERROR = 4;
    end
    
    methods (Access = private)
        function obj = StageScanner
            obj@EventSender(StageScanner.NAME);
            obj@Savable(StageScanner.NAME);
            obj@EventListener({SaveLoadCatImage.NAME, StageControlEvents.NAME});     % (stage control: Halt)
            obj.clear();
        end
    end
    
    methods
        function sendEventScanFinished(obj)
            obj.sendEvent(struct( ...
                obj.EVENT_SCAN_FINISHED, true, ...
                obj.EVENT_SCAN_STOPPED_MANUALLY, false, ...
                obj.PROPERTY_SCAN_PARAMS, obj.mStageScanParams));
        end
        function sendEventScanStopped(obj)
            obj.sendEvent(struct(...
                obj.EVENT_SCAN_FINISHED, true, ...
                obj.EVENT_SCAN_STOPPED_MANUALLY, true, ...
                obj.PROPERTY_SCAN_PARAMS, obj.mStageScanParams));
        end
        function sendEventScanStarting(obj)
            obj.sendEvent(struct(obj.EVENT_SCAN_STARTED, true));
        end
        function sendEventScanUpdated(obj, scanResults, scanExtraInfo)
            params = obj.mStageScanParams;
            if params.nDims == 3
                % A volume: the image axes (in xyz order), then the slice axis
                [imgAxes, sliceAxis] = StageScanner.volumeAxes(params);
                phAxes = {params.getScanAxisVector(imgAxes(1)), params.getScanAxisVector(imgAxes(2)), ...
                    params.getScanAxisVector(sliceAxis)};
                scanAxes = ClassStage.SCAN_AXES([imgAxes, sliceAxis]);
            else
                axis1 = params.getFirstScanAxisVector;
                axis2 = params.getSecondScanAxisVector;
                phAxes = {axis1, axis2};
                scanAxes = params.getScanAxes;
            end
            stageName = obj.mStageName;
            botLabel = obj.getBottomScanLabel;
            leftLabel = obj.getLeftScanLabel;
            if ~exist('scanExtraInfo', 'var')
                scanExtraInfo = struct();
            end
            extra = EventExtraScanUpdated(scanResults, phAxes, scanAxes, stageName, botLabel, leftLabel, scanExtraInfo);
            obj.sendEvent(struct(obj.EVENT_SCAN_UPDATED, extra));
        end
        
        function startScan(obj, advancedScanParams)
            if isnan(obj.mStageName)
                obj.sendError(['Can''t start scan: unknown stage name!' ...
                    'Please call stageScanner.switchTo(your_stage_name)' ...
                    'before calling stageScanner.startScan(). Exiting']);
            end
            
            stage = getObjByName(obj.mStageName);
            if isempty(stage); throwBaseObjException(obj.mStageName); end
            if ~exist('advancedScanParams', 'var')
                obj.mStageScanParams = stage.scanParams.copy;
                obj.mStageScanParams.scanBlocks = {};
            else
                obj.mStageScanParams = advancedScanParams;
            end
            isFastScan = obj.mStageScanParams.fastScan;
            stage.FastScan(isFastScan);

            % Same direction on every line, or zigzag -- for stages that can do both
            if isprop(stage, 'scanBidirectional')
                stage.scanBidirectional = logical(obj.mStageScanParams.zigzag);
            end

            % Where the stage is now: it goes back here when the scan ends (returnToStart)
            startPos = stage.Pos(stage.availableAxes);
            obj.mHalted = false;

            % From now on, changes in the GUI (from\to\fixed\numPoints)
            % won't affect the StageScanParmas object stored here.

            obj.mCurrentlyScanning = true;
            obj.mVolume = [];
            obj.mLive = [];
            obj.progressStart(obj.mStageScanParams, 1);
            obj.sendEventScanStarting();
            
            spcm = getObjByName(Spcm.NAME);
            if isempty(spcm); throwBaseObjException(Spcm.NAME); end
            spcm.setSPCMEnable(true);
            
            
            try
                stage.sanityCheckForScanRange(obj.mStageScanParams);
            catch err
                obj.mCurrentlyScanning = false;
                obj.sendEventScanStopped;
                rethrow(err)
            end
            
            if obj.mStageScanParams.isMWContrastScan
                %%% Set Frequency Generator
                fg = getObjByName(FrequencyGenerator.getDefaultFgName); 
                fg.amplitude = obj.mStageScanParams.MWAmplitude;
                fg.frequency = obj.mStageScanParams.MWFrequency;
                fg.output = 1;
            end
            timerVal = tic;
            disp('Initiating scan...');
            try
                [kcpsScanMatrix, extraInfo] = obj.scan(stage, spcm, obj.mStageScanParams);
                % kcps = kilo counts per second

                pass = 1;
                while (stage.scanParams.continuous && obj.mCurrentlyScanning)
                    pass = pass + 1;
                    obj.progressStart(obj.mStageScanParams, pass);
                    [kcpsScanMatrix, extraInfo] = obj.scan(stage, spcm, obj.mStageScanParams, kcpsScanMatrix, extraInfo);
                    drawnow; % todo check what happens if removing this line
                end
            catch err
                % We couldn't scan. Wrap it up nicely
                spcm.setSPCMEnable(false);
                obj.mCurrentlyScanning = false;
                obj.mVolume = [];
                obj.mLive = [];
                stage.sendEventStageAvailabilityChanged;
                if obj.mStageScanParams.isMWContrastScan
                    fg = getObjByName(FrequencyGenerator.getDefaultFgName);
                    fg.output = 0;
                end
                try
                    obj.returnToStart(stage, startPos);
                catch
                end

                rethrow(err)
            end
            if obj.mStageScanParams.isMWContrastScan
                fg = getObjByName(FrequencyGenerator.getDefaultFgName);
                fg.output = 0;
            end
            
            % Maybe we didn't encounter an error, but the scan still did not happen
            if isempty(kcpsScanMatrix)
                spcm.setSPCMEnable(false);
                obj.mCurrentlyScanning = false;
                obj.returnToStart(stage, startPos);
                stage.sendEventStageAvailabilityChanged;
                return
            end
            
            toc(timerVal)
            
            spcm.setSPCMEnable(false);
            obj.mScan = kcpsScanMatrix;
            obj.mScanExtraInfo = extraInfo;
            obj.mVolume = [];
            obj.mLive = [];

            scanStoppedManually = ~obj.mCurrentlyScanning;  % maybe someone has changed this boolean meanwhile
            obj.mCurrentlyScanning = false;                 % So when events are sent, they will know we are done.
            obj.returnToStart(stage, startPos);             % (before "finished": then the stage is back)
            obj.progressSend(true);
            
            if scanStoppedManually
                obj.sendEventScanStopped();
            else
                obj.sendEventScanFinished();
            end
            % (At least) two things should happen by this event:
            % 1. ImageScanResult will update
            % 2. SaveLoad will get the new scan, save it into local
            %    struct, and (if needed) will autosave it.
            
            stage.sendEventPositionChanged;
        end
        
        function [kcpsScanMatrix, extraInfo] = scan(obj, stage, spcm, scanParams, kcpsScanMatrixOptional, optionalExtraInfo)
            % Scan the stage.
            % stage - an object deriving from ClassStage
            % spcm - an object deriving from Spcm
            % scanParams - the scan parameters. an object deriving from StageScanParams
            % kcpsScanMatrixOptional - the last scanned matrix, if exists
            skipCreatingMatrix = exist('kcpsScanMatrixOptional', 'var');
            
            kcpsScanMatrix = [];        % Return value for nonvalid dinemsion number
            extraInfo = struct();
            nDimensions = sum(~scanParams.isFixed);
            switch nDimensions
                case 0
                    EventStation.anonymousWarning('Nothing to scan, dude!');
                    return;
                case 1
                    % One dimensional scanning
                    [kcpsScanMatrix, extraInfo] = obj.scan1D(stage, spcm, scanParams);
                case 2
                    % 2D scan
                    if skipCreatingMatrix
                        [kcpsScanMatrix, extraInfo] = obj.scan2D(stage, spcm, scanParams, kcpsScanMatrixOptional, optionalExtraInfo);
                    else
                        [kcpsScanMatrix, extraInfo] = obj.scan2D(stage, spcm, scanParams);
                    end
                case 3
                    % 3D scan: one 2D scan per slice
                    if skipCreatingMatrix
                        [kcpsScanMatrix, extraInfo] = obj.scan3D(stage, spcm, scanParams, kcpsScanMatrixOptional);
                    else
                        [kcpsScanMatrix, extraInfo] = obj.scan3D(stage, spcm, scanParams);
                    end
                otherwise
                    EventStation.anonymousWarning('%d-dimensional scan requested. String thoery is not yet implemented!\n', nDimensions);
            end
        end
        
        function [kcpsValue, steValue] = scanPoint(obj, stage, spcm, scanParams)
            % scan "0D"
            % input arg's:
            %    stage - object deriving from ClassStage
            %    spcm - object deriving from Spcm
            %    scanParams - object deriving from StageScanParams
            % output:
            %    kcpsValue - double. Value read by the spcm
            
            % To be used by TrackablePosition. Should NOT be used for
            % actual scanning.
            
            %%%% move to location %%%%
            pos = scanParams.fixedPos;
            phAxes = stage.availableAxes;
            stage.move(phAxes, pos);
            
            %%%% try to scan %%%%
            scanOk = false;
            for trial = 1:StageScanner.TRIALS_AMOUNT_ON_ERROR
                try
                    spcm.prepareReadByTime(scanParams.pixelTime);
                    [kcps, kcps_ste] = spcm.readFromTime();
                    spcm.clearTimeRead;
                    if kcps == 0
                        obj.sendError('No Signal Detected!')
                    end
                    scanOk = true;
                    break;
                catch err
                    %                    rethrow(err) % Uncomment to debug
                    warning(err.message);
                    fprintf('Reading from SPCM failed at trial %d, attempting to rescan.\n', trial);
                end
            end
            
            % This line will be reached when break()ing out of the for,
            % or when all trials went without success
            if ~scanOk
                obj.sendError(sprintf('Reading from SPCM failed after %d trials', StageScanner.TRIALS_AMOUNT_ON_ERROR));
            end
            
            kcpsValue = kcps;
            steValue = kcps_ste;
        end
        
        function [kcpsScanVector, extraInfo] = scan1D(obj, stage, spcm, scanParams)
            % scan 1D
            % stage - an object deriving from ClassStage
            % spcm - an object deriving from Spcm
            % scanParams - the scan parameters. an object deriving from StageScanParams
            % returns - a vector of the scan
            
            % If the area to scan is bigger than the maximum scanning-area
            % of the stage, divide it to smallers chunks to scan
            
            % ~~~~ preparing variables: ~~~~
            isFastScan = scanParams.fastScan;
            nFlat = 0;      % A flat section at the start of ramp. parameter not needed genrally, a stage can overwrite if needed
            nOverRun = 0;   % Let the waveform over run the start and end. Not needed genrally, a stage can overwrite if needed
            tPixel = scanParams.pixelTime;
            maxScanSize = stage.ReturnMaxScanSize(1);
            nPixels = length(scanParams.getFirstScanAxisVector());
            isMWContrastScan = scanParams.isMWContrastScan;
            kcpsScanVector = zeros(nPixels, 1+isMWContrastScan);
            axisToScan = scanParams.getScanAxes; % string of size 1
            % vector between [min, max] (in the order it is scanned) or the fixed position if exist
            x = scanParams.getScanAxisVectorInScanOrder(1);
            y = scanParams.getScanAxisVectorInScanOrder(2);
            z = scanParams.getScanAxisVectorInScanOrder(3);
            scannedAxis = find(~scanParams.isFixed, 1);
            isReversed = numel(scanParams.reverse) >= scannedAxis && scanParams.reverse(scannedAxis);
            
            bHasLiftime = spcm.hasLifetime();
            extraInfo = struct();
            extraInfo.sterr = nan(nPixels, 1+isMWContrastScan);  % per-point standard error, if the detector gives one
            spcm = getObjByName(Spcm.NAME);
            if bHasLiftime && ~isMWContrastScan
                extraInfo.hist = zeros(nPixels, spcm.nBins);
            end
            if isMWContrastScan
                pg = getObjByName(PulseGenerator.NAME);
            end
            
            if ~obj.mCurrentlyScanning
                return
            end
            
            % ~~~~ checks on size of scan ~~~~
            if (nPixels > maxScanSize)
                fprintf('Max number of points to scan is %d, %d were requested for %s axis. Scanning in parts.\n', maxScanSize, nPixels, axisToScan);
            end
            
            % ~~~~ continue scanning until finishing all the scan ~~~~
            
            vectorStartIndex = 1;
            while (vectorStartIndex <= nPixels)
                % every iteration in the "while" scans one chunk
                pixelsLeftToScan = nPixels - vectorStartIndex + 1;
                curPixelsAmountToScan = min(pixelsLeftToScan, maxScanSize);
                vectorEndIndex = vectorStartIndex + curPixelsAmountToScan - 1;
                nPoints = curPixelsAmountToScan + 2*(nFlat + nOverRun);
                timeout = 2*nPoints*tPixel;
                scanOk = false;
                
                % Prepare Scan
                
                prepareScanfuncName = sprintf('PrepareScan%s', upper(axisToScan));
                feval(prepareScanfuncName, stage, x, y, z, nFlat, nOverRun, tPixel);
                
                
                % try to scan
                kcps = zeros(curPixelsAmountToScan , size(kcpsScanVector, 2));
                kcpsErr = nan(size(kcps));
                for i = 1:size(kcpsScanVector, 2) % MW Contrast
                    for trial = 1:StageScanner.TRIALS_AMOUNT_ON_ERROR
                        if i==2
                            pg.on('MW');
                        end
                        try
                            if ~obj.mCurrentlyScanning
                                stage.AbortScan();
                                spcm.clearScanRead();  % todo - added to try resolving the problem. Wasn't here in the first place!
                                return
                            end

                            if isa(spcm, "CameraControlled") || isa(spcm, "PhotoDiodeDigitizerNiDaqControlled")
                                spcm.prepareCountByStage(stage.name, curPixelsAmountToScan, timeout, isFastScan, tPixel);
                            else
                                spcm.prepareCountByStage(stage.name, curPixelsAmountToScan, timeout, isFastScan);
                            end
                            spcm.startScanCount();
                            
                            % scan stage (drawing the points as they come: unmeasured = NaN)
                            liveVector = kcpsScanVector;
                            liveVector(vectorStartIndex:vectorEndIndex, i) = NaN;
                            obj.liveStart(liveVector, extraInfo, '1D', 1, vectorStartIndex, vectorEndIndex, ...
                                false, isReversed, i, true);
                            scanfuncName = sprintf('Scan%s', upper(axisToScan));
                            feval(scanfuncName, stage, x,y,z, nFlat, nOverRun, tPixel);
                            obj.liveStop();
                            
                            % read counter
                            [lineValues, lineErr] = StageScanner.readScanLine(spcm);
                            kcps(:,i) = lineValues;
                            kcpsErr(:,i) = lineErr;
                            if bHasLiftime && ~isMWContrastScan
                                hist = spcm.lastScanHist;
                            end
                            spcm.clearScanRead();
                            if ~nnz(kcps(:,i))
                                obj.sendError('No Signal Detected!')
                            end
                            
                            scanOk = true;
                            break;
                        catch err
                            obj.liveStop();
%                             rethrow(err) % Uncomment to debug
                            warning(err.message);
                            
                            if ~obj.mCurrentlyScanning
                                stage.AbortScan();
                                return
                            end
                            fprintf('Scan failed at trial %d, attempting to rescan.\n', trial);
                        end
                    end
                    if i==2
                        pg.off('MW');
                    end
                end % MW Contrast
                
                % This line will be reached when break()ing out of the for,
                % or when all trials went without success
                stage.AbortScan();
                if ~scanOk; obj.sendError(sprintf('Scan failed after %d trials', StageScanner.TRIALS_AMOUNT_ON_ERROR));end
                
                % A reversed axis was scanned to -> from: store it from -> to
                if isReversed
                    kcps = flipud(kcps);
                    kcpsErr = flipud(kcpsErr);
                    if bHasLiftime && ~isMWContrastScan
                        hist = flipud(hist);
                    end
                end

                % Update the scan results in the returned vector
                if isMWContrastScan % Contrast
                    kcpsScanVector(vectorStartIndex:vectorEndIndex, :) = kcps;
                else
                    kcpsScanVector(vectorStartIndex:vectorEndIndex) = kcps;
                    if bHasLiftime
                        extraInfo.hist(vectorStartIndex:vectorEndIndex, :) = hist;
                    end
                end
                extraInfo.sterr(vectorStartIndex:vectorEndIndex, :) = kcpsErr;
                obj.progressLineDone(curPixelsAmountToScan * size(kcpsScanVector, 2));

                % Go tell everybody
                obj.sendEventScanUpdated(kcpsScanVector, extraInfo);
                
                % Prepare the next scan chunk
                vectorStartIndex = vectorEndIndex + 1;
            end
        end
        
        function [kcpsScanMatrix, extraInfo] = scan2D(obj, stage, spcm, scanParams, optionalKcpsScanMatrix, optionalExtraInfo)
            % Scans a 2D image.
            % stage - an object deriving from ClassStage
            % spcm - an object deriving from Spcm
            % scanParams - the scan parameters. an object deriving from StageScanParams
            % optionalKcpsScanMatrix - if exists, the previous scan results
            % returns - a matrix of the scan
            
            
            % Method:
            % First, check what would be the best way to scan (minimize
            % amount of scans needed)
            % Then, call scan2dChunk for each chunk to get scan results and combine them together
            
            % The movement between pixels in each line is called "axis a",
            % while movement between lines is called "axis b"
            
            if ~obj.mCurrentlyScanning
                return
            end
            
            isMWContrastScan = scanParams.isMWContrastScan;
            spcm = getObjByName(Spcm.NAME);
            
            if isempty(scanParams.scanBlocks) % Normal scan, look for chucks due to stages limits
                % ~~~~ preparing variables: ~~~~
                maxScanSize = stage.ReturnMaxScanSize(2);
                [axisAIndex, axisBIndex] = StageScanner.lineAndStepAxes(maxScanSize, scanParams);
                axisCIndex = setdiff(ClassStage.getAxis(ClassStage.SCAN_AXES), [axisAIndex, axisBIndex]);
                % axis c is where the zero point is not moving throgh the scan
                reverseA = numel(scanParams.reverse) >= axisAIndex && scanParams.reverse(axisAIndex);
                reverseB = numel(scanParams.reverse) >= axisBIndex && scanParams.reverse(axisBIndex);
                isFlipped = axisAIndex > axisBIndex;  % (for example, if axis y is before axis x, than "isFlipped" == true)
                isFastScan = scanParams.fastScan;    % a boolean
                tPixel = scanParams.pixelTime;   % time for each pixel
                numPointsAxisA = scanParams.numPoints(axisAIndex);
                numLinesAxisB = scanParams.numPoints(axisBIndex);
                vectorAxisA = scanParams.getScanAxisVector(axisAIndex);
                vectorAxisB = scanParams.getScanAxisVector(axisBIndex);
                pointAxisC = scanParams.getScanAxisVector(axisCIndex);
                chunksPerLine = ceil(numPointsAxisA / maxScanSize);
                
                bHasLiftime = spcm.hasLifetime();
                
                if exist('optionalKcpsScanMatrix', 'var')
                    kcpsScanMatrix = optionalKcpsScanMatrix;
                    extraInfo = optionalExtraInfo;
                else
                    extraInfo = struct();
                    if isFlipped
                        kcpsScanMatrix = zeros(numPointsAxisA, numLinesAxisB, 1+scanParams.isMWContrastScan);
                        if bHasLiftime && ~isMWContrastScan
                            extraInfo.hist = zeros(numPointsAxisA, numLinesAxisB, spcm.nBins);
                        end
                    else
                        kcpsScanMatrix = zeros(numLinesAxisB, numPointsAxisA, 1+scanParams.isMWContrastScan);
                        if bHasLiftime && ~isMWContrastScan
                            extraInfo.hist = zeros(numLinesAxisB, numPointsAxisA, spcm.nBins);
                        end
                    end
                end
                if ~isstruct(extraInfo); extraInfo = struct(); end
                if ~isfield(extraInfo, 'sterr') || ~isequal(size(extraInfo.sterr), size(kcpsScanMatrix))
                    extraInfo.sterr = nan(size(kcpsScanMatrix));   % per-pixel standard error, if the detector gives one
                end
                obj.mScan = kcpsScanMatrix;
                obj.mScanExtraInfo = extraInfo;
                
                %%%%%% iterate through the chunks and scan each chunk %%%%%%%
                startIndexAxisA = 1;
                for chunkIndex = 1 : chunksPerLine
                    if ~obj.mCurrentlyScanning
                        break;
                    end
                    
                    pixelsLeftAxisA = numPointsAxisA - startIndexAxisA + 1;
                    chunkSizeAxisA = min(maxScanSize, pixelsLeftAxisA);
                    endIndexAxisA = startIndexAxisA + chunkSizeAxisA - 1;
                    chunkAxisA = vectorAxisA(startIndexAxisA:endIndexAxisA);
                    chunkAxisB = vectorAxisB;  % no restriction on number of lines!
                    lineIndexMap = 1:length(chunkAxisB);
                    % A reversed axis is visited to -> from; the matrix keeps from -> to
                    if reverseA; chunkAxisA = fliplr(chunkAxisA); end
                    if reverseB
                        chunkAxisB = fliplr(chunkAxisB);
                        lineIndexMap = fliplr(lineIndexMap);
                    end
                    [kcpsScanMatrix, extraInfo] = obj.scan2dChunk(...
                        kcpsScanMatrix, spcm, stage, tPixel, chunkAxisA, chunkAxisB, ...
                        axisAIndex, axisBIndex, pointAxisC, isFastScan, isFlipped, ...
                        1, length(chunkAxisB), ...
                        startIndexAxisA, endIndexAxisA, ...
                        extraInfo, lineIndexMap, reverseA);
                    obj.mScan = kcpsScanMatrix;
                    obj.mScanExtraInfo = extraInfo;
                    obj.publishScan(kcpsScanMatrix, extraInfo);
                    
                    startIndexAxisA = endIndexAxisA + 1;     % Updating for next chunk
                end
            else % Advanced scan, look for chucks defined in scan params
                % ~~~~ preparing variables: ~~~~
                maxScanSize = stage.ReturnMaxScanSize(2);
                isFastScan = scanParams.fastScan;    % a boolean
                tPixel = scanParams.pixelTime;   % time for each pixel
                [axisAIndex, axisBIndex] = StageScanner.optimize2dScanDirections(maxScanSize, scanParams.copy);
                isFlipped = axisAIndex > axisBIndex;  % (for example, if axis y is before axis x, than "isFlipped" == true)
                fullVectorAxisA = scanParams.getScanAxisVector(axisAIndex);
                fullVectorAxisB = scanParams.getScanAxisVector(axisBIndex);
                
                bHasLiftime = spcm.hasLifetime();
                
                if exist('optionalKcpsScanMatrix', 'var')
                    kcpsScanMatrix = optionalKcpsScanMatrix;
                    extraInfo = optionalExtraInfo;
                else
                    extraInfo = struct();
                    if isFlipped
                        kcpsScanMatrix = zeros(length(fullVectorAxisA), length(fullVectorAxisB), 1+scanParams.isMWContrastScan);
                        if bHasLiftime && ~isMWContrastScan
                            extraInfo.hist = zeros(length(fullVectorAxisA), length(fullVectorAxisB), spcm.nBins);
                        end
                    else
                        kcpsScanMatrix = zeros(length(fullVectorAxisB), length(fullVectorAxisA), 1+scanParams.isMWContrastScan);
                        if bHasLiftime && ~isMWContrastScan
                            extraInfo.hist = zeros(length(fullVectorAxisB), length(fullVectorAxisA), spcm.nBins);
                        end
                    end
                end
                obj.mScan = kcpsScanMatrix;
                obj.mScanExtraInfo = extraInfo;
                
                %%%%%% iterate through the chunks and scan each chunk %%%%%%%
                for i = 1:length(scanParams.scanBlocks)
                    if ~obj.mCurrentlyScanning
                        break;
                    end
                    if isFlipped
                        scanBlock = scanParams.scanBlocks{i}([3:4, 1:2, 5]);
                    else
                        scanBlock = scanParams.scanBlocks{i};
                    end
                    fromA = scanBlock(1);
                    toA = scanBlock(2);
                    fromB = scanBlock(3);
                    toB = scanBlock(4);
                    pointAxisC = scanBlock(5);
                    
                    indicesToScanA = find(fullVectorAxisA >= fromA & fullVectorAxisA <= toA);
                    indicesToScanB = find(fullVectorAxisB >= fromB & fullVectorAxisB <= toB);
                    vectorAxisA = fullVectorAxisA(indicesToScanA);
                    vectorAxisB = fullVectorAxisB(indicesToScanB);
                    
                    [kcpsScanMatrix, extraInfo] = obj.scan2dChunk(...
                        kcpsScanMatrix, spcm, stage, tPixel, vectorAxisA, vectorAxisB, ...
                        axisAIndex, axisBIndex, pointAxisC, isFastScan, isFlipped, ...
                        indicesToScanB(1), indicesToScanB(end), ...
                        indicesToScanA(1), indicesToScanA(end), ...
                        extraInfo);
                    obj.mScan = kcpsScanMatrix;
                    obj.mScanExtraInfo = extraInfo;
                    obj.publishScan(kcpsScanMatrix, extraInfo);
                end
            end
        end
        
        function [kcpsMatrix, extraInfo] = scan2dChunk(...
                obj, ...  StageScanner object
                kcpsMatrix, ... scan matrix (maybe partly filled)
                spcm, ... Spcm object
                stage, ... Stage object. to scan
                tPixel, ... double. Time per pixel (in sec)
                axisAPixelsPerLine, ... vector of double
                axisBLinesPerScan, ... vector of double
                axisADirectionIndex, ... integer. Index in {1, 2, 3} for "xyz"
                axisBDirectionIndex, ... integer. Index in {1, 2, 3} for "xyz"
                axisCPoint0, ... scalar. Zero point in the 3rd axis
                isFastScan, ... logical
                isFlipped, ... logical. Is the matrix flipped or not
                matrixIndexLineStart, ... integer. Index at which to start inserting lines to the matrix
                matrixIndexLineEnd, ... integer. Index at which to stop inserting lines to the matrix
                matrixIndexPixelInLineStart, ... integer. Index at which to start inserting pixels to the line in the matrix
                matrixIndexPixelInLineEnd, ... integer. Index in which to stop inserting pixels to the line in the matrix
                extraInfo, ... extra info matrix (maybe partly filled)
                lineIndexMap, ... (optional) matrix line of each scanned line, in scan order. Default: start:end
                flipPixels ... (optional) logical. axisAPixelsPerLine runs to -> from; store it reversed
                )
            % This function scans a chunk from the scan-matrix
            if ~exist('lineIndexMap', 'var') || isempty(lineIndexMap)
                lineIndexMap = matrixIndexLineStart : matrixIndexLineEnd;
            end
            if ~exist('flipPixels', 'var') || isempty(flipPixels)
                flipPixels = false;
            end
            if length(axisAPixelsPerLine) ~= matrixIndexPixelInLineEnd - matrixIndexPixelInLineStart + 1
                obj.sendError('Can''t scan - mismatch size of vector!')
            end
            if length(axisBLinesPerScan) ~= matrixIndexLineEnd - matrixIndexLineStart + 1
                obj.sendError('Can''t scan - mismatch size of vector!')
            end
            
            if ~obj.mCurrentlyScanning; return; end %Yoav: Will create error because no output. Nadav: or will it?
            if size(kcpsMatrix, 3) == 2 % MW Contrast
                pg = getObjByName(PulseGenerator.NAME);
            end
            
            % for each in {x, y, z}, it could be one of:
            % the axisA vector, the axisB vector, or the axisC point
            [x, y, z] = obj.getXYZfor2dScanChunk(axisAPixelsPerLine, axisBLinesPerScan, axisADirectionIndex, axisBDirectionIndex, axisCPoint0);
            nPixels = length(axisAPixelsPerLine);
            timeout = 2*nPixels*tPixel;
            
            % prepare scan
            nFlat = 0;      % A flat section at the start of ramp. parameter not needed genrally, a stage can overwrite if needed. BACKWARD_COPITABILITY
            nOverRun = 0;   % Let the waveform over run the start and end. Not needed genrally, a stage can overwrite if needed. BACKWARD_COPITABILITY
            axesLettersUpper = upper(ClassStage.SCAN_AXES([axisADirectionIndex, axisBDirectionIndex]));
            prapareScanFuncHndl = sprintf('PrepareScan%s', axesLettersUpper);
            feval(prapareScanFuncHndl, stage, x, y, z, nFlat, nOverRun, tPixel);
            spcm.prepareCountByStage(stage.name, nPixels, timeout, isFastScan);
            
            bHasLiftime = spcm.hasLifetime();
            isMWContrastScan = (size(kcpsMatrix, 3) == 2);
            
            % do the scan
            spcm.startScanCount();
            obj.mProgress.nLines = numel(lineIndexMap);
            lineCount = 0;      % lines of this chunk so far (the stage counts the same way)
            isZigzag = isprop(stage, 'scanBidirectional') && stage.scanBidirectional;
            for lineIndex = lineIndexMap
                if ~obj.mCurrentlyScanning; break; end
                obj.mProgress.line = obj.mProgress.line + 1;
                lineCount = lineCount + 1;
                success = false;
                for i = 1:size(kcpsMatrix, 3) % MW Contrast
                    if i==2
                        stage.PrepareRescanLine();
                        pg.on('MW');
                    end
                    for trial = 1 : StageScanner.TRIALS_AMOUNT_ON_ERROR    
                        try
                            if isa(spcm, 'PhotoDiodeDigitizerNiDaqControlled')
                                spcm.startScanCount();
                            end
                            obj.liveStart(kcpsMatrix, extraInfo, '2D', lineIndex, ...
                                matrixIndexPixelInLineStart, matrixIndexPixelInLineEnd, isFlipped, flipPixels, i, ...
                                ~isZigzag || mod(lineCount, 2) == 1);
                            wasScannedForward = stage.ScanNextLine();
                            obj.liveStop();
                            % forwards - if true, the line was scanned normally,
                            %            if false - should flip the results
                            [kcpsVector, errVector] = StageScanner.readScanLine(spcm);
                            kcpsVector = kcpsVector(:)';
                            errVector = errVector(:)';
                            if ~wasScannedForward
                                kcpsVector = fliplr(kcpsVector);
                                errVector = fliplr(errVector);
                            end
                            if flipPixels
                                kcpsVector = fliplr(kcpsVector);    % visited to -> from, stored from -> to
                                errVector = fliplr(errVector);
                            end
                            hasErr = isfield(extraInfo, 'sterr') && isequal(size(extraInfo.sterr), size(kcpsMatrix));
                            if isFlipped
                                kcpsMatrix(matrixIndexPixelInLineStart:matrixIndexPixelInLineEnd, lineIndex, i) = kcpsVector;
                                if hasErr; extraInfo.sterr(matrixIndexPixelInLineStart:matrixIndexPixelInLineEnd, lineIndex, i) = errVector; end
                            else
                                kcpsMatrix(lineIndex, matrixIndexPixelInLineStart:matrixIndexPixelInLineEnd, i) = kcpsVector;
                                if hasErr; extraInfo.sterr(lineIndex, matrixIndexPixelInLineStart:matrixIndexPixelInLineEnd, i) = errVector; end
                            end
                            
                            if bHasLiftime && ~isMWContrastScan
                                hist = spcm.lastScanHist;
                                hist = BooleanHelper.ifTrueElse(wasScannedForward, hist, fliplr(hist));
                                if flipPixels; hist = fliplr(hist); end
                                if isFlipped
                                    extraInfo.hist(matrixIndexPixelInLineStart:matrixIndexPixelInLineEnd, lineIndex, :) = hist;
                                else
                                    extraInfo.hist(lineIndex, matrixIndexPixelInLineStart:matrixIndexPixelInLineEnd, :) = hist;
                                end
                            end
                            
                            success = true;
                            break;
                        catch err
                            obj.liveStop();
%                             rethrow(err);  % Uncomment to debug
                            warning(err.message);
                            
                            if ~obj.mCurrentlyScanning
                                stage.AbortScan();
                                spcm.clearScanRead();
                                return;
                            end
                            
                            fprintf('Line %d failed at trial %d, attempting to rescan line.\n', lineIndex, trial);
                            
                            try
                                stage.PrepareRescanLine(); % Prepare to rescan the line
                            catch err2
                                stage.AbortScan();
                                spcm.clearScanRead();
                                rethrow(err2)
                            end
                        end % try catch
                    end % for trial = 1 : StageScanner.TRIALS_AMOUNT_ON_ERROR
                    if i==2
                        pg.off('MW');
                    end
                    obj.progressLineDone(nPixels);
                end % MW Contrast
                obj.publishScan(kcpsMatrix, extraInfo);
                
                if ~success
                    % We failed at reading the line, so there is probably
                    % no point in scanning next line
                    obj.mCurrentlyScanning = false;     % will abort the scan
                end
            end % lineIndex = 1 : length(axisBLinesPerScan)
            stage.AbortScan();
            spcm.clearScanRead();
        end
        
        function [x,y,z] = getXYZfor2dScanChunk(obj, axisAPointsPerLine, axisBLinesPerScan, axisADirectionIndex, axisBDirectionIndex, axisCPoint0) %#ok<INUSD,STOUT,INUSL>
            % Converts from "axisA", "axisB", to xyz: calculates the vectors for x,y,z to be used in scan2dChunk()
            for letter = ClassStage.SCAN_AXES
                % Iterate over the string "xyz" letter by letter
                ll = lower(letter);
                switch ClassStage.getAxis(letter)
                    case axisADirectionIndex
                        eval(sprintf('%s = axisAPointsPerLine;', ll));
                    case axisBDirectionIndex
                        eval(sprintf('%s = axisBLinesPerScan;', ll));
                    otherwise
                        eval(sprintf('%s = axisCPoint0;', ll));
                end
            end
            % todo: this needs to be redone. No need to use eval here.
        end
        
        
        function [volume, extraInfo] = scan3D(obj, stage, spcm, scanParams, optionalVolume)
            % A 3D scan: one 2D scan per slice. scanParams.scanOrder decides the axes: the
            % first runs along lines, the second steps between lines, the third steps
            % between slices (default: xy slices, one per z).
            %
            % volume: [rows x cols x slices (x 2 with MW contrast)]. Rows and columns are the
            % two image axes in xyz order, exactly as a 2D scan of those axes stores them, so
            % every slice is a normal 2D image. Slice k is at the k-th point of the slice
            % axis (from -> to, whatever order they were visited in).
            extraInfo = struct();
            [imgAxes, sliceAxis] = StageScanner.volumeAxes(scanParams);
            order = scanParams.orderedScanAxes;
            sliceVector = scanParams.getScanAxisVector(sliceAxis);
            nSlices = numel(sliceVector);
            nRows = scanParams.numPoints(imgAxes(2));
            nCols = scanParams.numPoints(imgAxes(1));
            nLayers = 1 + logical(scanParams.isMWContrastScan);
            if exist('optionalVolume', 'var') && isequal(size(optionalVolume, 1:4), [nRows nCols nSlices nLayers])
                volume = optionalVolume;
            else
                volume = nan(nRows, nCols, nSlices, nLayers);
            end
            sliceOrder = 1:nSlices;
            if numel(scanParams.reverse) >= sliceAxis && scanParams.reverse(sliceAxis)
                sliceOrder = fliplr(sliceOrder);
            end

            obj.mProgress.nSlices = nSlices;
            for n = 1:nSlices
                if ~obj.mCurrentlyScanning; break; end
                k = sliceOrder(n);
                obj.mProgress.slice = n;
                obj.mProgress.line = 0;

                % This slice is a 2D scan with the slice axis held at its k-th point
                sliceParams = scanParams.copy;
                sliceParams.scanBlocks = {};
                sliceParams.isFixed(sliceAxis) = true;
                sliceParams.fixedPos(sliceAxis) = sliceVector(k);
                sliceParams.scanOrder = ClassStage.SCAN_AXES(order(1:2));

                obj.mVolume = struct('vol', volume, 'slice', k);
                obj.publishScan(volume(:, :, k, :), struct());    % show the slice being scanned
                slice = obj.scan2D(stage, spcm, sliceParams);
                if isempty(obj.mVolume)     % (cleared if the scan was stopped from outside)
                    break
                end
                volume = obj.mVolume.vol;
                if ~isempty(slice)
                    volume(:, :, k, :) = reshape(slice, nRows, nCols, 1, []);
                end
            end
            obj.mVolume = [];
            obj.mLive = [];
            obj.mScan = volume;
            obj.mScanExtraInfo = extraInfo;
        end

        function publishScan(obj, matrix, extraInfo)
            % Tell everybody about a partly scanned image. During a 3D scan, the matrix is
            % the slice being scanned: it goes into the volume, and the volume is sent.
            if isempty(obj.mVolume)
                obj.sendEventScanUpdated(matrix, extraInfo);
                return
            end
            k = obj.mVolume.slice;
            vol = obj.mVolume.vol;
            if ~isempty(matrix)
                vol(:, :, k, :) = reshape(matrix, size(vol, 1), size(vol, 2), 1, []);
            end
            obj.mVolume.vol = vol;
            extra3D = struct('sliceIndex', k);
            if isstruct(extraInfo) && isfield(extraInfo, 'livePixel')
                extra3D.livePixel = true;
            end
            obj.sendEventScanUpdated(vol, extra3D);
        end

        %% Progress
        function progressStart(obj, scanParams, pass)
            n = scanParams.numPoints(~scanParams.isFixed);
            obj.mProgress = struct('done', 0, 'inLine', 0, ...
                'total', prod(n) * (1 + logical(scanParams.isMWContrastScan)), ...
                'tic', tic, 'lastSent', [], 'line', 0, 'nLines', 0, 'slice', 0, 'nSlices', 0, ...
                'pass', pass, 'dims', numel(n));
            obj.progressSend(true);
        end

        function progressLineDone(obj, nPixels)
            % A line (or a 1D chunk) is complete
            obj.mProgress.done = min(obj.mProgress.done + nPixels, obj.mProgress.total);
            obj.mProgress.inLine = 0;
            obj.progressSend(false);
        end

        function progressPixelDone(obj)
            % One pixel of the current line is measured (the stage calls notePixel)
            obj.mProgress.inLine = obj.mProgress.inLine + 1;
            obj.progressSend(false);
            obj.livePreview();
        end

        %% Drawing a line while it is scanned
        % The image used to change once per line (a 1D scan: once at the end). Now, while a
        % line runs, the pixels measured so far are drawn too, at most every
        % LIVE_MIN_INTERVAL seconds -- for detectors that can say what they have measured so
        % far (partialScanLine: the photodiode in its default 'software' mode). It costs
        % little: the image only swaps its pixel data (ImageScanResult.updateImageInPlace),
        % and the redraw does not run button callbacks (they wait for the end of the line).
        function liveStart(obj, matrix, extraInfo, kind, row, pixelStart, pixelEnd, isFlipped, flipPixels, layer, forwards)
            % matrix: the scan so far. kind: '1D' or '2D'. row: the matrix line (2D).
            % pixelStart..pixelEnd: where this line goes. isFlipped: the matrix is stored
            % [pixels x lines]. flipPixels: the line was given to the stage to -> from.
            % layer: 1, or 2 for the MW-on pass. forwards: the stage runs this line in the
            % order it was given (false: a zigzag line running back).
            obj.mLive = struct('matrix', matrix, 'extra', extraInfo, 'kind', kind, 'row', row, ...
                'p0', pixelStart, 'p1', pixelEnd, 'isFlipped', isFlipped, 'flip', logical(flipPixels), ...
                'layer', layer, 'forwards', logical(forwards), 'lastDraw', tic, 'shown', 0);
        end

        function liveStop(obj)
            obj.mLive = [];
        end

        function livePreview(obj)
            L = obj.mLive;
            if isempty(L) || toc(L.lastDraw) < obj.LIVE_MIN_INTERVAL; return; end
            spcm = getObjByName(Spcm.NAME);
            if isempty(spcm) || ~ismethod(spcm, 'partialScanLine'); return; end
            v = spcm.partialScanLine();
            n = L.p1 - L.p0 + 1;
            k = numel(v);
            if k == 0 || k > n || k == L.shown; return; end
            % pixel p of the line, as measured -> its place in the matrix
            idx = 1:k;
            if ~L.forwards; idx = n - idx + 1; end
            if L.flip; idx = n - idx + 1; end
            pos = L.p0 - 1 + idx;
            M = L.matrix;
            if strcmp(L.kind, '1D')
                M(pos, L.layer) = v(:);
            elseif L.isFlipped
                M(pos, L.row, L.layer) = v(:);
            else
                M(L.row, pos, L.layer) = v(:)';
            end
            obj.mLive.lastDraw = tic;
            obj.mLive.shown = k;
            extra = L.extra;
            extra.livePixel = true;
            obj.publishScan(M, extra);
        end

        function info = progressInfo(obj)
            % fraction, done, total, elapsed (s), remaining (s), line, nLines, slice,
            % nSlices, pass, dims, running
            p = obj.mProgress;
            done = min(p.done + p.inLine, p.total);
            info = struct('fraction', 0, 'done', done, 'total', p.total, 'elapsed', 0, ...
                'remaining', NaN, 'line', p.line, 'nLines', p.nLines, 'slice', p.slice, ...
                'nSlices', p.nSlices, 'pass', p.pass, 'dims', p.dims, ...
                'running', obj.mCurrentlyScanning);
            if p.total > 0
                info.fraction = done / p.total;
            end
            if ~isempty(p.tic)
                info.elapsed = toc(p.tic);
                if done > 0
                    info.remaining = info.elapsed / done * (p.total - done);
                end
            end
        end

        function progressSend(obj, force)
            % Progress event, at most every PROGRESS_MIN_INTERVAL seconds unless forced
            p = obj.mProgress;
            if ~force && ~isempty(p.lastSent) && toc(p.lastSent) < obj.PROGRESS_MIN_INTERVAL
                return
            end
            obj.mProgress.lastSent = tic;
            obj.sendEvent(struct(obj.EVENT_SCAN_PROGRESS, obj.progressInfo));
        end

        function returnToStart(obj, stage, startPos)
            % After a scan (finished, stopped, or failed): move the stage back to where it
            % was when the scan started -- unless scanParams.returnToStart is off, or the
            % stages were halted (Halt means "stop moving").
            p = obj.mStageScanParams;
            if obj.mHalted || isempty(startPos) || ~isa(p, 'StageScanParams') || ~p.returnToStart
                return
            end
            if all(abs(stage.Pos(stage.availableAxes) - startPos) < 1e-6); return; end
            stage.move(stage.availableAxes, startPos);  % (sends an error event itself if it fails)
        end

        function stopScan(obj)
            obj.mCurrentlyScanning = false;
            % todo: check maybe we need also
            %   stage = getObjByName(obj.mStageName);
            %   stage.scanRunning = false;
        end
        
        function switchTo(obj, newStageName)
            if obj.mCurrentlyScanning && ~strcmp(obj.mStageName, newStageName)
                obj.sendError('Can''t switch stage when scan running!');
            end
            obj.mStageName = newStageName;
        end
        
        function autosaveAfterScan(obj) %#ok<MANU>
            saveLoad = SaveLoad.getInstance(Savable.CATEGORY_IMAGE);
            saveLoad.autoSave();
        end
        
        function number = getScanDimensions(obj)
            number = sum(~obj.mStageScanParams.isFixed);
        end
        
        function [globalFirstAxisPos, globalSecondAxisPos] = getGlobalScanPosition(obj, stagesCell)
            % Returns the location of the scan-vector(s), relative to a
            % global zero
            pos = [0, 0, 0];
            % First get data from scanning stage
            
            params = obj.mStageScanParams;
            fixedInd = find(params.isFixed);
            pos(fixedInd) = pos(fixedInd) + params.fixedPos(fixedInd);
            for i = 1:length(stagesCell)
                stage = stagesCell{i};
                if ~strcmp(stage.name, obj.mStageName)
                    phAxes = stage.availableAxes;
                    ind = stage.getAxis(phAxes);
                    pos(ind) = pos(ind) + stage.Pos(ind);
                end
            end
            
            firstAxisInd = params.getFirstScanAxisIndex;
            globalFirstAxisPos = pos(firstAxisInd);
            
            secondAxisInd = params.getSecondScanAxisIndex;
            if secondAxisInd == -1
                globalSecondAxisPos = 0;
            else
                globalSecondAxisPos = pos(secondAxisInd);
            end
        end
    end
    
    methods
        function clear(obj)
            if (obj.mCurrentlyScanning)
                obj.stopScan();
            end
            
            obj.mScan = nan;
            obj.mScanExtraInfo = nan;
            obj.mStageName = nan;
            obj.mStageScanParams = nan;
            obj.mCurrentlyScanning = false;     % probably redundant; appears in obj.stopScan
        end
        
        function dummyScan(obj)
            obj.sendEventScanStarting();
            obj.mStageName = ClassStage.getStages{1}.name;
            data = 5 + peaks;
            obj.mScan = data;
            obj.mScanExtraInfo = struct();
            obj.mStageScanParams = StageScanParams([0,0,0], [48,48,48], [49,49,49], [0,0,0], [false false true], 1, 0, 0, 0);
            obj.sendEventScanUpdated(data, struct());
            obj.sendEventScanFinished();
            obj.autosaveAfterScan();
        end
        
        function value = dummyScanGaussian(obj,scanParams)
            pos = scanParams.fixedPos;
            stage = getObjByName(obj.mStageName);
            if isempty(stage); throwBaseObjException(obj.mStageName); end
            
            phAxes = stage.availableAxes;
            stage.move(phAxes, pos);
            
            X = scanParams.getScanAxisVector(1);
            Y = scanParams.getScanAxisVector(2);
            Z = scanParams.getScanAxisVector(3);
            f = @(x,y,z) 100*exp(-((z+5).^2+x.^2+(y-7).^2)/100);   % some test function
            value = f(X,Y,Z);
        end
        
        function boolean = isScanReady(obj)
            boolean = ~isnan(obj.mScan);
        end
        
        function string = getBottomScanLabel(obj)
            axisLetters = obj.mStageScanParams.getScanAxes;
            switch obj.getScanDimensions
                case 1
                    string = sprintf('%s (%s)', axisLetters, StringHelper.MICRON); % for example, 'x (?m)'
                case 2
                    string = sprintf('%s (%s)', axisLetters(1), StringHelper.MICRON); % (as above)
                case 3
                    imgAxes = StageScanner.volumeAxes(obj.mStageScanParams);
                    string = sprintf('%s (%s)', ClassStage.SCAN_AXES(imgAxes(1)), StringHelper.MICRON);
                otherwise
                    string = '';
            end
        end
        
        function string = getLeftScanLabel(obj)
            switch obj.getScanDimensions
                case 1
                    string = Spcm.signalUnits();    % 'V' with the photodiode, 'kcps' with a counter
                case 2
                    axisLetters = obj.mStageScanParams.getScanAxes;
                    string = sprintf('%s (%s)', axisLetters(2), StringHelper.MICRON); % for example, 'y (?m)'
                case 3
                    imgAxes = StageScanner.volumeAxes(obj.mStageScanParams);
                    string = sprintf('%s (%s)', ClassStage.SCAN_AXES(imgAxes(2)), StringHelper.MICRON);
                otherwise
                    string = '';
            end
        end
    end
    
    methods (Static)
        function obj = init
            obj = getObjByName(StageScanner.NAME);
            if isempty(obj)
                obj = StageScanner;
                addBaseObject(obj);
            end
        end
        
        function [axisAIndex, axisBIndex] = lineAndStepAxes(stageMaxScanSizeInt, scanParams)
            % The axis each line runs along (A) and the axis stepped between lines (B), for
            % a 2D scan: as scanParams.scanOrder says, or chosen by
            % optimize2dScanDirections when it is automatic.
            if scanParams.isOrderAutomatic
                [axisAIndex, axisBIndex] = StageScanner.optimize2dScanDirections(stageMaxScanSizeInt, scanParams.copy);
            else
                order = scanParams.orderedScanAxes;
                axisAIndex = order(1);
                axisBIndex = order(2);
            end
        end

        function [imgAxes, sliceAxis] = volumeAxes(scanParams)
            % For a 3D scan: the two axes of each slice image (in xyz order) and the axis
            % stepped between slices (the last in the scan order).
            order = scanParams.orderedScanAxes;
            sliceAxis = order(end);
            imgAxes = sort(order(1:end-1));
        end

        function [values, sterr] = readScanLine(spcm)
            % The scanned line, and each pixel's standard error when the detector gives one
            % (the photodiodes: the spread of the digitizer samples in the pixel, / sqrt(N)).
            % Photon counters give values only: sterr is then NaN. Which kind a detector is,
            % is read from readFromScan's signature -- calling it twice would read twice.
            persistent cache
            if isempty(cache); cache = containers.Map; end
            cls = class(spcm);
            if ~isKey(cache, cls)
                m = findobj(metaclass(spcm).MethodList, 'Name', 'readFromScan');
                cache(cls) = ~isempty(m) && numel(m(1).OutputNames) >= 2;
            end
            if cache(cls)
                [values, sterr] = spcm.readFromScan();
                if numel(sterr) ~= numel(values); sterr = nan(size(values)); end
            else
                values = spcm.readFromScan();
                sterr = nan(size(values));
            end
        end

        function notePixel()
            % Called by a stage after each pixel it measures during a scan (Standa, dummy),
            % so the progress bar moves within a line. Cheap; does nothing outside a scan.
            scanner = getObjByName(StageScanner.NAME);
            if ~isempty(scanner) && scanner.mCurrentlyScanning
                scanner.progressPixelDone();
            end
        end

        function [axisAIndex, axisBIndex] = optimize2dScanDirections(stageMaxScanSizeInt, scanParams)
            % calculates the best way to scan a 2d scan (given by the scanParams)
            % the movement between pixels in each line is called "axis a",
            % while movement between lines is called "axis b"
            %
            % returns:
            % axisAIndex - for each line in the 2d-scan, scan the pixels in the line in this direction
            % axisBIndex - move between lines in this direction
            %
            % for example, in a 2D scan (XY) where x is [0...30] and y is
            % [0...70] and maxScanSize is 100, the returned values will be:
            % axisAIndex  = 2 (Y)
            % axisBIndex  = 1 (X)
            % as the optimized results will be 30 scans
            
%             %%% Ty debug setup 1 10/04/21
%             firstAxisIndex = scanParams.getFirstScanAxisIndex();
%             secondAxisIndex = scanParams.getSecondScanAxisIndex();
%             axisAIndex = secondAxisIndex;
%             axisBIndex = firstAxisIndex;

            firstAxisIndex = scanParams.getFirstScanAxisIndex();
            secondAxisIndex = scanParams.getSecondScanAxisIndex();
            
            firstAxisNumPoints = scanParams.numPoints(firstAxisIndex);
            secondAxisNumPoints = scanParams.numPoints(secondAxisIndex);
            
            % Option ONE: scan each line by axis 1, move between lines in axis 2
            % Calculate how many steps would be needed then:
            timesMaxInFirst = ceil(firstAxisNumPoints / stageMaxScanSizeInt);
            totalTimesOptionOne = timesMaxInFirst * secondAxisNumPoints;
            
            % Option TWO: scan each line by axis 2, move between lines in axis 1
            % Calculate how many steps would be needed then:
            timesMaxInSecond = ceil(secondAxisNumPoints / stageMaxScanSizeInt);
            totalTimesOptionTwo = timesMaxInSecond * firstAxisNumPoints;
            
            if totalTimesOptionOne <= totalTimesOptionTwo
                axisAIndex = firstAxisIndex;
                axisBIndex = secondAxisIndex;
            else
                axisAIndex = secondAxisIndex;
                axisBIndex = firstAxisIndex;
            end
        end
    end
    
    %% overridden from EventListener
    methods
        % When events happen, this function jumps.
        % event is the event sent from the EventSender
        function onEvent(obj, event)
            % Halt: the scan is not followed by a return to the start position
            if isfield(event.extraInfo, StageControlEvents.HALT)
                obj.mHalted = true;
                return
            end
            
            % Check if event is "loaded file to SaveLoad" and need to show the image
            if strcmp(event.creator.name, SaveLoadCatImage.NAME) ...
                    && isfield(event.extraInfo, SaveLoad.EVENT_LOAD_SUCCESS_FILE_TO_LOCAL)
                % Need to load the image!
                category = Savable.CATEGORY_IMAGE;
                subcat = Savable.SUB_CATEGORY_DEFAULT;
                saveLoad = event.creator;
                struct = saveLoad.getStructToSavable(obj);
                if ~isempty(struct)
                    obj.loadStateFromStruct(struct, category, subcat);
                end
            end
        end
    end
    
    %% overriding from Savable
    methods (Access = protected)
        function outStruct = saveStateAsStruct(obj, category, type)
            % Saves the state as struct. if you want to save stuff, make
            % (outStruct = struct;) and put stuff inside. If you dont
            % want to save, make (outStruct = NaN;)
            %
            % category - string. Some objects saves themself only with
            %                    specific category (image/experiments/etc.)
            % type - string.     Whether the objects saves at the beginning
            %                    of the run (parameter) or at its end (result)
            if ~strcmp(category, Savable.CATEGORY_IMAGE)
                outStruct = nan;
                return
            end
            
            outStruct = struct;
            switch type
                case Savable.TYPE_PARAMS
                    outStruct.scanParams = obj.mStageScanParams.asStruct();
                    outStruct.stageName = obj.mStageName;
                case Savable.TYPE_RESULTS
                    if isnan(obj.mScan)
                        outStruct = NaN;
                    else
                        outStruct.scan = obj.mScan;
                        outStruct.scanExtraInfo = obj.mScanExtraInfo;
                    end
            end
        end
        
        function loadStateFromStruct(obj, savedStruct, category, subCategory)
            % Loads the state from a struct.
            % to support older versions, always check for a value in the
            % struct before using it. View example in the first line.
            if ~strcmp(category, Savable.CATEGORY_IMAGE); return; end
            if ~any(strcmp(subCategory, {Savable.CATEGORY_IMAGE_SUBCAT_STAGE})); return; end
            
            
            % for field = {'scan', 'scanParams', 'stageName'} - removed 'scan', for when scan has not yet been saved.
            for field = {'scanParams', 'stageName'}
                if ~isfield(savedStruct, field{:})
                    return
                end
            end
            if ~isfield(savedStruct, 'scan')
                savedStruct.scan = [];
                obj.sendWarning('No scan results found. Loading only scan parameters');
            end
            stage = getObjByName(savedStruct.stageName);
            if isempty(stage)
                obj.sendError(sprintf('Can''t load stage! No stage with name "%s"', savedStruct.stageName));
            end
            obj.mStageScanParams = StageScanParams.fromStruct(savedStruct.scanParams);
            stage.scanParams = obj.mStageScanParams; % This will send an event that the scan-parameters have changed
            
            obj.mScan = savedStruct.scan;
            obj.mStageName = savedStruct.stageName;
            if isfield(savedStruct, 'extraInfo')
                obj.mScanExtraInfo = savedStruct.extraInfo;
            else
                obj.mScanExtraInfo = struct();
            end
            obj.sendEventScanUpdated(savedStruct.scan, obj.mScanExtraInfo);
        end
        
        function string = returnReadableString(~, savedStruct)
            % Return a readable string to be shown. If this object
            % doesn't need a readable string, make (string = NaN;) or
            % (string = '');
            
            string = NaN;
            % for field = {'scan', 'scanParams', 'stageName'} - removed 'scan', for when scan has not yet been saved.
            for field = {'scanParams', 'stageName'}
                if ~isfield(savedStruct, field{:})
                    return
                end
            end
            
            % Getting initial information from stage ...
            stageName = savedStruct.stageName;
            stage = getObjByName(stageName);
            if isempty(stage); throwBaseObjException(obj.mStageName); end
            availAxes = stage.getAxis(stage.availableAxes);
            
            % and from scan parameters, ...
            params = savedStruct.scanParams;
            scanAxes = find(~params.isFixed);
            axesLetters = upper(ClassStage.GetLetterFromAxis(scanAxes));
            fixAxes = find(params.isFixed);
            fixAxes = intersect(fixAxes, availAxes);    % Some of the "fixed axes" might not actually be real
            
            % so we can create the output string
            string = sprintf('Scanning %s %s:', stageName, axesLetters);
            
            indentation = 5;
            for i = 1:length(fixAxes)
                index = fixAxes(i);
                ax = ClassStage.GetLetterFromAxis(index);
                position = params.fixedPos(index);
                axisString = sprintf('%s position: %.3f', upper(ax), position);
                string = sprintf('%s\n%s', string, ...
                    StringHelper.indent(axisString, indentation));
            end
            for i = 1:length(scanAxes)
                index = scanAxes(i);
                ax = ClassStage.GetLetterFromAxis(index);
                from = params.from(index);
                to = params.to(index);
                numPoints = params.numPoints(index);
                axisString = sprintf('%s: %d points from %.3f to %.3f', upper(ax), numPoints, from, to);
                string = sprintf('%s\n%s', string, ...
                    StringHelper.indent(axisString, indentation));
            end
            
        end
    end
end