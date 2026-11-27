classdef StageScanParams < handle
    %STAGESCANPARAMS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        DEFAULT_POINT_NUM = 51;
        DEFAULT_PIXEL_TIME = 0.015;
    end
    
    properties
        % all below are number-arrays with size 1x3
        from
        to
        numPoints
        fixedPos
        
        % boolean(logicals \ 1\0)-with size 1x3
        isFixed
        
        % double
        pixelTime
        MWAmplitude
        MWFrequency
        
        % all below are logicals (boolean)
        continuous
        fastScan
        autoSave
        isMWContrastScan
        
        % Advanced map indicating where to scan:
        % Should be a cell array of vectors of length 5 [x1, x2, y1, y2, z].
        % Each cell is a block from (x1, y1) to (x2, y2), scanned at "z".
        scanBlocks = {}; % empty means that it is a normal scan.

        % ---- how the scan runs (main window, "Scan" section) ----
        % scanOrder: the scanned axes in the order they move. The first runs along each
        % line, the second steps between lines, the third (3D) steps between slices.
        % E.g. 'yxz' = lines along y, stepping x, one xy slice per z. Letters of axes that
        % are not scanned are ignored. '' = automatic (lines along the first scanned axis,
        % unless the stage's maximum line length makes the other way shorter).
        scanOrder = ''
        % zigzag: true = every other line runs backwards (serpentine, no return trip).
        % false = every line runs in the same direction (immune to backlash). Only for
        % stages that can do both (they have a scanBidirectional property).
        zigzag = false
        % reverse(i): true = axis i is scanned from "to" towards "from". The image is
        % always stored and shown with the axis in its from -> to order.
        reverse = [false false false]
        % returnToStart: when the scan ends (or is stopped), move the stage back to where it
        % was when the scan started. Not after Halt.
        returnToStart = true
    end

    properties (Constant, Hidden)
        % Added after files were already saved without them: load silently with defaults.
        LATER_FIELDS = {'scanOrder', 'zigzag', 'reverse', 'returnToStart', 'scanBlocks'};
    end
    
    methods
        function obj = StageScanParams(from, to, numPoints, fixedPos, isFixed, pixelTime, continuous, fastScan, autoSave, isMWContrastScan, MWAmplitude, MWFrequency)
            
            %%%% default no-args constructor if needed %%%%
            if nargin == 0
                obj.from = zeros(1,length(ClassStage.SCAN_AXES));
                obj.to = zeros(1,length(ClassStage.SCAN_AXES));
                obj.fixedPos = zeros(1,length(ClassStage.SCAN_AXES));
                obj.numPoints = obj.DEFAULT_POINT_NUM * ones(1,length(ClassStage.SCAN_AXES));
                obj.isFixed = false(1,length(ClassStage.SCAN_AXES));
                obj.pixelTime = obj.DEFAULT_PIXEL_TIME;
                obj.continuous = false;
                obj.fastScan = false;
                obj.autoSave = true;
                obj.isMWContrastScan = false;
                obj.MWAmplitude = -10;
                obj.MWFrequency = 2870;
                return
            end
            
            %%%% input params checks %%%%
            objectsToCheckNumeric1x3Vector = {'from', 'to', 'numPoints', 'fixedPos'};
            
            for i = 1 : length(objectsToCheckNumeric1x3Vector)
                funcParamName = objectsToCheckNumeric1x3Vector{i};
                funcParamValue = eval(funcParamName);
                if ~isnumeric(funcParamValue) || length(funcParamValue) ~= 3
                    EventStation.anonymousError('"%s" parameter should be a 1x3 number array!', funcParamName)
                end
            end
            
            % check if "pixelTime" is a positive double
            if ~ValidationHelper.isValuePositive(pixelTime)
                EventStation.anonymousError('pixelTime parameter should be a positive double!')
            end
            
            % check if "isFixed" is logical array of size 3
            isFixedIsBooleanArray = all(or(isFixed==0, isFixed==1));
            if ~isFixedIsBooleanArray || length(isFixed) ~= 3
                EventStation.anonymousError('"isFixed" parameter should be a 1x3 boolean array!')
            end
            
            % check if the booleans are actually booleans
            objects_to_check_boolean = {'continuous', 'fastScan', 'autoSave', 'isMWContrastScan'};
            for i = 1 : length(objects_to_check_boolean)
                funcParamName = objects_to_check_boolean{i};
                funcParamValue = eval(funcParamName);
                if ~ValidationHelper.isTrueOrFalse(funcParamValue)
                    EventStation.anonymousError('"%s" should be logical!', funcParamName);
                end
            end
            
            
            
            %%%% create object %%%%
            obj.from = from;
            obj.to = to;
            obj.numPoints = numPoints;
            obj.fixedPos = fixedPos;
            obj.isFixed = isFixed;
            obj.pixelTime = pixelTime;
            obj.continuous = continuous;
            obj.fastScan = fastScan;
            obj.autoSave = autoSave;
            obj.isMWContrastScan = isMWContrastScan;
            obj.MWAmplitude = MWAmplitude;
            obj.MWFrequency = MWFrequency;
        end
        
        function outStruct = asStruct(obj)
            outStruct = jsondecode(jsonencode(obj));
        end
        
        function wasChanged = updateByLimit(obj, phAxis, newLimNeg, newLimPos)
            % phAxis - the axis to change
            % newLimNeg, newLimPos - the new limits to consider
            % supports vectorial axis!
            % NO INPUT CHECKS!
            phAxis = ClassStage.getAxis(phAxis);
            oldFrom = obj.from;
            oldTo = obj.to;
            oldFixedPos = obj.fixedPos;
            
            from_ = obj.from(phAxis);
            from_(from_ > newLimPos) = newLimPos(from_ > newLimPos);
            from_(from_ < newLimNeg) = newLimNeg(from_ < newLimNeg);
            obj.from(phAxis) = from_;
            
            to_ = obj.to(phAxis);
            to_(to_ > newLimPos) = newLimPos(to_ > newLimPos);
            to_(to_ < newLimNeg) = newLimNeg(to_ < newLimNeg);
            obj.to(phAxis) = to_;
            
            fixedPos_ = obj.fixedPos(phAxis);
            fixedPos_(fixedPos_ > newLimPos) = newLimPos(fixedPos_ > newLimPos);
            fixedPos_(fixedPos_ < newLimNeg) = newLimNeg(fixedPos_ < newLimNeg);
            obj.fixedPos(phAxis) = fixedPos_;
            
            sameAsBefore = all(oldFrom == obj.from) ...
                && all(oldTo == obj.to) ...
                && all(oldFixedPos == obj.fixedPos);
            wasChanged = ~sameAsBefore;
        end
        
        function newScanParams = copy(obj)
            newScanParams = StageScanParams.fromStruct(obj.asStruct());
        end
        
        function plotScanMap(obj)
            if isempty(obj.scanBlocks)
                return
            end
            firstScanVector = obj.getFirstScanAxisVector;
            secondScanVector = obj.getSecondScanAxisVector;
            scanMap = zeros(length(firstScanVector), length(secondScanVector));
            for i = 1:length(obj.scanBlocks)
                fromA = obj.scanBlocks{i}(1);
                toA = obj.scanBlocks{i}(2);
                fromB = obj.scanBlocks{i}(3);
                toB = obj.scanBlocks{i}(4);
                scanMap(firstScanVector >= fromA & firstScanVector <= toA, ...
                        secondScanVector >= fromB & secondScanVector <= toB) = 1;
            end
            figure;
            imagesc(firstScanVector, secondScanVector, scanMap');
            xlabel(obj.getFirstScanAxisLetter)
            ylabel(obj.getSecondScanAxisLetter)
            axis xy
        end
    end
    
    % scan order, direction
    methods
        function n = nDims(obj)
            % Number of scanned axes (0-3)
            n = sum(~obj.isFixed);
        end

        function order = orderedScanAxes(obj)
            % Indices of the scanned axes, in scanOrder (line axis first). Axes scanOrder
            % does not name follow in xyz order. Empty scanOrder = xyz order.
            scanned = find(~obj.isFixed);
            order = [];
            if ischar(obj.scanOrder)
                for letter = lower(obj.scanOrder)
                    k = strfind(ClassStage.SCAN_AXES, letter);
                    if ~isempty(k) && any(scanned == k) && ~any(order == k)
                        order(end + 1) = k; %#ok<AGROW>
                    end
                end
            end
            order = [order, setdiff(scanned, order, 'stable')];
        end

        function tf = isOrderAutomatic(obj)
            tf = isempty(obj.scanOrder);
        end

        function v = getScanAxisVectorInScanOrder(obj, axisIndex)
            % The points of an axis in the order they are visited (to -> from if reversed).
            v = obj.getScanAxisVector(axisIndex);
            if ~obj.isFixed(axisIndex) && numel(obj.reverse) >= axisIndex && obj.reverse(axisIndex)
                v = fliplr(v);
            end
        end

        function s = describe(obj, maxLineLength)
            % One line saying what the scan will do, for the GUI:
            % '2D · lines along X, stepping Y · 51 × 41 = 2,091 pixels'
            if nargin < 2; maxLineLength = inf; end
            order = obj.orderedScanAxes;
            L = upper(ClassStage.SCAN_AXES);
            n = obj.numPoints(order);
            switch numel(order)
                case 0
                    s = 'Nothing to scan: tick Scan for at least one axis.';
                    return
                case 1
                    s = sprintf('1D along %s', L(order));
                case 2
                    if obj.isOrderAutomatic && maxLineLength < max(n)
                        s = sprintf('2D %s%s (line axis: automatic)', L(sort(order)));
                    else
                        s = sprintf('2D · lines along %s, stepping %s', L(order(1)), L(order(2)));
                    end
                case 3
                    s = sprintf('3D · lines along %s, stepping %s, slices along %s', L(order(1)), L(order(2)), L(order(3)));
            end
            total = prod(n) * (1 + logical(obj.isMWContrastScan));
            s = sprintf('%s · %s = %s pixels', s, strjoin(arrayfun(@num2str, n, 'UniformOutput', false), ' × '), ...
                StageScanParams.thousands(total));
        end
    end

    methods (Static)
        function s = thousands(n)
            s = sprintf('%d', round(n));
            s = regexprep(s, '(\d)(?=(\d{3})+$)', '$1,');
        end
    end

    % getters
    methods
        function string = getScanAxes(obj)
            string = '';
            for i = 1 : ClassStage.SCAN_AXES_SIZE
                if ~obj.isFixed(i)
                    phAxis = ClassStage.SCAN_AXES(i);
                    string = [string phAxis]; %#ok<AGROW>
                end
            end
        end
        
        function string = getFixedAxes(obj)
            obj.sendWarningIfNotScannable();
            string = '';
            for i = 1 : ClassStage.SCAN_AXES_SIZE
                if obj.isFixed(i)
                    phAxis = ClassStage.SCAN_AXES(i);
                    string = [string phAxis]; %#ok<AGROW>
                end
            end
        end
        
        function index = getFirstScanAxisIndex(obj)
            obj.sendWarningIfNotScannable();
            for i = 1 : ClassStage.SCAN_AXES_SIZE
                if ~obj.isFixed(i)
                    index = i;
                    return
                end
            end
            index = -1;
        end
        
        function index = getSecondScanAxisIndex(obj)
            isScannable = obj.sendWarningIfNotScannable;
            if ~isScannable
                index = -1;
                return
            end
            
            for i = obj.getFirstScanAxisIndex + 1 : ClassStage.SCAN_AXES_SIZE
                if ~obj.isFixed(i)
                    index = i;
                    return
                end
            end
            index = -1;  % if not found
        end

        function phAxis = getFirstScanAxisVector(obj)
            obj.sendWarningIfNotScannable();
            phAxis = obj.getScanAxisVector(obj.getFirstScanAxisIndex);
        end
        
        function phAxis = getSecondScanAxisVector(obj)
            phAxis = obj.getScanAxisVector(obj.getSecondScanAxisIndex);
        end
        
        function axisVector = getScanAxisVector(obj, axisIndex)
            % gets a vector to scan, by the axis-index
            % or, if it's a fixed position, gets the position
            if axisIndex < 0 || axisIndex > length(obj.from)
                axisVector = [];
                return
            end
            
            if obj.isFixed(axisIndex)
                axisVector = obj.fixedPos(axisIndex);
                return
            end
            
            first = obj.from(axisIndex);
            last = obj.to(axisIndex);
            jumps = obj.numPoints(axisIndex) -1;
            dist = (last - first) / jumps;
            axisVector = first : dist : last;
        end
        
        function string1Letter = getFirstScanAxisLetter(obj)
            string1Letter = ClassStage.SCAN_AXES(obj.getFirstScanAxisIndex);
        end
        function string1Letter = getSecondScanAxisLetter(obj)
            % returns a 1char string, or empty string if this scan
            % parameters object doesn't request a 2 dim scan
            index = obj.getSecondScanAxisIndex;
            if index > ClassStage.SCAN_AXES_SIZE || index <= 0
                string1Letter = '';
                return
            end
            string1Letter = ClassStage.SCAN_AXES(obj.getSecondScanAxisIndex);
        end
    end
    
    % Helper methods
    methods (Access = protected)
        function isOk = sendWarningIfNotScannable(obj)
            if all(obj.isFixed)
                EventStation.anonymousWarning('Attention - StageScanParams object is not scan-friendly! All the axes are fixed!');
                isOk = false;
            else
                isOk = true;
            end
        end
        
    end
    
    methods (Static)
        function propNames = varProps
            % This should probably have been implemebted as a method of
            % some superclass, if anybody knows how to
            mc = metaclass(StageScanParams);
            allProps = mc.PropertyList;
            mProps = allProps(not([allProps.Constant]));
            propNames = {mProps.Name};
        end
        
        function obj = fromStruct(inputStruct)
            neededFields = StageScanParams.varProps;
            
            obj = StageScanParams;
            for fieldIndex = 1 : length(neededFields)
                fieldName = neededFields{fieldIndex};
                if isfield(inputStruct, fieldName)
                    value = inputStruct.(fieldName);
                    if ischar(value)
                        obj.(fieldName) = value(:)';    % a string stays a row (transposing made 'xyz' a column)
                    elseif isempty(value) && iscell(obj.(fieldName))
                        obj.(fieldName) = {};           % jsondecode turns {} into []
                    else
                        obj.(fieldName) = value';
                    end
                elseif ~any(strcmp(fieldName, StageScanParams.LATER_FIELDS))
                    warning('Couldn''t find field "%s" in struct. Using default value.', fieldName);
                end
            end
            obj.reverse = logical(obj.reverse(:)');
            obj.zigzag = logical(obj.zigzag);
            obj.returnToStart = logical(obj.returnToStart);
            if isempty(obj.scanOrder) || isnumeric(obj.scanOrder); obj.scanOrder = ''; end
        end
    end
end