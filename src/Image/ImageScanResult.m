classdef ImageScanResult < Savable & EventSender & EventListener
    %IMAGESCANRESULT stores the last scan-image
    %   it is savable, so the last image can be saved and opened
    
    % Hopefully, no other figure will need such functionality. Otherwise,
    % this class should be abstracted, and a child class will need to
    % implement it for scan result.
    
    properties (GetAccess = public, SetAccess = private)
        mData = [];     % double. Scan results
        mDimNumber = 2; % intgeger. Number of dimensions of scan. Has default value;
        mFirstAxis      % vector
        mSecondAxis     % vector or point
        mStageName      % string
        mAxesString     % string
        mLabelBot       % string
        mLabelLeft      % string
        mHist = []      % double. Scan lifetime results.
        mThirdAxis = [] % vector. The slice axis of a 3D scan (mAxesString(3)); [] otherwise
        mSterr = []     % double, like mData: each point's standard error, when the detector gives one (photodiode)
    end

    properties
        % 3D scans: which slice is shown, and how
        sliceIndex = 1          % integer. Slice shown (index into mThirdAxis)
        sliceFollow = true      % logical. Show the slice being scanned (set false when the user picks one)
        sliceView = 1           % integer. Index in SLICE_VIEW_OPTIONS
    end

    properties (Access = private, Transient)
        lastLabel = ''          % colour-bar label of the image on screen (see updateImageInPlace)
        lastAxesLabels = {}     % its axis labels
    end

    properties (SetAccess = private, Transient, Hidden)
        lastUpdateFull = true   % the last update redrew everything (which also resets the mouse mode)
    end
    
    properties (Access = private) % For plotting over Image
        gAxes
        gLimits = [];           % graphic handle; stores limits drawn onto image (for zoom)
        crosshairs = struct;    % struct of graphic handles, for the elements of the crosshairs\arrow for current position
        cursor                  % datacursor handle. Stores the information about recent mouse actions
    end
    
    properties
        plotStyle               % integer. Index for value in obj.PLOT_STYLE_OPTIONS (either _1D or _2D)
        colormapType            % integer. Index for value in COLORMAP_OPTIONS
        colormapLimits = [0 1]; % 2x1 double
        colormapAuto = true;    % logical
        cursorType              % integer. Index for value in CURSOR_OPTIONS
        contrastType            % integer. Index for value in CONTRAST_OPTIONS
        imageType               % integer. Index for value in IMAGE_OPTIONS
    end
    
    properties (Constant)
        NAME = 'imageScanResult'
        IMAGE_FILE_SUFFIX = 'png'
        FIGURE_FILE_SUFFIX = 'fig'
        
        % Events
        EVENT_IMAGE_UPDATED = 'imageUpdated';

        % Figure Options
        PLOT_STYLE_OPTIONS_1D = {'Normal', 'Square'};
        PLOT_STYLE_OPTIONS_2D = {'Normal', 'Equal', 'Square'};
        COLORMAP_OPTIONS = {'Pink', 'Parula', 'HSV', 'Hot', 'Cool', ...
            'Spring', 'Summer', 'Autumn', 'Winter', ...
            'Gray', 'Bone', 'Copper', 'Lines'};
        CURSOR_OPTIONS = {'Marker', 'Zoom', 'Location'};
        CONTRAST_OPTIONS = {'Difference', 'Ratio', 'Without MW', 'With MW'};
        IMAGE_OPTIONS = {'kcps', 'Lifetime'};
        SLICE_VIEW_OPTIONS = {'One slice', 'Maximum over slices', 'Mean over slices'};
    end
    
    methods
        function obj = ImageScanResult
            obj@Savable(ImageScanResult.NAME);
            obj@EventSender(ImageScanResult.NAME);
            obj@EventListener({SaveLoadCatImage.NAME, StageScanner.NAME});
            
            obj.plotStyle = 2;      % Equal
            obj.colormapType = 2;   % Parula
            obj.cursorType = 1;     % Marker
            obj.contrastType = 1;   % Difference
            obj.imageType = 1;      % kcps
        end
        
        %% Updating Image
        function update(obj, newStruct)
            % This function could be called in two cases:
            % 1. When new data arrives. We then have newStruct.
            % 2. When ViewImageResultImage starts, and wants data to plot.
            %    Then, no structExtra will be available, but if there is
            %    available data, we still want to plot it.
            if isempty(obj.gAxes) || ~ishandle(obj.gAxes)
                return
            end
            if exist('newStruct', 'var')
                % Get EVERYTHING from the struct
                obj.mData = newStruct.mData;
                obj.mDimNumber = newStruct.mDimNumber;
                obj.mFirstAxis = newStruct.mFirstAxis;
                obj.mSecondAxis = newStruct.mSecondAxis;
                obj.mLabelBot = newStruct.mLabelBot;
                obj.mLabelLeft = newStruct.mLabelLeft;
                
                if isfield(newStruct, 'mHist') % Only available in new scans in setup 4 (29.9.20)
                    obj.mHist = newStruct.mHist;
                else
                    obj.mHist = [];
                end

                % Per-point standard errors (2026-11), when the detector gives them
                if isfield(newStruct, 'mSterr') && isequal(size(newStruct.mSterr), size(obj.mData))
                    obj.mSterr = newStruct.mSterr;
                else
                    obj.mSterr = [];
                end

                % 3D scans (2026-11): the slice axis, and the slice being scanned
                if isfield(newStruct, 'mThirdAxis') && obj.mDimNumber == 3
                    obj.mThirdAxis = newStruct.mThirdAxis;
                else
                    obj.mThirdAxis = [];
                end
                if obj.mDimNumber == 3
                    if isfield(newStruct, 'scanningSlice') && ~isempty(newStruct.scanningSlice) && obj.sliceFollow
                        obj.sliceIndex = newStruct.scanningSlice;
                    elseif isfield(newStruct, 'sliceIndex') && ~isempty(newStruct.sliceIndex)
                        obj.sliceIndex = newStruct.sliceIndex;     % a saved file: the slice it was saved at
                    end
                    obj.sliceIndex = min(max(1, round(obj.sliceIndex)), max(1, numel(obj.mThirdAxis)));
                end

                % In a previous version we did not save mStageName and
                % mAxesString. For compatability, we need:
                if isfield(newStruct, 'mStageName')
                    obj.mStageName = newStruct.mStageName;
                    if isfield(newStruct, 'mAxesString')
                        obj.mAxesString = newStruct.mAxesString;
                    else
                        obj.sendWarning('This should not have happenned...')
                    end
                else
                    % We are loading an old version, so maybe we can get
                    % missing parameters from current system
                    obj.sendWarning(['Save file does not contain stage name.\n', ...
                        'We keep using the old stage, and hope for the best.'])
                end
            elseif ~obj.isDataAvailable
                % No data is available, neither externally nor internally.
                return
            end
            
            [data, label] = obj.getDataForPlotting;

            % We need to calculate this before sending event (so as to update the header):
            if obj.colormapAuto
                obj.colormapLimits = obj.calcColormapLimits(obj.dataForColormapLimits(data));
            end

            % Same image as on screen (a scan's next line, another slice, a re-scan of the
            % same area): only swap the pixels, so a zoom or pan the user made stays.
            if obj.updateImageInPlace(data, label)
                obj.lastUpdateFull = false;
                obj.styleAxes();
                if obj.mDimNumber > 1
                    obj.applyColormapLimits();
                end
                obj.sendEventImageUpdated();
                return
            end
            obj.lastUpdateFull = true;

            switch obj.mDimNumber
                case {2, 3}     % 3D: one slice (or a projection) is shown as a 2D image
                    AxesHelper.fill(obj.gAxes, data, 2, obj.mFirstAxis, obj.mSecondAxis, obj.mLabelBot, obj.mLabelLeft, [], label);
                case 1
                    % with error bars when the detector gave them (fitStyle's look: same hue, no caps)
                    yLabel = label;
                    if strcmp(label, Spcm.signalUnits())
                        yLabel = sprintf('signal (%s)', label);     % not just 'V'
                    end
                    AxesHelper.fill(obj.gAxes, data, obj.mDimNumber, obj.mFirstAxis, obj.mSecondAxis, obj.mLabelBot, yLabel, ...
                        obj.errorsForPlotting(data));
                    ExpPlotStyle.apply(obj.gAxes, 1, obj.mFirstAxis);   % the experiment plots' look
            end
            obj.lastLabel = label;
            obj.lastAxesLabels = {obj.mLabelBot, obj.mLabelLeft};
            obj.styleAxes();
            MainGuiStyle.makeInteractive(obj.gAxes);    % (a full redraw replaces the axes toolbar)
            obj.sendEventImageUpdated();
            obj.imagePostProcessing;
        end

        function tf = updateImageInPlace(obj, data, label)
            % Put new pixels (points) into the image (trace) already drawn, if it is the
            % same one: same size, same axes, same labels
            tf = false;
            if isempty(obj.gAxes) || ~isgraphics(obj.gAxes) || ~strcmp(label, obj.lastLabel) ...
                    || ~isequal(obj.lastAxesLabels, {obj.mLabelBot, obj.mLabelLeft})
                return
            end
            if obj.mDimNumber == 1
                x = obj.mFirstAxis(:)';
                series = findobj(obj.gAxes, 'Type', 'errorbar', '-or', 'Type', 'line');
                series = series(arrayfun(@(g) numel(g.XData) == numel(x), series));
                if numel(series) ~= 1 || numel(x) < 3 || numel(data) ~= numel(x) ...
                        || any(abs(series.XData - x) > 1e-9)
                    return
                end
                err = obj.errorsForPlotting(data);
                if isempty(err) == strcmp(series.Type, 'errorbar'); return; end    % bars appeared / went: redraw
                series.YData = data(:)';
                if ~isempty(err)
                    set(series, 'YNegativeDelta', err(:)', 'YPositiveDelta', err(:)');
                end
                tf = true;
                return
            end
            if ~any(obj.mDimNumber == [2 3]); return; end
            img = findobj(obj.gAxes, 'Type', 'image');
            if numel(img) ~= 1 || ~isequal(size(img.CData), size(data)); return; end
            x = img.XData; y = img.YData;
            if isempty(obj.mFirstAxis) || isempty(obj.mSecondAxis) ...
                    || abs(x(1) - obj.mFirstAxis(1)) > 1e-9 || abs(x(end) - obj.mFirstAxis(end)) > 1e-9 ...
                    || abs(y(1) - obj.mSecondAxis(1)) > 1e-9 || abs(y(end) - obj.mSecondAxis(end)) > 1e-9
                return
            end
            img.CData = data;
            tf = true;
        end

        function applyColormapLimits(obj)
            try
                clim(obj.gAxes, obj.colormapLimits);
            catch
                clim(obj.gAxes, 'auto');
            end
        end

        function showSlice(obj, k)
            % 3D scans: show slice k (the user picked it, so stop following the scan)
            if obj.mDimNumber ~= 3 || isempty(obj.mThirdAxis); return; end
            obj.sliceIndex = min(max(1, round(k)), numel(obj.mThirdAxis));
            obj.sliceFollow = false;
            obj.update;
        end

        function setSliceView(obj, viewIndex)
            obj.sliceView = viewIndex;
            obj.update;
        end

        function styleAxes(obj)
            % The main window's look for the scan image: secondary ink, ticks out
            if isempty(obj.gAxes) || ~isgraphics(obj.gAxes); return; end
            C = ExpPlotStyle.col();
            set(obj.gAxes, 'FontSize', 9, 'XColor', C.ink2, 'YColor', C.ink2, 'TickDir', 'out', 'Box', 'off');
            try
                cb = obj.gAxes.Colorbar;
                if ~isempty(cb)
                    set(cb, 'Color', C.ink2, 'FontSize', 9, 'TickDirection', 'out', 'Box', 'off');
                end
            catch
            end
            if obj.mDimNumber == 3 && ~isempty(obj.mThirdAxis)
                t = sprintf('%s = %s %s   (slice %d of %d)', obj.mAxesString(3), ...
                    MainGuiStyle.num(obj.mThirdAxis(obj.sliceIndex), 3), MainGuiStyle.micron, ...
                    obj.sliceIndex, numel(obj.mThirdAxis));
                if obj.sliceView > 1
                    t = sprintf('%s over %d %s slices', obj.SLICE_VIEW_OPTIONS{obj.sliceView}, ...
                        numel(obj.mThirdAxis), obj.mAxesString(3));
                end
                title(obj.gAxes, t, 'FontWeight', 'normal', 'FontSize', 9, 'Color', C.ink);
            else
                title(obj.gAxes, '');
            end
        end

        function data = dataForColormapLimits(obj, shown)
            % Auto colour limits. A 3D scan uses all its slices, so that stepping through
            % the slices compares like with like.
            data = shown;
            if obj.mDimNumber == 3 && obj.sliceView == 1 && obj.imageType == 1
                data = obj.contrastOf(obj.mData, 3);
            end
        end
        
        function imagePostProcessing(obj)
            %%%% After plotting is done, we can make some slight changes
            
            % Plot style
            switch obj.mDimNumber
                case 1
                    styleOptions = obj.PLOT_STYLE_OPTIONS_1D;
                    if obj.plotStyle == 3
                        % Not available anymore
                        obj.plotStyle = 2;
                    end
                case {2, 3}
                    styleOptions = obj.PLOT_STYLE_OPTIONS_2D;
                    % Switching from 1D to 2D will get us in Equal axes.
                    % Not great, but to worth fixing, for now.
            end
            
            styleName = lower(styleOptions{obj.plotStyle});
            axis(obj.gAxes, styleName)
            
            % Colormap
            colormapName = obj.COLORMAP_OPTIONS{obj.colormapType};
            colormap(obj.gAxes, colormapName)
            
            data = obj.getDataForPlotting;

            if obj.colormapAuto
                obj.colormapLimits = obj.calcColormapLimits(obj.dataForColormapLimits(data));
            end
            
            try
                caxis(obj.gAxes, obj.colormapLimits);
            catch % There is a problem with the limits we tried to enforce
                caxis(obj.gAxes, 'auto');
                obj.sendWarning(['Assigned colorbar limits were problematic.\n' ...
                    'Limits were set automatically, instead.']);
            end
            
            % Get current stage position
            try     % This will not happen without an operating stage
                stage = getObjByName(obj.mStageName);
                    if isempty(stage); throwBaseObjException(obj.mStageName); end
                pos = stage.Pos(obj.imageAxesString);
                % Draw
                if isempty(obj.cursor)
                    % If current cursor does not exist (for some reason)
                    warning('This shouldn''t have happenned')
                    fig = ancestor(obj.gAxes, 'figure');
                    obj.cursor = datacursormode(fig);
                    set(obj.cursor, 'UpdateFcn', @obj.cursorMarkerDisplay);
                end
                gLimitsVector = axis(obj.gAxes);   % vector of [x_min x_max y_min y_max]
                obj.drawCrosshairs(gLimitsVector, pos)
            catch err
                obj.sendWarning('Unable to set crosshairs');
                err2warning(err)
            end
            
            obj.clearCursorData;    % left outside of updateDataCursor(obj), so that zoom bar is not deleted
            obj.updateDataCursor;
        end
        
        function [data, label] = getDataForPlotting(obj)
            switch obj.imageType
                case 1 % Normal intensity image (kcps from a counter, V from the photodiode)
                    if obj.mDimNumber == 3
                        [data, label] = obj.contrastOf(obj.sliceOrProjection, 2);
                    else
                        [data, label] = obj.contrastOf(obj.mData, obj.mDimNumber);
                    end
                case 2 % Lifetime image
                    label = 'ns';
                    if isempty(obj.mHist) || obj.mDimNumber == 3 % Just incase it is some leftover, switch back to normal kcps image.
                        obj.imageType = 1;
                        [data, label] = getDataForPlotting(obj);
                        return
                    end
                    % Uses RLD method from https://en.wikipedia.org/wiki/Fluorescence-lifetime_imaging_microscopy#Analysis                    
                    
                    % Find maximum and get rid of the data before it (get
                    % rid of IRF)
                    switch obj.mDimNumber
                        case 2
                            [~, maxI] = max(obj.mHist(:,:,:), [], 3);
                            maxI = round(mean(mean(maxI)));
                        case 1
                            [~, maxI] = max(obj.mHist(:,:), [], 2);
                            maxI = round(mean(maxI));
                    end
                    
                    spcm = getObjByName(Spcm.NAME);
                    % Total time in ps:
                    T = spcm.binWidth*(spcm.nBins - maxI);
                    
                    dt = (T/2) / 1000; % Half of the time, in ns.
                    midI = round((spcm.nBins + maxI)/2);
                    
                    switch obj.mDimNumber
                        case 2
                            D0 = sum(obj.mHist(:,:,1:midI), 3) * dt;
                            D1 = sum(obj.mHist(:,:,midI+1:end), 3) * dt;
                        case 1
                            D0 = sum(obj.mHist(:,1:midI), 2) * dt;
                            D1 = sum(obj.mHist(:,midI+1:end), 2) * dt;
                    end
                    data = dt ./ log(D0./D1);
            end
        end
        
        function [data, label] = contrastOf(obj, raw, nDim)
            % The quantity shown, from raw data whose dimension nDim+1 holds [without MW,
            % with MW] in a MW contrast scan (and is a singleton otherwise).
            units = Spcm.signalUnits();
            label = units;
            if size(raw, nDim + 1) ~= 2
                data = raw;
                return
            end
            idx = repmat({':'}, 1, nDim);
            withoutMW = raw(idx{:}, 1);
            withMW = raw(idx{:}, 2);
            switch obj.contrastType
                case 1 % Difference
                    data = withoutMW - withMW;
                    label = [units ' difference'];
                case 2 % Ratio
                    data = withMW./withoutMW;
                    label = 'ratio (with MW / without MW)';     % dimensionless
                case 3 % Without MW
                    data = withoutMW;
                otherwise % With MW
                    data = withMW;
            end
        end

        function err = errorsForPlotting(obj, data)
            % Error bars for a 1D scan: each point's standard error, when the detector gave
            % them (the photodiode does; photon counters do not) and the plot shows the
            % measured values themselves (not a MW-contrast difference or ratio). [] = none.
            err = [];
            if obj.mDimNumber ~= 1 || obj.imageType ~= 1 || size(obj.mData, 2) ~= 1 ...
                    || numel(obj.mSterr) ~= numel(data) || ~any(isfinite(obj.mSterr(:)))
                return
            end
            err = obj.mSterr(:);
        end

        function e = errorAt(obj, row, col)
            % A 2D pixel's standard error (NaN when unknown), for its data tip
            e = NaN;
            if obj.mDimNumber == 2 && obj.imageType == 1 && size(obj.mData, 3) == 1 ...
                    && isequal(size(obj.mSterr), size(obj.mData)) && row <= size(obj.mSterr, 1) && col <= size(obj.mSterr, 2)
                e = obj.mSterr(row, col);
            end
        end

        function raw = sliceOrProjection(obj)
            % 3D scans: the slice shown, or the maximum / mean over all slices, as
            % [rows x cols (x 2 with MW contrast)]
            vol = obj.mData;
            k = min(max(1, obj.sliceIndex), size(vol, 3));
            switch obj.sliceView
                case 2
                    raw = max(vol, [], 3, 'omitnan');
                case 3
                    raw = mean(vol, 3, 'omitnan');
                otherwise
                    raw = vol(:, :, k, :);
            end
            raw = reshape(raw, size(vol, 1), size(vol, 2), []);
        end

        function s = imageAxesString(obj)
            % The axes drawn in the image: the scan axes, without the slice axis of a 3D scan
            s = obj.mAxesString;
            if obj.mDimNumber == 3 && numel(s) >= 2
                s = s(1:2);
            end
        end

        %%%% Data Cursor methods %%%%
        function clearCursorData(obj)
            % Clear the cursor data between passing from one cursor type to
            % another. To avoid collision.
            
            % Stop zoom if active
            % (uses undocumented feature)
            global GETRECT_H1
            if ~isempty(GETRECT_H1) && ishghandle(GETRECT_H1)
                set(GETRECT_H1, 'UserData', 'Completed');
            end
            
            if isempty(obj.cursor)
                % If current cursor does not exist (for some reason)
                warning('This shouldn''t have happenned')
                fig = ancestor(obj.gAxes, 'figure');
                obj.cursor = datacursormode(fig);
            end
            obj.cursor.removeAllDataCursors;
            set(obj.cursor, 'Enable', 'off');
            
            % Remove drawn limits, if exists
            if ~isempty(obj.gLimits)
                delete(obj.gLimits);
            end
            
            % Disable button press
            set([obj.gAxes; obj.gAxes.Children], 'ButtonDownFcn', '');
        end
        
        function drawRectangle(obj, pos)
            % Draw rectangle
            if HandleHelper.isType(obj.gLimits, 'rectangle')
                obj.gLimits.Position = pos;
            else
                obj.gLimits = rectangle(obj.gAxes, ...
                    'Position', pos, ...
                    'EdgeColor', 'g', ...
                    'LineWidth', 1, ...
                    'LineStyle', '-.', ...
                    'HitTest', 'Off');
            end
        end
        
        function drawLimitBar(obj,pos)
            % Draw limit bar
            if ~isempty(obj.gLimits)
                delete(obj.gLimits);
            end
            dx = pos(3)/2;
            xPos = pos(1) + dx;         % the center of the bar
            yPos = pos(2) + pos(4)/2;	% position of bar in the middle of selected rectangle
            
            hold(obj.gAxes,'on');
            obj.gLimits = errorbar(obj.gAxes, ...
                xPos,yPos, ...
                dx, 'horizontal', ...
                'Color', 'r', ...
                'LineWidth', 1);
            hold(obj.gAxes,'off');
        end
        
        function updateDataCursor(obj, action)
            % Needs to happen only when everything else finished, whether
            % by scanning or loading
            %
            % action - string. One of obj.CURSOR_OPTIONS

            if isempty(obj.cursor)
                % Nothing to do here
                return
            end
            
            if ~exist('action', 'var')
                actionIndex = 1;    % Marker is default
            else
                actionIndex = find(strcmp(action, obj.CURSOR_OPTIONS));
            end
            switch actionIndex
                case 1      % Inspect: zoom, pan, data tips (hover, or click to pin)
                    % No figure mode: the data-cursor MODE (on here before) turns zoom and
                    % pan off in every plot of the window. Data tips come from the plots'
                    % own interactions and the axes toolbar, with this text.
                    set(obj.cursor, 'Enable', 'off')
                    obj.cursor.UpdateFcn = @obj.cursorMarkerDisplay;
                    set([obj.gAxes; obj.gAxes.Children], 'ButtonDownFcn', '');
                    MainGuiStyle.makeInteractive(obj.gAxes);
                case 2      % Create a rectangle on the selected area, and update the GUI scanParams's min and max values accordingly
                    obj.clearCursorData;
                    obj.updataDataByZoom;
                case 3      % Move the stage to selected location
                    set(obj.cursor, 'Enable', 'off')
                    try; disableDefaultInteractivity(obj.gAxes); catch; end   % a click moves the stage, not a data tip
                    set([obj.gAxes; obj.gAxes.Children], 'ButtonDownFcn', @obj.setLocationFromCursor);
            end
        end
        
        function txt = cursorMarkerDisplay(obj, ~, event_obj)
            % Displays the location of the cursor on the plot and the kcps
            % (the color level from the colormap)
            
            % Displays the location of the cursor on the plot and the kcps
            % (the color level from the colormap)
            
            imageScanResult = getObjByName(ImageScanResult.NAME);
            if isempty(imageScanResult); throwBaseObjException(ImageScanResult.NAME); end

            % The data-tip text is set for the whole window: a tip on another plot (the
            % photodiode trace) gets plain x / y
            target = ancestor(get(event_obj, 'Target'), 'axes');
            if ~isequal(target, obj.gAxes)
                pos = event_obj.Position;
                txt = {sprintf('x: %.4g', pos(1)), sprintf('y: %.4g', pos(2))};
                return
            end

            if ~imageScanResult.isDataAvailable
                txt = '';
                EventStation.anonymousWarning('Image is empty');
                return
            end
            
            % Get scanning parameters
            dim = imageScanResult.mDimNumber;
            scanAxes = imageScanResult.imageAxesString;
            firstAxis = scanAxes(1);	% char (either 'x', 'y' or 'z')
            secondAxis = scanAxes(end); % char (either 'x', 'y' or 'z')
            
            % Customizes text of data tips
            [data, label] = imageScanResult.getDataForPlotting;
            pos = event_obj.Position;

            % %.4g, not %.1f: a photodiode reads tens of mV, and %.1f printed 0.05 V as "0.1".
            % "± error" when the point's standard error is known
            pm = @(e) BooleanHelper.ifTrueElse(isfinite(e), sprintf(' %s %.2g', char(177), e), '');
            if dim == 1
                err = imageScanResult.errorsForPlotting(data);
                e = NaN;
                if ~isempty(err)
                    [~, k] = min(abs(imageScanResult.mFirstAxis - pos(1)));
                    e = err(k);
                end
                txt = sprintf('%s = %.3f\n%.4g%s %s', firstAxis, pos(1), pos(2), pm(e), label);
            else
                % The pixel under the tip (getimage needed the Image Processing Toolbox)
                [~, col] = min(abs(imageScanResult.mFirstAxis - pos(1)));
                [~, row] = min(abs(imageScanResult.mSecondAxis - pos(2)));
                level = NaN;
                if row <= size(data, 1) && col <= size(data, 2)
                    level = data(row, col);
                end
                txt = {sprintf('(%s,%s) = (%.3f, %.3f)\n%.4g%s %s', ...
                    firstAxis, secondAxis,...
                    pos(1), pos(2), level, pm(imageScanResult.errorAt(row, col)), label)};
                if dim == 3 && imageScanResult.sliceView == 1 && ~isempty(imageScanResult.mThirdAxis)
                    txt = {sprintf('%s\n%s = %.3f (slice %d)', txt{1}, imageScanResult.mAxesString(3), ...
                        imageScanResult.mThirdAxis(imageScanResult.sliceIndex), imageScanResult.sliceIndex)};
                end
            end
        end
        
        function updataDataByZoom(obj)
            % Draw the rectangle on the selected area on the plot,
            % and update the GUI with the max and min values
            
            if ~obj.isDataAvailable      % nothing to zoom to
                EventStation.anonymousWarning('Image is empty');
                return
            end
            
            % "try" getting user input
            warning('off','all');
            rect = getrect(obj.gAxes);
            warning('on','all');
            if rect(3) == 0; return; end  % Selection has no width. No use in continuing
            
            % Draw, according to the dimensions of the image
            dim = min(obj.mDimNumber, 2);   % a 3D scan shows a 2D slice
            switch dim
                case 1
                    obj.drawLimitBar(rect);
                case 2
                    if rect(4) == 0     % Selection has no height
                        return
                    end
                    obj.drawRectangle(rect);
                otherwise
                    EventStation.anonymousError('Can''t fetch limits in higher dimensions');
            end
            
            stage = getObjByName(obj.mStageName);
            if isempty(stage)
                EventStation.anonymousError('Can''t set new scan parameters, since there is no relevant stage');
            end
            assert(stage.isScannable)
            
            for i = 1:dim
                axisIndex = stage.getAxis(obj.imageAxesString(i));
                
                stage.scanParams.from(axisIndex) = rect(i);            % rect(1)==horizontal position, rect(2)==vertical position
                stage.scanParams.to(axisIndex) = rect(i)+rect(i+2);    % rect(3)==width;	rect(4)==height
            end
            stage.sendEventScanParamsChanged;
        end
        
        function drawCrosshairs(obj, limits, pos)
            % Draws an arrow or a cross on the plot on a selected position
            %       limits: limits of the view in the axes,
            %               [x_min x_max y_min y_max]
            %       pos:    current stage position, [x] or [x y]
            %
            %       Note: this function assumes the input arguments have
            %       been validated.
            
            % Get properties of axes
            xLimits = limits(1:2);
            yLimits = limits(3:4);
            dim = length(pos);
            
            % If 1D, arbitrarily set the position in the middle
            if dim == 1
                pos(2) = (yLimits(2) + yLimits(1))/2;
            end
            
            %%%% Deciding on which area the position is located %%%%
            % Trinary encoding of the position yields:
            %-----------------------%
            %  2,0  |  2,1  |  2,2  |
            %-----------------------|
            %  1,0  |  1,1  |  1,2  |
            %-----------------------|
            %  0,0  |  0,1  |  0,2  |
            %-----------------------%
            % trinary receives:
            %       x, scalar
            %       I, vector representing interval
            % returns:
            %   0, if x < I,
            %   1, if x in I (including edges)
            %   2, if x > I
            quadrant = [trinary(pos(1),xLimits) , trinary(pos(2),yLimits)];
            
            % Delete previous and recreate correct handles
            hold(obj.gAxes,'on'); % Keep image while creating crosshair/arrow
            
            if all(quadrant == [1 1]) % Inside
                if isfield(obj.crosshairs, 'arrowHandle')
                    delete(obj.crosshairs.arrowHandle);
                end
                if ~isfield(obj.crosshairs, 'xLineHandle') || ~isvalid(obj.crosshairs.xLineHandle)
                    obj.crosshairs.xLineHandle = line(obj.gAxes, ...
                        [0 0], [0 0], 'Color', ExpPlotStyle.col().fit, 'LineWidth', 1.2, 'HitTest', 'Off');
                end
                if ~isfield(obj.crosshairs, 'yLineHandle') || ~isvalid(obj.crosshairs.yLineHandle)
                    obj.crosshairs.yLineHandle = line(obj.gAxes, ...
                        [0 0], [0 0], 'Color', ExpPlotStyle.col().fit, 'LineWidth', 1.2, 'HitTest', 'Off');
                end
                
                hold(obj.gAxes,'off'); % Allow picture to change when needed
                
                set(obj.crosshairs.xLineHandle, 'XData', xLimits);
                set(obj.crosshairs.xLineHandle, 'YData', [pos(2) pos(2)]);
                set(obj.crosshairs.yLineHandle, 'XData', [pos(1) pos(1)]);
                set(obj.crosshairs.yLineHandle, 'YData', yLimits);
                
                % If 1D, don't draw horizontal line and change color to black
                if dim == 1
                    delete(obj.crosshairs.xLineHandle);
                    set(obj.crosshairs.yLineHandle, 'Color', 'k') 
                end
                
            else % Outside
                if isfield(obj.crosshairs, 'xLineHandle')
                    delete(obj.crosshairs.xLineHandle);
                end
                if isfield(obj.crosshairs,'yLineHandle')
                    delete(obj.crosshairs.yLineHandle);
                end
                if isfield(obj.crosshairs, 'arrowHandle')
                    if (dim == 1 && ~isa(obj.crosshairs.arrowHandle,'matlab.graphics.primitive.Line')) ...
                            || (dim == 2 && ~isa(obj.crosshairs.arrowHandle,'matlab.graphics.chart.primitive.Quiver'))
                        delete(obj.crosshairs.arrowHandle);
                    end
                end
                if ~isfield(obj.crosshairs, 'arrowHandle') || ~ishandle(obj.crosshairs.arrowHandle)
                    if dim == 1
                        obj.crosshairs.arrowHandle = line(obj.gAxes, ...
                        [0 0], [0 0], 'Color', 'k', 'LineWidth', 1, 'HitTest', 'Off', 'MarkerSize', 10);
                    else
                        obj.crosshairs.arrowHandle = quiver(obj.gAxes, 0, 0, 0, 0, 'Color', ExpPlotStyle.col().fit, ...
                        'LineWidth', 1, 'MaxHeadSize', 10, 'HitTest', 'Off');
                    end
                end
                
                hold(obj.gAxes,'off'); % Allow picture to change when needed
                
                % Draw arrow according to position
                if dim == 1
                    switch quadrant(1)
                        case 0 % left column
                            set(obj.crosshairs.arrowHandle, 'XData', [xLimits(1), xLimits(1)]);
                            set(obj.crosshairs.arrowHandle, 'YData', [pos(2), pos(2)]);
                            set(obj.crosshairs.arrowHandle, 'Marker', '<')
                        case 2 % right column
                            set(obj.crosshairs.arrowHandle, 'XData', [xLimits(2), xLimits(2)]);
                            set(obj.crosshairs.arrowHandle, 'YData', [pos(2), pos(2)]);
                            set(obj.crosshairs.arrowHandle, 'Marker', '>')
                    end
                else
                    row = quadrant(1);
                    switch row
                        case 0 % left column
                            set(obj.crosshairs.arrowHandle, 'XData', xLimits(1)+(xLimits(2)-xLimits(1))/10);
                            set(obj.crosshairs.arrowHandle, 'UData', (xLimits(1)-xLimits(2))/10);
                        case 1 % middle column
                            set(obj.crosshairs.arrowHandle, 'XData', pos(1));
                            set(obj.crosshairs.arrowHandle, 'UData', 0);
                        case 2 % right column
                            set(obj.crosshairs.arrowHandle, 'XData', xLimits(2)+(xLimits(1)-xLimits(2))/10);
                            set(obj.crosshairs.arrowHandle, 'UData', (xLimits(2)-xLimits(1))/10);
                    end
                    
                    column = quadrant(2);
                    switch column
                        case 0 % in bottom row
                            set(obj.crosshairs.arrowHandle, 'YData', yLimits(1)+(yLimits(2)-yLimits(1))/10);
                            set(obj.crosshairs.arrowHandle, 'VData', (yLimits(1)-yLimits(2))/10);
                        case 1 % middle row
                            set(obj.crosshairs.arrowHandle, 'YData', pos(2));
                            set(obj.crosshairs.arrowHandle, 'VData', 0);
                        case 2 % top row
                            set(obj.crosshairs.arrowHandle, 'YData', yLimits(2)+(yLimits(1)-yLimits(2))/10);
                            set(obj.crosshairs.arrowHandle, 'VData', (yLimits(2)-yLimits(1))/10);
                    end
                end
            end
        end
        
        function setLocationFromCursor(obj, ~, ~)
            if isempty(obj.gAxes.Children)
                EventStation.anonymousWarning('Image is empty');
                return
            end
            
            % Get everything we need
            stage = getObjByName(obj.mStageName);
            if isempty(stage); throwBaseObjException(obj.mStageName); end
            pos = obj.gAxes.CurrentPoint(1, 1:min(obj.mDimNumber, 2));
            axesToMove = obj.imageAxesString;
            if obj.mDimNumber == 3 && obj.sliceView == 1 && ~isempty(obj.mThirdAxis)
                % A slice of a 3D scan: go to that slice too
                axesToMove = obj.mAxesString;
                pos(3) = obj.mThirdAxis(obj.sliceIndex);
            end

            % Now move to the corresponding location
            stage.move(axesToMove, pos);    % This should update the datacursor, by relevant event
        end
        
        %% Helper methods
        function sendEventImageUpdated(obj)
            obj.sendQueuedEvent(struct(obj.EVENT_IMAGE_UPDATED, true));
        end
        
        function tf = isDataAvailable(obj)
            tf = ~isempty(obj.mData);
        end
        
        function addGraphicAxes(obj, gAxes)
            % "Setter" for the axes, when they are created in the GUI
            if ~(isgraphics(gAxes) && isvalid(gAxes))
                obj.sendWarning('Graphic Axes were not created. Plotting is unavailable');
                return
            end
            
            obj.gAxes = gAxes;
            fig = ancestor(obj.gAxes,'figure');
            obj.cursor = datacursormode(fig);
            set(obj.cursor,'UpdateFcn',@obj.cursorMarkerDisplay);
        end

        function checkGraphicAxes(obj)
            if exist('obj.gAxes', 'var') && ~(isgraphics(obj.gAxes) && isvalid(obj.gAxes))
                % gAxes are no longer available, so we discard them
                obj.gAxes = [];
                obj.cursor = [];
            end
        end
        
        function fig = copyToFigure(obj, isVisible)
            % isVisible - logical. Should we create the figure as visible,
            % to begin with
            if ~exist('isVisible', 'var') || isVisible
                fig = figure;
            else
                fig = figure('Visible', 'off');
            end

            newAxes = copyobj(obj.gAxes, fig);
            if obj.mDimNumber >= 2
                c = colorbar(newAxes);
                [~, label] = obj.getDataForPlotting;
                xlabel(c, label)
            end
            
            try
                notes = SaveLoad.getInstance(Savable.CATEGORY_IMAGE).mNotes;
                title(notes); % Set the notes as the figure title
            catch err
                EventStation.anonymousWarning(err.message)
            end
        end
        
        function [stageName, axesString] = compatibilityHelper1(obj, struct) %#ok<INUSL>
            % Because of changes in saved struct, we need to [try to]
            % recover these variables
            scanner = getObjByName(StageScanner.NAME);
            stageName = scanner.mStageName;
            switch struct.mDimNumber
                case 1
                    axesString = struct.mLabelBot(1);
                case 2
                    axesString = [struct.mLabelBot(1), struct.mLabelLeft(1)];
            end
        end
        
        %% Saving
        function fullpath = savePlottingImage(obj, folder, filename)
            % Create a new invisible figure, and than save it with a same
            % filename.
            %
            % Input:
            %   folder - string. the path.
            %   filename - string. the file that was saved by the SaveLoad
            %
            % Output:
            % fullpath - the fullpath of the image file that was saved
            
            if isempty(obj.mData); return; end
            
            isVisible = false;
            figureInvis = obj.copyToFigure(isVisible);
            
            filename = PathHelper.removeDotSuffix(filename);
            fullpath = PathHelper.joinToFullPath(folder, filename);
            
            %%% Save image (.png)
            fullPathImage = [fullpath '.' ImageScanResult.IMAGE_FILE_SUFFIX];
            saveas(figureInvis, fullPathImage);
            
            %%% Save figure (.fig)
            % The figure is saved as invisible, but we set its creation
            % function to set it as visible
            set(figureInvis, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')'); % No other methods of specifying the function seemed to work...
            savefig(figureInvis, fullpath)
            
            %%% close the figure
            close(figureInvis);
        end
    end
    
    methods (Static)
        function init
            replaceBaseObject(ImageScanResult);  % in base object map
        end
        
        function limits = calcColormapLimits(data)
            % Calculate auto-limits from data.
            % (NaN: a 3D scan's slices not scanned yet)
            data = data(isfinite(data) & data ~= 0);
            if isempty(data)
                limits = [0 1];     % No information.
                return
            end
            limits = [min(data) max(data)];
            if limits(2) <= limits(1)
                limits = limits(1) + [-0.5 0.5] * max(abs(limits(1)), 1) * 1e-3;
            end
        end
        
        function newStruct = ScanStructToInternal(scanStruct)
            % Converts a struct, as output by StageScanner to the way it is
            % represented within this class (ImageScanResult)
            newStruct.mData = scanStruct.scan;
            newStruct.mDimNumber = scanStruct.dimNumber;
            newStruct.mFirstAxis = scanStruct.getFirstAxis;
            newStruct.mSecondAxis = scanStruct.getSecondAxis;
            newStruct.mStageName = scanStruct.stageName;
            newStruct.mAxesString = scanStruct.axesString;
            newStruct.mLabelBot = scanStruct.botLabel;
            newStruct.mLabelLeft = scanStruct.leftLabel;
            newStruct.mThirdAxis = scanStruct.getThirdAxis;
            newStruct.scanningSlice = [];
            if isprop(scanStruct, 'scanExtraInfo')
                if isfield(scanStruct.scanExtraInfo, 'hist')
                    newStruct.mHist = scanStruct.scanExtraInfo.hist;
                end
                if isfield(scanStruct.scanExtraInfo, 'sliceIndex')
                    newStruct.scanningSlice = scanStruct.scanExtraInfo.sliceIndex;
                end
                if isfield(scanStruct.scanExtraInfo, 'livePixel')
                    newStruct.livePixel = true;
                end
                if isfield(scanStruct.scanExtraInfo, 'sterr')
                    newStruct.mSterr = scanStruct.scanExtraInfo.sterr;
                end
            end
        end
    end
    
    %% Setters
    methods
        function set.mStageName(obj, stageName)
            % The stage might have changed, and we need to start listening to it
            
            % We might not need to do anything, since we are using
            % the same stage.
            if strcmp(stageName, obj.mStageName)
                return
            end
            
            % Check that required stage is available
            newStage = getObjByName(stageName);
            if isempty(newStage)
                obj.sendWarning('New stage is not available. Some things might not work properly.');
            else
                % Switch to new stage
                obj.stopListeningTo(obj.mStageName);
                obj.mStageName = stageName;
                obj.startListeningTo(obj.mStageName);
            end 
        end
    end
    
    %% Getters (for backward compatibility)
    methods 
        function limits = get.colormapLimits(obj)
            if isnumeric(obj.colormapLimits) && (length(obj.colormapLimits) == 2)
                limits = obj.colormapLimits;
            elseif obj.isDataAvailable
                limits = obj.calcColormapLimits(obj.mData);
            else
                limits = [0 1];
            end
        end
        
        function tf = get.colormapAuto(obj)
            % The second condition in this 'if' statement is needed since
            % true ~= 1, but we do want it to be accepted.
            if islogical(obj.colormapAuto) || obj.colormapAuto == 1
                tf = obj.colormapAuto;
            else
                tf = false;
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
            if ~strcmp(category, Savable.CATEGORY_IMAGE) || ~strcmp(type, Savable.TYPE_RESULTS)
                outStruct = NaN;
                return
            end
            
            if isempty(obj.mData)
                outStruct = NaN;
                return
            end
            
            outStruct = struct();
            propNameCell = obj.getAllNonConstProperties();
            for i = 1: length(propNameCell)
                propName = propNameCell{i};
                outStruct.(propName) = obj.(propName);
            end
            
        end
        
        function loadStateFromStruct(obj, savedStruct, category, subCategory)
            % Loads the state from a struct.
            % To support older versoins, always check for a value in the
            % struct before using it. view example in the first line.
            % subCategory - string. could be empty string

            if ~strcmp(category, Savable.CATEGORY_IMAGE); return; end
            if ~any(strcmp(subCategory, {Savable.SUB_CATEGORY_DEFAULT})); return; end
            
            hasChanged = false;     % initialize
            for propNameCell = obj.getAllPropertiesThisClassDefined()
                propName = propNameCell{:};
                if isfield(savedStruct, propName)
                    obj.(propName) = savedStruct.(propName);
                    hasChanged = true;
                    break       % If even one has changed, that's enough
                end
            end
            if hasChanged
                obj.update(savedStruct);
            end
        end
        
        function string = returnReadableString(obj, savedStruct) %#ok<INUSD>
            % Return a readable string to be shown. If this object
            % doesn't need a readable string, make (string = NaN;) or
            % (string = '');
            string = NaN;
        end
    end
    
    %% overridden from EventListener
    methods
        % When events happen, this function jumps.
        % event is the event sent from the EventSender
        function onEvent(obj, event)
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
            
            % Check if event is "SaveLoad wants to save a file" and need to
            % save an image file of the figure
            if strcmp(event.creator.name, SaveLoadCatImage.NAME) ...
                    && isfield(event.extraInfo, SaveLoad.EVENT_SAVE_SUCCESS_LOCAL_TO_FILE) ...
                    && ~isempty(obj.mData)
                
                folder = event.extraInfo.(SaveLoad.EVENT_FOLDER);
                filename = event.extraInfo.(SaveLoad.EVENT_FILENAME);
                obj.savePlottingImage(folder, filename);
            end
            
            % Check if event is "scanner started a new scan"
            if strcmp(event.creator.name, StageScanner.NAME) ...
                    && isfield(event.extraInfo, StageScanner.EVENT_SCAN_STARTED)
                % We might need to change our reference stage
                obj.mStageName = event.creator.mStageName;
                obj.sliceFollow = true;     % a new scan: show the slice being scanned
                % Make sure the cursor is in marker mode
                obj.updateDataCursor;
            end
                    
            
            % Check if event is "scanner has a new line scanned" and need
            % to updated the image
            if strcmp(event.creator.name, StageScanner.NAME) ...
                    && isfield(event.extraInfo, StageScanner.EVENT_SCAN_UPDATED)
                
                extra = event.extraInfo.(StageScanner.EVENT_SCAN_UPDATED);
                % "extra" now points to an object of class EventExtraScanUpdated,
                % but we want it in the in-house format
                extraInternal = obj.ScanStructToInternal(extra);
                obj.update(extraInternal);
                if isfield(extraInternal, 'livePixel') && extraInternal.livePixel
                    drawnow limitrate nocallbacks   % mid-line: clicks wait for the end of the line
                else
                    drawnow
                end
                return % To avoid drawing crosshairs twice
            end
            
            % Check if event is "Stage moved" (not by scanning,
            % which would have been caught by previous condition),
            % and redraw crosshairs
            if isfield(event.extraInfo, ClassStage.EVENT_POSITION_CHANGED) ...
                    && strcmp(event.creator.name, obj.mStageName)
                
                % If scanner is running, and there is no point in updating.
                % Crosshairs would really delay it
                scanner = getObjByName(StageScanner.NAME);
                if scanner.mCurrentlyScanning; return; end
                
                stage = getObjByName(obj.mStageName);
                if isempty(stage)
                    % Probably, there is nothing to draw on. Moving on!
                else
                    physAxes = obj.imageAxesString;
                    pos = stage.Pos(physAxes);
                    limits = axis(obj.gAxes);   % vector of [x_min x_max y_min y_max]
                    obj.drawCrosshairs(limits, pos)
                end
            end
        end
    end

end