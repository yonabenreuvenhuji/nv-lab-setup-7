classdef ViewStagePanelScanParams < GuiComponent & EventListener & EventSender
    %VIEWSTAGEPANELSCANPARAMS The "Scan area" section: what each axis does in a scan.
    %
    %   Axis Scan  From (µm)  To (µm)  Dir  Points  Step (µm)  Centre (µm)  Span (µm)
    %    X   [x]   [-20   ]   [ 20  ]  [→]  [ 41 ]  [  1   ]   [   0     ]  [  40   ]
    %    Y   [x]   ...
    %    Z   [ ]   (greyed: not scanned)                       [  2.5    ]  (grey)
    %   [Centre ← stage]  [Range = centre ± span/2]  [Full soft range]  [Go to centre]
    %   2D · lines along X, stepping Y · 41 × 41 = 1,681 pixels · at least 25 s
    %
    % WHAT USED TO BE "Fix" / "Fixed". "Fix" ticked meant "this axis is NOT scanned", and
    % "Fixed" was the position such an axis is held at -- which is ALSO the centre that "Scan
    % Around" built the range from. So the two were not the same thing, and greying out
    % "Fixed" when "Fix" was off hid a number that was still in use. Now:
    %   Scan     ticked = this axis is scanned (the opposite of the old Fix).
    %   Centre   always editable: where an unscanned axis stays during the scan, and the
    %            middle of the range when you set a Span. (Stored as scanParams.fixedPos.)
    %   From / To / Dir / Points / Step / Span only matter for a scanned axis, so they are
    %            greyed when Scan is off (they keep their values).
    %   Step     editable: changing it changes Points to fit the range.
    %   Span     editable: From/To = Centre ± Span/2.
    %   Dir      → = from "From" to "To"; ← = the other way. The image is the same either
    %            way (it is always stored From -> To); only the order the stage moves in
    %            changes, which matters with backlash or drift.

    properties
        cbxScan             % 1xN checkbox: scan this axis
        edtFrom             % 1xN
        edtTo               % 1xN
        btnDir              % 1xN toggle: direction
        edtNumPoints        % 1xN
        edtStep             % 1xN
        edtCentre           % 1xN (scanParams.fixedPos)
        edtSpan             % 1xN
        tvSummary           % one line: what the scan will do
        btnCentreHere
        btnRangeAround
        btnFullRange
        btnGoCentre

        stageName   % string. the name of the relevant stage
        stageAxes   % string. valid axes of the stage. example: "zx"
    end

    properties (Constant)
        COL_WIDTHS = [20 38 66 66 30 48 62 70 58]
    end

    methods
        function obj = ViewStagePanelScanParams(parent, controller, stage)
            obj@GuiComponent(parent, controller);
            obj@EventListener(stage.name);
            obj@EventSender(sprintf('%s%s', stage.name, ' _ panel Scan params'));
            obj.stageName = stage.name;
            obj.stageAxes = stage.availableAxes;
            n = length(stage.availableAxes);
            S = MainGuiStyle;
            um = MainGuiStyle.micron;

            [panel, inner] = MainGuiStyle.section(parent.component, 'Scan area');
            obj.component = panel;

            gridRow = MainGuiStyle.hbox(inner, 0);
            grid = MainGuiStyle.grid(gridRow, 4);
            MainGuiStyle.empty(gridRow);

            % columns: axis | scan | from | to | dir | points | step | centre | span
            MainGuiStyle.empty(grid);
            for i = 1:n; MainGuiStyle.axisTag(grid, obj.stageAxes(i)); end

            MainGuiStyle.header(grid, 'Scan');
            obj.cbxScan = gobjects(1, n);
            for i = 1:n
                box = MainGuiStyle.hbox(grid, 0);
                MainGuiStyle.empty(box);
                obj.cbxScan(i) = MainGuiStyle.checkbox(box, '', false, 'TooltipString', ...
                    sprintf('Scan along %s. Off: %s stays at its Centre during the scan.', upper(obj.stageAxes(i)), upper(obj.stageAxes(i))));
                box.Widths = [-1 18];
            end

            obj.edtFrom = obj.column(grid, sprintf('From (%s)', um), n, 'First point of the range');
            obj.edtTo = obj.column(grid, sprintf('To (%s)', um), n, 'Last point of the range');

            MainGuiStyle.header(grid, 'Dir');
            obj.btnDir = gobjects(1, n);
            for i = 1:n
                obj.btnDir(i) = MainGuiStyle.button(grid, char(8594), 'normal', 'FontSize', 11);
            end

            obj.edtNumPoints = obj.column(grid, 'Points', n, 'Number of points along this axis');
            obj.edtStep = obj.column(grid, sprintf('Step (%s)', um), n, ...
                'Distance between points. Typing a step changes Points to fit From-To.');
            obj.edtCentre = obj.column(grid, sprintf('Centre (%s)', um), n, ...
                'Where this axis stays when it is not scanned, and the middle of the range when you type a Span');
            obj.edtSpan = obj.column(grid, sprintf('Span (%s)', um), n, ...
                'Length of the range. Typing a span sets From/To = Centre ± Span/2.');

            set(grid, 'Widths', obj.COL_WIDTHS, 'Heights', [S.HEAD_H, S.ROW_H * ones(1, n)]);
            gridRow.Widths = [sum(obj.COL_WIDTHS) + 4*numel(obj.COL_WIDTHS), -1];

            buttons = MainGuiStyle.hbox(inner, S.GAP);
            obj.btnCentreHere = MainGuiStyle.button(buttons, [sprintf('Centre %s stage', char(8592))], 'normal', ...
                'TooltipString', 'Set every Centre to where the stage is now');
            obj.btnRangeAround = MainGuiStyle.button(buttons, sprintf('Range = centre %s span/2', char(177)), 'normal', ...
                'TooltipString', 'For every scanned axis: From/To = Centre ± Span/2');
            obj.btnFullRange = MainGuiStyle.button(buttons, 'Full soft range', 'normal', ...
                'TooltipString', 'For every scanned axis: From/To = the soft limits');
            obj.btnGoCentre = MainGuiStyle.button(buttons, 'Go to centre', 'normal', ...
                'TooltipString', 'Move the stage to the Centre of every axis');
            buttons.Widths = [-1.1 -1.4 -1 -0.9];

            obj.tvSummary = MainGuiStyle.hint(inner, '');

            inner.Heights = [S.HEAD_H + n*S.ROW_H + 4*n, S.ROW_H, 16];
            obj.height = MainGuiStyle.sectionHeight(MainGuiStyle.rowsHeight(inner.Heights));
            obj.width = sum(obj.COL_WIDTHS) + 60;

            %%%% Callbacks %%%%
            for i = 1:n
                set(obj.cbxScan(i), 'Callback', @(h,e) obj.cbxScanCallback(i));
                set(obj.edtFrom(i), 'Callback', @(h,e) obj.edtFromChangedCallback(i));
                set(obj.edtTo(i), 'Callback', @(h,e) obj.edtToChangedCallback(i));
                set(obj.btnDir(i), 'Callback', @(h,e) obj.btnDirCallback(i));
                set(obj.edtNumPoints(i), 'Callback', @(h,e) obj.edtNumPointsChangedCallback(i));
                set(obj.edtStep(i), 'Callback', @(h,e) obj.edtStepCallback(i));
                set(obj.edtCentre(i), 'Callback', @(h,e) obj.edtFixedChangedCallback(i));
                set(obj.edtSpan(i), 'Callback', @(h,e) obj.edtSpanCallback(i));
            end
            obj.btnCentreHere.Callback = @(h,e) obj.btnCentreHereCallback();
            obj.btnRangeAround.Callback = @(h,e) obj.btnRangeAroundCallback();
            obj.btnFullRange.Callback = @(h,e) obj.btnFullRangeCallback();
            obj.btnGoCentre.Callback = @(h,e) obj.btnGoCentreCallback();

            obj.refresh();
        end

        function refresh(obj)
            [stage, scanParams] = obj.getStageAndParams();
            um = MainGuiStyle.micron;
            for i = 1 : length(obj.stageAxes)
                k = ClassStage.getAxis(obj.stageAxes(i));     % index in scan params
                isScanned = ~scanParams.isFixed(k);
                obj.cbxScan(i).Value = isScanned;
                obj.edtFrom(i).String = MainGuiStyle.num(scanParams.from(k));
                obj.edtTo(i).String = MainGuiStyle.num(scanParams.to(k));
                obj.edtNumPoints(i).String = sprintf('%d', scanParams.numPoints(k));
                obj.edtStep(i).String = MainGuiStyle.num(obj.stepOf(scanParams, k), 4);
                obj.edtCentre(i).String = MainGuiStyle.num(scanParams.fixedPos(k));
                obj.edtSpan(i).String = MainGuiStyle.num(abs(scanParams.to(k) - scanParams.from(k)));

                reversed = numel(scanParams.reverse) >= k && scanParams.reverse(k);
                if reversed
                    set(obj.btnDir(i), 'String', char(8592), 'TooltipString', ...
                        sprintf('%s is scanned from To to From (%s %s %s). Click: From to To.', upper(obj.stageAxes(i)), ...
                        MainGuiStyle.num(scanParams.to(k)), char(8594), MainGuiStyle.num(scanParams.from(k))));
                else
                    set(obj.btnDir(i), 'String', char(8594), 'TooltipString', ...
                        sprintf('%s is scanned from From to To (%s %s %s %s). Click: the other way.', upper(obj.stageAxes(i)), ...
                        MainGuiStyle.num(scanParams.from(k)), char(8594), MainGuiStyle.num(scanParams.to(k)), um));
                end

                MainGuiStyle.setFieldEnabled([obj.edtFrom(i), obj.edtTo(i), obj.edtNumPoints(i), ...
                    obj.edtStep(i), obj.edtSpan(i)], isScanned);
                obj.btnDir(i).Enable = BooleanHelper.boolToOnOff(isScanned);
                MainGuiStyle.setFieldEnabled(obj.edtCentre(i), true);
            end

            % Summary line, with the shortest it can take (pixel time only)
            summary = scanParams.describe(stage.ReturnMaxScanSize(2));
            if scanParams.nDims > 0
                nPixels = prod(scanParams.numPoints(~scanParams.isFixed)) * (1 + logical(scanParams.isMWContrastScan));
                summary = sprintf('%s · at least %s', summary, MainGuiStyle.duration(nPixels * scanParams.pixelTime));
            end
            obj.tvSummary.String = summary;
            obj.tvSummary.TooltipString = 'At least = pixel time × pixels; moving between pixels adds to it.';
        end
    end

    %% Callbacks
    methods
        function cbxScanCallback(obj, index)
            [stage, scanParams, k] = obj.getStageAndParams(index);
            scanParams.isFixed(k) = ~obj.cbxScan(index).Value;
            stage.sendEventScanParamsChanged();
        end

        function btnDirCallback(obj, index)
            [stage, scanParams, k] = obj.getStageAndParams(index);
            if numel(scanParams.reverse) < 3; scanParams.reverse(end+1:3) = false; end
            scanParams.reverse(k) = ~scanParams.reverse(k);
            stage.sendEventScanParamsChanged();
        end

        function edtNumPointsChangedCallback(obj, index)
            [stage, scanParams, k] = obj.getStageAndParams(index);
            str = obj.edtNumPoints(index).String;
            if ValidationHelper.isValuePositiveInteger(str) && str2double(str) >= 2
                scanParams.numPoints(k) = str2double(str);
            else
                obj.sendWarning('Points must be a whole number of at least 2. Reverting.');
            end
            stage.sendEventScanParamsChanged();
        end

        function edtStepCallback(obj, index)
            % A step: Points = the number that fits From-To with (about) that step
            [stage, scanParams, k] = obj.getStageAndParams(index);
            str = obj.edtStep(index).String;
            span = abs(scanParams.to(k) - scanParams.from(k));
            if ~ValidationHelper.isValuePositive(str) || span == 0
                obj.sendWarning('The step must be a positive number (and From must differ from To). Reverting.');
            else
                scanParams.numPoints(k) = max(2, round(span / str2double(str)) + 1);
            end
            stage.sendEventScanParamsChanged();
        end

        function edtSpanCallback(obj, index)
            % A span: From/To = Centre ± span/2, clipped to the soft limits
            [stage, scanParams, k] = obj.getStageAndParams(index);
            str = obj.edtSpan(index).String;
            if ~ValidationHelper.isValueNonNegative(str)
                obj.sendWarning('The span must be a non-negative number. Reverting.');
            else
                obj.setRangeAround(stage, scanParams, k, str2double(str));
            end
            stage.sendEventScanParamsChanged();
        end

        function edtFromChangedCallback(obj, index)
            obj.setEnd(index, 'from');
        end

        function edtToChangedCallback(obj, index)
            obj.setEnd(index, 'to');
        end

        function setEnd(obj, index, which)
            [stage, scanParams, k] = obj.getStageAndParams(index);
            if strcmp(which, 'from'); field = obj.edtFrom(index); else; field = obj.edtTo(index); end
            str = field.String;
            [lowerBound, upperBound] = stage.ReturnLimits(k);
            if ~ValidationHelper.isStringValueANumber(str)
                obj.sendWarning('Only numbers can be accepted! Reverting.');
            elseif ~ValidationHelper.isInBorders(str2double(str), lowerBound, upperBound)
                obj.sendWarning(sprintf('"%s" of %s must be within the soft limits [%s, %s] %s. Reverting.', ...
                    which, upper(obj.stageAxes(index)), MainGuiStyle.num(lowerBound), MainGuiStyle.num(upperBound), MainGuiStyle.micron));
            else
                scanParams.(which)(k) = str2double(str);
                if scanParams.from(k) == scanParams.to(k)
                    obj.sendWarning(sprintf('The scan over %s starts and ends at the same point!', upper(obj.stageAxes(index))));
                end
            end
            stage.sendEventScanParamsChanged();
        end

        function edtFixedChangedCallback(obj, index)
            % Centre (scanParams.fixedPos)
            [stage, scanParams, k] = obj.getStageAndParams(index);
            str = obj.edtCentre(index).String;
            if ~ValidationHelper.isStringValueANumber(str)
                obj.sendWarning('Only numbers can be accepted! Reverting.');
                stage.sendEventScanParamsChanged();
                return
            end
            centre = str2double(str);
            [lowerBound, upperBound] = stage.ReturnLimits(k);
            if ~ValidationHelper.isInBorders(centre, lowerBound, upperBound)
                obj.sendWarning(sprintf('The centre of %s must be within the soft limits [%s, %s] %s. Set to the limit.', ...
                    upper(obj.stageAxes(index)), MainGuiStyle.num(lowerBound), MainGuiStyle.num(upperBound), MainGuiStyle.micron));
                centre = min(max(centre, lowerBound), upperBound);
            end
            scanParams.fixedPos(k) = centre;
            stage.sendEventScanParamsChanged();
        end

        function btnCentreHereCallback(obj)
            stage = obj.getStageAndParams();
            stage.sendPosToScanParams();    % centre = current position, all axes (sends an event)
        end

        function btnRangeAroundCallback(obj)
            [stage, scanParams] = obj.getStageAndParams();
            for i = 1:length(obj.stageAxes)
                k = ClassStage.getAxis(obj.stageAxes(i));
                if ~scanParams.isFixed(k)
                    obj.setRangeAround(stage, scanParams, k, abs(scanParams.to(k) - scanParams.from(k)));
                end
            end
            stage.sendEventScanParamsChanged();
        end

        function btnFullRangeCallback(obj)
            [stage, scanParams] = obj.getStageAndParams();
            for i = 1:length(obj.stageAxes)
                k = ClassStage.getAxis(obj.stageAxes(i));
                if ~scanParams.isFixed(k)
                    [scanParams.from(k), scanParams.to(k)] = stage.ReturnLimits(k);
                end
            end
            stage.sendEventScanParamsChanged();
        end

        function btnGoCentreCallback(obj)
            [stage, scanParams] = obj.getStageAndParams();
            k = ClassStage.getAxis(obj.stageAxes);
            stage.move(obj.stageAxes, scanParams.fixedPos(k));
        end
    end

    methods (Access = private)
        function edits = column(obj, grid, title, n, tooltip) %#ok<INUSL>
            MainGuiStyle.header(grid, title);
            edits = gobjects(1, n);
            for i = 1:n
                edits(i) = MainGuiStyle.edit(grid, '', 'TooltipString', tooltip);
            end
        end

        function [stage, scanParams, k] = getStageAndParams(obj, index)
            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
            scanParams = stage.scanParams;
            k = [];
            if nargin > 1
                k = ClassStage.getAxis(obj.stageAxes(index));
            end
        end

        function setRangeAround(obj, stage, scanParams, k, span) %#ok<INUSL>
            [lowerBound, upperBound] = stage.ReturnLimits(k);
            c = scanParams.fixedPos(k);
            scanParams.from(k) = max(c - span/2, lowerBound);
            scanParams.to(k) = min(c + span/2, upperBound);
        end
    end

    methods (Static)
        function s = stepOf(scanParams, k)
            if scanParams.numPoints(k) < 2
                s = NaN;
            else
                s = abs(scanParams.to(k) - scanParams.from(k)) / (scanParams.numPoints(k) - 1);
            end
        end
    end

    %% overridden from EventListener
    methods
        function onEvent(obj, event)
            if event.isError || isfield(event.extraInfo, ClassStage.EVENT_SCAN_PARAMS_CHANGED) ...
                    || isfield(event.extraInfo, ClassStage.EVENT_LIM_CHANGED)
                obj.refresh()
            end
        end
    end
end
