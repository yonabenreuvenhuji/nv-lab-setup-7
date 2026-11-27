classdef ViewStagePanelScan < GuiComponent & EventListener
    %VIEWSTAGEPANELSCAN The "Scan" section: start/stop, progress, and how the scan runs.
    %
    %   [  Scan  ]  [██████████░░░░░░░░░░░░]  42 %
    %               Line 17 of 41 · 0:48 min left
    %   Pixel time (ms) [ 15 ]   [ ] Repeat   [x] Auto-save   [ ] Fast scan
    %   Line order [ Lines along X, step Y ▾]    Pattern [ Same direction ▾]
    %
    % Line order: which scanned axis runs along each line, which steps between lines (and,
    % in 3D, which steps between slices). Pattern: every line in the same direction (raster:
    % immune to backlash) or every other line backwards (zigzag: no return trip). The
    % pattern needs a stage that can do both; it is greyed out otherwise.
    %
    % Pixel time is shown in ms; the scan parameters keep it in seconds.

    properties
        btnScan             % ButtonStartStop
        progressBox         % uix.HBox holding the two halves of the bar
        progressFill        % filled part of the bar
        progressRest        % empty part
        tvPercent           % '42 %'
        tvProgress          % 'Line 17 of 41 · 0:48 min left'
        edtPixelTime        % edit, ms
        cbxContinuous       % checkbox
        cbxFastScan         % checkbox
        cbxAutoSave         % checkbox
        cbxReturn           % checkbox: back to the start position after the scan
        popOrder            % line order
        popPattern          % same direction / zigzag

        stageName           % string
    end

    properties (Access = private)
        orderCodes = {''}   % scanOrder string for each entry of popOrder
    end

    properties (Constant)
        PATTERN_OPTIONS = {'Same direction (raster)', 'Zigzag (serpentine)'}
    end

    methods
        function obj = ViewStagePanelScan(parent, controller, stageName)
            obj@GuiComponent(parent, controller);
            expNames = Experiment.getExperimentNames();
            obj@EventListener([{stageName}, expNames(:)', {StageScanner.NAME}]);
            obj.stageName = stageName;
            S = MainGuiStyle;
            C = MainGuiStyle.col();

            [panel, inner] = MainGuiStyle.section(parent.component, 'Scan');
            obj.component = panel;

            %%%% Row 1: Scan button + progress %%%%
            top = MainGuiStyle.hbox(inner, 10);
            obj.btnScan = ButtonStartStop(top, 'Scan', 'Stop scan');
            obj.btnScan.colors = {C.accent, C.stop};
            obj.btnScan.isRunning = false;
            set(obj.btnScan.component, 'FontSize', 13, 'TooltipString', 'Start a scan of the Scan area');
            progress = MainGuiStyle.vbox(top, 3);
            barRow = MainGuiStyle.hbox(progress, 8);
            obj.progressBox = uix.HBox('Parent', barRow, 'Spacing', 0, 'Padding', 1, 'BackgroundColor', C.line);
            obj.progressFill = uicontrol('Parent', obj.progressBox, 'Style', 'text', 'String', '', ...
                'BackgroundColor', C.accent, 'Enable', 'inactive');
            obj.progressRest = uicontrol('Parent', obj.progressBox, 'Style', 'text', 'String', '', ...
                'BackgroundColor', [1 1 1], 'Enable', 'inactive');
            obj.progressBox.Widths = [0 -1];
            obj.tvPercent = uicontrol('Parent', barRow, 'Style', 'text', 'String', '', ...
                'FontSize', S.FONT + 1, 'FontWeight', 'bold', 'ForegroundColor', C.ink, ...
                'BackgroundColor', C.panel, 'HorizontalAlignment', 'right');
            barRow.Widths = [-1 44];
            obj.tvProgress = MainGuiStyle.hint(progress, 'Ready');
            progress.Heights = [16 -1];
            top.Widths = [130 -1];

            %%%% Row 2: settings %%%%
            row2 = MainGuiStyle.hbox(inner, 10);
            ptBox = MainGuiStyle.hbox(row2, 4);
            MainGuiStyle.label(ptBox, 'Pixel time (ms)', 'right');
            obj.edtPixelTime = MainGuiStyle.edit(ptBox, '', 'TooltipString', ...
                'How long the detector measures at each pixel, in milliseconds');
            ptBox.Widths = [90 60];
            obj.cbxContinuous = MainGuiStyle.checkbox(row2, 'Repeat', false, 'TooltipString', ...
                'Scan again and again until Stop');
            obj.cbxAutoSave = MainGuiStyle.checkbox(row2, 'Auto-save', false, 'TooltipString', ...
                'Save every finished scan to the autosave folder');

            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
            fastScannable = stage.hasFastScan;
            slowScannable = stage.hasSlowScan;
            obj.cbxFastScan = MainGuiStyle.checkbox(row2, 'Fast scan', fastScannable, ...
                'Enable', BooleanHelper.boolToOnOff(fastScannable && slowScannable));
            if fastScannable && slowScannable
                obj.cbxFastScan.TooltipString = 'Continuous (on-the-fly) scanning instead of stopping at every pixel';
            elseif fastScannable
                obj.cbxFastScan.TooltipString = 'This stage only scans on the fly';
            else
                obj.cbxFastScan.TooltipString = 'This stage stops at every pixel; it has no fast (on-the-fly) scan';
            end
            obj.cbxReturn = MainGuiStyle.checkbox(row2, 'Return to start', true, 'TooltipString', ...
                ['When the scan ends or is stopped, move the stage back to where it was ' ...
                'when Scan was pressed. (Not after Halt.)']);
            row2.Widths = [150 70 80 80 -1];

            %%%% Row 3: order and pattern %%%%
            row3 = MainGuiStyle.hbox(inner, 6);
            MainGuiStyle.label(row3, 'Line order', 'left');
            obj.popOrder = MainGuiStyle.popup(row3, {'Automatic'}, 'TooltipString', ...
                'Which axis runs along each line, which steps between lines (and between slices, in 3D)');
            MainGuiStyle.label(row3, 'Pattern', 'right');
            obj.popPattern = MainGuiStyle.popup(row3, obj.PATTERN_OPTIONS);
            if isprop(stage, 'scanBidirectional')
                obj.popPattern.TooltipString = sprintf(['Same direction: every line runs From%sTo (no backlash between lines).\n' ...
                    'Zigzag: every other line runs backwards (no return trip, faster).'], char(8594));
            else
                obj.popPattern.Enable = 'off';
                obj.popPattern.TooltipString = 'This stage decides the line direction itself';
            end
            row3.Widths = [64 -1.4 50 -1];

            inner.Heights = [40 S.ROW_H S.ROW_H];
            obj.height = MainGuiStyle.sectionHeight(MainGuiStyle.rowsHeight(inner.Heights));
            obj.width = 560;

            %%%% callbacks %%%%
            obj.cbxAutoSave.Callback = @(h,e) obj.setParam('autoSave', logical(obj.cbxAutoSave.Value));
            obj.cbxReturn.Callback = @(h,e) obj.setParam('returnToStart', logical(obj.cbxReturn.Value));
            obj.cbxContinuous.Callback = @(h,e) obj.setParam('continuous', logical(obj.cbxContinuous.Value));
            obj.cbxFastScan.Callback = @(h,e) obj.setParam('fastScan', logical(obj.cbxFastScan.Value));
            obj.edtPixelTime.Callback = @(h,e) obj.edtPixelTimeCallback;
            obj.popOrder.Callback = @(h,e) obj.setParam('scanOrder', obj.orderCodes{obj.popOrder.Value});
            obj.popPattern.Callback = @(h,e) obj.setParam('zigzag', obj.popPattern.Value == 2);
            obj.btnScan.startCallback = @(h,e) obj.btnScanCallback;
            obj.btnScan.stopCallback = @(h,e) obj.btnStopScanCallback;

            obj.refresh();  % init values
        end

        function refresh(obj)
            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
            scanParams = stage.scanParams;

            obj.cbxAutoSave.Value = scanParams.autoSave;
            obj.cbxContinuous.Value = scanParams.continuous;
            obj.cbxFastScan.Value = scanParams.fastScan;
            obj.cbxReturn.Value = scanParams.returnToStart;
            obj.edtPixelTime.String = MainGuiStyle.num(scanParams.pixelTime * 1e3, 3);
            obj.popPattern.Value = 1 + logical(scanParams.zigzag);
            obj.refreshOrder(scanParams);
        end

        function refreshOrder(obj, scanParams)
            % The line-order choices for the axes being scanned now
            L = upper(ClassStage.SCAN_AXES);
            scanned = find(~scanParams.isFixed);
            switch numel(scanned)
                case 2
                    a = scanned(1); b = scanned(2);
                    items = {'Automatic', ...
                        sprintf('Lines along %s, step %s', L(a), L(b)), ...
                        sprintf('Lines along %s, step %s', L(b), L(a))};
                    codes = {'', lower([L(a) L(b)]), lower([L(b) L(a)])};
                case 3
                    p = perms(1:3);
                    p = sortrows(p);    % xyz, xzy, yxz, ...
                    items = {sprintf('Automatic (lines %s, step %s, slices %s)', L(1), L(2), L(3))};
                    codes = {''};
                    for r = 1:size(p, 1)
                        items{end + 1} = sprintf('Lines %s, step %s, slices %s', L(p(r, 1)), L(p(r, 2)), L(p(r, 3))); %#ok<AGROW>
                        codes{end + 1} = lower(L(p(r, :))); %#ok<AGROW>
                    end
                otherwise
                    items = {BooleanHelper.ifTrueElse(numel(scanned) == 1, ...
                        sprintf('Along %s (1D)', L(scanned)), '—')};
                    codes = {''};
            end
            obj.orderCodes = codes;
            % The choice that matches the parameters ('' = automatic)
            current = 1;
            if ~scanParams.isOrderAutomatic
                order = lower(L(scanParams.orderedScanAxes));
                current = find(strcmp(codes, order), 1);
                if isempty(current); current = 1; end
            end
            set(obj.popOrder, 'String', items, 'Value', current, ...
                'Enable', BooleanHelper.boolToOnOff(numel(scanned) >= 2));
            obj.popPattern.Enable = BooleanHelper.boolToOnOff(numel(scanned) >= 2 && ...
                isprop(getObjByName(obj.stageName), 'scanBidirectional'));
        end

        function setParam(obj, name, value)
            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
            stage.scanParams.(name) = value;
            stage.sendEventScanParamsChanged();
        end

        function edtPixelTimeCallback(obj)
            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
            if ~ValidationHelper.isValuePositive(obj.edtPixelTime.String)
                EventStation.anonymousWarning('Pixel time has to be a positive number (in ms)! Reverting.');
            else
                stage.scanParams.pixelTime = str2double(obj.edtPixelTime.String) / 1e3;
            end
            stage.sendEventScanParamsChanged();
        end

        function btnScanCallback(obj)
            scanner = getObjByName(StageScanner.NAME);
            scanner.switchTo(obj.stageName);
            scanner.startScan();
            obj.btnScan.isRunning = false;  % Scan is now finished
        end

        function btnStopScanCallback(obj) %#ok<MANU>
            scanner = getObjByName(StageScanner.NAME);
            scanner.stopScan();
        end

        function showProgress(obj, p)
            % p: StageScanner.progressInfo
            C = MainGuiStyle.col();
            f = min(max(p.fraction, 0), 1);
            if f <= 0
                obj.progressBox.Widths = [0 -1];
            elseif f >= 1
                obj.progressBox.Widths = [-1 0];
            else
                obj.progressBox.Widths = [-f -(1 - f)];
            end
            obj.tvPercent.String = sprintf('%d %%', floor(100 * f));
            parts = {};
            if p.pass > 1; parts{end + 1} = sprintf('Pass %d', p.pass); end
            if p.nSlices > 0; parts{end + 1} = sprintf('Slice %d of %d', p.slice, p.nSlices); end
            if p.dims >= 2 && p.nLines > 0; parts{end + 1} = sprintf('Line %d of %d', min(p.line, p.nLines), p.nLines); end
            if p.running
                obj.progressFill.BackgroundColor = C.accent;
                if isfinite(p.remaining)
                    parts{end + 1} = sprintf('%s left', MainGuiStyle.duration(p.remaining));
                end
            else
                if f >= 1
                    obj.progressFill.BackgroundColor = C.ok;
                    parts = {sprintf('Done in %s', MainGuiStyle.duration(p.elapsed))};
                else
                    obj.progressFill.BackgroundColor = C.ink3;
                    parts{end + 1} = sprintf('stopped after %s', MainGuiStyle.duration(p.elapsed));
                end
            end
            obj.tvProgress.String = strjoin(parts, ' · ');
            drawnow limitrate nocallbacks   % (progress only: button clicks wait for the end of the line)
        end
    end

    %% overridden from EventListener
    methods
        function onEvent(obj, event)
            if isfield(event.extraInfo, StageScanner.EVENT_SCAN_PROGRESS)
                obj.showProgress(event.extraInfo.(StageScanner.EVENT_SCAN_PROGRESS));
                return
            end

            if event.isError ...
                    || isfield(event.extraInfo, StageScanner.EVENT_SCAN_FINISHED)
                % Scan stopped, either intentionally or not
                obj.refresh();
                obj.btnScan.Enable = 'on';
                obj.btnScan.isRunning = false;
                scanner = getObjByName(StageScanner.NAME);
                if ~isempty(scanner) && isfield(event.extraInfo, StageScanner.EVENT_SCAN_FINISHED)
                    obj.showProgress(scanner.progressInfo);
                end
            end

            if isfield(event.extraInfo, StageScanner.EVENT_SCAN_STARTED)
                obj.btnScan.isRunning = true;
            end

            if isfield(event.extraInfo, ClassStage.EVENT_SCAN_PARAMS_CHANGED)
                obj.refresh();
            end

            if isfield(event.extraInfo, Experiment.EVENT_EXP_RESUMED) ...
                    || isfield(event.extraInfo, Experiment.EVENT_EXP_PAUSED)
                % When experiments run, we can't scan
                exp = event.creator;
                isScanPossible = ~exp.isRunning;
                obj.btnScan.Enable = BooleanHelper.boolToOnOff(isScanPossible);
            end
        end
    end
end
