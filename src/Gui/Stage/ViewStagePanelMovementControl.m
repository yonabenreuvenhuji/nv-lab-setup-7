classdef ViewStagePanelMovementControl < GuiComponent & EventListener
    %VIEWSTAGEPANELMOVEMENTCONTROL The "Position" section: where the stage is, and moving it.
    %
    %   Axis   −   Position (µm)   +   Target (µm)  |  Step (µm) [ 1 ]  0.1  1  10  100
    %    X    [−]  [  12.345  ]   [+]  [        ]   |  [ Go to target ]  [ Clear ]
    %    Y    ...                                   |  [ ↶ Back ]        [ HALT  ]
    %   Saved positions [ NV1 (12.3, -8.1, 0) ▾]  [Go]  [Save current…]  [Delete]
    %
    % TARGETS. Typing a position in an axis' Target field (and Enter) makes it pending: the
    % field turns amber, and "Go to target" lights up. Go moves every pending axis at once,
    % then clears the targets. Clear (or emptying the field) drops them. This replaces the
    % old "click the position, type into the blue box, press Move To Blue", where the
    % target was a pale blue field that only appeared after a click and was easy to miss.
    %
    % BACK returns to where the stage was before the last move made from this section
    % (arrows, Go, saved positions). SAVED POSITIONS are kept between MATLAB sessions
    % (MATLAB preferences, per stage).

    properties
        btnMinus            % 1xN button: one step down
        btnPlus             % 1xN button: one step up
        tvPos               % 1xN read-only position
        edtTarget           % 1xN target field
        edtStep             % edit
        btnStepPreset       % buttons
        btnGo               % move to the targets
        btnClear            % clear the targets
        btnBack             % back to the previous position
        btnHaltStage        % halt
        popSaved            % saved positions
        btnSavedGo
        btnSavedAdd
        btnSavedDelete
        cbxJoystick         % checkbox (stages with a joystick)

        stageName           % string
        stageAxes           % string, e.g. 'xyz'
    end

    properties (Access = private)
        targets             % 1xN double, NaN = no target
        history = zeros(0, 3)   % previous positions (rows), most recent last
    end

    properties (Constant)
        STEP_PRESETS = [0.1 1 10 100]   % um
        HISTORY_LENGTH = 30
        PREF_GROUP = 'nvLabMainGui'
    end

    methods
        function obj = ViewStagePanelMovementControl(parent, controller, stage)
            expNames = Experiment.getExperimentNames();
            namesToListenTo = [{stage.name}, expNames(:)', {StageScanner.NAME}];
            obj@GuiComponent(parent, controller);
            obj@EventListener(namesToListenTo);
            obj.stageName = stage.name;
            obj.stageAxes = stage.availableAxes;
            n = length(obj.stageAxes);
            obj.targets = nan(1, n);
            S = MainGuiStyle;
            um = MainGuiStyle.micron;

            [panel, inner] = MainGuiStyle.section(parent.component, 'Position');
            obj.component = panel;

            %%%% Top: axes (left) | step and actions (right) %%%%
            top = MainGuiStyle.hbox(inner, 14);
            grid = MainGuiStyle.grid(top, 4);
            % column 1: axis
            MainGuiStyle.empty(grid);
            for i = 1:n; MainGuiStyle.axisTag(grid, obj.stageAxes(i)); end
            % column 2: minus
            MainGuiStyle.empty(grid);
            obj.btnMinus = gobjects(1, n);
            for i = 1:n
                obj.btnMinus(i) = MainGuiStyle.button(grid, char(8722), 'normal', 'FontSize', 12, ...
                    'TooltipString', sprintf('Move %s down by one step', upper(obj.stageAxes(i))));
            end
            % column 3: position
            MainGuiStyle.header(grid, sprintf('Position (%s)', um));
            obj.tvPos = gobjects(1, n);
            for i = 1:n
                obj.tvPos(i) = MainGuiStyle.readout(grid, '', ...
                    'TooltipString', sprintf('Where the stage is now (%s axis)', upper(obj.stageAxes(i))));
            end
            % column 4: plus
            MainGuiStyle.empty(grid);
            obj.btnPlus = gobjects(1, n);
            for i = 1:n
                obj.btnPlus(i) = MainGuiStyle.button(grid, '+', 'normal', 'FontSize', 12, ...
                    'TooltipString', sprintf('Move %s up by one step', upper(obj.stageAxes(i))));
            end
            % column 5: target
            MainGuiStyle.header(grid, sprintf('Target (%s)', um));
            obj.edtTarget = gobjects(1, n);
            for i = 1:n
                obj.edtTarget(i) = MainGuiStyle.edit(grid, '', 'TooltipString', ...
                    sprintf(['Type a position for %s and press Enter. It turns amber until you press ' ...
                    '"Go to target"; empty = %s does not move.'], upper(obj.stageAxes(i)), upper(obj.stageAxes(i))));
            end
            set(grid, 'Widths', [20 30 92 30 92], 'Heights', [S.HEAD_H, S.ROW_H * ones(1, n)]);
            gridHeight = S.HEAD_H + n*S.ROW_H + 4*n;

            right = MainGuiStyle.vbox(top, S.GAP);
            stepRow = MainGuiStyle.hbox(right, 4);
            MainGuiStyle.label(stepRow, sprintf('Step (%s)', um), 'right');
            obj.edtStep = MainGuiStyle.edit(stepRow, '', 'TooltipString', ...
                sprintf('How far − and + move (minimum %s %s)', MainGuiStyle.num(stage.STEP_MINIMUM_SIZE), um));
            obj.btnStepPreset = gobjects(1, numel(obj.STEP_PRESETS));
            for k = 1:numel(obj.STEP_PRESETS)
                obj.btnStepPreset(k) = MainGuiStyle.button(stepRow, MainGuiStyle.num(obj.STEP_PRESETS(k)), 'quiet', ...
                    'TooltipString', sprintf('Step = %s %s', MainGuiStyle.num(obj.STEP_PRESETS(k)), um));
            end
            stepRow.Widths = [56 56 -1 -1 -1 -1];
            goRow = MainGuiStyle.hbox(right, S.GAP);
            obj.btnGo = MainGuiStyle.button(goRow, 'Go to target', 'normal', 'Enable', 'off', ...
                'TooltipString', 'Move every axis that has a target (amber) to it, all at once');
            obj.btnClear = MainGuiStyle.button(goRow, 'Clear', 'quiet', 'Enable', 'off', ...
                'TooltipString', 'Forget the targets');
            goRow.Widths = [-2 -1];
            backRow = MainGuiStyle.hbox(right, S.GAP);
            obj.btnBack = MainGuiStyle.button(backRow, [char(8630) ' Back'], 'normal', 'Enable', 'off', ...
                'TooltipString', 'Return to where the stage was before the last move from this section');
            obj.btnHaltStage = MainGuiStyle.button(backRow, 'Halt', 'stop', ...
                'TooltipString', 'Stop every stage movement now (also stops a scan)');
            backRow.Widths = [-1 -1];
            MainGuiStyle.empty(right);
            right.Heights = [S.ROW_H S.ROW_H S.ROW_H -1];
            top.Widths = [sum([20 30 92 30 92]) + 16, -1];

            %%%% Saved positions %%%%
            savedRow = MainGuiStyle.hbox(inner, S.GAP);
            MainGuiStyle.label(savedRow, 'Saved positions', 'left');
            obj.popSaved = MainGuiStyle.popup(savedRow, {'(none saved)'}, ...
                'TooltipString', 'Positions you saved, kept between sessions');
            obj.btnSavedGo = MainGuiStyle.button(savedRow, 'Go', 'normal', ...
                'TooltipString', 'Move the stage to the selected saved position');
            obj.btnSavedAdd = MainGuiStyle.button(savedRow, 'Save current...', 'normal', ...
                'TooltipString', 'Save where the stage is now, under a name');
            obj.btnSavedDelete = MainGuiStyle.button(savedRow, 'Delete', 'quiet', ...
                'TooltipString', 'Delete the selected saved position');
            widths = [100 -1 40 100 56];
            if stage.hasJoystick
                obj.cbxJoystick = MainGuiStyle.checkbox(savedRow, 'Joystick', false, ...
                    'Callback', @(h,e) obj.cbxJoystickCallback);
                widths(end + 1) = 72;
            end
            savedRow.Widths = widths;

            inner.Heights = [max(gridHeight, 3*S.ROW_H + 2*S.GAP), S.ROW_H];
            obj.height = MainGuiStyle.sectionHeight(MainGuiStyle.rowsHeight(inner.Heights));
            obj.width = 600;

            %%%% Callbacks %%%%
            obj.edtStep.Callback = @(h,e) obj.edtStepSizeCallback();
            for k = 1:numel(obj.STEP_PRESETS)
                obj.btnStepPreset(k).Callback = @(h,e) obj.setStep(obj.STEP_PRESETS(k));
            end
            for i = 1:n
                obj.btnMinus(i).Callback = @(h,e) obj.btnMoveCallback(i, -1);
                obj.btnPlus(i).Callback = @(h,e) obj.btnMoveCallback(i, +1);
                obj.edtTarget(i).Callback = @(h,e) obj.edtTargetCallback(i);
            end
            obj.btnGo.Callback = @(h,e) obj.btnGoCallback();
            obj.btnClear.Callback = @(h,e) obj.clearTargets();
            obj.btnBack.Callback = @(h,e) obj.btnBackCallback();
            obj.btnHaltStage.Callback = @(h,e) StageControlEvents.sendHalt;
            obj.btnSavedGo.Callback = @(h,e) obj.btnSavedGoCallback();
            obj.btnSavedAdd.Callback = @(h,e) obj.btnSavedAddCallback();
            obj.btnSavedDelete.Callback = @(h,e) obj.btnSavedDeleteCallback();

            obj.refreshSaved();
            obj.refresh();
        end

        function refresh(obj)
            stage = obj.getStage();
            obj.edtStep.String = MainGuiStyle.num(stage.stepSize, 4);
            currentPosition = stage.Pos(obj.stageAxes);
            for i = 1: length(obj.stageAxes)
                obj.tvPos(i).String = sprintf('%.3f', currentPosition(i));
            end
            obj.checkMovementEnabled;
            drawnow % If we are in mid-operation
        end
    end

    %% Moving
    methods
        function edtStepSizeCallback(obj)
            stage = obj.getStage();
            stepSizeString = obj.edtStep.String;
            if ~ValidationHelper.isStringValueANumber(stepSizeString) ...
                    || ~ValidationHelper.isStringValueInBorders(stepSizeString, stage.STEP_MINIMUM_SIZE, inf)
                obj.edtStep.String = MainGuiStyle.num(stage.stepSize, 4);
                EventStation.anonymousWarning('The step must be a number of at least %s %s. Reverting.', ...
                    MainGuiStyle.num(stage.STEP_MINIMUM_SIZE), MainGuiStyle.micron);
                return
            end
            stage.stepSize = str2double(stepSizeString);
        end

        function setStep(obj, step)
            stage = obj.getStage();
            stage.stepSize = max(step, stage.STEP_MINIMUM_SIZE);
        end

        function btnMoveCallback(obj, index, direction)
            stage = obj.getStage();
            phAxis = ClassStage.getAxis(obj.stageAxes(index));
            obj.remember();
            stage.relativeMove(phAxis, direction * stage.stepSize);
        end

        function edtTargetCallback(obj, index)
            % Enter in a target field: check it, mark it pending (or clear it)
            stage = obj.getStage();
            str = strtrim(obj.edtTarget(index).String);
            if isempty(str)
                obj.targets(index) = NaN;
            else
                phAxis = ClassStage.getAxis(obj.stageAxes(index));
                [limNeg, limPos] = stage.ReturnLimits(phAxis);
                if ValidationHelper.isStringValueInBorders(str, limNeg, limPos)
                    obj.targets(index) = str2double(str);
                    obj.edtTarget(index).String = MainGuiStyle.num(obj.targets(index), 4);
                else
                    obj.targets(index) = NaN;
                    obj.edtTarget(index).String = '';
                    EventStation.anonymousWarning('%s target must be a number within the soft limits [%s, %s] %s.', ...
                        upper(obj.stageAxes(index)), MainGuiStyle.num(limNeg), MainGuiStyle.num(limPos), MainGuiStyle.micron);
                end
            end
            obj.showTargets();
        end

        function btnGoCallback(obj)
            idx = find(~isnan(obj.targets));
            if isempty(idx); return; end
            stage = obj.getStage();
            axesToMove = obj.stageAxes(idx);
            positions = obj.targets(idx);
            obj.remember();
            obj.clearTargets();
            stage.move(axesToMove, positions);
        end

        function clearTargets(obj)
            obj.targets(:) = NaN;
            set(obj.edtTarget, 'String', '');
            obj.showTargets();
        end

        function showTargets(obj)
            % Amber for pending targets; "Go to target" lit when there is one
            for i = 1:numel(obj.targets)
                MainGuiStyle.setPending(obj.edtTarget(i), ~isnan(obj.targets(i)));
            end
            obj.checkMovementEnabled();
        end

        function btnBackCallback(obj)
            if isempty(obj.history); return; end
            stage = obj.getStage();
            previous = obj.history(end, 1:numel(obj.stageAxes));
            obj.history(end, :) = [];
            stage.move(obj.stageAxes, previous);
            obj.checkMovementEnabled();
        end

        function remember(obj)
            % Keep the current position, for Back
            stage = obj.getStage();
            pos = zeros(1, 3);
            pos(1:numel(obj.stageAxes)) = stage.Pos(obj.stageAxes);
            if isempty(obj.history) || any(abs(obj.history(end, :) - pos) > 1e-6)
                obj.history(end + 1, :) = pos;
                if size(obj.history, 1) > obj.HISTORY_LENGTH
                    obj.history(1, :) = [];
                end
            end
        end

        function cbxJoystickCallback(obj)
            % Enables or disable the joystick according to the GUI settings.
            % While the joystick is on, updates the location every ~300ms.
            stage = obj.getStage();
            jstick = getObjByName(Joystick.NAME);
                if isempty(jstick); throwBaseObjException(Joystick.NAME); end
            stageScanner = getObjByName(StageScanner.NAME);

            isEnabled = obj.cbxJoystick.Value;
            % Don't let other things interrupt
            if strcmp(stageScanner.mStageName, obj.stageName) ...
                    && stageScanner.mCurrentlyScanning
                stageScanner.stopScan;
            end
            jstick.isEnabled = isEnabled;
            obj.checkMovementEnabled;
            stage.JoystickControl(isEnabled);
            stage.sendEventPositionChanged;
            while jstick.isEnabled
                jstick.buttonAction;
                stage.sendEventPositionChanged;
                drawnow();
            end
        end

        function checkMovementEnabled(obj)
            % Movement is possible when: closed loop, no scan running, joystick off, and
            % no experiment running.
            stage = getObjByName(obj.stageName);
            isEnabled = ~isempty(stage) && strcmp(stage.loopMode, 'Closed');
            scanner = getObjByName(StageScanner.NAME);
            isEnabled = isEnabled && ~scanner.mCurrentlyScanning;
            if ~isempty(stage) && stage.hasJoystick
                jstick = getObjByName(Joystick.NAME);
                if ~isempty(jstick)
                    isEnabled = isEnabled && ~jstick.isEnabled;
                end
            end
            exp = getObjByName(Experiment.current);
            if ~isempty(exp) && ~isa(exp, 'SpcmCounter')
                isEnabled = isEnabled && ~exp.isRunning;
            end

            onOff = BooleanHelper.boolToOnOff(isEnabled);
            set([obj.btnMinus, obj.btnPlus, obj.btnSavedGo], 'Enable', onOff);
            hasTarget = any(~isnan(obj.targets));
            obj.btnGo.Enable = BooleanHelper.boolToOnOff(isEnabled && hasTarget);
            obj.btnClear.Enable = BooleanHelper.boolToOnOff(hasTarget);
            MainGuiStyle.paintButton(obj.btnGo, BooleanHelper.ifTrueElse(isEnabled && hasTarget, 'pending', 'normal'));
            obj.btnBack.Enable = BooleanHelper.boolToOnOff(isEnabled && ~isempty(obj.history));
            if ~isempty(obj.history)
                obj.btnBack.TooltipString = sprintf('Return to (%s) %s, where the stage was before the last move from this section', ...
                    strjoin(arrayfun(@(v) MainGuiStyle.num(v), obj.history(end, 1:numel(obj.stageAxes)), 'UniformOutput', false), ', '), ...
                    MainGuiStyle.micron);
            end
        end
    end

    %% Saved positions
    methods
        function list = savedList(obj)
            list = struct('name', {}, 'pos', {});
            try
                if ispref(obj.PREF_GROUP, obj.prefName)
                    list = getpref(obj.PREF_GROUP, obj.prefName);
                end
            catch
            end
        end

        function saveList(obj, list)
            setpref(obj.PREF_GROUP, obj.prefName, list);
            obj.refreshSaved();
        end

        function name = prefName(obj)
            name = ['savedPositions_' matlab.lang.makeValidName(obj.stageName)];
        end

        function refreshSaved(obj)
            list = obj.savedList;
            if isempty(list)
                set(obj.popSaved, 'String', {'(none saved)'}, 'Value', 1);
                set([obj.btnSavedDelete], 'Enable', 'off');
                obj.btnSavedGo.Visible = 'on';
                obj.btnSavedGo.Enable = 'off';
                return
            end
            items = cell(1, numel(list));
            for k = 1:numel(list)
                items{k} = sprintf('%s   (%s)', list(k).name, strjoin(arrayfun(@(v) MainGuiStyle.num(v, 2), ...
                    list(k).pos(:)', 'UniformOutput', false), ', '));
            end
            set(obj.popSaved, 'String', items, 'Value', min(max(1, obj.popSaved.Value), numel(items)));
            obj.btnSavedDelete.Enable = 'on';
            obj.checkMovementEnabled();
        end

        function btnSavedAddCallback(obj)
            stage = obj.getStage();
            pos = stage.Pos(obj.stageAxes);
            list = obj.savedList;
            default = sprintf('P%d', numel(list) + 1);
            answer = inputdlg(sprintf('Name for the position (%s) %s:', ...
                strjoin(arrayfun(@(v) MainGuiStyle.num(v), pos, 'UniformOutput', false), ', '), MainGuiStyle.micron), ...
                'Save position', 1, {default});
            if isempty(answer) || isempty(strtrim(answer{1})); return; end
            list(end + 1) = struct('name', strtrim(answer{1}), 'pos', pos);
            obj.saveList(list);
            obj.popSaved.Value = numel(list);
        end

        function btnSavedGoCallback(obj)
            list = obj.savedList;
            k = obj.popSaved.Value;
            if isempty(list) || k > numel(list); return; end
            stage = obj.getStage();
            obj.remember();
            stage.move(obj.stageAxes, list(k).pos(1:numel(obj.stageAxes)));
        end

        function btnSavedDeleteCallback(obj)
            list = obj.savedList;
            k = obj.popSaved.Value;
            if isempty(list) || k > numel(list); return; end
            if ~QuestionUserYesNo('Delete saved position', sprintf('Delete the saved position "%s"?', list(k).name))
                return
            end
            list(k) = [];
            obj.saveList(list);
        end
    end

    methods (Access = private)
        function stage = getStage(obj)
            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
        end
    end

    %% overridden from EventListener
    methods
        function onEvent(obj, event)
            if event.isError ...
                    || isfield(event.extraInfo, ClassStage.EVENT_STEP_SIZE_CHANGED) ...
                    || isfield(event.extraInfo, ClassStage.EVENT_POSITION_CHANGED)
                obj.refresh();
            elseif isfield(event.extraInfo, ClassStage.EVENT_STAGE_AVAILABLITY_CHANGED) ...
                    || isfield(event.extraInfo, StageScanner.EVENT_SCAN_STARTED) ...
                    || isfield(event.extraInfo, StageScanner.EVENT_SCAN_FINISHED) ...
                    || isfield(event.extraInfo, Experiment.EVENT_EXP_PAUSED) ...
                    || isfield(event.extraInfo, Experiment.EVENT_EXP_RESUMED)
                obj.checkMovementEnabled;
            end
        end
    end
end
