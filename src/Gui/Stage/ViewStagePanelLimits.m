classdef ViewStagePanelLimits < GuiComponent & EventListener
    %VIEWSTAGEPANELLIMITS The "Soft limits" section: how far the stage may move.
    %
    %   Axis  Lower (µm)  Upper (µm)   Hard range (µm)
    %    X    [ -6500 ]   [ 6500  ]    -6500 … 6500
    %   Set [x] lower [x] upper  of  [x] X [x] Y [x] Z  to:
    %   [Hard limits]  [Stage position]  [Stage position ± [10] µm]
    %
    % Typing a limit sets it at once. The buttons set the ticked limits of the ticked axes.
    % Moving a limit past the stage position asks first, because the stage will move.

    properties
        edtLower            % 1xN edit-text
        edtUpper            % 1xN edit-text
        tvHard              % 1xN text
        cbxAxis             % 1xN checkbox
        cbxUpperLim         % checkbox
        cbxLowerLim         % checkbox
        btnSetToMax         % button
        btnSetToCurPos      % button
        btnSetToAround      % button
        edtSetAround        % edit-text

        stageName           % string
        stageAxes           % string. example: "xy"
    end

    properties (Constant)
        SET_AROUND_DEFAULT_STRING = '10';
        LIM_LOWER = 0;
        LIM_UPPER = 1;
    end

    methods
        function obj = ViewStagePanelLimits(parent, controller, stage)
            obj@GuiComponent(parent, controller);
            obj@EventListener(stage.name);
            obj.stageName = stage.name;
            obj.stageAxes = stage.availableAxes;
            n = length(obj.stageAxes);
            S = MainGuiStyle;
            um = MainGuiStyle.micron;

            [panel, inner] = MainGuiStyle.section(parent.component, 'Soft limits');
            obj.component = panel;

            gridRow = MainGuiStyle.hbox(inner, 0);
            grid = MainGuiStyle.grid(gridRow, 4);
            MainGuiStyle.empty(gridRow);
            MainGuiStyle.empty(grid);
            for i = 1:n; MainGuiStyle.axisTag(grid, obj.stageAxes(i)); end
            MainGuiStyle.header(grid, sprintf('Lower (%s)', um));
            obj.edtLower = gobjects(1, n);
            for i = 1:n
                obj.edtLower(i) = MainGuiStyle.edit(grid, '', 'TooltipString', 'Lowest position the stage may go to');
            end
            MainGuiStyle.header(grid, sprintf('Upper (%s)', um));
            obj.edtUpper = gobjects(1, n);
            for i = 1:n
                obj.edtUpper(i) = MainGuiStyle.edit(grid, '', 'TooltipString', 'Highest position the stage may go to');
            end
            MainGuiStyle.header(grid, sprintf('Hard range (%s)', um), 'left');
            obj.tvHard = gobjects(1, n);
            for i = 1:n
                phAxis = ClassStage.getAxis(obj.stageAxes(i));
                [hardNeg, hardPos] = stage.ReturnHardLimits(phAxis);
                obj.tvHard(i) = MainGuiStyle.hint(grid, sprintf('%s %s %s', MainGuiStyle.num(hardNeg), char(8230), ...
                    MainGuiStyle.num(hardPos)), 'left', 'TooltipString', 'The physical range of the axis');
            end
            widths = [20 76 76 130];
            set(grid, 'Widths', widths, 'Heights', [S.HEAD_H, S.ROW_H * ones(1, n)]);
            gridRow.Widths = [sum(widths) + 4*numel(widths), -1];

            %%%% which limits %%%%
            pick = MainGuiStyle.hbox(inner, 6);
            MainGuiStyle.label(pick, 'Set', 'left');
            obj.cbxLowerLim = MainGuiStyle.checkbox(pick, 'lower', false);
            obj.cbxUpperLim = MainGuiStyle.checkbox(pick, 'upper', false);
            MainGuiStyle.label(pick, 'of', 'center');
            obj.cbxAxis = gobjects(1, n);
            for i = 1:n
                obj.cbxAxis(i) = MainGuiStyle.checkbox(pick, upper(obj.stageAxes(i)), false);
            end
            MainGuiStyle.label(pick, 'to:', 'left');
            pick.Widths = [26 56 58 18 40*ones(1, n) -1];

            buttons = MainGuiStyle.hbox(inner, 6);
            obj.btnSetToMax = MainGuiStyle.button(buttons, 'Hard limits', 'normal', ...
                'TooltipString', 'The whole physical range');
            obj.btnSetToCurPos = MainGuiStyle.button(buttons, 'Stage position', 'normal', ...
                'TooltipString', 'Where the stage is now');
            obj.btnSetToAround = MainGuiStyle.button(buttons, sprintf('Stage position %s', char(177)), 'normal', ...
                'TooltipString', 'Where the stage is now, minus (lower) or plus (upper) the distance on the right');
            obj.edtSetAround = MainGuiStyle.edit(buttons, obj.SET_AROUND_DEFAULT_STRING);
            MainGuiStyle.label(buttons, um, 'left');
            buttons.Widths = [-1 -1 -1.1 56 24];

            inner.Heights = [S.HEAD_H + n*S.ROW_H + 4*n, S.ROW_H, S.ROW_H];
            obj.height = MainGuiStyle.sectionHeight(MainGuiStyle.rowsHeight(inner.Heights));
            obj.width = 520;

            for i = 1 : n
                obj.edtLower(i).Callback = @(h,e) obj.edtLowerLimitCallback(i);
                obj.edtUpper(i).Callback = @(h,e) obj.edtUpperLimitCallback(i);
            end
            obj.btnSetToAround.Callback = @(h,e) obj.btnSetToAroundCallback;
            obj.btnSetToCurPos.Callback = @(h,e) obj.btnSetToCurPosCallback;
            obj.btnSetToMax.Callback = @(h,e) obj.btnSetToMaxCallback;

            obj.refresh();
        end

        function refresh(obj)
            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
            for i = 1 : length(stage.availableAxes)
                phAxis = ClassStage.getAxis(stage.availableAxes(i));
                [limNeg, limPos] = stage.ReturnLimits(phAxis);
                obj.edtLower(i).String = MainGuiStyle.num(limNeg);
                obj.edtUpper(i).String = MainGuiStyle.num(limPos);
            end
        end

        function edtLowerLimitCallback(obj, index)
            obj.typedLimit(index, obj.LIM_LOWER, obj.edtLower(index));
        end

        function edtUpperLimitCallback(obj, index)
            obj.typedLimit(index, obj.LIM_UPPER, obj.edtUpper(index));
        end

        function typedLimit(obj, index, which, field)
            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
            phAxis = ClassStage.getAxis(obj.stageAxes(index));
            [limNeg, limPos] = stage.ReturnLimits(phAxis);
            [negHardLimit, posHardLimit] = stage.ReturnHardLimits(phAxis);
            value = str2double(field.String);
            msg = '';
            if ~ValidationHelper.isStringValueANumber(field.String)
                msg = 'A limit must be a number.';
            elseif ~ValidationHelper.isInBorders(value, negHardLimit, posHardLimit)
                msg = sprintf('A limit must be within the hard range [%s, %s].', MainGuiStyle.num(negHardLimit), MainGuiStyle.num(posHardLimit));
            elseif which == obj.LIM_LOWER && value > limPos
                msg = 'The lower limit cannot be above the upper one.';
            elseif which == obj.LIM_UPPER && value < limNeg
                msg = 'The upper limit cannot be below the lower one.';
            end
            if ~isempty(msg)
                obj.refresh();
                EventStation.anonymousWarning('%s Reverting.', msg);
                return
            end
            stage.setLim(value, phAxis, which);     % if it succeeds, an event refreshes us
        end

        function btnSetToAroundCallback(obj)
            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
            around = obj.edtSetAround.String;
            if ~ValidationHelper.isValueNonNegative(around)
                obj.edtSetAround.String = obj.SET_AROUND_DEFAULT_STRING;
                EventStation.anonymousWarning('The distance must be a non-negative number! Reverting.');
                return
            end
            around = str2double(around);
            obj.setTicked(@(phAxis, which) stage.Pos(phAxis) + BooleanHelper.ifTrueElse(which == obj.LIM_LOWER, -around, around));
        end

        function btnSetToCurPosCallback(obj)
            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
            obj.setTicked(@(phAxis, which) stage.Pos(phAxis));
        end

        function btnSetToMaxCallback(obj)
            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
            obj.setTicked(@(phAxis, which) obj.hardLimit(stage, phAxis, which));
        end
    end

    methods (Access = private)
        function setTicked(obj, valueFcn)
            % Set the ticked limits (lower/upper) of the ticked axes to valueFcn(axis, which)
            stage = getObjByName(obj.stageName);
            if ~any([obj.cbxAxis.Value]) || ~(obj.cbxLowerLim.Value || obj.cbxUpperLim.Value)
                EventStation.anonymousWarning('Tick which limits (lower / upper) and which axes to set first.');
                return
            end
            for i = 1 : length(obj.stageAxes)
                if ~obj.cbxAxis(i).Value; continue; end
                phAxis = ClassStage.getAxis(obj.stageAxes(i));
                if obj.cbxLowerLim.Value
                    stage.setLim(valueFcn(phAxis, obj.LIM_LOWER), phAxis, obj.LIM_LOWER);
                end
                if obj.cbxUpperLim.Value
                    stage.setLim(valueFcn(phAxis, obj.LIM_UPPER), phAxis, obj.LIM_UPPER);
                end
            end
        end

        function v = hardLimit(obj, stage, phAxis, which)
            [hardNeg, hardPos] = stage.ReturnHardLimits(phAxis);
            v = BooleanHelper.ifTrueElse(which == obj.LIM_LOWER, hardNeg, hardPos);
        end
    end

    %% overridden from EventListener
    methods
        function onEvent(obj, event)
            if event.isError || isfield(event.extraInfo, ClassStage.EVENT_LIM_CHANGED)
                obj.refresh();
            end
        end
    end
end
