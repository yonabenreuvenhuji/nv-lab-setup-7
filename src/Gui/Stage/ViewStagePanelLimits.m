classdef ViewStagePanelLimits < GuiComponent & EventListener
    %VIEWSTAGEPANELLIMITS The "Soft limits" section: how far the stage may move.
    %
    %          Lower (µm)        Upper (µm)        Hard range (µm)     Set both limits
    %   X   [ -6500 ] [Here]   [ 6500 ] [Here]    -6500 … 6500       [Full] [± 10 µm]
    %   Y   ...
    %   ± distance [ 10 ] µm      Here = where the stage is now
    %
    %   typing a limit   sets it
    %   Here             sets that one limit to the stage's position now
    %   Full             both limits = the hard range (the whole travel of the axis)
    %   ± d              lower = stage − d, upper = stage + d
    % (They replace "Set [lower][upper] of [X][Y][Z] to: Hard limits / Stage position /
    % Stage position ±", which needed boxes ticked before a button did anything.)
    % Moving a limit past the stage position asks first, because the stage will move.

    properties
        edtLower            % 1xN edit
        edtUpper            % 1xN edit
        btnLowerHere        % 1xN
        btnUpperHere        % 1xN
        tvHard              % 1xN text
        btnFull             % 1xN
        btnAround           % 1xN
        edtSetAround        % edit: the ± distance

        stageName           % string
        stageAxes           % string. example: "xy"
    end

    properties (Constant)
        SET_AROUND_DEFAULT_STRING = '10';
        LIM_LOWER = 0;
        LIM_UPPER = 1;
        % axis | lower | here | gap | upper | here | gap | hard | gap | full | ±
        COL_WIDTHS = [20 70 46 8 70 46 8 104 8 46 70]
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
            sp = 4;
            W = obj.COL_WIDTHS;

            [panel, inner] = MainGuiStyle.section(parent.component, 'Soft limits');
            obj.component = panel;

            MainGuiStyle.groupTitles(inner, W, sp, {[2 3], [5 6], 8, [10 11]}, ...
                {sprintf('Lower (%s)', um), sprintf('Upper (%s)', um), sprintf('Hard range (%s)', um), 'Set both limits'});

            gridRow = MainGuiStyle.hbox(inner, 0);
            grid = MainGuiStyle.grid(gridRow, sp);
            MainGuiStyle.empty(gridRow);
            obj.edtLower = gobjects(1, n); obj.edtUpper = gobjects(1, n);
            obj.btnLowerHere = gobjects(1, n); obj.btnUpperHere = gobjects(1, n);
            obj.tvHard = gobjects(1, n); obj.btnFull = gobjects(1, n); obj.btnAround = gobjects(1, n);
            % (a uix.Grid fills column by column)
            for i = 1:n; MainGuiStyle.axisTag(grid, obj.stageAxes(i)); end
            for i = 1:n
                obj.edtLower(i) = MainGuiStyle.edit(grid, '', 'TooltipString', ...
                    sprintf('Lowest position %s may go to. Type to set it.', upper(obj.stageAxes(i))));
            end
            for i = 1:n
                obj.btnLowerHere(i) = MainGuiStyle.button(grid, 'Here', 'quiet', 'TooltipString', ...
                    sprintf('%s lower limit = where the stage is now', upper(obj.stageAxes(i))));
            end
            for i = 1:n; MainGuiStyle.empty(grid); end
            for i = 1:n
                obj.edtUpper(i) = MainGuiStyle.edit(grid, '', 'TooltipString', ...
                    sprintf('Highest position %s may go to. Type to set it.', upper(obj.stageAxes(i))));
            end
            for i = 1:n
                obj.btnUpperHere(i) = MainGuiStyle.button(grid, 'Here', 'quiet', 'TooltipString', ...
                    sprintf('%s upper limit = where the stage is now', upper(obj.stageAxes(i))));
            end
            for i = 1:n; MainGuiStyle.empty(grid); end
            for i = 1:n
                phAxis = ClassStage.getAxis(obj.stageAxes(i));
                [hardNeg, hardPos] = stage.ReturnHardLimits(phAxis);
                obj.tvHard(i) = MainGuiStyle.hint(grid, sprintf('%s %s %s', MainGuiStyle.num(hardNeg), char(8230), ...
                    MainGuiStyle.num(hardPos)), 'center', 'TooltipString', 'The whole travel of the axis');
            end
            for i = 1:n; MainGuiStyle.empty(grid); end
            for i = 1:n
                obj.btnFull(i) = MainGuiStyle.button(grid, 'Full', 'normal', 'TooltipString', ...
                    sprintf('%s: both limits = the hard range (the whole travel)', upper(obj.stageAxes(i))));
            end
            for i = 1:n
                obj.btnAround(i) = MainGuiStyle.button(grid, '', 'normal');
            end
            set(grid, 'Widths', W, 'Heights', S.ROW_H * ones(1, n));
            gridRow.Widths = [sum(W) + sp*numel(W), -1];

            % The distance for the ± buttons
            bottom = MainGuiStyle.hbox(inner, 6);
            MainGuiStyle.label(bottom, sprintf('%s distance', char(177)), 'left');
            obj.edtSetAround = MainGuiStyle.edit(bottom, obj.SET_AROUND_DEFAULT_STRING, ...
                'TooltipString', sprintf('How far below and above the stage position the %s buttons put the limits', char(177)));
            MainGuiStyle.label(bottom, um, 'left');
            MainGuiStyle.hint(bottom, 'Here = where the stage is now', 'right');
            bottom.Widths = [62 56 24 -1];

            inner.Heights = [15, n*S.ROW_H + sp*(n-1), S.ROW_H];
            obj.height = MainGuiStyle.sectionHeight(MainGuiStyle.rowsHeight(inner.Heights));
            obj.width = sum(W) + sp*numel(W) + 20;

            for i = 1 : n
                obj.edtLower(i).Callback = @(h,e) obj.edtLowerLimitCallback(i);
                obj.edtUpper(i).Callback = @(h,e) obj.edtUpperLimitCallback(i);
                obj.btnLowerHere(i).Callback = @(h,e) obj.setHere(i, obj.LIM_LOWER);
                obj.btnUpperHere(i).Callback = @(h,e) obj.setHere(i, obj.LIM_UPPER);
                obj.btnFull(i).Callback = @(h,e) obj.setFull(i);
                obj.btnAround(i).Callback = @(h,e) obj.setAround(i);
            end
            obj.edtSetAround.Callback = @(h,e) obj.edtSetAroundCallback();

            obj.refresh();
            obj.labelAroundButtons();
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

        function setHere(obj, index, which)
            % One limit of one axis = the stage position now
            stage = getObjByName(obj.stageName);
            phAxis = ClassStage.getAxis(obj.stageAxes(index));
            stage.setLim(stage.Pos(phAxis), phAxis, which);
        end

        function setFull(obj, index)
            % Both limits = the hard range. (Upper first when widening, so the lower limit
            % never has to pass the upper one.)
            stage = getObjByName(obj.stageName);
            phAxis = ClassStage.getAxis(obj.stageAxes(index));
            [hardNeg, hardPos] = stage.ReturnHardLimits(phAxis);
            stage.setLim(hardPos, phAxis, obj.LIM_UPPER);
            stage.setLim(hardNeg, phAxis, obj.LIM_LOWER);
        end

        function setAround(obj, index)
            % lower = stage − d, upper = stage + d (clipped to the hard range)
            stage = getObjByName(obj.stageName);
            phAxis = ClassStage.getAxis(obj.stageAxes(index));
            d = obj.aroundDistance();
            pos = stage.Pos(phAxis);
            [hardNeg, hardPos] = stage.ReturnHardLimits(phAxis);
            lower = max(pos - d, hardNeg);
            upper = min(pos + d, hardPos);
            [limNeg, ~] = stage.ReturnLimits(phAxis);
            if upper >= limNeg      % (order so the lower limit never passes the upper one)
                stage.setLim(upper, phAxis, obj.LIM_UPPER);
                stage.setLim(lower, phAxis, obj.LIM_LOWER);
            else
                stage.setLim(lower, phAxis, obj.LIM_LOWER);
                stage.setLim(upper, phAxis, obj.LIM_UPPER);
            end
        end

        function edtSetAroundCallback(obj)
            if ~ValidationHelper.isValueNonNegative(obj.edtSetAround.String)
                obj.edtSetAround.String = obj.SET_AROUND_DEFAULT_STRING;
                EventStation.anonymousWarning('The distance must be a non-negative number! Reverting.');
            end
            obj.labelAroundButtons();
        end
    end

    methods (Access = private)
        function d = aroundDistance(obj)
            d = str2double(obj.edtSetAround.String);
            if ~isfinite(d) || d < 0; d = str2double(obj.SET_AROUND_DEFAULT_STRING); end
        end

        function labelAroundButtons(obj)
            d = MainGuiStyle.num(obj.aroundDistance());
            for i = 1:numel(obj.btnAround)
                set(obj.btnAround(i), 'String', sprintf('%s %s', char(177), d), 'TooltipString', ...
                    sprintf('%s: lower = stage %s %s, upper = stage + %s %s', upper(obj.stageAxes(i)), ...
                    char(8722), d, d, MainGuiStyle.micron));
            end
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
