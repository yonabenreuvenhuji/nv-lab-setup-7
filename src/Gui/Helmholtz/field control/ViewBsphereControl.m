classdef ViewBsphereControl < GuiComponent & EventListener
    %VIEWBSPHERECONTROL |B|, theta and phi, with the limits each can move between.
    %
    % theta is the POLAR angle from +z, phi the AZIMUTH from +x (physics convention); an
    % NV axis is theta 54.7356, phi 45 (+k*90) once the coils are calibrated. The limits
    % are recomputed after every change: each is how far that one coordinate can go with
    % the other two held (Helmholtz.limits), which is not a fixed box.

    properties (Constant)
        AXES = {'Brho', 'Btheta', 'Bphi'};
        LABELS = {'|B|', 'theta', 'phi'};
        TIPS = {'|B|, Gauss', ...
                'Btheta: POLAR angle from +z, 0..180 deg', ...
                'Bphi: AZIMUTH from +x toward +y, -180..180 deg'};
    end

    properties
        btnMoveLeft         % 1x3 button
        btnMoveRight        % 1x3 button
        btnMoveToBlue       % button
        isBlueMovingAvailable   % logical
        edtSteps            % edit-text
        tvCurB              % 1x3 text-view for the present field
        edtCurB             % 1x3 edit-input for a new value
        hboxCurB            % 1x3 hbox holding those two
        lowerLimLbl         % 1x3
        upperLimLbl         % 1x3
        edtCurBValue        % 1x3: the value typed in edtCurB, inf when none
        panelMain
        helm
    end

    methods
        function obj = ViewBsphereControl(parent, controller)
            expNames = Experiment.getExperimentNames();
            namesToListenTo = {expNames{:}, Helmholtz.NAME}; %#ok<CCAT>

            obj@GuiComponent(parent, controller);
            obj@EventListener(namesToListenTo);
            obj.helm = getObjByName('Helmholtz');

            %%%% panel init %%%%
            obj.panelMain = uix.Panel('Parent', parent.component, 'Title', ...
                'Spherical (G, deg): theta = polar from z, phi = azimuth from x', 'Padding', 5);
            hboxMain = uix.HBox('Parent', obj.panelMain, 'Spacing', 25, 'Padding', 0);
            vboxLeft = uix.VBox('Parent',hboxMain, 'Spacing', 6);

            %%%% step %%%%
            hboxStep = uix.HBox('Parent', vboxLeft, 'Spacing', 5);
            uicontrol(obj.PROP_LABEL{:}, 'Parent', hboxStep, 'String', 'Step:');
            obj.edtSteps = uicontrol(obj.PROP_EDIT{:}, 'Parent', hboxStep);
            hboxStep.Widths = [-1 -1];
            stepHeightNoSpacing = 25;

            %%%% the grid (filled column by column) %%%%
            gridLeftSide = uix.Grid('Parent', vboxLeft, 'Spacing', 5);
            for i = 1:3
                uicontrol(obj.PROP_LABEL{:}, 'Parent', gridLeftSide, 'String', obj.LABELS{i}, ...
                    'TooltipString', obj.TIPS{i});
            end
            obj.lowerLimLbl = gobjects(1, 3);
            for i = 1:3
                obj.lowerLimLbl(i) = uicontrol(obj.PROP_CHANGED_LABEL{:}, 'Parent', gridLeftSide, ...
                    'String', '', 'TooltipString', 'lowest value with the other two held');
            end
            obj.btnMoveLeft = gobjects(1, 3);
            for i = 1:3
                obj.btnMoveLeft(i) = uicontrol(obj.PROP_BUTTON{:}, 'Parent', gridLeftSide, ...
                    'String', StringHelper.LEFT_ARROW, 'FontSize', 20);
            end
            obj.tvCurB = gobjects(1, 3);
            obj.edtCurB = gobjects(1, 3);
            obj.hboxCurB = gobjects(1, 3);
            obj.edtCurBValue = inf * ones(1, 3);
            for i = 1:3
                obj.hboxCurB(i) = uix.HBox('Parent', gridLeftSide, 'Padding', 0, 'Spacing', 0);
                obj.tvCurB(i) = uicontrol(obj.PROP_EDIT{:}, 'Parent', obj.hboxCurB(i), ...
                    'TooltipString', obj.TIPS{i});
                obj.edtCurB(i) = uicontrol(obj.PROP_EDIT{:}, 'Parent', obj.hboxCurB(i), ...
                    'ForegroundColor', 'blue', 'BackgroundColor', [0.8 0.9 1]);
                obj.showEdtCurB(i, false);
            end
            obj.btnMoveRight = gobjects(1, 3);
            for i = 1:3
                obj.btnMoveRight(i) = uicontrol(obj.PROP_BUTTON{:}, 'Parent', gridLeftSide, ...
                    'String', StringHelper.RIGHT_ARROW, 'FontSize', 20);
            end
            obj.upperLimLbl = gobjects(1, 3);
            for i = 1:3
                obj.upperLimLbl(i) = uicontrol(obj.PROP_CHANGED_LABEL{:}, 'Parent', gridLeftSide, ...
                    'String', '', 'TooltipString', 'highest value with the other two held');
            end

            gridWidths = [55 60 35 70 35 60];
            heightGridTotalNoSpacing = 90;
            lineHeight = heightGridTotalNoSpacing / 3;
            set(gridLeftSide, 'Widths', -gridWidths, 'Heights', lineHeight * ones(1, 3));
            gridLeftHeightTotal = heightGridTotalNoSpacing + 5 * 3;
            gridLeftWidthTotal = sum(gridWidths) + 4 * length(gridWidths);
            vboxLeft.Heights = [stepHeightNoSpacing, heightGridTotalNoSpacing+10];
            vboxLeftHeightTotal = gridLeftHeightTotal + stepHeightNoSpacing + 42;

            %%%% 'move to' button %%%%
            vboxRight = uix.VBox('Parent',hboxMain, 'Spacing', 6);
            obj.btnMoveToBlue = uicontrol(obj.PROP_BUTTON{:}, 'Parent', vboxRight, ...
                'String', 'Change To Blue', 'Enable', 'off', ...
                'TooltipString', 'Set the field to the values in blue (all at once)');
            vboxRight.Heights = -2;
            rightSideTotalWidth = 120;
            widthsMain = [gridLeftWidthTotal, rightSideTotalWidth];
            hboxMain.Widths = -widthsMain;

            %%%% Callbacks %%%%
            obj.edtSteps.Callback = @(h,e) obj.edtStepSizeCallback();
            obj.btnMoveToBlue.Callback = @(h,e) obj.btnMoveToBlueCallback();
            obj.isBlueMovingAvailable = false;
            for i = 1:3
                obj.tvCurB(i).Callback = @(h,e) obj.tvCurBCallback(i);
                obj.edtCurB(i).Callback = @(h,e) obj.checkEdtCurBValue(i);
                obj.btnMoveLeft(i).Callback = @(h,e) obj.btnMoveCallback(i, true);
                obj.btnMoveRight(i).Callback = @(h,e) obj.btnMoveCallback(i, false);
                obj.lowerLimLbl(i).Callback = @(h,e) obj.limEdtCallback();
                obj.upperLimLbl(i).Callback = @(h,e) obj.limEdtCallback();
            end

            obj.width = sum(widthsMain) + (hboxMain.Spacing -3) * length(widthsMain);
            obj.height = vboxLeftHeightTotal;
            obj.refresh();
        end

        function refresh(obj)
            obj.edtSteps.String = StringHelper.formatNumber(obj.helm.stepSize);
            currentB = obj.helm.Bsphere;
            for i = 1:3
                obj.tvCurB(i).String = sprintf('%.3f', currentB(i));
                [lo, hi] = obj.helm.limits(obj.AXES{i});
                obj.lowerLimLbl(i).String = ViewBsphereControl.fmtLim(lo);
                obj.upperLimLbl(i).String = ViewBsphereControl.fmtLim(hi);
            end
            obj.panelMain.Title = sprintf(['Spherical (G, deg): theta = polar from z, ' ...
                'phi = azimuth from x  [%s]'], obj.helm.frameName);
            obj.btnMoveToBlue.Enable = BooleanHelper.boolToOnOff(obj.isBlueMovingAvailable);
            drawnow
        end

        function edtStepSizeCallback(obj)
            stepSizeString = obj.edtSteps.String;
            if ~ValidationHelper.isStringValueANumber(stepSizeString)
                obj.edtSteps.String = StringHelper.formatNumber(obj.helm.stepSize);
                EventStation.anonymousError('Step size must be a valid number! Reverting.');
            end
            obj.helm.stepSize = str2double(stepSizeString);
        end

        function btnMoveCallback(obj, index, isLeft)
            step = BooleanHelper.ifTrueElse(isLeft, -1, 1) * obj.helm.stepSize;
            obj.helm.relativeBsphereStep(index, step);
        end

        function btnMoveToBlueCallback(obj)
            % All typed values at once, so the field goes straight to the target rather than
            % through intermediate settings that might be out of reach.
            target = obj.helm.Bsphere;
            typed = obj.edtCurBValue ~= inf;
            target(typed) = obj.edtCurBValue(typed);
            for i = find(typed); obj.clearEdtCurB(i); end
            obj.isBlueMovingAvailable = false;
            obj.helm.setBsphere(target(1), target(2), target(3));
        end

        function tvCurBCallback(obj, index)
            str = obj.tvCurB(index).String;
            obj.tvCurB(index).String = sprintf('%.3f', obj.helm.Bsphere(index));
            if obj.isSettable(index, str)
                obj.edtCurB(index).String = str;
                obj.edtCurBValue(index) = str2double(str);
                obj.showEdtCurB(index, true);
                obj.isBlueMovingAvailable = true;
                obj.refresh
            end
        end

        function limEdtCallback(obj) %#ok<MANU>
            EventStation.anonymousWarning(['The limits follow from the supplies'' current ' ...
                'range and the coil calibration; they cannot be typed.']);
        end

        function showEdtCurB(obj, index, shouldShow)
            if shouldShow
                obj.hboxCurB(index).Widths = [-1 -1];
            else
                obj.hboxCurB(index).Widths = [-1 0];
            end
        end

        function checkEdtCurBValue(obj, index)
            str = obj.edtCurB(index).String;
            if isempty(str)
                obj.clearEdtCurB(index);
                if all(obj.edtCurBValue == inf)
                    obj.isBlueMovingAvailable = false;
                    obj.btnMoveToBlue.Enable = 'off';
                end
                return
            end
            if obj.isSettable(index, str)
                obj.edtCurBValue(index) = str2double(str);
            end
            obj.edtCurB(index).String = StringHelper.formatNumber(obj.edtCurBValue(index));
        end

        function ok = isSettable(obj, index, str)
            ok = ValidationHelper.isStringValueANumber(str);
            if ok; ok = obj.helm.checkLimits(obj.AXES{index}, str2double(str), 0); end
            if ~ok
                [lo, hi] = obj.helm.limits(obj.AXES{index});
                EventStation.anonymousWarning('%s', sprintf(['%s must be a number the coils can reach ' ...
                    'with the other two held: %s to %s.'], obj.AXES{index}, ...
                    ViewBsphereControl.fmtLim(lo), ViewBsphereControl.fmtLim(hi)));
            end
        end

        function clearEdtCurB(obj, index)
            obj.edtCurB(index).String = '';
            obj.edtCurBValue(index) = inf;
            obj.showEdtCurB(index, false);
            set(gcf, 'CurrentObject', gcf);
        end
    end

    methods (Static)
        function s = fmtLim(v)
            if isnan(v); s = '--'; else; s = sprintf('%.2f', v); end
        end
    end

    %% overridden from EventListener
    methods
        function onEvent(obj, event)
            if event.isError ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_STEP_SIZE_CHANGED) ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_MAGNETIC_FIELD_CHANGED) ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_CALIBRATION_CHANGED)
                obj.refresh();
            end
        end
    end
end
