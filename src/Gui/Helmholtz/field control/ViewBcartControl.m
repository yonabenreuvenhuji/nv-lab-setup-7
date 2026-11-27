classdef ViewBcartControl < GuiComponent & EventListener
    %VIEWBCARTCONTROL Bx, By, Bz, with the limits each can move between.
    %
    % In the NV crystal frame once the coils are calibrated (nominal coil frame before).
    % Each limit is how far that one component can go with the other two held -- the
    % reachable fields are a parallelepiped, not a box, so these move with the field.

    properties (Constant)
        AXES = {'Bx', 'By', 'Bz'};
        TIPS = {'x component of B, in gauss (in the frame named in the title)', ...
                'y component of B, in gauss (in the frame named in the title)', ...
                'z component of B, in gauss (in the frame named in the title)'};
    end

    properties
        btnMoveLeft         % 1x3 button
        btnMoveRight        % 1x3 button
        btnMoveToBlue       % button: "Set to target"
        btnClearTargets     % button
        isBlueMovingAvailable   % logical: some target is pending
        edtSteps            % edit-text
        tvCurB              % 1x3 text-view for the present field
        edtCurB             % 1x3 edit-input for a new value
        hboxCurB            % 1x3
        lowerLimLbl         % 1x3
        upperLimLbl         % 1x3
        edtCurBValue        % 1x3: the value typed in edtCurB, inf when none
        panelMain
        titleText           % the section's title (it names the frame)
        helm
    end

    methods
        function obj = ViewBcartControl(parent, controller)
            expNames = Experiment.getExperimentNames();
            namesToListenTo = {expNames{:}, Helmholtz.NAME}; %#ok<CCAT>

            obj@GuiComponent(parent, controller);
            obj@EventListener(namesToListenTo);
            obj.helm = getObjByName('Helmholtz');

            %%%% layout (shared with the spherical and current sections) %%%%
            L = HelmholtzAxesLayout.build(parent.component, sprintf('Field, cartesian   %s   %s', char(183), obj.helm.frameName), ...
                {'Bx (G)', 'By (G)', 'Bz (G)'}, obj.TIPS, 'G', ...
                'How far − and + move a component (G)');
            obj.panelMain = L.panel;
            obj.titleText = L.title;
            obj.component = L.panel;
            obj.edtSteps = L.edtSteps;
            obj.lowerLimLbl = L.lowerLimLbl;
            obj.upperLimLbl = L.upperLimLbl;
            obj.btnMoveLeft = L.btnMoveLeft;
            obj.btnMoveRight = L.btnMoveRight;
            obj.tvCurB = L.tvCurB;
            obj.edtCurB = L.edtCurB;
            obj.hboxCurB = L.hboxCurB;
            obj.btnMoveToBlue = L.btnMoveToBlue;
            obj.btnClearTargets = L.btnClearTargets;
            obj.btnClearTargets.Callback = @(h,e) obj.clearAllTargets();
            obj.edtCurBValue = inf * ones(1, 3);
            for i = 1:3
                TargetEntry.styleNow(obj.tvCurB(i));
                set(obj.tvCurB(i), 'TooltipString', obj.TIPS{i});
                TargetEntry.styleTarget(obj.edtCurB(i));
            end

            obj.edtSteps.Callback = @(h,e) obj.edtStepSizeCallback();
            obj.btnMoveToBlue.Callback = @(h,e) obj.btnMoveToBlueCallback();
            obj.isBlueMovingAvailable = false;
            for i = 1:3
                obj.edtCurB(i).Callback = @(h,e) obj.checkEdtCurBValue(i);
                obj.btnMoveLeft(i).Callback = @(h,e) obj.btnMoveCallback(i, true);
                obj.btnMoveRight(i).Callback = @(h,e) obj.btnMoveCallback(i, false);
            end

            obj.width = L.width;
            obj.height = L.height;
            obj.refresh();
        end

        function refresh(obj)
            obj.edtSteps.String = StringHelper.formatNumber(obj.helm.stepSize);
            currentB = obj.helm.B;
            for i = 1:3
                obj.tvCurB(i).String = sprintf('%.3f', currentB(i));
                [lo, hi] = obj.helm.limits(obj.AXES{i});
                obj.lowerLimLbl(i).String = ViewBsphereControl.fmtLim(lo);
                obj.upperLimLbl(i).String = ViewBsphereControl.fmtLim(hi);
            end
            obj.titleText.String = sprintf('Field, cartesian   %s   %s', char(183), obj.helm.frameName);
            TargetEntry.paintApply(obj.btnMoveToBlue, obj.isBlueMovingAvailable);
            obj.btnClearTargets.Enable = BooleanHelper.boolToOnOff(obj.isBlueMovingAvailable);
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
            obj.helm.relativeBcartStep(index, step);
        end

        function btnMoveToBlueCallback(obj)
            target = obj.helm.B;
            typed = obj.edtCurBValue ~= inf;
            target(typed) = obj.edtCurBValue(typed);
            for i = find(typed); obj.clearEdtCurB(i); end
            obj.isBlueMovingAvailable = false;
            obj.helm.setB(target);
        end

        function tvCurBCallback(obj, index)
            str = obj.tvCurB(index).String;
            obj.tvCurB(index).String = sprintf('%.3f', obj.helm.B(index));
            if obj.isSettable(index, str)
                obj.edtCurB(index).String = str;
                obj.edtCurBValue(index) = str2double(str);
                obj.showEdtCurB(index, true);
                obj.isBlueMovingAvailable = true;
                obj.refresh
            end
        end

        function showEdtCurB(~, ~, ~)
            % now | target, both always shown
            % (now and target are separate columns: nothing to show or hide)
        end

        function checkEdtCurBValue(obj, index)
            % Enter in a target: pending (amber) if it is reachable; empty = no target
            str = strtrim(obj.edtCurB(index).String);
            if isempty(str)
                obj.clearEdtCurB(index);
            elseif obj.isSettable(index, str)
                obj.edtCurBValue(index) = str2double(str);
                obj.edtCurB(index).String = StringHelper.formatNumber(obj.edtCurBValue(index));
                TargetEntry.showPending(obj.edtCurB(index), true);
            elseif isinf(obj.edtCurBValue(index))
                obj.clearEdtCurB(index);
            else
                obj.edtCurB(index).String = StringHelper.formatNumber(obj.edtCurBValue(index));
            end
            obj.isBlueMovingAvailable = any(~isinf(obj.edtCurBValue));
            obj.refresh();
        end

        function clearAllTargets(obj)
            for i = 1:3; obj.clearEdtCurB(i); end
            obj.isBlueMovingAvailable = false;
            obj.refresh();
        end

        function ok = isSettable(obj, index, str)
            ok = ValidationHelper.isStringValueANumber(str);
            if ok; ok = obj.helm.checkLimits(obj.AXES{index}, str2double(str), 0); end
            if ~ok
                [lo, hi] = obj.helm.limits(obj.AXES{index});
                EventStation.anonymousWarning('%s', sprintf(['%s must be a number the coils can ' ...
                    'reach with the other two held: %s to %s G.'], obj.AXES{index}, ...
                    ViewBsphereControl.fmtLim(lo), ViewBsphereControl.fmtLim(hi)));
            end
        end

        function clearEdtCurB(obj, index)
            obj.edtCurB(index).String = '';
            obj.edtCurBValue(index) = inf;
            TargetEntry.showPending(obj.edtCurB(index), false);
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
