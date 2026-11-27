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
        TIPS = {'|B|: field magnitude, in gauss', ...
                [char(952) ' (Btheta): POLAR angle from +z, 0 to 180' char(176)], ...
                [char(966) ' (Bphi): AZIMUTH from +x toward +y, -180 to 180' char(176)]};
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
        hboxCurB            % 1x3 hbox holding those two
        lowerLimLbl         % 1x3
        upperLimLbl         % 1x3
        edtCurBValue        % 1x3: the value typed in edtCurB, inf when none
        panelMain
        titleText           % the section's title (it names the frame)
        helm
    end

    methods
        function obj = ViewBsphereControl(parent, controller)
            expNames = Experiment.getExperimentNames();
            namesToListenTo = {expNames{:}, Helmholtz.NAME}; %#ok<CCAT>

            obj@GuiComponent(parent, controller);
            obj@EventListener(namesToListenTo);
            obj.helm = getObjByName('Helmholtz');

            %%%% layout (shared with the cartesian and current sections) %%%%
            deg = [' (' char(176) ')'];
            L = HelmholtzAxesLayout.build(parent.component, obj.sectionTitle(), ...
                {'|B| (G)', [char(952) deg], [char(966) deg]}, obj.TIPS, ['G, ' char(176)], ...
                'How far − and + move |B| (G) and the angles (degrees)');
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

            %%%% Callbacks %%%%
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
            currentB = obj.helm.Bsphere;
            for i = 1:3
                obj.tvCurB(i).String = sprintf('%.3f', currentB(i));
                [lo, hi] = obj.helm.limits(obj.AXES{i});
                obj.lowerLimLbl(i).String = ViewBsphereControl.fmtLim(lo);
                obj.upperLimLbl(i).String = ViewBsphereControl.fmtLim(hi);
            end
            obj.titleText.String = obj.sectionTitle();
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
                EventStation.anonymousWarning('%s', sprintf(['%s must be a number the coils can reach ' ...
                    'with the other two held: %s to %s.'], obj.AXES{index}, ...
                    ViewBsphereControl.fmtLim(lo), ViewBsphereControl.fmtLim(hi)));
            end
        end

        function clearEdtCurB(obj, index)
            obj.edtCurB(index).String = '';
            obj.edtCurBValue(index) = inf;
            TargetEntry.showPending(obj.edtCurB(index), false);
        end
    end

    methods (Access = private)
        function s = sectionTitle(obj)
            % 'Field, spherical  ·  θ from z, φ from x  ·  nominal coil frame'
            s = sprintf('Field, spherical   %s   %s from z, %s from x   %s   %s', char(183), ...
                char(952), char(966), char(183), obj.helm.frameName);
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
