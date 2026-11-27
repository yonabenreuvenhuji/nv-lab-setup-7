classdef ViewCurrentControl < GuiComponent & EventListener

    properties (Constant)
        AXES = {'Ix', 'Iy', 'Iz'};
    end
    
    properties
        btnMoveLeft         % 1x3 button
        btnMoveRight        % 1x3 button
        
        btnMoveToBlue       % button: "Set to target"
        btnClearTargets     % button
        isBlueMovingAvailable   % logical: some target is pending
        
        edtSteps            % edit-text 

        tvCurB            % 1x3 text-view for current fiekd
        edtCurB           % 1x3 edit-input for current field input
        hboxCurB          % 1x3 vbox that holds those ^ 2 views
        
        edtCurBValue      % 1x3 double that holds the values in edtCurB. 
                            % when an edt is visible, this always holds the
                            % numeric value shown. if it's invisible, this
                            % always stores Infinity
        helm
    end
    
    methods
        function obj = ViewCurrentControl(parent, controller)
            expNames = Experiment.getExperimentNames();
            namesToListenTo = {expNames{:}, Helmholtz.NAME}; %#ok<CCAT>
            
            obj@GuiComponent(parent, controller);
            obj@EventListener(namesToListenTo);
            obj.helm = getObjByName('Helmholtz');
            
            %%%% layout (shared with the spherical and cartesian sections) %%%%
            L = HelmholtzAxesLayout.build(parent.component, 'Coil currents', ...
                {'Ix (A)', 'Iy (A)', 'Iz (A)'}, ...
                {'Current in the x coil pair, in amperes', 'Current in the y coil pair, in amperes', ...
                'Current in the z coil pair, in amperes'}, 'A', 'How far − and + move a current (A)');
            obj.component = L.panel;
            obj.edtSteps = L.edtSteps;
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
                % the supply's range (fixed)
                L.lowerLimLbl(i).String = StringHelper.formatNumber(obj.helm.minCurrent(i));
                L.upperLimLbl(i).String = StringHelper.formatNumber(obj.helm.maxCurrent(i));
                set([L.lowerLimLbl(i), L.upperLimLbl(i)], 'TooltipString', 'The supply''s current range');
                TargetEntry.styleNow(obj.tvCurB(i));
                TargetEntry.styleTarget(obj.edtCurB(i));
            end
            
            %%%% Callbacks %%%%
            obj.edtSteps.Callback = @(h,e) obj.edtStepSizeCallback();
            obj.btnMoveToBlue.Callback = @(h,e) obj.btnMoveToBlueCallback();
            obj.isBlueMovingAvailable = false; % By default
            
            for i = 1:3
                obj.edtCurB(i).Callback = @(h,e) obj.checkEdtCurBValue(i);

                isLeft = true;
                obj.btnMoveLeft(i).Callback = @(h,e) obj.btnMoveCallback(i, isLeft);
                obj.btnMoveRight(i).Callback = @(h,e) obj.btnMoveCallback(i, ~isLeft);
            end
            
            
            obj.width = L.width;
            obj.height = L.height;
            
            obj.refresh();
        end
        
        function refresh(obj)
            obj.edtSteps.String = StringHelper.formatNumber(obj.helm.stepSizeI);
            currentI = obj.helm.I;   
            for i = 1:3
                obj.tvCurB(i).String = StringHelper.formatNumber(currentI(i));
            end
            TargetEntry.paintApply(obj.btnMoveToBlue, obj.isBlueMovingAvailable);
            obj.btnClearTargets.Enable = BooleanHelper.boolToOnOff(obj.isBlueMovingAvailable);
            drawnow % If we are in mid-operation
        end
        
        function edtStepSizeCallback(obj)
            stepSizeString = obj.edtSteps.String;
            try
                if ~ValidationHelper.isStringValueANumber(stepSizeString)
                    error('Step size must be a valid number! Reverting.');
                end
            catch err
                obj.edtSteps.String = StringHelper.formatNumber(obj.helm.stepSizeI);
                EventStation.anonymousError(err.message);
            end
            obj.helm.stepSizeI = str2double(stepSizeString);
        end
        
        function btnMoveCallback(obj, index, trueForLeftFalseForRight)
            step = BooleanHelper.ifTrueElse(trueForLeftFalseForRight,-1,1)*obj.helm.stepSizeI;
            obj.helm.relativeCurrentStep(index, step);
        end
        
        function btnMoveToBlueCallback(obj)
            isBeingGrayed = true;
            for i = 1:3
                obj.btnMoveLeft(i).Enable = 'off';
                obj.recolor(obj.tvCurB(i),isBeingGrayed)
                obj.btnMoveRight(i).Enable = 'off';
            end
            newValue = obj.helm.I;
            for i = 1:3
                Bi = obj.edtCurBValue(i);
                if Bi ~= inf
                    newValue(i) = Bi;
                    obj.clearEdtCurB(i);
                end
            end
            obj.helm.setI(newValue);
            for i = 1:3
                obj.btnMoveLeft(i).Enable = 'on';
                TargetEntry.styleNow(obj.tvCurB(i));
                obj.btnMoveRight(i).Enable = 'on';
            end
            % We are now finished moving. Waiting for new Input.
            obj.isBlueMovingAvailable = false;
            obj.refresh();
        end
        
        function tvCurBCallback(obj, index)
            limNeg = obj.helm.minCurrent(index);
            limPos = obj.helm.maxCurrent(index);
            if ValidationHelper.isStringValueInBorders(obj.tvCurB(index).String, limNeg, limPos)
                obj.edtCurB(index).String = obj.tvCurB(index).String;
                obj.showEdtCurB(index, true);
                obj.checkEdtCurBValue(index);
                obj.isBlueMovingAvailable = true;
                obj.refresh
            else
                EventStation.anonymousWarning('Current must be a number within [%.3f, %.3f] A', limNeg, limPos)
            end
            obj.tvCurB(index).String = StringHelper.formatNumber(obj.helm.I(index));
        end
        
        function showEdtCurB(~, ~, ~)
            % now | target, both always shown
            % (now and target are separate columns: nothing to show or hide)
        end

        function checkEdtCurBValue(obj, index)
            % Enter in a target: pending (amber) if it is within the supply's range;
            % empty = no target
            limNeg = obj.helm.minCurrent(index);
            limPos = obj.helm.maxCurrent(index);
            str = strtrim(obj.edtCurB(index).String);
            if isempty(str)
                obj.clearEdtCurB(index);
            elseif ValidationHelper.isStringValueInBorders(str, limNeg, limPos)
                obj.edtCurBValue(index) = str2double(str);
                obj.edtCurB(index).String = StringHelper.formatNumber(obj.edtCurBValue(index));
                TargetEntry.showPending(obj.edtCurB(index), true);
            else
                EventStation.anonymousWarning('Current must be a number between %.3f and %.3f A! Reverting.', limNeg, limPos)
                if isinf(obj.edtCurBValue(index))
                    obj.clearEdtCurB(index);
                else
                    obj.edtCurB(index).String = StringHelper.formatNumber(obj.edtCurBValue(index));
                end
            end
            obj.isBlueMovingAvailable = any(~isinf(obj.edtCurBValue));
            obj.refresh();
        end

        function clearAllTargets(obj)
            for i = 1:3; obj.clearEdtCurB(i); end
            obj.isBlueMovingAvailable = false;
            obj.refresh();
        end

        function clearEdtCurB(obj, index)
            % Clears the target (and its stored value)
            obj.edtCurB(index).String = '';
            obj.edtCurBValue(index) = inf;
            TargetEntry.showPending(obj.edtCurB(index), false);
        end
        
    end

    %% overridden from EventListener
    methods
        % When events happen, this function jumps.
        % event is the event sent from the EventSender
        function onEvent(obj, event)
            if event.isError ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_STEP_SIZE_CHANGED) ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_MAGNETIC_FIELD_CHANGED)
                obj.refresh();
            end
        end
    end
    
    
end