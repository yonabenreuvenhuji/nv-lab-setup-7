classdef ViewCoilCalibration < GuiComponent & EventListener
    %VIEWCOILCALIBRATION What the field is computed from: B = M*I + b0, and where M came from.

    properties
        helm
        txtInfo
        btnClear
    end

    methods
        function obj = ViewCoilCalibration(parent, controller)
            expNames = Experiment.getExperimentNames();
            namesToListenTo = {expNames{:}, Helmholtz.NAME}; %#ok<CCAT>
            obj@GuiComponent(parent, controller);
            obj@EventListener(namesToListenTo);
            obj.helm = getObjByName('Helmholtz');

            panelMain = uix.Panel('Parent', parent.component, 'Title', ...
                'Coil model: B = M*I + b0', 'Padding', 5);
            vbox = uix.VBox('Parent', panelMain, 'Spacing', 4);
            obj.txtInfo = uicontrol(obj.PROP_TEXT_NORMAL{:}, 'Parent', vbox, ...
                'HorizontalAlignment', 'left', 'FontName', 'FixedWidth', 'Max', 2);
            obj.btnClear = uicontrol(obj.PROP_BUTTON_SMALL{:}, 'Parent', vbox, ...
                'String', 'Use the nominal model', 'TooltipString', ...
                ['Drop the coil calibration and go back to the ratios in setupInfo.json ' ...
                 '(the currents do not change). Load a json to get a calibration back.']);
            vbox.Heights = [-1, 24];
            obj.btnClear.Callback = @(h,e) obj.helm.clearCoilCalibration();
            obj.height = 175;
            obj.width = 400;
            obj.refresh();
        end

        function refresh(obj)
            M = obj.helm.coilMatrix; b = obj.helm.coilOffset;
            lines = {obj.helm.calibrationSummary(), ''};
            lines{end+1} = sprintf('M (G/A)    coil x    coil y    coil z');
            ax = 'xyz';
            for i = 1:3
                lines{end+1} = sprintf('  B%s     %8.4f  %8.4f  %8.4f', ax(i), M(i, :)); %#ok<AGROW>
            end
            lines{end+1} = sprintf('b0 (G)   %+8.4f  %+8.4f  %+8.4f', b);
            lines{end+1} = sprintf('I (A)     %8.4f  %8.4f  %8.4f', obj.helm.I);
            w = obj.helm.extrapolationWarning();
            if ~isempty(w); lines{end+1} = w; end
            obj.txtInfo.String = lines;
            if obj.helm.isCalibrated
                obj.txtInfo.ForegroundColor = [0 0 0];
                obj.btnClear.Enable = 'on';
            else
                obj.txtInfo.ForegroundColor = [0.8 0 0];
                obj.btnClear.Enable = 'off';
            end
        end
    end

    methods
        function onEvent(obj, event)
            if event.isError ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_MAGNETIC_FIELD_CHANGED) ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_CALIBRATION_CHANGED)
                obj.refresh();
            end
        end
    end
end
