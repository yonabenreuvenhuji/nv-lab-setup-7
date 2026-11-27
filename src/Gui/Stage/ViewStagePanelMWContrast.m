classdef ViewStagePanelMWContrast < GuiComponent & EventListener & EventSender
    %VIEWSTAGEPANELMWCONTRAST The "MW contrast" section: scan every line twice, without and
    %with microwaves, to image the ODMR contrast at one frequency.

    properties
        cbxMWContrast       % checkbox
        edtFrequency        % edit text, MHz
        edtAmplitude        % edit text, dBm

        stageName           % string
    end
    methods
        function obj = ViewStagePanelMWContrast(parent, controller, stage)
            obj@GuiComponent(parent, controller);
            obj@EventListener(stage.name);
            obj@EventSender(sprintf('%s%s', stage.name, ' _ panel Scan Contrast'));
            obj.stageName = stage.name;
            S = MainGuiStyle;

            [panel, inner] = MainGuiStyle.section(parent.component, 'MW contrast');
            obj.component = panel;
            row = MainGuiStyle.hbox(inner, 6);
            obj.cbxMWContrast = MainGuiStyle.checkbox(row, 'Scan with and without MW', false, 'TooltipString', ...
                'Scan every line twice, without and with microwaves; the image shows the contrast');
            MainGuiStyle.label(row, 'Frequency (MHz)', 'right');
            obj.edtFrequency = MainGuiStyle.edit(row, '');
            MainGuiStyle.label(row, 'Power (dBm)', 'right');
            obj.edtAmplitude = MainGuiStyle.edit(row, '');
            row.Widths = [-1 100 70 76 56];
            inner.Heights = S.ROW_H;

            %%%% callbacks %%%%
            obj.cbxMWContrast.Callback = @(h,e) obj.cbxMWContrastCallback;
            obj.edtFrequency.Callback = @(h,e) obj.edtFrequencyCallback;
            obj.edtAmplitude.Callback = @(h,e) obj.edtAmplitudeCallback;

            obj.height = MainGuiStyle.sectionHeight(S.ROW_H);
            obj.width = 520;
            obj.refresh();  % init values
        end

        function refresh(obj)
            stage = getObjByName(obj.stageName);
            scanParams = stage.scanParams;
            obj.cbxMWContrast.Value = scanParams.isMWContrastScan;
            obj.edtFrequency.String = MainGuiStyle.num(scanParams.MWFrequency, 4);
            obj.edtAmplitude.String = MainGuiStyle.num(scanParams.MWAmplitude, 2);
            MainGuiStyle.setFieldEnabled([obj.edtFrequency, obj.edtAmplitude], logical(scanParams.isMWContrastScan));
        end

        function cbxMWContrastCallback(obj)
            stage = getObjByName(obj.stageName);
            stage.scanParams.isMWContrastScan = logical(obj.cbxMWContrast.Value);
            stage.sendEventScanParamsChanged();
        end

        function edtFrequencyCallback(obj)
            stage = getObjByName(obj.stageName);
            if ~ValidationHelper.isValuePositive(obj.edtFrequency.String)
                EventStation.anonymousWarning('Frequency has to be a positive number (MHz)! Reverting.');
            else
                stage.scanParams.MWFrequency = str2double(obj.edtFrequency.String);
            end
            stage.sendEventScanParamsChanged();
        end

        function edtAmplitudeCallback(obj)
            stage = getObjByName(obj.stageName);
            if ~ValidationHelper.isStringValueANumber(obj.edtAmplitude.String)
                EventStation.anonymousWarning('Power has to be a number (dBm)! Reverting.');
            else
                stage.scanParams.MWAmplitude = str2double(obj.edtAmplitude.String);
            end
            stage.sendEventScanParamsChanged();
        end
    end

    %% overridden from EventListener
    methods
        function onEvent(obj, event)
            if event.isError || isfield(event.extraInfo, ClassStage.EVENT_SCAN_PARAMS_CHANGED)
                obj.refresh()
            end
        end
    end
end
