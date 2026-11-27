classdef ViewImageResultPanelColormap < GuiComponent
    %VIEWIMAGERESULTPANELCOLORMAP "Colormap [Parula ▾] [x] Auto  Min [ ]  Max [ ]" in the
    %image toolbar. Typing a limit switches Auto off.

    properties
        popupTypes  % popup/dropdown-menu. Chooses names of the available colormaps
        cbxAuto     % checkbox. Wheteher the colormap is auto-scaled
        edtMin      % edit-input. Minimum value of colormap
        edtMax      % edit-input. Maximum value of colormap
    end

    methods
        function obj = ViewImageResultPanelColormap(parent, controller)
            obj@GuiComponent(parent, controller);
            box = MainGuiStyle.hbox(parent.component, 4);
            obj.component = box;
            MainGuiStyle.label(box, 'Colormap', 'right');
            obj.popupTypes = MainGuiStyle.popup(box, ImageScanResult.COLORMAP_OPTIONS, ...
                'Callback', @obj.popupTypesCallback);
            obj.cbxAuto = MainGuiStyle.checkbox(box, 'Auto', true, 'Callback', @obj.cbxAutoCallback, ...
                'TooltipString', 'Colour range from the data (all slices, in 3D)');
            MainGuiStyle.label(box, 'Min', 'right');
            obj.edtMin = MainGuiStyle.edit(box, '0', 'Callback', @obj.edtMinCallback);
            MainGuiStyle.label(box, 'Max', 'right');
            obj.edtMax = MainGuiStyle.edit(box, '1', 'Callback', @obj.edtMaxCallback);
            box.Widths = [60 80 50 26 64 28 64];
            obj.height = MainGuiStyle.ROW_H;
            obj.width = sum(box.Widths) + 4*6;

            try
                obj.update;     % It might not succeed if there is no image
            catch
            end
        end

        function update(obj)
            imageScanResult = obj.getImageScanResult();
            obj.popupTypes.Value = imageScanResult.colormapType;
            obj.cbxAuto.Value = imageScanResult.colormapAuto;
            obj.edtMin.String = obj.fmt(imageScanResult.colormapLimits(1));
            obj.edtMax.String = obj.fmt(imageScanResult.colormapLimits(2));
            MainGuiStyle.setPending([obj.edtMin obj.edtMax], false);
            if imageScanResult.colormapAuto
                set([obj.edtMin obj.edtMax], 'ForegroundColor', MainGuiStyle.col().ink3);
            end
        end

        function popupTypesCallback(obj, ~, ~)
            imageScanResult = obj.getImageScanResult();
            imageScanResult.colormapType = obj.popupTypes.Value;
            imageScanResult.imagePostProcessing;    % which now updates added layer (including colormap)
        end

        function cbxAutoCallback(obj, ~, ~)
            imageScanResult = obj.getImageScanResult();
            imageScanResult.colormapAuto = logical(obj.cbxAuto.Value);
            if imageScanResult.isDataAvailable
                imageScanResult.imagePostProcessing;
            end
            obj.update;
        end

        function edtMinCallback(obj, ~, ~)
            obj.setLimit(1, obj.edtMin);
        end

        function edtMaxCallback(obj, ~, ~)
            obj.setLimit(2, obj.edtMax);
        end
    end

    methods (Access = private)
        function setLimit(obj, which, field)
            imageScanResult = obj.getImageScanResult();
            limits = imageScanResult.colormapLimits;     % = [minVal maxVal]
            value = str2double(field.String);
            if ~ValidationHelper.isStringValueANumber(field.String)
                EventStation.anonymousWarning('Colour limits must be numbers! Reverting.')
            elseif (which == 1 && value > limits(2)) || (which == 2 && value < limits(1))
                EventStation.anonymousWarning('The colour minimum must be below the maximum! Reverting.')
            else
                imageScanResult.colormapAuto = false;
                limits(which) = value;
                imageScanResult.colormapLimits = limits;
                imageScanResult.imagePostProcessing;    % Updates added layer (including colormap)
            end
            obj.update;
        end
    end

    methods (Static)
        function isr = getImageScanResult
            isr = getObjByName(ImageScanResult.NAME);
            if isempty(isr)
                throwBaseObjException(ImageScanResult.NAME);
            end
        end

        function s = fmt(v)
            s = sprintf('%.4g', v);
        end
    end
end
