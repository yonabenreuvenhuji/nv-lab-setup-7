classdef ViewImageResultPanelPlot < GuiComponent
    %VIEWIMAGERESULTPANELPLOT "Aspect [Equal ▾]" in the image toolbar (and the
    %"Open in figure" action).

    properties
        popupStyle
    end

    methods
        function obj = ViewImageResultPanelPlot(parent, controller)
            obj@GuiComponent(parent, controller);
            box = MainGuiStyle.hbox(parent.component, 4);
            obj.component = box;
            MainGuiStyle.label(box, 'Aspect', 'right');
            obj.popupStyle = MainGuiStyle.popup(box, ImageScanResult.PLOT_STYLE_OPTIONS_2D, ...
                'Callback', @obj.popupStyleCallback, 'TooltipString', ...
                'Normal: fill the space. Equal: 1 µm is as long on both axes. Square: a square plot.');
            box.Widths = [46 80];
            obj.height = MainGuiStyle.ROW_H;
            obj.width = 130;
            obj.update;
        end

        function update(obj)
            imageScanResult = getObjByName(ImageScanResult.NAME);
            if isempty(imageScanResult); throwBaseObjException(ImageScanResult.NAME); end
            switch imageScanResult.mDimNumber
                case 1
                    obj.popupStyle.String = ImageScanResult.PLOT_STYLE_OPTIONS_1D;
                otherwise
                    obj.popupStyle.String = ImageScanResult.PLOT_STYLE_OPTIONS_2D;
            end
            obj.popupStyle.Value = min(imageScanResult.plotStyle, numel(obj.popupStyle.String));
        end

        function btnOpenInFigureCallback(obj, ~, ~) %#ok<INUSD>
            imageScanResult = getObjByName(ImageScanResult.NAME);
            if isempty(imageScanResult); throwBaseObjException(ImageScanResult.NAME); end
            isVisible = true;
            imageScanResult.copyToFigure(isVisible);
        end
    end

    methods (Access = private)
        function popupStyleCallback(obj, ~, ~)
            imageScanResult = getObjByName(ImageScanResult.NAME);
            if isempty(imageScanResult); throwBaseObjException(ImageScanResult.NAME); end
            imageScanResult.plotStyle = obj.popupStyle.Value;
            imageScanResult.imagePostProcessing;    % Update added layer (including plot style)
        end
    end
end
