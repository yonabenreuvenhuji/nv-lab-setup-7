classdef ViewImageResultPanelContrast < GuiComponent
    %VIEWIMAGERESULTPANELCONTRAST "MW contrast [Difference ▾]" in the image toolbar: what a
    %MW contrast scan shows. Greyed out for ordinary scans.

    properties
        popup
    end

    properties (Constant)
        ITEMS = {'Difference (without - with)', 'Ratio (with / without)', 'Without MW only', 'With MW only'}
    end

    methods
        function obj = ViewImageResultPanelContrast(parent, controller)
            obj@GuiComponent(parent, controller);
            box = MainGuiStyle.hbox(parent.component, 4);
            obj.component = box;
            MainGuiStyle.label(box, 'MW contrast', 'right');
            obj.popup = MainGuiStyle.popup(box, obj.ITEMS, 'Callback', @obj.popupCallback, ...
                'TooltipString', 'For MW contrast scans: what the image shows');
            box.Widths = [74 170];
            obj.height = MainGuiStyle.ROW_H;
            obj.width = 250;
        end

        function update(obj)
            imageScanResult = getObjByName(ImageScanResult.NAME);
            d = imageScanResult.mDimNumber;
            isContrast = size(imageScanResult.mData, d + 1 + (d == 3)) == 2;
            obj.popup.Value = imageScanResult.contrastType;
            obj.popup.Enable = BooleanHelper.boolToOnOff(isContrast);
        end

        function popupCallback(obj, ~, ~)
            imageScanResult = getObjByName(ImageScanResult.NAME);
            imageScanResult.contrastType = obj.popup.Value;
            imageScanResult.update;
        end
    end
end
