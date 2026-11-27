classdef ViewImageResultPanelBinning < GuiComponent
    %VIEWIMAGERESULTPANELBINNING "Image [kcps ▾ | Lifetime]" in the image toolbar, for
    %detectors that record a time histogram.

    properties
        popup
    end

    methods
        function obj = ViewImageResultPanelBinning(parent, controller)
            obj@GuiComponent(parent, controller);
            box = MainGuiStyle.hbox(parent.component, 4);
            obj.component = box;
            MainGuiStyle.label(box, 'Image', 'right');
            obj.popup = MainGuiStyle.popup(box, {Spcm.signalUnits(), 'Lifetime'}, ...
                'Callback', @obj.popupCallback);
            box.Widths = [40 90];
            obj.height = MainGuiStyle.ROW_H;
            obj.width = 136;
        end

        function update(obj)
            imageScanResult = getObjByName(ImageScanResult.NAME);
            hasHist = ~isempty(imageScanResult.mHist);
            if ~hasHist
                obj.popup.String = {Spcm.signalUnits()};
            else
                obj.popup.String = {Spcm.signalUnits(), 'Lifetime'};
            end
            obj.popup.Value = min(imageScanResult.imageType, numel(obj.popup.String));
        end

        function popupCallback(obj, ~, ~)
            imageScanResult = getObjByName(ImageScanResult.NAME);
            imageScanResult.imageType = obj.popup.Value;
            imageScanResult.update;
        end
    end
end
