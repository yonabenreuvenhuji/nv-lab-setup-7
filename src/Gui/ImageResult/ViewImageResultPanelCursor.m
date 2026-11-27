classdef ViewImageResultPanelCursor < GuiComponent
    %VIEWIMAGERESULTPANELCURSOR "Mouse [Inspect][Select scan area][Move stage here]" in the
    %image toolbar: what a click on the image does.
    %
    %   Inspect           a data tip with the position and the value (was "Marker")
    %   Select scan area  drag a rectangle: it becomes the Scan area's From/To (was "Zoom");
    %                     back to Inspect when done
    %   Move stage here   a click moves the stage there (was "Move to"); in a 3D scan's
    %                     slice, also to that slice

    properties
        btnModes     % 1x3 toggle buttons, in the order of ImageScanResult.CURSOR_OPTIONS
    end

    properties (Constant)
        LABELS = {'Inspect', 'Select scan area', 'Move stage here'}
        TIPS = {'Click a pixel to see where it is and its value', ...
            'Drag a rectangle on the image: it becomes the From/To of the Scan area', ...
            'Click on the image to move the stage there'}
    end

    methods
        function obj = ViewImageResultPanelCursor(parent, controller)
            obj@GuiComponent(parent, controller);
            box = MainGuiStyle.hbox(parent.component, 2);
            obj.component = box;
            MainGuiStyle.label(box, 'Mouse', 'right');
            obj.btnModes = gobjects(1, 3);
            for k = 1:3
                obj.btnModes(k) = MainGuiStyle.toggle(box, obj.LABELS{k}, ...
                    'TooltipString', obj.TIPS{k}, 'Callback', @(h,e) obj.select(k));
            end
            box.Widths = [44 62 110 110];
            obj.height = MainGuiStyle.ROW_H;
            obj.width = sum(box.Widths) + 6;
            obj.show(1);
        end

        function update(obj)
            % Executes when image updates. A full redraw puts the image back in Inspect mode;
            % a scan's line-by-line (or pixel-by-pixel) updates do not.
            imageScanResult = getObjByName(ImageScanResult.NAME);
            if isempty(imageScanResult) || imageScanResult.lastUpdateFull
                obj.show(1);
            end
        end

        function select(obj, k)
            imageScanResult = getObjByName(ImageScanResult.NAME);
            if isempty(imageScanResult); throwBaseObjException(ImageScanResult.NAME); end
            obj.show(k);
            imageScanResult.updateDataCursor(ImageScanResult.CURSOR_OPTIONS{k});
            if k == 2
                % The rectangle has been drawn (or cancelled): back to Inspect
                obj.show(1);
                imageScanResult.updateDataCursor(ImageScanResult.CURSOR_OPTIONS{1});
            end
        end
    end

    methods (Access = private)
        function show(obj, k)
            C = MainGuiStyle.col();
            for j = 1:3
                if j == k
                    set(obj.btnModes(j), 'Value', 1, 'BackgroundColor', C.accent, 'ForegroundColor', [1 1 1]);
                else
                    set(obj.btnModes(j), 'Value', 0, 'BackgroundColor', [0.94 0.94 0.94], 'ForegroundColor', C.ink);
                end
            end
        end
    end
end
