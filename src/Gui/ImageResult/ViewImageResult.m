classdef ViewImageResult < ViewVBox
    %VIEWIMAGERESULT The scan image: a toolbar, the image, and (for 3D scans) a slice bar.

    properties
        vHeader % the toolbar
        vImage  % the image view
        vSlices % slice bar (3D scans)
    end

    methods
        function obj = ViewImageResult(parent, controller, minImageSize, camera)
            if ~exist('minImageSize', 'var'); minImageSize = []; end
            if ~exist('camera', 'var'); camera = []; end
            obj@ViewVBox(parent, controller, 0, MainGuiStyle.GAP);
            obj.component.BackgroundColor = MainGuiStyle.col().window;
            obj.vHeader = ViewImageResultHeader(obj, controller, camera);
            obj.vHeader.startListeningTo(ImageScanResult.NAME);     % This will ensure that header actions will occur after image has updated
            obj.vImage = ViewImageResultImage(obj, controller, minImageSize);
            obj.vSlices = ViewImageResultSlices(obj, controller);

            obj.height = obj.vHeader.height + obj.vImage.height + 10;
            obj.width = max([obj.vHeader.width, obj.vImage.width]) + 10;
            if ~isempty(minImageSize); obj.width = 600; end

            obj.setHeights([obj.vHeader.height, -1, 0]);
            obj.vSlices.onShow = @(tf) obj.showSliceBar(tf);
            obj.vSlices.update();
        end

        function showSliceBar(obj, tf)
            h = obj.component.Heights;
            h(3) = tf * obj.vSlices.height;
            obj.component.Heights = h;
        end
    end
end
