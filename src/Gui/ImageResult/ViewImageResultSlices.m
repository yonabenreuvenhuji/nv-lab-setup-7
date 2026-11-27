classdef ViewImageResultSlices < GuiComponent & EventListener
    %VIEWIMAGERESULTSLICES The bar under the image for 3D scans: which slice is shown.
    %
    %   Slice  [◀] [=========|===========] [▶]   z = 2.5 µm (3 of 11)   Show [One slice ▾]
    %
    % Hidden unless the image is a 3D scan. While a 3D scan runs, the slice being scanned is
    % shown, until you pick one yourself. "Show" can also be the maximum or the mean over all
    % slices (a projection along the slice axis).

    properties
        btnPrev
        btnNext
        slider
        tvSlice
        popView
        onShow = []     % function handle @(tf): the parent shows / hides this bar
    end

    properties (Access = private)
        isShown = false
        isUpdating = false
    end

    methods
        function obj = ViewImageResultSlices(parent, controller)
            obj@GuiComponent(parent, controller);
            obj@EventListener(ImageScanResult.NAME);
            C = MainGuiStyle.col();

            box = uix.HBox('Parent', parent.component, 'Spacing', 6, 'Padding', 4, 'BackgroundColor', C.panel);
            obj.component = box;
            MainGuiStyle.label(box, 'Slice', 'right', 'FontWeight', 'bold');
            obj.btnPrev = MainGuiStyle.button(box, char(9664), 'normal', 'Callback', @(h,e) obj.step(-1), ...
                'TooltipString', 'Previous slice');
            obj.slider = MainGuiStyle.slider(box, 'Min', 1, 'Max', 2, 'Value', 1, 'SliderStep', [1 1]);
            addlistener(obj.slider, 'ContinuousValueChange', @(h,e) obj.sliderMoved());
            obj.btnNext = MainGuiStyle.button(box, char(9654), 'normal', 'Callback', @(h,e) obj.step(+1), ...
                'TooltipString', 'Next slice');
            obj.tvSlice = MainGuiStyle.label(box, '', 'left', 'FontSize', MainGuiStyle.FONT + 1);
            MainGuiStyle.label(box, 'Show', 'right');
            obj.popView = MainGuiStyle.popup(box, ImageScanResult.SLICE_VIEW_OPTIONS, ...
                'Callback', @(h,e) obj.viewChanged(), 'TooltipString', ...
                'One slice, or the maximum / mean over all slices (a projection along the slice axis)');
            box.Widths = [36 28 -1 28 190 36 150];
            obj.height = 34;
            obj.width = 600;
        end

        function update(obj)
            isr = getObjByName(ImageScanResult.NAME);
            is3D = ~isempty(isr) && isr.mDimNumber == 3 && numel(isr.mThirdAxis) >= 1;
            if is3D ~= obj.isShown && ~isempty(obj.onShow)
                obj.onShow(is3D);
            end
            obj.isShown = is3D;
            if ~is3D; return; end

            obj.isUpdating = true;
            n = numel(isr.mThirdAxis);
            k = isr.sliceIndex;
            if n > 1
                set(obj.slider, 'Min', 1, 'Max', n, 'Value', k, 'SliderStep', [1 max(1, round(n/5))] / (n - 1), 'Enable', 'on');
            else
                set(obj.slider, 'Min', 1, 'Max', 2, 'Value', 1, 'Enable', 'off');
            end
            single = isr.sliceView == 1;
            set([obj.btnPrev obj.btnNext], 'Enable', BooleanHelper.boolToOnOff(single && n > 1));
            if ~single; obj.slider.Enable = 'off'; end
            obj.popView.Value = isr.sliceView;
            follow = '';
            scanner = getObjByName(StageScanner.NAME);
            if isr.sliceFollow && ~isempty(scanner) && scanner.mCurrentlyScanning
                follow = '  · following the scan';
            end
            obj.tvSlice.String = sprintf('%s = %s %s  (%d of %d)%s', isr.mAxesString(3), ...
                MainGuiStyle.num(isr.mThirdAxis(k), 3), MainGuiStyle.micron, k, n, follow);
            obj.isUpdating = false;
        end

        function step(obj, d)
            isr = getObjByName(ImageScanResult.NAME);
            isr.showSlice(isr.sliceIndex + d);
        end

        function sliderMoved(obj)
            if obj.isUpdating; return; end
            isr = getObjByName(ImageScanResult.NAME);
            k = round(obj.slider.Value);
            if k ~= isr.sliceIndex
                isr.showSlice(k);
            end
        end

        function viewChanged(obj)
            isr = getObjByName(ImageScanResult.NAME);
            isr.setSliceView(obj.popView.Value);
        end
    end

    methods
        function onEvent(obj, event)
            if isfield(event.extraInfo, ImageScanResult.EVENT_IMAGE_UPDATED)
                obj.update();
            end
        end
    end
end
