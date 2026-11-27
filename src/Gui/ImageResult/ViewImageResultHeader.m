classdef ViewImageResultHeader < GuiComponent & EventListener
    %VIEWIMAGERESULTHEADER The toolbar above the scan image, two rows:
    %
    %   Aspect [Equal ▾]   Colormap [Parula ▾] [x] Auto  Min [ ]  Max [ ]   [Open in figure]
    %   Mouse [Inspect][Select scan area][Move stage here]   MW contrast [Difference ▾]  (Image [V ▾])

    properties
        vPlotOptions
        vColorMap
        vCursor
        vBinning
        vContrast
        vCameraScan
        views   % all of them.
    end

    methods
        function obj = ViewImageResultHeader(parent, controller, camera)
            try
                spcm = getObjByName(Spcm.NAME);
                hasAdvanced = spcm.hasG2 || spcm.hasLifetime;
            catch
                hasAdvanced = 0;
            end
            obj@GuiComponent(parent, controller);
            obj@EventListener;
            S = MainGuiStyle;

            [panel, inner] = MainGuiStyle.section(parent.component, 'Scan image');
            obj.component = panel;
            innerView = struct('component', inner);     % (a parent for the row views)
            row1 = ViewBox(innerView, controller, @uix.HBox, 'Spacing', 6, 'BackgroundColor', MainGuiStyle.col().panel);
            row2 = ViewBox(innerView, controller, @uix.HBox, 'Spacing', 6, 'BackgroundColor', MainGuiStyle.col().panel);

            % Row 1: aspect, colours, open in figure
            obj.vPlotOptions = ViewImageResultPanelPlot(row1, controller);
            obj.vColorMap = ViewImageResultPanelColormap(row1, controller);
            MainGuiStyle.empty(row1.component);
            figureButton = MainGuiStyle.button(row1.component, 'Open in figure', 'normal', ...
                'TooltipString', 'Copy the image into its own MATLAB figure', ...
                'Callback', @(h,e) obj.vPlotOptions.btnOpenInFigureCallback());
            row1.component.Widths = [obj.vPlotOptions.width, obj.vColorMap.width, -1, 110];

            % Row 2: mouse mode, contrast, (lifetime image), (camera z scan)
            obj.vCursor = ViewImageResultPanelCursor(row2, controller);
            obj.vContrast = ViewImageResultPanelContrast(row2, controller);
            widths = [obj.vCursor.width, obj.vContrast.width];
            obj.views = {obj.vPlotOptions, obj.vColorMap, obj.vCursor, obj.vContrast};
            if hasAdvanced
                obj.vBinning = ViewImageResultPanelBinning(row2, controller);
                widths(end + 1) = obj.vBinning.width;
                obj.views{end + 1} = obj.vBinning;
            end
            if ~isempty(camera)
                obj.vCameraScan = ViewCameraZScan(row2, controller);
                widths(end + 1) = obj.vCameraScan.width;
                obj.views{end + 1} = obj.vCameraScan;
            end
            MainGuiStyle.empty(row2.component);
            row2.component.Widths = [widths, -1];
            set(figureButton, 'Tag', 'openInFigure');

            inner.Heights = [S.ROW_H S.ROW_H];
            obj.height = MainGuiStyle.sectionHeight(MainGuiStyle.rowsHeight(inner.Heights));
            obj.width = 700;
        end

        function updateAxes(obj)
            cellfun(@(v) v.update, obj.views);
        end
    end

    %% overridden from EventListener
    methods
        function onEvent(obj, event)
            if strcmp(event.creator.name, ImageScanResult.NAME)
                if isfield(event.extraInfo, ImageScanResult.EVENT_IMAGE_UPDATED) || ...
                        event.isError
                    obj.updateAxes;
                end
            end
        end
    end
end
