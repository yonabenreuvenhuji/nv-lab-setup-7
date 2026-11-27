classdef ViewImageResultImage < GuiComponent
    %VIEWIMAGERESULTIMAGE view that shows the scan results
    %   it is being used by other GUI components (such as the various
    %   options above it) as well as the StageScanner when it needs to
    %   duplicate the axes() object (which is obj.vAxes)

    properties
        vAxes       % the axes view to use for the plotting
    end
    properties (Constant)
        NAME = 'ViewImageResultImage';
    end

    methods
        function obj = ViewImageResultImage(parent, controller, minImageSize)
            obj@GuiComponent(parent, controller);
            imageScanResult = getObjByName(ImageScanResult.NAME);
            if isempty(imageScanResult)
                throwBaseObjException(ImageScanResult.NAME)
            end
            C = MainGuiStyle.col();

            obj.component = uipanel('Parent', parent.component, 'BorderType', 'line', ...
                'BackgroundColor', C.surface);
            MainGuiStyle.borderColor(obj.component, C.line);
            obj.vAxes = axes('Parent', obj.component, 'ActivePositionProperty', 'outerposition', ...
                'Color', C.surface, 'XColor', C.ink2, 'YColor', C.ink2, 'FontSize', 9, ...
                'TickDir', 'out', 'Box', 'off');
            imageScanResult.addGraphicAxes(obj.vAxes);
            MainGuiStyle.makeInteractive(obj.vAxes);

            % Set colormap and colorbar
            cMap = imageScanResult.COLORMAP_OPTIONS{imageScanResult.colormapType};
            colormap(obj.vAxes, cMap);

            % Creating floating axes() so that default calls to axes (such
            % as image() surf() etc.) won't reach this view but rather the
            % invisible floating one
            axes('Visible', 'off', 'HandleVisibility', 'on', 'Position', [0 0 eps eps]);
            imageScanResult.update;

            if isempty(minImageSize)
                obj.height = 600;   % minimum
                obj.width = 600;    % minimum
            else
                obj.height = minImageSize(1);   % minimum
                obj.width = minImageSize(2);    % minimum
            end
        end
    end
end
