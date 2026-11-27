classdef ViewMainImage < ViewVBox
    %VIEWMAINIMAGE The main setup window.
    %
    %   +-- stage(s) ----------+---- scan image ---------------+-- laser ------------+
    %   | Position             |  view · colour · contrast ·   |  photodiode:        |
    %   | Scan area            |  cursor toolbar               |   live reading      |
    %   | Scan  [progress]     |                               |   live trace        |
    %   | MW contrast          |        image                  |   controls          |
    %   | Tracking | Stage     |                               |  save & load        |
    %   | Soft limits          |  (3D: slice slider)           |                     |
    %   +----------------------+-------------------------------+---------------------+
    %   status bar: the last error or warning, and where it came from
    %
    % The window opens maximized (GuiControllerImage). The stage column scrolls if the
    % screen is not tall enough for it.

    properties (Constant)
        LEFT_WIDTH = 640    % stage column, px
        RIGHT_WIDTH = 450   % detector column, px
    end

    methods
        % constructor - parent of a main view is the figure itself!
        function obj = ViewMainImage(parent, controller)
            obj@ViewVBox(parent, controller, 0, 0);
            C = MainGuiStyle.col();
            fig = ancestor(obj.component, 'figure');
            set(fig, 'Color', C.window);
            obj.component.BackgroundColor = C.window;

            columns = ViewBox(obj, controller, @uix.HBox, 'Spacing', 10, 'Padding', 8, ...
                'BackgroundColor', C.window);

            % Left: the stage(s)
            left = ViewBox(columns, controller, @uix.ScrollingPanel, 'BackgroundColor', C.window);
            stages = ViewContainerStage(left, controller);
            set(left.component, 'Height', -1, 'MinimumHeight', stages.height, ...
                'Width', -1, 'MinimumWidth', obj.LEFT_WIDTH - 20);

            % Centre: the scan image
            ViewImageResult(columns, controller);

            % Right: laser, detector, save & load
            right = ViewVBox(columns, controller, 0, MainGuiStyle.GAP);
            right.component.BackgroundColor = C.window;
            viewLaserContainer = ViewLasersContainer(right, controller);
            ViewSpcm(right, controller);
            viewSaveLoad = ViewSaveLoad(right, controller, Savable.CATEGORY_IMAGE);
            right.setHeights([viewLaserContainer.height, -1, viewSaveLoad.height]);

            columns.component.Widths = [obj.LEFT_WIDTH, -1, obj.RIGHT_WIDTH];

            status = ViewStatusBar(obj, controller);
            obj.component.Heights = [-1, status.height];

            obj.width = obj.LEFT_WIDTH + obj.RIGHT_WIDTH + 600;
            obj.height = 860;
        end
    end
end
