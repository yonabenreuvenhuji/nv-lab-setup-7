classdef ViewLasersContainer < ViewVBox
    % ViewLasersContainer One ViewLaser section per laser of the setup, stacked.

    methods
        function obj = ViewLasersContainer(parent, controller)
            obj@ViewVBox(parent, controller, 0, MainGuiStyle.GAP);
            obj.component.BackgroundColor = MainGuiStyle.col().window;

            lasers = LaserGate.getLasers();
            heights = [];
            width = 0;
            for i = 1 : length(lasers)
                viewLaser = ViewLaser(obj, controller, lasers{i});
                width = max(width, viewLaser.width);
                heights(i) = viewLaser.height; %#ok<AGROW>
            end
            if ~isempty(heights)
                set(obj.component, 'Heights', heights);
            end
            obj.height = MainGuiStyle.rowsHeight(heights);
            obj.width = width;
        end
    end
end
