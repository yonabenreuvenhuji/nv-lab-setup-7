classdef ViewFieldMap < ViewVBox
    %VIEWFIELDMAP The fields the coils can make, the NV axes, and B -- in one 3D plot.
    %
    % Replaces the diamond-placement drawing: with the coils calibrated against the NVs
    % there is no hand-entered diamond orientation left to draw. What is worth seeing is
    % the parallelepiped of reachable fields (the current box mapped through M) against
    % the NV axes, which shows at a glance which orientations can be used and how far.
    % Rotate it with the axes toolbar; the view is kept across refreshes.

    methods
        function obj = ViewFieldMap(parent, controller)
            vMainPanel = ViewExpandablePanel(parent, controller, 'Reachable fields and the NV axes');
            obj@ViewVBox(vMainPanel, controller);
            obj.component.Spacing = 5;
            obj.component.Padding = 5;
            vPlot = ViewFieldMapPlot(obj, controller);
            obj.height = vPlot.height + 10;
            obj.width = vPlot.width + 10;
        end
    end
end
