classdef ViewExpectedESR < ViewVBox
    methods
        function obj = ViewExpectedESR(parent, controller)
            vMainPanel = ViewExpandablePanel(parent, controller, ...
                'Expected ESR  (zoom, pan and data tips: the toolbar at the top right of the plot)');
            obj@ViewVBox(vMainPanel, controller);
            obj.component.Spacing = 5;
            obj.component.Padding = 5;

            vImage = ViewESRPlot(obj, controller);
            obj.height = vImage.height + 10;
            obj.width = vImage.width + 10;
        end
    end
end
