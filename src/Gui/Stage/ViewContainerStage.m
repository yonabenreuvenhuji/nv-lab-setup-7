classdef ViewContainerStage < ViewVBox
    %VIEWCONTAINERSTAGE One ViewStage per stage of the setup, stacked.

    methods
        function obj = ViewContainerStage(parent, controller, camera)
            if ~exist('camera', 'var'); camera = []; end
            obj@ViewVBox(parent, controller, 0, 3*MainGuiStyle.GAP);
            obj.component.BackgroundColor = MainGuiStyle.col().window;
            stages = ClassStage.getStages;

            heights = zeros(1, length(stages));
            width = 0;
            for k = 1 : length(stages)
                stageView = ViewStage(obj, controller, stages{k}, camera);
                heights(k) = stageView.height;
                width = max(width, stageView.width);
            end

            obj.component.Heights = heights;
            obj.height = sum(heights) + obj.component.Spacing * max(numel(heights) - 1, 0);
            obj.width = width;
        end
    end
end
