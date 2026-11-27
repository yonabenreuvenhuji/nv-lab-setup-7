classdef ViewBox < GuiComponent
    %VIEWBOX Any GUI Layout Toolbox container as a GuiComponent, so views can be put in it.
    %
    %   col = ViewBox(parent, controller, @uix.ScrollingPanel, 'BackgroundColor', c);
    %   ViewContainerStage(col, controller);     % the stage views go inside it
    %
    % (ViewVBox / ViewHBox do the same for one container type each.)

    methods
        function obj = ViewBox(parent, controller, constructor, varargin)
            obj@GuiComponent(parent, controller);
            obj.component = constructor('Parent', parent.component, varargin{:});
        end
    end
end
