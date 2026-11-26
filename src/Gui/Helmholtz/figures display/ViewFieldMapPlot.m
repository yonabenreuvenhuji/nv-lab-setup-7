classdef ViewFieldMapPlot < GuiComponent & EventListener

    properties
        helm
        vAxes
    end

    methods
        function obj = ViewFieldMapPlot(parent, controller)
            expNames = Experiment.getExperimentNames();
            namesToListenTo = {expNames{:}, Helmholtz.NAME}; %#ok<CCAT>
            obj@GuiComponent(parent, controller);
            obj@EventListener(namesToListenTo);
            obj.helm = getObjByName('Helmholtz');
            obj.component = uicontainer('parent', parent.component);
            obj.vAxes = axes('Parent', obj.component, 'ActivePositionProperty', 'outerposition');
            obj.helm.plotFieldMap(obj.vAxes);
            obj.height = 420;
            obj.width = 450;
        end

        function refresh(obj)
            [az, el] = view(obj.vAxes);
            obj.helm.plotFieldMap(obj.vAxes);
            view(obj.vAxes, az, el);
        end
    end

    methods
        function onEvent(obj, event)
            if event.isError ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_MAGNETIC_FIELD_CHANGED) ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_CALIBRATION_CHANGED)
                obj.refresh();
            end
        end
    end
end
