classdef ViewESRPlot < GuiComponent & EventListener
    %VIEWESRPLOT The expected ESR at the present field, and at the field being considered.
    %
    % Drawn by Helmholtz.plotExpectedESR: a smooth spectrum from the exact Hamiltonian, every
    % dip marked and labelled (frequency, and x3 where three orientations share it), the
    % preview from the alignment panel dashed in purple, D, and the operating frequency from
    % the json. The axes have their own toolbar (zoom, pan, data tips, restore, export) --
    % the window has no figure toolbar -- and a zoom survives the redraws that every field
    % change causes, until "restore view".

    properties
        helm
        vAxes
        autoXLim = []       % the limits the last redraw chose; anything else is a user zoom
        operatingMHz = NaN
        operatingFile = ''
    end

    methods
        function obj = ViewESRPlot(parent, controller)
            expNames = Experiment.getExperimentNames();
            namesToListenTo = {expNames{:}, Helmholtz.NAME}; %#ok<CCAT>
            obj@GuiComponent(parent, controller);
            obj@EventListener(namesToListenTo);
            obj.helm = getObjByName('Helmholtz');

            obj.component = uicontainer('parent', parent.component, ...
                'BackgroundColor', Helmholtz.plotColours().surface);
            obj.vAxes = axes('Parent', obj.component, 'ActivePositionProperty', 'outerposition');
            try
                axtoolbar(obj.vAxes, {'export', 'datacursor', 'pan', 'zoomin', 'zoomout', ...
                                      'restoreview'});
                enableDefaultInteractivity(obj.vAxes);
            catch
            end
            obj.refresh();
            obj.height = 470;   % minimum
            obj.width = 620;    % minimum
        end

        function refresh(obj)
            ax = obj.vAxes;
            curX = ax.XLim; curY = ax.YLim;
            userZoom = ~isempty(obj.autoXLim) && any(abs(curX - obj.autoXLim) > 1e-9);
            obj.readOperatingFrequency();
            obj.autoXLim = obj.helm.plotExpectedESR(ax, struct('operatingMHz', obj.operatingMHz));
            % the auto view is what "restore view" goes back to
            try; resetplotview(ax, 'SaveCurrentView'); catch; end
            if userZoom
                ax.XLim = curX; ax.YLim = curY;
            end
        end

        function readOperatingFrequency(obj)
            % The json's top-level "frequency" -- the MW the experiments run at. Re-read when
            % the state file changes, which is when it can have changed.
            f = obj.helm.stateFile;
            if strcmp(f, obj.operatingFile); return; end
            obj.operatingFile = f;
            obj.operatingMHz = NaN;
            try
                J = jsondecode(fileread(f));
                if isfield(J, 'frequency') && isnumeric(J.frequency)
                    obj.operatingMHz = J.frequency(1);
                end
            catch
            end
        end
    end

    %% overridden from EventListener
    methods
        function onEvent(obj, event)
            if event.isError ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_MAGNETIC_FIELD_CHANGED) ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_CALIBRATION_CHANGED) ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_PREVIEW_CHANGED)
                obj.refresh();
            end
        end
    end
end
