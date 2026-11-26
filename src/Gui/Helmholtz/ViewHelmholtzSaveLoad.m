classdef ViewHelmholtzSaveLoad < GuiComponent & EventListener
    %VIEWHELMHOLTZSAVELOAD Save and load the Helmholtz state: the "Helmholtz" block of
    % expProperties.json (currents, field, steps and the coil calibration).
    %
    %   Save     to the json the state was last loaded from or saved to -- normally
    %            setup7Paths().json, the newest calibration folder
    %   Save As  choose a json (created if it does not exist; other keys are left alone)
    %   Load     a json with a "Helmholtz" block, or a legacy helmholtzState.mat

    properties
        btnSave
        btnSaveAs
        btnLoad
        tvFileName
        helm
    end

    methods
        function obj = ViewHelmholtzSaveLoad(parent, controller)
            panel = ViewExpandablePanel(parent, controller, 'Save & Load (expProperties.json)');
            obj@GuiComponent(panel, controller);
            obj@EventListener({Helmholtz.NAME});
            obj.helm = getObjByName('Helmholtz');

            vboxSave = uix.VBox('Parent', panel.component, 'Spacing', 5, 'Padding', 5);
            saveLoadRow = uix.HBox('Parent', vboxSave, 'Spacing', 5);
            obj.btnSave = uicontrol(obj.PROP_BUTTON{:}, 'Parent', saveLoadRow, 'string', 'Save', ...
                'TooltipString', 'write the state into the "Helmholtz" block of the json below');
            obj.btnSaveAs = uicontrol(obj.PROP_BUTTON{:}, 'Parent', saveLoadRow, 'string', 'Save As');
            uix.Empty('Parent', saveLoadRow);
            obj.btnLoad = uicontrol(obj.PROP_BUTTON{:}, 'Parent', saveLoadRow, 'string', 'Load', ...
                'TooltipString', 'a json with a "Helmholtz" block, or an old helmholtzState.mat');
            saveLoadRow.Widths = [-4 -4 -2 -4];
            obj.tvFileName = uicontrol(obj.PROP_TEXT_NORMAL{:}, 'Parent', vboxSave, 'String', '', ...
                'HorizontalAlignment', 'left');
            vboxSave.Heights = [40, 32];

            obj.btnSave.Callback = @obj.btnSaveCallback;
            obj.btnSaveAs.Callback = @obj.btnSaveAsCallback;
            obj.btnLoad.Callback = @obj.btnLoadCallback;

            obj.height = sum(vboxSave.Heights) + 50;
            obj.width = 300;
            obj.refresh;
        end

        function refresh(obj)
            f = obj.defaultFile();
            if isempty(f)
                obj.tvFileName.String = 'no state file yet: use Save As';
            else
                obj.tvFileName.String = {'state file:', f};
            end
        end

        function f = defaultFile(obj)
            f = obj.helm.stateFile;
            if isempty(f) && exist('setup7Paths', 'file') == 2
                try; f = setup7Paths().json; catch; f = ''; end
            end
        end

        function btnSaveCallback(obj, ~, ~)
            f = obj.defaultFile();
            if isempty(f); obj.btnSaveAsCallback(); return; end
            obj.helm.save(f);
            obj.refresh;
        end

        function btnSaveAsCallback(obj, ~, ~)
            [fileName, folder] = uiputfile('*.json', 'Save the Helmholtz state into...', ...
                                           obj.defaultFile());
            if isequal(fileName, 0); return; end
            obj.helm.save(fullfile(folder, fileName));
            obj.refresh;
        end

        function btnLoadCallback(obj, ~, ~)
            [fileName, folder] = uigetfile({'*.json;*.mat', 'expProperties.json or helmholtzState.mat'}, ...
                                           'Load the Helmholtz state from...', obj.defaultFile());
            if isequal(fileName, 0); return; end
            obj.helm.load(fullfile(folder, fileName));
            obj.refresh;
        end
    end

    methods
        function onEvent(obj, event)
            if ~event.isError
                obj.refresh();
            end
        end
    end
end
