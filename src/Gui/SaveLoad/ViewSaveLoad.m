classdef ViewSaveLoad < GuiComponent & EventListener
    %VIEWSAVELOAD View to save and load
    %   Contains ViewSave and ViewLoad
    
    properties
        saveView    % the child view
        tvFileName  % shows the file name with status
        loadView    % the child view
        
        categoryToSave  % string. the category to save ('image', 'experiments', etc...)
    end
    
    properties (Constant)
        DELAY_BG_COLOR_BACK_SEC = 0.2;
        GREEN_BG_COLOR = [0.1 1 0.1];
        RED_BG_COLOR = [1 0.1 0.1];
        
        COLOR_DEFAULT = [0.32 0.32 0.31];   % secondary ink
        COLOR_UNSAVED = [0.70 0.12 0.10];   % red: not saved
        COLOR_SAVED = [0.23 0.44 0.85];     % blue: saved
        COLOR_AUTOSAVED = [0.20 0.55 0.30]; % green: autosaved
    end
    
    methods
        % constructor
        function obj = ViewSaveLoad(parent, controller, categoryToSave)
            % categoryToSave - string. Which category you should save
            % (mostly 'image' or 'experiments')
            
            obj@GuiComponent(parent, controller);
            
            saveLoad = SaveLoad.getInstance(categoryToSave);
            obj@EventListener(saveLoad.name);
            
            % A titled section; the save and load views go into its inner box
            [~, vbox] = MainGuiStyle.section(parent.component, 'Save & load scan');
            obj.component = vbox;
            obj.categoryToSave = categoryToSave;
            
            obj.saveView = ViewSave(obj, controller, categoryToSave);
            obj.tvFileName = uicontrol('Parent', obj.component, 'Style', 'text', 'String', 'No file', ...
                'FontSize', MainGuiStyle.FONT, 'BackgroundColor', MainGuiStyle.col().panel, ...
                'HorizontalAlignment', 'left');
            obj.loadView = ViewLoad(obj, controller, categoryToSave);
            heights = [obj.saveView.height, 16, obj.loadView.height];
            widths = [obj.saveView.width, obj.loadView.width];
            set(vbox, 'Heights', heights);
            
            obj.height = MainGuiStyle.sectionHeight(MainGuiStyle.rowsHeight(heights));
            obj.width = max(widths);
            
            % Dispaly current file name and status
            saveload = SaveLoad.getInstance('image');
            fileName = saveload.mLoadedFileName;
            status = saveload.mLocalStructStatus;
            obj.refresh(fileName, status);
        end
        
        function refresh(obj, fileName, status)
            if isnan(fileName)
                % No new file could be loaded
                toDisplay = 'No file';
            elseif ischar(fileName) && ischar(status)
                toDisplay = sprintf('File: %s   (%s)', fileName, status);
            elseif ischar(fileName)
                toDisplay = sprintf('File: %s', fileName);
            elseif ischar(status)
                toDisplay = sprintf('File: %s', status);
            end
            obj.tvFileName.String = toDisplay;
            obj.tvFileName.ForegroundColor = obj.statusColor(status);
        end
    end
    
    
    methods
        %% When events happen, this function jumps.
        % event is the event sent from the EventSender
        function onEvent(obj, event)
            
            if event.isError
                return
            end
            
            % Get filename and status
            if isfield(event.extraInfo, SaveLoad.EVENT_STATUS_LOCAL_STRUCT)
                status = event.extraInfo.(SaveLoad.EVENT_STATUS_LOCAL_STRUCT);
            else
                status = NaN;
            end
            
            if isfield(event.extraInfo, SaveLoad.EVENT_FILENAME)
                fileName = PathHelper.getFileNameFromFullPathFile(event.extraInfo.(SaveLoad.EVENT_FILENAME));
            else
                fileName = NaN;
            end
            
            noNewFile = isfield(event.extraInfo, SaveLoad.EVENT_DELETE_FILE_SUCCESS);
            
            % Change display, if needed
            if ischar(fileName) || ischar(status) || noNewFile
                obj.refresh(fileName, status)
            end
            
        end
    end
    
    methods (Static)
        function RGB = statusColor(status)
            % Returns appropriate color for status
            RGB =  ViewSaveLoad.COLOR_DEFAULT;
            if ~ischar(status); return; end
            switch status
                case SaveLoad.STRUCT_STATUS_SAVED
                    RGB = ViewSaveLoad.COLOR_SAVED;
                case SaveLoad.STRUCT_STATUS_AUTO_SAVED
                    RGB = ViewSaveLoad.COLOR_AUTOSAVED;
                case SaveLoad.STRUCT_STATUS_NOT_SAVED
                    RGB = ViewSaveLoad.COLOR_UNSAVED;
            end
        end
    end
end

