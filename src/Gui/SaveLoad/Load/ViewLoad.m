classdef ViewLoad < ViewVBox & EventListener
    %VIEWLOAD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        btnLoad         % to load file
        btnPrev         % to load prev file in the folder
        btnNext         % to load next file in the folder
        btnLast         % to load last file in the folder
        btnDelete       % to delete the currently loaded file
        radioLoadAuto   % to load from the autosave folder
        radioLoadManual % to load from the folder that "save" saves into
        tvLoadedInfo    % shows the loaded struct in the saveLoad object
        btnLoadToSystem % vector of buttons - to load from the local struct to the system components
        
        category        % string. the category of the saveLoad
        subCategories   % cell of strings. the sub-categories needed for the various "load to system" buttons
    end
    
    methods
        function obj = ViewLoad(parent, controller, saveLoadCategory)
            % saveLoadCategory - string. the category to work with
            
            obj@ViewVBox(parent, controller);
            
            saveLoad = SaveLoad.getInstance(saveLoadCategory);
            obj@EventListener(saveLoad.name);
            obj.category = saveLoadCategory;
            
            %%%% UI elements init   %%%%
            greyOut = {'Enable', 'off'};
            S = MainGuiStyle;
            obj.component.Spacing = S.GAP;
            obj.component.BackgroundColor = MainGuiStyle.col().panel;
            
            % Row 1: load, browse the folder, delete
            row1 = MainGuiStyle.hbox(obj.component, 4);
            obj.btnLoad = MainGuiStyle.button(row1, 'Load...', 'normal', 'TooltipString', 'Choose a saved scan to load');
            obj.btnPrev = MainGuiStyle.button(row1, [char(9664) ' Previous'], 'normal', greyOut{:}, 'TooltipString', 'The file before this one in the folder');
            obj.btnNext = MainGuiStyle.button(row1, ['Next ' char(9654)], 'normal', greyOut{:}, 'TooltipString', 'The file after this one in the folder');
            obj.btnLast = MainGuiStyle.button(row1, 'Latest', 'normal', 'TooltipString', 'The newest file in the folder');
            obj.btnDelete = MainGuiStyle.button(row1, 'Delete', 'quiet', greyOut{:}, 'TooltipString', 'Delete the loaded file (asks first)');
            row1.Widths = [-1.1 -1.2 -1 -1 -0.9];
            
            % Row 2: which folder
            row2 = MainGuiStyle.hbox(obj.component, 6);
            MainGuiStyle.label(row2, 'Browse', 'left');
            obj.radioLoadAuto = uicontrol('Parent', row2, 'Style', 'radiobutton', 'String', 'Autosave folder', ...
                'FontSize', S.FONT, 'BackgroundColor', MainGuiStyle.col().panel);
            obj.radioLoadManual = uicontrol('Parent', row2, 'Style', 'radiobutton', 'String', 'Last folder saved to', ...
                'FontSize', S.FONT, 'BackgroundColor', MainGuiStyle.col().panel);
            row2.Widths = [50 -1 -1.2];
            
            % Row 3: what is in the loaded file. (An inactive multi-line edit, so it scrolls.)
            obj.tvLoadedInfo = uicontrol('Parent', obj.component, 'Style', 'edit', 'String', 'No file loaded', ...
                'Enable', 'inactive', 'Min', 0, 'Max', 2, 'FontSize', S.FONT - 1, ...
                'HorizontalAlignment', 'left', 'BackgroundColor', MainGuiStyle.col().readonly, ...
                'ForegroundColor', MainGuiStyle.col().ink2);
            
            % Row 4: apply the loaded file's settings
            row4 = MainGuiStyle.hbox(obj.component, 4);
            MainGuiStyle.label(row4, 'Apply its', 'left');
            labels = obj.getDisplayLoadToSystemButtons(saveLoadCategory);
            subCats = obj.getSubCategoriesForLoadToSystemButtons(saveLoadCategory);
            if length(labels) ~= length(subCats)
                EventStation.anonymousError('Lengths must be same! Coding work must be done here, please. Category: %s', saveLoadCategory)
            end
            len = length(labels);
            obj.btnLoadToSystem = gobjects(1, len);
            for i = 1 : len
                obj.btnLoadToSystem(i) = MainGuiStyle.button(row4, labels{i}, 'normal', greyOut{:});
                obj.btnLoadToSystem(i).Callback = @(h,e) obj.handleButtonLoadToSystem(subCats{i});
            end
            row4.Widths = [56, -ones(1, len)];
            
            %%%% Ui components set callbacks  %%%%
            set(obj.radioLoadAuto, 'Callback', @(h,e)obj.handleRadioAuto);
            set(obj.radioLoadManual, 'Callback', @(h,e)obj.handleRadioManual);
            set(obj.btnLoad, 'Callback', @(h,e)obj.handleButtonLoad);
            set(obj.btnDelete, 'Callback', @(h,e)obj.handleButtonDelete);
            set(obj.btnPrev, 'Callback', @(h,e)obj.handleButtonPrev);
            set(obj.btnNext, 'Callback', @(h,e)obj.handleButtonNext);
            set(obj.btnLast, 'Callback', @(h,e)obj.handleButtonLast);
            
            %%%% Internal properties set values %%%%
            heights = [S.ROW_H, 20, 74, S.ROW_H];
            obj.setHeights(heights);
            obj.height = MainGuiStyle.rowsHeight(heights);
            obj.width = 400;
            
            %%%% Set values from the saveLoad %%%%
            obj.refresh();
        end
        
        function cellOfStrings = getDisplayLoadToSystemButtons(~, category)
            % Strings to be displayed in the "load to" buttons, depending
            % on the relevant category.
            % Implement here other buttons for other categories!
            switch category
                case Savable.CATEGORY_IMAGE
                    cellOfStrings = {'scan settings', 'laser settings'};
                otherwise
                    cellOfStrings = {'settings'};
            end
        end
        
        function cellOfStrings = getSubCategoriesForLoadToSystemButtons(~, category)
            % Sub-categories matching the "load to" buttons.
            % should be a cell of same length as
            % obj.getDisplayLoadToSystemButtons with the same category
            switch category
                case Savable.CATEGORY_IMAGE
                    cellOfStrings = {Savable.CATEGORY_IMAGE_SUBCAT_STAGE, Savable.CATEGORY_IMAGE_SUBCAT_LASER};
                otherwise
                    cellOfStrings = {Savable.SUB_CATEGORY_DEFAULT};
            end
        end
        
        
        function setAvailableLoads(obj)
            saveLoad = SaveLoad.getInstance(obj.category);
            [allFiles, fileIndex] = saveLoad.availableFiles;
            
            isAnyAvailable = iscell(allFiles) || ischar(allFiles);
            isPrevAvailable = isAnyAvailable && ~isnan(fileIndex) && fileIndex ~= 1;
            isNextAvailable = isAnyAvailable && ~isnan(fileIndex) && fileIndex ~= length(allFiles);
            
            obj.btnLast.Enable = BooleanHelper.boolToOnOff(isAnyAvailable);
            obj.btnPrev.Enable = BooleanHelper.boolToOnOff(isPrevAvailable);
            obj.btnNext.Enable = BooleanHelper.boolToOnOff(isNextAvailable);
            
        end
        
        function handleButtonLoad(obj)
            % Sets the 'Enable' property of the quick loading buttons
            % ('prev', 'next') according to files existing in the folder
            saveLoad = SaveLoad.getInstance(obj.category);
            loadingFolder = saveLoad.mLoadingFolder;
            [fileName,folderName,~] = uigetfile('*.mat', 'Choose file to load...', loadingFolder);
            
            if ~isequal(fileName,0)  % user has not canceled
                saveLoad = SaveLoad.getInstance(obj.category);
                saveLoad.setLoadingFolder(folderName);
                fullPath = [folderName fileName];
                saveLoad.loadFileToLocal(fullPath);
            end
            
            obj.setAvailableLoads;
        end
        
        function handleButtonDelete(obj)
            saveLoad = SaveLoad.getInstance(obj.category);
            filename = saveLoad.mLoadedFileName;
            title = 'Delete file...';
            msg = sprintf('Are you sure you want to delete %s?', filename);
            if QuestionUserOkCancel(title, msg)
                saveLoad.deleteCurrentFile();
            end
            
            obj.setAvailableLoads;
        end
        
        function handleButtonLoadToSystem(obj, subCategory)
            saveLoad = SaveLoad.getInstance(obj.category);
            saveLoad.loadLocalToSystem(subCategory);
        end
        
        function handleButtonPrev(obj)
            saveLoad = SaveLoad.getInstance(obj.category);
            try 
                saveLoad.loadPreviousFile();
            catch e
                err2warning(e)
            end
            
            obj.setAvailableLoads;
        end
        
        function handleButtonNext(obj, ~, ~)
            saveLoad = SaveLoad.getInstance(obj.category);
            try
                saveLoad.loadNextFile();
            catch e
                err2warning(e)
            end
            
            obj.setAvailableLoads;
        end
        
        function handleButtonLast(obj)
            saveLoad = SaveLoad.getInstance(obj.category);
            try
                saveLoad.loadLastFile();
            catch e
                err2warning(e)
            end
            
            obj.setAvailableLoads;
        end
        
        function handleRadioAuto(obj)
            if obj.radioLoadAuto.Value
                saveLoad = SaveLoad.getInstance(obj.category);
                saveLoad.setLoadingFolder('autosave');
            else
                % ignore user pressing off
                obj.radioLoadAuto.Value = 1;
            end
        end
        function handleRadioManual(obj)
            if obj.radioLoadManual.Value
                saveLoad = SaveLoad.getInstance(obj.category);
                saveLoad.setLoadingFolder('manual_saves');
            else
                % ignore user pressing off
                obj.radioLoadManual.Value = 1;
            end
        end
        
        function refresh(obj)
            saveLoad = SaveLoad.getInstance(obj.category);
            if strcmp(saveLoad.mLoadingFolder, saveLoad.PATH_DEFAULT_AUTO_SAVE)
                obj.radioLoadAuto.Value = true;
                obj.radioLoadManual.Value = false;
            else
                obj.radioLoadAuto.Value = false;
                obj.radioLoadManual.Value = true;
            end
            
            % If there is a saved struct, we want to be able to see its string
            loadedStruct = saveLoad.mLocalSaveStruct;
            if isstruct(loadedStruct)
                obj.tvLoadedInfo.String = saveLoad.onlyReadableStrings(loadedStruct);
            end
            
            obj.setAvailableLoads;
        end
        
    end
        
    %% overridden from EventListener
    methods
        %% When events happen, this function jumps.
        % event is the event sent from the EventSender
        function onEvent(obj, event)
            
            %%%% handle showing the local loaded struct %%%%
            if isfield(event.extraInfo, SaveLoad.EVENT_DELETE_FILE_SUCCESS)
                obj.tvLoadedInfo.String = 'No file loaded';
            elseif isfield(event.extraInfo, SaveLoad.EVENT_LOCAL_STRUCT)
                localSavedStruct = event.extraInfo.(SaveLoad.EVENT_LOCAL_STRUCT);
                stringToDisplay = Savable.onlyReadableStrings(localSavedStruct);
                obj.tvLoadedInfo.String = stringToDisplay;
            end
            
            %%%% handle the radio buttons %%%%
            if ~event.isError
                obj.refresh();
            end
            
            %%%% handle the visibility (based on whether a struct loaded) ~~~~
            if isfield(event.extraInfo, SaveLoad.EVENT_STATUS_LOCAL_STRUCT)
                structStatus = event.extraInfo.(SaveLoad.EVENT_STATUS_LOCAL_STRUCT);
                statusesNoFile = {SaveLoad.STRUCT_STATUS_NOT_SAVED, SaveLoad.STRUCT_STATUS_NOTHING};
                if any(strcmp(structStatus, statusesNoFile))
                    obj.btnDelete.Enable = 'off';
                    obj.btnPrev.Enable = 'off';
                    obj.btnNext.Enable = 'off';
                    for btn = obj.btnLoadToSystem; btn.Enable = 'off'; end
                else
                    obj.btnDelete.Enable = 'on';
                    obj.btnPrev.Enable = 'on';
                    obj.btnNext.Enable = 'on';
                    for btn = obj.btnLoadToSystem; btn.Enable = 'on'; end
                end
            end
        end
    end
    
end