classdef ViewSave < GuiComponent & EventListener
    %VIEWSAVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        btnSave;        % button
        btnSaveAs;      % button
        btnUpdate       % button
        edtNotes;       % edit-text. Contains "notes"
        
        category;       % category to save to. Mostly 'image' or 'experiment'
    end
    
    methods
        % constructor
        function obj = ViewSave(parent, controller, category)
            %%%%% Get object to work with %%%%
            saveLoadPhysical = SaveLoad.getInstance(category);
            
            %%%% Constructors %%%%        
            obj@EventListener(saveLoadPhysical.name);
            obj@GuiComponent(parent, controller);
            
            %%%% Ui components init %%%%
            S = MainGuiStyle;
            vboxSaveImage = MainGuiStyle.vbox(parent.component, S.GAP);
            saveImage1stRow = MainGuiStyle.hbox(vboxSaveImage, 4);
            saveImage2stRow = MainGuiStyle.hbox(vboxSaveImage, 4);
            
            obj.btnSave = MainGuiStyle.button(saveImage1stRow, 'Save', 'normal', 'FontWeight', 'bold', ...
                'TooltipString', 'Save the scan (and its settings) to the save folder');
            obj.btnSaveAs = MainGuiStyle.button(saveImage1stRow, 'Save as...', 'normal', ...
                'TooltipString', 'Save the scan under a name and folder you choose');
            obj.btnUpdate = MainGuiStyle.button(saveImage1stRow, 'Update notes in file', 'normal', ...
                'TooltipString', 'Write the notes below into the file that is loaded/saved');
            saveImage1stRow.Widths = [-1 -1 -1.5];
            
            MainGuiStyle.label(saveImage2stRow, 'Notes', 'left');
            obj.edtNotes = MainGuiStyle.edit(saveImage2stRow, '', 'HorizontalAlignment', 'left', ...
                'FontSize', S.FONT, 'TooltipString', 'Saved with the scan, and shown as the title of its figure');
            saveImage2stRow.Widths = [40 -1];
            
            vboxSaveImage.Heights = [S.ROW_H S.ROW_H];
            obj.component = vboxSaveImage;
            %%%% UI components set callbacks %%%%
            set(obj.edtNotes, 'Callback', @obj.edtNotesCallback);
            set(obj.btnSave, 'Callback', @obj.btnSaveCallback);
            set(obj.btnSaveAs, 'Callback', @obj.btnSaveAsCallback);
            set(obj.btnUpdate, 'Callback', @obj.btnUpdateCallback);
            
            %%%% Set internal values %%%%
            obj.height = 2*MainGuiStyle.ROW_H + MainGuiStyle.GAP;
            obj.width = 400;
            obj.category = category;
            
            %%%% Init. values from the saveLoad object %%%%
            obj.refresh(saveLoadPhysical);
        end
        
        function respondToInputNotes(obj, ~, ~)
            saveLoad = SaveLoad.getInstance(obj.category);
            saveLoad.setNotes(obj.edtNotes.String);
        end
        
        %%%% Callbacks %%%%
        function edtNotesCallback(obj, ~, ~)
            obj.respondToInputNotes;
            saveLoad = SaveLoad.getInstance(obj.category);
            obj.refresh(saveLoad);
        end
        
        function btnSaveCallback(obj, ~, ~)
            obj.respondToInputNotes(); % Let the notes get into the saveload obj.
            saveLoad = SaveLoad.getInstance(obj.category);
            saveLoad.save();
        end
                
        function btnSaveAsCallback(obj, ~, ~)
            obj.respondToInputNotes(); % Let the notes get into the saveload obj.
            
            [fileName, fullPathFolder, ~] = uiputfile('.mat', 'Save As...');
            if ~ischar(fileName) || ~ischar(fullPathFolder)
                % User pressed cancel
                return
            end
            fullPathFolder = PathHelper.appendBackslashIfNeeded(fullPathFolder);
            
            saveLoad = SaveLoad.getInstance(obj.category);
            saveLoad.saveAs([fullPathFolder fileName]);
        end
        
        function btnUpdateCallback(obj, ~, ~)
            saveLoad = SaveLoad.getInstance(obj.category);
            saveLoad.saveNotesToFile;
        end
        
        %%%% end (Callbacks) %%%%
        
        function refresh(obj, saveLoadObject)
            % Refresh the view on the screen
            %
            % saveLoadObject - saveLoad object to take info from
            
            obj.edtNotes.String = saveLoadObject.mNotes;
            status = saveLoadObject.notesStatus;
            obj.edtNotes.ForegroundColor = ViewSaveLoad.statusColor(status);
        end
    end
        
        
    % overridden from EventListener
    methods
        % When events happen, this function jumps.
        %
        % event is the event sent from the EventSender
        function onEvent(obj, event) %#ok<INUSD>
            saveLoad = SaveLoad.getInstance(obj.category);
            obj.refresh(saveLoad);
        end
    end
    
end

