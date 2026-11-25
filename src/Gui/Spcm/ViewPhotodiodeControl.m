classdef ViewPhotodiodeControl < GuiComponent
    %ViewSpcmControl
    %   Detailed explanation goes here
    
    properties
        photodiodeSelection     % we have 3 photodiode: #1: fluorescence, #2: laser sampling, #3: laser block (optional)
        radioNormalInput        % #1
        radioRefInput           % #2
        radioBothInput          % #1 and #2
        radioDifferential       % #1 - <#2>
        radioReferenced         % #1 / <#2>
        radioLaserBlock         % #3
        
        cbxLifetime
        edtFrom
        edtTo
    end

    methods
        function obj = ViewPhotodiodeControl(parent, controller)
            obj@GuiComponent(parent, controller);
            hboxMain = uix.HBox('Parent', parent, 'Spacing', 5, 'Padding', 0);
            obj.component = hboxMain;
            spcm = getObjByName(Spcm.NAME);
            
            obj.photodiodeSelection = uibuttongroup(...
                'Parent', hboxMain, ...
                'Title', 'Photodiode Control', ...
                'SelectionChangedFcn',@obj.callbackPhotodiodeRadioSelection);
            
            rbHeight = 15; % "rb" stands for "radio button"
            rbWidth = 150;
            paddingFromLeft = 10;
            
            obj.radioNormalInput = uicontrol(obj.PROP_RADIO{:}, 'Parent', obj.photodiodeSelection, ...
                'String', 'signal', ...
                'Position', [paddingFromLeft 60 rbWidth rbHeight], ...  % [fromLeft, fromBottom, width, height]
                'Tag', PhotoDiodeDigitizerNiDaqControlled.PHOTODIODE_OPTIONS{1});
            obj.radioRefInput = uicontrol(obj.PROP_RADIO{:}, 'Parent', obj.photodiodeSelection, ...
                'String', 'reference', ...
                'Position', [paddingFromLeft 35 rbWidth rbHeight], ...  % [fromLeft, fromBottom, width, height]
                'Tag', PhotoDiodeDigitizerNiDaqControlled.PHOTODIODE_OPTIONS{2});
            obj.radioBothInput = uicontrol(obj.PROP_RADIO{:}, 'Parent', obj.photodiodeSelection, ...
                'String', 'signal & reference', ...
                'Position', [paddingFromLeft 10 rbWidth+30 rbHeight], ...  % [fromLeft, fromBottom, width, height]
                'Tag', PhotoDiodeDigitizerNiDaqControlled.PHOTODIODE_OPTIONS{3});
            obj.radioDifferential = uicontrol(obj.PROP_RADIO{:}, 'Parent', obj.photodiodeSelection, ...
                'String', 'differential (sig-ref)', ...
                'Position', [paddingFromLeft+rbWidth-5 60 rbWidth+40 rbHeight], ...  % [fromLeft, fromBottom, width, height]
                'Tag', PhotoDiodeDigitizerNiDaqControlled.PHOTODIODE_OPTIONS{4});
            obj.radioReferenced = uicontrol(obj.PROP_RADIO{:}, 'Parent', obj.photodiodeSelection, ...
                'String', 'referenced (sig/<ref>)', ...
                'Position', [paddingFromLeft+rbWidth-5 35 rbWidth+40 rbHeight], ...  % [fromLeft, fromBottom, width, height]
                'Tag', PhotoDiodeDigitizerNiDaqControlled.PHOTODIODE_OPTIONS{5});
            obj.radioLaserBlock = uicontrol(obj.PROP_RADIO{:}, 'Parent', obj.photodiodeSelection, ...
                'String', 'laser block', ...
                'Position', [paddingFromLeft+rbWidth-5 10 rbWidth+40 rbHeight], ...  % [fromLeft, fromBottom, width, height]
                'Tag', PhotoDiodeDigitizerNiDaqControlled.PHOTODIODE_OPTIONS{6});
            
            hboxMain.Widths = [-1];
            obj.height = 90;
            obj.width = -1;
            
            obj.update()
        end
        
        function update(obj)
            % Executes when spcm updates
            spcm = getObjByName(Spcm.NAME);
            obj.radioNormalInput.Enable = 'on';
            obj.radioReferenced.Enable = 'on';
            
            obj.radioNormalInput.Value = 0;
            obj.radioReferenced.Value = 0;
            obj.radioRefInput.Value = 0;
            obj.radioLaserBlock.Value = 0;
            obj.radioBothInput.Value = 0;
            switch spcm.currentChannelGUI
                case 1
                    obj.radioNormalInput.Value = 1;
                    obj.photodiodeSelection.SelectedObject = obj.radioNormalInput;
                case 2
                    obj.radioRefInput.Value = 1;
                    obj.photodiodeSelection.SelectedObject = obj.radioRefInput;
                case 3
                    obj.radioBothInput.Value = 1;
                    obj.photodiodeSelection.SelectedObject = obj.radioBothInput;
                case 4
                    obj.radioDifferential.Value = 1;
                    obj.photodiodeSelection.SelectedObject = obj.radioDifferential;
                case 5
                    obj.radioReferenced.Value = 1;
                    obj.photodiodeSelection.SelectedObject = obj.radioReferenced;
                case 6
                    obj.radioLaserBlock.Value = 1;
                    obj.photodiodeSelection.SelectedObject = obj.radioLaserBlock;
            end

        end
        
        %%%% Callbacks %%%%
        function callbackPhotodiodeRadioSelection(obj, ~, event)
            spcm = getObjByName(Spcm.NAME);
            action = event.NewValue.Tag;
            actionIndex = find(strcmp(action, PhotoDiodeDigitizerNiDaqControlled.PHOTODIODE_OPTIONS));
            spcm.currentChannelGUI = actionIndex;
            obj.update();
        end
        
        function callbackTypeRadioSelection(obj, ~, event)
            spcm = getObjByName(Spcm.NAME);
            action = event.NewValue.Tag;
            actionIndex = find(strcmp(action, PhotoDiodeDigitizerNiDaqControlled.TYPE_OPTIONS));
            spcm.diffrentialInputGUI = actionIndex;
            obj.update();
        end
        
        function callbackStartEndSelection(obj, ~, ~)
            spcm = getObjByName(Spcm.NAME);
            if ~ValidationHelper.isValueInteger(obj.edtFrom.String)
                obj.edtFrom.String = spcm.startRead;
                EventStation.anonymousError('Only integer numbers can be accepted! Reverting.');
            end
            if ~ValidationHelper.isValueInteger(obj.edtTo.String)
                obj.edtTo.String = spcm.endRead;
                EventStation.anonymousError('Only integer numbers can be accepted! Reverting.');
            end
            if ~ValidationHelper.isInBorders(str2double(obj.edtFrom.String), 1, spcm.endRead)
                obj.edtFrom.String = spcm.startRead;
                EventStation.anonymousError('Start can''t be lower than 1 or higher than end! Reverting.');
            end
            if ~ValidationHelper.isInBorders(str2double(obj.edtTo.String), spcm.startRead, spcm.nBins)
                obj.edtTo.String = spcm.endRead;
                EventStation.anonymousError('End can''t be lower than start or too high! Reverting.');
            end
            spcm.startRead = round(str2double(obj.edtFrom.String));
            spcm.endRead = round(str2double(obj.edtTo.String));
        end
    end
end