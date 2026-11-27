classdef GuiControllerImage < GuiController
    %GUICONTROLLERIMAGE Gui Controller for the image GUI
    %   
    properties
        dummuMode
    end
    
    methods
        function obj = GuiControllerImage(dummyMode)
            if ~exist("dummyMode", "var")
                dummyMode = 0;
            end
            shouldConfirmOnExit = true;
            openOnlyOne = true;
            % The title bar (was 'ImageNVC_touch_new'). Also how a second start finds the
            % window that is already open.
            windowName = 'NV lab';
            try
                windowName = sprintf('NV lab %s Setup %s', char(8212), JsonInfoReader.setupNumber);
            catch
            end
            obj = obj@GuiController(windowName, shouldConfirmOnExit, openOnlyOne);
            obj.dummuMode = dummyMode;
        end
        
        function view = getMainView(obj, figureWindowParent)
            % This function should get the main View of this GUI.
            % It can call any view constructor with the params:
            % parent=figureWindowParent, controller=obj
            if ~obj.dummuMode
                view = ViewMainImage(figureWindowParent, obj);
            else
                view = ViewMainImageCamera(figureWindowParent, obj);
            end
        end
        
        function onStarted(obj)
            obj.windowMinHeight = 1;
            obj.windowMinWidth = 1;
            obj.maximize();     % the window is always used maximized
            MainGuiStyle.relayout(obj.figureWindow);    % (section titles covering their content)
        end
        
        function onClose(obj)
            % Callback. Things to run when need to close the GUI.
            imageScanResult = getObjByName(ImageScanResult.NAME);
            if ~isempty(imageScanResult)
                imageScanResult.checkGraphicAxes;       % Tell it that vAxes are no longer available (without being EventSender)
            else
                % Could not find ImageScanResult, so nothing needs updating.
            end
%             StageControlEvents.sendCloseConnection;
                        % requires GUI for closing connection demand
        end
    end
    
end

