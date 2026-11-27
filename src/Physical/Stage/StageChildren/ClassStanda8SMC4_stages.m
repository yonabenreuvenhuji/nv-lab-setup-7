classdef ClassStanda8SMC4_stages < ClassStage
    % Standa 8SMC4 controllers (libximc), Setup 7.
    %
    % Scan fix, 2026-11-23:
    %  * A move now waits on the MVCMD_RUNNING bit of MvCmdSts. The old
    %    WaitFor waited on MoveSts (MOVE_STATE_MOVING), which the libximc
    %    manual says must not be used to wait for the end of a move. A move
    %    could therefore "finish" while the stage was still travelling, so
    %    the detector integrated on the way to the pixel: the first pixels
    %    of a scan saw the start position, and every line was shifted in
    %    its scan direction (about half a pixel each way).
    %  * The position is checked against the command after every move, and
    %    the scan waits scanSettleTime before each pixel is measured.
    %  * Each line is entered from scanLeadIn um before its first pixel, so
    %    every pixel is approached from the scan direction. By default all
    %    lines run in the same direction (scanBidirectional = false), so
    %    backlash cannot shift alternate lines either.
    %  * The stage no longer prepares the detector or writes into its data.
    %    It asks the detector for one pixel with spcm.acquireScanPixel(t)
    %    when that method exists, or gates it with the trigger otherwise.
   
    properties (Constant)
        controllerModel = '8SMC4';
        NAME = 'Stage - 8SMC4'
        UNITS = 'um'
        
        NEEDED_FIELDS = {'address'}
        
        COMM_DELAY = 0.005; % 5ms delay needed between consecutive commands sent to the controllers.
        STEP_MINIMUM_SIZE = 0.001                   % in um
        STEP_DEFAULT_SIZE = 0.02                    % in um
        MAX_SCAN_SIZE = 1e6;
        
        WARNING_PREVIOUS_SCAN_CANCELLED = '2D Scan is in progress! Previous scan was cancelled';

        % libximc status flags (ximc.h)
        MVCMD_RUNNING = 128     % 0x80, MvCmdSts: move command is being executed
        MVCMD_ERROR = 64        % 0x40, MvCmdSts: last move command finished with an error
        ENGINE_ANTIPLAY = 8     % 0x08, EngineFlags: controller backlash compensation on
        USTEPS = 256            % microsteps per step, as used everywhere in this class
    end
       
    properties
        % ---- scan behaviour, can be changed from the command line ----
        scanBidirectional = false   % false: every line runs in the same direction (immune to backlash). true: serpentine (saves the return trip)
        scanLeadIn = 20             % um. Each line starts this far before its first pixel, so every pixel is approached from the scan direction. 0 = off
        scanSettleTime = 0.01       % s. Wait after a move has finished and before the pixel is measured
        positionTolerance = []      % um. Allowed |reported - commanded| position after a move. [] = one full motor step
        moveTimeout = 60            % s. Longest wait for a single move
    end

    properties (Access = protected)
        physicalAddress
        id
        
        leadScrewPitch
        travelRange
        countsPerTurn     % encoder counts per turn
        stepsPerTurn      % without encoder we use motor steps per turn
        maxSpeed
        isEncoder           % beelean. if there is encoder in the stage

        posSoftRangeLimit
        negSoftRangeLimit
        
        defaultVel = 1e3;   % in um/s (== 1mm/s)
        curPos
        curVel
        forceStop
        scanRunning = false
        paramsOfRunningScan = []    % StageScanParams of the scan now running (see refuseNestedScan)
        antiplayOn         % logical per axis, read from EngineFlags at start-up (information only)
        macroScanAxis
        macroScanVector
        macroNormalScanAxis
        macroNormalScanVector
        macroScanStartEndVector
        macroScanPixelSizeInum
        macrotPixel
        macroPixelTime
        macroIndex
        macroNumberOfPixels
        fastScan
            
            
        % NiDaq
        digitalPulseTask = -1; % Digital pulse task for scanning
        
        % Tilt
        tiltCorrectionEnable = false
        tiltThetaXZ = 0
        tiltThetaYZ = 0
        
        VALID_AXES
        POSITIVE_HARD_LIMITS
        NEGATIVE_HARD_LIMITS
    end
    
    methods (Static, Access = public)
        function obj = create(stageStruct)
            % Get instance constructor
            missingField = FactoryHelper.usualChecks(stageStruct, ClassStanda8SMC4_stages.NEEDED_FIELDS);
            if ~isnan(missingField)
                EventStation.anonymousError(...
                    'Trying to create the Standa 8SMC4 stage, needed field "%s" was missing. Aborting',...
                    missingField);
            end
            
            if isfield(stageStruct, 'pg_controlled') && stageStruct.pg_controlled
                triggerChannel = 'pulseGenerator';
            else
                triggerChannel = stageStruct.niDaqChannel;
            end
            removeObjIfExists(ClassStanda8SMC4_stages.NAME);

            addressStruct = stageStruct.address;
            axes = fieldnames(addressStruct);
            validAxes = '';
            address = cell(1, length(axes));
            for k = 1:length(axes)
                validAxes = [validAxes, axes{k}];
                address{k} = addressStruct.(axes{k});
            end
            obj = ClassStanda8SMC4_stages(lower(validAxes), address, triggerChannel);
        end
        
        function axis = GetAxisInternal(axisName)
            % Gives the axis number (0 for x, 1 for y, 2 for z) when the
            % user enters the axis name (x,y,z or 1 for x, 2 for y and 3
            % for z).
            axis = zeros(size(axisName));
            for i = 1:length(axisName)
                if ((strcmpi(axisName(i),'x')) || (axisName(i) == 1))
                    axis(i) = 1;
                elseif (axisName(i) == 'y') || (axisName(i) == 'Y') || (axisName(i) == 2)
                    axis(i) = 2;
                elseif (axisName(i) == 'z') || (axisName(i) == 'Z') || (axisName(i) == 3)
                    axis(i) = 3;
                else
                    error(['Unknown axis: ' axisName]);
                end
            end
        end
    end
    
    methods (Access = protected) % Protected Functions
        function obj = ClassStanda8SMC4_stages(availableAxes, address, triggerChannel)
            % name - string
            % availableAxes - string. example: "xyz"
            name = ClassStanda8SMC4_stages.NAME;
            
            obj@ClassStage(name, availableAxes)
            
            obj.availableAxes = availableAxes;
            obj.VALID_AXES = availableAxes;
            obj.availableProperties.(obj.HAS_FAST_SCAN) = false;
            obj.availableProperties.(obj.HAS_SLOW_SCAN) = true;
            
            obj.physicalAddress = address;
            obj.Connect;
            
            obj.availableProperties.(obj.HAS_FAST_SCAN) = false;
            obj.fastScan = 0;
            obj.tiltCorrectionEnable = 0;
            obj.tiltThetaXZ = 0;
            obj.tiltThetaYZ = 0;
            
            try
                obj.readStagesParameters;
                
                if ~strcmp(triggerChannel, 'pulseGenerator')
                    nidaq = getObjByName(NiDaq.NAME);
                    nidaq.registerChannel(triggerChannel, obj.name);
                else
                    obj.pgControlled = 1;
                end
                obj.triggerChannel = triggerChannel;
                Initialization(obj);
            catch err
                obj.CloseConnection;
                rethrow(err);
            end
        end
        
        function Connect(obj)
            %load librery if neded
            if ~libisloaded('libximc')
                % libximc comes from the repo's drivers folder (git-ignored, so
                % it has to be copied into each clone by hand):
                %   <repo>\drivers\ximc-2.10.5\ximc\win64\                  libximc.dll + its dlls
                %   <repo>\drivers\ximc-2.10.5\ximc\win64\wrappers\matlab\  ximcm.m prototype + thunk
                ximcWin64 = fullfile(PathHelper.getPathToRepo(), 'drivers', 'ximc-2.10.5', 'ximc', 'win64');
                ximcMatlab = fullfile(ximcWin64, 'wrappers', 'matlab');
                if ~isfile(fullfile(ximcWin64, 'libximc.dll')) || ~isfile(fullfile(ximcMatlab, 'ximcm.m'))
                    error(['Standa stages: libximc driver not found. Expected\n' ...
                           '    libximc.dll in %s\n    ximcm.m     in %s\n' ...
                           'Copy the ximc-2.10.5 folder into <repo>\\drivers\\ (it is git-ignored).'], ...
                          ximcWin64, ximcMatlab);
                end
                addpath(ximcMatlab);
                addpath(ximcWin64);
                loadlibrary('libximc.dll',@ximcm);
                disp('libximc is now loaded')
            else
                disp('libximc is already loaded')
            end
            % connect to device
            obj.id = cell(1, length(obj.VALID_AXES));
            for k = 1:length(obj.VALID_AXES)
                device_name = sprintf('xi-com:\\\\.\\%s', obj.physicalAddress{k});
                disp("connecting to device name " + device_name)
                obj.id{k} = calllib('libximc','open_device', device_name);
                if isempty(obj.id{k})
                    error('Could not connect to stage at %s', obj.physicalAddress{k});
                end
                fprintf('connected to stage in %s\n', obj.physicalAddress{k});
            end
        end
        
        function Close(obj)
            for k = 1:length(obj.VALID_AXES)
                device_id_ptr = libpointer('int32Ptr', obj.id{k});
                calllib('libximc','close_device', device_id_ptr);
                obj.id{k} = [];
                fprintf('STANDA stage %s disconnected - remove me when time comes\n', obj.VALID_AXES(k))
            end
        end
        
        function readStagesParameters(obj)
            % This function read the relevant parameters of the specific
            % stages used in the setup.
            
            N = length(obj.VALID_AXES);
            obj.countsPerTurn = zeros(1, N);
            obj.leadScrewPitch = zeros(1, N);
            obj.travelRange = zeros(1, N);
            obj.curVel = zeros(1, N);
            obj.curPos = zeros(1, N);
            obj.maxSpeed = zeros(1, N);
            obj.countsPerTurn = zeros(1, N);     % encoder counts per turn
            obj.stepsPerTurn = zeros(1, N);      % without encoder we use motor steps per turn
            obj.isEncoder = zeros(1, N);         % beelean. if there is encoder in the stage
            obj.negSoftRangeLimit = zeros(1, N);
            obj.posSoftRangeLimit = zeros(1, N);
            obj.antiplayOn = false(1, N);
            for k = 1:N
                try
                    % Engine settings:
                    engine_settings_t = struct('NomVoltage',1,'NomCurrent',1,'NomSpeed',1,...
                        'uNomSpeed',1,'EngineFlags',1,'Antiplay',1,'MicrostepMode',1','StepsPerRev',1);
                    engine_settings = obj.SendCommand('get_engine_settings', k, engine_settings_t);
                    obj.stepsPerTurn(k) = engine_settings.StepsPerRev;
                    obj.antiplayOn(k) = bitand(double(engine_settings.EngineFlags), obj.ENGINE_ANTIPLAY) > 0;
                    
                    % Stage settings:
                    stage_settings_t = struct('LeadScrewPitch',2,'Units',100,'MaxSpeed',50,'TravelRange',360,...
                        'SupplyVoltageMin',12,'SupplyVoltageMax',36,'MaxCurrentConsumption',0,...
                        'HorizontalLoadCapacity',0,'VerticalLoadCapacity',0);
                    stage_settings = obj.SendCommand('get_stage_settings', k, stage_settings_t);
                    obj.leadScrewPitch(k) = stage_settings.LeadScrewPitch;
                    obj.travelRange(k) = stage_settings.TravelRange*1e3;
                    
                    % Feedback settings (in presence of encoder):
%                     encoder_information_t = struct('Manufacturer', 0, 'PartNumber', 0);
%                     encoder_information = obj.SendCommand('get_encoder_information', k, encoder_information_t);
%                     obj.isEncoder(k) = encoder_information.Manufacturer ~= 0;
%                     if obj.isEncoder(k)
%                         feedback_settings_t = struct('IPS',0,'FeedbackType',7,'FeedbackFlags',7);
%                         feedback_settings = obj.SendCommand('get_feedback_settings', k, feedback_settings_t);
%                         obj.countsPerTurn(k) = feedback_settings.IPS;
%                     end
                    feedback_settings_t = struct('IPS',0,'FeedbackType',7,'FeedbackFlags',7);
                    feedback_settings = obj.SendCommand('get_feedback_settings', k, feedback_settings_t);
                    obj.countsPerTurn(k) = feedback_settings.IPS;
                    obj.isEncoder(k) = feedback_settings.FeedbackType ~= 5;
                    
                    % Borders:
                    stage_edges_t = struct('BorderFlags',1,'EnderFlags',11,'LeftBorder',0,'uLeftBorder',0,...
                        'RightBorder',10,'uRightBorder',100);
                    stage_edges = obj.SendCommand('get_edges_settings', k, stage_edges_t);
                    obj.maxSpeed(k) = stage_settings.MaxSpeed*1e3;
                    obj.posSoftRangeLimit(k) = obj.stepsToValue(stage_edges.RightBorder, k)*1e3;
                    obj.negSoftRangeLimit(k) = obj.stepsToValue(stage_edges.LeftBorder, k)*1e3;
                    
                    %%% spscific stage with a problem - should be delete after fixing!! Yachel 05.09.21
%                     units = char(stage_settings.Units);
%                     if units~='m' && units~='d'
%                         obj.countsPerTurn(k) = 200;
%                         obj.leadScrewPitch(k) = 0.25;
%                         obj.travelRange(k) = 30*1e3;
%                         obj.isEncoder(k) = 0;
% %                         obj.speed = 1.25;
%                         obj.maxSpeed(k) = 5*1e3;
%                         obj.posSoftRangeLimit(k) = 14*1e3;
%                         obj.negSoftRangeLimit(k) = -14*1e3;
%                     end
                    %%%


                    %% y stage (index 2) is problematic adding code to match x axis - Yona and Yarden 11.06.26
                    units = char(stage_settings.Units);
                    if units~='m' && units~='d'
                        obj.countsPerTurn(k) = obj.countsPerTurn(1);
                        obj.leadScrewPitch(k) = obj.leadScrewPitch(1);
                        obj.travelRange(k) = obj.travelRange(1);
                        obj.isEncoder(k) = obj.isEncoder(1);
                        obj.maxSpeed(k) = obj.maxSpeed(1);
                        obj.posSoftRangeLimit(k) = obj.posSoftRangeLimit(1);
                        obj.negSoftRangeLimit(k) = obj.negSoftRangeLimit(1);
                    end
                    %%
                    
                catch errr
                    obj.Close();
                    rethrow(errr)
                end

            end
            
            obj.QueryVel;
            obj.QueryPos;
            obj.NEGATIVE_HARD_LIMITS = -obj.travelRange/2;
            obj.POSITIVE_HARD_LIMITS =  obj.travelRange/2;
        end
        
        function res_struct = GetStageStatus(obj, k)
            % here is a trick.
            % we need to init a struct with any real field from the header.
            dummy_struct = struct('Flags',999);
            parg_struct = libpointer('status_t', dummy_struct);
            [err, res_struct] = calllib('libximc','get_status', obj.id{k}, parg_struct);
            obj.testErr(err);
        end
        
        function full_res_struct = GetStatus(obj)
            full_res_struct = cell(1, length(obj.VALID_AXES));
            for k = 1:length(obj.VALID_AXES)
                full_res_struct{k} = obj.GetStageStatus(k);
            end
        end
        
        function Delay(obj, seconds) %#ok<INUSL>
            % Pauses for the given seconds.
            delay = tic;
            while toc(delay) < seconds
            end
        end
        
        function CommunicationDelay(obj)
            % Pauses for obj.commDelay seconds.
            obj.Delay(obj.COMM_DELAY);
        end
        
        function value = stepsToValue(obj, steps, k)
            if obj.isEncoder(k)
                value = steps*obj.leadScrewPitch(k)/obj.countsPerTurn(k);  
            else
                value = steps*obj.leadScrewPitch(k)/obj.stepsPerTurn(k);  
            end                
        end
        
        function steps = valueToSteps(obj, value, k)
            if obj.isEncoder(k)
                steps = value*obj.countsPerTurn(k)/obj.leadScrewPitch(k);  
            else
                steps = value*obj.stepsPerTurn(k)/obj.leadScrewPitch(k);  
            end                
        end

        function boardersCalibration(obj)
            % This function takes all the stages to the physical limits and
            % zero them in the middle of the travel range.
            prompt = 'The operation will take the stages to their limits, make sure the space is free. To continue write ''yes''\n';
            inputStr = input(prompt);
            if strcmp(inputStr, 'yes')
                for k = 1:length(obj.availableAxes)
                    fprintf('Calibrating axis %s... ', obj.availableAxes(k));
                    obj.Calibrate(k);
                    fprintf('Done!\n');
                end
            end
        end
        
        function structOut = SendCommand(obj, command, axisNumber, structIn)
            % Send the command to the controller and returns the output.
            [err, structOut] = calllib('libximc',command, obj.id{axisNumber}, structIn);
            obj.testErr(err)
        end
                
        function axisIndex = GetAxisIndex(obj, phAxis)
            % Converts x,y,z into the corresponding index for this
            % controller; if stage only has 'z' then z is 1.
            phAxis = GetAxis(obj, phAxis);
            CheckAxis(obj, phAxis)
            axisIndex = zeros(size(phAxis));
            for i = 1:length(phAxis)
                idx = phAxis(i);
                if idx > length(obj.availableAxes)
                    obj.sendError('Invalid axis')
                end
                axisIndex(i) = idx;
            end
        end
        
        function CheckAxis(obj, phAxis)
            % Checks that the given axis matches the connected stage.
            % 'phAxis' can be either a specific axis (x,y,z or 1 for x,
            % 2 for y, and 3 for z) or any vectorial combination of them.
            phAxis = obj.GetAxis(phAxis);
            if ~isempty(setdiff(phAxis, GetAxis(obj, obj.VALID_AXES)))
                if isscalar(phAxis)
                    string = 'axis is';
                else
                    string = 'axes are';
                end
                obj.sendError(sprintf('%s %s invalid for this setup.', upper(obj.SCAN_AXES(phAxis)), string));
            end
        end
        
        
        
        
        
        
        
%         function CheckReference(obj, phAxis)
%             % Checks whether the given (physical) axis is referenced, 
%             % and if not, asks for confirmation to reference it.
%             phAxis = GetAxis(obj, phAxis);
%             referenced = IsRefernced(obj, phAxis);
%             if all(referenced)
%                 % We are done
%                 return
%             end
%             
%             % Otherwise
%             unreferncedAxesNames = obj.VALID_AXES(phAxis(~referenced));
%             if isscalar(unreferncedAxesNames)
%                 questionStringPart1 = sprintf('WARNING!\n%s axis is unreferenced.\n', unreferncedAxesNames);
%             else
%                 questionStringPart1 = sprintf('WARNING!\n%s axes are unreferenced.\n', unreferncedAxesNames);
%             end
%             % Ask for user confirmation
%             questionStringPart2 = sprintf('Stages must be referenced before use.\nThis will move the stages.\nPlease make sure the movement will not cause damage to the equipment!');
%             questionString = [questionStringPart1 questionStringPart2];
%             titleString = 'Referencing Confirmation';
%             referenceString = 'Reference';
%             referenceCancelString = 'Cancel';
%             confirm = questdlg(questionString, titleString, referenceString, referenceCancelString, referenceCancelString);
%             switch confirm
%                 case referenceString
%                     Refernce(obj, unreferncedAxesNames)
%                 case referenceCancelString
%                     obj.sendError(sprintf('Referencing cancelled for controller MCS2: %s', ...
%                         unreferncedAxesNames));
%                 otherwise
%                     obj.sendError(sprintf('Referencing failed for controller MCS2: %s - No user confirmation was given', ...
%                         unreferncedAxesNames));
%             end
%         end
%         
%         function tf = IsRefernced(obj, phAxis)
%             % Check reference status for the given axis.
%             % 'phAxis' can be either a specific axis (x,y,z or 1 for x, 2 for y
%             % and 3 for z) or any vectorial combination of them.
%             len = length(phAxis);
%             tf = false(1, len);
%             realAxis = obj.GetAxisIndex(phAxis);
%             
%             for i = 1:length(realAxis)
%                 chState = obj.getChannelState(realAxis(i));
%                 tf(i) = bitget(chState, 8);     % Bit 7, if we count from 0
%             end
%             tf = logical(tf);
%         end
%         
%         function Refernce(obj, phAxis)
%             % Reference the given axis.
%             % 'phAxis' can be either a specific axis (x,y,z or 1 for x, 2 for y
%             % and 3 for z) or any vectorial combination of them.
%             CheckAxis(obj, phAxis)
%             realAxis = GetAxisIndex(obj, phAxis);
%             
%             len = length(realAxis);
%             for i = 1:len
%                 SendCommand(obj, 'SA_CTL_Reference', obj.id, realAxis(i), 0);
%             end
%             for i = 1:len
%                 WaitFor(obj, 'ReferencingDone', phAxis(i))
%             end
%             isRef = IsRefernced(obj, phAxis);
%             
%             % Check if ready & if referenced succeeded
%             if (~all(isRef))
%                 obj.sendError(sprintf('Referencing failed for controller MCS2 with ID %d: Reason unknown.', ...
%                     obj.id));
%             end
%         end
%         
%         function CheckCalibration(obj, phAxis)
%             % Checks whether the given (physical) axis is calibrated, 
%             % and if not, asks for confirmation to calibrate it.
%             phAxis = GetAxis(obj, phAxis);
%             calibrated = IsCalibrated(obj, phAxis);
%             if all(calibrated)
%                 % We are done
%                 return
%             end
%             
%             % Otherwise
%             uncalibratedAxesNames = obj.VALID_AXES(phAxis(~calibrated));
%             if isscalar(uncalibratedAxesNames)
%                 questionStringPart1 = sprintf('WARNING!\n%s axis is uncalibrated.\n', uncalibratedAxesNames);
%             else
%                 questionStringPart1 = sprintf('WARNING!\n%s axes are uncalibrated.\n', uncalibratedAxesNames);
%             end
%             % Ask for user confirmation
%             questionStringPart2 = sprintf('Stages should be calibrated before use.\nThis will move the stages.\nPlease make sure the movement will not cause damage to the equipment!');
%             questionString = [questionStringPart1 questionStringPart2];
%             titleString = 'Calibration Confirmation';
%             calibrateString = 'Calibrate';
%             calibrateCancelString = 'Cancel';
%             confirm = questdlg(questionString, titleString, calibrateString, calibrateCancelString, calibrateCancelString);
%             switch confirm
%                 case calibrateString
%                     Calibrate(obj, uncalibratedAxesNames)
%                 case calibrateCancelString
%                     obj.sendWarning(sprintf('Calibration cancelled for controller MCS2: %s', ...
%                         uncalibratedAxesNames));
%                 otherwise
%                     obj.sendWarning(sprintf('Calibration failed for controller MCS2: %s - No user confirmation was given', ...
%                         uncalibratedAxesNames));
%             end
%         end
%         
%         function tf = IsCalibrated(obj, phAxis)
%             % Check reference status for the given axis.
%             % 'phAxis' can be either a specific axis (x,y,z or 1 for x, 2 for y
%             % and 3 for z) or any vectorial combination of them.
%             len = length(phAxis);
%             tf = false(1, len);
%             realAxis = obj.GetAxisIndex(phAxis);
%             
%             for i = 1:length(realAxis)
%                 chState = obj.getChannelState(realAxis(i));
%                 tf(i) = bitget(chState, 7);     % Bit 6, if we count from 0
%             end
%             tf = logical(tf);
%         end
%         
        function Calibrate(obj, phAxis)
            % Reference the given axis.
            % 'phAxis' can be either a specific axis (x,y,z or 1 for x, 2 for y
            % and 3 for z) or any vectorial combination of them.
            CheckAxis(obj, phAxis)
            realAxis = GetAxisIndex(obj, phAxis);
            
%                 SendCommand(obj, 'SA_CTL_Calibrate', obj.id, realAxis(i), 0);
%                 WaitFor(obj, 'CalibratingDone', phAxis(i))
%             isCalib = IsCalibrated(obj, phAxis);
%             
%             % Check if ready & if calibration succeeded
%             if (~all(isCalib))
%                 obj.sendError(sprintf('Calibration failed for controller MCS2 with ID %d: Reason unknown.', ...
%                     obj.id));
%             end
        end
        
        function Initialization(obj)
            % Checks
%             CheckCalibration(obj, obj.VALID_AXES)
%             CheckReference(obj, obj.VALID_AXES)
            
            % Set velocity
            obj.SetVelocity(obj.VALID_AXES, zeros(size(obj.VALID_AXES))+obj.defaultVel);
            
            % Update position and velocity
            QueryPos(obj);
            QueryVel(obj);
            for i = 1:length(obj.VALID_AXES)
                fprintf('%s axis - Position: %.4f%s, Velocity: %d%s/s.\n', upper(obj.VALID_AXES(i)), obj.curPos(i), obj.UNITS, obj.curVel(i), obj.UNITS);
            end
            if ~isempty(obj.antiplayOn)
                onOff = {'OFF', 'ON'};
                fprintf('Standa backlash compensation (controller antiplay): %s\n', ...
                    strjoin(arrayfun(@(i) sprintf('%s %s', upper(obj.VALID_AXES(i)), onOff{1+obj.antiplayOn(i)}), ...
                    1:length(obj.antiplayOn), 'UniformOutput', false), ', '));
            end
        end
        
        function QueryPos(obj)
            % Queries the position and updates the internal variable.
            for k = 1:length(obj.VALID_AXES)
                get_position_calt_t = struct('Position',1,'uPosition', 1, 'EncPosition', 1);
                get_position_calt = obj.SendCommand('get_position', k, get_position_calt_t);
                pos = obj.stepsToValue(get_position_calt.Position + get_position_calt.uPosition/256, k);
                obj.curPos(k) = pos * 1e3; % Convert from mm to um
            end
        end
        
        function QueryVel(obj)
            % Queries the velocity and updates the internal variable.
            for k = 1:length(obj.VALID_AXES)
                move_settings_calb_t = struct('Speed', 0, 'Accel', 0, 'Decel', 0, 'AntiplaySpeed', 0);
                move_settings_calb = obj.SendCommand('get_move_settings', k, move_settings_calb_t);
                vel = obj.stepsToValue(move_settings_calb.Speed, k);
                obj.curVel(k) = vel * 1e3; % Convert from mm/s to um/s
            end
        end
        
        function WaitFor(obj, what, phAxis, targetMicro)
            % Waits until the move command on one axis has finished.
            % 'phAxis' must a specific (physical) axis (x,y,z or 
            % 1 for x, 2 for y and 3 for z).
            % Options for 'what' (all wait for the move command to end):
            % MovementDone, CalibratingDone, ReferencingDone.
            % targetMicro (optional): the commanded position in
            % microsteps (Position*256 + uPosition). When given, the
            % reported position is checked against it.
            %
            % Waits on MVCMD_RUNNING in MvCmdSts. The libximc manual on
            % MOVE_STATE_MOVING (MoveSts, used here before): "Don't use
            % this flag to wait for the completion of the movement
            % command. Use the MVCMD_RUNNING flag from the MvCmdSts field
            % instead."
            CheckAxis(obj, phAxis);
            realAxis = GetAxisIndex(obj, phAxis);
            if ~any(strcmp(what, {'MovementDone', 'CalibratingDone', 'ReferencingDone'}))
                obj.sendError(sprintf('Wrong Input %s', what));
            end
                
            timer = tic;
            while true
                drawnow % Needed in order to get input from GUI
                if obj.forceStop % Checks if the user pressed the Halt Button
                    HaltPrivate(obj, phAxis);
                    return
                end
                chState = GetStageStatus(obj, realAxis);
                if ~bitand(double(chState.MvCmdSts), obj.MVCMD_RUNNING)
                    break   % the controller has finished the command
                end
                if toc(timer) > obj.moveTimeout
                    error('Standa %s axis: move did not finish within %g s.', upper(obj.VALID_AXES(realAxis)), obj.moveTimeout);
                end
            end

            if bitand(double(chState.MvCmdSts), obj.MVCMD_ERROR)
                error(['Standa %s axis: the controller reports that the move finished with an error ', ...
                    '(limit switch, alarm or backlash compensation failure).'], upper(obj.VALID_AXES(realAxis)));
            end
            if exist('targetMicro', 'var') && ~isempty(targetMicro)
                obj.checkReachedTarget(realAxis, chState, targetMicro);
            end
        end
        
        function checkReachedTarget(obj, k, chState, targetMicro)
            % Errors if axis k did not end at the commanded position.
            % Reported and commanded positions are compared in
            % microsteps, so no unit conversion enters the check.
            if isempty(obj.positionTolerance)
                tolMicro = obj.USTEPS;     % one full motor step
            else
                tolMicro = abs(obj.valueToSteps(obj.positionTolerance * 1e-3, k)) * obj.USTEPS;
            end
            current = double(chState.CurPosition) * obj.USTEPS + double(chState.uCurPosition);
            t0 = tic;
            while abs(current - targetMicro) > tolMicro && toc(t0) < 0.2
                pause(0.005);
                chState = GetStageStatus(obj, k);
                current = double(chState.CurPosition) * obj.USTEPS + double(chState.uCurPosition);
            end
            if abs(current - targetMicro) > tolMicro
                offUm = obj.stepsToValue((current - targetMicro) / obj.USTEPS, k) * 1e3;
                error('Standa %s axis stopped %.3f um away from the commanded position.', upper(obj.VALID_AXES(k)), offUm);
            end
        end

        function MovePrivate(obj, phAxis, pos)
            % Absolute change in position (pos) of axis (x,y,z or 1 for x,
            % 2 for y and 3 for z).
            % Vectorial axis is possible.
            % Does not check if scan is running.
            % Does not move if HaltStage was triggered.
            % This function is the one used by all internal functions.
            % Returns only when every axis has finished its move and is
            % at the commanded position (see WaitFor).
            CheckAxis(obj, phAxis)
            
            if obj.forceStop % Doesn't move if forceStop is enabled
                return
            end
            
            if obj.tiltCorrectionEnable
                [phAxis, pos] = TiltCorrection(phAxis, pos);
            end
            
            if ~PointIsInRange(obj, phAxis, pos) % Check that point is in limits
                obj.sendError('Move Command is outside the soft limits');
            end
            
%             CheckReference(obj, phAxis)
            
            % Convert position to mm, which are the units the stage recieves
            pos = pos * 1e-3;

            % Send the move command
            for i = 1:length(phAxis)
                if obj.forceStop % Halt pressed while a previous axis was moving
                    return
                end
                realAxes = obj.GetAxisInternal(phAxis(i));
                setPositionPoints = floor(obj.valueToSteps(pos(i), realAxes));
                set_uPosition = round((obj.valueToSteps(pos(i), realAxes) - setPositionPoints) * 256);
                if set_uPosition == 256
                    setPositionPoints = setPositionPoints + 1;
                    set_uPosition = 0;
                end
                [err] = calllib('libximc','command_move', obj.id{realAxes}, setPositionPoints, set_uPosition);
                obj.testErr(err);
                WaitFor(obj, 'MovementDone', realAxes, setPositionPoints * obj.USTEPS + set_uPosition)
                % target_mm = pos(i);

                % disp("ffffffffffffffffffffffffffffff")
                % disp(stringAxis(i))
                % 
                % if stringAxis(i) == "1" || stringAxis(i) == "x"
                %     hid = 1;
                %     Step_multiplier = 200;
                %     Unit_multiplier = 0.25;
                %     microsteps_per_step = 256;
                % end
                % 
                % if stringAxis(i) == "2" || stringAxis(i) == "y"
                %     hid = 2;
                %     Step_multiplier = 200;
                %     Unit_multiplier = 0.25;
                %     microsteps_per_step = 256;
                % end
                % 
                % if stringAxis(i) == "3" || stringAxis(i) == "z"
                %     hid = 3;
                %     Step_multiplier = 200;
                %     Unit_multiplier = 1;
                %     microsteps_per_step = 256;
                % end
                % 
                % mm_per_full_step = Unit_multiplier / Step_multiplier;
                % total_usteps = target_mm / mm_per_full_step * microsteps_per_step;
                % 
                % % Split into integer steps and microsteps remainder
                % full_steps = floor(total_usteps / microsteps_per_step);
                % usteps     = round(mod(total_usteps, microsteps_per_step));
                % 
                % disp([hid,full_steps,usteps])
                % 
                % % Move and wait
                % calllib('libximc','command_move',hid,int32(full_steps),int32(usteps));
                % disp("oooooooooooooooooooooooooooooooooooo")
                % disp([hid,int32(100)])
                % calllib('libximc','command_wait_for_stop',hid,int32(100));


            end
        end

        function HaltPrivate(obj, phAxis)
            % Halts the stage.
            realAxes = GetAxisIndex(obj, phAxis);
            for i = 1:length(realAxes)
                calllib('libximc', 'command_stop', obj.id{realAxes(i)});  % was SendCommand(..., realAxis(i), []): undefined variable
            end
            obj.AbortScan;
            obj.sendWarning('Stage Halted!');
        end
        
        function SetVelocityPrivate(obj, phAxis, vel)
            % Absolute change in velocity (vel) of axis (x,y,z or 1 for x,
            % 2 for y and 3 for z).
            % Does not check if scan is running.
            % Vectorial axis is possible.
            CheckAxis(obj, phAxis)
            vel = vel * 1e-3; % Convert from um/s to mm/s
            for i = 1:length(phAxis)
                realAxis = GetAxisIndex(obj, phAxis(i));
                move_settings_calb_t = struct('Speed', 0, 'Accel', 0, 'Decel', 0, 'AntiplaySpeed', 0);
                move_settings_calb = obj.SendCommand('get_move_settings', realAxis, move_settings_calb_t);
                move_settings_calb.Speed = obj.valueToSteps(vel(i), realAxis);
                obj.SendCommand('set_move_settings', realAxis, move_settings_calb);
            end
        end
        
        function [phAxis, pos] = TiltCorrection(obj, phAxis, pos)
            % Corrects the movement axis & pos vectors according to the 
            % tilt angles: If x and/or y axes are given, then z also moves
            % according the tilt angles.
            % However, if also or only z axis is given, then no changes
            % occurs.
            % Assumes the stage has all three xyz axes.
            phAxis = GetAxis(obj, phAxis);
            if ~contains(obj.SCAN_AXES(phAxis), 'z') % Only do something if there is no z axis
                QueryPos(obj);
                pos = [pos, obj.curPos(3)]; % Adds the z position command, start by writing the current position (as the base)
                for i=1:length(phAxis)
                    switch obj.SCAN_AXES(phAxis)
                        case 'x'
                            dx = pos(i) - obj.curPos(1);
                            pos(end) = pos(end) + dx*tan(obj.tiltThetaXZ*pi/180); % Adds movements according to the angles
                        case 'y'
                            dy = pos(i) - obj.curPos(2);
                            pos(end) = pos(end) + dy*tan(obj.tiltThetaYZ*pi/180); % Adds movements according to the angles
                    end
                end
                phAxis = [phAxis, 3]; % Add Z axis at the end
            end
        end
    end
    
    methods (Access = public)
        function CloseConnection(obj)
            % Closes the connection to the controllers.
            for k = 1:length(obj.VALID_AXES)
                device_id_ptr = libpointer('int32Ptr', obj.id{k});
                calllib('libximc','close_device', device_id_ptr);
                obj.id{k} = [];
                fprintf('STANDA stage %s disconnected\n', obj.VALID_AXES(k))
            end
        end
        
        function delete(obj)
            obj.CloseConnection;
        end
        
        function Reconnect(obj)
            % Reconnects the controller.
            CloseConnection(obj);
            Connect(obj);
            Initialization(obj);
        end

        function report = checkMoveWait(obj, phAxis, stepUm, nSteps)
            % Hardware check of the move-completion criterion. No detector.
            % Makes nSteps relative moves of stepUm on one axis, then goes
            % back to the start. For each move it polls get_status and
            % records when the OLD wait (MoveSts == 0) would have returned,
            % when the NEW wait (MVCMD_RUNNING cleared) returns, and how far
            % the stage still was from the target when the old wait said
            % "done".
            %   s = getObjByName('Stage - 8SMC4');
            %   r = s.checkMoveWait('x', 25, 10);   % defaults
            if ~exist('phAxis', 'var') || isempty(phAxis); phAxis = 'x'; end
            if ~exist('stepUm', 'var') || isempty(stepUm); stepUm = 25; end
            if ~exist('nSteps', 'var') || isempty(nSteps); nSteps = 10; end
            if obj.scanRunning
                error('A scan is running.');
            end
            CheckAxis(obj, phAxis);
            k = GetAxisIndex(obj, phAxis);
            QueryPos(obj);
            startUm = obj.curPos(k);
            targetsUm = startUm + (1:nSteps) * stepUm;
            if ~all(arrayfun(@(p) obj.PointIsInRange(k, p), targetsUm))
                error('The test moves would leave the soft limits.');
            end

            report = struct('targetUm', num2cell(targetsUm), 'tOldDone', NaN, 'tNewDone', NaN, ...
                'shortAtOldDoneUm', NaN, 'nPolls', 0);
            for i = 1:nSteps
                stepsF = obj.valueToSteps(targetsUm(i) * 1e-3, k);
                P = floor(stepsF);
                U = round((stepsF - P) * obj.USTEPS);
                if U == obj.USTEPS; P = P + 1; U = 0; end
                targetMicro = P * obj.USTEPS + U;
                err = calllib('libximc', 'command_move', obj.id{k}, P, U);
                obj.testErr(err);
                t0 = tic;
                tOld = NaN; posAtOld = NaN; n = 0;
                while true
                    st = GetStageStatus(obj, k);
                    t = toc(t0); n = n + 1;
                    pos = double(st.CurPosition) * obj.USTEPS + double(st.uCurPosition);
                    if isnan(tOld) && double(st.MoveSts) == 0
                        tOld = t; posAtOld = pos;
                    end
                    if ~bitand(double(st.MvCmdSts), obj.MVCMD_RUNNING)
                        break
                    end
                    if t > obj.moveTimeout
                        error('Move did not finish within %g s.', obj.moveTimeout);
                    end
                end
                report(i).tOldDone = tOld;
                report(i).tNewDone = t;
                report(i).nPolls = n;
                report(i).shortAtOldDoneUm = abs(obj.stepsToValue((targetMicro - posAtOld) / obj.USTEPS, k)) * 1e3;
            end
            MovePrivate(obj, k, startUm);

            tOld = [report.tOldDone] * 1e3;
            tNew = [report.tNewDone] * 1e3;
            short = [report.shortAtOldDoneUm];
            early = short > 0.5;
            fprintf('\ncheckMoveWait, %s axis, %d moves of %g um (antiplay %s):\n', upper(obj.VALID_AXES(k)), nSteps, stepUm, ...
                BooleanHelper.ifTrueElse(obj.antiplayOn(k), 'ON', 'OFF'));
            fprintf('  old wait (MoveSts == 0) said done after   %6.1f ms (median), range %.1f-%.1f ms\n', median(tOld, 'omitnan'), min(tOld), max(tOld));
            fprintf('  new wait (MVCMD_RUNNING clear) done after %6.1f ms (median), range %.1f-%.1f ms\n', median(tNew), min(tNew), max(tNew));
            fprintf('  stage still > 0.5 um from target when the old wait said done: %d of %d moves (worst %.1f um)\n', ...
                sum(early), nSteps, max(short));
        end
        
        function ok = PointIsInRange(obj, phAxis, point)
            % Checks if the given point is within the soft (and hard)
            % limits of the given axis (x,y,z or 1 for x, 2 for y and 3 for z).
            % Vectorial axis is possible.
            CheckAxis(obj, phAxis)
            axisIndex = GetAxis(obj, phAxis);
            ok = all((point >= obj.negSoftRangeLimit(axisIndex)) & (point <= obj.posSoftRangeLimit(axisIndex)));
        end
        
        function [negSoftLimit, posSoftLimit] = ReturnLimits(obj, phAxis)
            % Return the soft limits of the given axis (x,y,z or 1 for x,
            % 2 for y and 3 for z).
            % Vectorial axis is possible.
            CheckAxis(obj, phAxis)
            axisIndex = GetAxis(obj, phAxis);
            negSoftLimit = obj.negSoftRangeLimit(axisIndex);
            posSoftLimit = obj.posSoftRangeLimit(axisIndex);
        end
        
        function [negHardLimit, posHardLimit] = ReturnHardLimits(obj, phAxis)
            % Return the hard limits of the given axis (x,y,z or 1 for x,
            % 2 for y and 3 for z).
            % Vectorial axis is possible.
            CheckAxis(obj, phAxis)
            axisIndex = GetAxis(obj, phAxis);
            negHardLimit = obj.NEGATIVE_HARD_LIMITS(axisIndex);
            posHardLimit = obj.POSITIVE_HARD_LIMITS(axisIndex);
        end
        
        function SetSoftLimits(obj, phAxis, softLimit, negOrPos)
            % Set the new soft limits:
            % if negOrPos = 0 -> then softLimit = lower soft limit
            % if negOrPos = 1 -> then softLimit = higher soft limit
            % This is because each time this function is called only one of
            % the limits updates
            CheckAxis(obj, phAxis)
            axisIndex = GetAxis(obj, phAxis);
            if ((softLimit >= obj.NEGATIVE_HARD_LIMITS(axisIndex)) && (softLimit <= obj.POSITIVE_HARD_LIMITS(axisIndex)))
                if negOrPos == 0
                    obj.negSoftRangeLimit(axisIndex) = softLimit;
                else
                    obj.posSoftRangeLimit(axisIndex) = softLimit;
                end
            else
                obj.sendError(sprintf('Soft limit %.4f is outside of the hard limits %.4f - %.4f', ...
                    softLimit, obj.NEGATIVE_HARD_LIMITS(axisIndex), obj.POSITIVE_HARD_LIMITS(axisIndex)))
            end
        end
        
        function pos = Pos(obj, phAxis)
            % Query and return position of axis (x,y,z or 1 for x, 2 for y
            % and 3 for z)
            % Vectorial axis is possible.
            CheckAxis(obj, phAxis)
            QueryPos(obj);
            axisIndex = GetAxis(obj, phAxis);
            pos = obj.curPos(axisIndex);
        end
        
        function vel = Vel(obj, phAxis)
            % Query and return velocity of axis (x,y,z or 1 for x, 2 for y
            % and 3 for z)
            % Vectorial axis is possible.
            CheckAxis(obj, phAxis)
            QueryVel(obj);
            axisIndex = GetAxis(obj, phAxis);
            vel = obj.curVel(axisIndex);
        end
        
        function Move(obj, phAxis, pos)
            % Absolute change in position (pos) of axis (x,y,z or 1 for x,
            % 2 for y and 3 for z).
            % Vectorial axis is possible.
            % Checks that a scan is not currently running and whether
            % HaltStage was triggered.
            if obj.scanRunning
                obj.sendWarning(obj.WARNING_PREVIOUS_SCAN_CANCELLED);
                obj.AbortScan;
            end
            
            if obj.forceStop % Ask for user confirmation if forcestop was triggered
                questionString = sprintf('Stages were forcefully halted!\nAre you sure you want to move?');
                yesString = 'Yes';
                noString = 'No';
                confirm = questdlg(questionString, 'Movement Confirmation', yesString, noString, yesString);
                switch confirm
                    case yesString
                        obj.forceStop = 0;
                    case noString
                        obj.sendWarning('Movement aborted!')
                        return;
                    otherwise
                        obj.sendWarning('Movement aborted!')
                        return;
                end
            end
            
            obj.MovePrivate(phAxis, pos);
        end
        
        function RelativeMove(obj, phAxis, change)
            % Relative change in position (pos) of axis (x,y,z or 1 for x,
            % 2 for y and 3 for z).
            % Vectorial axis is possible.
            if obj.scanRunning
                obj.sendWarning(obj.WARNING_PREVIOUS_SCAN_CANCELLED);
                AbortScan(obj);
            end
            CheckAxis(obj, phAxis)
            QueryPos(obj);
            axisIndex = GetAxis(obj, phAxis);
            Move(obj, phAxis, obj.curPos(axisIndex) + change);
        end
        
        function Halt(obj)
            % Halts all stage movements.
            % This works by setting the parameter below to 1, which is
            % checked inside the "WaitFor" function. When the WaitFor is
            % triggered, is calls an internal function, "HaltPrivate", which
            % immediately sends a halt command to the controller.
            % Afterwards it also tries to abort scan. The reason abort scan
            % happens afterwards is to minimize the time it takes to
            % send the halt command to the controller.
            % This parameters also denies the "MovePrivate" command from
            % running.
            % It is reset by a normal/relative "Move Command", which will
            % be triggered whenever a new external move or scan command is
            % sent to the stage.
            obj.forceStop = 1;
        end
        
        function SetVelocity(obj, phAxis, vel)
            % Absolute change in velocity (vel) of axis (x,y,z or 1 for x,
            % 2 for y and 3 for z).
            % Vectorial axis is possible.
            if obj.scanRunning
                obj.sendWarning('2D Scan is in progress, previous scan canceled');
                AbortScan(obj);
            end
            
            CheckAxis(obj, phAxis)
            SetVelocityPrivate(obj, phAxis, vel);
        end
        
        function ScanX(obj, x, y, z, nFlat, nOverRun, tPixel)
            %%%%%%%%%%%%%% ONE DIMENSIONAL X SCAN MACRO %%%%%%%%%%%%%%
            % Does a macro scan for x axis.
            % x - A vector with the points to scan, points should have
            % equal distance between them.
            % y/z - The starting points for the other axes.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.scanRunning
                obj.sendWarning('2D Scan is in progress, previous scan canceled');
                AbortScan(obj);
            end
            Move(obj, 'yz', [y z]);
            ScanOneDimension(obj, x, nFlat, nOverRun, tPixel, 'x');
            QueryPos(obj);
        end
        
        function ScanY(obj, x, y, z, nFlat, nOverRun, tPixel)
            %%%%%%%%%%%%%% ONE DIMENSIONAL Y SCAN MACRO %%%%%%%%%%%%%%
            % Does a macro scan for y axis.
            % y - A vector with the points to scan, points should have
            % equal distance between them.
            % x/z - The starting points for the other axes.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.scanRunning
                obj.sendWarning('2D Scan is in progress, previous scan canceled');
                AbortScan(obj);
            end
            Move(obj, 'xz', [x z]);
            ScanOneDimension(obj, y, nFlat, nOverRun, tPixel, 'y');
            QueryPos(obj);
        end
        
        function ScanZ(obj, x, y, z, nFlat, nOverRun, tPixel)
            %%%%%%%%%%%%%% ONE DIMENSIONAL Z SCAN MACRO %%%%%%%%%%%%%%
            % Does a macro scan for z axis.
            % z - A vector with the points to scan, points should have
            % equal distance between them.
            % x/y - The starting points for the other axes.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.scanRunning
                obj.sendWarning('2D Scan is in progress, previous scan canceled');
                AbortScan(obj);
            end
            Move(obj, 'xy', [x y]);
            ScanOneDimension(obj, z, nFlat, nOverRun, tPixel, 'z');
            QueryPos(obj);
        end
        
        function PrepareScanX(obj, x, y, z, nFlat, nOverRun, tPixel)
            %%%%%%%%%%%%%%%%% ONE DIMENSIONAL X SCAN %%%%%%%%%%%%%%%%%
            % Prepares a scan for X axis, writes the waveform and 
            % initializes the DDL function.
            % x - A vector with the points to scan, points should have
            % equal distance between them.
            % y/z - The starting points for the other axes.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.scanRunning
                warning(obj.WARNING_PREVIOUS_SCAN_CANCELLED);
                AbortScan(obj);
            end
            Move(obj, 'yz', [y z]);
            PrepareScanInOneDimension(obj, x, nFlat, nOverRun, tPixel, 'x');
            QueryPos(obj);
        end
        
        function PrepareScanY(obj, x, y, z, nFlat, nOverRun, tPixel)
            %%%%%%%%%%%%%%%%% ONE DIMENSIONAL Y SCAN %%%%%%%%%%%%%%%%%
            % Prepares a scan for Y axis, writes the waveform and 
            % initializes the DDL function.
            % y - A vector with the points to scan, points should have
            % equal distance between them.
            % x/z - The starting points for the other axes.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.scanRunning
                warning(obj.WARNING_PREVIOUS_SCAN_CANCELLED);
                AbortScan(obj);
            end
            Move(obj, 'xz', [x z]);
            PrepareScanInOneDimension(obj, y, nFlat, nOverRun, tPixel, 'y');
            QueryPos(obj);
        end
        
        function PrepareScanZ(obj, x, y, z, nFlat, nOverRun, tPixel)
            %%%%%%%%%%%%%%%%% ONE DIMENSIONAL Z SCAN %%%%%%%%%%%%%%%%%
            % Prepares a scan for Z axis, writes the waveform and 
            % initializes the DDL function.
            % z - A vector with the points to scan, points should have
            % equal distance between them.
            % x/y - The starting points for the other axes.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.scanRunning
                warning(obj.WARNING_PREVIOUS_SCAN_CANCELLED);
                AbortScan(obj);
            end
            Move(obj, 'xy', [x y]);
            PrepareScanInOneDimension(obj, z, nFlat, nOverRun, tPixel, 'z');
            QueryPos(obj);
        end
        
        
        
        
        
        
        
        
        
        
        function PrepareScanInOneDimension(obj, scanAxisVector, nFlat, nOverRun, tPixel, scanAxis) %#ok<INUSL>
            % Prepares a one dimensional scan, writes the waveform and 
            % initializes the DDL function.
            % scanAxisVector - A vector with the points to scan, points
            % should increase with equal distances between them.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel (in seconds).
            % scanAxis - The axis to scan (x,y,z or 1 for x, 2 for y and 3
            % for z).
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Get parameters
            scanAxis = GetAxis(obj, scanAxis);
            
%             if (nOverRun < 10); nOverRun = 10; end % In order to be centered around the pixel we need at least one extra point from each side, and it's 2 because negative direction has a bug
            
            
            numberOfPixels = length(scanAxisVector); % This is the number of pixels
            pixelSizeInum = (scanAxisVector(end) - scanAxisVector(1))/(numberOfPixels-1); % Will be negative if scanning in reverse
            nOverRunInum = pixelSizeInum*nOverRun; % Can be negative if scanning in reverse
            startPointInum = scanAxisVector(1)-nOverRunInum;
            endPointInum = scanAxisVector(end)+nOverRunInum;
            
            maxPointInum = max(endPointInum, startPointInum);    
            minPointInum = min(endPointInum, startPointInum);
            if maxPointInum > obj.posSoftRangeLimit(scanAxis) || minPointInum < obj.negSoftRangeLimit(scanAxis)
                obj.sendError(sprintf('Scan is outside the limits: scan overhead requires adding an additional pixel from each side which increase the scan from %.3f %s to %.3f %s\n', startPointInum, StringHelper.MICRON, endPointInum, StringHelper.MICRON));
            end
            
            obj.macroScanAxis = scanAxis;
            obj.macrotPixel = tPixel;
            obj.macroScanVector = scanAxisVector;
            obj.macroScanStartEndVector = [startPointInum, endPointInum];
            obj.macroScanPixelSizeInum = pixelSizeInum;
            
            % Move to starting position, might be redundant sometimes, but
            % is used in some cases.
            MovePrivate(obj, scanAxis, startPointInum);
        end
        
        function ScanOneDimension(obj, scanAxisVector, nFlat, nOverRun, tPixel, scanAxis)  %#ok<INUSL>
            %%%%%%%%%%%%%% ONE DIMENSIONAL SCAN %%%%%%%%%%%%%%
            % Stop-and-go scan of one line, in the order of scanAxisVector.
            % Last 2 variables are for 2D scans.
            % scanAxisVector - A vector with the points to scan, in the
            % order they are visited.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel (in seconds).
            % scanAxis - The axis to scan (x,y,z or 1 for x, 2 for y and 3
            % for z).
            %
            % For every pixel: move, wait until the controller has
            % finished the move and reports the commanded position, wait
            % scanSettleTime, then measure. The line is entered from
            % scanLeadIn um before its first pixel, so every pixel is
            % approached from the scan direction.
            %
            % The detector does the measuring and keeps the data. If it
            % has acquireScanPixel (PhotoDiodeDigitizerNiDaqControlled),
            % it is called once per pixel and readFromScan returns the
            % line. Other detectors are gated with the stage trigger for
            % tPixel, as in the other stage classes.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            scanAxis = GetAxis(obj, scanAxis);
            spcm = getObjByName(Spcm.NAME);
            detectorTakesPixels = ~isempty(spcm) && ismethod(spcm, 'acquireScanPixel');
            if ~detectorTakesPixels
                obj.prepareWriteTrigger;
            end

            obj.moveToLineStart(scanAxis, scanAxisVector);
            for i = 1:length(scanAxisVector)
                if obj.forceStop
                    error('Stage halted - scan line aborted.');
                end
                MovePrivate(obj, scanAxis, scanAxisVector(i));
                if obj.scanSettleTime > 0
                    pause(obj.scanSettleTime);
                end
                if detectorTakesPixels
                    spcm.acquireScanPixel(tPixel);
                else
                    obj.writeTrigger(1);
                    pause(tPixel);
                    obj.writeTrigger(0);
                end
                StageScanner.notePixel();   % the GUI's progress bar
            end
        end

        function moveToLineStart(obj, scanAxis, scanAxisVector)
            % Moves to scanLeadIn um before the first pixel of the line, on
            % the side the line comes from, so the first pixel is reached
            % moving in the scan direction like all the others. The point
            % is clipped to the soft limits.
            leadIn = abs(obj.scanLeadIn);
            if isempty(leadIn) || leadIn == 0 || length(scanAxisVector) < 2
                return
            end
            direction = sign(scanAxisVector(end) - scanAxisVector(1));
            if direction == 0
                return
            end
            target = scanAxisVector(1) - direction * leadIn;
            [negLim, posLim] = obj.ReturnLimits(scanAxis);
            target = min(max(target, negLim), posLim);
            MovePrivate(obj, scanAxis, target);
        end


        % function ScanOneDimension(obj, scanAxisVector, nFlat, nOverRun, tPixel, scanAxis)  %#ok<INUSL>
        %     %%%%%%%%%%%%%% ONE DIMENSIONAL SCAN %%%%%%%%%%%%%%
        %     % Does a scan for the given axis.
        %     % Last 2 variables are for 2D scans.
        %     % scanAxisVector - A vector with the points to scan, points
        %     % should increase with equal distances between them.
        %     % nFlat - How many flat points should be in the beginning of the scan.
        %     % nOveRun - How many extra points should be taken from each.
        %     % tPixel - Scan time for each pixel (in seconds).
        %     % scanAxis - The axis to scan (x,y,z or 1 for x, 2 for y and 3
        %     % for z).
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     disp("do you love me?")
        %     scanAxis = GetAxis(obj, scanAxis);
        %     obj.prepareWriteTrigger;
        %     for i = 1:length(scanAxisVector)
        %         Move(obj, scanAxis, scanAxisVector(i));
        %         obj.writeTrigger(1);
        %         pause(obj.macrotPixel);
        %         obj.writeTrigger(0);
        %     end
        %     obj.SetVelocity(obj.macroScanAxis, obj.defaultVel);
        % end
        
        function PrepareScanXY(obj, x, y, z, nFlat, nOverRun, tPixel)
            %%%%%%%%%%%%%% TWO DIMENSIONAL XY SCAN MACRO %%%%%%%%%%%%%%
            % Prepare a macro scan for xy axes!
            % Scanning is done by calling 'ScanNextLine'.
            % Aborting via 'AbortScan'.
            % x/y - Vectors with the points to scan, points should have
            % equal distance between them.
            % z - The starting points for the other axis.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.scanRunning
                obj.sendWarning('2D Scan is in progress, previous scan canceled');
                AbortScan(obj);
            end
            Move(obj, 'z', z);
            PrepareScanInTwoDimensions(obj, x, y, nFlat, nOverRun, tPixel, 'x', 'y');
            QueryPos(obj);
        end
        
        function PrepareScanXZ(obj, x, y, z, nFlat, nOverRun, tPixel)
            %%%%%%%%%%%%%% TWO DIMENSIONAL XZ SCAN MACRO %%%%%%%%%%%%%%
            % Prepare a macro scan for xz axes!
            % Scanning is done by calling 'ScanNextLine'.
            % Aborting via 'AbortScan'.
            % x/z - Vectors with the points to scan, points should have
            % equal distance between them.
            % y - The starting points for the other axis.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.scanRunning
                obj.sendWarning('2D Scan is in progress, previous scan canceled');
                AbortScan(obj);
            end
            Move(obj, 'y', y);
            PrepareScanInTwoDimensions(obj, x, z, nFlat, nOverRun, tPixel, 'x', 'z');
            QueryPos(obj);
        end
        
        function PrepareScanYX(obj, x, y, z, nFlat, nOverRun, tPixel)
            %%%%%%%%%%%%%% TWO DIMENSIONAL XY SCAN MACRO %%%%%%%%%%%%%%
            % Prepare a macro scan for xy axes!
            % Scanning is done by calling 'ScanNextLine'.
            % Aborting via 'AbortScan'.
            % x/y - Vectors with the points to scan, points should have
            % equal distance between them.
            % z - The starting points for the other axis.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.scanRunning
                obj.sendWarning('2D Scan is in progress, previous scan canceled');
                AbortScan(obj);
            end
            Move(obj, 'z', z);
            PrepareScanInTwoDimensions(obj, y, x, nFlat, nOverRun, tPixel, 'y', 'x');
            QueryPos(obj);
        end
        
        function PrepareScanYZ(obj, x, y, z, nFlat, nOverRun, tPixel)
            %%%%%%%%%%%%%% TWO DIMENSIONAL YZ SCAN MACRO %%%%%%%%%%%%%%
            % Prepare a macro scan for yz axes!
            % Scanning is done by calling 'ScanNextLine'.
            % Aborting via 'AbortScan'.
            % y/z - Vectors with the points to scan, points should have
            % equal distance between them.
            % x - The starting points for the other axis.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.scanRunning
                obj.sendWarning('2D Scan is in progress, previous scan canceled');
                AbortScan(obj);
            end
            Move(obj, 'x', x);
            PrepareScanInTwoDimensions(obj, y, z, nFlat, nOverRun, tPixel, 'y', 'z');
            QueryPos(obj);
        end
        
        function PrepareScanZX(obj, x, y, z, nFlat, nOverRun, tPixel)
            %%%%%%%%%%%%%% TWO DIMENSIONAL XZ SCAN MACRO %%%%%%%%%%%%%%
            % Prepare a macro scan for xz axes!
            % Scanning is done by calling 'ScanNextLine'.
            % Aborting via 'AbortScan'.
            % x/z - Vectors with the points to scan, points should have
            % equal distance between them.
            % y - The starting points for the other axis.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.scanRunning
                obj.sendWarning('2D Scan is in progress, previous scan canceled');
                AbortScan(obj);
            end
            Move(obj, 'y', y);
            PrepareScanInTwoDimensions(obj, z, x, nFlat, nOverRun, tPixel, 'z', 'x');
            QueryPos(obj);
        end
        
        function PrepareScanZY(obj, x, y, z, nFlat, nOverRun, tPixel)
            %%%%%%%%%%%%%% TWO DIMENSIONAL YZ SCAN MACRO %%%%%%%%%%%%%%
            % Prepare a macro scan for yz axes!
            % Scanning is done by calling 'ScanNextLine'.
            % Aborting via 'AbortScan'.
            % y/z - Vectors with the points to scan, points should have
            % equal distance between them.
            % x - The starting points for the other axis.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.scanRunning
                obj.sendWarning('2D Scan is in progress, previous scan canceled');
                AbortScan(obj);
            end
            Move(obj, 'x', x);
            PrepareScanInTwoDimensions(obj, z, y, nFlat, nOverRun, tPixel, 'z', 'y');
            QueryPos(obj);
        end
        
        function PrepareScanInTwoDimensions(obj, macroScanAxisVector, normalScanAxisVector, nFlat, nOverRun, tPixel, macroScanAxis, normalScanAxis) 
            %%%%%%%%%%%%%% TWO DIMENSIONAL SCAN MACRO %%%%%%%%%%%%%%
            % Does a macro scan for given axes!
            % scanAxisVector1/2 - Vectors with the points to scan, points
            % should increase with equal distances between them.
            % nFlat - How many flat points should be in the beginning of the scan.
            % nOveRun - How many extra points should be taken from each.
            % tPixel - Scan time for each pixel.
            % scanAxis1/2 - The axes to scan (x,y,z or 1 for x, 2 for y and
            % 3 for z).
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            numberOfMacroPixels = length(macroScanAxisVector);
            numberOfNormalPixels = length(normalScanAxisVector);
            
            % Prepare Scan
            obj.macroPixelTime = tPixel;
            obj.macroNumberOfPixels = numberOfNormalPixels;
            obj.macroScanVector = macroScanAxisVector;
            obj.macroScanAxis = macroScanAxis;
            obj.macroNormalScanAxis = normalScanAxis;
            obj.macroNormalScanVector = normalScanAxisVector;
            obj.macroIndex = 1;
            obj.macroScanStartEndVector = macroScanAxisVector([1, end]);
            
            obj.PrepareScanInOneDimension(macroScanAxisVector, nFlat, nOverRun, tPixel, macroScanAxis);
            obj.MovePrivate([normalScanAxis, macroScanAxis], [normalScanAxisVector(1), macroScanAxisVector(1)]);
            % From here until AbortScan the scan owns the stage: a Move from
            % outside (GUI) cancels it. This flag was never set before, so
            % that check never worked. Scan moves use MovePrivate.
            obj.scanRunning = true;
        end
        
        function [forwards, done] = ScanNextLine(obj)
            % Scans the next line for the 2D scan, to be used after
            % 'PrepareScanXX'.
            % No other commands should be used between 'PrepareScanXX' or
            % until 'AbortScan' has been called.
            % forwards is 1 when the line ran in the order of the scan
            % vector and 0 when it ran in reverse (the caller flips it).
            % scanBidirectional = false: every line forwards.
            % scanBidirectional = true: odd lines forwards, even reversed.
            if ~obj.scanRunning || isempty(obj.macroIndex) || obj.macroIndex < 1 || obj.macroIndex > obj.macroNumberOfPixels
                obj.sendError(sprintf('No scan detected.\nFunction can only be called after ''PrepareScanXX!'''));
            end
            
            % Advance the line counter before scanning, so that if this line
            % fails part-way, PrepareRescanLine steps back to THIS line (as
            % it does when the failure comes later, in readFromScan).
            % Counting at the end made a mid-line failure rescan the line
            % before, which stopped the scan on line 1 and shifted every
            % later line by one.
            line = obj.macroIndex;
            obj.macroIndex = line + 1;
            
            % Scan
            MovePrivate(obj, obj.macroNormalScanAxis, obj.macroNormalScanVector(line));
            forwards = obj.isLineForwards(line);
            if forwards
                ScanOneDimension(obj, obj.macroScanVector, 0, 0, obj.macroPixelTime, obj.macroScanAxis)
            else
                ScanOneDimension(obj, flip(obj.macroScanVector), 0, 0, obj.macroPixelTime, obj.macroScanAxis)
            end
            done = (line == obj.macroNumberOfPixels);
        end
        
        function forwards = isLineForwards(obj, lineIndex)
            % Direction of line number lineIndex (1-based).
            forwards = double(~obj.scanBidirectional || mod(lineIndex, 2) == 1);
        end

        function PrepareRescanLine(obj)
            % Prepares the previous line for rescanning.
            % Scanning is done with "ScanNextLine"
            if ~obj.scanRunning
                obj.sendWarning(sprintf('No scan detected. Function can only be called after ''PrepareScanXX!''.\nThis will work if attempting to rescan a 1D scan, but macro will be rewritten.\n'));
                return
            elseif (obj.macroIndex <= 1)
                obj.sendError('Scan did not start yet. Function can only be called after ''ScanNextLine!''');
            end
            
            % Decrease index. No move is needed: ScanNextLine enters every
            % line from its lead-in point, wherever the stage is now.
            obj.macroIndex = obj.macroIndex - 1;
        end
        
        function AbortScan(obj)
            % Aborts the 2D scan defined by 'PrepareScanXX';
            if obj.scanRunning
                obj.macroIndex = -1;
                obj.scanRunning = 0;
            end
        end
        
         function [tiltEnabled, thetaXZ, thetaYZ] = GetTiltStatus(obj)
            % Return the status of the tilt control.
            tiltEnabled = obj.tiltCorrectionEnable;
            thetaXZ = obj.tiltThetaXZ;
            thetaYZ = obj.tiltThetaYZ;
          end 
        
        function JoystickControl(obj, enable) %#ok<INUSD>
            % Changes the joystick state for all axes to the value of
            % 'enable' - 1 to turn Joystick on, 0 to turn it off.
            obj.sendWarning(sprintf('No joystick support for the %s controller.\n', obj.controllerModel));
        end
        
        function binaryButtonState = ReturnJoystickButtonState(obj)
            % Returns the state of the buttons in 3 bit decimal format.
            % 1 for first button, 2 for second and 4 for the 3rd.
            obj.sendWarning(sprintf('No joystick support for the %s controller.\n', obj.controllerModel));
            binaryButtonState = 0;
        end
        
        function maxScanSize = ReturnMaxScanSize(obj, nDimensions)
            % Returns the maximum number of points allowed for an
            % 'nDimensions' scan.
            maxScanSize = obj.MAX_SCAN_SIZE;
        end
        
        function ChangeLoopMode(obj, mode)
            % Changes between closed and open loop.
            % 'mode' should be either 'Open' or 'Closed'.
            % Stage will auto-lock when in open mode, which should increase
            % stability.
            error('Standa stage cannot chnage loop mode!')
        end
        
        function FastScan(obj, enable)
            % Changes the scan between fast & slow mode
            % 'enable' - 1 for fast scan, 0 for slow scan.
            % StageScanner.startScan calls this first in every scan, so it
            % is also where a scan started INSIDE another one is refused
            % (see refuseNestedScan).
            obj.refuseNestedScan();
            if enable
                obj.sendWarning('No fast scan for this stage!');
            end
            obj.fastScan = 0;
        end

        function refuseNestedScan(obj)
            % MATLAB runs a button click whenever it processes events: in
            % the drawnow while a move is waited for, in the drawnow after
            % each image update, and inside saveas/print while the previous
            % scan's image is being auto-saved. A Scan click handled there
            % starts a second scan INSIDE the first one (or inside its save).
            % The two share the StageScanner object and the stage: the
            % inner scan cancels the outer one, the outer one resumes later
            % at the wrong position ("stopped ... um away"), and the image
            % updates run inside the save's graphics update ("Attempt to
            % modify the tree during an update traversal"), so the picture
            % disappears. (Seen 2026-11-23.) Such a click is refused here,
            % before the second scan moves the stage, clears the image or
            % marks the scanner as scanning.
            % StageScanner.startScan has already copied the GUI scan settings
            % into scanner.mStageScanParams before it calls FastScan; the
            % running scan uses that property for its image updates and its
            % saved file, so its own settings are put back here first.
            stack = dbstack;                % (MATLAB: no dot-indexing on a call without parentheses)
            names = {stack.name};
            isStartScan = ~cellfun(@isempty, regexp(names, 'startScan$', 'once'));
            if sum(isStartScan) == 0; return; end   % not called by a scan start
            scanner = getObjByName(StageScanner.NAME);
            if sum(isStartScan) == 1
                if ~isempty(scanner)
                    obj.paramsOfRunningScan = scanner.mStageScanParams;
                end
                return
            end
            if ~isempty(scanner) && ~isempty(obj.paramsOfRunningScan)
                scanner.mStageScanParams = obj.paramsOfRunningScan;
            end
            msg = ['Scan not started: the previous scan is still running or still saving its image. ', ...
                'Wait until the scan has finished and its image is saved, then press Scan.'];
            obj.sendWarning(msg);
            error('ClassStanda8SMC4_stages:nestedScan', '%s', msg);
        end
        
        function success = EnableTiltCorrection(obj, enable)
            % Enables the tilt correction according to the angles.
            if ~strcmp(obj.VALID_AXES, obj.SCAN_AXES)
                string = BooleanHelper.ifTrueElse(isscalar(obj.VALID_AXES), 'axis', 'axes');
                obj.sendWarning(sprintf('Controller %s has only %s %s, and can''t do tilt correction.', ...
                    obj.controllerModel, obj.VALID_AXES, string));
                success = 0;
                return;
            end
            obj.tiltCorrectionEnable = enable;
            success = 1;
        end
        
        function success = SetTiltAngle(obj, thetaXZ, thetaYZ)
            % Sets the tilt angles between Z axis and XY axes.
            % Angles should be in degrees, valid angles are between -5 and 5
            % degrees.
            if (thetaXZ < -5 || thetaXZ > 5) || (thetaYZ < -5 || thetaYZ > 5)
                obj.sendWarning(sprintf('Angles are outside the limits (-5 to 5 degrees).\nAngles were not set.'));
                success = 0;
                return
            end
            obj.tiltThetaXZ = thetaXZ;
            obj.tiltThetaYZ = thetaYZ;
            success = 1;
        end
        
        function testErr(obj, err)
            if err ~= 0
                error('error %d was found on connection to standa stage', err)                
            end            
        end
        
    end
end