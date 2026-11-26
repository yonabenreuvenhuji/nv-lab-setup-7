classdef FrequencyGeneratorSRS < FrequencyGenerator
    %FREQUENCYGENERATORSRS SRS frequency genarator class
    
    properties (Constant, Hidden)
        TYPE = 'srs';
%         NAME = 'srsFrequencyGenerator';
        
        NEEDED_FIELDS = {'address', 'port', 'serialNumber', 'minFrequency', 'maxFrequency', 'minAmplitude', 'maxAmplitude'}
        OPTIONAL_FIELDS = {'keepOn'};
    end
       
    properties (Access = private)
        t           % tcpclient object; empty while disconnected
        tcpAddress  % char. IP address of the SRS
        tcpPort     % double. TCP port of the SRS
    end

    methods (Access = private)
        function obj = FrequencyGeneratorSRS(name, address, port, frequencyLimits, amplitudeLimits, keepOn)
            % All models are the same in regards to controlling them, but
            % the limitations on the amplitude and on the allowed frequencies may vary.
            obj@FrequencyGenerator(name, frequencyLimits, amplitudeLimits, keepOn);
            % tcpip() and its fopen/fprintf/fscanf throw an error from MATLAB R2026b,
            % so this talks through tcpclient (R2014b+). A tcpclient connects as soon
            % as it is created, so it is created in connect(), keeping the old
            % open-on-first-use / close-in-disconnect timing.
            obj.tcpAddress = address;
            obj.tcpPort = port;

            obj.initialize;
        end
    end

    methods
        function connect(obj)
            if isempty(obj.t) || ~isvalid(obj.t)
                % Defaults match the old tcpip object: Terminator "LF", Timeout 10 s.
                obj.t = tcpclient(obj.tcpAddress, obj.tcpPort);
            end
        end

        function disconnect(obj)
            if ~isempty(obj.t) && isvalid(obj.t)
                delete(obj.t);      % deleting a tcpclient closes the connection
            end
            obj.t = [];
        end

        function delete(obj)
            obj.disconnect;
        end

        function sendCommand(obj, command)
            % Actually sends command to hardware
            obj.connect() % Only happens if it is closed
            writeline(obj.t, command);      % appends the "LF" terminator, like fprintf did
        end

        function value = readOutput(obj, what)
            % Get value returned from object. readline drops the terminator;
            % strtrim also drops the CR the SRS sends, as fscanf('%s') did.
            value = strtrim(char(readline(obj.t)));
            switch what
                case {'frequency', 'freq', 'f'}
                    value = num2str(str2double(value)/1e6);  % convert Hz to MHz
            end
        end
    end
    
    %%
    methods (Static)
        function obj = getInstance(struct)
            missingField = FactoryHelper.usualChecks(struct, ...
                FrequencyGeneratorSRS.NEEDED_FIELDS);
            if ~isnan(missingField)
                EventStation.anonymousError(...
                    'Trying to create an SRS frequency generator, encountered missing field - "%s". Aborting',...
                    missingField);
            end
            struct = FactoryHelper.supplementStruct(struct, FrequencyGeneratorSRS.OPTIONAL_FIELDS);
            
            name = [FrequencyGeneratorSRS.TYPE, '-', struct.serialNumber];
            frequencyLimits = [struct.minFrequency, struct.maxFrequency];
            amplitudeLimits = [struct.minAmplitude, struct.maxAmplitude];
            keepOn = struct.keepOn;
            obj = FrequencyGeneratorSRS(name, struct.address, struct.port, frequencyLimits, amplitudeLimits, keepOn);

            addBaseObject(obj);
        end
        
        function command = createCommand(what, value, channel) %#ok<INUSD>
           switch lower(what)
               case {'enableoutput', 'output', 'enabled', 'enable'}
                   name = 'ENBR';
               case {'frequency', 'freq'}
                   name = 'FREQ';
                   if ~strcmp(value, '?')
                       value = value*1e6; % convert MHz to Hz
                   end
               case {'amplitude', 'ampl', 'amp'}
                   name = 'AMPR';
               case {'phase'}
                   name = 'PHAS';   % phase is in deg
               otherwise
                   error('Unknown command type %s', what)
           end
           
           if isnumeric(value)
               value = num2str(value);
           end
           command = [name, value];
       end
    end
end