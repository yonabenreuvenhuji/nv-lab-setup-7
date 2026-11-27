classdef powerSupply_dummy < powerSupply
    %POWERSUPPLY_DUMMY A power supply with no hardware behind it.
    %
    % Behaves like a single-channel E3632A (15 V / 7 A and 30 V / 4 A ranges) and remembers
    % what it is told, so the Helmholtz class and its GUI can be built and exercised offline
    % -- testHelmholtzOffline does exactly that. Use "type": "dummy" in a setup json.

    methods
        function obj = powerSupply_dummy(name, address)
            obj@powerSupply(name, address);
            obj.channels = [];
            obj.ranges = ["P15V", "P30V"];
            obj.minCurrent = [0, 0];
            obj.maxCurrent = [7, 4];
            obj.minVoltage = [0, 0];
            obj.maxVoltage = [15, 30];
            obj.current = 0;
            obj.voltage = 0;
            obj.range = "P15V";
            obj.mode = "current";
        end

        function initilize(obj, address) %#ok<INUSD>
            obj.id = 'dummy';
        end

        function connect(obj) %#ok<MANU>
        end

        function closeConnection(obj)
            obj.id = '';
        end

        function reset(obj)
            obj.current = 0; obj.voltage = 0;
        end

        function setCurrent(obj, I, channel)
            if ~exist('channel', 'var'); channel = []; end
            obj.checkCurrentLim(I, channel);
            obj.current = I;
        end

        function I = getCurrent(obj, channel) %#ok<INUSD>
            I = obj.current;
        end

        function setVoltage(obj, V, channel)
            if ~exist('channel', 'var'); channel = []; end
            obj.checkVoltageLim(V, channel);
            obj.voltage = V;
        end

        function V = getVoltage(obj, channel) %#ok<INUSD>
            V = obj.voltage;
        end

        function setMode(obj, mode, channel) %#ok<INUSD>
            obj.checkMode(mode);
            obj.mode = mode;
        end

        function m = getMode(obj, channel) %#ok<INUSD>
            m = obj.mode;
        end

        function setRange(obj, range, channel)
            if ~exist('channel', 'var'); channel = []; end
            obj.checkRange(range, channel);
            obj.range(1) = range;
        end

        function r = getRange(obj, channel) %#ok<INUSD>
            r = obj.range(1);
        end

        function output(obj, state, channel) %#ok<INUSD>
            obj.checkOutput(state);
        end

        function o = getOutput(obj, channel) %#ok<INUSD,STOUT>
            o = 0;
        end

        function answer = sendCommand(obj, what) %#ok<INUSD>
            answer = '';
        end

        function checkError(obj) %#ok<MANU>
        end
    end

    methods (Static)
        function powerSupplyObj = create(powerSupplyName, powerSupplyStruct)
            powerSupplyObj = powerSupply_dummy(powerSupplyName, powerSupplyStruct.address);
        end
    end
end
