classdef SpcmPhotodiodeDummy < SpcmDummy
    %SPCMPHOTODIODEDUMMY A photodiode on no hardware, for running the GUI offline.
    %
    %   "spcm": {"classname": "photodiodeDummy", "Switch": {...}}
    %
    % Behaves like PhotoDiodeDigitizerNiDaqControlled as far as the GUI can tell: it has a
    % photodiode, two inputs (signal, reference) and the same six input options
    % (currentChannelGUI). The signal is a few bright spots on a dim background, read at the
    % stage's position, so a scan of a dummy stage draws a recognisable image, and 'signal &
    % reference' gives two values per reading. Values are in volts (tens of mV), like the
    % real photodiode.
    %
    % Scans: the stage calls acquireScanPixel(tPixel) at every pixel (ClassDummyStage and
    % the Standa stage both do), and readFromScan returns the line.

    properties
        currentChannelGUI = 1   % as in PhotoDiodeDigitizerNiDaqControlled.PHOTODIODE_OPTIONS
        noise = 0.02            % relative noise of one reading
        pixelPause = 0          % s per scan pixel (0 = as fast as possible)
    end

    properties (Access = private)
        lineValues = []
        nLine = 0
    end

    methods
        function obj = SpcmPhotodiodeDummy
            obj@SpcmDummy;
            obj.availableProperties.(obj.HAS_PHOTODIODE) = true;
            obj.availableProperties.(obj.HAS_DUAL_INPUT) = true;
        end

        %%% By time %%%
        function [value, sterr] = readFromTime(obj, ~)
            if obj.integrationTime <= 0
                obj.sendError('Can''t call readFromTime() without calling obj.prepareReadByTime() first!');
            end
            pause(obj.integrationTime)
            [sig, ref] = obj.signalHere();
            [value, sterr] = obj.reduce(sig, ref);
        end

        %%% By stage %%%
        function prepareCountByStage(obj, ~, nPixels, ~, ~, ~)
            if ~ValidationHelper.isValuePositiveInteger(nPixels)
                obj.sendError(sprintf('Can''t prepare for reading %d points, only positive integers allowed!', nPixels));
            end
            if obj.currentChannelGUI == 3
                obj.sendError('The ''both inputs'' photodiode option gives two values per pixel and cannot be scanned. Choose another option.');
            end
            obj.timesToRead = nPixels;
            obj.lineValues = [];
            obj.nLine = 0;
        end

        function startScanCount(obj)
            obj.calledScanStart = true;
            obj.lineValues = [];
            obj.nLine = 0;
        end

        function acquireScanPixel(obj, tPixel)
            if obj.pixelPause > 0
                pause(obj.pixelPause);
            elseif exist('tPixel', 'var') && ~isempty(tPixel)
                pause(min(tPixel, 0.002));
            end
            [sig, ref] = obj.signalHere();
            v = obj.reduce(sig, ref);
            obj.nLine = obj.nLine + 1;
            obj.lineValues(obj.nLine) = v(1);
        end

        function [voltage, sterr] = readFromScan(obj)
            if obj.nLine ~= obj.timesToRead
                n = obj.nLine;
                obj.nLine = 0;
                error('Scan line incomplete: %d of %d pixels were measured.', n, obj.timesToRead);
            end
            voltage = obj.lineValues;
            sterr = obj.noise * abs(voltage);
            obj.nLine = 0;
        end

        function clearScanRead(obj)
            obj.timesToRead = 0;
            obj.nLine = 0;
        end

        function values = partialScanLine(obj)
            % The pixels of the current line measured so far (drawn live)
            values = obj.lineValues(1:min(obj.nLine, numel(obj.lineValues)));
        end
    end

    methods (Access = private)
        function [sig, ref] = signalHere(obj)
            % Fluorescence at the dummy stage's position: three spots (FWHM ~ 6 um) and a
            % slow background, in volts.
            pos = [0 0 0];
            try
                stages = ClassStage.getStages;
                if ~isempty(stages)
                    pos = stages{1}.Pos('xyz');
                end
            catch
            end
            spots = [ 12  -8   0; -20  15   5;  30  25 -10];
            amp   = [0.060 0.035 0.045];
            r2 = sum((spots - pos).^2 .* [1 1 0.15], 2)';
            sig = 0.012 + 0.002*sin(pos(1)/40) + sum(amp .* exp(-r2/(2*2.5^2)));
            ref = 0.080 * (1 + 0.01*sin(2*pi*rem(posixtime(datetime('now')), 7)/7));
            sig = sig * (1 + obj.noise*randn);
            ref = ref * (1 + 0.5*obj.noise*randn);
        end

        function [value, sterr] = reduce(obj, sig, ref)
            switch obj.currentChannelGUI
                case 1; value = sig;
                case 2; value = ref;
                case 3; value = [sig; ref];
                case 4; value = sig - ref;
                case 5; value = sig / (ref / 0.080);
                case 6; value = 0.0005 * (1 + obj.noise*randn);
                otherwise; error('Unknown photodiode option %d', obj.currentChannelGUI);
            end
            sterr = obj.noise * abs(value) / 3;
        end
    end
end
