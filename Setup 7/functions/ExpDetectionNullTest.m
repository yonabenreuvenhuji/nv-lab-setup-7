classdef ExpDetectionNullTest < Experiment
    %EXPDETECTIONNULLTEST Laser on all the time, detection windows opened and closed. No MW.
    %
    % Every repeat is the CW-ESR timing with the microwave taken out:
    %
    %   [0.1 us] [window 1: D - 0.1 us, detector] [gap] [0.1 us] [window 2: D - 0.1 us, detector] [gap]
    %
    % with the green laser ON in every event, so the light never changes. Nothing on the
    % sample changes either: no AWG, no frequency generator, no trigger. Whatever differs
    % between window 1 and window 2, or inside a window, is the detection chain itself --
    % the detector gate, the digitizer arming and capturing a segment, the photodiode
    % reference normalisation, or electrical pickup from channels switched during window 1.
    %
    % toggleChannels lists pulse-generator channels to switch ON during window 1 only,
    % exactly where the ESR switches its MW switches. {'MW','MW2'} with the AWG stopped
    % toggles the switches with no microwave behind them: a test for TTL pickup.
    %
    % The digitizer records each window as a full trace. The reduction here is done
    % directly on the raw segments -- NOT through Experiment.processData, whose photodiode
    % normalisation would divide neighbouring samples of a full trace.
    %
    % Usually run through detectionNullTest, which runs all the conditions and reports.
    % On its own:
    %     e = ExpDetectionNullTest; e.gap = 10; e.run; r = e.analyse();

    properties (Constant)
        NAME = 'Detection_NullTest';
    end

    properties
        gap = 10;               % us between the windows, = ESR delayBetweenDetections
        toggleChannels = {};    % channels switched ON during window 1 only (no MW behind them)
        dropRepeats = [];       % repeats discarded at the start of each pg.run, while the
                                % sample settles after the laser comes on. [] = auto (>= 500 us).
    end

    properties (SetAccess = protected)
        perAverage = struct('ratio', {}, 'ratioSterr', {}, 'level1', {}, 'level2', {}, ...
                            'start1', {}, 'start2', {});
                                % one entry per completed average, see reduce()
    end

    properties (Dependent)
        time                    % us, time inside a window
    end

    properties (Constant, Hidden)
        SAMPLE_RATE = 100;      % MS/s, the digitizer's rate in experiments
        TRIG = 0.1;             % us, the leading event of each window (ESR: TRIG_DURATION)
        START_WINDOW = 1;       % us, "start of window" for the transient check
    end

    %%
    methods
        function obj = ExpDetectionNullTest()
            obj@Experiment(ExpDetectionNullTest.NAME);
            obj.parameterName = 'time';
            obj.repeats = 200;
            obj.averages = 20;
            obj.isTracking = 0;
            obj.trackThreshhold = 0.7;
            obj.detectionDuration = 125;            % the ESR window
            obj.referenceDetectionDuration = 125;
            obj.digitizerFullDataAcquisition = true;
            obj.isPlotAlternateAvailable = true;
            obj.displayType1 = 'Normal';
            obj.displayType2 = 'Alt (One Color)';
            obj.mCurrentXAxisParam = ExpParamDoubleVector('Time in window', [], [], 'us', obj.NAME);
            obj.signalParam  = ExpResultDoubleVector('FL', [], [], 'V', obj.NAME, 'Window 1');
            obj.signalParam2 = ExpResultDoubleVector('FL', [], [], 'V', obj.NAME, 'Window 2');
        end

        function set.gap(obj, v)
            if ~(isnumeric(v) && isscalar(v) && v >= 0.5)
                obj.sendError('gap must be a number of us, at least 0.5')
            end
            obj.gap = v;
            obj.changeFlag = true;
        end

        function set.toggleChannels(obj, v)
            if ischar(v); v = {v}; end
            if ~iscell(v)
                obj.sendError('toggleChannels must be a cell array of channel names, e.g. {''MW'',''MW2''}')
            end
            if any(cellfun(@(c) any(strcmpi(c, {'greenLaser', 'detector'})), v))
                obj.sendError('toggleChannels must not include greenLaser or detector')
            end
            obj.toggleChannels = v(:)';
            obj.changeFlag = true;
        end

        function t = get.time(obj)
            dt = 1 / obj.SAMPLE_RATE;
            t = 0 : dt : obj.detectionDuration - dt;
        end

        function n = getTotalNumberOfParams(obj)
            n = length(obj.time);
        end

        function r = analyse(obj, quiet)
            % Summary over all completed averages. Returns a struct; prints unless quiet.
            if ~exist('quiet', 'var'); quiet = false; end
            % Only the averages of THIS run: reset() zeroes signal but not perAverage, so
            % a re-run of the same object can leave older entries beyond currIter.
            K = min(obj.currIter, numel(obj.perAverage));
            if K == 0
                obj.sendError('No completed averages to analyse.');
            end
            P = obj.perAverage(1:K);
            ratios = [P.ratio];
            r.ratio      = mean(ratios);
            if K > 1
                r.ratioSterr = std(ratios) / sqrt(K);
            else
                r.ratioSterr = P(1).ratioSterr;
            end
            r.offset     = r.ratio - 1;                 % window 1 relative to window 2
            r.start1     = mean([P.start1]);   % (first 1 us)/(rest) - 1, window 1
            r.start2     = mean([P.start2]);
            r.level1     = mean([P.level1]);
            r.level2     = mean([P.level2]);
            r.nAverages  = K;
            r.gap        = obj.gap;
            r.window     = obj.detectionDuration;
            r.toggled    = obj.toggleChannels;
            S1 = squeeze(obj.signal(1, :, 1:K)); S2 = squeeze(obj.signal(2, :, 1:K));
            if K > 1; S1 = mean(S1, 2); S2 = mean(S2, 2); end
            r.t     = obj.time(:);
            r.trace1 = S1(:);
            r.trace2 = S2(:);
            if ~quiet
                fprintf(['\nwindow 1 / window 2 = %.6f +/- %.6f  (offset %+.2e)   start transients: ' ...
                         'w1 %+.2e, w2 %+.2e   [gap %g us, toggled: %s, %d averages]\n'], ...
                        r.ratio, r.ratioSterr, r.offset, r.start1, r.start2, r.gap, ...
                        strjoin([{''}, r.toggled], ' '), K);
            end
        end
    end

    %% Overridden from Experiment
    methods (Access = protected)
        function prepare(obj)
            obj.referenceDetectionDuration = obj.detectionDuration;
            obj.detectionPeriodsPerRepeat = 2;
            obj.runsPerPerform = 1;
            obj.digitizerFullDataAcquisition = true;
            D = obj.detectionDuration;
            T = obj.TRIG;
            extra = obj.toggleChannels;

            % The ESR sequence (AWGExpESR.prepare) with MW and trigger removed; any
            % toggleChannels sit exactly where the ESR has its MW switches.
            S = Sequence;
            S.addEvent(T,       [{'greenLaser'}, extra]);
            S.addEvent(D - T,   [{'greenLaser'}, extra, {'detector'}]);
            S.addEvent(obj.gap, {'greenLaser'});
            S.addEvent(T,       {'greenLaser'});
            S.addEvent(D - T,   {'greenLaser', 'detector'});
            S.addEvent(obj.gap, {'greenLaser'});
            obj.prepareInternal(S);

            obj.mCurrentXAxisParam.value = obj.time;
        end

        function perform(obj)
            pg = getObjByName(PulseGenerator.NAME);
            spcm = getObjByName(Spcm.NAME);
            success = false;
            for trial = 1:5
                if obj.checkEmergencyStop(); return; end
                try
                    data = obj.getRawData(pg, spcm);
                    [w1, w2, e1, e2, res] = obj.reduce(data);
                    obj.signal(1, :, obj.currIter) = w1;
                    obj.signal(2, :, obj.currIter) = w2;
                    obj.sterr(1, :, obj.currIter) = e1;
                    obj.sterr(2, :, obj.currIter) = e2;
                    obj.perAverage(obj.currIter) = res;
                    success = true;
                    break
                catch err
                    err2warning(err);
                    fprintf('Experiment failed at trial %d, attempting again.\n', trial);
                    try; spcm.stopExperimentCount; catch; end
                end
            end
            if ~success
                obj.emergencyStop;
                return
            end

            % Plot the running mean of both windows.
            K = obj.currIter;
            S1 = reshape(obj.signal(1, :, 1:K), [], K);
            S2 = reshape(obj.signal(2, :, 1:K), [], K);
            obj.signalParam.value  = mean(S1, 2)';
            obj.signalParam2.value = mean(S2, 2)';
            obj.signalParam.sterr  = (std(S1, 0, 2) / sqrt(max(K-1, 1)))';
            obj.signalParam2.sterr = (std(S2, 0, 2) / sqrt(max(K-1, 1)))';
            obj.currParamIter = obj.getTotalNumberOfParams;
            sendEventParamIterationDone(obj);
        end

        function [w1, w2, e1, e2, res] = reduce(obj, data)
            % data: [2*repeats x nSamples] raw segments, in acquisition order, so odd
            % rows are window 1 and even rows window 2 of the same repeat.
            nSamp = length(obj.time);
            if size(data, 2) ~= nSamp && size(data, 1) == nSamp
                data = data';
            end
            if size(data, 1) ~= 2 * obj.repeats
                error('ExpDetectionNullTest: got %d segments, expected %d', size(data, 1), 2 * obj.repeats);
            end
            if size(data, 2) > nSamp
                % The digitizer may add samples for its channel delay; keep the window length.
                data = data(:, 1:nSamp);
            elseif size(data, 2) < nSamp
                error('ExpDetectionNullTest: segments have %d samples, expected %d', size(data, 2), nSamp);
            end
            W1 = data(1:2:end, :);
            W2 = data(2:2:end, :);
            drop = obj.dropRepeats;
            if isempty(drop)
                repeatUs = 2 * (obj.detectionDuration + obj.gap);
                drop = max(1, ceil(500 / repeatUs));
            end
            drop = min(drop, obj.repeats - 2);
            W1 = W1(drop+1:end, :);
            W2 = W2(drop+1:end, :);
            n = size(W1, 1);

            w1 = mean(W1, 1);
            w2 = mean(W2, 1);
            e1 = std(W1, 0, 1) / sqrt(n);
            e2 = std(W2, 0, 1) / sqrt(n);

            % One ratio per repeat: the two windows are 135 us apart, so laser drift cancels.
            m1 = mean(W1, 2);
            m2 = mean(W2, 2);
            rr = m1 ./ m2;

            % Start of window against the rest of it.
            t = obj.time;
            st = t < obj.START_WINDOW;
            rest = t >= 10 * obj.START_WINDOW;
            res = struct('ratio', mean(rr), 'ratioSterr', std(rr) / sqrt(n), ...
                         'level1', mean(m1), 'level2', mean(m2), ...
                         'start1', mean(w1(st)) / mean(w1(rest)) - 1, ...
                         'start2', mean(w2(st)) / mean(w2(rest)) - 1);
        end

        function wrapUp(obj)
            obj.wrapUpInternal();
            if ~isempty(obj.perAverage)
                obj.analyse();
            end
        end

        function params = alternateSignal(obj)
            % Alternate view: (window 1 - window 2) / window 2, point by point.
            if obj.currIter == 0 || isempty(obj.signalParam.value)
                params = ExpResultDoubleVector('Difference', [], [], 'relative', obj.NAME);
                return
            end
            d = (obj.signalParam.value - obj.signalParam2.value) ./ mean(obj.signalParam2.value);
            params = ExpResultDoubleVector('Difference', d, [], 'relative', obj.NAME, 'W1 - W2');
        end
    end
end
