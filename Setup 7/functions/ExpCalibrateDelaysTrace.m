classdef ExpCalibrateDelaysTrace < ExpAWGDigitizerCalibrateDelays
    %EXPCALIBRATEDELAYSTRACE The level-1 delay staircase, averaged as a plain trace.
    %
    % Identical to ExpAWGDigitizerCalibrateDelays -- same sequence, same AWG load, same
    % perform -- except for ONE step: how the raw data of one pg.run becomes a trace.
    % Used by calibChannelDelays; there is no reason to run it by hand.
    %
    % WHY IT EXISTS
    %
    % The inherited Experiment.processData was written for experiments that read one
    % number per detection window. On a full trace (digitizerFullDataAcquisition) it does
    % two things that break this calibration:
    %
    %  1. The glitch filter. expProperties.json sets checkGlitchInData = true, and on a
    %     photodiode that throws away a WHOLE repeat if any ONE of its samples is more
    %     than 4 sigma from the mean of that sample. A trace has thousands of samples, so
    %     the chance that a repeat survives falls exponentially with the trace length:
    %     52 us of trace (channelDelay 10 us, three channels) is 5200 chances per repeat.
    %     Sooner or later an average comes back with every repeat thrown away, i.e. all
    %     NaN, and the running mean across averages in perform() does not skip NaN -- so
    %     from that average on the displayed trace is NaN everywhere, the GUI draws
    %     nothing, and at the end summaryData.y is the plot's placeholder instead of the
    %     trace. That is the "no data at channelDelay > 5 us" failure and the
    %     out-of-bounds error at calibChannelDelays line 223.
    %
    %  2. The photodiode normalisation divides odd columns by even ones. On a trace the
    %     columns are SAMPLES, so every odd sample is divided by its 10 ns neighbour.
    %     Harmless on a flat stretch, noisy and slightly biased on the edges and in the
    %     dark -- which is exactly where the delays are read.
    %
    % Here: the mean over repeats (first repeat dropped, as processData does), skipping
    % NaN samples, and nothing else. A bad repeat costs one repeat's worth of averaging
    % instead of the whole run.

    methods (Access = protected)
        function [signal, sterr] = processData(obj, rawData)
            n = obj.repeats;
            m = numel(rawData) / n;
            if m ~= round(m)
                error('ExpCalibrateDelaysTrace: %d samples do not split into %d repeats.', ...
                      numel(rawData), n);
            end
            s = reshape(rawData, m, n)';            % one row per repeat
            if n > 10
                s = s(2:end, :);                    % first repeat follows no initialisation
            end

            ok = ~isnan(s);
            k  = sum(ok, 1);
            s0 = s;  s0(~ok) = 0;
            signal = sum(s0, 1) ./ k;               % NaN only where every repeat was NaN
            d = s - signal;  d(~ok) = 0;
            sterr = sqrt(sum(d.^2, 1) ./ max(k - 1, 1)) ./ sqrt(max(k, 1));

            % The class preallocated getTotalNumberOfParams samples per average. Match it,
            % so a one-sample rounding difference cannot fail the assignment in perform.
            M = obj.getTotalNumberOfParams;
            if m > M
                signal = signal(1:M);  sterr = sterr(1:M);
            elseif m < M
                signal(end+1:M) = NaN; sterr(end+1:M) = NaN;
            end
        end
    end
end
