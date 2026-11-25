classdef AWGExpMwQuadrature < AWGExpRabi
    %AWGEXPMWQUADRATURE The two MW channels' balance, measured with BOTH channels on.
    %
    % Used by calibMwQuadrature (level 4, cell 4.2b). Not a general experiment.
    %
    % THE IDEA. Drive with both channels at once, but 90 degrees apart instead of in phase.
    % The rotation axis of that pulse is the direction of Omega1 + i*Omega2, i.e. at
    % atan(Omega2/Omega1): exactly 45 degrees when the two channels deliver the same field
    % AT THE NVs, while both are playing. A second pulse, the in-phase pair as a pi/2 about
    % 135 or -45 degrees (perpendicular to 45), turns a tilt of that axis into population:
    % the two readouts differ by ~sin(tilt). Doing the same with channel 1 leading instead of
    % channel 2 tilts the axis the other way, and a channelPhase error moves both the same
    % way, so
    %
    %       Q = [A(135) - A(-45)] - [B(135) - B(-45)]     A: ch2 leads, B: ch1 leads
    %
    % is ~ +/- 2*sin(imbalance/2) with the phase error cancelled to first order. The drive is
    % strong throughout -- unlike the antiphase null, the spins never sit under a weak
    % residual drive, so the 14N hyperfine detuning does not wash it out.
    %
    % Every configuration is played at several channel-2 levels (quadOffsetsDb, relative to
    % the amplitude set), so the sign convention, the slope and the contrast scale are all
    % MEASURED: the balance is where the fitted Q(offset) line crosses zero, in dB.
    %
    % Two references close the list: 'bright' (no pulse) and 'dark' (an in-phase pi), for
    % the contrast scale.
    %
    % IMPLEMENTATION. A subclass of AWGExpRabi, so the pulse-generator sequence, walk mode,
    % the load cache, the readout and processData are the Rabi's, unchanged and proven. Only
    % the waveforms are different: one per configuration and channel, [zero pad, pulse], the
    % same block structure as a short-tau Rabi point. Every "tau" is the same MW gate length
    % (the longest configuration), which is what walk mode needs anyway; the AWG waveform
    % defines the pulse, as in the Rabi.
    %
    % It inherits the Rabi's NAME ('AWG Rabi'), so constructing one replaces any Rabi object
    % in the workspace, exactly as calibMwPair's own Rabi runs do.
    %
    % WAVEFORM. Channel c plays, for t from 0 after its zero pad,
    %       A(c) * k(c) * cos(w*t - phase(c) - extra)
    % with phase(c) from channelsDifference + channelPhase (the Rabi's own convention, so
    % channelPhase means the same here as in calibMwPair), k(2) = 10^(offset/20), k(1) = 1,
    % and extra = 90 deg on the leading channel during the first pulse, refPhaseDeg on both
    % during the second.

    properties
        quadT1 = 0.1            % us, the quadrature pulse (pi/2 for |Omega1 + i Omega2|)
        quadT2 = 0.07           % us, the in-phase pi/2 that reads the axis; 'dark' is 2*quadT2
        quadOffsetsDb = [-1 0 1]  % dB on channel 2, relative to its set amplitude
    end

    properties (SetAccess = protected)
        configs = []            % one struct per point: name, kind, offsetDb, lead, refPhaseDeg
    end

    methods
        function obj = AWGExpMwQuadrature()
            obj@AWGExpRabi();
            obj.shouldAutosave = false;
            obj.useSequenceWalk = true;
            obj.rebuildConfigs();
        end

        function set.quadT1(obj, v)
            if ~(isnumeric(v) && isscalar(v) && v > 0 && v < 10)
                error('AWGExpMwQuadrature: quadT1 must be a duration in us, 0 < quadT1 < 10.');
            end
            obj.quadT1 = v;
            obj.changeFlag = true;
            obj.rebuildConfigs();
        end

        function set.quadT2(obj, v)
            if ~(isnumeric(v) && isscalar(v) && v > 0 && v < 10)
                error('AWGExpMwQuadrature: quadT2 must be a duration in us, 0 < quadT2 < 10.');
            end
            obj.quadT2 = v;
            obj.changeFlag = true;
            obj.rebuildConfigs();
        end

        function set.quadOffsetsDb(obj, v)
            if ~(isnumeric(v) && isvector(v) && numel(v) >= 2 && all(isfinite(v)) && all(abs(v) <= 6))
                error('AWGExpMwQuadrature: quadOffsetsDb must be >= 2 offsets, each within +/-6 dB.');
            end
            obj.quadOffsetsDb = v(:)';
            obj.changeFlag = true;
            obj.rebuildConfigs();
        end

        function rebuildConfigs(obj)
            % The configuration list and the MW gate ("tau") that goes with it. tau is only
            % re-assigned when it changes, so a prepare on an unchanged object does not raise
            % changeFlag (and with it the framework's "parameters changed" question).
            c = AWGExpMwQuadrature.configList(obj.quadOffsetsDb);
            gate = max(obj.quadT1 + obj.quadT2, 2 * obj.quadT2);
            gate = ceil(gate * 1e3 - 1e-6) / 1e3;          % whole PulseStreamer ns, never shorter
            newTau = gate * ones(1, numel(c));
            if ~isequal(obj.configs, c)
                obj.configs = c;
            end
            if ~isequal(obj.tau, newTau)
                obj.tau = newTau;
            end
        end
    end

    methods (Static)
        function c = configList(offsets)
            % Fixed order: for each offset A(135) A(-45) B(135) B(-45), then bright, dark.
            item = @AWGExpMwQuadrature.cfgItem;
            c = item('', '', 0, 0, 0);
            c = c([]);
            for d = offsets(:)'
                c(end+1) = item(sprintf('A+ %+.2f', d), 'quad', d, 2, 135);   %#ok<AGROW>
                c(end+1) = item(sprintf('A- %+.2f', d), 'quad', d, 2, -45);   %#ok<AGROW>
                c(end+1) = item(sprintf('B+ %+.2f', d), 'quad', d, 1, 135);   %#ok<AGROW>
                c(end+1) = item(sprintf('B- %+.2f', d), 'quad', d, 1, -45);   %#ok<AGROW>
            end
            c(end+1) = item('bright', 'bright', 0, 0, 0);
            c(end+1) = item('dark', 'dark', 0, 0, 0);
        end

        function s = cfgItem(name, kind, d, lead, ref)
            s = struct('name', name, 'kind', kind, 'offsetDb', d, 'lead', lead, ...
                       'refPhaseDeg', ref);
        end

        function segs = segmentsFor(cfg, chan, T1, T2)
            % Rows [tStart tEnd amplitudeFactor extraPhaseRad] for one channel of one point.
            k = 1;
            if chan == 2; k = 10^(cfg.offsetDb / 20); end
            switch cfg.kind
                case 'quad'
                    lead = (chan == cfg.lead) * pi/2;
                    ref = deg2rad(cfg.refPhaseDeg);
                    segs = [0, T1, k, lead; T1, T1 + T2, k, ref];
                case 'dark'
                    segs = [0, 2 * T2, k, 0];
                case 'bright'
                    segs = zeros(0, 4);
                otherwise
                    error('AWGExpMwQuadrature: unknown configuration kind ''%s''.', cfg.kind);
            end
        end

        function y = quadWave(t, segs, Ach, w, ph)
            % The waveform of one channel: Ach * k * cos(w t - ph - extra) inside each segment.
            y = zeros(size(t));
            for s = 1:size(segs, 1)
                m = t >= segs(s, 1) & t < segs(s, 2);
                y(m) = Ach * segs(s, 3) * cos(w * t(m) - ph - segs(s, 4));
            end
        end
    end

    methods (Access = protected)
        function prepare(obj)
            obj.rebuildConfigs();
            prepare@AWGExpRabi(obj);
            % The x axis is the configuration index; every "tau" is the same gate length.
            obj.mCurrentXAxisParam.value = 1:numel(obj.configs);
        end

        function [fp, parts] = awgLoadFingerprint(obj)
            % The Rabi's parts plus what makes these waveforms different from a Rabi's with
            % the same (constant) tau list -- without them a Rabi load could be taken for
            % this one, or one offset list for another.
            [~, parts] = awgLoadFingerprint@AWGExpRabi(obj);
            parts.className     = class(obj);
            parts.quadT1        = obj.quadT1;
            parts.quadT2        = obj.quadT2;
            parts.quadOffsetsDb = obj.quadOffsetsDb;
            parts.configNames   = strjoin({obj.configs.name}, ',');
            fp = DataHash(parts);
        end

        function LoadAWG(obj)
            [hit, fpAtEntry] = obj.tryLoadFromCache();
            if hit; return; end
            obj.walkParts = [];

            % ---- waveforms: same set-up steps as the Rabi's LoadAWG ----------
            baseWFtime = obj.RoundDurationByFrequencies(0.2, obj.frequency);
            obj.AWG.ResetStoredData;
            obj.AWG.workChannels = obj.activeChannels;
            obj.AWG.setAmplitude(obj.amplitude);
            obj.AWG.ChangeSampleRate(baseWFtime);
            minWaveformTime = obj.AWG.minWaveformDuration;
            if obj.nChannels ~= 2
                error('AWGExpMwQuadrature: needs nChannels = 2 (it measures the pair).');
            end
            A = obj.AWG.MaxNormPowerAllowed;
            if numel(A) < 2; A = A(1) * [1 1]; end
            kMax = 10^(max(obj.quadOffsetsDb) / 20);
            if A(2) * kMax > 1
                error(['AWGExpMwQuadrature: channel 2 at +%g dB would exceed the AWG''s full ' ...
                       'scale (%.3f x %.3f > 1). Lower the amplitude or the offsets.'], ...
                      max(obj.quadOffsetsDb), A(2), kMax);
            end
            w = 2*pi*obj.frequency;
            [prePulseTime, phase] = obj.channelsDifference;
            phase = phase + obj.channelPhaseRadians();

            N = numel(obj.configs);
            initZeroIndex = zeros(1, 2);
            wfIndex = zeros(2, N);
            for chan = 1:2
                initZeroIndex(chan) = obj.AWG.AddWaveformByFunction(@(t) 0 * t, prePulseTime(chan), 0);
            end
            fprintf('Quadrature load: %d configurations x 2 channels...\n', N);
            for i = 1:N
                for chan = 1:2
                    segs = AWGExpMwQuadrature.segmentsFor(obj.configs(i), chan, obj.quadT1, obj.quadT2);
                    dur = 0;
                    if ~isempty(segs); dur = max(segs(:, 2)); end
                    Ac = A(chan);  ph = phase(chan);
                    f = @(t) AWGExpMwQuadrature.quadWave(t, segs, Ac, w, ph);
                    wfIndex(chan, i) = obj.AWG.AddWaveformByFunction(f, dur + minWaveformTime, 0, 1);
                end
            end
            obj.AWG.LoadWaveform();

            if obj.useSequenceWalk
                % One sequence per channel, one [zero pad, pulse] block per configuration:
                % exactly the Rabi's walk blocks for a tau shorter than its base step, so the
                % Rabi's uploadWalkSequence and performWalk run it unchanged.
                p = struct();
                p.initZeroIndex       = initZeroIndex;
                p.constStepIndex      = initZeroIndex;      % never played: no base steps
                p.waveformsIndex      = wfIndex;
                p.numZeroStepsInPoint = zeros(1, N);
                p.nWarm               = obj.walkWarmupTriggers;
                obj.walkParts = p;
                obj.uploadWalkSequence(randperm(N));
            else
                obj.indexSeq = zeros(2, N);
                for i = 1:N
                    for chan = 1:2
                        obj.indexSeq(chan, i) = obj.AWG.AddSequence( ...
                            [initZeroIndex(chan), wfIndex(chan, i)], [1 1], [1 0], 0);
                    end
                end
                obj.AWG.LoadSequence();
            end
            obj.stampLoaded(fpAtEntry);
        end
    end
end
