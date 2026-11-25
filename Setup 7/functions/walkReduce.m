function [signal, sterrOut] = walkReduce(raw, nDet, nFreq, opts)
%WALKREDUCE Reduce interleaved sweep-walk data to per-frequency signal and error.
%
%   [signal, sterrOut] = walkReduce(raw, nDet, nFreq, opts)
%
% For the useSequenceWalk acquisition mode, where the AWG sequencer advances one
% frequency per trigger and a single pg.run therefore returns a whole pass over
% the sweep. Each frequency's "repeats" are the PASSES, which arrive spread out
% in time rather than back to back.
%
% This is a faithful transcription of Experiment.processData applied per
% frequency, so the normalisation, the discard-the-first-repeat rule and the
% photodiode branch all behave identically to the classic per-point path.
% testWalkReduce checks that equivalence against a replica of processData.
%
% INPUTS
%   raw   [nDet x nFreq x nPasses]. raw(1,k,p) is the SIGNAL window (MW on) for
%         frequency k on pass p; raw(2,k,p) is the REFERENCE window (MW off).
%         Comes from the digitizer stream, which arrives as
%         [s1 r1 s2 r2 ... sN rN] per pass.
%   nDet  detection periods per repeat. Must be EVEN: the windows are consumed as
%         consecutive (signal, reference) PAIRS, so
%             nDet = 2   [signal reference]                    ESR, Rabi
%             nDet = 4   [sig ref sig ref]                      Ramsey, Echo, T1 with
%                                                               doubleMeasurement on
%         Each pair is reduced exactly as a lone pair would be -- including the
%         photodiode normalisation, which pairs a signal window with ITS OWN
%         reference rather than with a global one. With nDet = 2 this file behaves
%         identically to the version that only accepted 2.
%   nFreq number of frequency points.
%
%         nDet and nFreq are MANDATORY rather than inferred from size(raw),
%         because MATLAB drops trailing singleton dimensions: raw(:,:,1:1) comes
%         through as [2 x nFreq], not [2 x nFreq x 1], so a shape check on the
%         incoming array fails on the very first pass when passesPerRun is 1.
%         Reshaping against explicit dimensions is immune to that.
%   opts  struct with fields:
%           hasPhotodiode              logical, spcm.hasPhotodiode
%           detectionWithGI            logical, spcm.detectionWithGI
%           GIgain                     double,  spcm.GIgain (only used with GI)
%           detectionDuration          us
%           referenceDetectionDuration us
%
% OUTPUTS
%   signal, sterrOut   [nDet x nFreq]. Rows alternate signal, reference, signal,
%                      reference..., matching the layout of obj.signal(:, :, currIter)
%                      that processData produces.

    if ~ValidationHelper.isValuePositiveInteger(nDet) || ~ValidationHelper.isValuePositiveInteger(nFreq)
        error('walkReduce needs positive integer nDet and nFreq');
    end
    nPasses = numel(raw) / (nDet * nFreq);
    if mod(nPasses, 1) ~= 0 || nPasses < 1
        error(['walkReduce got %d elements, which is not a whole number of passes ' ...
               'for nDet = %d and nFreq = %d'], numel(raw), nDet, nFreq);
    end
    raw = reshape(raw, nDet, nFreq, nPasses);
    if mod(nDet, 2) ~= 0
        error(['walkReduce needs an EVEN number of detection periods per repeat -- they ' ...
               'are consumed as (signal, reference) pairs -- got %d'], nDet);
    end
    nPair = nDet / 2;

    signal   = zeros(nDet, nFreq);
    sterrOut = zeros(nDet, nFreq);

    % processData drops the first repeat when there are more than 10 of them,
    % because the first detection is not preceded by a full init period. The
    % passes play the role of repeats here, so the same rule applies to them.
    start = 1 + double(nPasses > 10);

    musec = 1e-6;
    kc    = 1e3;
    timeNormalization = [opts.detectionDuration, opts.referenceDetectionDuration] * musec;

    for k = 1:nFreq
      for q = 1:nPair
        % Two columns, one row per pass -- the same shape processData builds with
        % s = reshape(rawData, m, n)' for m = 2 reads per repeat. For nDet > 2 this
        % runs once per PAIR, which is what makes the multi-window experiments work:
        % Ramsey's rows 1-2 are the |1> projection and rows 3-4 the |0> one, and they
        % are independent measurements that must not be mixed.
        dS = 2*q - 1; dR = 2*q;
        s = [squeeze(raw(dS, k, :)), squeeze(raw(dR, k, :))];

        if opts.hasPhotodiode
            % Normalise per pass before averaging, to suppress low-frequency
            % noise, then park ALL the noise in column 1 and make column 2 a
            % constant. That is deliberate upstream: it makes the ratio's error
            % exact and stops the reference contributing twice.
            refSignal = s(:,1) ./ s(:,2);
            ref = nanMeanCol(s(start:end, 2));
            s(:,1) = ref .* refSignal;
            s(:,2) = ref + 0 * refSignal;
        end

        sig = nanMeanCol(s(start:end, :));
        err = ste(s(start:end, :));

        if ~opts.hasPhotodiode
            sig = sig ./ timeNormalization / kc;    % kcounts per second
            err = err ./ timeNormalization / kc;
        elseif opts.detectionWithGI
            tn  = timeNormalization / musec;
            sig = sig ./ tn ./ opts.GIgain;
            err = err ./ tn ./ opts.GIgain;
        end

        signal([dS dR], k)   = sig(:);
        sterrOut([dS dR], k) = err(:);
      end
    end
end

function m = nanMeanCol(x)
    % Column-wise mean ignoring NaN. Spelled out rather than using
    % mean(x,1,'omitnan') so the behaviour is identical under Octave, where the
    % offline test runs.
    m = zeros(1, size(x, 2));
    for c = 1:size(x, 2)
        col = x(:, c);
        col = col(~isnan(col));
        if isempty(col)
            m(c) = NaN;
        else
            m(c) = sum(col) / numel(col);
        end
    end
end
