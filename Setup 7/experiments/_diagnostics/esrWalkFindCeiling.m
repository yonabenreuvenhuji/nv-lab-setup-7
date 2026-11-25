function results = esrWalkFindCeiling(a, centerFreq, stepMHz, pointCounts)
%ESRWALKFINDCEILING How many frequencies fit in one walk-mode sequence.
%
%   results = esrWalkFindCeiling(a)
%   results = esrWalkFindCeiling(a, 2845, 2, [12 24 48 96])
%
% WHY THIS AND NOT ANOTHER PROBE. Synthetic probes have now shown that the AWG
% arms a 288-step, 96-trigger, two-channel sequence in under a second -- but those
% sequences referenced ONE shared zero waveform with rcount = 1. The real walk
% sequence references 96 DISTINCT 16000-sample bursts with rcount ~110 each, and
% that refuses to arm. So instead of guessing which of those differences matters,
% this exercises the actual LoadAWGConcatenated path with a shrinking sweep.
%
% Two things come out of it:
%   1. The largest frequency count that arms -- which is exactly the chunk size a
%      chunked sweep would use.
%   2. Whether walk mode's DATA is right at that size, by comparing against a
%      classic run of the same sweep. Arming is necessary but not sufficient; the
%      jump mode fiasco armed fine and produced scrambled frequencies.
%
% Stops at the first size that fails to arm, and always leaves the experiment on
% the classic path.
%
% 'a' must be a configured AWGExpESR. centerFreq should sit on your resonance so
% that even a short sweep has contrast to correlate -- from the classic runs in
% this session the dip is around 2840-2850 MHz.

    if ~exist('centerFreq', 'var') || isempty(centerFreq);  centerFreq = 2845;  end
    if ~exist('stepMHz', 'var')    || isempty(stepMHz);     stepMHz = 2;        end
    if ~exist('pointCounts', 'var')|| isempty(pointCounts); pointCounts = [12 24 48 96]; end

    if ~isa(a, 'AWGExpESR')
        error('esrWalkFindCeiling expects an AWGExpESR, got a %s', class(a));
    end

    fprintf('\n=== esrWalkFindCeiling: centre %g MHz, %g MHz steps ===\n', centerFreq, stepMHz);
    fprintf('%-8s %-8s %-12s %-14s %-12s %s\n', ...
            'points', 'steps', 'walk armed', 'correlation', 'best shift', 'verdict');

    results = struct('points',{},'steps',{},'armed',{},'correlation',{}, ...
                     'bestShift',{},'dataOk',{});

    % Steps per frequency depends on the layout: walk mode with walkCompactSteps
    % drops the trailing zero, and walkWarmupTriggers adds one block up front.
    % Arm time is linear in step count, so this is the number that matters.
    stepsPerBlock = BooleanHelper.ifTrueElse(a.walkCompactSteps, 2, 3);
    stepsFor = @(n) stepsPerBlock * (n + a.walkWarmupTriggers);

    for nPts = pointCounts(:).'
        entry = struct('points', nPts, 'steps', stepsFor(nPts), 'armed', false, ...
                       'correlation', NaN, 'bestShift', NaN, 'dataOk', false);

        a.useSequenceJump = false;
        a.useSequenceWalk = false;
        a.frequency = centerFreq + ((0:nPts-1) - (nPts-1)/2) * stepMHz;

        % ---- classic reference at this sweep length -----------------------
        ref = [];
        try
            a.run;
            ref  = a.signalParam.value(:);
            refE = a.signalParam.sterr(:);
        catch err
            fprintf('%-8d %-8d %-12s %-14s %-12s classic reference FAILED - %s\n', ...
                    nPts, stepsFor(nPts), '-', '-', '-', shortMsg(err.message));
            results(end+1) = entry;   %#ok<AGROW>
            continue
        end

        % ---- walk mode at the same sweep length --------------------------
        walk = [];
        try
            a.useSequenceWalk = true;
            a.run;
            walk  = a.signalParam.value(:);
            walkE = a.signalParam.sterr(:);
            entry.armed = true;
        catch err
            a.useSequenceWalk = false;
            msg = shortMsg(err.message);
            fprintf('%-8d %-8d %-12s %-14s %-12s %s\n', ...
                    nPts, stepsFor(nPts), 'NO', '-', '-', msg);
            results(end+1) = entry;   %#ok<AGROW>
            fprintf(['\n  Stopped: %d frequencies (%d steps) is past the ceiling.\n' ...
                     '  The largest size that armed is the chunk size to build against.\n'], ...
                    nPts, stepsFor(nPts));
            break
        end
        a.useSequenceWalk = false;

        % ---- is the data actually right? ---------------------------------
        if isempty(walk) || numel(walk) ~= numel(ref)
            fprintf('%-8d %-8d %-12s %-14s %-12s armed but produced %d of %d values\n', ...
                    nPts, stepsFor(nPts), 'yes', '-', '-', numel(walk), numel(ref));
            results(end+1) = entry;   %#ok<AGROW>
            continue
        end

        % Judge agreement against the NOISE-LIMITED ceiling, not against 1.
        % Two correct but noisy measurements cannot correlate at 0.9 when the
        % sweep's own SNR caps them lower -- that mistake failed walk mode
        % repeatedly while its graph visibly matched the reference.
        mm = walkCompareMetrics(ref, refE, walk, walkE);
        entry.correlation = mm.correlation;
        entry.bestShift   = mm.bestShift;
        entry.dataOk      = mm.pass;

        if mm.pass
            verdict = sprintf('MATCHES (r/rMax = %.2f, %.1f sigma)', mm.rRatio, mm.sigmaDiff);
        elseif mm.bestShift ~= 0
            verdict = sprintf('shifted by %+d (r there %.2f)', mm.bestShift, mm.bestShiftCorr);
        else
            verdict = sprintf('does NOT match (r/rMax = %.2f, %.1f sigma)', mm.rRatio, mm.sigmaDiff);
        end
        fprintf('%-8d %-8d %-12s %-+14.4f %-12d %s\n', ...
                nPts, stepsFor(nPts), 'yes', entry.correlation, entry.bestShift, verdict);

        results(end+1) = entry;   %#ok<AGROW>
    end

    % ---- what to do with it ----------------------------------------------
    a.useSequenceWalk = false;
    a.useSequenceJump = false;

    good = results([results.dataOk]);
    fprintf('\n');
    if isempty(good)
        armedOnly = results([results.armed]);
        if isempty(armedOnly)
            fprintf('  Walk mode did not arm at any size, not even %d frequencies.\n', min(pointCounts));
            fprintf('  Nothing about sequence length explains that -- stop here.\n');
        else
            fprintf('  Walk mode armed up to %d frequencies but the data never matched.\n', ...
                    max([armedOnly.points]));
            fprintf('  Arming is not the problem; the acquisition mapping is. Send me the\n');
            fprintf('  correlations and best shifts above.\n');
        end
    else
        best = max([good.points]);
        fprintf('  Largest sweep that arms AND gives correct data: %d frequencies (%d steps).\n', ...
                best, stepsFor(best));
        nChunks = ceil(96 / best);
        fprintf('  For a 96-point sweep that is %d chunk(s) -- %d sequence assignments per\n', ...
                nChunks, nChunks);
        fprintf('  pass instead of 96, so roughly %.0f s of arming instead of ~48 s.\n', ...
                nChunks * 0.5);
        if best >= max(pointCounts)
            fprintf('  NOTE that was the largest size tried, so the real ceiling may be higher.\n');
        end
    end
    fprintf('  Left useSequenceWalk = false.\n\n');
end

function m = shortMsg(msg)
    m = strtrim(msg);
    m = strrep(m, newline, ' ');
    if length(m) > 58
        m = [m(1:55) '...'];
    end
end

function r = corr2col(x, y)
    x = x(:) - mean(x(:));
    y = y(:) - mean(y(:));
    d = sqrt(sum(x.^2) * sum(y.^2));
    if d == 0
        r = 0;
    else
        r = sum(x .* y) / d;
    end
end
