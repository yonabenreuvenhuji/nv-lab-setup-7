function results = esrWalkCacheCheck(a)
%ESRWALKCACHECHECK Time repeated walk runs against a warm AWG.
%
%   results = esrWalkCacheCheck(a)
%
% WHAT THE FIRST VERSION OF THIS FILE GOT WRONG, twice over.
%
% 1. It called a.run twice. Experiment.run resumes from currIter, so once
%    currIter has reached averages the second call returns immediately with the
%    OLD data -- 4.7 s and a byte-identical contrast, which is not a fast run, it
%    is no run at all. The reuse measurement has to go through a.restart, which
%    resets the average counter while leaving the loaded waveforms alone.
%
% 2. Its verdict logic then read that 4.7 s no-op as "the setup cost is being
%    paid again", which is backwards. What the numbers actually said is that a
%    warm walk run takes 28.9 s.
%
% AND THAT 28.9 s IS THE RESULT THAT MATTERS. I had attributed the 161.7 s
% validation run as 73 s setup + 89 s acquisition, and concluded from that
% acquisition figure that the sweep was trigger-bound and walk mode could never
% be much faster than classic. A cache-hit walk run finishing the same 96 x 100
% sweep in 28.9 s says that was wrong: the setup in that run was ~130 s, not
% 73 s, because the waveform upload is not in either of the printed settle times.
%
% So acquisition is ~29 s on both paths, and the ~94 s difference between that
% and classic's 123.2 s was never acquisition -- it was 96 asset assignments and
% 96 arms. Walk mode removes essentially all of it. That makes it about 4x
% faster than classic on a warm AWG, on top of having no misarmed outliers.
%
% Leaves the experiment on the classic path.

    if ~isa(a, 'AWGExpESR')
        error('esrWalkCacheCheck expects an AWGExpESR, got a %s', class(a));
    end

    fprintf('\n=== esrWalkCacheCheck: %d points, repeats = %d ===\n', ...
            length(a.frequency), a.repeats);

    results = struct('seconds', [], 'contrast', [], 'ranProperly', []);

    try
        a.useSequenceJump = false;
        a.useSequenceWalk = true;

        for k = 1:2
            if k == 1
                fprintf('\n--- walk run 1 (reload if the sweep changed) ---\n');
            else
                fprintf('\n--- walk run 2, via restart (must reuse the loaded sequence) ---\n');
            end

            before = a.currIter;
            t0 = tic;               % local tic: the AWG class runs its own tic/toc,
            if k == 1               % so a bare toc here would read the wrong one
                a.run;
            else
                a.restart;          % NOT a.run -- see note 1 above
            end
            results.seconds(k) = toc(t0);

            v = a.signalParam.value(:);
            results.contrast(k) = 100*(1 - min(v));
            % A run that did not advance the average counter did not acquire
            % anything, whatever its wall time says.
            results.ranProperly(k) = (a.currIter > 0) && ~(k == 2 && a.currIter == before);
            fprintf('  run %d: %.1f s, contrast %.2f%%%s\n', k, results.seconds(k), ...
                    results.contrast(k), BooleanHelper.ifTrueElse(results.ranProperly(k), ...
                    '', '   <-- NO ACQUISITION, this number is meaningless'));
        end
    catch err
        a.useSequenceWalk = false;
        fprintf('\n  aborted: %s\n', err.message);
        rethrow(err);
    end

    a.useSequenceWalk = false;

    fprintf('\n');
    if ~all(results.ranProperly)
        fprintf(['  One of the runs did not acquire -- ignore the timings. Check a.averages\n' ...
                 '  and a.currIter; run returns immediately once currIter has reached\n' ...
                 '  averages.\n']);
    else
        warm = min(results.seconds);
        fprintf('  Warm walk run: %.1f s for %d points x %d repeats.\n', ...
                warm, length(a.frequency), a.repeats);
        fprintf(['  Compare your classic run of the same sweep. If warm walk is well under\n' ...
                 '  it, the per-point asset assignment was the whole cost and walk mode is\n' ...
                 '  both the faster and the cleaner path once the sequence is loaded.\n']);
    end
    fprintf('  Left useSequenceWalk = false.\n\n');
end
