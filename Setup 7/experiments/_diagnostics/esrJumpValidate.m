function results = esrJumpValidate(a, settleTimes)
%ESRJUMPVALIDATE Check that jump mode reproduces the classic spectrum.
%
%   results = esrJumpValidate(a)
%   results = esrJumpValidate(a, [0 0.002 0.005 0.010 0.020 0.050])
%
% Runs the CLASSIC path once as a reference, then the JUMP path at each
% candidate AWG.jumpSettleTime, and reports whether the jump spectrum actually
% matches. Leaves the experiment on the classic path when it finishes.
%
% 'a' must be a configured AWGExpESR -- frequency, repeats, amplitude and
% detection durations already set. Use enough points that the resonance is
% resolved; a 2-point check is not sensitive enough to catch misassignment
% (learned the hard way).
%
% THE METRICS, and why these three:
%
%   correlation  Pearson r between the jump and classic spectra. This is the one
%                that matters. A frequency misassignment leaves the value
%                DISTRIBUTION intact -- same contrast, same median -- while
%                destroying the frequency mapping, so contrast alone cannot
%                detect it. Measured r was +0.04 with no settle time.
%
%   scatterRatio point-to-point scatter divided by what the error bars predict.
%                A real ESR line is smooth, so this should be ~1.4. It was 5.5
%                when every point was reporting some other point's frequency.
%
%   sigmaDiff    median |classic - jump| in units of the combined error bar.
%                ~1 means the two paths agree within statistics.
%
% Verdict is PASS when correlation > 0.9, scatterRatio < 2.5 and sigmaDiff < 3.
%
% Example:
%   a = AWGExpESR; a.loadJsonParams('all', jsonLocation);
%   a.frequency = 2770:2:2960; a.repeats = 100; a.averages = 1;
%   a.amplitude = -20; a.isTracking = 0;
%   a.detectionDuration = 100; a.referenceDetectionDuration = a.detectionDuration;
%   esrJumpValidate(a)

    if ~exist('settleTimes', 'var') || isempty(settleTimes)
        settleTimes = [0 0.002 0.005 0.010 0.020 0.050];
    end
    if ~isa(a, 'AWGExpESR')
        error('esrJumpValidate expects an AWGExpESR, got a %s', class(a));
    end
    if length(a.frequency) < 20
        error(['Use at least 20 frequency points -- a short sweep cannot distinguish ' ...
               'a scrambled spectrum from a correct one.']);
    end

    awg = a.AWG;
    settleOriginal = awg.jumpSettleTime;

    % Save and clear the OTHER sequencer flag. This file predates walk mode and set
    % useSequenceJump without touching useSequenceWalk, so arriving here with walk
    % enabled -- which is now the normal state, and what a preceding diagnostic leaves
    % behind -- aborts on LoadAWG's mutual-exclusion check before a single point is
    % taken. Exactly the same omission esrFgenValidate had.
    savedWalk = a.useSequenceWalk;
    savedJump = a.useSequenceJump;
    savedFgen = a.useFgen;

    % Also worth stating plainly: armCheckEvery must NOT be 1 during this sweep. The
    % per-point AWGControl:RSTate? query is itself a command barrier -- it forces the
    % instrument to finish its queued jumps before replying -- so with it enabled every
    % settle time passes, including zero, and the sweep measures nothing. Auto ([])
    % gives 25 for the jump path, which leaves the settle time as the only synchroniser.
    % Plain fprintf, NOT a.sendWarning: sendWarning and sendError are protected methods
    % of EventSender, callable only as obj.method from inside the class. From a helper
    % script they raise "Unrecognized method, property, or field", which is a confusing
    % way to fail and has nothing to do with the thing being warned about.
    if isequal(a.armCheckEvery, 1)
        fprintf(['\nWARNING armCheckEvery = 1 during a settle-time sweep. The per-point\n' ...
                 '        AWGControl:RSTate? query is a command barrier in its own right, so\n' ...
                 '        every settle time would pass for the wrong reason, including zero.\n' ...
                 '        Setting armCheckEvery = [] for this sweep.\n']);
    end
    savedCheck = a.armCheckEvery;
    a.armCheckEvery = [];

    fprintf('\n=== esrJumpValidate: %d points, repeats = %d, averages = %d ===\n', ...
            length(a.frequency), a.repeats, a.averages);
    fprintf('jumpSettleTime is applied ONCE PER CHANNEL, so a %d-channel run pays it %d times per point.\n', ...
            a.nChannels, a.nChannels);

    try

    % ---- reference ------------------------------------------------------
    fprintf('\nreference run (classic path)...\n');
    a.useFgen         = false;
    a.useSequenceJump = false;
    a.useSequenceWalk = false;          % classic: neither sequencer path
    a.run;
    ref  = a.signalParam.value(:);
    refE = a.signalParam.sterr(:);
    refScatter = median(abs(diff(ref))) / max(1.4*median(refE), eps);
    fprintf('  contrast %.2f%%, dip %.1f sigma, scatter ratio %.2f\n', ...
            100*(1-min(ref)), (1-min(ref))/median(refE), refScatter);
    if refScatter > 2.5
        fprintf('  WARNING the CLASSIC reference is itself noisy. Any comparison below\n');
        fprintf('          is against a shaky baseline -- fix that first.\n');
    end

    % ---- candidates -----------------------------------------------------
    results = struct('settleTime',{},'correlation',{},'scatterRatio',{}, ...
                     'sigmaDiff',{},'contrast',{},'pass',{});

    fprintf('\n%-12s %-12s %-13s %-11s %-10s %s\n', ...
            'settle (ms)', 'correlation', 'scatterRatio', 'sigmaDiff', 'contrast', 'verdict');
    for st = settleTimes(:).'
        awg.jumpSettleTime = st;
        a.useSequenceWalk  = false;     % must be cleared BEFORE the jump flag is set
        a.useSequenceJump  = true;
        a.run;
        v  = a.signalParam.value(:);
        e  = a.signalParam.sterr(:);

        r       = corr2col(v, ref);
        scatter = median(abs(diff(v))) / max(1.4*median(e), eps);
        sigmaD  = median(abs(v - ref) ./ max(sqrt(e.^2 + refE.^2), eps));
        passed  = (r > 0.9) && (scatter < 2.5) && (sigmaD < 3);

        fprintf('%-12.1f %-+12.4f %-13.2f %-11.2f %-10.2f %s\n', ...
                1e3*st, r, scatter, sigmaD, 100*(1-min(v)), ...
                BooleanHelper.ifTrueElse(passed, 'PASS', 'FAIL'));

        results(end+1) = struct('settleTime', st, 'correlation', r, ...
                                'scatterRatio', scatter, 'sigmaDiff', sigmaD, ...
                                'contrast', 100*(1-min(v)), 'pass', passed); %#ok<AGROW>
    end

    % ---- pick and report -------------------------------------------------
    ok = find([results.pass], 1);       % smallest passing settle time
    fprintf('\n');
    if isempty(ok)
        awg.jumpSettleTime = settleOriginal;
        fprintf('No settle time passed. Jump mode is NOT safe to use -- left\n');
        fprintf('useSequenceJump = false and jumpSettleTime = %g.\n', settleOriginal);
        fprintf('If correlation stayed near zero at every settle time, the problem is\n');
        fprintf('not timing and the step arithmetic needs looking at instead.\n');
    else
        best = results(ok).settleTime;
        awg.jumpSettleTime = best;
        fprintf('Smallest settle time that reproduces the classic spectrum: %.0f ms\n', 1e3*best);
        fprintf('Set awg.jumpSettleTime = %g in your speed-knobs cell.\n', best);
        fprintf('Sweep cost of that delay: %.2f s per average over %d points.\n', ...
                best*length(a.frequency)*a.nChannels, length(a.frequency));
        if ok == length(results)
            fprintf('NOTE that was the largest value tried -- extend settleTimes upward\n');
            fprintf('     to check whether the true minimum is even higher.\n');
        end
    end

    catch err
        % Restore everything on the failure path too, so an abort does not leave the
        % setup on a half-configured sequencer path for the next experiment to inherit.
        a.useSequenceJump = savedJump;
        a.useSequenceWalk = savedWalk;
        a.useFgen         = savedFgen;
        a.armCheckEvery   = savedCheck;
        awg.jumpSettleTime = settleOriginal;
        rethrow(err);
    end

    a.useSequenceJump = false;          % never leave the setup on the unproven path
    a.useSequenceWalk = savedWalk;      % back to whatever path was in use before
    a.useFgen         = savedFgen;
    a.armCheckEvery   = savedCheck;
    fprintf('Left useSequenceJump = false, useSequenceWalk = %d, armCheckEvery = %s.\n\n', ...
            a.useSequenceWalk, mat2str(savedCheck));
end

function r = corr2col(x, y)
    % Pearson correlation without needing the Statistics Toolbox
    x = x(:) - mean(x(:));
    y = y(:) - mean(y(:));
    d = sqrt(sum(x.^2) * sum(y.^2));
    if d == 0
        r = 0;
    else
        r = sum(x .* y) / d;
    end
end
