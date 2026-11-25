function results = awgArmLengthProbe(awg, stepCounts, skipReset)
%AWGARMLENGTHPROBE Find what sequence size the AWG will actually arm.
%
%   results = awgArmLengthProbe(awg)
%   results = awgArmLengthProbe(awg, [3 12 24 48 96 192 288])
%
% Established by experiment: a 3-step sequence arms in ~0.5 s (with a tail to
% 4.5 s), and a 288-step sequence never arms at all -- 30 s with an empty error
% queue, so the instrument accepts it and then refuses to run it. Somewhere
% between those two is a ceiling, and where it sits decides whether a CHUNKED
% sweep is viable: eight sequences of twelve frequencies instead of ninety-six of
% one would cut re-assignments 12x, but only if a ~36-step sequence still arms.
%
% Each size is tried TWICE, which is the point of the probe:
%
%   1 trigger    only the first step waits for a trigger. Tests whether raw STEP
%                COUNT is the limit.
%   N/3 triggers every third step waits, exactly like the concatenated sweep.
%                Tests whether the number of WAIT-FOR-TRIGGER steps is the limit.
%
% If "many steps, one trigger" arms but "many steps, many triggers" does not, the
% ceiling is on trigger-wait steps and chunking must be sized by trigger count.
% If both fail at the same size, it is plain step count.
%
% SAFETY: builds its sequences entirely from a ZERO-amplitude waveform, so no
% microwaves are produced at any point. It does overwrite the AWG's waveform and
% sequence lists, so whatever was loaded is gone -- run it between experiments,
% not during one. Invalidates the load cache on the way out.

    if ~exist('stepCounts', 'var') || isempty(stepCounts)
        stepCounts = [3 12 24 48 96 192 288];
    end
    if ~exist('skipReset', 'var') || isempty(skipReset)
        skipReset = false;
    end

    % Start from a known state. The first version of this probe ran straight after
    % a failed load and reported that even a 3-step sequence would not arm -- which
    % cannot be true, since that is the shape the classic sweep uses every point.
    % So the instrument was in a bad state, and every result was meaningless.
    if ~skipReset
        fprintf('Resetting the AWG first (*RST is slow, ~10-30 s)...\n');
        awg.Reset();
        fprintf('Done.\n');
    end

    % ---- CONTROL -----------------------------------------------------------
    % 3 steps, one trigger wait: exactly a classic ESR sequence. If this does not
    % arm, the probe or the instrument is broken and no other row means anything.
    fprintf('\nControl: 3 steps, 1 trigger (the shape the classic sweep uses)...\n');
    [ctrlArmed, ctrlSeconds, ctrlMsg] = tryArm(awg, 3, 4);
    if ctrlArmed
        fprintf('  control ARMED in %.1f s -- results below are meaningful.\n\n', ctrlSeconds);
    else
        fprintf('  control REFUSED after %.1f s: %s\n', ctrlSeconds, ctrlMsg);
        fprintf(['\n  STOP. A 3-step sequence is what the classic sweep arms at every point,\n' ...
                 '  so if it will not arm here the instrument is in a bad state, not at a\n' ...
                 '  size limit. Power-cycle the AWG, re-run main, and try again. Do not draw\n' ...
                 '  any conclusion about sequence length from this run.\n\n']);
        results = struct('steps',{},'nTriggers',{},'buildSeconds',{},'armSeconds',{},'armed',{});
        return
    end

    fprintf('\n=== awgArmLengthProbe ===\n');
    fprintf('%-8s %-14s %-12s %-12s %s\n', 'steps', 'trigger steps', 'build (s)', 'arm (s)', 'result');

    results = struct('steps',{},'nTriggers',{},'buildSeconds',{},'armSeconds',{},'armed',{});

    for n = stepCounts(:).'
        for variant = 1:2
            if variant == 1
                trigEvery = n + 1;          % effectively: only step 1 waits
            else
                trigEvery = 3;              % same shape as the concatenated sweep
            end

            entry = struct('steps', n, 'nTriggers', NaN, 'buildSeconds', NaN, ...
                           'armSeconds', NaN, 'armed', false);
            try
                [armed, armSeconds, msg, buildSeconds, nTrig] = tryArm(awg, n, trigEvery);
                entry.nTriggers    = nTrig;
                entry.buildSeconds = buildSeconds;
                entry.armSeconds   = armSeconds;
                entry.armed        = armed;
                if armed
                    verdict = 'ARMED';
                else
                    verdict = sprintf('REFUSED - %s', msg);
                end
                fprintf('%-8d %-14d %-12.1f %-12.1f %s\n', ...
                        n, entry.nTriggers, entry.buildSeconds, entry.armSeconds, verdict);
            catch err
                fprintf('%-8d %-14s %-12s %-12s FAILED - %s\n', ...
                        n, '?', '?', '?', err.message);
            end
            results(end+1) = entry;     %#ok<AGROW>
        end
    end

    % --- where is the ceiling? --------------------------------------------
    fprintf('\n');
    for variant = 1:2
        sel = results(variant:2:end);
        armed = [sel.armed];
        label = BooleanHelper.ifTrueElse(variant == 1, 'one trigger step', 'trigger every 3rd step');
        if all(armed)
            fprintf('  %-24s : armed at every size tried, up to %d steps\n', ...
                    label, max([sel.steps]));
        elseif ~any(armed)
            fprintf('  %-24s : never armed, even at %d steps\n', label, min([sel.steps]));
        else
            lastOk = max([sel(armed).steps]);
            firstBad = min([sel(~armed).steps]);
            fprintf('  %-24s : armed up to %d steps, refused from %d\n', label, lastOk, firstBad);
        end
    end

    okTrig = results(2:2:end);
    armedTrig = [okTrig.armed];
    if any(armedTrig)
        best = max([okTrig(armedTrig).steps]);
        fprintf('\n  Largest sweep-shaped sequence that arms: %d steps = %d frequencies.\n', ...
                best, floor(best/3));
        fprintf('  A chunked sweep would need ceil(nFreq / %d) assignments per pass\n', floor(best/3));
        fprintf('  instead of nFreq -- for 96 points that is %d instead of 96.\n', ...
                ceil(96 / max(floor(best/3), 1)));
    else
        fprintf('\n  No sweep-shaped sequence armed at any size. Chunking will not help.\n');
    end

    awg.ResetStoredData;
    awg.invalidateLoaded();
    fprintf('\n  Cleared the AWG waveform/sequence lists. Reload before your next run.\n\n');
end

function [armed, armSeconds, msg, buildSeconds, nTrig] = tryArm(awg, nSteps, trigEvery)
% Build an nSteps sequence from a single zero waveform and try to arm it.
% Used for both the control and every probe row, so they cannot differ.
    awg.ResetStoredData;

    % The real LoadAWG does BOTH of these before arming, and the AWG will not
    % enter Run with its outputs disabled -- which is what invalidated the first
    % version of this probe.
    awg.workChannels = 1;

    zeroIdx = awg.AddWaveformByFunction(@(t) 0*t, awg.minWaveformDuration, 0, 1);

    wf  = repmat(zeroIdx, 1, nSteps);
    rpt = ones(1, nSteps);
    trg = zeros(1, nSteps);
    trg(1:trigEvery:nSteps) = 1;        % trigger-wait steps; rcount 1 as required
    nTrig = sum(trg);

    tBuild = tic;
    awg.AddSequence(wf, rpt, trg, 0);
    awg.LoadWaveform();
    awg.LoadSequence();
    awg.assignSequence(1, 1);
    buildSeconds = toc(tBuild);

    msg = '';
    tArm = tic;
    try
        awg.Run;
        armSeconds = toc(tArm);
        armed = true;
    catch armErr
        armSeconds = toc(tArm);
        armed = false;
        msg = armErr.message;
        if length(msg) > 70
            msg = msg(1:70);
        end
    end
end
