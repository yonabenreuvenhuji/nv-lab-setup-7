function results = awgTwoChannelArmProbe(awg, stepCounts)
%AWGTWOCHANNELARMPROBE Why a long 2-channel sequence refuses to arm.
%
%   results = awgTwoChannelArmProbe(awg)
%   results = awgTwoChannelArmProbe(awg, [3 48 96 288])
%
% CONTEXT. awgArmLengthProbe showed that a 288-step sequence with 96 trigger-wait
% steps arms fine in 3.5 s -- on ONE channel. But LoadAWGConcatenated builds TWO
% SEPARATE single-track sequences and assigns one to each channel, and that
% refuses to arm after 30 s with an empty error queue. One channel works, two do
% not, so the channel arrangement is the remaining suspect.
%
% On the AWG70002A both channels share one clock and one sequencer. The intended
% way to drive them together is ONE sequence with TWO TRACKS -- which is what the
% 'tasset1' ("track asset 1") mnemonic already in WriteSequenceToAWG is for --
% rather than two independent sequences running independently.
%
% So each size is tried two ways:
%
%   A  two sequences, one track each   assigned "1"->ch1, "2"->ch2.
%      Exactly what LoadAWGConcatenated does now.
%
%   B  one sequence, two tracks        assigned ("T",1)->ch1, ("T",2)->ch2.
%      The arrangement the hardware is built around.
%
% If A fails at a length where B succeeds, the fix is to teach WriteSequenceToAWG
% about tracks, and walk mode comes back. If both fail at two channels, the
% concatenated approach cannot work on this instrument and we stop.
%
% SAFETY: every step plays a ZERO-amplitude waveform, so no microwaves are
% produced. It overwrites the waveform and sequence lists -- run between
% experiments. Resets the AWG first and clears the lists afterwards.

    if ~exist('stepCounts', 'var') || isempty(stepCounts)
        stepCounts = [3 48 96 288];
    end

    fprintf('Resetting the AWG first (*RST is slow)...\n');
    awg.Reset();
    fprintf('Done.\n');

    % ---- CONTROL: 3 steps, two channels, arrangement A ---------------------
    % This is what the classic sweep does at every point with nChannels = 2, so
    % it MUST arm. If it does not, the instrument is in a bad state and nothing
    % below means anything.
    fprintf('\nControl: 3 steps, 2 channels, two single-track sequences...\n');
    [ok, secs, msg] = tryTwoSeq(awg, 3, 3);
    if ok
        fprintf('  control ARMED in %.1f s -- results below are meaningful.\n', secs);
    else
        fprintf('  control REFUSED after %.1f s: %s\n', secs, msg);
        fprintf(['\n  STOP. This is the classic sweep''s own configuration, so a refusal here\n' ...
                 '  means the instrument is in a bad state. Power-cycle it, re-run main, and\n' ...
                 '  try again. Draw no conclusions from this run.\n\n']);
        results = struct('steps',{},'arrangement',{},'armSeconds',{},'armed',{});
        return
    end

    fprintf('\n=== awgTwoChannelArmProbe ===\n');
    fprintf('%-8s %-34s %-12s %s\n', 'steps', 'arrangement', 'arm (s)', 'result');

    results = struct('steps',{},'arrangement',{},'armSeconds',{},'armed',{});

    for n = stepCounts(:).'
        for variant = 1:2
            if variant == 1
                label = 'A: 2 sequences, 1 track each';
                [ok, secs, msg] = tryTwoSeq(awg, n, 3);
            else
                label = 'B: 1 sequence, 2 tracks';
                [ok, secs, msg] = tryTwoTrack(awg, n, 3);
            end
            verdict = BooleanHelper.ifTrueElse(ok, 'ARMED', sprintf('REFUSED - %s', msg));
            fprintf('%-8d %-34s %-12.1f %s\n', n, label, secs, verdict);
            results(end+1) = struct('steps', n, 'arrangement', label, ...
                                    'armSeconds', secs, 'armed', ok);   %#ok<AGROW>
        end
    end

    % ---- verdict -----------------------------------------------------------
    A = results(1:2:end); B = results(2:2:end);
    fprintf('\n');
    fprintf('  A (what the code does now) : %s\n', ceilingText(A));
    fprintf('  B (one sequence, 2 tracks) : %s\n', ceilingText(B));

    aMax = maxArmed(A); bMax = maxArmed(B);
    fprintf('\n');
    if bMax > aMax
        fprintf('  => Tracks are the answer. B reaches %d steps where A stops at %d.\n', bMax, aMax);
        fprintf('     WriteSequenceToAWG needs track support and walk mode is viable again.\n');
    elseif bMax == aMax && bMax >= max(stepCounts)
        fprintf('  => Both arrangements arm at every size tried, so the arm failure in\n');
        fprintf('     LoadAWGConcatenated is NOT about channels or length. The remaining\n');
        fprintf('     difference is the waveforms: 195 real bursts of 16000 samples each,\n');
        fprintf('     versus one shared zero waveform here.\n');
    else
        fprintf('  => Neither arrangement scales on two channels (A to %d, B to %d steps).\n', aMax, bMax);
        fprintf('     The concatenated approach cannot work on this instrument.\n');
    end

    awg.ResetStoredData;
    awg.invalidateLoaded();
    fprintf('\n  Cleared the AWG lists. Reload before your next run.\n\n');
end


function [armed, armSeconds, msg] = tryTwoSeq(awg, nSteps, trigEvery)
% Arrangement A: two separate single-track sequences, one per channel.
% Identical to what LoadAWGConcatenated does today.
    zeroIdx = prepareZeroWaveform(awg, 2);
    wf  = repmat(zeroIdx, 1, nSteps);
    rpt = ones(1, nSteps);
    trg = zeros(1, nSteps); trg(1:trigEvery:nSteps) = 1;

    awg.AddSequence(wf, rpt, trg, 0);       % sequence "1"
    awg.AddSequence(wf, rpt, trg, 0);       % sequence "2"
    awg.LoadSequence();
    awg.assignSequence(1, 1);
    awg.assignSequence(2, 2);
    [armed, armSeconds, msg] = timedArm(awg);
end

function [armed, armSeconds, msg] = tryTwoTrack(awg, nSteps, trigEvery)
% Arrangement B: ONE sequence with two tracks, one track per channel. Built with
% raw SCPI because WriteSequenceToAWG only knows about single-track sequences.
    prepareZeroWaveform(awg, 2);
    name = 'T';

    awg.sendCommand('slist:sequence:delete ALL');
    awg.sendCommand(sprintf('slist:seq:new "%s", %d, 2', name, nSteps));   % 2 TRACKS

    cmds = cell(1, 0);
    for r = 1:nSteps
        cmds{end+1} = sprintf('slist:seq:step%d:tasset1:wav "%s", "1"', r, name);  %#ok<AGROW>
        cmds{end+1} = sprintf('slist:seq:step%d:tasset2:wav "%s", "1"', r, name);  %#ok<AGROW>
        cmds{end+1} = sprintf('slist:seq:step%d:rcount "%s", 1', r, name);         %#ok<AGROW>
        if mod(r - 1, trigEvery) == 0
            cmds{end+1} = sprintf('slist:seq:step%d:winput "%s", atrigger', r, name); %#ok<AGROW>
        else
            cmds{end+1} = sprintf('slist:seq:step%d:winput "%s", off', r, name);      %#ok<AGROW>
        end
    end
    cmds{end+1} = sprintf('slist:seq:step%d:goto "%s", first', nSteps, name);
    awg.sendCommandBatch(cmds);
    awg.sendCommand('*WAI');

    awg.sendCommand(sprintf('source1:casset:sequence "%s", 1', name));
    awg.sendCommand(sprintf('source2:casset:sequence "%s", 2', name));
    [armed, armSeconds, msg] = timedArm(awg);
end

function zeroIdx = prepareZeroWaveform(awg, nChannels)
% One shared zero-amplitude waveform, outputs enabled. The AWG will not enter Run
% with its outputs off, which invalidated an earlier probe entirely.
    awg.ResetStoredData;
    awg.workChannels = 1:nChannels;
    zeroIdx = awg.AddWaveformByFunction(@(t) 0*t, awg.minWaveformDuration, 0, 1);
    awg.LoadWaveform();
end

function [armed, armSeconds, msg] = timedArm(awg)
    msg = '';
    t = tic;
    try
        awg.Run;
        armSeconds = toc(t);
        armed = true;
    catch err
        armSeconds = toc(t);
        armed = false;
        msg = err.message;
        if length(msg) > 60
            msg = msg(1:60);
        end
    end
end

function txt = ceilingText(sel)
    armed = [sel.armed];
    if all(armed)
        txt = sprintf('armed at every size tried, up to %d steps', max([sel.steps]));
    elseif ~any(armed)
        txt = sprintf('never armed, even at %d steps', min([sel.steps]));
    else
        txt = sprintf('armed up to %d steps, refused from %d', ...
                      max([sel(armed).steps]), min([sel(~armed).steps]));
    end
end

function m = maxArmed(sel)
    armed = [sel.armed];
    if any(armed)
        m = max([sel(armed).steps]);
    else
        m = 0;
    end
end
