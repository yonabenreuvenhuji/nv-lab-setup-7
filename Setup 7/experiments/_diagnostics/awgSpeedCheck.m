function awgSpeedCheck(expObj)
%AWGSPEEDCHECK Report where an AWG experiment spends its time.
%
%   awgSpeedCheck(expObj)
%
% Run this AFTER the experiment has been run (or at least prepared) once, so
% that the waveforms/sequences are loaded and expObj.indexSeq is populated:
%
%   a = AWGExpESR; a.loadJsonParams('all', jsonLocation);
%   a.frequency = 2870 + (-30:0.5:30); a.repeats = 150; a.averages = 1;
%   a.run;
%   awgSpeedCheck(a)
%
% It reports the SCPI command budget, the measured per-frequency-point
% overhead, and what fraction of the wall clock is actually collecting
% photons. Use it before and after changing awg.commDelayTime /
% awg.skipRedundantDeletes so you can see whether a change did anything
% instead of guessing.
%
% Rough expectation for a 121-point ESR sweep:
%   legacy (commDelayTime = 0.005, skipRedundantDeletes = false)
%       ~2190 commands, ~11 s of busy-wait alone, ~12 ms/point overhead
%   patched (commDelayTime = 0, skipRedundantDeletes = true)
%       ~1950 commands, 0 s of busy-wait,        ~2-4 ms/point overhead
%
% NOTE this briefly reassigns sequences on the AWG to time the assignment,
% so do not call it while an experiment is running.

    if ~isa(expObj, 'ExperimentAWG')
        error('awgSpeedCheck expects an ExperimentAWG subclass, got a %s', class(expObj));
    end
    if isempty(expObj.indexSeq) || all(expObj.indexSeq(:) == 0)
        error(['expObj.indexSeq is empty -- the AWG has not been loaded yet. ' ...
               'Run the experiment once (averages = 1 is enough) and call this again.']);
    end

    awg = expObj.AWG;

    fprintf('\n=== awgSpeedCheck: %s ===\n', expObj.NAME);
    fprintf('  commDelayTime        = %g s%s\n', awg.commDelayTime, ...
            repmat(' <-- legacy value, this is the 5 ms/command busy-wait', ...
                   1, awg.commDelayTime >= 0.005));
    fprintf('  skipRedundantDeletes = %d\n',   awg.skipRedundantDeletes);
    fprintf('  SampleRate           = %g MHz\n', awg.SampleRate);
    fprintf('  frequency points     = %d\n',  length(expObj.frequency));
    fprintf('  nChannels            = %d\n',  expObj.nChannels);

    %% ---- what is loaded, and what it cost to get there ----
    nWf  = length(awg.waveforms);
    nSeq = length(awg.sequences.waveforms);
    if nWf > 0
        wfLen = length(awg.waveforms{1});
        fprintf('\n  waveforms loaded      : %d  (%d samples each, %.1f MB transferred)\n', ...
                nWf, wfLen, nWf*wfLen*4/1e6);
    else
        fprintf('\n  waveforms loaded      : 0\n');
    end
    fprintf('  sequences loaded      : %d\n', nSeq);

    % Command budget per item, given the current knob settings.
    %   waveform: [delete] + new + data-header + data + *WAI
    %   sequence: [delete] + new + 3 per step + goto + *WAI
    stepsPerSeq = 0;
    if nSeq > 0
        stepsPerSeq = length(awg.sequences.waveforms{1});
    end
    cmdPerWf  = 4 + double(~awg.skipRedundantDeletes);
    cmdPerSeq = 3 + 3*stepsPerSeq + double(~awg.skipRedundantDeletes);
    nExpected = nWf*cmdPerWf + nSeq*cmdPerSeq;

    fprintf('  load command budget   : %d waveform + %d sequence = ~%d commands\n', ...
            nWf*cmdPerWf, nSeq*cmdPerSeq, nExpected);
    fprintf('    of which busy-wait  : %.1f s at commDelayTime = %g\n', ...
            nExpected*awg.commDelayTime, awg.commDelayTime);
    if ~awg.skipRedundantDeletes
        fprintf('    (%d of those are redundant deletes -- set skipRedundantDeletes = true)\n', ...
                nWf + nSeq);
    end
    fprintf('  commandCount now      : %d  (since the last resetCommandCount)\n', ...
            awg.commandCount);
    if awg.commandCount > 1.15*nExpected
        fprintf('    NOTE: that is %d more than the load budget above. Some of the excess is\n', ...
                awg.commandCount - nExpected);
        fprintf('    the per-point assignSequence during the sweep (%d points x %d chan x %d avg\n', ...
                length(expObj.frequency), expObj.nChannels, expObj.averages);
        fprintf('    = %d), but if a large remainder is left over, check whether LoadAWG ran\n', ...
                length(expObj.frequency)*expObj.nChannels*expObj.averages);
        fprintf('    twice -- prepareInternal calls obj.LoadAWG() again when an SGT is\n');
        fprintf('    registered and obj.IQ.useIQ is true.\n');
    end

    %% ---- measured per-point sweep overhead ----
    %
    % CAVEAT, read this before trusting the number below. This probe runs with
    % the pulse generator IDLE -- no triggers are arriving. An AWG armed in
    % triggered mode normally sits in "Waiting for Trigger" (RSTate = 1) during
    % a real sweep, which EnsureRunning accepts without re-arming. With no
    % triggers it can instead sit at Stopped, so this loop may pay a full
    % re-arm on every point that the real sweep does not.
    %
    % Watch the "found STOPPED" fraction printed below. If it is high here,
    % confirm against a REAL sweep before acting on it:
    %
    %     awg.resetCommandCount;  a.run;  awg.printInstrumentation
    %
    % which accumulates the same timings while the PG is actually firing.
    nProbe = min(20, size(expObj.indexSeq, 2));
    awg.resetCommandCount;
    tSweep = tic;
    for i = 1:nProbe
        awg.assignSequence(expObj.indexSeq(1,i), 1);
        if expObj.nChannels == 2
            awg.assignSequence(expObj.indexSeq(2,i), 2);
        end
        awg.EnsureRunning;
    end
    tPerPoint = toc(tSweep)/nProbe;
    nCmdSweep = awg.commandCount;
    tProbe    = awg.timing;

    fprintf('\n  per-point overhead PROBE (assignSequence + EnsureRunning, PG idle):\n');
    fprintf('    %.2f ms per frequency point  (%d commands over %d points)\n', ...
            1e3*tPerPoint, nCmdSweep, nProbe);
    fprintf('    assignSequence %.2f ms/call, re-arm %.0f ms/call, found STOPPED %d/%d points\n', ...
            1e3*tProbe.assign/max(tProbe.assignCalls,1), ...
            1e3*tProbe.run/max(tProbe.runCalls,1), ...
            tProbe.ensureStopped, tProbe.ensureCalls);
    if tProbe.ensureCalls > 0 && tProbe.ensureStopped/tProbe.ensureCalls > 0.5
        fprintf('    ^ mostly re-arm cost. This may be a probe artefact (PG idle) -- verify\n');
        fprintf('      with printInstrumentation on a real run before changing anything.\n');
    end

    %% ---- compare against the actual acquisition ----
    seqPerRepeat = 2*expObj.detectionDuration;      % signal window + reference window
    if isprop(expObj, 'delayBetweenDetections') && ~isempty(expObj.delayBetweenDetections)
        seqPerRepeat = seqPerRepeat + 2*expObj.delayBetweenDetections;
    end
    tAcq = expObj.repeats * seqPerRepeat * 1e-6;    % us -> s

    fprintf('\n  acquisition: %.2f ms per point (repeats = %d, %.0f us per repeat)\n', ...
            1e3*tAcq, expObj.repeats, seqPerRepeat);
    fprintf('  efficiency : %.0f%% of the sweep wall clock is photons\n', ...
            100*tAcq/(tAcq + tPerPoint));

    if tPerPoint > tAcq
        repeatsFor90 = ceil(9*tPerPoint/(seqPerRepeat*1e-6));
        fprintf('  *** overhead exceeds acquisition ***\n');
        reArmDominated = tProbe.ensureCalls > 0 && ...
                         tProbe.run > 0.5*(tProbe.assign + tProbe.run);
        if reArmDominated
            % Do NOT tell the user to raise repeats here. When the cost is the
            % AWG re-arming itself, "repeats" would have to reach absurd values
            % (tens of thousands) to amortise it, which exceeds the digitizer
            % buffer and would leave the tracker blind for minutes at a time.
            fprintf('      The overhead is the AWG RE-ARMING (%.0f ms per point), not SCPI\n', ...
                    1e3*tProbe.run/max(tProbe.runCalls,1));
            fprintf('      chatter. Raising repeats cannot sensibly fix that -- it would take\n');
            fprintf('      repeats = %d, far past the digitizer buffer and long enough to\n', repeatsFor90);
            fprintf('      leave the tracker blind. The per-point sequence re-assign has to go\n');
            fprintf('      instead: one sequence holding every frequency, advanced by the\n');
            fprintf('      hardware on successive triggers, so the AWG is armed once per sweep.\n');
        else
            fprintf('      Overhead is per-command, so raising repeats amortises it.\n');
            fprintf('      repeats = %d would put you at ~90%% efficiency', repeatsFor90);
            if repeatsFor90 > 5000
                fprintf(' -- but that is a lot;\n      check the digitizer buffer and your tracking interval first');
            end
            fprintf('.\n');
        end
    end

    tSweepTotal = length(expObj.frequency) * (tAcq + tPerPoint) * expObj.averages;
    fprintf('\n  projected sweep time  : %.1f s  (%d points x %d averages)\n', ...
            tSweepTotal, length(expObj.frequency), expObj.averages);
    fprintf('\n');
end
