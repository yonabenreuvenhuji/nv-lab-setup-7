function results = awgFgenBringUp(awg, freqMHz, dBm)
%AWGFGENBRINGUP Find out why the FGEN output never goes live.
%
%   results = awgFgenBringUp(awg)
%   results = awgFgenBringUp(awg, 2500, -20)
%
% WHERE THIS SITS. Everything else about FGEN mode has now been eliminated by
% measurement rather than argument:
%
%   frequency    sine accepted to at least 6 GHz            (awgFgenProbe)
%   amplitude    0.063246 Vpp in BOTH modes for -20 dBm     (awgFgenLevelProbe)
%   gating       driveMWSwitches = false gives a FLAT spectrum, so the switches
%                are normally closed, pass only when driven, and were gating
%                correctly all along
%   chain        AWG -> switch -> amplifier, so the amplifier sees the same
%                pulsed drive in either mode
%
% So an FGEN ESR had the right frequency, the right power and a working gate, and
% still measured nothing. The one reading that never added up is AWGCONTROL:RSTATE?
% returning 0 after AWGCONTROL:RUN:IMMEDIATE. I dismissed that once as "the run
% state means nothing in FGEN mode", which was a guess, and then spent four rounds
% on theories -- leakage, a missing RUN, a missing switch, channel phase -- each
% of which was ruled out by something I could have checked first.
%
% This checks it properly: every command's error queue, the output and run state
% after each step, and the orderings I did not try.
%
% SAFETY. This does enable the AWG output, which the earlier probes did not. That
% is safe here precisely because of the switch finding: with no experiment running
% the pulse generator is not driving channels 5 or 6, so the switches are closed
% and nothing reaches the sample. The test frequency also defaults to 2500 MHz,
% far off resonance, so even a leak could not drive the NVs. Restores the starting
% mode, turns the outputs back off, and invalidates the load cache.

    if ~exist('freqMHz', 'var') || isempty(freqMHz); freqMHz = 2500; end
    if ~exist('dBm', 'var')     || isempty(dBm);     dBm = -20;      end

    fprintf('\n=== awgFgenBringUp (%g MHz, %g dBm) ===\n', freqMHz, dBm);
    fprintf('Outputs WILL be enabled. Safe: the MW switches are undriven and closed.\n');

    results = struct();
    awg.sendCommand('*CLS');
    startMode = strtrim(awg.queryRaw('INST:MODE?'));
    startChannels = awg.workChannels;
    results.startMode = startMode;

    try
        awg.workChannels = [];
        step(awg, 'awgcontrol:stop:immediate');
        drain(awg, 20);

        % ---- variant A: what SetupFgen does now -----------------------------
        fprintf('\n--- A: outputs on, then run (what SetupFgen does) ---\n');
        results.A = tryVariant(awg, freqMHz, dBm, {'OUTPUT1:STATE ON', 'AWGCONTROL:RUN:IMMEDIATE'});

        % ---- variant B: run first, then outputs -----------------------------
        % Some instruments will not enter run with the outputs already switching.
        fprintf('\n--- B: run first, then outputs on ---\n');
        results.B = tryVariant(awg, freqMHz, dBm, {'AWGCONTROL:RUN:IMMEDIATE', 'OUTPUT1:STATE ON'});

        % ---- variant C: the class's own compound form -----------------------
        % setValue('enableOutput','1') sends exactly this. It is the only FGEN
        % start sequence in this codebase that anyone may once have used.
        fprintf('\n--- C: the compound form from setValue(''enableOutput'',''1'') ---\n');
        results.C = tryVariant(awg, freqMHz, dBm, {'OUTP1 ON;:AWGCONTROL:RUN:IMMEDIATE'});

        % ---- variant D: RUN without :IMMEDIATE ------------------------------
        fprintf('\n--- D: AWGCONTROL:RUN, no :IMMEDIATE ---\n');
        results.D = tryVariant(awg, freqMHz, dBm, {'OUTPUT1:STATE ON', 'AWGCONTROL:RUN'});

        % ---- what else might be blocking it ---------------------------------
        % If none of the above runs, the instrument may be refusing for a reason
        % it will name. These are the usual suspects on this family.
        fprintf('\n--- extra state ---\n');
        for q = {'AWGCONTROL:RSTATE?', 'OUTPUT1:STATE?', 'OUTP1?', ...
                 'INST:MODE?', 'FGEN:CHAN1:TYPE?', 'FGEN:CHAN1:FREQ?', ...
                 'FGEN:CHAN1:AMPL:POW?', 'CLOCK:SRATE?', 'CLOCK:SOURCE?', ...
                 '*ESR?', 'SYSTEM:ERROR:COUNT?'}
            v = '(unsupported)';
            try
                v = strtrim(awg.queryRaw(q{1}, 5));
            catch
            end
            fprintf('  %-26s %s\n', q{1}, v);
        end
        drain(awg, 20);

    catch err
        restore(awg, startMode, startChannels);
        fprintf('\n  bring-up aborted: %s\n', err.message);
        rethrow(err);
    end

    restore(awg, startMode, startChannels);

    % ---- verdict ---------------------------------------------------------
    fprintf('\n--- verdict ---\n');
    names = {'A', 'B', 'C', 'D'};
    ranAny = false;
    for k = 1:4
        v = results.(names{k});
        if v.running
            fprintf('  %s DID reach run state %s.\n', names{k}, v.finalRunState);
            ranAny = true;
        end
    end
    if ranAny
        fprintf(['  Use whichever variant ran, in SetupFgen, and re-run\n' ...
                 '  esrFgenValidate. Everything else about FGEN already checks out.\n']);
    else
        fprintf(['  NOTHING put the instrument into run state in FGEN mode. If the error\n' ...
                 '  queue above is empty too, then this firmware does not use the run state\n' ...
                 '  for FGEN and the output should already have been live -- in which case\n' ...
                 '  the microwave is leaving the AWG but not arriving, and that is a bench\n' ...
                 '  question: put a power meter or scope on the AWG output in FGEN mode and\n' ...
                 '  compare it with AWG mode. No further code will settle it.\n']);
    end
    fprintf('  Restored: mode %s, outputs off.\n\n', startMode);
end

function out = tryVariant(awg, freqMHz, dBm, startCmds)
    % Fresh start each time, then the variant's own ordering.
    step(awg, 'awgcontrol:stop:immediate');
    step(awg, 'OUTPUT1:STATE OFF');
    step(awg, 'OUTPUT2:STATE OFF');
    awg.setValue('mode', 'FGEN');
    step(awg, 'FGEN:CHAN1:TYPE SINE');
    step(awg, sprintf('FGEN:CHAN1:FREQ %.6f', freqMHz*1e6));
    step(awg, sprintf('FGEN:CHAN1:AMPL:POW %g', dBm));
    drain(awg, 10);

    out = struct('commands', {startCmds}, 'errors', {{}}, 'runStates', {{}});
    for c = 1:length(startCmds)
        e = step(awg, startCmds{c});
        out.errors{c} = e;
        fprintf('  %-42s error: %s\n', startCmds{c}, e);
    end

    % Read the run state immediately and again after a wait: entering run may
    % simply not be instantaneous, which would make every earlier reading of
    % "run state = 0" an artefact of asking too soon.
    delays = [0 0.2 1.0];
    for d = 1:length(delays)
        if delays(d) > 0; pause(delays(d)); end
        rs = strtrim(awg.queryRaw('AWGCONTROL:RSTATE?'));
        out.runStates{d} = rs;
        fprintf('  run state after %.1f s: %s\n', sum(delays(1:d)), rs);
    end
    out.finalRunState = out.runStates{end};
    out.outputState = strtrim(awg.queryRaw('OUTPUT1:STATE?'));
    out.running = ~strcmp(out.finalRunState, '0');
    fprintf('  output1 state: %s   -> %s\n', out.outputState, ...
            BooleanHelper.ifTrueElse(out.running, 'RUNNING', 'not running'));
end

function e = step(awg, cmd)
    awg.sendCommand(cmd);
    e = drain(awg, 3);
end

function txt = drain(awg, maxEntries)
    entries = {};
    for k = 1:maxEntries
        try
            r = strtrim(awg.queryRaw('SYSTEM:ERROR?', 5));
        catch
            break
        end
        if isempty(r); break; end
        if ~isempty(regexp(r, '^\s*[+-]?0\s*,', 'once')); break; end
        entries{end+1} = r;    %#ok<AGROW>
    end
    if isempty(entries)
        txt = 'none';
    else
        txt = strjoin(entries, ' | ');
    end
end

function restore(awg, startMode, startChannels)
    try
        awg.sendCommand('awgcontrol:stop:immediate');
        awg.sendCommand('OUTPUT1:STATE OFF');
        awg.sendCommand('OUTPUT2:STATE OFF');
        awg.setValue('mode', startMode);
        awg.workChannels = startChannels;
        awg.invalidateLoaded;
        drain(awg, 20);
    catch err
        warning('awgFgenBringUp could not fully restore the AWG: %s', err.message);
    end
end
