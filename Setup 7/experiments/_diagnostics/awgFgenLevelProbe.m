function results = awgFgenLevelProbe(awg, dBm)
%AWGFGENLEVELPROBE Does FGEN mode emit the same power as AWG mode?
%
%   results = awgFgenLevelProbe(awg)
%   results = awgFgenLevelProbe(awg, -20)
%
% WHY. An FGEN ESR run gave 0.07% contrast where AWG mode gives 3.9%. The setup
% json settles what that is NOT: two microwave switches are declared on
% pulseStreamer channels 5 and 6 ("MW" and "MW2"), the AWG is the only frequency
% generator on this config, and those channels are asserted during the signal
% window -- so the gate exists and is driven correctly. Zero contrast therefore
% is not a gating failure.
%
% What is left is LEVEL, and that is my code's assumption. The two modes set
% amplitude by completely different routes:
%
%   AWG mode   setAmplitude(P) puts the ANALOG output at MAX_VPP (full scale) and
%              encodes P as MaxNormPowerAllowed, a normalised scaling applied to
%              the waveform SAMPLES. Emitted Vpp = dBmToVpp(P).
%
%   FGEN mode  there are no samples to scale, so SetupFgen sends
%              FGEN:CHAN<n>:AMPL:POW P and assumes that means the same thing in
%              dBm. That assumption was never tested. If the node is in volts, or
%              if FGEN takes its level from somewhere else entirely, the emitted
%              power can be far below what the NVs need -- which looks exactly
%              like no microwave at all.
%
% This reads what each mode reports for the output level, in both modes, with the
% OUTPUTS OFF the whole time. Nothing reaches the sample; only settings are read.
% Restores the starting mode and channel state.
%
% WHAT TO DO WITH THE ANSWER. Compare the last two lines. If AWG mode's analog
% amplitude is full scale (MAX_VPP) while FGEN's is a small fraction of it, the
% level is the bug and it is fixable in software. If they match, the level is not
% the bug and the next step is a power measurement at the connector.

    if ~exist('dBm', 'var') || isempty(dBm); dBm = -20; end

    fprintf('\n=== awgFgenLevelProbe (test level %g dBm) ===\n', dBm);
    results = struct('askedDBm', dBm);

    awg.sendCommand('*CLS');
    startMode = strtrim(awg.queryRaw('INST:MODE?'));
    startChannels = awg.workChannels;
    fprintf('starting mode %s, MAX_VPP = %g V\n', startMode, awg.MAX_VPP);
    results.startMode = startMode;
    results.maxVpp = awg.MAX_VPP;

    try
        awg.workChannels = [];                          % outputs off throughout
        awg.sendCommand('awgcontrol:stop:immediate');

        % ---- AWG mode -------------------------------------------------------
        awg.setValue('mode', 'AWG');
        awg.setAmplitude(dBm * ones(1, 2));
        results.awgAnalogVpp   = str2double(awg.queryRaw('SOURCE1:VOLTAGE:AMPLITUDE?'));
        results.maxNormPower   = awg.MaxNormPowerAllowed;
        results.impliedAwgVpp  = results.awgAnalogVpp * results.maxNormPower(1);
        fprintf(['\nAWG mode\n' ...
                 '  analog amplitude          %.4f Vpp   (SOURCE1:VOLTAGE:AMPLITUDE?)\n' ...
                 '  MaxNormPowerAllowed       %s          (waveform sample scaling)\n' ...
                 '  => microwave actually at  %.5f Vpp  (product of the two)\n'], ...
                results.awgAnalogVpp, mat2str(results.maxNormPower, 4), results.impliedAwgVpp);

        % ---- FGEN mode ------------------------------------------------------
        awg.setValue('mode', 'FGEN');
        awg.sendCommand('FGEN:CHAN1:TYPE SINE');
        awg.sendCommand(sprintf('FGEN:CHAN1:AMPL:POW %g', dBm));
        results.fgenAmplPow = str2double(awg.queryRaw('FGEN:CHAN1:AMPL:POW?'));

        % Ask for the level through every node that might report it. Unknown
        % mnemonics simply error into the queue; a NaN here means "not supported",
        % which is itself informative.
        candidates = {'FGEN:CHAN1:AMPL?', 'FGEN:CHAN1:AMPLITUDE?', ...
                      'FGEN:CHAN1:HIGH?', 'FGEN:CHAN1:LOW?', ...
                      'FGEN:CHAN1:OFFSET?', 'SOURCE1:VOLTAGE:AMPLITUDE?'};
        results.fgenReadings = struct('command', {}, 'value', {});
        fprintf('\nFGEN mode\n  AMPL:POW readback         %.4f (the value I set)\n', ...
                results.fgenAmplPow);
        for c = 1:length(candidates)
            v = NaN;
            try
                v = str2double(awg.queryRaw(candidates{c}, 5));
            catch
                % no such node on this firmware
            end
            drainErrors(awg, 3);
            results.fgenReadings(end+1) = struct('command', candidates{c}, 'value', v);  %#ok<AGROW>
            if ~isnan(v)
                fprintf('  %-28s %.6f\n', candidates{c}, v);
            else
                fprintf('  %-28s (not supported)\n', candidates{c});
            end
        end

        % Whichever of those is a voltage lets us compare like with like.
        fgenVpp = NaN;
        for c = 1:length(results.fgenReadings)
            nm = results.fgenReadings(c).command;
            if any(strcmp(nm, {'FGEN:CHAN1:AMPL?', 'FGEN:CHAN1:AMPLITUDE?'})) && ...
                    ~isnan(results.fgenReadings(c).value)
                fgenVpp = results.fgenReadings(c).value;
                break
            end
        end
        results.fgenVpp = fgenVpp;

        fprintf('\n--- comparison ---\n');
        fprintf('  AWG mode delivers   %.5f Vpp for %g dBm\n', results.impliedAwgVpp, dBm);
        if isnan(fgenVpp)
            fprintf(['  FGEN mode reports no voltage node, so the two cannot be compared\n' ...
                     '  directly from here. Read the amplitude off the front panel in FGEN\n' ...
                     '  mode instead and compare it with the number above.\n']);
        else
            fprintf('  FGEN mode delivers  %.5f Vpp for the same request\n', fgenVpp);
            ratio = fgenVpp / max(results.impliedAwgVpp, eps);
            results.vppRatio = ratio;
            fprintf('  ratio %.3f  (%.1f dB)\n', ratio, 20*log10(max(ratio, eps)));
            if ratio < 0.5 || ratio > 2
                fprintf(['  THAT IS THE BUG. The FGEN amplitude is not the same power as AWG\n' ...
                         '  mode, so the NVs were driven far off resonance-saturation and the\n' ...
                         '  contrast vanished. Fixable: set the FGEN level to match, in\n' ...
                         '  whatever unit this node actually wants.\n']);
            else
                fprintf(['  The levels agree, so amplitude is NOT the explanation and the\n' ...
                         '  next step is measuring actual power at the connector.\n']);
            end
        end

    catch err
        restore(awg, startMode, startChannels, dBm);
        fprintf('\n  probe aborted: %s\n', err.message);
        rethrow(err);
    end

    restore(awg, startMode, startChannels, dBm);
    fprintf('\n  restored: mode %s, amplitude re-applied, outputs as before.\n\n', startMode);
end

function restore(awg, startMode, startChannels, dBm)
    try
        awg.setValue('mode', startMode);
        awg.setAmplitude(dBm * ones(1, 2));   % leave MaxNormPowerAllowed consistent
        awg.workChannels = startChannels;
        awg.invalidateLoaded;                 % mode changes stopped playback
        drainErrors(awg, 20);
    catch err
        warning('awgFgenLevelProbe could not fully restore the AWG: %s', err.message);
    end
end

function txt = drainErrors(awg, maxEntries)
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
        txt = 'no error';
    else
        txt = strjoin(entries, ' | ');
    end
end
