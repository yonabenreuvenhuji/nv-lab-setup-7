function info = awgFgenTone(varargin)
%AWGFGENTONE Put the AWG into function-generator mode and leave a waveform running.
%
%   awgFgenTone                                  % 1 kHz sine, 0.5 Vpp, channel 1
%   awgFgenTone(1e3)
%   awgFgenTone(1e3, 0.5, [1 2])                 % both channels, same settings
%   awgFgenTone(1e3, 0.5, 1, 'SQUARE')
%   awgFgenTone(struct('freq', [1e3 2e3], 'vpp', [0.5 0.2], ...
%                      'type', {{'SINE','SQUARE'}}, 'channels', [1 2]))
%
% A bench helper, not a calibration. Puts a continuous waveform on the AWG output
% so there is something to look at on a scope, and reports what the instrument says
% it is actually doing.
%
% ===================================================================
% THE ORDER OF OPERATIONS IS THE WHOLE POINT
% ===================================================================
%
% Enabling the output and then telling the generator to run DOES NOT WORK on this
% instrument, and it fails silently -- no SCPI error, the error queue stays clean,
% and the output sits there producing nothing. awgFgenBringUp tried four orderings:
%
%   outputs on -> RUN:IMMEDIATE     run state 0 at 0.0, 0.2 and 1.2 s   FAILS
%   RUN:IMMEDIATE -> outputs on     0, 0, then 2                        works
%   OUTP1 ON;:AWGCONTROL:RUN:IMM    0, 0, then 2                        works
%   outputs on -> RUN               0, 0, then 2                        works
%
% Note the second column: even the working orderings read 0 for about a second
% after the command. So checking the run state immediately tells you nothing, and
% this polls until it is non-zero rather than asking once. That pair of mistakes --
% wrong order, and believing an immediate readback -- is why the FGEN experiments
% produced a flat spectrum for four rounds.
%
% ===================================================================
% FINDING THE INSTRUMENT
% ===================================================================
%
% getObjByName('AWG') does NOT work: the object is named after
% FrequencyGeneratorTektronixAWG.TYPE plus its serial number, e.g.
% 'TektronixAWG-SN-NA'. So the lookup here mirrors what ExperimentAWG's constructor
% does -- walk FrequencyGenerator.getFG() and match on the name -- rather than
% guessing a name. It also accepts 'AWGdummy', for the same reason ExperimentAWG
% does.
%
% ===================================================================
% AMPLITUDE, AND WHY IT IS ASKED FOR IN VOLTS
% ===================================================================
%
% Given in VOLTS PEAK-TO-PEAK, because that is what a scope shows. Two commands can
% set it and they are not equivalent:
%
%   FGEN:CHANn:AMPLitude    volts peak-to-peak, shape-independent
%   FGEN:CHANn:AMPL:POW     dBm, which is a SINE-power convention
%
% This tries the volts form first and verifies by readback. That matters as soon as
% you ask for anything other than a sine: converting Vpp to dBm assumes a sine's
% crest factor, so a square wave requested as "0.5 Vpp" through the dBm path would
% come out at a different peak-to-peak than asked for. If the volts form is not
% accepted the code falls back to dBm and says so, so the discrepancy is visible
% rather than silent.
%
% For reference, into 50 ohms: 500 mVpp sine is -2.0 dBm, and the +10 dBm ceiling
% in setupInfo.json is 2 Vpp. Scope input and AWG output are both 50 ohms, so the
% scope reads the set amplitude with no factor of two.
%
% ===================================================================
% OPTIONS (struct form; every field may be scalar or per-channel)
%   channels    which AWG outputs (default 1)
%   freq        Hz
%   vpp         volts peak-to-peak
%   type        SINE SQUARE TRIANGLE RAMP PULSE DC GAUSSIAN SINC LORENTZ
%               ERISE EDECAY HAVERSINE   (char for all channels, or cellstr per)
%   phase       degrees (default 0)
%   offset      volts (default 0)
%   dutyCycle   percent, SQUARE and PULSE only
%   symmetry    percent, RAMP and TRIANGLE only
%
% AFTERWARDS, to hand the instrument back to the experiments:
%
%   awg = awgFgenTone('findonly');
%   awg.setValue('mode', 'AWG');
%   awg.invalidateLoaded;      % the upload cache is stale after a mode switch
%
% The invalidateLoaded is not optional: a mode switch stops playback, so the assets
% are no longer being played and the next experiment must re-upload rather than
% trust a cache hit that would otherwise look valid.

    % ---- 'findonly' shortcut, so the restore snippet needs no duplicate lookup
    if nargin == 1 && ischar(varargin{1}) && strcmpi(varargin{1}, 'findonly')
        info = findAwg();
        return
    end

    % ---- arguments: struct form or positional -------------------------------
    if nargin >= 1 && isstruct(varargin{1})
        o = varargin{1};
    else
        o = struct();
        if nargin >= 1 && ~isempty(varargin{1}); o.freq     = varargin{1}; end
        if nargin >= 2 && ~isempty(varargin{2}); o.vpp      = varargin{2}; end
        if nargin >= 3 && ~isempty(varargin{3}); o.channels = varargin{3}; end
        if nargin >= 4 && ~isempty(varargin{4}); o.type     = varargin{4}; end
    end
    o = withDefaults(o, struct('channels', 1, 'freq', 1e3, 'vpp', 0.5, ...
                               'type', 'SINE', 'phase', 0, 'offset', 0, ...
                               'dutyCycle', [], 'symmetry', []));

    ch = o.channels(:)';
    nC = numel(ch);
    freq   = perChannel(o.freq,   nC, 'freq');
    vpp    = perChannel(o.vpp,    nC, 'vpp');
    phase  = perChannel(o.phase,  nC, 'phase');
    offset = perChannel(o.offset, nC, 'offset');
    duty   = perChannel(o.dutyCycle, nC, 'dutyCycle');
    symm   = perChannel(o.symmetry,  nC, 'symmetry');
    types  = perChannelStr(o.type, nC);

    awg = findAwg();

    fprintf('\n=== awgFgenTone ===\n');
    fprintf('instrument: %s\n', awg.name);
    for i = 1:nC
        fprintf('  ch%d: %s, %g Hz, %g Vpp, phase %g deg, offset %g V\n', ...
                ch(i), types{i}, freq(i), vpp(i), phase(i), offset(i));
    end

    % ---- mode --------------------------------------------------------------
    current = strtrim(awg.queryRaw('INST:MODE?'));
    if ~strcmpi(current, 'FGEN')
        fprintf('instrument mode: %s -> FGEN\n', current);
        awg.setValue('mode', 'FGEN');
        try; awg.invalidateLoaded; catch; end
    else
        fprintf('instrument mode: already FGEN\n');
    end

    % ---- configure, BEFORE running ----------------------------------------
    for i = 1:nC
        c = ch(i);
        cmds = { sprintf('FGEN:CHAN%d:TYPE %s', c, types{i}), ...
                 sprintf('FGEN:CHAN%d:FREQ %.9g', c, freq(i)), ...
                 sprintf('FGEN:CHAN%d:OFFSET %g', c, offset(i)), ...
                 sprintf('FGEN:CHAN%d:PHASE %g', c, phase(i)) };
        % Shape parameters only where they mean something. Sending a duty cycle to a
        % sine is the kind of command this instrument rejects into the error queue
        % and nowhere else.
        if ~isnan(duty(i)) && any(strcmpi(types{i}, {'SQUARE','SQU','PULSE','PULS'}))
            cmds{end+1} = sprintf('FGEN:CHAN%d:DCYCLE %g', c, duty(i));           %#ok<AGROW>
        end
        if ~isnan(symm(i)) && any(strcmpi(types{i}, {'RAMP','TRIANGLE','TRI'}))
            cmds{end+1} = sprintf('FGEN:CHAN%d:SYMMETRY %g', c, symm(i));         %#ok<AGROW>
        end
        awg.sendCommandBatch(cmds);
        setAmplitude(awg, c, vpp(i));
    end

    % ---- RUN, then WAIT, then outputs -------------------------------------
    awg.sendCommand('awgcontrol:run:immediate');
    t0 = tic; runState = '0';
    while toc(t0) < 10
        runState = strtrim(awg.queryRaw('AWGCONTROL:RSTATE?'));
        if ~strcmp(runState, '0'); break; end
        pause(0.05);
    end
    fprintf('run state reached %s after %.1f s\n', runState, toc(t0));
    if strcmp(runState, '0')
        drainErrors(awg);
        error(['The generator never started: run state still 0 after 10 s. Nothing ' ...
               'will appear on the scope. Run awgFgenBringUp to see which start ' ...
               'ordering this firmware accepts.']);
    end

    awg.workChannels = ch;              % outputs on, AFTER the run

    % ---- readback ----------------------------------------------------------
    % Ask the instrument what it thinks it is doing rather than assuming the commands
    % took. FGEN silently coerces frequencies and amplitudes it cannot produce, and a
    % coerced value is far easier to spot here than on a scope.
    info = struct('awg', awg, 'runState', runState, 'channels', ch);
    fprintf('\nreadback:\n');
    for i = 1:nC
        c = ch(i);
        f  = str2double(awg.queryRaw(sprintf('FGEN:CHAN%d:FREQ?', c)));
        ty = strtrim(awg.queryRaw(sprintf('FGEN:CHAN%d:TYPE?', c)));
        ou = strtrim(awg.queryRaw(sprintf('OUTPUT%d:STATE?', c)));
        vb = readAmplitude(awg, c);
        fprintf('  ch%d: output %s, %s, %.9g Hz, %.4g Vpp\n', c, ou, ty, f, vb);
        info.(sprintf('ch%d', c)) = struct('freqHz', f, 'type', ty, 'vpp', vb, ...
                                           'output', ou);
        if abs(f - freq(i)) > max(1e-6*freq(i), 1e-3)
            fprintf('       COERCED frequency: asked %g Hz, got %.9g Hz\n', freq(i), f);
        end
        if isfinite(vb) && abs(vb - vpp(i)) > 0.1*vpp(i)
            fprintf('       COERCED amplitude: asked %g Vpp, got %.4g Vpp\n', vpp(i), vb);
        end
        if ~any(strcmp(ou, {'1', 'ON'}))
            fprintf('       output reads back OFF -- nothing is on the cable\n');
        end
    end

    drainErrors(awg);

    % ---- what to expect on the scope --------------------------------------
    % Spelled out rather than left to be worked out, because the commonest way to see
    % nothing on a 4 GHz scope is a timebase that is right for the instrument and
    % wrong for the signal by six orders of magnitude.
    period = 1/max(freq);
    % The EXR404A inputs are SWITCHABLE 50 ohm / 1 Mohm, unlike the DSO90404A this text
    % was written for, which had no choice to get wrong. The 1 Mohm path rolls off a few
    % hundred MHz up, so a GHz tone into it reads as very nearly nothing -- worth saying
    % out loud every time, because it is the failure that looks like a dead cable.
    fprintf('\nOn the scope:\n');
    fprintf('  input      50 ohm, NOT 1 Mohm (the EXR can do either; 1 Mohm kills a GHz tone)\n');
    fprintf('  timebase   ~%s/div  (fastest period is %s, so that shows ~2 cycles)\n', ...
            engStr(period/2), engStr(period));
    fprintf('  vertical   ~%.0f mV/div, offset 0\n', 1000*max(vpp)/4);
    fprintf('  trigger    Auto sweep, edge on that channel at 0 V\n');
    fprintf('  or Default Setup then Autoscale, which fixes all of the above at once.\n');
    fprintf('\nTo give the instrument back to the experiments:\n');
    fprintf('  awg = awgFgenTone(''findonly''); awg.setValue(''mode'', ''AWG''); awg.invalidateLoaded;\n\n');
end

% ======================================================================

function awg = findAwg()
% Mirror what ExperimentAWG's constructor does. getObjByName('AWG') fails because
% the object is named TYPE plus serial number, e.g. 'TektronixAWG-SN-NA'.
    awg = [];
    try
        fgs = FrequencyGenerator.getFG();
        for i = 1:numel(fgs)
            if contains(fgs{i}.name, 'TektronixAWG') || contains(fgs{i}.name, 'AWGdummy')
                awg = fgs{i}; return
            end
        end
    catch
    end
    % Last resort, in case the registry is not populated but the object exists.
    for nm = {'TektronixAWG', 'AWG'}
        try
            candidate = getObjByName(nm{1});
            if ~isempty(candidate); awg = candidate; return; end
        catch
        end
    end
    error(['No AWG found. FrequencyGenerator.getFG() returned nothing matching ' ...
           '''TektronixAWG'' or ''AWGdummy''. Has the setup been initialised?']);
end

function setAmplitude(awg, c, vpp)
% Prefer the volts command and verify it. Fall back to dBm, which assumes a sine
% crest factor and so is wrong for any other shape -- but say so rather than
% quietly delivering a different amplitude than was asked for.
    awg.sendCommand(sprintf('FGEN:CHAN%d:AMPLITUDE %.6g', c, vpp));
    back = readAmplitude(awg, c);
    if isfinite(back) && abs(back - vpp) <= 0.05*vpp
        return
    end
    dbm = 10*log10( (vpp/(2*sqrt(2)))^2 / 50 / 1e-3 );
    awg.sendCommand(sprintf('FGEN:CHAN%d:AMPL:POW %g', c, dbm));
    fprintf(['  ch%d: FGEN:AMPLITUDE (volts) not accepted, used AMPL:POW %.2f dBm\n' ...
             '        instead. That conversion assumes a SINE, so for another shape the\n' ...
             '        peak-to-peak will not be the %g V asked for.\n'], c, dbm, vpp);
end

function vpp = readAmplitude(awg, c)
    vpp = NaN;
    try
        v = str2double(awg.queryRaw(sprintf('FGEN:CHAN%d:AMPLITUDE?', c), 2));
        if isfinite(v) && v > 0; vpp = v; return; end
    catch
    end
    try
        a = str2double(awg.queryRaw(sprintf('FGEN:CHAN%d:AMPL:POW?', c), 2));
        if isfinite(a); vpp = 2*sqrt(2)*sqrt(10^(a/10)*1e-3*50); end
    catch
    end
end

function drainErrors(awg)
% Empty the SCPI error queue and print anything in it. Execution errors on this
% instrument are invisible to fwrite -- they sit in the queue and nowhere else -- so
% a rejected command leaves no trace at the MATLAB prompt.
    for k = 1:20
        r = strtrim(awg.queryRaw('SYSTEM:ERROR?'));
        if isempty(r) || startsWith(r, '0,'); break; end
        fprintf('  SCPI error: %s\n', r);
    end
end

function v = perChannel(x, nC, name)
    if isempty(x); v = nan(1, nC); return; end
    x = x(:)';
    if isscalar(x); v = repmat(x, 1, nC); return; end
    if numel(x) ~= nC
        error('%s has %d entries for %d channels.', name, numel(x), nC);
    end
    v = x;
end

function c = perChannelStr(t, nC)
    if ischar(t) || isstring(t)
        c = repmat({upper(char(t))}, 1, nC);
    elseif iscell(t)
        if numel(t) == 1; t = repmat(t, 1, nC); end
        if numel(t) ~= nC
            error('type has %d entries for %d channels.', numel(t), nC);
        end
        c = cellfun(@(s) upper(char(s)), t, 'UniformOutput', false);
    else
        error('type must be a char or a cell array of chars.');
    end
    known = {'SINE','SIN','SQUARE','SQU','TRIANGLE','TRI','RAMP','PULSE','PULS', ...
             'DC','GAUSSIAN','GAUS','SINC','LORENTZ','LOR','ERISE','EDECAY', ...
             'HAVERSINE','HAV'};
    for i = 1:numel(c)
        if ~any(strcmpi(c{i}, known))
            error(['Unknown waveform type "%s". Known: %s. (Unrecognised types are ' ...
                   'rejected into the SCPI error queue and nowhere else, so they are ' ...
                   'checked here instead.)'], c{i}, strjoin(known, ' '));
        end
    end
end

function o = withDefaults(o, d)
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(o, f{i}); o.(f{i}) = d.(f{i}); end
    end
end

function s = engStr(x)
    if x >= 1;        s = sprintf('%.3g s', x);
    elseif x >= 1e-3; s = sprintf('%.3g ms', x*1e3);
    elseif x >= 1e-6; s = sprintf('%.3g us', x*1e6);
    else;             s = sprintf('%.3g ns', x*1e9);
    end
end
