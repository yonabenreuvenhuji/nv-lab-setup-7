function out = scopeProbe(freqHz, opts)
%SCOPEPROBE One frequency, one acquisition, and a full account of what came back.
%
%   out = scopeProbe                      % 2.5 GHz, channel 1
%   out = scopeProbe(2.87e9)
%   out = scopeProbe(2.5e9, struct('autoscale', true, 'vpp', 0.1))
%
% A diagnostic, not a measurement. Turns on the microwave at one frequency, takes a
% single acquisition, and prints every quantity needed to tell whether the trace
% actually contains the tone -- plus a plot of the trace and its full spectrum.
%
% ===================================================================
% WHY THIS EXISTS
% ===================================================================
%
% calibResonator's first run reported 12 to 107 MICROvolts where the same signal
% read 141 mV peak on the scope face: a factor of about 4000. Everything downstream
% was then garbage -- the auto-range ratcheted to its 16 mV floor because the
% amplitude looked like 0.2% of full scale, the linearity check compared noise with
% noise, and the flattening table produced +155 dB entries.
%
% One number would have caught it immediately, and calibResonator never asked for it:
% the SAMPLE RATE. A Fourier amplitude evaluated at 2.5 GHz is only meaningful if
% 2.5 GHz is below Nyquist. If the scope is acquiring at, say, 2 GS/s -- because
% memory depth was pinned, or it was left in segmented mode, or the timebase forced
% a decimation -- then the tone aliases somewhere else entirely and the DFT at the
% requested frequency correctly returns approximately nothing.
%
% So this prints the sample rate, the samples per cycle, and where the energy in the
% spectrum actually is. If the peak is not at the drive frequency, that is the answer.
%
% ===================================================================
% WHAT IT CHECKS
% ===================================================================
%
%   sample rate and samples per cycle    below ~4 samples per cycle the amplitude
%                                        estimate degrades; below 2 it is aliased
%   raw ADC extremes                     near full scale means clipping, which
%                                        flattens the sine and biases the amplitude
%   volts from two independent paths      the binary WORD transfer and an ASCII
%                                        transfer of the same trace. If those two
%                                        disagree the fault is in the preamble
%                                        scaling, not the acquisition
%   spectral peak                        where the energy is, against where it
%                                        was asked to be
%   Vrms and Vpp                         directly comparable with what the scope's
%                                        own Measure menu shows, which is the one
%                                        number already known to be right
%
% ===================================================================
% OPTIONS (struct fields, all optional)
%   channel       AWG channel (default 1)
%   scopeChannel  scope input (default = channel)
%   scopeResource VISA resource (default: the EXR/MXR/DSO found by visadevlist)
%   vpp           AWG amplitude (default 0.1)
%   recordTime    seconds per acquisition (default 200e-9)
%   range         scope full-scale volts, [] to leave alone (default 0.5)
%   autoscale     issue :AUToscale before acquiring (default false)
%   settle        seconds to wait after setting the tone (default 1)
%   leaveOn       leave microwave on when finished (default false)
%   mwChannel     pulse-generator channel (default 'MW' / 'MW2')
%
% Writes nothing. Restores the bench unless leaveOn is true.

    if ~exist('freqHz', 'var') || isempty(freqHz); freqHz = 2.5e9; end
    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('channel', 1, 'scopeChannel', [], ...
        'scopeResource', '', 'vpp', 0.1, 'recordTime', 200e-9, 'range', 0.5, ...
        'autoscale', false, 'settle', 1, 'leaveOn', false, 'mwChannel', ''));
    if isempty(opts.scopeChannel); opts.scopeChannel = opts.channel; end
    if isempty(opts.mwChannel)
        if opts.channel == 1; opts.mwChannel = 'MW'; else; opts.mwChannel = 'MW2'; end
    end
    ch = opts.scopeChannel;

    fprintf('\n=== scopeProbe ===\n');
    fprintf('%.4f GHz, %g Vpp, AWG ch%d -> scope ch%d\n', ...
            freqHz/1e9, opts.vpp, opts.channel, ch);

    awg = findAwg();
    sc  = openScope(opts.scopeResource);
    pg  = getObjByName(PulseGenerator.NAME);

    % ---- microwave on ------------------------------------------------------
    pg.on(opts.mwChannel);
    awgFgenTone(struct('channels', opts.channel, 'freq', freqHz, ...
                       'vpp', opts.vpp, 'type', 'SINE'));
    fprintf('\nsettling %g s\n', opts.settle);
    pause(opts.settle);

    try
        % ---- acquisition setup ---------------------------------------------
        % Every setup command is sent CHECKED: written, then the error queue is read
        % immediately so a rejected mnemonic is attributed to the command that caused
        % it. The first EXR run ended with five -113 "Undefined header" and two -420
        % "Query UNTERMINATED" in one undifferentiated heap at the end of the output,
        % which says something was refused but not what. Draining per command turns
        % that heap into a list of names.
        writeline(sc, ':SYSTem:HEADer OFF');
        drainErrors(sc);
        refused = {};
        setup = { ':ACQuire:MODE RTIMe'
                  ':ACQuire:AVERage OFF'
                  ':ACQuire:SRATe:ANALog MAX'
                  ':ACQuire:POINts:ANALog AUTO'
                  % Interpolation OFF. The first EXR run reported 16 GSa/s from
                  % :ACQuire:SRATe:ANALog? and 32 GSa/s from the preamble's Xincrement,
                  % with 3200 acquired points returned as 6399 -- the scope was
                  % sinc-interpolating 2x and scopeProbe, unlike calibResonator, never
                  % asked it to stop. Interpolated points carry no information: they
                  % double the transfer time, they make "samples per cycle" report 11.31
                  % where the hardware delivered 5.65, and they hide how much real
                  % Nyquist margin there is.
                  ':ACQuire:INTerpolate OFF'
                  ':TRIGger:SWEep AUTO'
                  sprintf(':TIMebase:SCALe %g', opts.recordTime/10)
                  ':TIMebase:POSition 0'
                  sprintf(':CHANnel%d:DISPlay ON', ch) };
        % :WAVeform:UNSigned is NOT sent. The EXR rejects it with -113 "Undefined
        % header" -- it has no such setting -- and it turns out not to be needed: the
        % ASCII cross-check returned 70.66 mV against the WORD path's 70.66 mV, four
        % significant figures identical. That is precisely the test that proves the
        % int16 interpretation and the preamble scaling are right, so the signedness
        % question is settled by measurement rather than by a command.
        for i = 1:numel(setup)
            if ~sendChecked(sc, setup{i}); refused{end+1} = setup{i}; end  %#ok<AGROW>
        end
        setInputFifty(sc, ch);
        if ~isempty(opts.range)
            if ~sendChecked(sc, sprintf(':CHANnel%d:RANGe %g', ch, opts.range))
                refused{end+1} = 'RANGe';
            end
            sendChecked(sc, sprintf(':CHANnel%d:OFFSet 0', ch));
        end
        if ~isempty(refused)
            fprintf('\nthe EXR refused these (harmless unless noted, but now named):\n');
            for i = 1:numel(refused); fprintf('  %s\n', refused{i}); end
        end
        out.refused = refused;
        if opts.autoscale
            fprintf('issuing :AUToscale\n');
            writeline(sc, ':AUToscale');
            pause(3);
        end

        % ---- what the scope says its state is -------------------------------
        fprintf('\n--- scope state ---\n');
        st = struct();
        items = { 'acq mode',      ':ACQuire:MODE?'
                  'sample rate',   ':ACQuire:SRATe:ANALog?'
                  'points',        ':ACQuire:POINts:ANALog?'
                  'averaging',     ':ACQuire:AVERage?'
                  'timebase s/div',':TIMebase:SCALe?'
                  'trigger sweep', ':TRIGger:SWEep?'
                  'ch range',      sprintf(':CHANnel%d:RANGe?', ch)
                  'ch offset',     sprintf(':CHANnel%d:OFFSet?', ch)
                  'ch input',      sprintf(':CHANnel%d:INPut?', ch)
                  'ch display',    sprintf(':CHANnel%d:DISPlay?', ch) };
        % Every state query goes through queryOpt, not query1. This is a diagnostic
        % whose entire job is to survive an instrument it does not fully know: one
        % unsupported mnemonic must print "(not supported)" and move on, not stall for
        % the full VISA timeout and then abort before any waveform is read. The EXR
        % rejects :CHANnel<n>:IMPedance? -- it spells the same setting :INPut DC50 --
        % and that one query used to kill the whole run.
        for i = 1:size(items,1)
            v = queryOpt(sc, items{i,2});
            fprintf('  %-16s %s\n', items{i,1}, v);
            st.(matlab.lang.makeValidName(items{i,1})) = v;
        end
        out.state = st;

        % ---- acquire, binary -----------------------------------------------
        writeline(sc, sprintf(':DIGitize CHANnel%d', ch));
        writeline(sc, sprintf(':WAVeform:SOURce CHANnel%d', ch));
        writeline(sc, ':WAVeform:FORMat WORD');
        writeline(sc, ':WAVeform:BYTeorder LSBFirst');
        if ~sendChecked(sc, ':WAVeform:STReaming OFF')
            out.refused{end+1} = ':WAVeform:STReaming OFF';
        end

        preStr = query1(sc, ':WAVeform:PREamble?');
        pre = str2double(strsplit(strtrim(preStr), ','));
        fprintf('\n--- preamble ---\n');
        names = {'format','type','points','count','Xinc','Xorg','Xref', ...
                 'Yinc','Yorg','Yref'};
        for i = 1:min(10, numel(pre))
            fprintf('  %-8s %g\n', names{i}, pre(i));
        end
        out.preamble = pre;

        writeline(sc, ':WAVeform:DATA?');
        raw = double(readbinblock(sc, 'int16'));
        try; readline(sc); catch; end        % swallow the terminator

        xinc = pre(5); xorg = pre(6); xref = pre(7);
        yinc = pre(8); yorg = pre(9); yref = pre(10);
        v = (raw - yref)*yinc + yorg;
        n = numel(v);
        t = ((0:n-1) - xref)*xinc + xorg;

        % ---- the numbers that matter ---------------------------------------
        sr  = 1/xinc;
        spc = sr/freqHz;
        fprintf('\n--- acquisition ---\n');
        fprintf('  points returned      %d\n', n);
        fprintf('  Xincrement           %.4g s  ->  sample rate %.4g GS/s\n', xinc, sr/1e9);
        % The preamble rate is what the returned array actually implies. The SRATe query
        % is what the hardware digitised at. They differ if the scope is interpolating,
        % and it is the HARDWARE number that sets the real Nyquist margin -- so print
        % both and say which is which rather than quietly believing the flattering one.
        srHw = NaN;
        try; srHw = str2double(queryOpt(sc, ':ACQuire:SRATe:ANALog?')); catch; end
        if isfinite(srHw) && srHw > 0
            fprintf('  hardware sample rate %.4g GS/s (:ACQuire:SRATe:ANALog?)\n', srHw/1e9);
            if abs(srHw - sr)/srHw > 0.05
                fprintf(['     interpolating %.2fx: %d returned points are not %d\n' ...
                         '     independent measurements. Real samples per cycle is %.2f,\n' ...
                         '     not %.2f. Harmless for a single-tone amplitude, but the\n' ...
                         '     Nyquist margin is the smaller number.\n'], ...
                        sr/srHw, n, round(n*srHw/sr), srHw/freqHz, spc);
            end
        end
        out.hardwareSampleRate = srHw;
        fprintf('  record length        %.4g ns  (%.0f cycles at %.4f GHz)\n', ...
                n*xinc*1e9, n*xinc*freqHz, freqHz/1e9);
        fprintf('  samples per cycle    %.2f\n', spc);
        if spc < 2
            fprintf(['     *** ALIASED. Below two samples per cycle the tone is folded to a\n' ...
                     '     different frequency, and a DFT at %.4f GHz will read essentially\n' ...
                     '     zero however strong the real signal is. THIS IS THE BUG if the\n' ...
                     '     spectrum below peaks somewhere else.\n'], freqHz/1e9);
        elseif spc < 4
            fprintf('     marginal: under 4 samples per cycle the amplitude estimate degrades.\n');
        end
        fprintf('  raw ADC min/max      %d / %d\n', min(raw), max(raw));
        if max(abs(raw)) > 30000
            fprintf('     *** near full scale: the sine is clipping, raise the range\n');
        end
        fprintf('  volts min/max        %+.4g / %+.4g V\n', min(v), max(v));
        fprintf('  Vpp (min to max)     %.4g mV\n', 1e3*(max(v)-min(v)));
        fprintf('  Vrms                 %.4g mV   <-- compare with the scope''s Measure menu\n', ...
                1e3*sqrt(mean(v.^2)));

        % ---- amplitude at the drive frequency, and where the energy really is
        w = 0.5*(1 - cos(2*pi*(0:n-1)'/(n-1)));
        A = 2*abs(sum(v(:).*w.*exp(-2i*pi*freqHz*t(:))))/sum(w);
        fprintf('  DFT amplitude at %.4f GHz   %.4g mV pk\n', freqHz/1e9, 1e3*A);

        V = abs(fft(v(:).*w));
        V = V(1:floor(n/2));
        fAx = (0:numel(V)-1)' * (sr/n);
        [~, ip] = max(V(2:end)); ip = ip + 1;
        fPeak = fAx(ip);
        Apeak = 2*V(ip)/sum(w);
        fprintf('  spectral peak at     %.4f GHz, %.4g mV pk\n', fPeak/1e9, 1e3*Apeak);
        if abs(fPeak - freqHz) > 0.02*freqHz
            fprintf(['     *** THE ENERGY IS NOT AT THE DRIVE FREQUENCY. Asked for %.4f GHz,\n' ...
                     '     found %.4f GHz. With %.2f samples per cycle, an alias would appear\n' ...
                     '     at %.4f GHz -- check whether that matches.\n'], ...
                    freqHz/1e9, fPeak/1e9, spc, abs(freqHz - sr*round(freqHz/sr))/1e9);
        else
            fprintf('     the peak is at the drive frequency, so the trace is genuine.\n');
        end

        out.freqHz = freqHz; out.t = t(:); out.v = v(:); out.raw = raw(:);
        out.sampleRate = sr; out.samplesPerCycle = spc;
        out.amplitudeAtDrive = A; out.peakFreq = fPeak; out.peakAmplitude = Apeak;
        out.vrms = sqrt(mean(v.^2));

        % ---- independent cross-check: the same trace as ASCII ---------------
        % If the ASCII and WORD amplitudes disagree, the fault is in the preamble
        % scaling rather than in the acquisition -- a different problem with a
        % different fix, and worth separating rather than guessing.
        try
            % Read the ASCII block with readbinblock-style patience: a single readline
            % cannot swallow 32000 comma-separated values, and whatever it leaves
            % pending makes the NEXT query fail with -410 "Query INTERRUPTED". That is
            % where the -410 in the first run came from -- confined to this
            % cross-check, but worth not emitting at all.
            % NumBytesAvailable is gone from visadev in current MATLAB ("no longer
            % supported"), which is what skipped this check on the EXR -- a MATLAB API
            % change, nothing to do with the scope. The loop now terminates on having
            % collected enough VALUES rather than on an emptied input buffer, which is
            % the honest condition anyway: we know how many points the preamble promised.
            writeline(sc, ':WAVeform:FORMat ASCii');
            writeline(sc, ':WAVeform:DATA?');
            a = "";
            tRead = tic;
            for k = 1:200
                if toc(tRead) > 20; break; end
                try
                    chunk = readline(sc);
                catch
                    break                      % nothing more coming
                end
                a = a + chunk;
                if count(a, ',') + 1 >= n; break; end
            end
            writeline(sc, ':WAVeform:FORMat WORD');   % leave it as we found it
            va = str2double(strsplit(strtrim(char(a)), ','));
            va = va(isfinite(va));
            if numel(va) > 16
                na = numel(va);
                wa = 0.5*(1 - cos(2*pi*(0:na-1)'/(na-1)));
                ta = (0:na-1)'*xinc;
                Aa = 2*abs(sum(va(:).*wa.*exp(-2i*pi*freqHz*ta)))/sum(wa);
                fprintf('\n  ASCII cross-check: %d points, amplitude %.4g mV pk\n', na, 1e3*Aa);
                out.asciiAmplitude = Aa;
                if A > 0 && (Aa/A > 1.2 || Aa/A < 0.8)
                    fprintf(['     ASCII and WORD disagree by %.2fx, so the binary scaling from\n' ...
                             '     the preamble is wrong, not the acquisition.\n'], Aa/A);
                end
            end
        catch e
            fprintf('\n  ASCII cross-check skipped: %s\n', e.message);
        end

        drainScopeErrors(sc);

        % ---- plots ----------------------------------------------------------
        figure
        subplot(2,1,1)
        nShow = min(n, round(20*max(spc,2)));       % ~20 cycles, so a sine is visible
        plot(t(1:nShow)*1e9, 1e3*v(1:nShow), '.-');
        xlabel('time (ns)'); ylabel('mV');
        title(sprintf('trace, first %d points at %.3g GS/s (%.2f samples/cycle)', ...
              nShow, sr/1e9, spc));
        subplot(2,1,2)
        plot(fAx/1e9, 1e3*2*V/sum(w));
        hold on; line([freqHz freqHz]/1e9, ylim, 'Color', 'r');
        xlabel('GHz'); ylabel('mV pk'); set(gca, 'YScale', 'log');
        title('spectrum -- red is the drive frequency; the peak should sit on it');

    catch err
        cleanUp(awg, pg, opts);
        rethrow(err);
    end

    if ~opts.leaveOn
        cleanUp(awg, pg, opts);
    else
        fprintf('\n(microwave LEFT ON at %.4f GHz -- call pg.off when done)\n', freqHz/1e9);
    end
end

% ======================================================================

function cleanUp(awg, pg, ~)
    try; pg.off; catch; end
    try; awg.workChannels = []; catch; end
    try; awg.setValue('mode', 'AWG'); awg.invalidateLoaded; catch; end
    fprintf('\n(bench restored: switch off, AWG outputs off, AWG mode)\n');
end

function trySend(sc, cmd)
% Some of these mnemonics vary between Infiniium firmware versions. A rejected
% command lands in the error queue and nowhere else, so send it, then clear.
    try
        writeline(sc, cmd);
    catch
        fprintf('  (command not accepted: %s)\n', cmd);
    end
end

function drainScopeErrors(sc)
    for k = 1:20
        r = strtrim(query1(sc, ':SYSTem:ERRor? STRing'));
        if isempty(r) || startsWith(r, '0') || startsWith(r, '+0'); break; end
        fprintf('  scope error: %s\n', r);
    end
end

function s = query1(sc, cmd)
    writeline(sc, cmd);
    s = char(readline(sc));
end

function ok = sendChecked(sc, cmd)
% Write a command and immediately ask whether the instrument objected. Returns false if
% it did, and prints the error alongside the command that produced it. The point is
% attribution: an error queue drained at the end of a run tells you something was
% refused, not which line to fix.
    ok = true;
    try
        writeline(sc, cmd);
    catch
        fprintf('  (could not send: %s)\n', cmd);
        ok = false; return
    end
    e = queryOpt(sc, ':SYSTem:ERRor? STRing');
    if ~isempty(e) && ~startsWith(e, '0') && ~startsWith(e, '+0') && ...
       ~strcmp(e, '(not supported)') && ~strcmp(e, '(empty response)')
        fprintf('  REFUSED %s\n     -> %s\n', cmd, e);
        ok = false;
        drainErrors(sc);
    end
end

function s = queryOpt(sc, cmd)
% A query the instrument is allowed to refuse. An unsupported mnemonic on the Infiniium
% family produces no response at all -- the command lands in the error queue and the
% read blocks -- so the timeout is dropped to 2 s for the attempt and the connection is
% cleaned up afterwards. Without the flush, the next legitimate read would collect this
% query's late-arriving answer and every reading after it would be off by one.
    old = sc.Timeout;
    s = '(not supported)';
    try
        sc.Timeout = 2;
        writeline(sc, cmd);
        s = strtrim(char(readline(sc)));
        if isempty(s); s = '(empty response)'; end
    catch
        try; flush(sc); catch; end
        drainErrors(sc);
    end
    try; sc.Timeout = old; catch; end
end

function drainErrors(sc)
% Empty the error queue quietly. A refused command leaves an entry behind, and a stale
% entry read later would be blamed on whatever command happened to be running then.
    old = sc.Timeout;
    try
        sc.Timeout = 2;
        for k = 1:20
            writeline(sc, ':SYSTem:ERRor? STRing');
            r = strtrim(char(readline(sc)));
            if isempty(r) || startsWith(r, '0') || startsWith(r, '+0'); break; end
        end
    catch
    end
    try; sc.Timeout = old; catch; end
end

function sc = openScope(resource)
% Keysight EXR404A / MXR, or the older DSO9000 Infiniium. All three run the Infiniium
% software so the SCPI below is unchanged; only the model match had to widen. USB is
% preferred when the same instrument also appears over LAN.
    if isempty(resource)
        L = visadevlist;
        M = upper(string(L.Model));
        hit = startsWith(M, "EXR") | startsWith(M, "MXR") | contains(M, "DSO");
        if ~any(hit)
            error(['No Keysight scope found in visadevlist (looked for EXR*, MXR*, ' ...
                   '*DSO*). Pass opts.scopeResource explicitly.']);
        end
        idx = find(hit);
        isUsb = startsWith(upper(string(L.ResourceName(idx))), "USB");
        if any(isUsb); idx = idx(find(isUsb, 1)); else; idx = idx(1); end
        resource = L.ResourceName(idx);
        fprintf('scope resource: %s (%s)\n', resource, L.Model(idx));
    end
    sc = visadev(resource);
    sc.Timeout = 20;
    fprintf('scope: %s\n', strtrim(query1(sc, '*IDN?')));
end

function ok = setInputFifty(sc, c)
% 50 ohm termination, VERIFIED by readback. The DSO90404A was 50 ohm only, so this
% never had to be set. The EXR inputs are switchable, and the 1 Mohm path is bandwidth
% limited to a few hundred MHz -- a 2.87 GHz tone into a 1 Mohm input reads as almost
% nothing, which is indistinguishable from a dead signal path. Two command spellings
% are tried because the Infiniium and InfiniiVision families disagree on the mnemonic.
    ok = false;
    attempts = { sprintf(':CHANnel%d:INPut DC50', c),      sprintf(':CHANnel%d:INPut?', c); ...
                 sprintf(':CHANnel%d:IMPedance FIFTy', c), sprintf(':CHANnel%d:IMPedance?', c) };
    for k = 1:size(attempts, 1)
        % queryOpt, not query1: on an instrument that rejects the first spelling the
        % readback is a 20 s stall per attempt otherwise. The EXR takes :INPut DC50 on
        % the first try and never reaches the second.
        try; writeline(sc, attempts{k,1}); catch; end
        r = upper(queryOpt(sc, attempts{k,2}));
        if contains(r, "DC50") || contains(r, "DCF") || contains(r, "FIFT") || contains(r, "50")
            fprintf('  scope ch%d input: %s (50 ohm confirmed)\n', c, r);
            ok = true; return
        end
    end
    warning(['scopeProbe:inputImpedance', newline, ...
             'Could not confirm 50 ohm on scope channel %d. If this input is on 1 Mohm ' ...
             'the tone is being rolled off by the input path and every amplitude below ' ...
             'is meaningless. Set it on the front panel before trusting anything.'], c);
end

function awg = findAwg()
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
    error('No AWG found via FrequencyGenerator.getFG().');
end

function o = withDefaults(o, d)
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(o, f{i}) || isempty(o.(f{i})); o.(f{i}) = d.(f{i}); end
    end
end
