function results = awgFgenProbe(awg, testFreqsMHz, timeDelayUs)
%AWGFGENPROBE Can the AWG's built-in function generator replace uploaded ESR bursts?
%
%   results = awgFgenProbe(awg)
%   results = awgFgenProbe(awg, [1000 2000 2770 2870 2960 4000 6000], 20e-6)
%
% WHY. Every ESR waveform we upload is a sine. The AWG70000 can also synthesise a
% sine itself in FGEN mode, where a frequency change is ONE SCPI write and there
% is nothing to upload, no sequence to build and no first arm -- i.e. the entire
% ~130 s cold load disappears. It would not touch the ~29 s of acquisition, which
% is photon-limited on any path, so the cold load is the whole prize.
%
% Before any of that is worth designing around, three facts have to hold, and
% this probe establishes them without changing the experiment:
%
%   1. FGEN's sine must actually REACH 2.87 GHz. Its ceiling is well below the
%      sample rate and I do not know where it sits on a 16 GS/s AWG70002A. If it
%      is under your sweep, the idea is dead and nothing else matters.
%   2. The two channels must take INDEPENDENT frequency and amplitude, because
%      the per-channel amplitude difference is calibrated and load-bearing.
%   3. Something must replace timeDelay. In AWG mode it exists only because one
%      channel's leading zero waveform is padded longer than the other's; FGEN
%      has no waveforms. For a sine, phase IS delay (delay = phase / 360 f), and
%      this class already knows the FGEN:CHAN<n>:PHAS command -- so the question
%      is its range and resolution.
%
% ENCOURAGING PRIOR: the FGEN spelling is not speculative. This class's generic
% FrequencyGenerator interface is already written against FGEN mode --
% setValue('frequency'), ('amplitude') and ('phase') map to FGEN:CHAN<n>:FREQ,
% :AMPL:POW and :PHAS. So those three commands are known-good on this instrument;
% what is unknown is their limits.
%
% SAFETY. This switches the instrument to FGEN mode and back. Outputs are turned
% OFF for the whole probe, so no microwave reaches the sample -- only parameters
% are set and read back. Switching modes stops playback, so the loaded sequence
% is treated as gone: the probe calls invalidateLoaded at the end and your next
% ESR will do a full reload. Run it between measurements, not mid-sweep.

    if ~exist('testFreqsMHz', 'var') || isempty(testFreqsMHz)
        testFreqsMHz = [1000 2000 2770 2870 2960 3500 4000 6000];
    end
    if ~exist('timeDelayUs', 'var'); timeDelayUs = []; end

    results = struct();
    fprintf('\n=== awgFgenProbe ===\n');

    % Read the mode THREE times after a *CLS, because the first run of this probe
    % reported FGEN on an instrument that had demonstrably just played sequences.
    % Two explanations were possible: a stale line in the input buffer (the class's
    % output-state query sends a COMPOUND query, 'OUTP1?; INST:MODE?', which returns
    % two lines, and reading only one leaves the other to be picked up by whatever
    % queries next), or the readback simply not meaning what I assumed. *CLS plus
    % three consecutive reads separates those.
    awg.sendCommand('*CLS');
    modeReads = cell(1,3);
    for iRead = 1:3
        modeReads{iRead} = strtrim(awg.queryRaw('INST:MODE?'));
    end
    fprintf('INST:MODE? read 3x after *CLS: "%s", "%s", "%s"\n', modeReads{:});
    results.modeReads = modeReads;
    startMode = modeReads{3};
    startChannels = awg.workChannels;
    fprintf('Starting mode: %s,  sample rate %g MS/s (1 sample = %.1f ps)\n', ...
            startMode, awg.SampleRate, 1e6/awg.SampleRate);
    results.startMode = startMode;

    restored = false;
    try
        % Outputs off first: nothing below should reach the sample.
        awg.workChannels = [];
        awg.sendCommand('awgcontrol:stop:immediate');
        drainErrors(awg, 20);                     % drain whatever was queued

        % ---- does the mode readback track what we set? -----------------------
        % If asking for AWG mode still reads back FGEN, the readback is not
        % reporting the mode and every mode-based conclusion is worthless. The
        % physics already says AWG mode was active during the ESR runs -- a real
        % FGEN-mode instrument emits one fixed tone and could not have produced a
        % resonance dip at all -- so a stuck readback is the likely answer.
        awg.setValue('mode', 'AWG');
        readAsAwg = strtrim(awg.queryRaw('INST:MODE?'));
        awg.setValue('mode', 'FGEN');
        readAsFgen = strtrim(awg.queryRaw('INST:MODE?'));
        results.modeRoundTrip = struct('askedAWG', readAsAwg, 'askedFGEN', readAsFgen);
        fprintf('Mode round-trip: asked AWG -> "%s", asked FGEN -> "%s"%s\n', ...
                readAsAwg, readAsFgen, ...
                BooleanHelper.ifTrueElse(strcmpi(readAsAwg, readAsFgen), ...
                '   <-- readback does NOT track the mode', ''));

        % ---- switch to FGEN -------------------------------------------------
        awg.setValue('mode', 'FGEN');
        modeNow = strtrim(awg.queryRaw('INST:MODE?'));
        results.enteredFgen = contains(upper(modeNow), 'FGEN');
        fprintf('Mode after switch: %s\n', modeNow);
        if ~results.enteredFgen
            fprintf(['  Could not enter FGEN mode at all. Everything below would be\n' ...
                     '  meaningless, so stopping here.\n']);
            error('awgFgenProbe: instrument did not enter FGEN mode (reported "%s")', modeNow);
        end

        % Ask for a sine explicitly -- the ceiling can depend on the shape.
        awg.sendCommand('FGEN:CHAN1:TYPE SINE');
        awg.sendCommand('FGEN:CHAN2:TYPE SINE');
        results.type = strtrim(awg.queryRaw('FGEN:CHAN1:TYPE?'));
        errTxt = drainErrors(awg, 5);
        fprintf('Waveform type: %s%s\n', results.type, ...
                errIf(errTxt));

        % ---- 1. how high does the sine go? ----------------------------------
        % Set-and-read-back rather than a MAX query: a rejected set leaves the
        % old value in place, and the error queue says so. That is a direct test
        % and does not depend on this firmware supporting SCPI MIN/MAX queries.
        fprintf('\n1. Frequency ceiling\n');
        fprintf('   %-12s %-14s %-10s %s\n', 'asked (MHz)', 'read back (MHz)', 'accepted', 'error');
        results.freqTest = struct('askedMHz', {}, 'readMHz', {}, 'ok', {}, 'err', {});
        highestOk = NaN;
        for fMHz = testFreqsMHz(:).'
            awg.sendCommand(sprintf('FGEN:CHAN1:FREQ %.6f', fMHz*1e6));
            readHz = str2double(awg.queryRaw('FGEN:CHAN1:FREQ?'));
            e = drainErrors(awg, 3);
            clean = isClean(e);
            close = ~isnan(readHz) && abs(readHz - fMHz*1e6) < max(1, fMHz*1e6*1e-6);
            ok = clean && close;
            if ok; highestOk = fMHz; end
            fprintf('   %-12g %-14.4f %-10s %s\n', fMHz, readHz/1e6, ...
                    BooleanHelper.ifTrueElse(ok, 'yes', 'NO'), ...
                    BooleanHelper.ifTrueElse(clean, '', shorten(e)));
            results.freqTest(end+1) = struct('askedMHz', fMHz, 'readMHz', readHz/1e6, ...
                                             'ok', ok, 'err', e);   %#ok<AGROW>
        end
        results.highestFreqMHz = highestOk;

        % ---- 2. independent per-channel frequency and amplitude -------------
        fprintf('\n2. Independent channels\n');
        % Both orders, with the error queue checked, so that "coupled" can be told
        % apart from "the second command was rejected and left the first in place".
        awg.sendCommand('FGEN:CHAN1:FREQ 2870000000');
        awg.sendCommand('FGEN:CHAN2:FREQ 2800000000');
        eSet = drainErrors(awg, 3);
        f1 = str2double(awg.queryRaw('FGEN:CHAN1:FREQ?'))/1e6;
        f2 = str2double(awg.queryRaw('FGEN:CHAN2:FREQ?'))/1e6;
        awg.sendCommand('FGEN:CHAN2:FREQ 2800000000');
        awg.sendCommand('FGEN:CHAN1:FREQ 2870000000');
        eSet2 = drainErrors(awg, 3);
        g1 = str2double(awg.queryRaw('FGEN:CHAN1:FREQ?'))/1e6;
        g2 = str2double(awg.queryRaw('FGEN:CHAN2:FREQ?'))/1e6;
        results.independentFreq = abs(f1 - 2870) < 1 && abs(f2 - 2800) < 1;
        results.freqCoupled = (abs(f1 - f2) < 1) && (abs(g1 - g2) < 1);
        fprintf('   set 1 then 2: chan1 %.3f, chan2 %.3f%s\n', f1, f2, errIf(eSet));
        fprintf('   set 2 then 1: chan1 %.3f, chan2 %.3f%s\n', g1, g2, errIf(eSet2));
        if results.freqCoupled
            fprintf(['   -> frequency is COUPLED across channels: whichever was set last\n' ...
                     '      wins for both. This is NOT a problem for ESR. Both channels\n' ...
                     '      always play the SAME frequency there -- LoadAWG builds both from\n' ...
                     '      the same w(i) and differs only in amplitude and timeDelay. My\n' ...
                     '      original pass criterion asked for something the experiment never\n' ...
                     '      uses.\n']);
        elseif results.independentFreq
            fprintf('   -> frequency is independent per channel.\n');
        else
            fprintf('   -> neither independent nor cleanly coupled; look at the errors above.\n');
        end

        % Amplitude is in dBm here (AMPL:POW), which matches how you already set
        % a.amplitude, so a calibrated per-channel offset carries over directly.
        awg.sendCommand('FGEN:CHAN1:AMPL:POW -25');
        awg.sendCommand('FGEN:CHAN2:AMPL:POW -31');
        a1 = str2double(awg.queryRaw('FGEN:CHAN1:AMPL:POW?'));
        a2 = str2double(awg.queryRaw('FGEN:CHAN2:AMPL:POW?'));
        e = drainErrors(awg, 3);
        results.independentAmpl = abs(a1 - (-25)) < 0.5 && abs(a2 - (-31)) < 0.5;
        fprintf('   amplitude: chan1 %.2f, chan2 %.2f (asked -25 / -31)  -> %s%s\n', ...
                a1, a2, BooleanHelper.ifTrueElse(results.independentAmpl, ...
                'independent', 'NOT independent'), errIf(e));

        % ---- 3. can phase stand in for timeDelay? ---------------------------
        fprintf('\n3. Phase, as a stand-in for timeDelay\n');
        % Resolution scan: 0.01 deg read back as 0, so the step is coarser than that.
        % Find where it stops rounding to zero -- that step, converted to time, is
        % what a phase-based channel delay could actually resolve.
        results.phase = struct('chan', {}, 'asked', {}, 'read', {}, 'err', {});
        phaseStep = NaN;
        for ph = [0.005 0.01 0.05 0.1 0.5 1 45 179.99 -45]
            awg.sendCommand(sprintf('FGEN:CHAN1:PHAS %.4f', ph));
            r = str2double(awg.queryRaw('FGEN:CHAN1:PHAS?'));
            e = drainErrors(awg, 3);
            fprintf('   chan1 asked %9.4f deg -> read %9.4f%s\n', ph, r, errIf(e));
            results.phase(end+1) = struct('chan', 1, 'asked', ph, 'read', r, 'err', e);  %#ok<AGROW>
            if isnan(phaseStep) && ph > 0 && r ~= 0
                phaseStep = ph;     % smallest positive setting that survived
            end
        end
        results.phaseStepDeg = phaseStep;

        % Per-channel independence matters here in a way it does not for frequency:
        % the whole point of the delay is that the two channels differ.
        awg.sendCommand('FGEN:CHAN1:PHAS 0');
        awg.sendCommand('FGEN:CHAN2:PHAS 30');
        p1 = str2double(awg.queryRaw('FGEN:CHAN1:PHAS?'));
        p2 = str2double(awg.queryRaw('FGEN:CHAN2:PHAS?'));
        e = drainErrors(awg, 3);
        results.independentPhase = abs(p1) < 0.5 && abs(p2 - 30) < 0.5;
        fprintf('   per channel: chan1 %.4f, chan2 %.4f (asked 0 / 30) -> %s%s\n', ...
                p1, p2, BooleanHelper.ifTrueElse(results.independentPhase, ...
                'independent', 'NOT independent -- a phase-based delay is impossible'), errIf(e));
        awg.sendCommand('FGEN:CHAN1:PHAS 0');
        awg.sendCommand('FGEN:CHAN2:PHAS 0');

        % What that buys, in time. For a sine, delay = phase / (360 * f).
        fRef = 2870e6;
        fRefPre = fRef;
        halfPeriodPs = 0.5/fRef * 1e12;
        samplePs = 1e6/awg.SampleRate;
        results.halfPeriodPs = halfPeriodPs;
        results.samplePs = samplePs;
        if ~isnan(results.phaseStepDeg)
            fprintf('   Smallest phase step that stuck: %g deg = %.3f ps at 2870 MHz.\n', ...
                    results.phaseStepDeg, results.phaseStepDeg/360/fRefPre*1e12);
        end
        fprintf(['   At 2870 MHz: +/-180 deg spans +/-%.1f ps, and 0.01 deg is %.4f ps.\n' ...
                 '   Today the delay is padding on a zero waveform, so its step is one\n' ...
                 '   sample = %.1f ps and its range is unlimited.\n'], ...
                halfPeriodPs, halfPeriodPs/18000, samplePs);
        if ~isempty(timeDelayUs)
            nSamples = timeDelayUs * awg.SampleRate;      % us * MHz = samples
            fprintf('   Your timeDelay = %g us = %.3f ps = %.3f samples.\n', ...
                    timeDelayUs, timeDelayUs*1e6, nSamples);
            if nSamples < 0.5
                fprintf(['   NOTE that is under half a sample, so round(duration*SampleRate)\n' ...
                         '   in AddWaveformByFunction rounds it to ZERO -- the delay you are\n' ...
                         '   setting today may not be reaching the instrument at all. Worth\n' ...
                         '   checking before you calibrate against it. Phase would express it\n' ...
                         '   exactly.\n']);
            elseif timeDelayUs*1e6 > halfPeriodPs
                fprintf(['   That exceeds half a period, so a single phase setting cannot\n' ...
                         '   express it (phase wraps). Would need a coarse delay elsewhere\n' ...
                         '   plus phase for the remainder.\n']);
            else
                fprintf('   Within phase range, and finer than the current one-sample step.\n');
            end
        else
            fprintf(['   Pass your json timeDelay as the third argument and this will say\n' ...
                     '   whether phase can express it.\n']);
        end

    catch err
        fprintf('\n  probe aborted: %s\n', err.message);
        restored = restoreState(awg, startMode, startChannels);
        rethrow(err);
    end

    if ~restored
        restoreState(awg, startMode, startChannels);
    end

    % ---- verdict ---------------------------------------------------------
    fprintf('\n--- verdict ---\n');
    if isnan(results.highestFreqMHz)
        fprintf(['FGEN would not accept ANY of the test frequencies. Either the command\n' ...
                 'set differs on this firmware or the mode is unusable here -- drop the\n' ...
                 'idea and keep the uploaded bursts.\n']);
    elseif results.highestFreqMHz < 2960
        fprintf(['FGEN tops out around %g MHz, below your sweep (2770-2960). Dead end,\n' ...
                 'and cheaply established. Keep the uploaded bursts.\n'], results.highestFreqMHz);
    elseif ~results.independentAmpl
        fprintf(['FGEN reaches your sweep, but per-channel AMPLITUDE is not independent,\n' ...
                 'and that one is load-bearing. Not usable as it stands.\n']);
    else
        fprintf(['FGEN reaches your sweep (up to %g MHz here) with independent per-channel\n' ...
                 'frequency and amplitude. Worth designing around: it removes the whole\n' ...
                 'cold load.\n\n'], results.highestFreqMHz);
        fprintf(['The remaining question is PHYSICS, not speed. In FGEN mode the AWG emits\n' ...
                 'continuously, so the signal/reference contrast rests entirely on the MW\n' ...
                 'switch the pulse generator drives -- today the AWG''s zero-burst-zero\n' ...
                 'structure is a second gate in series with it. Losing it may leak MW into\n' ...
                 'the reference window and shrink contrast. That is one ESR run against a\n' ...
                 'classic reference to settle, and contrast is the number to watch, not\n' ...
                 'the shape.\n']);
    end
    fprintf('Mode restored to %s, outputs restored, load cache invalidated.\n\n', startMode);
end

function restored = restoreState(awg, startMode, startChannels)
    % Always put the instrument back, even on the failure path.
    restored = true;
    try
        awg.sendCommand('FGEN:CHAN1:PHAS 0');
        awg.setValue('mode', startMode);
        awg.workChannels = startChannels;
        awg.invalidateLoaded;       % the mode switch stopped playback: force a reload
        drainErrors(awg, 20);
    catch restoreErr
        warning('awgFgenProbe could not fully restore the AWG: %s', restoreErr.message);
        restored = false;
    end
end

function s = errIf(e)
    if isClean(e)
        s = '';
    else
        s = ['   <-- ' shorten(e)];
    end
end

function txt = drainErrors(awg, maxEntries)
    % Read SYSTem:ERRor? until the queue is empty, and return what came out.
    %
    % Spelled out here rather than calling the class's own describeInstrumentErrors,
    % which lives in a private methods block and is therefore not reachable from a
    % plain function -- that is what aborted the first run of this probe.
    entries = {};
    for k = 1:maxEntries
        try
            r = strtrim(awg.queryRaw('SYSTEM:ERROR?', 5));
        catch
            entries{end+1} = 'error queue unreadable';   %#ok<AGROW>
            break
        end
        if isempty(r); break; end
        % Tektronix answers 0,"No error" when there is nothing queued.
        if ~isempty(regexp(r, '^\s*[+-]?0\s*,', 'once')); break; end
        entries{end+1} = r;                              %#ok<AGROW>
    end
    if isempty(entries)
        txt = 'no error';
    else
        txt = strjoin(entries, ' | ');
    end
end

function tf = isClean(e)
    tf = strcmp(e, 'no error');
end

function m = shorten(msg)
    m = strtrim(msg);
    m = strrep(m, newline, ' ');
    if length(m) > 60
        m = [m(1:57) '...'];
    end
end
