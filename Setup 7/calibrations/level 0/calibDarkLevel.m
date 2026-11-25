function results = calibDarkLevel(integrationTime, nRepeats, neutraliseOffset, nCycles)
%CALIBDARKLEVEL Measure the detector's zero, decomposed into its parts.
%
%   results = calibDarkLevel()
%   results = calibDarkLevel(0.1, 20)     % 0.1 s per read, 20 reads per level
%   results = calibDarkLevel([], [], false)  % leave voltage_offset as configured
%
% LEVEL 0 CALIBRATION. Depends on nothing. Nothing depends on it except that every
% other measurement's zero point comes from here, which is why it goes first.
%
% This is the photodiode equivalent of dark counts. On a counting setup you would
% measure a rate; here the meaningful quantity is a voltage, and it has more than
% one contribution. Measuring only "laser off" lumps them together, so this
% separates them by taking three states:
%
%   A  photodiode UNPOWERED (setSPCMEnable false)
%      -> digitizer input offset alone. Pure electronics, no light path at all.
%
%   B  photodiode powered, laser gate CLOSED
%      -> A + photodiode dark current + any green leaking through the closed AOM.
%
%   C  photodiode powered, laser gate closed, room light blocked or not
%      -> whatever ambient light reaches the detector, if you run it twice with
%         the enclosure open and closed. Optional, controlled by the prompt.
%
% B minus A is the detector's own contribution. Anything unexpectedly large in
% that difference is either a warm photodiode or real light leakage, and those
% call for different fixes, which is the reason to separate them.
%
% WHAT IT DOES NOT MEASURE. Background fluorescence and scattered green with the
% laser ON are power-dependent, so they belong to the saturation curve
% (calibLaserSaturation), which fits them as a linear-in-power term. Do not treat
% the number from here as the background for a laser-on measurement.
%
% OUTPUT. A struct with, per digitizer input, the mean level and an uncertainty
% that accounts for drift (see pdReadRepeated).
%
% WHERE THIS NUMBER DOES *NOT* GO. An earlier version of this file suggested
% writing the measured offset into spcm.voltage_offset. That was wrong and the
% suggestion is gone. Following the value through the code:
%
%   PhotoDiodeDigitizerNiDaqControlled -> digi.voltageOffset
%   Digitizer.prepareAcquisition       -> channelInfo(i).DcOffset
%
% so voltage_offset is a DC offset the digitizer APPLIES AT ITS INPUT to centre a
% signal inside the input range. It is an acquisition setting, not a correction
% constant. Writing a measured dark level there would shift what gets digitized
% rather than correct any reading, and would then be double-counted by any
% analysis that also subtracts a dark level.
%
% The measured dark level belongs in a calibration record of its own, consumed by
% whatever subtracts it. This function therefore writes nothing anywhere.
%
% SELF-REFERENCE, AND WHY voltage_offset IS ZEROED WHILE MEASURING
%
% A calibration is only reproducible if measuring the quantity does not depend on
% the currently stored value of that quantity. voltage_offset breaks that: the
% digitizer applies it at the input, so a non-zero stored value shifts exactly the
% level being measured here. Measure, store, re-measure, and you would get a
% different answer the second time -- and a third value the time after.
%
% So the offset is set to ZERO for the duration of the measurement (restored
% afterwards), which makes the result absolute and independent of history. Pass
% neutraliseOffset = false to measure under the configured offset instead, which
% is what you want if you are checking that a *deliberate* input offset is doing
% what you think.
%
% The opposite rule applies to max_voltage, input_impedance and the integration
% time: those must NOT be neutralised, because the dark level is only valid for
% subtraction if it was measured under the same acquisition settings the
% experiment uses. They are recorded instead. The distinction is: neutralise what
% you are measuring, match what the experiment will use.

    if ~exist('integrationTime', 'var') || isempty(integrationTime); integrationTime = 0.1; end
    if ~exist('nRepeats', 'var') || isempty(nRepeats); nRepeats = 20; end
    if ~exist('neutraliseOffset', 'var') || isempty(neutraliseOffset); neutraliseOffset = true; end
    % Interleaving blocks: more cycles cancels slower drift, at the cost of toggling
    % the photodiode supply more often. 5 is a reasonable compromise.
    if ~exist('nCycles', 'var') || isempty(nCycles); nCycles = 5; end
    settleTime = 0.5;

    spcm  = getObjByName(Spcm.NAME);
    laser = getObjByName(LaserGate.GREEN_LASER_NAME);
    if isempty(spcm);  throwBaseObjException(Spcm.NAME); end
    if isempty(laser); throwBaseObjException(LaserGate.GREEN_LASER_NAME); end
    if ~spcm.hasPhotodiode
        error(['calibDarkLevel is written for the photodiode path. This setup reports ' ...
               'hasPhotodiode = false, so the quantity to measure is a count rate, not ' ...
               'a voltage, and the decomposition below does not apply.']);
    end

    fprintf('\n=== calibDarkLevel: %g s per read, %d reads per level ===\n', ...
            integrationTime, nRepeats);

    results = struct('integrationTime', integrationTime, 'nRepeats', nRepeats);

    % Save everything we are about to touch.
    savedLaserOn = laser.isOn;
    savedConfig  = spcm.currentChannelGUI;
    digi = getObjByName(Digitizer.NAME);
    savedOffset = [];
    if neutraliseOffset && ~isempty(digi)
        savedOffset = digi.voltageOffset;
        if any(savedOffset ~= 0)
            fprintf(['voltage_offset is %s and is applied at the digitizer input, so it\n' ...
                     'would shift the very level being measured. Zeroing it for this\n' ...
                     'measurement and restoring afterwards.\n'], mat2str(savedOffset));
        end
        digi.voltageOffset = zeros(size(savedOffset));
    end

    try
        % ---- A and B, INTERLEAVED -------------------------------------------
        % A: photodiode unpowered -- setSPCMEnable(false) drops the niDaq gate line
        %    that powers it, so this is the digitizer's own offset.
        % B: photodiode powered, laser gate still closed.
        %
        % Measured in ALTERNATING short blocks rather than one long A then one long
        % B, and B-A formed as the mean of PAIRED differences. Once the range change
        % brought the resolution to ~5 uV, sequential A-then-B started returning
        % B-A = -53 uV at 7.7 sigma -- negative, which no dark current can be -- and
        % the same figure came out +1 uV on the previous run. It was resolving
        % baseline wander over the ~3 s between the two blocks, not the detector.
        % Pairing cancels any drift slower than one cycle, which is the same reason
        % calibLaserSaturation interleaves its reference.
        laser.isOn = false;
        perBlock = max(2, ceil(nRepeats / nCycles));
        Araw = []; Braw = []; dPairs = [];
        for cyc = 1:nCycles
            spcm.setSPCMEnable(false);
            pause(settleTime);
            a = pdReadRepeated(spcm, integrationTime, perBlock, 3);

            spcm.setSPCMEnable(true);
            pause(settleTime);
            b = pdReadRepeated(spcm, integrationTime, perBlock, 3);

            Araw = [Araw; a.raw];                       %#ok<AGROW>
            Braw = [Braw; b.raw];                       %#ok<AGROW>
            dPairs = [dPairs; b.mean - a.mean];         %#ok<AGROW>
        end

        A = summarise(Araw);
        B = summarise(Braw);
        results.digitizerOffset  = A;
        results.darkWithDetector = B;
        results.nCycles = nCycles;
        report(sprintf('A  digitizer offset (photodiode unpowered), %d blocks', nCycles), A);
        report(sprintf('B  detector powered, laser gate closed, %d blocks', nCycles), B);

        % ---- B - A, from paired differences ---------------------------------
        d = mean(dPairs, 1);
        if nCycles > 1
            dErr = std(dPairs, 0, 1) ./ sqrt(nCycles);
        else
            dErr = sqrt(B.sterr.^2 + A.sterr.^2);
        end
        results.detectorContribution      = d;
        results.detectorContributionSterr = dErr;
        results.pairedDifferences         = dPairs;
        fprintf('\nB - A  detector dark + any green leakage (paired, %d cycles)\n', nCycles);
        for c = 1:numel(d)
            fprintf('   input_%d  %+9.6f V  +/- %.6f   (%s)\n', c, d(c), dErr(c), ...
                    significance(d(c), dErr(c)));
        end
        if nCycles > 1
            fprintf('   per-cycle spread: %s uV\n', ...
                    mat2str(round(1e6*std(dPairs, 0, 1), 1)));
        end

    catch err
        restore(spcm, laser, savedLaserOn, savedConfig, digi, savedOffset);
        fprintf('\n  calibDarkLevel aborted: %s\n', err.message);
        rethrow(err);
    end

    restore(spcm, laser, savedLaserOn, savedConfig, digi, savedOffset);

    % ---- the acquisition settings this was measured UNDER -----------------
    % A calibration measures the system as configured, so the configuration is
    % part of the result. Recording it here is what makes the number reusable.
    fprintf('\n--- measured under these acquisition settings ---\n');
    fprintf('   voltage_offset was %s while measuring\n', ...
            BooleanHelper.ifTrueElse(neutraliseOffset, 'ZEROED (result is absolute)', ...
                                     'left as configured'));
    try
        j = JsonInfoReader.getJson().spcm;
        results.config = j;
        fprintf('   max_voltage      %s   -> digitizer input range +/- that\n', mat2str(j.max_voltage));
        fprintf('   voltage_offset   %s   (applied AT THE INPUT, not subtracted)\n', mat2str(j.voltage_offset));
        fprintf('   input_impedance  %s\n', mat2str(j.input_impedance));
        fprintf('   channel_delay    %s   (see the warning: zeroed at this sample rate)\n', mat2str(j.channel_delay));

        % Report the resolution actually achieved, and DO NOT recommend a range from
        % the dark level. An earlier version compared max_voltage against the dark
        % offset and concluded the range was 100x too big -- which was right at
        % +/-2 V by luck, and then kept saying it at +/-0.2 V where the range is
        % correctly matched to the ~20 mV SIGNAL. The dark level is not the signal;
        % judging headroom by it is comparing against the wrong number.
        % Two different numbers, and they were being conflated. The pooled scatter
        % across interleaved blocks includes drift BETWEEN blocks, so quoting it as
        % "resolution" overstated the noise once interleaving was added (5 uV became
        % 10.9 uV without anything getting worse). The paired spread is the honest
        % floor on a DIFFERENCE; the pooled scatter is the floor on an ABSOLUTE level
        % over the whole measurement.
        fprintf('\n   noise floor on an absolute level : %.1f uV  (includes drift between blocks)\n', ...
                1e6*max(A.sterr));
        if isfield(results, 'pairedDifferences') && size(results.pairedDifferences,1) > 1
            fprintf('   noise floor on a paired difference: %.1f uV  (drift cancelled)\n', ...
                    1e6*max(std(results.pairedDifferences, 0, 1)));
        end
        fprintf(['   Headroom is set by your largest SIGNAL, not by this dark level. On\n' ...
                 '   Setup 7 the fluorescence runs ~20 mV, so +/-0.2 V leaves 10x headroom\n' ...
                 '   and spends the resolution where it is wanted.\n']);
    catch
        fprintf('   could not read the spcm section of the setup json\n');
    end

    % ---- drift warnings --------------------------------------------------
    if any(A.drifted) || any(B.drifted)
        fprintf(['\n   NOTE the repeat-to-repeat scatter is more than 1%% of the level on at\n' ...
                 '   least one input, so it moved while measuring. On a dark level that\n' ...
                 '   usually means the photodiode is still warming up -- leave it powered\n' ...
                 '   for a few minutes and repeat.\n']);
    end

    % ---- the history file, and what to do with all of it ------------------
    results.history = appendHistory(results);
    verdict(results);
    fprintf('\n');
end

% ======================================================================

function h = appendHistory(results)
% Append this run to darkLevelHistory.mat in the current results folder.
%
% This exists for one reason: the ABSOLUTE dark level is not reproducible and the only way
% to see that is to compare runs. A single run reports +/- 8 uV and looks like a constant
% measured to five digits. It is not one, and a five-digit number that moves is worse than
% no number at all, because it invites being written down and reused.
%
% Soft on every failure. A calibration that cannot write its log has still measured the
% thing, and refusing to return a good measurement over a bookkeeping problem is the wrong
% trade.
    h = struct([]);
    try
        P = setup7Paths();
        f = fullfile(P.folder, 'darkLevelHistory.mat');
        rec = struct('when', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...      %#ok<TNOW1,DATST>
                     'digitizerOffset',  results.digitizerOffset.mean, ...
                     'darkWithDetector', results.darkWithDetector.mean, ...
                     'sterr',            results.darkWithDetector.sterr, ...
                     'contribution',     results.detectorContribution, ...
                     'contributionSterr',results.detectorContributionSterr);
        if isfile(f)
            S = load(f);
            if isfield(S, 'history') && isstruct(S.history) && ~isempty(S.history)
                h = S.history(:)';
            end
        end
        if isempty(h); h = rec; else; h(end+1) = rec; end
        history = h;                     %#ok<NASGU>  save() names the VARIABLE, not the field
        save(f, 'history', '-v7');
    catch err
        fprintf('\n   (could not update darkLevelHistory.mat: %s)\n', err.message);
    end
end

function verdict(r)
% Say plainly what to do, because "here are six voltages" is not an instruction.
    nCh = numel(r.darkWithDetector.mean);
    fprintf('\n============================================================\n');
    fprintf('WHAT TO DO WITH THIS\n');
    fprintf('============================================================\n');

    % ---- 1. the thing people expect and should not do ----
    fprintf(['\n1. NOTHING GOES IN expProperties.json. There is no dark-level field and\n' ...
             '   nothing downstream reads one. Do not add it. The dark level is consumed in\n' ...
             '   the SAME SESSION by calibLaserSaturation -- pass this struct straight in:\n\n' ...
             '       d = calibDarkLevel();\n' ...
             '       s = calibLaserSaturation([], [], [], d);\n\n' ...
             '   which is what level 0 does. Everything above level 0 is RATIOMETRIC -- ESR\n' ...
             '   and Rabi divide signal by reference from the same sequence -- so the dark\n' ...
             '   level cancels there and never needs to be known at all.\n']);

    % ---- 2. the pass/fail gate ----
    fprintf('\n2. THE GATE. This run either passes or it does not:\n\n');
    bad = false;
    for c = 1:nCh
        v = r.detectorContribution(c); e = r.detectorContributionSterr(c);
        if e > 0 && abs(v) > 3*e
            fprintf('   input_%d  FAIL  B - A is %.1f sigma from zero (%+.1f uV)\n', ...
                    c, abs(v)/e, 1e6*v);
            bad = true;
        else
            fprintf('   input_%d  pass  B - A consistent with zero\n', c);
        end
    end
    if bad
        fprintf(['\n   A non-zero B - A means the detector is contributing dark current or\n' ...
                 '   green is reaching it with the gate closed. POSITIVE and on input_1 is\n' ...
                 '   almost always a light leak: check the gate really closes, and cover the\n' ...
                 '   detector. Do not proceed to 0.2 -- a leak adds a power-independent\n' ...
                 '   offset to the fluorescence and inflates the fitted Fsat.\n']);
    else
        fprintf(['\n   Both pass, so the photodiode contributes nothing measurable and nothing\n' ...
                 '   is leaking past the gate. Go on to 0.2.\n']);
    end

    % ---- 3. the number that actually matters downstream ----
    % Over BOTH levels and both inputs, so this is the same figure the acquisition-settings
    % block above quotes rather than a second, slightly different one.
    nf = 1e6*max([r.digitizerOffset.sterr(:); r.darkWithDetector.sterr(:)]);
    fprintf(['\n3. THE NOISE FLOOR, %.1f uV, is the one number here worth remembering. It\n' ...
             '   bounds every SNR above this level: against Setup 7''s ~20 mV fluorescence it\n' ...
             '   is %.3f %%, so a 1 %% Rabi contrast is ~%.0fx the floor per read. If a later\n' ...
             '   fit is noise-limited, compare its scatter against this before blaming the\n' ...
             '   fit -- if they match, more averages is the only fix.\n'], ...
            nf, 100*nf/20e3, 0.01*20e3/nf);

    % ---- 4. run-to-run drift, the reason none of this is stored ----
    if isfield(r, 'history') && numel(r.history) > 1
        h = r.history;
        prev = h(end-1); this = h(end);
        d1 = 1e6*(this.digitizerOffset - prev.digitizerOffset);
        fprintf('\n4. AGAINST THE PREVIOUS RUN (%s):\n', prev.when);
        for c = 1:min(nCh, numel(d1))
            sig = 1e6*r.digitizerOffset.sterr(c);
            fprintf('   input_%d  digitizer zero moved %+.0f uV  (%.0fx the %.1f uV within-run error)\n', ...
                    c, d1(c), abs(d1(c))/max(sig,eps), sig);
        end
        if any(abs(d1) > 5*1e6*r.digitizerOffset.sterr(1:numel(d1)))
            fprintf(['\n   The ABSOLUTE zero moved by far more than the within-run error, which is\n' ...
                     '   the expected behaviour and the reason nothing stores it. A dark level\n' ...
                     '   measured an hour ago is wrong by roughly this much now. So:\n' ...
                     '     - re-measure it in the session you use it in, never carry it over;\n' ...
                     '     - trust PAIRED differences (B - A above) and ratios, which cancelled\n' ...
                     '       this drift completely -- B - A was consistent with zero both times\n' ...
                     '       while the absolute level shifted underneath it.\n']);
        end
    else
        fprintf(['\n4. No previous run on record, so there is nothing to compare against yet.\n' ...
                 '   This run has been appended to darkLevelHistory.mat in the results folder;\n' ...
                 '   the next one will report how far the absolute zero moved in between. That\n' ...
                 '   comparison is the point -- the absolute level drifts by far more than the\n' ...
                 '   quoted error, which is why it is never stored and always re-measured.\n']);
    end
    fprintf('\n============================================================\n');
end

function report(label, r)
    fprintf('\n%s\n', label);
    for c = 1:numel(r.mean)
        fprintf('   input_%d  %+9.6f V  +/- %.6f  (within %.6f, across %.6f)%s\n', ...
                c, r.mean(c), r.sterr(c), r.sterrWithin(c), r.sterrAcross(c), ...
                BooleanHelper.ifTrueElse(r.drifted(c), '  <-- drifting', ''));
    end
end

function s = significance(value, err)
    if err <= 0
        s = 'no error estimate';
    elseif abs(value) > 3*err
        s = sprintf('%.1f sigma from zero', abs(value)/err);
    else
        s = 'consistent with zero';
    end
end

function restore(spcm, laser, savedLaserOn, savedConfig, digi, savedOffset)
    try
        spcm.setSPCMEnable(true);
        spcm.currentChannelGUI = savedConfig;
        laser.isOn = savedLaserOn;
        if ~isempty(savedOffset) && ~isempty(digi)
            digi.voltageOffset = savedOffset;
        end
    catch err
        warning('calibDarkLevel could not fully restore the setup: %s', err.message);
    end
end

function out = summarise(raw)
    % Same summary pdReadRepeated produces, over reads pooled from several blocks.
    % The within-read error is not available here (it belongs to individual blocks),
    % so the quoted error is the scatter of the pooled reads.
    out.raw         = raw;
    out.nRepeats    = size(raw, 1);
    out.mean        = mean(raw, 1, 'omitnan');
    out.sterrWithin = zeros(1, size(raw, 2));
    out.sterrAcross = std(raw, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(raw), 1));
    out.sterr       = out.sterrAcross;
    out.hasWithin   = false(1, size(raw, 2));
    out.relScatter  = out.sterrAcross ./ max(abs(out.mean), eps);
    out.drifted     = out.relScatter > 0.01;
end
