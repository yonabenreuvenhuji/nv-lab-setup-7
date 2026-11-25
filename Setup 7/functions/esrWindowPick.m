function out = esrWindowPick(opts)
%ESRWINDOWPICK The CW-ESR detection window, derived and checked instead of scanned.
%
%   out = esrWindowPick                              % everything from expProperties.json
%   out = esrWindowPick(struct('tauPol', 7.97))      % override the json's tauPol
%   out = esrWindowPick(struct('detectionDuration', 100))   % check a value you chose
%
% NO HARDWARE. Arithmetic on the json and the AWG's own waveform rules. Replaces
% calibCwWindow, which scanned eight windows looking for a maximum that cannot exist.
%
% ===================================================================
% WHY THIS IS ARITHMETIC AND NOT A MEASUREMENT
% ===================================================================
%
% calibCwWindow scaled repeats as 1/window so that every point cost the same sequence
% time. That makes window x repeats -- the total integration time, the only thing shot
% noise cares about -- constant by construction, so shot-noise SNR per unit time does not
% depend on the window at all. There is no interior optimum to find, and a scan that
% appears to find one is finding an artefact. On the run that prompted this, it recommended
% 1600 us on the strength of two noise estimates, each a standard deviation from four
% off-resonance frequencies, that fell by a factor of five for no reason anyone could name.
%
% Only three things actually constrain the window, and all three are computable:
%
%   1. THE AWG BURST MUST FIT INSIDE ONE REPEAT.   <-- a hard limit, and the one that bites
%   2. The window must clear the start-of-window transient.
%   3. Drift eventually breaks the signal/reference cancellation. That ceiling sits far
%      above (1) and (2) -- the FoM was flat from 25 to 400 us with no turnover -- so it
%      does not constrain the choice, and confirming that is a two-point check rather than
%      an eight-point scan.
%
% ===================================================================
% (1) THE HARD LIMIT, AND WHY A SHORT WINDOW BREAKS
% ===================================================================
%
% This is the one that produced a zigzag spectrum -- alternate frequencies at ratio 1 with
% no dip at all -- at detectionDuration = 10 us while 100 us was clean.
%
% In CW ESR the AWG does not play one long tone. It plays a short waveform holding a whole
% number of MW cycles and LOOPS it:
%
%       segmentTime = nCycles / f,   nCycles = ceil(1.05 * MIN_WAVEFORM_LENGTH * f / S)
%       numSegments = round( (detectionDuration + laserInitializationDuration)
%                            / segmentTime )
%
% so the microwave plays for about  detectionDuration + laserInitializationDuration.
% Meanwhile the pulse sequence for one repeat is
%
%       2 * (detectionDuration + delayBetweenDetections)
%
% -- signal window, gap, reference window, gap. The AWG is triggered once per repeat and
% cannot accept a trigger while it is still playing. So:
%
%       detectionDuration + laserInitializationDuration
%                <=  2 * (detectionDuration + delayBetweenDetections)
%
%       i.e.   detectionDuration >= laserInitializationDuration
%                                   - 2 * delayBetweenDetections
%
% NOTE WHAT laserInitializationDuration IS DOING HERE. The CW ESR sequence has NO
% initialisation event -- the green laser is on continuously and the "initialisation" is
% just the reference window. The field enters ONLY through numSegments, where it is padding:
% how much extra microwave the AWG plays past the detection window so the burst comfortably
% covers it. AWGExpESR's own constructor sets it to 10 us, which satisfies the inequality
% for any window. expProperties.json holds a much larger value because the PULSED
% experiments need one, and loadJsonParams copies it in. One field, two unrelated jobs --
% the same collision detectionDuration itself has.
%
% With delayBetweenDetections = 10 us (the default when the detector is the NiDaq
% photodiode, i.e. Setup 7):
%
%   laserInit = 10.0 us  ->  any window is legal.        This is why it used to work.
%   laserInit = 26.5 us  ->  window >= 6.5 us.  At 10 us the burst is 36.5 us against a
%                            40 us repeat: legal, but 9% of headroom, and that is the run
%                            whose dip came back at half depth.
%   laserInit = 70.6 us  ->  window >= 50.5 us. At 10 us the burst is 80.5 us against a
%                            40 us repeat -- almost exactly TWO repeat periods, so the AWG
%                            is ready for every SECOND trigger and misses the ones between.
%                            Half the repeats get no microwave, and because the AWG is
%                            re-armed per frequency the accepted/missed pattern alternates
%                            with frequency index. That is the zigzag: every other point
%                            sits at ratio 1 because it genuinely had no drive.
%
% It fails silently. Nothing errors, no warning is printed, and the spectrum looks like a
% measurement.
%
% ===================================================================
% (2) THE TRANSIENT FLOOR
% ===================================================================
%
% Laser and microwave both start at the top of every window, so its first ~tau_pol carries
% less than the steady-state contrast and a short window dilutes the dip. The dilution is the
% AREA of that transient over the window, which is why tauPol is defined as the mean
% re-polarisation time (area under c/c0) -- see calibDetectionWindow. 10 x tauPol keeps it to
% about a tenth.
%
% tauPol is read from expProperties.json, where level 2 (cell 2.3) writes it. It is not an
% experiment property, so loadJsonParams ignores it; it is in the json so that this function
% and every script that calls it read one measured number instead of a hand-typed copy.
%
% OPTIONS
%   tauPol              us. Default [] = read "tauPol" from expProperties.json; if it is not
%                       there either, the transient floor is skipped (and said so).
%   detectionDuration   us. A value to CHECK rather than derive. Default [] = derive.
%   laserInit           us. Default [] = read expProperties.json.
%   delayBetweenDetections  us. Default [] = 10 if the detector is the NiDaq photodiode.
%   sampleRate          MHz. Default [] = ask the AWG, else 16000.
%   minWaveformLength   samples. Default [] = ask the AWG, else 2400.
%   jsonLocation        Default '' = setup7Paths.
%   frequency           MHz, for the segment arithmetic. Default 2870.
%
% OUTPUT
%   out.detectionDuration   what to use
%   out.minLegal            the hard limit from (1)
%   out.transientFloor      10*tauPol, or NaN
%   out.burstUs             the AWG burst length at the chosen window
%   out.repeatUs            one repeat at the chosen window
%   out.headroom            repeatUs / burstUs; must exceed 1
%   out.tauPol              the tauPol that was used (NaN if none)

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('tauPol', [], 'detectionDuration', [], ...
        'laserInit', [], 'delayBetweenDetections', [], 'sampleRate', [], ...
        'minWaveformLength', [], 'jsonLocation', '', 'frequency', 2870, ...
        'quiet', false), 'esrWindowPick');

    % ---- gather the inputs ----------------------------------------------
    [laserInit, srcInit] = resolveInit(opts);
    [opts.tauPol, srcTau] = resolveTauPol(opts);
    gap = opts.delayBetweenDetections;
    srcGap = 'given';
    if isempty(gap)
        gap = 10; srcGap = 'AWGExpESR default for the NiDaq photodiode detector';
        try
            spcm = getObjByName(Spcm.NAME);
            if ~isa(spcm, 'PhotoDiodeDigitizerNiDaqControlled')
                gap = 0; srcGap = 'AWGExpESR default for this detector class';
            end
        catch
        end
    end
    [S, minWf, srcAwg] = resolveAwg(opts);

    % quiet mode exists so this can be called as a library function -- calibFieldSurvey
    % resolves its ESR window through here, and does not want two screens of derivation
    % printed inside every field scan.
    say = @(varargin) qprintf(opts.quiet, varargin{:});

    say('\n=== esrWindowPick ===\n');
    say('   laserInitializationDuration  %8.3f us   (%s)\n', laserInit, srcInit);
    if ~isempty(opts.tauPol)
        say('   tauPol                       %8.3f us   (%s)\n', opts.tauPol, srcTau);
    end
    say('   delayBetweenDetections       %8.3f us   (%s)\n', gap, srcGap);
    say('   AWG %g MS/s, min waveform %d samples   (%s)\n', S, minWf, srcAwg);

    % ---- (1) the hard limit ---------------------------------------------
    minLegal = laserInit - 2*gap;
    say('\n1. HARD LIMIT: window >= laserInit - 2*gap = %.1f us.\n', minLegal);
    if minLegal <= 0
        say(['   Negative, so every window is legal. That is the healthy case and it is\n' ...
                 '   what a laserInitializationDuration near the class default of 10 us buys.\n']);
        minLegal = 0;
    else
        say(['   Below that the AWG burst is longer than one repeat and triggers are\n' ...
                 '   MISSED, silently. See the header: this is the zigzag spectrum.\n']);
    end

    % ---- (2) the transient floor ----------------------------------------
    floorT = NaN;
    if ~isempty(opts.tauPol) && isfinite(opts.tauPol)
        floorT = 10*opts.tauPol;
        say('\n2. TRANSIENT FLOOR: 10 x tau_pol = %.1f us.\n', floorT);
    else
        say(['\n2. TRANSIENT FLOOR: not computed -- no tauPol in the json. Run level 2 (2.2, then\n' ...
                 '   write it in 2.3), or pass struct(''tauPol'', ...).\n' ...
                 '   Without it the window is set by the hard limit alone, which is a safety\n' ...
                 '   condition rather than a good choice.\n']);
    end

    % ---- the pick --------------------------------------------------------
    if ~isempty(opts.detectionDuration)
        w = opts.detectionDuration;
        say('\n3. CHECKING THE VALUE YOU GAVE: %g us.\n', w);
    else
        w = nanMax([minLegal*1.5, floorT]);     % 50% margin on the hard limit
        if ~isfinite(w) || w <= 0; w = 100; end
        w = roundUpNice(w);
        say('\n3. PICKED %g us.\n', w);
        if isfinite(floorT) && floorT >= minLegal*1.5
            say('   Set by the transient floor. It clears the hard limit with room to spare.\n');
        else
            say(['   Set by the hard limit plus 50%% margin. Pass tauPol to have the\n' ...
                     '   transient floor considered too -- it is usually the larger of the two\n' ...
                     '   and it is the one with physics behind it.\n']);
        end
    end

    % ---- verify the arithmetic at the chosen window ----------------------
    [burstUs, segUs, nSeg] = awgBurst(w, laserInit, opts.frequency, S, minWf);
    repeatUs = 2*(w + gap);
    head = repeatUs / max(burstUs, eps);
    say('\n   AWG burst  %8.2f us   (%d loops of a %.4f us segment)\n', burstUs, nSeg, segUs);
    say('   one repeat %8.2f us   (2 x (%g + %g))\n', repeatUs, w, gap);
    say('   headroom   %8.2fx\n', head);
    if head < 1
        say(['\n   FAILS. The burst is %.2fx one repeat, so roughly %d in every %d triggers\n' ...
                 '   is missed and those repeats get no microwave at all. Expect a spectrum\n' ...
                 '   with points at ratio 1 interleaved with real ones -- the visible period\n' ...
                 '   depends on where each frequency''s re-arm lands, so it need not equal the\n' ...
                 '   ratio above.\n'], ...
                1/head, max(ceil(1/head)-1, 1), ceil(1/head));
        say('   FIX EITHER ONE:\n');
        say('       a.detectionDuration = %g;   %% raise the window\n', roundUpNice(minLegal*1.5));
        say('       a.laserInitializationDuration = %g;   %% or lower the padding\n', ...
                max(10, floor(w + 2*gap)));
    elseif head < 1.3
        say(['\n   LEGAL BUT TIGHT (%.0f%% headroom). Some triggers will land close to the\n' ...
                 '   end of a burst. Prefer more margin.\n'], 100*(head-1));
    else
        say('\n   Comfortable.\n');
    end

    out = struct('detectionDuration', w, 'minLegal', minLegal, 'transientFloor', floorT, ...
                 'burstUs', burstUs, 'repeatUs', repeatUs, 'headroom', head, ...
                 'laserInit', laserInit, 'gap', gap, 'tauPol', nanIfEmpty(opts.tauPol));

    if ~opts.quiet; verdict(out, opts); end
end

% ======================================================================

function verdict(o, opts)
    fprintf('\n============================================================\n');
    fprintf('WHAT TO DO WITH THIS\n');
    fprintf('============================================================\n');

    % If the window under test FAILED, do not print it back as the thing to set. An
    % earlier draft did exactly that -- "FAILS ... roughly 2 in 3 triggers missed" and then
    % "a.detectionDuration = 10" three lines later.
    w = o.detectionDuration;
    if o.headroom < 1
        w = roundUpNice(max(o.minLegal*1.5, nanMax([o.transientFloor, 0])));
        fprintf('\n1. THE WINDOW YOU CHECKED FAILS, so use %g us instead:\n\n', w);
    else
        fprintf('\n1. IN THE ESR CELL of the run script, before a.run:\n\n');
    end
    fprintf('       a.detectionDuration          = %g;\n', w);
    fprintf('       a.referenceDetectionDuration = %g;\n', w);
    if o.laserInit > w + 2*o.gap
        fprintf('       a.laserInitializationDuration = %g;   %% NOT the json value\n', ...
                max(10, floor(w + 2*o.gap)));
    end
    fprintf(['\n   None of these go in expProperties.json. detectionDuration there is the\n' ...
             '   PULSED window from calibDetectionWindow, a few microseconds, and\n' ...
             '   laserInitializationDuration there is a real initialisation for the pulsed\n' ...
             '   experiments. Both fields mean something different in CW ESR, and\n' ...
             '   loadJsonParams copies them in regardless -- which is exactly how a window\n' ...
             '   that used to work stopped working when the level-2 calibrations changed the\n' ...
             '   json.\n']);

    fprintf(['\n2. WHY THERE IS NO SCAN ANY MORE. calibCwWindow scanned eight windows for an\n' ...
             '   SNR-per-unit-time maximum. At fixed total integration time that maximum does\n' ...
             '   not exist: repeats scale as 1/window, so window x repeats is constant and\n' ...
             '   shot noise does not care how the time is divided. The scan cost minutes,\n' ...
             '   recommended 1600 us on the strength of two unexplained noise points, and its\n' ...
             '   shortest window was silently broken by the trigger problem above. Deleting it\n' ...
             '   loses nothing: the window is set by two computable floors and a ceiling that\n' ...
             '   is nowhere near either.\n']);

    fprintf(['\n3. IF YOU WANT THE CEILING CONFIRMED, it is two points, not eight. Run the ESR\n' ...
             '   at %g us and at %g us and compare the fitted dip depth over its error. If the\n' ...
             '   longer window is no better per unit time, drift is not limiting you and %g us\n' ...
             '   is right. If it IS better, the transient floor was underestimated -- raise\n' ...
             '   tau_pol, not the window.\n'], o.detectionDuration, 4*o.detectionDuration, ...
            o.detectionDuration);

    fprintf(['\n4. THE FAILURE TO REMEMBER. A CW-ESR window that is too short does not error,\n' ...
             '   does not warn, and produces a spectrum that looks real -- alternate points at\n' ...
             '   ratio 1 with the dip only on the others. If an ESR ever comes back looking\n' ...
             '   like a comb, run this function before anything else.\n']);
    fprintf('\n   Nothing was written.\n');
    fprintf('\n============================================================\n\n');
end

% ======================================================================

function [burstUs, segUs, nSeg] = awgBurst(det, laserInit, f, S, minWf)
% Exactly AWGExpESR's arithmetic, so the number here is the number the instrument plays.
    nCycles = max(ceil(1.05 * minWf * f / S), 1);
    segUs   = nCycles / f;
    nSeg    = round((det + laserInit) / segUs);
    burstUs = nSeg * segUs;
end

function [v, src] = resolveInit(opts)
    if ~isempty(opts.laserInit) && isfinite(opts.laserInit)
        v = opts.laserInit; src = 'given'; return
    end
    v = 10; src = 'AWGExpESR constructor default (json not read)';
    f = opts.jsonLocation;
    if isempty(f)
        try; f = setup7Paths().json; catch; f = ''; end
    end
    if ~isempty(f) && isfile(f)
        try
            j = jsondecode(fileread(f));
            if isfield(j, 'laserInitializationDuration')
                v = j.laserInitializationDuration;
                src = 'expProperties.json';
            end
        catch
        end
    end
end

function v = nanIfEmpty(v)
    if isempty(v); v = NaN; end
end

function [v, src] = resolveTauPol(opts)
% tauPol: given, else the json's "tauPol" (written by level 2), else empty.
    if ~isempty(opts.tauPol) && isfinite(opts.tauPol)
        v = opts.tauPol; src = 'given'; return
    end
    v = []; src = 'not in the json';
    f = opts.jsonLocation;
    if isempty(f)
        try; f = setup7Paths().json; catch; f = ''; end
    end
    if ~isempty(f) && exist(f, 'file')
        try
            j = jsondecode(fileread(f));
            if isfield(j, 'tauPol') && isnumeric(j.tauPol) && isscalar(j.tauPol) && j.tauPol > 0
                v = j.tauPol; src = 'expProperties.json';
            end
        catch
        end
    end
end

function [S, minWf, src] = resolveAwg(opts)
    S = opts.sampleRate; minWf = opts.minWaveformLength; src = 'given';
    if isempty(S) || isempty(minWf)
        src = 'defaults for FrequencyGeneratorTektronixAWG';
        if isempty(S);     S = 16000; end
        if isempty(minWf); minWf = 2400; end
        try
            awg = getObjByName('TektronixAWG');
            if ~isempty(awg)
                if ~isempty(awg.SampleRate); S = awg.SampleRate; end
                minWf = awg.MIN_WAVEFORM_LENGTH;
                src = 'read from the AWG object';
            end
        catch
        end
    end
end

function w = roundUpNice(w)
    steps = [10 15 20 25 30 40 50 60 75 100 125 150 200 250 300 400 500 600 800 1000 1600];
    j = find(steps >= w, 1);
    if isempty(j); w = ceil(w/100)*100; else; w = steps(j); end
end

function m = nanMax(v)
    v = v(isfinite(v));
    if isempty(v); m = NaN; else; m = max(v); end
end

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d); given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('%s: unknown option "%s". Valid: %s.', who, bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}); o.(known{i}) = d.(known{i}); end
    end
end

function qprintf(quiet, varargin)
    if ~quiet; fprintf(varargin{:}); end
end
