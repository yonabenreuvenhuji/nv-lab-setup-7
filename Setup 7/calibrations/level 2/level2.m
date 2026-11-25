%% ========================================================================
%  LEVEL 2  --  detection windows: when to count, for how long, and how long to re-polarise
%
%    2.1  calibLastDelay        fit: where the MW has actually stopped           -> lastDelay
%    2.1b calibDarkTime         diagnostic: the slow hardware rise after the step (writes nothing)
%    2.2  calibDetectionWindow  fit: the polarisation decay                       -> tauPol
%                               and, from it, window + init chosen together       -> detectionDuration,
%                                                                                    laserInitializationDuration
%    2.3  write                 the five numbers into expProperties.json (the ONLY cell that writes)
%    2.4  esrWindowPick         the CW-ESR window, arithmetic on the json. Not written there.
%
%  FIELD OFF for all of it. Offline check of both fits, no hardware: testLevel2Fits
%
%  ========================================================================
%  WHY lastDelay COMES FIRST
%  ========================================================================
%  calibLastDelay locates a STEP POSITION on the gap axis. A position does not move when the
%  contrast is halved, so it runs fine on a provisional window and a generous init. The decay in
%  2.2 is fragile in exactly the place lastDelay controls: if MW is still present when the window
%  opens, the first samples -- the ones the fast component comes from -- are wrong, silently.
%  So: the step first, on provisional numbers; then the decay with the MW known to be gone.
%
%  ========================================================================
%  WHAT CHANGED ON 2026-11-04, IN ONE PARAGRAPH EACH
%  ========================================================================
%  lastDelay is FITTED (erf step + background). The old "knee" rule gave 0.705 us because the
%    curve keeps creeping up after the step and the rule read the creep as "not yet flat". The
%    step itself is done by ~0.12-0.16 us. The creep is reported, not used; 2.1b tests what it is.
%  The polarisation decay is a SUM OF EXPONENTIALS (1-3, by BIC) plus a nuisance 60 kHz
%    oscillation that the MW burst triggers in the signal window. One exponential fitted to
%    log C (the old tau_pol = 11.3 us) was off by chi2_nu ~ 1300. tauPol is now the MEAN
%    re-polarisation time, which is what esrWindowPick's 10 x tauPol floor needs.
%  The window and the init are chosen TOGETHER. laserInitializationDuration already contains the
%    window in every pulsed class, so a longer window costs no sequence time; what it costs is
%    reference contamination from the slow component. The init is set by the carry-over into
%    the next point in walk mode (cap 5%, opts.carryOverMax). One map shows both.
%  tauPol lives in expProperties.json, so the run script and level 3 read it from there.
%  ========================================================================

%% 2.0  shared settings
h = getObjByName('Helmholtz');
h.control('off');

% setup7Paths picks the newest folder under calibrations\results\ and prints which one.
P = setup7Paths();
jsonLocation = P.json;

% Provisional values for 2.1 only, deliberately generous: 2.1 measures a position.
%   provWindow  a few us, inside the pulsed contrast where the step is easiest to see
%   provInit    100 us, so 2.1 does not need tauPol from 2.2
provWindow = 3;
provInit   = 100;


%% 2.1  last delay  --  the MW step, fitted on a fine scan
%
% The step is at ~0.1 us, so this scans only 0 to 0.5 us in 10 ns steps and fits only gaps up to
% 0.6 us. The scan is BALANCED (calibLastDelay's default since 2026-11-04): the reference window
% follows the same dark as the signal window. Unbalanced, the reference came straight after laser
% light, and the lower fluorescence after a dark period showed up as a slow "contrast" rise with
% the gap -- the long-scan puzzle of 2026-11-04. See 2.1b.
q = calibLastDelay(struct('detectionDuration', provWindow, ...
                          'laserInitializationDuration', provInit, ...
                          'jsonLocation', jsonLocation, ...
                          'gaps', 0.005:0.01:0.5, ...
                          'repeats', 1000, 'averages', 20));

% No usable step -> 2.2 runs on the stored value, and says so.
lastDelayUse = q.lastDelay;
if ~isfinite(lastDelayUse)
    j = jsondecode(fileread(jsonLocation));
    lastDelayUse = j.lastDelay;
    fprintf('\n2.1 gave no lastDelay, so 2.2 runs on the stored %g us. Fix 2.1 first.\n', lastDelayUse);
end


%% 2.1b  the slow rise in the long scans  (diagnostic, writes nothing)
%
% 2026-11-04: contrast rose ~0.6% over 0-150 us of gap, identically on and off resonance, with the
% MW off, the AWG never triggered and the switches mirrored. None of it was hardware: the scan ran
% the UNBALANCED sequence (Experiment's default), so the reference window followed laser light and
% the signal window followed dark, and fluorescence after dark is lower (0.14% after ~1 us, ~0.7%
% after 100 us). calibDarkTime's default now runs the old unbalanced scan next to the balanced one
% and the spin contrast balanced; expect the balanced ones to be flat. Not needed every time --
% once confirmed, 2.1 alone is enough. The other presets ('trigger', 'hardware') rerun the
% isolation tests if this ever comes back.
dt = calibDarkTime(struct('jsonLocation', jsonLocation, ...
                          'detectionDuration', provWindow, ...
                          'laserInitializationDuration', provInit));
outDir = fullfile(P.folder, 'level 2'); if ~exist(outDir, 'dir'); mkdir(outDir); end
save(fullfile(outDir, 'darkTime.mat'), 'dt');      % replot later: calibDarkTime(struct('data', dt))


%% 2.2  polarisation decay  --  tauPol, detectionDuration, laserInitializationDuration
%
% 50 us window: the slow component is ~37 us, and the carry-over that sets the init is read at
% twice the init, so a shorter window leaves the answer resting on extrapolation.
%
% The two knobs, if you need them:
%   carryOverMax  how much of a point may leak into the next one in WALK mode (default 0.05).
%                 The printed table shows the cost of 10/5/2/1%.
%   darkPerHalf   dark time per half-repeat in the experiments you care about (default 2 us:
%                 trigger + a short pulse or free evolution + lastDelay). For T1-like sequences
%                 with long dark times, re-run the analysis with it set larger -- no new data:
%                     r = calibDetectionWindow(struct('data', r, 'darkPerHalf', 500));
r = calibDetectionWindow(struct('window', 50, 'repeats', 100, 'averages', 100, ...
                                'lastDelay', lastDelayUse, ...
                                'jsonLocation', jsonLocation));


%% 2.3  write to expProperties.json
%
% The only cell in level 2 that changes the json. expPropsWrite edits the numbers in place,
% keeps a timestamped backup, and verifies the file afterwards. 'create' is needed once, the
% first time tauPol goes in.
%
% Look at both figures before running this. To see what would change without changing it:
%     expPropsWrite(toWrite, struct('jsonLocation', jsonLocation, 'create', true, 'dryRun', true));
toWrite = r.toWrite;
if isfield(q.toWrite, 'lastDelay'); toWrite.lastDelay = q.toWrite.lastDelay; end
res = expPropsWrite(toWrite, struct('jsonLocation', jsonLocation, 'create', true));

% Keep the raw results next to the figures, so either fit can be redone without hardware:
%     q = calibLastDelay(struct('data', q));   r = calibDetectionWindow(struct('data', r));
outDir = fullfile(P.folder, 'level 2');
if ~exist(outDir, 'dir'); mkdir(outDir); end
save(fullfile(outDir, 'level2_results.mat'), 'q', 'r');
fprintf('saved q and r to %s\n', fullfile(outDir, 'level2_results.mat'));

% Experiments built BEFORE this cell still hold the old values: loadJsonParams runs when an
% experiment is constructed. Rebuild them (b = AWGExpRabi; b.loadJsonParams(...)).


%% 2.4  CW-ESR window  --  computed from the json, not scanned, not written
%
% CW ESR needs a window far longer than the pulsed one, so it cannot live under the json's
% detectionDuration (that is the pulsed value, loaded into every experiment). It is computed
% where it is used, from two json numbers:
%   HARD LIMIT   det >= laserInitializationDuration - 2*delayBetweenDetections, or the AWG
%                misses triggers and the spectrum becomes a comb (no error, no warning)
%   TRANSIENT    det >= 10 * tauPol, which esrWindowPick now reads from the json
% The run script's ESR cell and level 3 call it the same way.
cw = esrWindowPick(struct('jsonLocation', jsonLocation));


%% 2.5  summary
fprintf('\n==================== LEVEL 2 RESULT ====================\n');
fprintf('IN expProperties.json  (loaded by Rabi, Ramsey, Echo, T1 ...)\n');
f = fieldnames(toWrite);
for i = 1:numel(f)
    fprintf('   %-28s %8.3f\n', f{i}, toWrite.(f{i}));
end
fprintf('   (tauPol is read by esrWindowPick, not by the experiments)\n');
fprintf('\nCW ESR, set in the ESR cells only  ->  never in the json\n');
fprintf('   detectionDuration           %8.0f us   (headroom %.2fx)\n', cw.detectionDuration, cw.headroom);
if cw.headroom < 1.2
    fprintf(['   WARNING: only %.0f%% headroom against the AWG burst -- override\n' ...
             '   laserInitializationDuration in the ESR cell, or the spectrum becomes a comb.\n'], ...
            100*(cw.headroom - 1));
end
fprintf('========================================================\n');
