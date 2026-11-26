%% ========================================================================
%  runScriptBasicExpAWG.m  --  Setup 7
%
%  Replaces runScriptBasicExpAWG.m, runScriptBasicExpAWG (1).m,
%  runScriptBasicExpAWG_ESR_updated.m and runScriptBasicExpAWGNew.m. Delete those four.
%
%  The dead experiments are GONE rather than commented out: the jump-mode probes, the FGEN
%  attempts, the mirror-fold plots, the loose percent/power arrays. A commented-out cell is
%  indistinguishable from one you have not got round to yet, and this file had a dozen.
%
%  ========================================================================
%  WHAT IS ON
%  ========================================================================
%
%  visadev transport         queries 1.2 s -> 2.3 ms, writes 2.72 -> 0.96 ms. Transport
%                            only; it cannot change what is measured.
%  commDelayTime = 0         drops a 5 ms busy-wait after every SCPI write.
%  SCPI batching             3084 writes -> 971 on a 96-point load. Checked offline for
%                            exact equivalence by testSequenceBatching.
%  load caching              an unchanged sweep skips the upload entirely, and verifies
%                            against the instrument with SLISt:SIZE? rather than trusting
%                            the hash. Keyed on the REQUESTED state, not the loaded state
%                            -- LoadAWG rewrites obj.tau, and keying on the result is why
%                            the Rabi cache never hit for a week.
%  short burst waveforms     the looped MW waveform is ~0.158 us / 2520 samples instead of
%                            ~1 us / 16000, about 4x less data per point. Automatic.
%  walk mode (ESR, Rabi,     one sequence holds every point and the AWG is armed ONCE.
%   Ramsey and Echo)         Echo is the cheapest of the four: its tau lives in the sequence
%                            repeat counts, so the waveform set never grows with the sweep.
%                            It also has the largest excess MW gate -- see the Echo cell.
%                            A re-arm costs ~0.6 s because assignSequence stops playback,
%                            so the whole speed question is how many re-arms a path needs:
%
%                                classic   one per point   (nPoints * averages)
%                                walk      one per pg.run  (ceil(repeats/P) * averages)
%
%                            Walk also removed the ~1.000 outlier points the per-point
%                            path shows, because those came from 96 separate arms.
%  passesPerRun derived      P = min(repeats, floor(maxDigitizerSegments/(2*(nTau+warmup))))
%                            rather than hand-picked. At 99 taus that gives 25 -- which is
%                            the number that used to be hardcoded, and it was right only
%                            for that one tau count.
%  channelPhase              per-channel phase written into the cosine. timeDelay cannot do
%                            this: it lengthens a padding waveform, which rounds to whole
%                            samples, and one sample at 16 GS/s is 64.6 degrees at
%                            2870 MHz. SET timeDelay = 0 -- the two are additive.
%  error bars                ratios formed per average, combined by inverse-variance
%                            weighting. averages > 1 no longer inflates them.
%  robust fits               esrFitRobust, rabiFitRobust and ramseyFitRobust all take ONE
%                            struct and return ONE struct:
%                                o = rabiFitRobust(struct('exp', b, 'nFreq', 2))
%                                o = ramseyFitRobust(struct('t', t, 'y', y, 'e', e))
%                            The data rides in the same struct as the options, so the vector
%                            and experiment forms are the same call. No positional outputs
%                            to mis-order -- the old struct-is-the-sixth-output convention
%                            cost three separate debugging sessions. The old positional
%                            values survive as named fields (o.pulses, o.f, o.T, o.f0,
%                            o.width, o.coeffError), so nothing that READ them changed.
%                            Every figure comes from fitStyle, so the three look identical;
%                            residual and spectrum panels are off unless asked for.
%
%  ========================================================================
%  WHAT IS OFF, AND STAYS OFF
%  ========================================================================
%
%  useSequenceJump = false   armed fine, returned frequency-scrambled data at every settle
%                            time tried (r = +0.04 throughout). Abandoned.
%  useFgen = false           the AWG's internal generator. It reaches the sweep,
%                            per-channel amplitude and phase work, the generator starts --
%                            and no ESR signal appears. Whether it emits usable power at
%                            2.87 GHz at all is a connector measurement, not a code
%                            question. Do not turn on.
%  keepOn = true             a stopgap. See the session cell.
%
%  ========================================================================
%  OPEN ITEMS
%  ========================================================================
%
%  1. RESOLVED. AWGExpRamsey and AWGExpT1 set obj.AWG.workChannels = obj.activeChannels in
%     prepare, and AWGExpEcho sets it in LoadAWG (line ~208) -- so all three are covered and
%     the note that used to say otherwise was wrong. keepOn = true is still on below and does
%     no harm; it can be turned off once someone wants to.
%  2. AWGExpESR uses workChannels = 1:obj.nChannels where activeChannels is correct. Only
%     matters at nChannels = 1, where it drives output 1 whatever MWChannel says.
%  3. Frequency quantization on the short burst: it holds a whole number of MW cycles but
%     the sample count is rounded, so each point's true frequency can be off by ~0.5 MHz,
%     varying across the sweep. A quarter of a 2 MHz step. Not yet fixed.
%  4. 'pulse streamer time out' on 1-4 of every 100 runs, from PulseStreamerNewClass. On
%     the walk path resyncEveryPass absorbs it -- one pass of data lost, mapping intact.
%     On the per-point path it is what produces the ~1.000 outliers.
%  5. AWGExpRabi and AWGExpRabiCont share NAME = 'AWG Rabi', so getExpByName is ambiguous
%     between them. Construct directly, as the cells below do.
%  ========================================================================


%% Start of session  --  run once
%
% No paths written down. setup7Paths derives the Setup 7 root from its own location in
% Setup 7\functions\ and picks the NEWEST dated folder under calibrations\results\,
% printing which one it chose. To pin a session to an older calibration on purpose:
%     P = setup7Paths('2026-08-14');
% and setup7Paths('list') shows what is available.
main;

P            = setup7Paths;
jsonLocation = P.json;
currentPath  = [P.folder filesep];

h = Helmholtz;
h.load(P.helmholtz);
h.control('auto');

% frequency, amplitudeDiff, channelPhase and timeDelay all live in expProperties.json now,
% so loadJsonParams supplies them and no cell below has to remember to set them. Note the
% field is called `frequency`, not `resonantFrequency` -- that is why the old
% resonant_frequency.mat and the `b.frequency = resonantFrequency` lines are gone.


%% AWG session settings  --  run once, right after main
%
% main rebuilds the device objects, so re-run this cell every time you re-run main. To
% revert any one thing, flip that one line.
%
% keepOn deserves its explanation. 'AWGCONTROL:STOP;OUTP1 OFF' was malformed for years: a
% bare ';' keeps the SCPI parser inside the AWGControl subtree, so the instrument rejected
% the second half with -113 and carried on silently -- disconnect never actually turned the
% MW outputs off. That is fixed. Now an experiment that re-runs WITHOUT reloading the AWG
% emits no MW and returns every point at ratio ~= 1. AWGExpESR and AWGExpRabi re-enable
% their outputs on every run; Echo, T1 and Ramsey do not have that line yet. keepOn = true
% leaves the outputs alone, which is the old behaviour, so they keep working.
%
% It is a stopgap that works by restoring a bug's side effect. Fix the three classes, then
% set this false.
awg = getObjByName('TektronixAWG');
if isempty(awg)
    fgs = FrequencyGenerator.getFG();
    for iFg = 1:numel(fgs)
        if contains(fgs{iFg}.name, 'TektronixAWG'); awg = fgs{iFg}; break; end
    end
end

awg.useVisadev           = true;
awg.commDelayTime        = 0;
awg.skipRedundantDeletes = true;
awg.armWithoutChecking   = false;   % leave off: with visadev the check costs ~2 ms
awg.keepOn               = true;    % see above
awg.resetCommandCount;

fprintf(['AWG session: visadev %d, commDelayTime %g s, skipRedundantDeletes %d, ' ...
         'keepOn %d\n'], awg.useVisadev, awg.commDelayTime, ...
        awg.skipRedundantDeletes, awg.keepOn);


%% ESR window  --  run this ONCE per session, before either ESR cell
%
% ===================================================================
% WHY THIS CELL EXISTS, AND WHY IT IS NOT OPTIONAL
% ===================================================================
%
% CW ESR needs a window two orders of magnitude longer than the pulsed one, and the number
% is not free to choose. In CW ESR the AWG loops a short microwave burst for
%
%       detectionDuration + laserInitializationDuration
%
% while one repeat of the sequence lasts
%
%       2 * (detectionDuration + delayBetweenDetections)
%
% and the AWG cannot accept the next trigger while it is still playing. So
%
%       detectionDuration  >=  laserInitializationDuration - 2*delayBetweenDetections
%
% Break that and repeats get NO MICROWAVE AT ALL, because their trigger arrived while the
% AWG was busy. The spectrum comes back as a comb -- alternate points sitting at exactly
% ratio 1, the dip only on the others. Nothing errors and nothing warns.
%
% The trap is that laserInitializationDuration is not an initialisation here. The CW
% sequence has no init event -- the green laser is on throughout -- and the field only sets
% how much extra microwave the AWG plays past the window. AWGExpESR's constructor uses 10 us
% for that, safe at any window; expProperties.json holds a much larger value because the
% PULSED experiments need a real initialisation, and loadJsonParams copies it in regardless.
% So editing the json for level 2 can silently break a CW window that worked last week. It
% did: detectionDuration = 10 was fine here until laserInitializationDuration went up, and
% then produced exactly that comb.
%
% esrWindowPick does the arithmetic: the larger of that hard limit and ~10*tauPol (the
% start-of-window transient), with margin. No hardware, no field, instant.
%
% tauPol comes from level 2.1. Set it by hand if you have not run level 2 this session --
% 13.6 us on this setup at the time of writing, but it is a per-sample number.

cw = esrWindowPick(struct('jsonLocation', jsonLocation));

% What the ESR cells below use. esrInit is the json value UNLESS the json value would break
% the inequality, in which case it is the largest value that does not.
esrWindow = cw.detectionDuration;
esrInit   = min(cw.laserInit, max(10, floor(esrWindow + 2*cw.gap)));
fprintf('\nESR cells: detectionDuration %g us, laserInitializationDuration %g us', ...
        esrWindow, esrInit);
if esrInit < cw.laserInit
    fprintf('  (json says %g -- overridden)', cw.laserInit);
end
fprintf('\n');


%% ESR, zero field  --  the single dip at D
%
% Field off, so all four orientations sit on top of each other near 2870 and there is one
% line. This is the cell for checking the sample is alive and the MW chain works.
%
% The first run at a given sweep uploads and prints a 'Cold load breakdown'. Every later
% run at the SAME frequency list, amplitude, detection durations and nChannels prints
% 'AWG already holds this waveform set' and skips the upload. Changing any of those
% invalidates it, and the printout names the field that changed.
%
% repeats is nearly free -- one pg.run per pass, no software in the loop. averages re-walks
% the whole sweep, so keep it small but above 1: at averages = 1 the reported sterr is not
% a standard error of anything.
h.control('off');

a = AWGExpESR;
a.loadJsonParams('all', jsonLocation);
a.frequency          = 2870 + (-100:1:100);
a.repeats            = 50;
a.averages           = 3;   
a.amplitude          = -16;
a.newOrderEveryAverage = true;
% From the window cell, NOT hardcoded. This line used to read 10, which was safe until
% level 2 raised laserInitializationDuration in the json and turned it into a comb.

a.detectionDuration  = esrWindow;
a.referenceDetectionDuration = a.detectionDuration;
a.laserInitializationDuration = esrInit;    % overrides the json -- see the window cell
a.isTracking         = 0;

a.useSequenceWalk = true;       % the fast, outlier-free path
a.useSequenceJump = false;      % never: scrambled frequencies
a.useFgen         = false;      % never yet: emits nothing usable
a.driveMWSwitches = true;       % normal operation; false is a switch diagnostic

tId = tic; a.run;
fprintf('ESR: %.1f s, %d SCPI commands\n', toc(tId), a.AWG.commandCount);
% esrFitRobust, not ESRfit. Same first five outputs so it is a drop-in, but the fitting is
% esrMultiLineFit -- N Lorentzians on a shared baseline, seeded from a robust envelope, with
% an explicit covariance -- and the extraction does not go through extract1Ddata, whose
% three successive cuts at 5, 3 and ONE standard deviation discard about a third of good
% Gaussian data and shrink every error bar downstream.
o = esrFitRobust(struct('exp', a, 'nLines', 1));

%% ESR, with field  --  the four-orientation spectrum
%
% Same sweep the level 3 calibrations use, so within a session they share the load cache --
% which means level 3 must use the SAME window, and it takes it from esrWindowPick too.
%
% This line used to read 50, described as "the CW value from calibCwWindow". Two things were
% wrong with that. calibCwWindow is gone -- it scanned for an optimum that does not exist at
% fixed integration time -- and 50 us is now ILLEGAL on this setup: with
% laserInitializationDuration at 70.6 and delayBetweenDetections at 10, the hard limit is
% 50.5 us, so 50 fails the inequality by half a microsecond and drops triggers.
h.control('auto');

a = AWGExpESR;
a.loadJsonParams('all', jsonLocation);
a.frequency          = 2750 : 2 : 3000;
a.repeats            = 100;
a.averages           = 3;
a.amplitude          = -35;
a.newOrderEveryAverage = true;

a.detectionDuration  = esrWindow;           % from the window cell, same as the zero-field one
a.referenceDetectionDuration = a.detectionDuration;
a.laserInitializationDuration = esrInit;
a.isTracking         = 0;

a.useSequenceWalk = true;
a.useSequenceJump = false;
a.useFgen         = false;
a.driveMWSwitches = true;

tId = tic; a.run;
fprintf('ESR: %.1f s, %d SCPI commands\n', toc(tId), a.AWG.commandCount);

% Two options here, and they answer different questions.
%
% Just the lines: esrFitRobust fits every line at once and gives centres with their
% covariance, which is what splittings need.
o = esrFitRobust(struct('exp', a, 'baselineOrder', 2));

% The lines AND the field vector: calibFieldSurvey solves [Bx By Bz D] from those centres
% via the exact Hamiltonian. Use it rather than eyeballing which dip is which -- at 15-20 G
% no single D pairs all four orientations, because the off-axis pair is pushed up 7-8 MHz by
% the second-order transverse shift while its splitting barely moves.
sv = calibFieldSurvey(struct('replay', a));
%% 
h.control('auto');

p = AWGExpESRPulsed;
p.loadJsonParams('all', jsonLocation);   % piTime, amplitude, windows, lastDelay, laser init
p.frequency          = 2817 + (-7:0.2:7);
p.repeats            = 50;
p.averages           = 3;   
p.amplitude          = -44;
p.mwPulseDuration = 0.556;
p.useSequenceWalk      = true;             % one sequence, one arm per pg.run
p.newOrderEveryAverage = true;
%p.plotPulsedSequence()
p.run;

%% 
out = esrPulsedVsCwCompare;

%% Rabi, COARSE  --  find the period
%
% Walk mode, and the standard tau list: 0.02:0.01:1.0, 99 points. 0.01 us steps rather than
% 0.005 for two measured reasons -- the residue waveform count is 0.2/step, so halving the
% resolution halves that part of the upload, and Nyquist at 0.01 us is 50 MHz against a
% Rabi of 13-20 MHz, so nothing is lost. Points per period affects how the plot looks, not
% how well the frequency is determined.
%
% passesPerRun = [] means the class derives it and prints what it chose, along with the
% walk-versus-classic re-arm count and a warning if walk is the slower path here.
h.control('auto');

b = AWGExpRabi;
b.loadJsonParams('all', jsonLocation);
b.nChannels = 2;                    % BEFORE amplitude: set.nChannels re-enters
b.amplitude = -20 + [0 0];          % set.amplitude, which would add amplitudeDiff twice
b.timeDelay = 0;                    % channelPhase is the phase knob; they are additive
b.repeats            = 300;
b.averages           = 2;           % never 1: sterr is not a standard error at 1
b.isTracking         = 0;
b.balancedSequence   = 1;
% detectionDuration is LEFT AS THE JSON SUPPLIES IT rather than hardcoded here. The
% current json says 3.88 us, which is neither the ~0.3 us that calibDetectionWindow
% reported nor the ~50 us CW value -- worth resolving, but this script should not silently
% override whichever number the calibration decided. To force one for a single run:
%     b.detectionDuration = 0.3; b.referenceDetectionDuration = b.detectionDuration;
b.shouldAutosave     = true;       % ~7 s per run writing a .mat nothing reads
b.tau                = round(0.02 : 0.06 : 6, 10);

b.useSequenceWalk    = true;
b.passesPerRun       = [];          % [] = derived; 25 at 99 taus
b.walkWarmupTriggers = 1;
b.resyncEveryPass    = true;

b.AWG.resetCommandCount; tId = tic;
b.run;
fprintf('coarse Rabi: %.1f s, %d SCPI commands\n', toc(tId), b.AWG.commandCount);

o = rabiFitRobust(struct('exp', b));        % one struct in, one struct out
fprintf('pi = %.1f +/- %.1f ns at %.3f +/- %.3f MHz\n', ...
        1e3*o.pulses(2), 1e3*o.pulsesErr(2), o.f, o.fErr);


%% Rabi, FINE  --  the measurement
%
% Tau range from the coarse fit: about 3.5 periods at >= 8 points per period. repeats is
% the cheap axis -- photons cost ~55 us per repeat, so 3000 repeats is ~16 s of a ~75 s run
% and buys sqrt(10) better SNR than 300. Contrast is what limits the fit, not tau count.
% o.pulses IS ALREADY IN MICROSECONDS. rabiFitRobust matches rabiFit's units, which is why
% the coarse cell above prints 1e3*o.pulses(2) to get nanoseconds. This line used to read
% t_pi = 1e-3*o.pulses(2), dividing by a thousand instead of leaving it alone, and the
% consequence was not a wrong tau range -- it was NO tau range:
%
%       t_pi   = 6.6e-5        instead of 0.066
%       tauMax = 4.6e-4        instead of 0.462
%       dtau   = 0.005         the floor, since round(1.3e-5, 3) is 0
%       tau    = 0.02 : 0.005 : 0.00046   ->  EMPTY, because start > stop
%
% and an empty vector reaches set.tau, which reports "Parameter must be a valid vector, and
% not too long!". The "too long" half is a red herring -- the same message covers both ends
% and the list was empty. Hence the explicit check below: a message that names the actual
% problem is worth four lines.
t_pi = o.pulses(2);                 % microseconds
if ~isfinite(t_pi) || t_pi <= 0
    error(['Fine Rabi: the coarse fit did not return a usable pi time (%g us). Re-run the ' ...
           'coarse cell and check its R^2 against the noise ceiling before continuing.'], t_pi);
end

tauMax = 7 * t_pi;                  % about 3.5 periods
dtau   = max(round((2*t_pi)/10, 3), 0.005);   % 10 points per period, floored at one AWG step

h.control('auto');
b = AWGExpRabi;
b.loadJsonParams('all', jsonLocation);
b.nChannels = 2;
b.amplitude = -16 + [0 0];
b.timeDelay = 0;
b.repeats            = 3000;
b.averages           = 2;
b.isTracking         = 0;
b.balancedSequence   = 1;
b.shouldAutosave     = false;
tauList = round(max(0.02, dtau) : dtau : tauMax, 3);
if numel(tauList) < 12
    error(['Fine Rabi: the tau list came out with %d points (%g to %g us in steps of %g), ' ...
           'which is not a sweep. t_pi = %g us from the coarse fit -- if that looks wrong, ' ...
           'the coarse fit is the problem, not this cell.'], ...
          numel(tauList), max(0.02, dtau), tauMax, dtau, t_pi);
end
fprintf('fine tau: %d points, %.3f to %.3f us, %.1f points per period, %.1f periods\n', ...
        numel(tauList), min(tauList), max(tauList), (2*t_pi)/dtau, ...
        (max(tauList) - min(tauList))/(2*t_pi));
b.tau                = tauList;

b.useSequenceWalk    = true;
b.passesPerRun       = [];
b.walkWarmupTriggers = 1;
b.resyncEveryPass    = true;

b.AWG.resetCommandCount; tId = tic;
b.run;
fprintf('fine Rabi: %.1f s, %d SCPI commands\n', toc(tId), b.AWG.commandCount);

o = rabiFitRobust(struct('exp', b));
fprintf('pi = %.1f +/- %.1f ns at %.3f +/- %.3f MHz\n', ...
        1e3*o.pulses(2), 1e3*o.pulsesErr(2), o.f, o.fErr);


%% Rabi, sanity check  --  on whichever Rabi you just ran
%
% The noise ceiling is 1 - sum(e^2)/sum((y-mean(y))^2): the highest R^2 any fit could reach
% given how much of the variance is signal. R^2 of 0.81 against a ceiling of 0.83 is as
% good as the data allow; the same 0.81 against 0.99 is a real failure. One number cannot
% distinguish those, which is why both are printed here and in calibB1Rabi's table.
fprintf('\ncontrast          %.2f %%\n', 100*o.contrast);
fprintf('R^2               %.4f   (noise ceiling %.4f)\n', ...
        o.r2, 1 - sum(o.e.^2)/sum((o.y - mean(o.y)).^2));
fprintf('points per period %.1f     (below 8 the plot looks worse, the fit does not)\n', ...
        o.pointsPerPeriod);
fprintf('periods spanned   %.1f     (>= 3 for a frequency, >= 6 for a decay)\n', ...
        o.periodsSpanned);
fprintf('resid autocorr    %+.2f     (|.| > 0.4 means one frequency is the wrong model)\n', ...
        o.residAutocorr);
if o.contrast < 0.01
    fprintf(['\nContrast under 1%%: check the MW is reaching the sample. If workChannels\n' ...
             'is missing from prepare AND keepOn is false, the outputs are OFF and every\n' ...
             'point returns ratio ~= 1 with no error at all.\n']);
end
if abs(o.residAutocorr) > 0.4
    fprintf(['\nTry rabiFitRobust(struct(''exp'', b, ''nFreq'', 2)): a second orientation may be\n' ...
             'driven off resonance, adding a component at sqrt(Omega^2 + Delta^2).\n']);
end


%% Ramsey  --  T2*, on resonance, walk mode
%
% ========================================================================
% WHAT CHANGED, AND WHY EACH ONE
% ========================================================================
%
% WALK IS ON. It was suspected of halving the Ramsey contrast; it was not. That suspicion
% came from a contrast statistic in ramseyWalkDiag whose own scatter was 32% and which
% responded more to grid-versus-phase alignment than to contrast -- see testContrastMetric,
% which reproduces the whole thing in two seconds. Measured properly since:
%
%       Rabi     98% +/- 1.7%    ESR     99% +/- 1.9%
%       Ramsey  101% +/- 5.7% at a 2.53 us gate, 106% +/- 5.3% at 4.05 us
%
% so walk reproduces classic on all three. ramseyGateTest re-runs that check in ~4 minutes.
%
% NO DETUNING. The old cell used `r.frequency = r.frequency - 10` to "detune, so there is a
% fringe to fit". A Ramsey on this sample does not need one: the 14-N nucleus splits the
% line into three, and the fringe at zero detuning is the 2.158 MHz beat between them. That
% is what you see in every on-resonance run here. Detuning by 10 MHz makes the fastest
% component 12.2 MHz instead of 2.16, so the tau step has to shrink by 5.6x for the SAME
% T2* -- five times the points, five times the run, no extra information.
%
% ON-RESONANCE TAU GRID. The step must resolve the FASTEST component, |detuning| + A, not
% the detuning. At zero detuning that is A = 2.158 MHz, so 13 points per 0.463 us beat is
% 0.036 us and 2.5 us needs about 70 points. That is the whole 5x saving.
%
% laserInitializationDuration = 16. About 81% of a Ramsey shot is laser initialisation --
% four events of (initialisation - detection) per repeat, 109 us of a 134 us shot at the
% json's 29.4. Level 2 put the SNR-per-unit-time optimum at 15-20 us. Delete the line to go
% back to whatever the json says.
%
% averages = 2 at minimum, never 1. At one average sterr is the within-average scatter with
% nothing to compare it against, and every error bar downstream is decoration.
%
% workChannels: the class DOES set it now (AWGExpRamsey.prepare, obj.AWG.workChannels =
% obj.activeChannels). The old comment here saying it does not is stale. keepOn = true is
% still on in the session cell and does no harm.
% ========================================================================
h.control('auto');

r = AWGExpRamsey;
r.loadJsonParams('all', jsonLocation);
r.nChannels = 2;                    % BEFORE amplitude: set.nChannels re-enters
r.amplitude = -16 + [0 0];          % set.amplitude, which would add amplitudeDiff twice
r.timeDelay = 0;                    % channelPhase is the phase knob; they are additive
r.halfPiTime = 0.0539/2;
r.repeats          = 3000;
r.averages         = 2;
r.balancedSequence = 1;
r.constantTime     = 1;             % or the laser duty cycle changes across the scan and
r.frequency = 2820.6;                                    % biases T2* SHORT
r.smallDelay       = 0.5;
r.isTracking       = 0;
r.shouldAutosave   = false;
r.laserInitializationDuration = 16; % see above; remove to use the json's value

hyperfineMHz = 2.158;               % 14-N. Use 3.030 and nHyperfine = 2 for a 15-N sample.
tauMax       = 2.5;                 % past T2*, or the decay is an extrapolation
perBeat      = 30;
dtau         = round(1/(hyperfineMHz*perBeat), 4);
r.tau        = round(0.005 : dtau : tauMax, 4);
fprintf('Ramsey tau: %d points, %.4f us step, %.2f us span = %.1f beats of %.3f us\n', ...
        numel(r.tau), dtau, max(r.tau), max(r.tau)*hyperfineMHz, 1/hyperfineMHz);

r.useSequenceWalk    = true;
r.passesPerRun       = [];          % [] = derived from the digitizer ceiling, and printed
r.walkWarmupTriggers = 1;
r.resyncEveryPass    = true;

r.AWG.resetCommandCount; tId = tic;
r.run;
fprintf('Ramsey: %.1f s, %d SCPI commands\n', toc(tId), r.AWG.commandCount);

% ONE STRUCT IN, ONE STRUCT OUT -- for all three robust fits now. No positional outputs to
% get in the wrong order, and the data rides in the same struct as the options, so
% ramseyFitRobust(struct('t', t, 'y', y, 'e', e)) fits vectors with the same call.
o = ramseyFitRobust(struct('exp', r));
fprintf('T2* = %.3f +/- %.3f us, detuning %.3f +/- %.3f MHz, contrast %.2f %%\n', ...
        o.T2star, o.T2starErr, o.delta, o.deltaErr, 100*o.contrast);


%% Ramsey, sanity check  --  on the fit you just ran
%
% Same idea as the Rabi check: one number cannot tell a bad fit from a noisy one, so the
% ceiling is printed next to R^2 and the diagnostics that decide whether the numbers are
% quotable are printed next to those.
fprintf('\nT2*               %.3f +/- %.3f us\n', o.T2star, o.T2starErr);
fprintf('detuning          %.3f +/- %.3f MHz\n', o.delta, o.deltaErr);
fprintf('contrast          %.3f +/- %.3f %%\n', 100*o.contrast, 100*o.contrastErr);
fprintf('stretch p         %.2f +/- %.2f    (2 = Gaussian bath, 1 = Lorentzian)\n', o.p, o.pErr);
fprintf('components        %s MHz\n', mat2str(round(o.freqs*1000)/1000));
fprintf('R^2               %.4f   (noise ceiling %.4f)\n', o.r2, o.r2Ceiling);
fprintf('reduced chi2      %.2f\n', o.chi2Reduced);
fprintf('resid acf         %+.2f    (|.| > 0.4 = the residual walks, model is missing\n', ...
        o.residAutocorr);
fprintf('                          something. Phase IS fitted, freely and per projection,\n');
fprintf('                          so a residual that LOOKS phase-shifted is a FREQUENCY\n');
fprintf('                          error -- try showResid = true, then fitHyperfine)\n');
if o.hyperfineFreed
    fprintf('hyperfine A       FREED by autoHyperfine: %.3f -> %.3f MHz\n', o.AFixed, o.A);
end
fprintf('span              %.2f us = %.1f beats, %.1f points per beat\n', ...
        o.spanUs, o.beatsSpanned, o.pointsPerBeat);
if o.traces == 2
    fprintf('|0>/|1> amplitude %+.3f     (should be -1; anything else and the two readout\n', ...
            o.ampRatio);
    fprintf('                            windows are not what the class documents)\n');
end
if o.T2star > 0.9*o.spanUs
    fprintf(['\nT2* is at or past the end of the tau list, so it is an extrapolation rather\n' ...
             'than a measurement. Raise tauMax to about %.1f us and re-run.\n'], 3*o.T2star);
end
if o.nDown > 0
    fprintf(['\n%d point(s) were down-weighted, at tau = %s us. On the CLASSIC path those are\n' ...
             'usually the per-point re-arm: repeats that fire before the AWG resumes see no\n' ...
             'microwave, so the point is diluted toward ratio = 1. Walk mode does not produce\n' ...
             'them, so on this cell they are worth a look rather than a shrug.\n'], ...
            o.nDown, mat2str(round(o.downTau*1000)/1000));
end
if o.contrast < 0.002
    fprintf(['\nContrast under 0.2%%: check halfPiTime in the json before blaming the sample.\n' ...
             'A wrong pi/2 costs contrast on every path at once and looks exactly like this.\n']);
end

% Useful variants, none of which need a different cell:
%   ramseyFitRobust(struct('exp', r, 'showResid', true))      % see the residual structure
%   ramseyFitRobust(struct('exp', r, 'showSpectrum', true))    % the three components
%   ramseyFitRobust(struct('exp', r, 'useProjection', 'one'))  % |1> alone, not the joint fit
%   ramseyFitRobust(struct('exp', r, 'nHyperfine', 2))         % 15-N sample
%   ramseyFitRobust(struct('exp', r, 'fitHyperfine', true))    % force A free
%   testRamseyFit                                              % grade the fitter, no hardware


%% Echo  --  T2, walk mode
%
% ========================================================================
% BEFORE THE FIRST WALKED ECHO, RUN THESE TWO
% ========================================================================
%
%   >> testEchoWalk       % ~1 s, no hardware. Block layout, trigger count, tau mapping.
%   >> echoGateTest       % ~4 min. Does the fixed gate leak? SEE BELOW -- this one matters.
%
% ECHO IS THE CHEAPEST OF THE THREE TO WALK, for a reason worth knowing: tau already lives in
% the AWG SEQUENCE, as the repeat count of a zero-amplitude step. The waveform set is fixed at
% one zero step plus five per channel and does not grow with the tau list at all, so walking
% Echo costs no extra upload whatsoever -- it only concatenates sequences that were being
% built anyway. Each tau's existing sequence is already a walk block: two trigger-gated steps
% among 2*(3*cycles+2), which is exactly one pulse-generator repeat.
%
% AND IT HAS THE WORST GATE PROBLEM OF THE THREE. Walk cannot rewrite the microwave gate per
% point, so it sits at the longest tau:
%
%       Ramsey, MEASURED:  4.05 us of excess gate cost nothing (106% +/- 5% of classic)
%       Echo, 1-50 us:     98 us of excess.  25x more.
%       Echo, XY8 to 200:  3.2 ms of excess. 800x more.
%
% Leakage through an open switch is a rate, so the dose scales with the excess, and the effect
% would be lost contrast at SHORT tau -- which on an Echo is indistinguishable from a shorter
% T2. A plausible number, not a crash. Nothing measured on this setup covers 98 us, which is
% why echoGateTest exists and why prepare prints a warning with the actual excess in it.
%
% If echoGateTest says the gate leaks: BAND THE TAU LIST. Walk two or three narrow ranges,
% each with its own gate, so no block sees much excess. One extra arm per band, most of the
% saving kept.
%
% maxRunSeconds is Echo's own knob and Ramsey has no equivalent. Echo's sequence duration is
% 2*max(tau)*cycles and the fixed gate makes every block pay it, so a long-tau XY8 run can
% otherwise disappear into one uninterruptible pg.run for minutes -- no plot updates, no
% emergency stop, and a digitizer timeout if the estimate is off. The class derives
% passesPerRun from BOTH ceilings and prints which one binds.
% ========================================================================
h.control('auto');

e = AWGExpEcho;
e.loadJsonParams('all', jsonLocation);
e.nChannels = 2;                    % BEFORE amplitude: set.nChannels re-enters
e.amplitude = -16 + [0 0];          % set.amplitude, which would add amplitudeDiff twice
e.timeDelay = 0;                    % channelPhase is the phase knob; they are additive

e.repeats          = 2000;      % 200 makes the whole decay only ~7 error bars tall; the
                                % revivals still show, but T2 gets a wide error bar
e.averages         = 2;             % never 1: sterr is not a standard error at 1
e.balancedSequence = 1;
e.constantTime     = 1;
e.smallDelay       = 1;
e.isTracking       = 0;
e.shouldAutosave   = false;
e.laserInitializationDuration = 16;
e.cycles           = 1;             % 1 = plain Hahn echo. See the XY8 cell below for 8.
e.pulsePhase       = 'Y';
e.frequency=2820.6;
e.halfPiTime = 0.0539/2;
e.piTime = 0.0539;


% Log spacing: the T2 knee is what is being measured, and a linear grid spends most of its
% points on the flat tail. round to 3 dp because LoadAWG snaps tau to the base waveform
% period anyway and then removes duplicates -- better to know what was asked for.
% Linear, not log: a log grid cannot sample a revival at the long end. See the note above on
% how far to run -- if the fit reports one revival only, double this.
e.tau = unique(round(linspace(1, 60, 60), 3));
fprintf('Echo tau: %d points, %.1f to %.1f us -> evolution %.1f to %.1f us at cycles = %d\n', ...
        numel(e.tau), min(e.tau), max(e.tau), ...
        2*min(e.tau)*e.cycles, 2*max(e.tau)*e.cycles, e.cycles);

e.useSequenceWalk    = true;
e.passesPerRun       = [];          % [] = derived from BOTH ceilings, and printed
e.walkWarmupTriggers = 1;
e.resyncEveryPass    = true;
e.maxRunSeconds      = 30;          % cap on ONE pg.run, so the experiment stays responsive

e.AWG.resetCommandCount; tId = tic;
e.run;
fprintf('Echo: %.1f s, %d SCPI commands\n', toc(tId), e.AWG.commandCount);

% echoFitRobust, same one-struct-in / one-struct-out convention as the other three. The old
% EchoFit is no longer called from here: it had a ten-argument positional signature and a
% fixed decay exponent.
%
%   >> testEchoFit        % ~1 min, no hardware. 22 checks. Run it after touching the fitter.
%
% WHAT IT FITS, AND THE TWO THINGS TO READ FIRST:
%
%   p, the stretch exponent, is FITTED rather than pinned at 1, 2 or 3, because which one is
%   right is a statement about the bath: 1 is a fast bath or a relaxation-limited echo, 2 is
%   slow Gaussian noise, 3 is the dipolar 13-C spin bath -- the textbook single-NV answer.
%   Because p moves, T2 alone is not comparable between runs, so compare out.coherenceArea
%   (= T2*gamma(1+1/p)) whenever p has changed.
%
%   THE 13-C REVIVALS are in the model, off unless the data pays for them. With a field on,
%   the echo collapses and revives at the 13-C Larmor period, and a monotone decay fitted
%   through that returns a SHORT T2 with a small error bar -- the same class of mistake as a
%   single cosine on an on-resonance Ramsey. The comb is accepted only on a 3-sigma joint test
%   for its two extra parameters, and the report says which test failed when it is refused.
%
% THE TAU GRID DECIDES WHETHER THE REVIVALS ARE FITTABLE, so it is not a free choice.
%
% Log spacing puts the points at the knee and leaves ~13 us between them at the long end,
% coarser than a revival; the fit says "under 6 points per revival" when that happens. A LINEAR
% grid is the right one whenever a field is on. It costs nothing on the knee at 60 points.
%
% HOW FAR TO GO is the part worth thinking about, and the report tells you after the first run.
% The comb is anchored at t = 0, so ONE revival is enough to fit the period -- but not enough to
% separate T2 from p, which then slide along a curve of equally good fits. If the report says
% "ONE revival only", double the sweep so a second lands in range and both become measurable.
% The report prints the tau to use. Until then, either quote T2 with pFixed and say which p you
% assumed, or quote out.coherenceArea, which moves less across that curve than T2 does.
o = echoFitRobust(struct('exp', e, 'T1us', []));      % T1us fills in the T2 <= 2*T1 check
fprintf('Echo: T2 = %.3f +/- %.3f us, p = %.2f, contrast = %.3f%%\n', ...
        o.T2, o.T2Err, o.p, 100*o.contrast);


%% Echo, XY8  --  T2 under decoupling. Same class, three settings.
%
% AWGExpEcho already IS an XY-N sequence generator; a Hahn echo is just the cycles = 1 case.
% Three lines turn it into XY8:
%
%       e.cycles     = 8;
%       e.pulsePhase = 'XYXYYXYX';
%       e.tau        = <divided by 8 relative to a Hahn echo over the same total time>
%
% HOW THE CLASS BUILDS IT. LoadAWG makes a piX waveform and a piY waveform -- identical except
% for the +pi/2 in the cosine argument -- and blockFor lays out
%
%       pi/2  --  [ tau  pi  tau ] x cycles  --  pi/2
%
% so the spacing between adjacent pi pulses is 2*tau and the first and last are tau from the
% ends. That is the standard XY8 timing. pulsePhase is a REPEATING pattern truncated to cycles:
% 'XYXYYXYX' at cycles = 8 gives that pattern once, and at cycles = 16 gives it twice, which is
% XY16. A shorter string still works -- 'XY' at cycles = 8 gives XYXYXYXY -- but that is CPMG-2
% repeated, not XY8, and it does not have XY8's cancellation of pulse-error accumulation.
% Uppercase only: checkPhases rejects anything that is not 'X' or 'Y'.
%
% TAU IS THE HALF-SPACING, NOT THE TOTAL. Total evolution = 2*tau*cycles, which is what
% mCurrentXAxisParam and echoFitRobust use. So at cycles = 8 a tau of 15 us is 240 us of
% evolution, and a tau list for XY8 covering the same total range as a Hahn echo is the Hahn
% list divided by 8.
%
% TWO HARD LIMITS ON TAU, both from the class rather than from physics:
%   * tau is quantised to baseWFtime (~0.15 us) and the minimum is twice that, because
%     zeroSegNum = tau/baseWFtime - 1 has to be a positive integer. LoadAWG rounds tau up and
%     removes duplicates, and it prints what it did.
%   * so the shortest total evolution at cycles = 8 is about 16 x 0.30 = 4.8 us, on a grid of
%     about 2.4 us. Do not bother asking for finer.
%
% THE GATE IS THE REAL CONSTRAINT, AND IT SCALES WITH CYCLES. The fixed MW gate is
% 2*max(tau)*cycles + 2*halfPiTime, so at cycles = 8 it is EIGHT times a Hahn echo's for the
% same tau. echoGateTest validated the gate as harmless out to 240 us of excess; the tau list
% below stops exactly there on purpose. To run XY8 further -- and you will want to, since
% decoupling is the point -- re-run echoGateTest with cycles = 8 first rather than extrapolating
% a bound measured at cycles = 1.
%
% The digitizer and maxRunSeconds ceilings both tighten at cycles = 8 as well, since the
% sequence is 8x longer. choosePassesPerRun works them out and prints which one binds.
e = AWGExpEcho;
e.loadJsonParams('all', jsonLocation);
e.nChannels = 2;                    % BEFORE amplitude: set.nChannels re-enters
e.amplitude = -16 + [0 0];          % set.amplitude, which would add amplitudeDiff twice
e.timeDelay = 0;                    % channelPhase is the phase knob; they are additive

e.repeats          = 2000;
e.averages         = 2;
e.balancedSequence = 1;
e.constantTime     = 1;
e.smallDelay       = 1;
e.isTracking       = 0;
e.shouldAutosave   = false;
e.laserInitializationDuration = 16;

e.cycles           = 8;
e.pulsePhase       = 'XYXYYXYX';    % XY8 = XY4 then XY4 reversed. NOTE: the class's own property
                                    % comment says 'XYYX' for XY4, which is not the standard
                                    % XY4 = 'XYXY'. The string here is the standard XY8; if the
                                    % lab means something else by XY4, that comment and this line
                                    % disagree and the comment is the one to check.

% 15 us of tau is 240 us of evolution at cycles = 8, which is exactly the gate excess
% echoGateTest has validated. Log spaced because the T2 knee is what is being measured.
e.tau = unique(round(linspace(1, 15, 120), 3));
fprintf('XY8 tau: %d points, %.2f to %.2f us -> evolution %.1f to %.1f us, gate %.1f us\n', ...
        numel(e.tau), min(e.tau), max(e.tau), ...
        2*min(e.tau)*e.cycles, 2*max(e.tau)*e.cycles, ...
        2*max(e.tau)*e.cycles + 2*e.halfPiTime);

e.useSequenceWalk    = true;
e.passesPerRun       = [];
e.walkWarmupTriggers = 1;
e.resyncEveryPass    = true;
e.maxRunSeconds      = 30;

e.AWG.resetCommandCount; tId = tic;
e.run;
fprintf('XY8: %.1f s, %d SCPI commands\n', toc(tId), e.AWG.commandCount);

% echoFitRobust reads cycles off the object, so the x axis and the revival period are handled.
% ONE CAVEAT ON THE IMPLIED FIELD: the period-to-Larmor conversion in that file assumes the
% Hahn-echo relation (periodic in tau with period 1/f_L, hence 2*cycles/Trev on a total-time
% axis). Under XY8 the 13-C features sit where the PULSE SPACING resonates, which is not the
% same mapping. The fitted PERIOD is still what was measured and T2 does not depend on the
% conversion -- but do not read the gauss number as a field measurement here.
o = echoFitRobust(struct('exp', e));


%%
t = ExpT1;
t.loadJsonParams('all', jsonLocation);
t.tau = linspace(1,3000,100);
t.constantTime = 1;
t.run;

%o = t1FitRobust(struct('exp', t));
%There is also ExpExtendedT1 where one can measure 1 relaxation and -1 but
%some variable needs to be defined for that
%% T1  --  relaxation. NO WALK MODE, AND IT DOES NOT NEED ONE.
% Doesnt work yet.



% Same keepOn caveat.
%
% WHY THERE IS NO useSequenceWalk HERE. Walk mode exists to remove the ~0.6 s AWG re-arm that
% the classic path pays at every point. T1 does not pay it: AWGExpT1.perform puts assignSequence
% and Run OUTSIDE the tau loop, so the AWG is armed ONCE for the whole sweep. It can do that
% because T1's waveform set is tau-independent by construction -- three waveforms per channel
% and one four-step sequence, none of which change with tau -- since the wait is a
% pulse-generator duration and not an AWG waveform. Inside the loop there is no AWG call at all,
% only pg.changeSequence and getRawData.
%
% So the thing walk mode deletes has already been deleted, and porting the machinery would add
% walkOrder, walkUnscramble, two passesPerRun ceilings and a warm-up trigger to save nothing.
% For scale: one average below is about 118 s of pure sequence time against a few milliseconds
% per point of pulse-streamer upload. The run is physics, not overhead.
%
% WHAT DOES MATTER HERE IS THE GATE. AWGExpT1 held the microwave switch open for the WHOLE of
% tau, where the AWG only plays for 0.335 us after each trigger. At tau = 4900 that is 4900 us of
% open switch against 0.3 us of need -- twenty times the largest excess echoGateTest ever
% bounded, and twelve hundred times Ramsey's measured-harmless case.
%
% And on a T1 the artefact has the SHAPE OF THE SIGNAL: anything leaking removes contrast as
% exp(-r*tau), which just adds to the relaxation rate. 1/T1_measured = 1/T1_true + r. No kink,
% no bad residual, only a T1 that is too short with a small error bar. Echo's fixed gate at
% least put its worst excess at short tau, against the decay.
%
%   >> testT1Gate      % ~30 s, no hardware. Sequence arithmetic and the statistics.
%   >> t1GateTest      % ~25 min on hardware. Measures the leak rate three ways.
%
% mwGateMode is left at 'tau' below -- unchanged -- until t1GateTest has run, because switching
% it would silently change every T1 this setup measures and there is no measurement yet saying
% which way is better. After the test: 'pulse' if it finds a leak, and 'pulse' anyway on the
% noise argument, since an open switch passes amplifier noise for milliseconds with no upside.
t = AWGExpT1;
t.loadJsonParams('all', jsonLocation);
t.nChannels = 2;                % BEFORE amplitude: set.nChannels re-enters
% -16, NOT -3. This was -3 and every other cell in this file drives -16: Rabi coarse, Rabi
% fine, Ramsey and Echo. The pulse times in the json were calibrated at that amplitude, so
% driving T1 13 dB harder makes the Rabi frequency 10^(13/20) = 4.47x higher and turns the
% calibrated pi pulse into a 4.47*pi one. Population transfer is sin^2(theta/2), which at
% 4.47*pi is 0.45 -- so under half the inversion, from the steepest part of the curve where any
% calibration error is amplified rather than averaged out.
%
% That is the leading explanation for the first real T1 run coming back with 0.13% contrast
% where Echo gets 0.3%, and for the fit reporting no measurable decay.
%
% CHECK IT RATHER THAN TRUSTING THIS COMMENT. loadJsonParams has just run, so t.amplitude at
% this moment IS the amplitude the pulse times were calibrated at:
%
%   >> t = AWGExpT1; t.loadJsonParams('all', jsonLocation); t.amplitude
%
% If that prints -16, the line below is right. If it prints something else, use what it prints.
t.amplitude = -16 + [0 0];      % set.amplitude, which would add amplitudeDiff twice
t.isTracking     = 0;
t.GIMeas         = 0;
t.AWGorSRSswitch = 0;
t.lastDelay      = 1;
% ---- SHOTS. This is the whole reason the first two runs measured nothing. ----------------
% The Rabi at these settings gives 0.88% PEAK-TO-PEAK, which is 0.44% peak-to-plateau -- the
% convention a T1 contrast is quoted in. So there is about 0.4% of signal available here, and
% resolving it needs error bars well under that.
%
%   Rabi:  repeats 3000 x averages 2 = 6000 shots per point
%   T1:    repeats  500 x averages 3 = 1500          -- 4x fewer, so 2x the sterr
%
% At 1500 shots the fitted decay spanned 1.9 error bars, against the 5 that t1FitRobust needs
% before it will call anything measured. Error bars go as 1/sqrt(shots), so 1.9 -> 6.6 wants
% 12x the shots: 3000 x 6 = 18000. Half the tau points pays for a third of that, and dropping
% constantTime pays for the rest.
t.repeats        = 3000;
t.averages       = 6;           % 18000 shots per point, 12x the first attempt
t.shouldAutosave = false;

% ---- constantTime OFF, deliberately, and watch the control ------------------------------
% With it ON every point costs max(tau), so the sweep costs nTau*maxTau. OFF, each point costs
% its own tau, which on the log grid below is 4.6x cheaper -- 10 minutes against 43.
%
% What it buys in time it risks in systematics: the laser duty cycle then varies across the
% sweep, so the steady-state temperature and charge state vary WITH TAU, which is a common-mode
% change that looks exactly like a decay. That is precisely what the |0> control measures, and
% it is why t1FitRobust reports flatnessSigma at all.
%
% So: run it off, and read the control FIRST. It came back at 3.0 sigma on the constantTime-ON
% run, which is already at the threshold. If turning this off pushes it past about 3, put it
% back to true and pay the 43 minutes -- a fast wrong answer is worth nothing here.
t.constantTime   = false;
t.mwGateMode     = 'tau';       % 'tau' (as shipped) | 'pulse' (closed during the wait) | 'all'

% 12 log-spaced points from 100 us to 6000. The first T1 attempt used 24 points of (1:3:70).^2,
% which spends a third of them below 500 us where nothing has happened yet. With T1 around 3 ms
% the sweep only needs to reach about 2 T1 -- past 1.6 T1 the quoted error bar starts tracking
% the real scatter (see t4 in testT1Fit) -- and 12 points is plenty to fit two parameters.
t.tau            = unique(round(logspace(log10(100), log10(6000), 12)));
t.detectionDuration = 0.5;
t.referenceDetectionDuration = 0.5;
t.run;

% t1FitRobust, the fifth fitter on the one-struct-in / one-struct-out convention.
%
%   >> testT1Fit       % ~40 s, no hardware. 25 checks.
%
% THREE NUMBERS TO READ FIRST, and none of them is T1.
%
%   |0> control flatness. The |0> trace was never excited, so it must not move. If it does,
%   something common-mode is changing with tau -- laser power, tracking, charge state -- and it
%   enters BOTH traces with the same shape, so the fit absorbs it into T1 and the residual stays
%   clean. This is the only thing that catches it.
%
%   span/T1. Below about 2 the plateau is an extrapolation and the ERROR BAR DOES NOT KNOW:
%   measured over ten noise realisations at 0.8 T1 of reach, the fitted T1 is only 4% low but the
%   scatter between seeds is three times the quoted error. Act on span/T1, not on the error bar.
%   The tau list above reaches 4900 us, so this is fine for a T1 up to ~2 ms and not beyond.
%
%   p, the stretch exponent, which is FITTED rather than pinned at 1. p below 1 means a
%   distribution of rates rather than one rate, and a single exponential forced onto that returns
%   a T1 that depends on where the sweep stopped. Compare out.ratePerUs between runs -- rates are
%   what add -- or out.coherenceArea when p has moved.
o = t1FitRobust(struct('exp', t));

% READ out.measured BEFORE out.T1. On the first real run of this cell the fit returned
% T1 = 9.8e4 us with no error bar, p on its upper bound and a contrast of 103% -- all of which
% are just where the optimiser stopped on a flat likelihood, because the decay spanned 1.2 error
% bars and the quoted sigma was larger than the entire spread of the trace. The fit says so now,
% loudly, but the habit worth having is to check that it CAN measure before reading what it
% measured.
%
% AND CHECK CONTRAST BEFORE T1. Contrast is the thing the pulse and the readout control, and it
% is diagnostic: Echo on this setup gets about 0.3%, so a T1 -- a straight population difference
% with no coherence to lose -- should get at least that. Materially less means the pi pulse or
% the readout window, not the sample.
if ~o.measured
    fprintf(['\nT1: NOT MEASURED -- read the reasons the fit printed. Before touching repeats,\n' ...
             '    check the contrast (%.3f%%): Echo gets ~0.3%% here, so a T1 should too. If it is\n' ...
             '    well under that, the pi pulse is the first suspect -- run the Rabi sanity cell\n' ...
             '    at THESE settings and confirm piTime against it.\n'], 100*o.contrast);
else
    fprintf('T1 = %.0f +/- %.0f us, p = %.2f, rate = %.4g /us, contrast %.3f%%\n', ...
            o.T1, o.T1Err, o.p, o.ratePerUs, 100*o.contrast);
    if o.controlTestable
        fprintf('     |0> control %.1f sigma from flat\n', o.flatnessSigma);
    else
        fprintf('     |0> control NOT TESTABLE -- see the report\n');
    end
end


%% ---- diagnostics ------------------------------------------------------
% Only when something looks wrong. These live in
% helperFunctions\setup 7 run scripts\_diagnostics\.
%
%   esrWalkValidate(a)        walk against the per-point path, judged against the noise
%                             ceiling rather than a flat correlation threshold
%   esrWalkCacheCheck(a)      warm vs cold, to see what the load cache is saving
%   rabiWalkCompare(...)      classic vs walk at identical settings, re-arm arithmetic
%                             printed alongside
%   rabiStageTimer(...)       where a Rabi run's wall clock actually goes
%   awgDrainErrors(awg)       the SCPI error queue. Failures there are silent by design --
%                             fwrite cannot see an instrument-side rejection
%   awgTransportBench(awg)    after changing useVisadev or commDelayTime
%   scopeProbe                what is actually on the cable
%
% Are the MW switches gating? Run the sweep with the switch channels undriven:
%   a.driveMWSwitches = false; a.restart; a.driveMWSwitches = true;
% A FLAT spectrum means they gate correctly. An unchanged spectrum means they do not, and
% the AWG's zeros are the only gate you have.
%
% Do NOT trust anything remembered from rabiTimingProbe before 13 August. It used bare
% tic/toc -- one global stopwatch that any tic inside b.run resets -- and reported 1.9 s
% for a 68 s run. Every conclusion drawn from it was voided; not all were re-measured.
