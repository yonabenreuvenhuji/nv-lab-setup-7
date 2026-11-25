%% ========================================================================
%  LEVEL 3  --  the magnetic field, and the operating frequency
%
%  Four steps, and each one needs the previous:
%
%    3.1  testFieldPipeline   offline grader, no hardware. Run it FIRST, every time.
%    3.2  calibFieldSurvey    one spectrum -> the field vector [Bx By Bz] and D
%    3.3  calibFieldAlign     turn theta and phi until one NV axis is aligned with B
%    3.4  calibB1Rabi         scan |B| so the aligned line sweeps the antenna response,
%                             and find where B1 is largest
%
%  WHY THE OFFLINE TEST COMES FIRST. Both of the worst bugs in this chain were invisible on
%  real data and obvious against a known answer. The multi-line fitter could not split a
%  blend, because its seeded amplitudes came from a movmax narrower than the lines, so every
%  seed started at the noise floor. And the survey threw away correctly aligned pairs,
%  because the mirror-partner test assumed a single D that pairs all four orientations --
%  which does not exist: the off-axis pair is pushed up 7-8 MHz at 45 G by the second-order
%  transverse shift while its splitting barely moves. Neither had a symptom you could see in
%  a spectrum. testFieldPipeline takes seconds and needs nothing switched on.
%
%  Writes to expProperties.json: frequency (from 3.4b), through expPropsWrite in 3.5.
%  NOT "resonantFrequency" -- that key does nothing at all; see 3.5.
%  Writes to <results>\helmholtzState.mat: Btheta, Bphi (3.3) and Brho (3.4),
%  where <results> is setup7Paths().folder -- the newest dated calibration folder.
%  ========================================================================

%% 3.0  the CW-ESR window  --  once per session, before any hardware
%
% Every ESR in this level goes through calibFieldSurvey, and all of them used to run at a
% 10 us window. That was wrong twice over:
%
%   TOO SHORT FOR THE PHYSICS. tau_pol is 13.6 us on this setup, so a 10 us window is less
%   than ONE polarisation time. Laser and microwave both start at the top of every CW
%   window, so its first ~tau_pol carries less than steady-state contrast -- at 10 us the
%   dip comes back at roughly a third of its available depth. That is a 3x SNR penalty on
%   every line centre in this level, for nothing.
%
%   AND POSSIBLY ILLEGAL. In CW ESR the AWG loops a burst for detectionDuration +
%   laserInitializationDuration while one repeat is 2*(detectionDuration +
%   delayBetweenDetections). If the burst is longer, triggers are dropped, those repeats get
%   no microwave, and the spectrum comes back as a COMB -- alternate points at ratio 1 with
%   the dip only on the others. It does not error and it does not warn.
%
% esrWindowPick settles both. calibFieldSurvey now also resolves the window through it when
% detectionMHz is left empty, and clamps laserInitializationDuration either way, so a bare
% call is safe -- but computing it here once means the number is visible in the log and the
% same value is used by 3.2, 3.3 and 3.4.
%
% tauPol is read from expProperties.json, where level 2 (cell 2.3) writes it, so this cell no
% longer needs level 2's r in the workspace -- which also removes the old trap of 3.4
% reassigning r and breaking a re-run of this cell. calibFieldSurvey calls esrWindowPick the
% same way when detectionMHz is empty, so it now gets the transient floor too; before, it
% never passed tauPol and the floor was silently skipped.
cw = esrWindowPick(struct('jsonLocation', setup7Paths().json));
esrWindow = cw.detectionDuration;


%% 3.1  offline grader  --  no hardware, run before anything else
%
% Three tests against synthetic spectra with known fields: does the chain recover |B| and
% the tilt, does the alignment scan find the vertex, and does it survive a spurious non-NV
% line and a clipped sweep. Expect ALL CHECKS PASSED.
testFieldPipeline


%% 3.2  field survey  --  where are the four orientations
%
% A wide sweep, wide enough to hold both lines of every orientation. The fit solves
% [Bx By Bz D] from the line centres using the exact 3x3 Hamiltonian rather than pairing
% lines by hand -- see the header note above for why pairing does not work.
h = getObjByName('Helmholtz');  h.control('auto');

% detectionMHz was ABSENT here, so this cell ran at the old 10 us default -- see 3.0. That
% is the single change most likely to fix a survey whose line centres look off: the dip was
% a third as deep as it should have been, so every centre had three times the error bar it
% needed, and a shallow blend is exactly where the multi-line fit gives up first.
sv = calibFieldSurvey(struct('freqs', 2795:0.5:2935, 'amplitude', -35, ...
                             'detectionMHz', esrWindow, 'repeats', 100));

% To re-fit the SAME data with different fit options, pass it back in rather than
% re-measuring. 'replay' also accepts an experiment object.
%   r = calibFieldSurvey(struct('replay', sv));


%% 3.3  alignment  --  turn theta and phi until one axis lines up with B
%
% Dry run first. It reports what it would scan and how long it would take without moving
% the coils, which is worth ten seconds given that the real scan takes a while and a wrong
% sweep window wastes all of it.
h = getObjByName('Helmholtz');  h.control('auto');

a = calibFieldAlign(struct('dryRun', true, ...
    'nPoints',     15, ... 
    'amplitude',   -35, ...
    'scanFreqs',   2800 : 0.25 : 2930, ...
    'averages',    1, ...
    'repeats',     100, ...
    'detectionMHz', esrWindow));    % from 3.0. Was 50, described as "the CW window from
                                    % level 2" -- that came from calibCwWindow, which is
                                    % gone, and 50 us is now below the trigger limit anyway.

%% 3.3b  the real alignment scan
%
% The metric is the squared excess span of the inner pair, not the tilt: below about 2.5
% degrees the tilt estimate collapses toward zero because the off-axis splittings differ by
% less than a linewidth, while the span keeps moving all the way to zero tilt.
h = getObjByName('Helmholtz');  h.control('auto');
% h.setBrho(25);

a = calibFieldAlign(struct(...
    'rounds',               2, ...
    'tolDeg',               0, ...
    'thetaSpan',            15, ...
    'phiSpan',              15, ...
    'nPoints',              13, ... 
    'amplitude',            -30, ...
    'scanFreqs',            2760:1:2980, ...
    'detectionMHz',         esrWindow, ...
    'repeats',              100, ...
    'newOrderEveryAverage', true));

%% 3.3c  write the alignment
h = getObjByName('Helmholtz');
h.setBtheta(a.theta);   h.setBphi(a.phi);
h.save(setup7Paths().helmholtz);
fprintf('field direction saved: theta %.4f, phi %.4f deg\n', a.theta, a.phi);


%% 3.4  offline grader for the B1 scan  --  again, before the hardware
%
% Five tests: recover a known resonator, refuse when the peak is not in the data, and quote
% error bars that cover the truth as often as they claim to. Expect ALL TESTS PASS.
%
% Named tb, not r: r is level 2's calibDetectionWindow result and 3.0 needs r.tauPol out of
% it. Overwriting it here worked only because 3.0 had already run.
tb = testB1RabiChain;


%% 3.4b  B1 maximum  --  the operating frequency and field
%
% Scans |B| so the aligned line moves across the antenna response, running an ESR at each
% field to find the line and a Rabi at that line to measure B1. The Rabi frequency IS B1 --
% no saturation model in between, which is why this replaced the linewidth method.
%
% rhoCoarse spans 12 G because the Lorentzian's WIDTH needs about 38 MHz of sweep before it
% is determined; narrower still finds the peak, via the parabola fallback, but returns a
% biased Q. minR2 is 0 by default and set here only if you want a floor -- R^2 is capped by
% the signal-to-noise ratio, so 0.9 would reject sound fits for being quiet.
h = getObjByName('Helmholtz');  h.control('auto');

out = calibB1Rabi(struct( ...
    'rabiAmplitude',        -16, ...
    'rhoCoarse',            28.5:1:38.5, ...
    'repeats',              300, ...
    'averages',             3, ...
    'esrFreqs',             2750:2:2992, ...      % must reach D + gamma*|B| at the TOP field
    'esrAmplitude',         -25, ...
    'esrRepeats',           50, ...
    'esrAverages',          3, ...
    'esrDetectionMHz',      esrWindow, ...
    'esrNewOrderEveryAverage', true, ...    % ESR steps: new frequency order every average
    'newOrderEveryAverage',    true));      % Rabi runs: new tau order every average
%% 3.4c  write the field magnitude
h = getObjByName('Helmholtz');
h.setBrho(out.rhoOptimum);
h.save(setup7Paths().helmholtz);
fprintf('Brho saved: %.3f G\n', out.rhoOptimum);


%% 3.5  write the operating frequency
%
% THE JSON KEY IS "frequency", NOT "resonantFrequency". Worth being exact about, because
% getting it wrong fails silently in both directions: the AWG experiment classes declare a
% property called `frequency`, and loadStateFromStruct loops over the OBJECT's properties
% asking whether the json has a field of that name. A json field with no matching property is
% skipped with no error and no warning -- so a "resonantFrequency" entry would sit in the
% file looking authoritative while every experiment ran at its constructor default (2760 MHz
% in AWGExpRabi). This cell used to tell you to write exactly that.
%
% It is also why the old resonant_frequency.mat and the `b.frequency = resonantFrequency`
% lines are gone from the run script: the json supplies it on construction now.
%
% Writes only when 3.4b settled (out.ok). An unsettled scan says "nothing downstream should
% use these numbers", and this cell would otherwise move every experiment's frequency to one.
% To see the diff first, add 'dryRun', true.
if out.ok
    w = expPropsWrite(struct('frequency', out.resonantFrequency), ...
                      struct('jsonLocation', setup7Paths().json));
else
    fprintf(['\n3.4b did NOT settle (see its hand-over block), so the frequency has not been\n' ...
             'written. Extend the field scan and re-run 3.4b. To write anyway:\n' ...
             '   expPropsWrite(struct(''frequency'', out.resonantFrequency), ...\n' ...
             '                 struct(''jsonLocation'', setup7Paths().json));\n']);
end


%% 3.6  check before moving up
fprintf('\nfrequency in expProperties.json:  %.2f MHz\n', out.resonantFrequency);
fprintf('Brho %.3f G, direction theta %.4f phi %.4f\n', ...
        out.rhoOptimum, out.theta0, out.phi0);
if isfield(out, 'piNsAtOptimum') && isfinite(out.piNsAtOptimum)
    fprintf('pi pulse there: %.1f ns (both channels)\n', out.piNsAtOptimum);
end
fprintf(['\nThe field DIRECTION and MAGNITUDE stay in helmholtzState.mat (3.3c, 3.4c) --\n' ...
         'they are coil settings, and no experiment class has a property for them, so the\n' ...
         'json could not deliver them even if they were in it. The frequency is the\n' ...
         'opposite: it belongs in the json, where loadJsonParams hands it to every\n' ...
         'experiment on construction.\n']);
