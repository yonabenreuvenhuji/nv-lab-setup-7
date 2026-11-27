%% ========================================================================
%  LEVEL 3  --  the magnetic field, and the operating frequency
%
%  Each step needs the previous:
%
%    3.1  offline graders     no hardware. Run them FIRST, every time.
%    3.2  calibCoilMatrix     ~20 wide ESRs at spread coil currents -> the 3x3 coil matrix:
%                             B = M*I + b0, with B in the NV CRYSTAL FRAME
%    3.3  align               B along an NV axis is now a SETTING (theta 54.7356, phi 45),
%                             not a scan; one survey checks it
%    3.4  calibB1Rabi         scan |B| so the aligned line sweeps the antenna response,
%                             and find where B1 is largest
%
%  WHAT CHANGED (2026-11). The old 3.3 scanned theta and phi until one NV axis lined up,
%  because nothing said where the crystal axes were relative to the coils -- the coil model
%  was three hand-entered ratios with no cross-talk, and the diamond's placement was typed
%  in. calibCoilMatrix measures both at once against the NVs: column k of M is the field coil
%  k makes per Ampere, expressed in the crystal frame the NVs define. With that, every NV
%  axis is at a known (theta, phi), every coordinate's limits are exact, and "align" is
%  h.alignToNV. calibFieldAlign is still here as an independent check, not a step.
%
%  ANGLES are the physics convention: Btheta POLAR from +z, Bphi AZIMUTH from +x. (The class
%  used the opposite naming before; old saved states held only B, so they are unaffected.)
%
%  WHERE THINGS GO. Everything lands in expProperties.json of the newest results folder
%  (setup7Paths().json): "frequency" (3.5, a top-level number, read by every experiment), and
%  the "Helmholtz" block -- currents, field, coil calibration -- written by h.save (3.2b,
%  3.3, 3.4c). Experiments ignore the block (no experiment has a property called Helmholtz);
%  the Helmholtz class reads it. helmholtzState.mat is no longer written; an older folder
%  that has one still loads.
%  ========================================================================


%% 3.1  offline graders  --  no hardware, run before anything else
%
% testFieldPipeline: the field survey chain against synthetic spectra with known fields.
% testCoilMatrixFit: the coil matrix against known coils (a turned diamond, skewed coils,
%   sensitivities off nominal), including the number that matters: aim at an NV axis with
%   the fitted model and see where the TRUE field points. fullChain adds synthetic spectra
%   through esrFitRobust and calibCoilMatrix's replay path (~1 min).
% testHelmholtzLimits: the field limits against brute force, and the angle convention.
% Expect ALL CHECKS PASSED from each.
testFieldPipeline
testCoilMatrixFit(struct('fullChain', true))
testHelmholtzLimits


%% 3.2  coil calibration  --  dry run first
%
% Prints the ~20 current settings (spread over the supplies' range, kept where the model in
% use predicts 12-40 G), the sweep, and what each point costs. Nothing moves.
h = getObjByName('Helmholtz');  h.control('auto');
cc = calibCoilMatrix(struct('dryRun', true));


%% 3.2a  the real coil calibration
%
% One wide CW ESR per current setting, lines fitted by esrFitRobust, then coilMatrixFit puts
% one linear model through every line of every spectrum. The raw spectra are saved after
% each point (results\<date>\level 3\coilMatrix_*.mat), so a crash costs one point, and
%     cc = calibCoilMatrix(struct('replay', cc));         % or the .mat path
% re-fits without the hardware. Restores the currents when done.
%
% Low amplitude on purpose: only line POSITIONS are used, and power broadening is what blends
% lines (below 2850 MHz they came out at 18-23 MHz at higher power).
h = getObjByName('Helmholtz');  h.control('auto');
cc = calibCoilMatrix(struct( ...
    'nPoints',      20, ...
    'Brange',       [12 40], ...
    'amplitude',    -30, ...
    'repeats',      100, ...
    'averages',     3));


%% 3.2b  use it, and write it into the json
%
% Only when cc.ok -- the fit refuses (ok = false) when more than a quarter of the lines
% disagree with every prediction, which means it did not find a consistent crystal
% orientation. The currents do not move; what changes is what the field is known to be.
if cc.ok
    h = getObjByName('Helmholtz');
    h.setCoilCalibration(cc.calibration);
    h.save(setup7Paths().json);
else
    fprintf('\nNot applied: %s\n', cc.why);
end


%% 3.3  align  --  a setting now, not a scan
%
% cc printed which NV directions the supplies can reach and over what |B| (they are
% unipolar, so usually only one or two). k indexes Helmholtz.nvDirections:
%   1 [+1 +1 +1]   2 [-1 +1 +1]   3 [-1 -1 +1]   4 [+1 -1 +1]   (theta 54.7356, phi 45 135 -135 -45)
%   5-8 the opposites (theta 125.2644)
% h.nvTable() lists the same with the reach at the present calibration. In the GUI, the
% "Align B" panel does this without arithmetic: pick a direction, then type |B| OR where the
% outer (or inner) dips should land -- the rest follow -- check the dashed preview on the
% expected-ESR plot, and press Go. Cube axes (+-x, +-y, +-z) put all four orientations
% equal; from a script: h.setDirection([0 0 1], 25) (closest reachable if just out of reach).
h = getObjByName('Helmholtz');  h.control('auto');
h.alignToNV(1, 30);           % 30 G along [+1 +1 +1]: theta 54.7356, phi 45
h.save(setup7Paths().json);


%% 3.3b  check the alignment with one survey
%
% An aligned field shows four lines: the aligned pair split by 2*gamma*|B| and the three
% off-axis orientations degenerate. The survey reports the tilt (meaningful above ~2.5 deg)
% and the inner-span excess, which keeps shrinking all the way to zero tilt. Its |B| should
% match h.Brho to the calibration's accuracy (~0.05 G).
sv = calibFieldSurvey(struct('freqs', 2760:0.5:2980, 'amplitude', -35, 'repeats', 100));
fprintf('\n|B|: survey %.3f G, coil model %.3f G\n', sv.B, h.Brho);

% ONLY if the survey disagrees (tilt > ~0.5 deg, or the inner span clearly open): the
% calibration is off, and the cure is to re-run 3.2 -- not to scan. calibFieldAlign still
% works as an independent check (theta/phi are now the polar/azimuth angles in the NV frame):
%   a = calibFieldAlign(struct('rounds', 1, 'thetaSpan', 2, 'phiSpan', 2, 'nPoints', 7));


%% 3.4  offline grader for the B1 scan  --  again, before the hardware
%
% Five tests: recover a known resonator, refuse when the peak is not in the data, and quote
% error bars that cover the truth as often as they claim to. Expect ALL TESTS PASS.
tb = testB1RabiChain;


%% 3.4b  B1 maximum  --  the operating frequency and field
%
% Scans |B| so the aligned line moves across the antenna response, running an ESR at each
% field to find the line and a Rabi at that line to measure B1. The Rabi frequency IS B1 --
% no saturation model in between, which is why this replaced the linewidth method. Brho is
% now the calibrated |B| and the direction stays on the NV axis, because setBrho scales
% along the present direction.
%
% rhoCoarse spans 12 G because the Lorentzian's WIDTH needs about 38 MHz of sweep before it
% is determined; narrower still finds the peak, via the parabola fallback, but returns a
% biased Q. minR2 is 0 by default and set here only if you want a floor -- R^2 is capped by
% the signal-to-noise ratio, so 0.9 would reject sound fits for being quiet.
%
% The second panel now also draws where the coil calibration EXPECTS the tracked line, and
% prints "measured - expected". The run is saved to results\<date>\level 3\calibB1Rabi_*.mat;
% to look at it again without the hardware (the field is not moved):
%     out = calibB1Rabi(struct('replay', out));      % or the saved .mat path
h = getObjByName('Helmholtz');  h.control('auto');

out = calibB1Rabi(struct( ...
    'rabiAmplitude',        -16, ...
    'rhoCoarse',            34:1:46, ...
    'repeats',              300, ...
    'averages',             3, ...
    'esrFreqs',             2720:3:3020, ...      % must reach D + gamma*|B| at the TOP field
    'esrAmplitude',         -25, ...
    'esrRepeats',           50, ...
    'esrAverages',          3));
%% 3.4c  write the field magnitude
h = getObjByName('Helmholtz');
h.setBrho(out.rhoOptimum);
h.save(setup7Paths().json);
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
h = getObjByName('Helmholtz');
fprintf('\nfrequency in expProperties.json:  %.2f MHz\n', out.resonantFrequency);
fprintf('Brho %.3f G, theta %.4f (polar), phi %.4f (azimuth) -- %s\n', ...
        h.Brho, h.Btheta, h.Bphi, h.frameName);
[kNV, cNV] = h.nearestNV();
fprintf('%.3f deg from NV direction %d\n', acosd(min(cNV, 1)), kNV);
if isfield(out, 'piNsAtOptimum') && isfinite(out.piNsAtOptimum)
    fprintf('pi pulse there: %.1f ns (both channels)\n', out.piNsAtOptimum);
end
fprintf(['\nThe field and the coil calibration are in the "Helmholtz" block of\n' ...
         'expProperties.json (3.2b, 3.3, 3.4c); the frequency is the top-level "frequency",\n' ...
         'which loadJsonParams hands to every experiment on construction.\n']);
