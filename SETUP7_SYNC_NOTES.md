# Setup 7: Drive → repo sync notes

Snapshot taken **2026-11-25, ~20:30–22:00**, from `G:\...\NV Lab` (read-only; nothing on the Drive
was changed). Someone was editing Setup 7 on the Drive while this ran
(`calibrations/level 3/level3.m` 20:24, `results/2026-11-23/helmholtzState.mat` 20:20), so
re-sync before the switch-over. Nothing here is committed yet.

**Decision taken: this branch runs the Drive's core ("option A").** `main`'s FrequencyGenerator →
SignalGenerator refactor (`add_AWG`) is replaced by the Drive's architecture on this branch, so
setup 7 runs the same code it runs today from the Drive. This branch therefore diverges from `main`
in the shared core. See section 5 for what porting to `main`'s format later would take.

---

## 1. What this branch now contains

### Setup-7 code (setup-7-only)

**New folder `Setup 7/` at the repo root**, a mirror of `NV Lab/Setup 7` (194 files):
`functions/` (incl. `classes/`), `experiments/` (incl. `_diagnostics/`), `calibrations/`
(`level 0`–`level 4`, `waist`, `CALIBRATION_PLAN.md`, and `results/<date>/` with the
`expProperties.json`, `helmholtzState.mat` and other `.mat` data the code reads at runtime, plus
the 48 `.svg` calibration plots, ~52 MB, largest 5.2 MB).

- It sits outside `src/` on purpose. `setup7Paths` / `setup7AddPath` derive the root from their
  own location and add exactly the folders they want. `functions/classes` and
  `calibrations/results` are deliberately *not* on the path, and an `addpath(genpath('src'))`
  would put them there.
- Not copied: `old files/` (20.6 GB archive, deliberately off-path), the 12
  `expProperties.json.bak.*` backups in `results/`, `*.asv`, `desktop.ini`.
- Re-checked against the live Drive after the SVG copy: no new, removed or changed files.

**`src/` setup-7 files at the Drive version:**

| Area | Files |
|---|---|
| AWG base + driver | `Experiment/ExperimentAWG.m`, `Physical/FrequencyGenerator/FrequencyGeneratorTektronixAWG.m` |
| AWG experiments (changed) | `AWGExpESR`, `AWGExpEcho`, `AWGExpHelmholtzAlignment`, `AWGExpRabi`, `AWGExpRabiCont`, `AWGExpRamsey`, `AWGExpT1`, `AWGsensingSingleDriveCont` |
| AWG experiments (new) | `AWGExpESRPulsed`, `AWGExpESR_gapAom`, `AWGExpWeakMeasCorr`, `AWGPulsedDynamicalDecoupling`, `AWGsensingACfield`, `AWGEXPDRDOIDS/` (m + mat), `awgesrblock.txt` |
| Digitizer calibration | `AWGExpDigitizerDetectionDuration`, `AWGExpDigitizerInitDuration`, `AWGExpLastDelayDuration` (new) |
| Digitizer / photodiode | `Physical/Spcm/Digitizer.m`, `Physical/Spcm/PhotoDiodeDigitizerNiDaqControlled.m`, `Gui/Spcm/ViewPhotodiodeControl.m` |
| Standa stages | `Physical/Stage/StageChildren/ClassStanda8SMC4_stages.m` |

**Removed:** `src/helperFunctions/setup 7 run scripts/` (680 files). The Drive moved
`runScriptBasicExpAWG.m` to `Setup 7/experiments/`, and the old copy would shadow it.
`alignmentScript.m`, `beamWaist_calibration.m`, `parameters_optimization.m` and the
`setup calibration` data are no longer in the Drive's active folders. They remain in git history.

### Drive core restored (shared files, so setup 7 runs as it does on the Drive)

| Group | Files (all = Drive version) |
|---|---|
| Generator architecture | `Setup/Setup.m`; `Physical/FrequencyGenerator/` restored (`FrequencyGenerator`, `…Dummy`, `…RigolAWG`, `…SGT100A`, `…SRS`, `…TektronixAWGDummy`, `…Windfreak`); **`Physical/SignalGenerator/` removed** (`SignalGenerator`, `AWG/AWG`, `AWG/Waveform`, `AWG/wfSegment` and the refactored FG classes) |
| Core it is coupled to | `Experiment/Experiment.m`, `Experiment/Sequence/Pulse.m`, `Experiment/Sequence/Sequence.m`, `Physical/PulseGenerator/PulseGenerator.m`, `Physical/PulseGenerator/PulseStreamer/old/PH.m`, `Tracker/Trackables/TrackablePosition.m` |
| Experiments the refactor had migrated | `ExpESR`, `ExpRabi`, `ExpRamsey`, `ExpEcho`, `ExpEchoNpulse`, `ExpExtendedT1`, `ExpMagneticFieldAlignment`, `Calibrations/ExpCalibrateDelays` |
| Required by setup-7 classes | `Physical/Spcm/Spcm.m` (`hasDualInput`), `Physical/Stage/ClassStage.m` (`pgControlled` / trigger methods), `Experiment/NonExperiments/Scan/StageScanner.m` (`tPixel`), `Physical/PulseGenerator/PulseStreamerNewClass.m` (timeout) |
| On setup 7's call path | `Experiment/Experiments/ExpT1.m` (the run script uses it), `Gui/Util/progressbar.m`, `Image/ImageScanResult.m`, `Physical/Laser/Switch/SwitchPgControlled.m`, `Helper/PathHelper.m`; fit helpers with Drive contents, kept in the repo's folders: `fits/ESRfit.m`, `fits/rabiFit.m`, `extractData/extract1Ddata.m`, `extractData/extract2Ddata.m` (the Drive `rabiFit` needs the new `lastAvg` outputs) |

**Checks run:**
- No new duplicate function/class names on the path. The 8 inherited duplicates are unchanged
  (e.g. `PulseStreamer/P.m` vs `PulseStreamer/old/P.m`).
- Every function/class that setup-7 code references matches the Drive version (dependency scan,
  3 levels deep).
- 573 files at the same path are byte-identical to the Drive (line endings ignored).
- Nothing outside `ExpExtendedCoherentControl.m` references the removed SignalGenerator API.
- **Not tested in MATLAB:** this PC isn't the lab PC.

---

## 2. Lab-PC / config to-dos (not files in the repo)

1. **`C:\lab\setupInfoLocation.json`** (read by `JsonInfoReader`) points to the hardware JSON,
   presumably the Drive's `Control code\Lab json\Setup 7\setupInfo.json`. It stays compatible
   (same `frequencyGenerators` / `microwave` format). Decide whether it stays on the Drive or gets a
   copy in the repo (e.g. `Setup 7/json/`), and repoint it. The README's "settings file is at
   C:\lab\setupInfo.json" is inaccurate: it's a pointer file.
2. **`C:\lab\setupParamsLocation.json`** is used by `Experiment.loadJsonParams` when no location is
   passed. The run scripts and `expPropsWrite` use `setup7Paths().json`, which will now be the
   repo's `Setup 7/calibrations/results/<newest>/expProperties.json`. Point this file at the same
   place, or two calibrations will be live at once.
3. **MATLAB path:** remove the Drive's `prod` and `NV Lab\Setup 7` folders from the saved path
   (`pathdef.m`), add the repo's `src`, then run `setup7AddPath` from the repo's
   `Setup 7\functions`. Run it with `struct('check', true)` first: it reports duplicate names.
4. **Keep `G:\My Drive` mounted.** `ClassStanda8SMC4_stages` loads `libximc` from
   `<NV Lab>\Control code\Drivers\ximc\ximc-2.10.5\...` via `PathHelper.getPathToNvLab`, and saves
   still go to `NV Lab\Control code\Saves\Setup 7`. Alternatively, vendor the driver.
5. **Results now live in two places.** From the switch-over, calibrations write into the repo's
   `Setup 7/calibrations/results/`. The existing `.svg` plots are included (up to ~5 MB each).
   Every calibration run adds more, so decide whether to keep committing them or `.gitignore` new
   ones.
6. `AWGEXPDRDOIDS/scrippt_weak_meas.m` hardcodes
   `G:\My Drive\NV Lab\Setup 7\calibrations\calibration results\`, a folder that no longer exists.

---

## 3. Still different from the Drive. Deliberately untouched: not on setup 7's path

These 29 same-path files differ, but setup 7 never runs them (other setups' hardware or other
people's experiments):

- Experiments: `ExpESR_Nuclear`, `ExpEchoNpulseWithACMagnetic`, `ExpInitDuration`,
  `ExpMagneticFieldAlignmentFL`.
- Camera stack: `CameraControlled`, `Camera`, `CameraBasler`, `Gui/Camera/*` (7 files),
  `CameraDisplay`, `CameraCapture`, `CameraImageParams`, `EventExtraImageUpdated`, `fmeasure`.
- Other hardware: `ClassANC350`, `ClassECC`, `ClassGalvo`, `ClassMCS2`, `SpcmNiDaqControlled`
  (keeps `main`'s DAQ-task-close fix), `SpcmTimeTaggerControlledNiDaqEnabled`,
  `ClassExternalFieldControl` (Pavel's timeout), `LaserSourcePhysicalFactory` (whitespace only).
- Other: `fourierTransform`, `setupInfoLocation.json.lnk`.
- Fit helpers not used by setup 7: `fits/EchoFit.m`, `extractData/extract2D_rabi.m` (Drive is newer).

**`main`-only, left in place:** `ExpExtendedCoherentControl.m` (uses the removed SignalGenerator API,
so it won't run on this branch), `Exp2DRabi_FT.m`, `Physical/MirrorStage/ClassMirrorStage.m`.

**New on the Drive, other people's work (not copied):** Bayesian ESR, weak-measurement experiments,
nuclear experiments, `ESRsecondChannel`, `ExpCalibrateCheckIQinScope` (Ty), `PiezoMirrors` /
`MirrorStage` additions (setup 4), `WavePlateStage`, `PDQ80Reader`,
`ArduinoGP8413Daq`/`DueDaq`/`NiDaqFather`, `LaserSourcePDL810`, `L1_homotopy_v2.0`,
`cardiac project`, `IttaiNDsearch`, `Tests/`, new fit helpers.

**Unclear owner, not copied:** `helperFunctions/analyze_AOM_saturation/*` + `aomSweep_*.mat`
(2026-11-22). Copy it if it's setup 7's.

**Backups not copied:** `Copy_of_AWGExpESR.m`, `Copy_of_FrequencyGeneratorTektronixAWG.m`,
`FrequencyGeneratorSGT100A_old.m`, `create_pulse2.m`, `windfreakconnection.m`,
`ClassStanda8SMC4_stages-1.m`, `PhotoDiodeDigitizerNiDaqControlled-1.m`,
`AWGExpHelmholtzAlignment.mat`.

---

## 4. What `main`'s refactor had that this branch drops

Only relevant if setups 4, 8 or 10 ever share this branch (they use `main`): the
`signalGenerators` JSON format, the SGT100A AWG/IQ pipeline (`AWG`, `Waveform`, `wfSegment`,
`Experiment.loadAWG`), and add_AWG fixes in the restored files, e.g. `Experiment.wrapUp` after
`prepare`, the `WindFreakChannels` validation fix, and the `PulseGenerator` single-precision delay
workaround.

---

## 5. Later: porting setup 7 to `main`'s SignalGenerator format (option B)

This doesn't require rewriting the AWG experiments. They drive the Tektronix directly through
`obj.AWG` (about 700 call sites). Needed:

1. **JSON:** a `signalGenerators` entry for the Tektronix with `fgChannels` MW (ch 5) / MW2 (ch 6).
2. **Tektronix driver:**
   - pass `numChannels` to the new base constructor;
   - the base declares `createCommand` Abstract Static, but the driver's needs `obj`, so rename it
     and add a static stub;
   - override `queryValue` (the new base doesn't send the query before reading);
   - make sure startup queries don't flip the AWG/FGEN mode.
3. **`ExperimentAWG`:**
   - find the AWG through `SignalGenerator` (`freqGenName` is gone);
   - override the protected `validateCellInput` / `setFGparams`, so the base doesn't turn
     frequency/amplitude into cells or push output to the AWG (which would switch it to FGEN).
4. **14 `FrequencyGenerator.getFG()` calls** in 8 `Setup 7` tool files.
5. **Shared, coordinate with the refactor's owners:**
   - `TrackablePosition` (`cell2mat(parent.amplitude)` fails for AWG experiments);
   - `ExpT1` not migrated in `main`;
   - merge the Drive's photodiode changes into `main`'s `Experiment.m`.
6. **Re-run the calibration ladder (levels 0–4) on the lab PC.** `main`'s
   `PulseGenerator.changeSequence` edits the live sequence in place instead of rebuilding it, and
   AWG experiments call it every tau point.

Caveat: `main`'s migration is itself unfinished. Only the SRS and SGT100A drivers are fully
ported, about 13 of 75 non-AWG experiments are migrated, and `SignalGenerator.getFG` is called but
never defined.
