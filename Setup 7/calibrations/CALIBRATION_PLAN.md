# Setup 7 calibration plan

Revision of 2026-08-14. Supersedes the 2026-07-31 plan, which predates the level 3/4
renumbering, the Rabi-based B1 and microwave calibrations, and the move out of
`Control code\prod`.

## The tree

```
Setup 7\
  functions\        library code. On the MATLAB path, and everything below may call it.
                    setup7Paths  rabiFitRobust  esrFitRobust  esrMultiLineFit
                    walkReduce  ratioAcrossAverages
  experiments\      runScriptBasicExpAWG.m, runScriptAdvancedExpAWG.m
    _diagnostics\   probes and offline tests. Nothing here decides a calibrated value.
  calibrations\     level 0 .. level 4, and results\
    results\
      2026-08-14\   expProperties.json, helmholtzState.mat
      <next one>\   a new dated folder per calibration campaign
  old files\        the archive. Nothing on the path.
```

**`walkReduce.m` must stay in `functions\`.** `AWGExpESR` and `AWGExpRabi` both call it at
runtime, so if it is not on the path walk mode breaks on both — and it used to live in a
scratch folder that happened to be on the path, which is a bad reason for something to work.

**No file writes a path down.** `setup7Paths` derives the Setup 7 root from its own location
in `functions\`, so there is no drive letter anywhere and it works whether the tree is
reached as `G:\My Drive\NV Lab\Setup 7` or through the `.shortcut-targets-by-id` path. It
picks the newest dated folder under `results\` and prints which one, so the log records what
a session was run against. `setup7Paths('2026-08-14')` pins an older one on purpose;
`setup7Paths('list')` shows what exists.

Before this, twelve files each carried their own copy of the old path as a default. That is
not twelve times the work to change — it is one opportunity per copy to change it wrongly,
and a missed copy loads the *old* calibration into one experiment while everything else uses
the new one. Two calibrations live at once, with no symptom except data that will not
reproduce.

## The order, and why it is an order

Each level's measurement is only meaningful once the ones below it are fixed. That is the
whole content of the numbering.

| Level | What it fixes | Field | Depends on |
|---|---|---|---|
| 0 | dark level, laser power | off | nothing |
| 1 | channel delays | off | 0 |
| 2 | detection windows, last delay | off | 0, 1 |
| 3 | field vector, alignment, operating frequency | **on** | 0–2 |
| 4 | two-channel balance and phase | **on** | 0–3 |

Run 0–2 after any optical or cabling change. Run 3 after moving the magnet or the sample.
Run 4 after anything in the microwave chain, and after 3, because it is measured at the field
3 chooses.

---

## Level 0 — the optical chain

`calibDarkLevel` then `calibLaserSaturation`, in that order, because the second takes the
first's result rather than re-measuring it. Nothing here touches the AWG or the coils.

Block the laser before `calibDarkLevel`. It does not check, and a dark level measured with
the laser on is a positive offset that quietly cancels part of every later power reading.

Writes nothing to `expProperties.json` — there is no dark-level or laser-power field and
nothing reads one. The AOM setting lives on the laser object; cell 0.2 prints the line to
set it.

## Level 1 — channel delays

`calibChannelDelays`, four times, with the channel list in different orders. Each call
measures against the **first** channel in the list, so `greenLaser` is the reference every
time.

Cells 1.2 and 1.3 pass the same three channels in two different orders. That is not
redundant — it is the cross-check. A delay that depends on the order it was measured in is
not a delay, it is an artefact, and neither number should be written.

Cell 1.4 does `MW2`. The two microwave channels have different path lengths, and their
relative *delay* is calibrated here. Their relative *phase at the carrier* is a level 4 job
and is not the same measurement: `timeDelay` quantises to whole AWG samples, which is 64.6°
of phase at 2870 MHz.

Writes nothing to `expProperties.json`. The per-channel delays live in `setupInfo.json`,
which is edited by hand on purpose — each call prints measured against stored.

## Level 2 — detection windows

Run in this order (see `level2.m`; offline check of both fits: `testLevel2Fits`):

1. `calibLastDelay` **fits** the gap scan with an erf step plus a slowly varying background.
   `lastDelay` is the 95% upper bound of the point where the step is 99% complete, i.e. where
   the MW has actually stopped. The background (a creep of a few %/µs, likely the MW-triggered
   60 kHz artefact) is reported, not used.
2. `calibDetectionWindow` fits the polarisation decay as a **sum of 1–3 exponentials** (BIC)
   plus a nuisance damped oscillation, corrects the amplitudes for the contaminated reference,
   and defines `tauPol` as the **mean** re-polarisation time (area under c/c0). From the same
   fit it chooses `detectionDuration` and `laserInitializationDuration` **together**: best SNR
   per unit time with carry-over into the next point ≤ 5% (walk mode). The window sits inside
   `laserInitializationDuration` in every pulsed class, so it costs no sequence time.
3. Cell 2.3 writes all of it with `expPropsWrite`. Cell 2.4 computes the **CW** window with
   `esrWindowPick`, which reads `tauPol` from the json.

**The two windows are different numbers for different measurements, two orders of magnitude
apart.** Pulsed experiments read out a prepared state and the contrast is gone within a
microsecond. CW ESR is a steady-state measurement with laser and microwave both on, where
there is no transient to catch and a long window simply collects more photons. Using the CW
number in a Rabi throws away all the contrast; using the pulsed number in a CW ESR throws
away all the photons. Neither raises an error.

Writes (cell 2.3, `expPropsWrite`): `lastDelay`, `detectionDuration`,
`referenceDetectionDuration`, `laserInitializationDuration`, `tauPol`. `tauPol` is not an experiment property, so
`loadJsonParams` skips it; it is in the json so that `esrWindowPick` (run script, level 3,
`calibFieldSurvey`) reads one measured number instead of a hand-typed copy. The CW window is
**not** written: it would collide with the pulsed `detectionDuration`.

## Level 3 — the field, and the operating frequency

Four steps: `testFieldPipeline` (offline), `calibFieldSurvey`, `calibFieldAlign`,
`calibB1Rabi`.

**Run the offline grader first, every time.** Both of the worst bugs in this chain were
invisible on real data and obvious against a known answer. The multi-line fitter could not
split a blend, because its seeded amplitudes came from a `movmax` narrower than the lines, so
every seed started at the noise floor. And the survey discarded correctly aligned pairs,
because a mirror-partner test assumed a single D pairing all four orientations — which does
not exist: the off-axis pair is pushed up 7–8 MHz at 45 G by the second-order transverse
shift while its splitting barely moves. Neither had a symptom you could see in a spectrum.

`calibFieldSurvey` solves `[Bx By Bz D]` from the line centres via the exact 3×3 Hamiltonian.
Lines are never paired by hand, for the reason just given.

`calibFieldAlign` scans θ and φ. Its metric is the **squared excess span of the inner pair,
not the tilt**: below about 2.5° the tilt estimate collapses toward zero because the off-axis
splittings differ by less than a linewidth, while the span keeps moving all the way to zero
tilt.

`calibB1Rabi` scans |B| so the aligned line sweeps the antenna response, running an ESR at
each field to find the line and a Rabi at that line to measure B1.

**The Rabi frequency *is* B1**, up to one fixed constant: Ω = γ·B1·⟨matrix element⟩, and the
matrix element moves by under one part in 5000 across the scan, since the mixing that would
change it goes as (γB⊥/D)² ≈ 4×10⁻⁴ even for a fully off-axis orientation. No saturation
model in between. That is what replaced the linewidth method, which needed one and also
carried inhomogeneous broadening and hyperfine structure that have nothing to do with drive
strength.

The response is fitted as **Ω² against measured line frequency**, not as a parabola in field:
a resonator passes power as a Lorentzian and B1 goes as the square root of power, so Ω² is
the Lorentzian. That yields f₀, the width, Q and the on-resonance Rabi frequency from the
points a parabola would have reduced to one number.

Two sizing facts, both measured on synthetic scans of a known Q = 94:

- The **coarse scan must span about ±7 G** (38 MHz). At ±5 G the fit was usable only 77% of
  the time and returned Q ≈ 101 against a true 94, because over a short arc a Lorentzian
  shoulder is nearly a parabola and (f₀, W, A) trade off. At ±7 G it is usable 100% of the
  time and Q comes back at 93.9.
- A **narrower scan still finds the peak.** The weighted parabola vertex recovers f₀ to
  better than 1 MHz down to an 8 MHz span, at spans where the Lorentzian is correctly
  refused. So a short scan is a legitimate choice when you want the field setting and not Q.

The response curve is fitted to **both passes together**. A 4 G fine pass spans ~11 MHz, far
too short for a width on its own; combined, the coarse points supply the wings and the fine
points the curvature.

`coilLineRelation` converts the target frequency back to a coil setting, and is refitted from
each scan rather than saved. **This replaces plan §4.4 (`calibFieldMagnitude`), which was
never written and should not be.** Over 15–20.6 G the slope and offset of a straight line are
0.995 correlated, so two honest fits gave `|B| = 0.974·Brho − 0.813` and `1.000·Brho − 1.176`.
Both fit their own data and agree near 18 G; neither slope is "the" coil ratio. The pair is an
interpolator over the range just scanned, and the caller checks whether the answer landed
inside it.

Writes: `frequency` in the json (cell 3.5, `expPropsWrite`, only when 3.4b settled; **not**
`resonantFrequency` — that field name never existed), and `Btheta`, `Bphi`, `Brho` in
`helmholtzState.mat` (cells 3.3c, 3.4c).

## Level 4 — the two microwave channels

`testMwPairChain` (offline), then `calibMwPair`, at the field level 3 chose. Four stages: each
channel alone, a phase scan with both driven, a closure test, an antiphase check.

**One channel at a time is valid.** Everything between the AWG outputs and the diamond is a
linear network, so the transfer function from output *n* to B1 is a fixed property of the
hardware and cannot depend on what the other port is doing. The per-channel circulators make
this stronger, not weaker — a circulator is exactly what stops channel 2's reflected power
reaching channel 1's amplifier output. And since each amplifier only ever sees its own
channel, turning the other on does not move its operating point.

**Stage 3 measures that rather than assuming it.** In a linear network at the optimum phase
Ω_both = Ω₁ + Ω₂ exactly, so R = Ω_both/(Ω₁+Ω₂) is a direct test of superposition. R = 1
means the one-at-a-time calibration is sound. A shortfall at the *fitted* maximum cannot be a
phase error, so it is compression — and a compression both channels share scales them
together, costing drive without invalidating the balance.

**The real caveat is level, not both-on.** Per-channel amplifiers compress at different
points, so `delta` measured at −30 dBm need not hold at −16. Calibrate at the amplitude you
run at; re-run stage 1 alone if you change it. The phase does not need redoing — it is set by
path lengths.

**The antiphase null is a check, not the measurement**, and the reverse would have been right
one calibration ago. With CW dip depth as the observable, comparing two channels meant
comparing two points on an unmeasured curve, so a null was the only option. The Rabi
frequency changes that: it gives the ratio to about 0.05 dB in four runs. The null cannot
match that, because a 1° phase error on perfectly balanced channels leaves the same residual
as a 1% amplitude imbalance — so the null's floor is set by how well the phase is known, not
by how well a Rabi frequency can be measured. And a slow residual cannot be seen on a sweep
shorter than its period.

The phase is fitted from the **whole** interference curve. Ω² = c₀ + c₁cos φ + c₂sin φ is
exact for two coherent drives and linear in three coefficients, so there is no start point
and no local optimum. Peak-finding would be poor: a 25° error costs only 2.4% of the drive,
so the peak is nearly flat while the null is sharp.

Writes (cell 4.3, `expPropsWrite`, only when 4.2 settled): `amplitudeDiff`, `channelPhase`,
and **`timeDelay = 0`**. Cell 4.5 then writes `calibPulseTimes`' `halfPiTime`, `piTime`,
`threeHalvesPiTime`, `amplitude` and `frequency`, unless its self-check says the drive moved.
That `timeDelay = 0` is not optional — `timeDelay` is the other phase mechanism and the two are *additive* inside the
waveform builder, so leaving both non-zero means the real phase is their sum and nothing
records it.

`calibResonator` is a **diagnostic, not a calibration.** It measures power returning from the
circulators' reverse ports, i.e. power failing to get into the antenna. That is not B1 at the
NVs, and a reflection minimum can sit several MHz from the field maximum. The operating
frequency comes from `calibB1Rabi`, which uses the NVs as the probe.

---

## Fitting: what we use and what we do not

| Use | Not | Why |
|---|---|---|
| `rabiFitRobust` | `rabiFit`, `ourRabiFit` | the old fit's amplitude start is 2× the truth with the truth as its lower bound, so `a` pins at the bound; the frequency starts from a guess. The new one scans frequency exhaustively — no start point to get wrong. |
| `esrFitRobust` → `esrMultiLineFit` | `ESRfit` | N Lorentzians on a **shared** baseline with the centre covariance exposed. Splittings are differences of centres, so var = C₁₁ + C₂₂ − 2C₁₂ and the cross term is not optional; `confint`, which the old path uses, does not expose covariance at all. |
| own extraction | `extract1Ddata`, `extractESRdata` | `extract1Ddata` cuts the ratio across averages at 5, then 3, then **one** standard deviation. A one-sigma cut discards about a third of good Gaussian data, after two harsher passes have shrunk the spread it measures against, and it biases survivors toward the mean — so every error bar downstream is too small. `extractESRdata` is a four-line wrapper around it. |

Both robust fits report R² **next to the ceiling the noise imposes**, `1 − Σe²/Σ(y−ȳ)²`.
R² = 0.82 against a ceiling of 0.84 is as good as the photons allow; the same 0.82 against
0.99 is a real failure. One number cannot distinguish those, which is why an absolute R² gate
was wrong: at the ~1.4% contrast this setup gives, the ceiling is around 0.8, so a 0.9
threshold rejects sound fits for being quiet.

Every fit has an offline grader that runs in seconds with no hardware: `testFieldPipeline`,
`testB1RabiChain`, `testMwPairChain`. Run them before the hardware, not after. They live
*with* the calibrations rather than in a test folder, because that is what keeps them being
run.

## Still open

1. `AWGExpEcho`, `AWGExpT1`, `AWGExpRamsey` need `obj.AWG.workChannels =
   obj.activeChannels;` in `prepare`. `awg.keepOn = true` covers it by leaving the outputs
   on. With neither, they run with the microwave off and return every point at ratio ≈ 1 with
   no error raised.
2. `AWGExpESR` uses `workChannels = 1:obj.nChannels` where `activeChannels` is correct. Only
   matters at `nChannels = 1`, where it drives output 1 whatever `MWChannel` says.
3. `detectionDuration = 3.88` µs in the json against level 2's ~0.3 µs. See level 2.
4. Frequency quantization on the short burst waveform: a whole number of MW cycles but a
   rounded sample count, so each point's true frequency can be off by ~0.5 MHz, varying
   across the sweep — a quarter of a 2 MHz step.
5. `pulse streamer time out` on 1–4 of every 100 runs. On the walk path `resyncEveryPass`
   absorbs it; on the per-point path it is what produces the ratio ≈ 1 outliers.
6. `AWGExpRabi` and `AWGExpRabiCont` share `NAME = 'AWG Rabi'`, so `getExpByName` is
   ambiguous. Construct them directly.
7. Two `extractESRdata` / `extractESRData` files differing only in one letter's case, with
   different contents, both in `prod`. Windows does not distinguish them. Nothing of ours
   calls either, but the collision is still there.
