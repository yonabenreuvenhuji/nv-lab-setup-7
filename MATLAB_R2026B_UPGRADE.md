# Setup 7 on MATLAB R2026b

Branch `setup-7-r2026b` (off `setup-7`). Goal: everything setup 7 runs works on R2026b. The lab
PC currently runs R2023b.

**Status:** code changes done and tested on R2025b (the newest release on the dev PC), with and
without Java, against local stand-ins for the hardware. **Not yet run on R2026b, on the hardware,
or with the GUI** (the dev PC has no GUI Layout Toolbox). Run `matlabUpgradeCheck` on the lab
PC (below) before the first hardware test.

---

## 1. What breaks between R2023b and R2026b, and what changed

| Problem | Since | Fix | On R2023b |
|---|---|---|---|
| `visa()` and its `fopen`/`fscanf`/`query`… **throw an error**. The Tektronix AWG driver used `visa()` by default; the constructor connects before any script can switch. | R2026b | `FrequencyGeneratorTektronixAWG.useVisadev` now defaults to `true`. The `visadev` path already existed, and `runScriptBasicExpAWG` already switched to it. | Same as running the run script today |
| `tcpip()` throws an error. Only affects setup 7's SRS configuration (`setupInfoSRS.json`). | R2026b | `FrequencyGeneratorSRS` uses `tcpclient`, still connecting on first use and closing in `disconnect`. | Byte-identical commands (A/B test) |
| **No Java runtime is installed.** `urlread2` (every Pulse Streamer command) and `DataHash` (sequence caching) need Java. | R2026b | Java paths are kept when a JRE is present. Without one: `PulseStreamerHttpNoJava` (HTTP over `tcpclient`) and `DataHashNoJavaEngine` (`keyHash` digest). | Unchanged (Java path) |
| `urlread2` called three undocumented `com.mathworks.*` Java classes from the retired Java desktop. | removal announced | Replaced with core `java.nio` and Java's default proxy handling. | Identical replies, slightly faster |
| The Pulse Streamer base64 encoder used `org.apache.commons` (bundled JAR, being removed). Its fallback only triggered if the whole JVM was missing. | removal announced | Vectorized pure-MATLAB encoder. | Byte-identical output; 1M pulses ~90 ms vs ~65 ms |
| `JavaFrame` was removed (`GuiController.maximize`; nothing calls it). | R2025a | `figure.WindowState = 'maximized'`. | — |
| HTML labels in `uicontrol` show as raw text: the "Open In Figure" button of every experiment window showed `<html><center…`. | R2025a | Two-line plain label on R2025a+, HTML kept on older releases. | Unchanged |

New files:
- `src/Physical/PulseGenerator/PulseStreamer/PulseStreamerHttpNoJava.m`
- `src/Util/DataHashNoJavaEngine.m`
- `Setup 7/functions/matlabUpgradeCheck.m`

### How the no-Java Pulse Streamer transport behaves
- It sends the same request line and headers as Java's `HttpURLConnection`, and keeps the
  connection open between calls. Opening a `tcpclient` costs ~60–100 ms; a request ~2 ms.
- It closes an idle connection after 2 s (Java: ~5 s). On R2026b+ it also checks
  `tcpclient.Connected` before reusing one.
- If a *reused* connection fails before any reply byte arrives, it retries once on a fresh
  connection. It never re-sends after reply bytes have arrived, so a command can't run twice.
- Handles `Content-Length`, chunked and close-delimited replies. A non-2xx status is returned,
  not thrown, as with `urlread2`.

### DataHash without Java
- The digest is two 64-bit `keyHash` values over the same byte stream DataHash always built,
  instead of MD5. It is still **content-based**: a `Sequence` edited in place gets a new key.
- `keyHash` values are only stable within one MATLAB session. Every caller uses them as
  in-memory cache keys (the pulse generator's delay-fixed sequences, the AWG "already loaded"
  fingerprints), so that is enough.

### Pre-existing bugs fixed on the way (not R2026b-specific)
- `AWGsensingACfield`: `properties (Depentent)` was a typo that stopped the class from loading in
  any release. It is now `Dependent`, which matches its `get.signalFrequency` method.
- `SpcmDummy.prepareExperimentCount` lacked the `nParameters` argument every real detector class
  takes since 2021, so no experiment could run in dummy mode.

---

## 2. What was verified (dev PC, R2025b, no hardware)

- **Code Compatibility Analyzer** over all 422 files setup 7 runs (the dependency closure of
  `main.m`, the AWG experiments, the digitizer calibrations, `ExpT1` and all of `Setup 7\`).
  Every remaining finding is in code setup 7 does not run.
- `checkcode`: no new messages in any edited file; zero in the new files.
- **Tektronix AWG:** every legacy-only access is gated by `useVisadev`; the `visadev` branch is
  complete.
- **SRS:** old (`tcpip`) and new (`tcpclient`) versions sent byte-identical command streams to a
  local stand-in SRS and parsed the replies identically.
- **Encoder:** byte-identical to the Java encoder for 0 to 123,457 pulses, all padding cases.
- **HTTP transports:**
  - Patched `urlread2` and `PulseStreamerHttpNoJava` return the same replies and status codes as
    the original `urlread2`. Tested on keep-alive, chunked, close-delimited and HTTP 500 replies,
    a 2 MB upload, a server that drops connections (retry), and a non-HTTP server (clean timeout).
  - Per call: 1.8 ms `urlread2`, 2.4 ms no-Java transport.
- **Smoke test**, with and **without Java** (`matlab -nojvm`, standing in for R2026b):
  - `Setup.init` with setup 7's devices in dummy mode, except the pulse generator, which is the
    real `PulseStreamerNewClass` talking to a local stand-in Pulse Streamer: all objects created.
  - Sequence upload and run: 13 device calls, the same with and without Java.
  - `DataHash` checks pass.
  - An experiment runs up to opening its plot window, which needs the GUI Layout Toolbox.
- **Rendering:** the four "Open In Figure" label variants were drawn on R2025b to pick the fix.

### Not verified here
- **R2026b itself:** not installed on the dev PC. `matlabUpgradeCheck` covers this on the lab PC.
- **The GUI:** needs the GUI Layout Toolbox. R2025a changed `uicontrol` fonts and sizes, so check
  for clipped labels.
- **The hardware:**
  - AWG over `visadev`.
  - NI-DAQmx: `daq.ni.NIDAQmx` is an **undocumented** class of the Data Acquisition Toolbox's
    NI-DAQmx support package; `matlabUpgradeCheck` verifies every member the wrappers use.
  - The CompuScope MEX, the Standa stages, and the real Pulse Streamer server.

---

## 3. Lab PC: steps to switch to R2026b

1. **Install R2026b** with these toolboxes. Each one is used by setup 7's own code:

   | Toolbox | Used for |
   |---|---|
   | Instrument Control | `visadev` (Tektronix AWG); `calibResonator`, `scopeProbe` |
   | Data Acquisition **+ its NI-DAQmx support package** | `NiDaq` (`daq.ni.NIDAQmx`) |
   | Statistics and Machine Learning | `range()` in the AWG driver; `esrLineFit`, `fitStd` |
   | Curve Fitting | `fit()` in `ESRfit`, `rabiFit`, `calibLaserSaturation`, `esrLineFit` |
   | Signal Processing | `findpeaks` in `rabiFit` |
   | Mapping | `extractfield` in `PathHelper`, which builds the GUI's experiment list at startup |
   | Phased Array System | `rotx`/`roty`/`rotz` in `Helmholtz` |
   | Image Processing | `getrect` in `ImageScanResult` (selecting a scan region) |

2. **Add-ons and drivers:**
   - GUI Layout Toolbox **2.4 or newer** (2.4 is the first for R2025a+ graphics).
   - VISA: NI-VISA or Keysight IO Libraries. TekVISA and Keysight VISA can't both be installed.
   - A C compiler for `LoadNIDAQmx`'s `loadlibrary`: MinGW-w64 or Visual C++, then `mex -setup C`.
   - The NI-DAQmx driver, with `nidaqmx.h` at the path in `LoadNIDAQmx.m`.
   - The CompuScope driver.
   - `drivers\ximc-2.10.5` in the repo (it is git-ignored, so copy it into each clone).
3. **Optional: Java.** Install "MATLAB Support for OpenJDK" and run `jenv -clear`. The Pulse
   Streamer and `DataHash` then use exactly today's Java paths, which is the lowest-risk first
   switch-over. Without it, the no-Java paths are used.
4. `git checkout setup-7-r2026b`. Put `src` on the path and run `setup7AddPath`. The `C:\lab`
   pointer files stay as they are.
5. **Run `matlabUpgradeCheck`** and fix every FAIL. Its "expected" entries are code setup 7 does
   not run.
6. **Hardware test, in order:**
   - Start the GUI (`main`) and check the panels.
   - Stages connect; laser AOM; Helmholtz.
   - Calibration levels 0 and 1.
   - A short ESR and Rabi with `runScriptBasicExpAWG`.
   - Compare with recent R2023b results; check autosaves land in `C:\lab\Saves\Setup 7\`.
7. **Rollback:** `git checkout setup-7` and use R2023b, which is untouched.

---

## 4. Known leftovers (not fixed, not on setup 7's path)

- `ExpESR_CS_scan` and `ExpFluorescenceMagnetScan` don't load (pre-existing class-definition
  errors). They only cause warnings when the GUI builds its experiment list.
- Other setups' code that will error on R2026b if used:
  - `SerialControlled` (`serial`/`instrfind`): Windfreak, laser sources, Arduino.
  - `FrequencyGeneratorRigolAWG` (`visa`).
  - The stage tilt panel (`javacomponent`), shown only for tiltable stages, not the Standa.
- Two different `JsonInfoReader.m` copies (`src\Helper\` and `src\Setup\`). Which one runs depends
  on path order; the lab PC resolves the `Setup\` one, whose `setParams` takes one argument.
- Duplicate `P.m`/`PH.m` in `PulseStreamer\` and `PulseStreamer\old\`.
- `datestr`/`now`/`clock`/`nanmean` are "not recommended", but MathWorks has no removal plans.
