%% runScriptBasicExpAWG  --  Setup 7
% Run the session cell first, and again after re-running main.
% Walk mode, newOrderEveryAverage and the CW-ESR window are class defaults: no need to set them.
% Fits take one struct and return one struct: o = rabiFitRobust(struct('exp', b)).

%% Start of session
[jsonLocation, h, awg] = setup7Session;     % setup7Session('2026-08-14') pins an older calibration


%% ESR, zero field  --  the single dip at D
h.control('off');

a = AWGExpESR;
a.loadJsonParams('all', jsonLocation);      % also sets the CW window (esrWindowPick)
a.frequency  = 2870 + (-100:1:100);
a.repeats    = 50;
a.averages   = 3;                           % never 1: sterr is not a standard error at 1
a.amplitude  = -16;
a.isTracking = 0;

tId = tic; a.run;
fprintf('ESR: %.1f s, %d SCPI commands\n', toc(tId), a.AWG.commandCount);
o = esrFitRobust(struct('exp', a, 'nLines', 1));    % the resonator's dip is fitted separately


%% ESR, with field  --  the four-orientation spectrum
h.control('auto');

a = AWGExpESR;
a.loadJsonParams('all', jsonLocation);
a.frequency  = 2750 : 2 : 3000;
a.repeats    = 100;
a.averages   = 3;
a.amplitude  = -35;
a.isTracking = 0;

tId = tic; a.run;
fprintf('ESR: %.1f s, %d SCPI commands\n', toc(tId), a.AWG.commandCount);
o = esrFitRobust(struct('exp', a));     % lines, resonator dip, and the field: o.B, o.tiltDeg, o.D


%% Pulsed ESR
h.control('auto');

p = AWGExpESRPulsed;
p.loadJsonParams('all', jsonLocation);
p.frequency       = 2762 + (-10:0.1:10);
p.repeats         = 300;
p.averages        = 4;
p.amplitude       = -44;
p.mwPulseDuration = 0.572;
p.run;
o = esrFitRobust(struct('exp', p));     % also tests for the 14N triplet (o.hyperfine)
%% 
h.control('auto');

a = AWGExpESR;
a.loadJsonParams('all', jsonLocation);
a.frequency       = 2762 + (-10:0.1:10);
a.repeats         = 300;
a.averages        = 4;
a.amplitude       = -44;
a.isTracking = 0;

tId = tic; a.run;
fprintf('ESR: %.1f s, %d SCPI commands\n', toc(tId), a.AWG.commandCount);
o = esrFitRobust(struct('exp', a));     % lines, resonator dip, and the field: o.B, o.tiltDeg, o.D



%% Pulsed vs CW ESR comparison
out = esrPulsedVsCwCompare;


%% Rabi, COARSE  --  find the period
h.control('auto');

b = AWGExpRabi;
b.loadJsonParams('all', jsonLocation);
b.nChannels = 2;                    % BEFORE amplitude: set.nChannels re-enters
b.amplitude = -20 + [0 0];          % set.amplitude, which would add amplitudeDiff twice
b.timeDelay = 0;                    % channelPhase is the phase knob; they are additive
b.repeats          = 300;
b.averages         = 2;
b.isTracking       = 0;
b.balancedSequence = 1;
b.shouldAutosave   = true;
b.tau = round(0.02 : 0.06 : 6, 10);

b.AWG.resetCommandCount; tId = tic;
b.run;
fprintf('coarse Rabi: %.1f s, %d SCPI commands\n', toc(tId), b.AWG.commandCount);

o = rabiFitRobust(struct('exp', b));
fprintf('pi = %.1f +/- %.1f ns at %.3f +/- %.3f MHz\n', ...
        1e3*o.pulses(2), 1e3*o.pulsesErr(2), o.f, o.fErr);


%% Rabi, FINE  --  the measurement, tau range from the coarse fit
t_pi = o.pulses(2);                 % o.pulses is already in us
if ~isfinite(t_pi) || t_pi <= 0
    error(['Fine Rabi: the coarse fit did not return a usable pi time (%g us). Re-run the ' ...
           'coarse cell and check its R^2 against the noise ceiling before continuing.'], t_pi);
end
tauMax = 7 * t_pi;                              % about 3.5 periods
dtau   = max(round((2*t_pi)/10, 3), 0.005);     % 10 points per period, floored at one AWG step

h.control('auto');
b = AWGExpRabi;
b.loadJsonParams('all', jsonLocation);
b.nChannels = 2;
b.amplitude = -16 + [0 0];
b.timeDelay = 0;
b.repeats          = 3000;
b.averages         = 2;
b.isTracking       = 0;
b.balancedSequence = 1;
b.shouldAutosave   = false;
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
b.tau = tauList;

b.AWG.resetCommandCount; tId = tic;
b.run;
fprintf('fine Rabi: %.1f s, %d SCPI commands\n', toc(tId), b.AWG.commandCount);

o = rabiFitRobust(struct('exp', b));
fprintf('pi = %.1f +/- %.1f ns at %.3f +/- %.3f MHz\n', ...
        1e3*o.pulses(2), 1e3*o.pulsesErr(2), o.f, o.fErr);


%% Rabi, sanity check  --  on whichever Rabi you just ran
% R^2 means little alone: compare it with the noise ceiling, the best any fit could reach.
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


%% Ramsey  --  T2*, on resonance
% No detuning needed: the 14N hyperfine splitting gives a 2.158 MHz beat, which is the fringe.
h.control('auto');

r = AWGExpRamsey;
r.loadJsonParams('all', jsonLocation);
r.nChannels = 2;
r.amplitude = -16 + [0 0];
r.timeDelay = 0;
r.halfPiTime = 0.0539/2;
r.frequency  = 2820.6;
r.repeats          = 3000;
r.averages         = 2;
r.balancedSequence = 1;
r.constantTime     = 1;             % else the laser duty cycle varies across the scan: T2* biased short
r.smallDelay       = 0.5;
r.isTracking       = 0;
r.shouldAutosave   = false;
r.laserInitializationDuration = 16; % level 2's SNR-per-time optimum; delete to use the json

hyperfineMHz = 2.158;               % 14N. Use 3.030 and nHyperfine = 2 for a 15N sample.
tauMax       = 2.5;                 % past T2*, or the decay is an extrapolation
perBeat      = 30;
dtau         = round(1/(hyperfineMHz*perBeat), 4);
r.tau        = round(0.005 : dtau : tauMax, 4);
fprintf('Ramsey tau: %d points, %.4f us step, %.2f us span = %.1f beats of %.3f us\n', ...
        numel(r.tau), dtau, max(r.tau), max(r.tau)*hyperfineMHz, 1/hyperfineMHz);

r.AWG.resetCommandCount; tId = tic;
r.run;
fprintf('Ramsey: %.1f s, %d SCPI commands\n', toc(tId), r.AWG.commandCount);

o = ramseyFitRobust(struct('exp', r));
fprintf('T2* = %.3f +/- %.3f us, detuning %.3f +/- %.3f MHz, contrast %.2f %%\n', ...
        o.T2star, o.T2starErr, o.delta, o.deltaErr, 100*o.contrast);
% Variants: add 'showResid', 'showSpectrum', 'useProjection','one', 'nHyperfine',2 or
% 'fitHyperfine',true to the struct. testRamseyFit grades the fitter without hardware.


%% Ramsey, sanity check  --  on the fit you just ran
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


%% Echo  --  T2
% First time: testEchoWalk (~1 s, no hardware) and echoGateTest (~4 min, does the fixed MW
% gate leak?). Walk holds the gate at max(tau), so a leak would look like a shorter T2.
h.control('auto');

e = AWGExpEcho;
e.loadJsonParams('all', jsonLocation);
e.nChannels = 2;
e.amplitude = -16 + [0 0];
e.timeDelay = 0;
e.repeats          = 2000;
e.averages         = 2;
e.balancedSequence = 1;
e.constantTime     = 1;
e.smallDelay       = 1;
e.isTracking       = 0;
e.shouldAutosave   = false;
e.laserInitializationDuration = 16;
e.cycles           = 1;             % 1 = plain Hahn echo
e.pulsePhase       = 'Y';
e.frequency        = 2820.6;
e.halfPiTime       = 0.0539/2;
e.piTime           = 0.0539;

% Linear, not log: a log grid cannot sample the 13C revivals. If the fit reports only one
% revival, double the span.
e.tau = unique(round(linspace(1, 60, 60), 3));
fprintf('Echo tau: %d points, %.1f to %.1f us -> evolution %.1f to %.1f us at cycles = %d\n', ...
        numel(e.tau), min(e.tau), max(e.tau), ...
        2*min(e.tau)*e.cycles, 2*max(e.tau)*e.cycles, e.cycles);

e.AWG.resetCommandCount; tId = tic;
e.run;
fprintf('Echo: %.1f s, %d SCPI commands\n', toc(tId), e.AWG.commandCount);

% p (stretch) is fitted, so compare o.coherenceArea between runs rather than T2.
o = echoFitRobust(struct('exp', e, 'T1us', []));      % T1us fills in the T2 <= 2*T1 check
fprintf('Echo: T2 = %.3f +/- %.3f us, p = %.2f, contrast = %.3f%%\n', ...
        o.T2, o.T2Err, o.p, 100*o.contrast);


%% Echo, XY8  --  T2 under decoupling
% tau is the half-spacing: total evolution = 2*tau*cycles. The MW gate grows 8x with
% cycles; echoGateTest validated up to 240 us of excess, which is where this tau list stops.
e = AWGExpEcho;
e.loadJsonParams('all', jsonLocation);
e.nChannels = 2;
e.amplitude = -16 + [0 0];
e.timeDelay = 0;
e.repeats          = 2000;
e.averages         = 2;
e.balancedSequence = 1;
e.constantTime     = 1;
e.smallDelay       = 1;
e.isTracking       = 0;
e.shouldAutosave   = false;
e.laserInitializationDuration = 16;
e.cycles           = 8;
e.pulsePhase       = 'XYXYYXYX';    % XY8; at cycles = 16 it repeats, giving XY16

e.tau = unique(round(linspace(1, 15, 120), 3));
fprintf('XY8 tau: %d points, %.2f to %.2f us -> evolution %.1f to %.1f us, gate %.1f us\n', ...
        numel(e.tau), min(e.tau), max(e.tau), ...
        2*min(e.tau)*e.cycles, 2*max(e.tau)*e.cycles, ...
        2*max(e.tau)*e.cycles + 2*e.halfPiTime);

e.AWG.resetCommandCount; tId = tic;
e.run;
fprintf('XY8: %.1f s, %d SCPI commands\n', toc(tId), e.AWG.commandCount);

% T2 and the revival period are valid here; the implied 13C field is not (Hahn-echo mapping).
o = echoFitRobust(struct('exp', e));


%% T1  --  ExpT1
% ExpExtendedT1 can also measure relaxation from -1, but needs extra variables defined.
t = ExpT1;
t.loadJsonParams('all', jsonLocation);
t.tau = linspace(1,3000,100);
t.constantTime = 1;
t.run;
%o = t1FitRobust(struct('exp', t));


%% T1  --  AWGExpT1. Doesn't work yet.
% No walk mode: the AWG is already armed once per sweep. mwGateMode stays 'tau' until
% t1GateTest (~25 min) says whether the open gate leaks; a leak would just shorten T1.
t = AWGExpT1;
t.loadJsonParams('all', jsonLocation);
t.nChannels = 2;
t.amplitude = -16 + [0 0];          % the amplitude the pulse times were calibrated at
t.isTracking     = 0;
t.GIMeas         = 0;
t.AWGorSRSswitch = 0;
t.lastDelay      = 1;
t.repeats        = 3000;
t.averages       = 6;               % ~0.4% signal needs ~18000 shots per point
t.shouldAutosave = false;
t.constantTime   = false;           % 4.6x faster; if the |0> control goes past ~3 sigma, set true
t.mwGateMode     = 'tau';           % 'tau' | 'pulse' (closed during the wait) | 'all'
t.tau            = unique(round(logspace(log10(100), log10(6000), 12)));
t.detectionDuration = 0.5;
t.referenceDetectionDuration = 0.5;
t.run;

% Read o.measured, the |0> control and span/T1 before T1 itself.
o = t1FitRobust(struct('exp', t));
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
% Only when something looks wrong (experiments\_diagnostics\ and functions\):
%   esrWalkValidate(a)   esrWalkCacheCheck(a)   rabiWalkCompare(...)   rabiStageTimer(...)
%   awgDrainErrors(awg)  awgTransportBench(awg)  scopeProbe   esrWindowPick
% Are the MW switches gating? A flat spectrum here means yes:
%   a.driveMWSwitches = false; a.restart; a.driveMWSwitches = true;
