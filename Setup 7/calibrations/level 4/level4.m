%% ========================================================================
%  LEVEL 4  --  the microwave pair, and the pulse times everything above runs on
%
%  Setup 7 drives the antenna from two AWG channels, each through its own amplifier and its
%  own circulator, combined at the sample. Two numbers describe the pair:
%
%    amplitudeDiff  [0 delta]  dB to add to channel 2 so both deliver the same B1
%    channelPhase   [0 phi]    degrees on channel 2 for maximum combined drive
%
%  and once the pair is balanced and in phase, one more measurement turns the drive into
%  pulse times:
%
%    halfPiTime, piTime, threeHalvesPiTime   us, from a Rabi at the operating amplitude
%
%  4.4 IS THE END OF THE CHAIN. Everything below it exists so that piTime means something:
%  level 1 says when the pulses actually arrive, level 2 says when to count photons, level 3
%  says what frequency and field to sit at, and 4.2 says how to combine the two channels.
%  Ramsey, Echo and T1 then read halfPiTime and piTime out of the json and measure nothing
%  themselves -- so if those two numbers are wrong, every measurement above this level is
%  wrong in a way that looks like physics rather than like a mistake.
%
%  calibMwPair measures the pair in one run, in four stages. Run level 3 first and leave the
%  coils where calibB1Rabi put them, so this is measured at the antenna's response peak
%  where the drive is strongest.
%
%  WHY ONE CHANNEL AT A TIME IS VALID, since it is the obvious thing to doubt. Everything
%  between the AWG outputs and the diamond is a linear network, so the transfer function
%  from output n to B1 is a fixed property of the hardware and cannot depend on what the
%  other port is doing. The circulators make that stronger, not weaker: they are exactly
%  what stops channel 2's reflected power reaching channel 1's amplifier output. And since
%  each amplifier only sees its own channel, turning the other one on does not move its
%  operating point.
%
%  Stage 3 measures this rather than assuming it. In a linear network at the optimum phase
%  Omega_both = Omega1 + Omega2 exactly, so R = Omega_both/(Omega1+Omega2) is a direct test
%  of superposition; R = 1 means the one-at-a-time calibration is sound.
%
%  EACH CHANNEL AT THE LEVEL IT PLAYS (2026-11-23). Stage 1 now measures channel 2 at
%  operating + delta and iterates until it matches channel 1 there, instead of measuring
%  both at the same dBm and extrapolating. At -16 dBm the old way left the pair ~1.7 dB
%  unbalanced (closure R = 0.91: channel 2 compressing). And 4.2b then measures the balance
%  with BOTH channels playing (calibMwQuadrature), which no single-channel run can.
%
%  THE ORDER, and why there is no amplitude <-> phase loop: the phase optimum does not
%  depend on the balance (Omega^2 = O1^2 + O2^2 + 2 O1 O2 cos(phi - phi0) peaks at phi0 for
%  any O1, O2), stage 1 uses no phase at all, and 4.2b cancels a phase error to first order.
%  The one coupling is a compressing amplifier's level-dependent phase, which is why the
%  amplitude comes FIRST and the phase is scanned at the final level. A second round is
%  needed only when 4.2b disagrees with 4.2 -- and then one is enough.
%
%  THE REAL CAVEAT IS LEVEL, NOT BOTH-ON. Per-channel amplifiers compress at different
%  points, so delta measured at -30 dBm need not hold at -16. Calibrate at the amplitude
%  you run at, and re-run stage 1 alone if you change it. The PHASE does not need redoing --
%  it is set by path lengths, not by level.
%
%  ONE AMPLITUDE THROUGHOUT THIS LEVEL. Set it once in 4.0 and let 4.2 and 4.4 both take it
%  from there. amplitudeDiff, channelPhase and the three pulse times are all properties of the
%  drive AT A STATED LEVEL; measured at different levels they are individually correct and
%  jointly meaningless.
%
%  Writes to expProperties.json, through 4.3 and 4.5 and never as a side effect:
%    amplitudeDiff, channelPhase, timeDelay = 0                     from 4.2, checked by 4.2b
%    amplitude, frequency, halfPiTime, piTime, threeHalvesPiTime    from 4.4
%  ========================================================================

%% 4.0  the operating amplitude, chosen once
%
% -16 dBm per channel is the Rabi standard on this setup and calibMwPair's own default. Every
% number in this level is measured here and is valid at no other level, so if you change it,
% re-run 4.2 (stage 1 is enough -- see 4.7) AND 4.4. Neither one warns you.
opAmplitude  = -16;
jsonLocation = setup7Paths().json;

% A new random tau order before every average in the walk-mode Rabi runs (4.2, the coarse
% pass of 4.4, 4.7). Without it the order is drawn once per waveform load, so anything tied
% to position in the pass sits on the same taus for the whole run. Costs one sequence upload
% + arm per average. false = the old behaviour.
newOrder = true;

% 4.3 writes the pair only when 4.2b (both channels on) agrees with 4.2. Set false to write
% on 4.2 alone -- e.g. when 4.2b came back 'inconclusive' and you accept that.
requireQuadrature = true;


%% 4.1  offline grader  --  no hardware, ~2 s
%
% This one is mostly about SIGNS. The physics is easy; what goes wrong is whether delta
% comes out positive when channel 2 is the weak one, whether dB is 20*log10 or 10*log10 of
% the thing being compared, and whether the fitted phase is the maximum or the null. None of
% those fail on the bench -- they show up as a setup that ran for a month with the two
% channels 180 degrees apart. Expect ALL TESTS PASS.
r = testMwPairChain;

% The both-on test of 4.2b, graded against a simulated NV ensemble (14N triplet,
% inhomogeneous width, opposite field gradients from the two PCBs) with a known imbalance,
% a phase error and pulse-time errors. It plays AWGExpMwQuadrature's own pulse definitions,
% so it also catches a change there. Expect ALL TESTS PASS.
rq = testMwQuadrature;


%% 4.2  the pair  --  amplitudeDiff and channelPhase
%
% Coils ON, at the field level 3 chose. About 19 Rabi runs, ~30 minutes: every point
% reloads the AWG because both amplitude and phase change the waveform samples. It never
% touches the Helmholtz.
%
% amplitude comes from 4.0 rather than from calibMwPair's default, which is the same -16 dBm
% today. Passed explicitly so that changing 4.0 actually changes what is measured here, rather
% than leaving 4.2 at the default while 4.4 moves.
out = calibMwPair(struct('amplitude', opAmplitude, 'phases', 0:15:315, ...
                         'newOrderEveryAverage', newOrder, 'tau', linspace(0.005, 1, 100), 'minContrast', 0));


%% 4.2b  the balance with BOTH channels on  --  calibMwQuadrature
%
% Same field, same amplitude. About 10 minutes at the defaults (14 configurations, 10000
% passes x 20 averages). Writes nothing.
%
% A pulse with the two channels 90 degrees apart has its rotation axis at atan(Omega2/Omega1):
% 45 degrees exactly when both deliver the same field at the NVs WHILE BOTH PLAY. An in-phase
% pi/2 about 135 / -45 reads the tilt, with channel 2 leading and with channel 1 leading, at
% three channel-2 levels around 4.2's delta. Where that line crosses zero is the both-on
% balance. It cancels a channelPhase error to first order, so it does not need 4.2's phase to
% be perfect, and it does not use the antiphase null, which the hyperfine washes out.
%
% It uses the pair's pi time from 4.2's closure run (both channels, at the maximum phase),
% since the json's may belong to another amplitude.
piPair = [];
if isfield(out, 'closure') && isstruct(out.closure) && isfield(out.closure, 'piNs') ...
        && isfinite(out.closure.piNs)
    piPair = out.closure.piNs / 1e3;
end
q = calibMwQuadrature(struct('amplitude', opAmplitude, 'delta', out.delta, ...
                             'deltaErr', out.deltaErr, 'phiMax', out.phiMax, ...
                             'piTime', piPair, 'newOrderEveryAverage', newOrder));


%% 4.3  write the pair
%
% timeDelay goes to zero in the SAME call, which is why it is in the same struct rather than
% left to be done by hand afterwards. It is the other phase mechanism and the two are ADDITIVE
% inside the waveform builder, so with both non-zero the real phase is their sum and nothing
% records that. The stored 327 ps is within 6% of a full period at 2870 MHz anyway, so
% dropping it costs nothing.
%
% expPropsWrite backs the file up, edits only these three values in place, then re-reads it
% and checks that nothing else moved. Add 'dryRun', true to see the diff without writing.
%
% Writes only when 4.2 settled (out.ok) AND 4.2b agreed (unless requireQuadrature = false):
% calibMwPair's hand-over says nothing should be written until a failing stage is understood,
% and a balance that does not hold with both channels on is a failing stage.
quadOk = exist('q', 'var') && isstruct(q) && isfield(q, 'verdict') && strcmp(q.verdict, 'agree');
if out.ok && (quadOk || ~requireQuadrature)
    w = expPropsWrite(struct('amplitudeDiff', [0 out.delta], ...
                             'channelPhase',  [0 out.phiMax], ...
                             'timeDelay',     0), ...
                      struct('jsonLocation', jsonLocation));
    if ~quadOk
        fprintf('\nWritten on 4.2 alone (requireQuadrature = false); 4.2b did not agree.\n');
    end
elseif ~out.ok
    fprintf(['\n4.2 did NOT settle (see its hand-over block), so the pair has not been written.\n' ...
             'Fix the failing stage and re-run 4.2 first.\n']);
elseif exist('q', 'var') && isstruct(q) && strcmp(q.verdict, 'disagree')
    % One more round at the both-on balance: the phase scan at q.deltaBoth (stage 1 skipped,
    % the order rule: amplitude first, then phase), then 4.2b again at the new pair.
    fprintf(['\n4.2b DISAGREES: with both channels on the balance is %+.3f dB, not %+.3f.\n' ...
             'Nothing written. One more round -- run these two, then this cell again:\n\n' ...
             '  out = calibMwPair(struct(''amplitude'', opAmplitude, ''deltaStart'', %.4f, ' ...
             '''phases'', 0:15:315, ...\n' ...
             '                           ''newOrderEveryAverage'', newOrder, ''tau'', ' ...
             'linspace(0.005, 1, 100), ''minContrast'', 0));\n' ...
             '  q = calibMwQuadrature(struct(''amplitude'', opAmplitude, ''delta'', out.delta, ' ...
             '''phiMax'', out.phiMax, ...\n' ...
             '                               ''piTime'', out.closure.piNs/1e3, ' ...
             '''newOrderEveryAverage'', newOrder));\n'], q.deltaBoth, out.delta, q.deltaBoth);
else
    fprintf(['\n4.2b was inconclusive or not run, so the pair has not been written. Re-run 4.2b ' ...
             'with more averages\n(struct field ''averages''), or set requireQuadrature = false ' ...
             'in 4.0 to write on 4.2 alone.\n']);
end

if isfield(out, 'closureRatio') && isfinite(out.closureRatio)
    fprintf('\nclosure R = %.4f +/- %.4f  (1 means the two channels add linearly)\n', ...
            out.closureRatio, out.closureRatioErr);
end


%% 4.4  the pulse times  --  measures and prints, writes nothing
%
% A Rabi at the operating amplitude with the pair already balanced and in phase, from which
% pi/2, pi and 3pi/2 follow. Coils where level 3 left them.
%
% Two passes, for two different reasons:
%
%   COARSE, walk mode, ~20 s. Only has to find the period, so it can be cheap -- walk re-arms
%   about a dozen times instead of once per tau. Its possible contrast bias does not matter
%   when all that is wanted is a frequency to size the next pass with.
%
%   FINE, classic mode, high repeats. This is the measurement, so it uses the transport that
%   has been validated for data. repeats is the cheap axis: 4000 repeats costs photons (~55 us
%   each) and no re-arms, while every extra tau point costs a 0.6 s re-arm per average. So
%   repeats are high, the spacing is dense (13 points per period), and the tau list is only as
%   long as the fit needs -- six periods. Expect 4-6 minutes.
%
% THEN IT CHECKS ITSELF, which is the part worth waiting for. The fine curve is fitted three
% times: whole, first half, second half. If the halves disagree by more than 3 sigma on pi,
% either the drive drifted during the run (amplifier warm-up, usually) or the rotation angle
% is not linear in the pulse duration -- and the second of those is exactly what makes the
% derived pi/2 wrong while pi still looks right. A single fit cannot see either one; it
% reports a small error bar regardless, because the model has no way to say "the frequency
% changed halfway through".
%
% Read the WHAT TO DO block. The one thing in it worth knowing in advance: pi/2 is NOT an
% independent measurement. A Rabi fit yields one frequency and one phase offset, so
% t_k = k/(4f) + x0 -- the three pulse times are equally spaced by a quarter period BY
% CONSTRUCTION, and pi/2's error is degenerate with x0. Treat halfPiTime as provisional until
% a Ramsey confirms it.
p = calibPulseTimes(struct('amplitude', opAmplitude, 'jsonLocation', jsonLocation, ...
                           'newOrderEveryAverage', newOrder));


%% 4.5  write the pulse times
%
% A separate cell on purpose: 4.4 is free to re-run and this is not, because it moves the
% numbers Ramsey, Echo and T1 will use. Read 4.4's verdict before running it.
%
% threeHalvesPiTime is probably NOT in the file yet, and expPropsWrite refuses unknown keys by
% default -- deliberately, since a typo'd key would otherwise be created silently and then
% ignored forever by loadStateFromStruct. 'create' true is right here and only here.
%
% To see the diff first:
%   expPropsWrite(p.toWrite, struct('jsonLocation', jsonLocation, 'dryRun', true, ...
%                                   'create', true));
if ~isequal(p.stable, false)
    v = expPropsWrite(p.toWrite, struct('jsonLocation', jsonLocation, 'create', true));
else
    fprintf(['\n4.4 reported the drive was NOT stable across the run, so nothing has been\n' ...
             'written. Leave the microwave chain powered for ten minutes and re-run 4.4.\n' ...
             '\nIf the halves still disagree after that, it is not warm-up, and the pulse\n' ...
             'times are limited by that spread however small the fit''s own error bar is.\n' ...
             'Write them anyway -- call expPropsWrite by hand, as in the commented line\n' ...
             'above without dryRun -- but treat the spread as the real uncertainty, and\n' ...
             'expect Echo and Ramsey contrast below what the Rabi contrast suggests.\n']);
end


%% 4.6  check before moving up
fprintf('\n==================== LEVEL 4 RESULT ====================\n');
fprintf('measured at %+.2f dBm per channel -- valid there and at no other level\n\n', ...
        opAmplitude);
fprintf('   amplitudeDiff     [0, %+.3f] dB\n', out.delta);
if exist('q', 'var') && isstruct(q) && isfield(q, 'verdict')
    fprintf('   both channels on  [0, %+.3f] +/- %.3f dB   (4.2b: %s)\n', q.deltaBoth, ...
            q.deltaBothErr, q.verdict);
end
fprintf('   channelPhase      [0, %.2f] deg\n', out.phiMax);
fprintf('   timeDelay         0\n');
if isfield(p, 'piTime') && isfinite(p.piTime)
    fprintf('\n   Rabi f            %.4f +/- %.4f MHz\n', p.fRabi, p.fRabiErr);
    fprintf('   halfPiTime        %8.2f ns   -> Ramsey, Echo\n', 1e3*p.halfPiTime);
    fprintf('   piTime            %8.2f ns   -> Echo, T1, ExperimentAWG, digitizer calibs\n', ...
            1e3*p.piTime);
    fprintf('   threeHalvesPiTime %8.2f ns   -> nothing reads it; kept for the record\n', ...
            1e3*p.threeHalvesPiTime);
    fprintf('   contrast          %8.2f %%\n', 100*p.contrast);

    % ---- the one check that tests the whole chain at once ----
    % Cheap, and it catches what no single step can: a wrong operating amplitude, a frequency
    % off resonance, and a write that landed in a different json than these cells read. If a
    % pi pulse really is a pi pulse, a Rabi swept narrowly around piTime has its contrast
    % minimum AT piTime. Anywhere else means one of those three.
    fprintf(['\nCONFIRM ON HARDWARE, ~1 minute: one Rabi over a narrow window around\n' ...
             '   tau = %.4f us,   e.g.   b.tau = %.4f : 0.001 : %.4f\n' ...
             'and check the minimum sits AT piTime rather than beside it. That single run\n' ...
             'catches a wrong amplitude, a frequency off resonance, and a write that went to\n' ...
             'a different json -- none of which any cell above can see.\n'], ...
            p.piTime, max(p.piTime - 0.008, 0.002), p.piTime + 0.008);
end
fprintf('========================================================\n');


%% 4.7  re-run stage 1 only, after changing the operating amplitude
%
% The balance is level dependent; the phase is not. This measures delta at a new amplitude
% without spending twenty minutes re-scanning the phase.
%
% THEN RE-RUN 4.4 AND 4.5. The pulse times moved too, and the only record of which amplitude
% they belong to is the `amplitude` field 4.5 writes beside them.
%
% Stage 1 iterates to the level channel 2 will play (iterateBalance), so this is right in
% compression too. Run 4.2b at the new level before trusting it.
%
% opAmplitude = -10;
% out2 = calibMwPair(struct('amplitude', opAmplitude, 'doPhase', false, ...
%                           'doClosure', false, 'doAntiphase', false, ...
%                           'newOrderEveryAverage', newOrder));
% expPropsWrite(struct('amplitudeDiff', [0 out2.delta]), ...
%               struct('jsonLocation', jsonLocation));
% p = calibPulseTimes(struct('amplitude', opAmplitude, 'jsonLocation', jsonLocation, ...
%                            'newOrderEveryAverage', newOrder));
% expPropsWrite(p.toWrite, struct('jsonLocation', jsonLocation, 'create', true));


%% 4.8  DIAGNOSTIC, not a calibration  --  antenna reflection
%
% calibResonator measures power coming BACK out of the circulators' reverse ports, i.e.
% power failing to get into the antenna. That is not B1 at the NVs, and the two differ by
% everything between the antenna terminals and the diamond -- coupling, standing waves, and
% losses that dissipate without radiating. A reflection minimum can sit several MHz from the
% field maximum.
%
% The operating frequency comes from calibB1Rabi (level 3), which uses the NVs as the probe.
% This is here to diagnose the antenna, not to set a number.
%
% freqs = 1e9*linspace(2.7, 3, 201);
% r = calibResonator(struct('freqs', freqs, 'nSweeps', 1, 'vpp', 0.01, 'settle', 1));
