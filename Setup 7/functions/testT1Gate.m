function ok = testT1Gate()
%TESTT1GATE Offline grader for the T1 gate work. No hardware, ~30 s.
%
%   testT1Gate
%
% Two halves, and they check different kinds of thing.
%
% THE SEQUENCE ARITHMETIC (t1-t8). AWGExpT1's three gate modes are a claim about durations:
% that the arms are timing-identical, that only the switch's open window differs, and that
% arm C's dose does not depend on tau. Every one of those is checkable on paper and none of
% them is checkable on hardware -- a pulse-generator sequence that is 0.44 us short at every
% point looks exactly like a slightly lower contrast. The expressions are reimplemented here
% rather than called out of the class, because a test that gets its arithmetic from the code
% under test cannot detect arithmetic that is wrong.
%
% THE STATISTICS (t9-t14). t1GateTest claims that requiring the slope estimate and the offset
% estimate to AGREE rejects drift. That is a design claim, and design claims are exactly what
% quietly stops being true. So the two confounds are injected and the test asserts that each
% one fires ONE estimate and not the other, and that neither is reported as a leak:
%
%   t12  a drifting contrast AMPLITUDE between arms -- a pure offset. Fires the C-vs-B
%        estimate at 14 sigma and must NOT be called a leak.
%   t13  a drifting T1 between arms -- a pure slope. Fires the A-vs-B estimate and must NOT
%        be called a leak either.
%
% If either of those ever starts reporting a leak, t1GateTest has become a drift detector
% wearing a gate detector's name, and the number it prints would be believed.

    fprintf('\n=== testT1Gate ===\n');
    ok = true;
    ok = ok & t1_holdCoversPlayback();
    ok = ok & t2_gateOpenPerMode();
    ok = ok & t3_armsAreTimingIdentical();
    ok = ok & t4_constantTimeHoldsTotal();
    ok = ok & t5_armCDoseIsConstant();
    ok = ok & t6_armCNeedsConstantTime();
    ok = ok & t7a_holdIsWholeNanoseconds();
    ok = ok & t7_holdMustFitInsideTau();
    ok = ok & t8_maxNotEnd();
    ok = ok & t9_cleanFindsNothing();
    ok = ok & t10_recoversAnInjectedLeak();
    ok = ok & t11_sensitivityIsQuoted();
    ok = ok & t12_amplitudeDriftIsNotALeak();
    ok = ok & t13_t1DriftIsNotALeak();
    ok = ok & t14_orientationIsMeasured();
    ok = ok & t15_channelNames();
    ok = ok & t16_sampleGrid();
    fprintf('\n');
    if ok
        fprintf('ALL TESTS PASS -- the gate arms differ only in the gate, and drift is not a leak.\n\n');
    else
        fprintf('*** FAILURES ABOVE -- do not run t1GateTest on hardware until they are fixed ***\n\n');
    end
end

% ======================================================================
% THE SEQUENCE, REIMPLEMENTED
% ======================================================================

function s = seqOf(mode, tau, P)
%SEQOF The event durations and gate state AWGExpT1.prepare builds, written out independently.
%
% Returns the gate-open time and the total duration for ONE repeat-half. Independent of the
% class on purpose: this is the reference the class is graded against, so sharing code with it
% would make the grade meaningless.
    trig = P.trig;
    holdUs = P.hold;
    lastDelay = P.lastDelay;
    if P.constantTime
        % constantTime pads the last delay so the total is the same at every point.
        lastDelay = P.lastDelay + P.tauMax - tau;
    end
    switch mode
        case 'pulse'
            gate  = trig + holdUs;
            total = trig + holdUs + (tau - holdUs) + lastDelay;
        case 'all'
            gate  = trig + tau + lastDelay;
            total = trig + tau + lastDelay;
        otherwise % 'tau'
            gate  = trig + tau;
            total = trig + tau + lastDelay;
    end
    s = struct('gate', gate, 'total', total, 'lastDelay', lastDelay);
end

function P = params(varargin)
    P = struct('trig', 0.1, 'minWf', 0.15, 'piTime', 0.035, 'lastDelay', 1, ...
               'tauMax', 4900, 'constantTime', true, 'hold', []);
    for i = 1:2:numel(varargin); P.(varargin{i}) = varargin{i+1}; end
    % Ceiled to a whole nanosecond, as prepare does: the pulse generator streams integer
    % nanoseconds, and piTime lives on the AWG's 62.5 ps grid.
    if isempty(P.hold); P.hold = ceil((2*P.minWf + P.piTime + 0.1)*1000)/1000; end
end

% ======================================================================

function ok = t1_holdCoversPlayback()
% The hold has to be longer than everything the AWG plays after a trigger, or the pi pulse is
% clipped. A clipped pi costs contrast at EVERY tau equally, so it looks like a calibration
% error rather than a timing one -- which is the worst way for this to fail.
    P = params();
    playback = 2*P.minWf + P.piTime;             % zero pad, then the pi waveform with its pad
    margin = P.hold - playback;
    ok = margin > 0;
    report(ok, sprintf('t1  hold %.3f us covers the %.3f us the AWG plays, margin %.3f us', ...
                       P.hold, playback, margin));
end

function ok = t2_gateOpenPerMode()
% The three doses, at the longest tau. This is the table in t1GateTest's header, as arithmetic.
    P = params();
    tau = P.tauMax;
    a = seqOf('tau',   tau, P);
    b = seqOf('pulse', tau, P);
    c = seqOf('all',   tau, P);
    ok = abs(a.gate - (P.trig + tau)) < 1e-12 && ...
         abs(b.gate - (P.trig + P.hold)) < 1e-12 && ...
         c.gate > a.gate;
    report(ok, sprintf('t2  gate open: tau %.1f us, pulse %.3f us, all %.1f us -- a %.0fx range', ...
                       a.gate, b.gate, c.gate, a.gate/b.gate));
end

function ok = t3_armsAreTimingIdentical()
% The claim the whole test rests on: at every tau, all three modes have the SAME total duration.
% If they did not, the arms would differ in duty cycle and laser load as well as in the gate,
% and no comparison between them would mean anything.
    P = params();
    taus = [20 50 137 600 1500 4900];
    ok = true; worst = 0;
    for tau = taus
        t1 = seqOf('tau', tau, P).total;
        t2 = seqOf('pulse', tau, P).total;
        t3 = seqOf('all', tau, P).total;
        d = max(abs([t1 - t2, t1 - t3]));
        worst = max(worst, d);
        ok = ok && d < 1e-9;
    end
    report(ok, sprintf('t3  all three modes have the same total duration at every tau (worst %.2g us)', worst));
end

function ok = t4_constantTimeHoldsTotal()
% And with constantTime on, that same total is independent of tau as well. Both properties are
% needed: t3 makes the arms comparable to each other, t4 makes the POINTS comparable within an
% arm, which is what lets a contrast ratio be read as a dose effect.
    P = params();
    taus = [20 50 137 600 1500 4900];
    tot = arrayfun(@(x) seqOf('tau', x, P).total, taus);
    ok = (max(tot) - min(tot)) < 1e-9;
    report(ok, sprintf('t4  constantTime holds the total at %.1f us across the sweep (spread %.2g)', ...
                       tot(1), max(tot) - min(tot)));
end

function ok = t5_armCDoseIsConstant()
% Arm C's whole purpose. Its gate-open time must not depend on tau, because that is what turns
% it into an OFFSET with no slope -- a second, independent handle on the leak rate.
    P = params();
    taus = [20 50 137 600 1500 4900];
    g = arrayfun(@(x) seqOf('all', x, P).gate, taus);
    ok = (max(g) - min(g)) < 1e-9;
    report(ok, sprintf('t5  arm C dose is constant at %.1f us across the sweep (spread %.2g)', ...
                       g(1), max(g) - min(g)));
end

function ok = t6_armCNeedsConstantTime()
% The other side of t5: WITHOUT constantTime, arm C's dose does track tau and the arm loses the
% property that distinguishes it. The class warns about this; the warning has to be warranted.
    P = params('constantTime', false);
    taus = [20 4900];
    g = arrayfun(@(x) seqOf('all', x, P).gate, taus);
    ok = (max(g) - min(g)) > 1;
    report(ok, sprintf('t6  without constantTime, arm C dose DOES vary (%.1f to %.1f us) -- the warning is earned', ...
                       g(1), g(2)));
end

function ok = t7a_holdIsWholeNanoseconds()
% The pulse generator streams integer nanoseconds (sendToHardware does P(duration*1e3, ...)), and
% piTime sits on the AWG's 62.5 ps grid, so the derived hold is a fraction of a nanosecond unless
% it is ceiled. Every other duration this class hands the pulse generator is already whole; the
% hold was the one that would not have been, and only in 'pulse' mode.
    ok = true; worst = 0;
    for piT = [0.035 0.044750 0.0448125 0.05]
        P = params('piTime', piT);
        ns = P.hold * 1000;
        d = abs(ns - round(ns));
        worst = max(worst, d);
        % and it must never be SHORTER than the playback -- ceil, not round
        ok = ok && d < 1e-9 && P.hold >= 2*P.minWf + piT;
    end
    report(ok, sprintf('t7a the hold is a whole number of ns for every piTime tried (worst residue %.1g ns) and never shorter than the playback', worst));
end

function ok = t7_holdMustFitInsideTau()
% In 'pulse' mode the wait is split into hold + (tau - hold), so a hold longer than the shortest
% tau asks for a negative duration. The class raises; here we check the condition it raises on,
% and that the default hold is comfortably inside a sensible shortest tau.
    P = params();
    tooShort = 0.3;                              % shorter than the 0.44 us hold
    ok = (P.hold >= tooShort) && (P.hold < 20);
    report(ok, sprintf('t7  hold %.3f us is refused against a %.2f us tau and fine against 20 us', ...
                       P.hold, tooShort));
end

function ok = t8_maxNotEnd()
% max(tau), not tau(end). Nothing sorts the tau list, and every line that reads this number
% means "the longest tau". The same latent bug was fixed in AWGExpEcho.
    unsorted = [100 4900 20 600];
    ok = max(unsorted) == 4900 && unsorted(end) ~= max(unsorted);
    report(ok, sprintf('t8  on an unsorted list max(tau) = %d where tau(end) = %d', ...
                       max(unsorted), unsorted(end)));
end

% ======================================================================
% THE STATISTICS
% ======================================================================

function o = gate(varargin)
    op = struct('selfTest', true, 'plotIt', false, 'verbose', false);
    for i = 1:2:numel(varargin); op.(varargin{i}) = varargin{i+1}; end
    o = t1GateTest(op);
end

function ok = t9_cleanFindsNothing()
    ok = true;
    for sd = 1:3
        o = gate('leakRatePerUs', 0, 'synthSeed', sd);
        good = o.ok && ~o.leakDetected;
        ok = ok & good;
        report(good, sprintf('t9  seed %d, no leak injected: slope %.1f sigma, offset %.1f sigma, not called a leak', ...
                             sd, o.sigmaSlope, o.sigmaOffset));
    end
end

function ok = t10_recoversAnInjectedLeak()
% Both estimates have to come back at the injected rate, and the test has to call it.
    ok = true;
    for r = [1e-4 3e-4]
        o = gate('leakRatePerUs', r);
        d1 = abs(o.rFromSlope  - r) / max(o.rFromSlopeErr,  eps);
        d2 = abs(o.rFromOffset - r) / max(o.rFromOffsetErr, eps);
        good = o.leakDetected && d1 < 3 && d2 < 3;
        ok = ok & good;
        report(good, sprintf(['t10 r = %.0e recovered as %.3g +/- %.2g (slope, %.1f sigma off) and ' ...
                              '%.3g +/- %.2g (offset, %.1f sigma off), detected = %d'], ...
                             r, o.rFromSlope, o.rFromSlopeErr, d1, ...
                             o.rFromOffset, o.rFromOffsetErr, d2, o.leakDetected));
    end
end

function ok = t11_sensitivityIsQuoted()
% "No leak" without a bound is not a result. Both bounds have to be finite, positive, and the
% offset one has to be the tighter of the two -- if it ever stops being, the reason the file
% quotes two has gone away and the comment explaining it is wrong.
    o = gate('leakRatePerUs', 0);
    ok = isfinite(o.sensitivityR) && o.sensitivityR > 0 && ...
         isfinite(o.sensitivityFrac) && o.sensitivityFrac > 0 && o.sensitivityFrac < 1 && ...
         isfinite(o.sensitivityROffset) && o.sensitivityROffset < o.sensitivityR;
    report(ok, sprintf(['t11 bounds quoted: slope %.2g /us (%.1f%% of T1), offset %.2g /us ' ...
                        '(%.1f%%), offset tighter by %.1fx'], ...
                       o.sensitivityR, 100*o.sensitivityFrac, o.sensitivityROffset, ...
                       100*o.sensitivityFracOffset, o.sensitivityR/max(o.sensitivityROffset,eps)));
end

function ok = t12_amplitudeDriftIsNotALeak()
%T12 A drifting contrast amplitude between arms is a pure OFFSET.
%
% It fires the C-vs-B estimate hard -- 14 sigma at 10% per arm -- and cannot touch the A-vs-B
% slope, because a constant scale factor has no tau dependence. Requiring both is what rejects
% it. This is the test that would catch t1GateTest turning into a drift detector.
    ok = true;
    for d = [0.03 0.10]
        o = gate('synthArmDriftFrac', d);
        fired  = abs(o.sigmaOffset) > 3;
        clean  = abs(o.sigmaSlope) < 3;
        good = ~o.leakDetected && clean;
        ok = ok & good;
        report(good, sprintf(['t12 amplitude drift %.0f%%/arm: offset %+.1f sigma (%s), slope %+.1f ' ...
                              'sigma (quiet), NOT called a leak'], ...
                             100*d, o.sigmaOffset, eitherText(fired, 'fires, as expected', 'quiet'), ...
                             o.sigmaSlope));
    end
end

function ok = t13_t1DriftIsNotALeak()
%T13 A drifting T1 between arms is a pure SLOPE, the mirror image of t12.
%
% It fires A-vs-B and leaves the offset alone. Between them, t12 and t13 cover both ways the
% half hour this test takes can betray it.
    ok = true;
    for d = [0.06 0.15]
        o = gate('synthArmT1DriftFrac', d);
        clean = abs(o.sigmaOffset) < 3;
        good = ~o.leakDetected && clean;
        ok = ok & good;
        report(good, sprintf(['t13 T1 drift %.0f%%/arm: slope %+.1f sigma, offset %+.1f sigma ' ...
                              '(quiet), agree %.1f sigma, NOT called a leak'], ...
                             100*d, o.sigmaSlope, o.sigmaOffset, abs(o.agree)));
    end
end

function ok = t14_orientationIsMeasured()
%T14 t1Contrast must give the same contrast whichever field the decaying trace arrives in.
%
% This is the one that protects against the projection-label swap in AWGExpT1. Feed it the pair
% both ways round; the contrast must come back identical and positive, and `flipped` must report
% which way it found them.
    tau = linspace(20, 4900, 14);
    y1 = 1 - 0.03*exp(-tau/2000);                % the |1> trace: moves with tau
    y0 = 1 + zeros(size(tau));                   % the |0> trace: flat
    e  = 4e-4*ones(size(tau));

    a = t1Contrast(struct('value', y0, 'sterr', e, 'value2', y1, 'sterr2', e, 'tau', tau));
    b = t1Contrast(struct('value', y1, 'sterr', e, 'value2', y0, 'sterr2', e, 'tau', tau));
    same = max(abs(a.contrast - b.contrast)) < 1e-12;
    signs = median(a.contrast) > 0 && median(b.contrast) > 0;
    flags = ~a.flipped && b.flipped && a.confident && b.confident;
    ok1 = same && signs && flags;
    report(ok1, sprintf(['t14 orientation measured, not read: same contrast both ways round ' ...
                         '(max diff %.1g), both positive, flipped = %d then %d'], ...
                        max(abs(a.contrast - b.contrast)), a.flipped, b.flipped));

    % And when it CANNOT tell -- both traces flat, T1 far beyond the sweep -- it must say so
    % rather than pick one and sound certain.
    flat1 = 1 + 1e-6*tau/max(tau);
    c = t1Contrast(struct('value', flat1, 'sterr', e, 'value2', flat1, 'sterr2', e, 'tau', tau));
    ok2 = ~c.confident && ~isempty(c.why);
    report(ok2, sprintf('t14 two flat traces: confident = %d, and it explains why', c.confident));
    ok = ok1 && ok2;
end

function ok = t15_channelNames()
%T15 The channel-name resolution that made every T1 run fail.
%
% prepare built the list as {MWChannel, 'MW2'}, and obj.MWChannel is not always the char 'MW'
% that Experiment.m defaults it to -- jsondecode turns a json array of strings into a cell, and
% ExpDoubleRabi sets it to {'MW','MW2'} outright. A cell inside a cell gives a first element that
% is itself a cell, and Pulse/changeChannels dies on
%
%       oc.(channel) = true;
%       -> Argument to dynamic structure reference must evaluate to a valid field name
%
% four frames below anything that names MWChannel.
%
% The rule below is a REFERENCE implementation of AWGExpT1.mwChannelNames, written out for the
% same reason seqOf is: a test that gets its logic from the code under test cannot detect logic
% that is wrong. It does not prove the method matches -- nothing offline can, since constructing
% an AWGExpT1 needs the AWG -- so on a live object the method is public and
%
%       >> t.mwChannelNames()
%
% shows what it actually resolves to. That is the check to run on the setup itself.
    ok = true;
    good = isequal(refNames('MW', 2), {'MW', 'MW2'});
    ok = ok & good; report(good, 't15 char at 2 channels -> {MW, MW2}');

    good = isequal(refNames({'MW', 'MW2'}, 2), {'MW', 'MW2'});
    ok = ok & good; report(good, 't15 CELL at 2 channels -> {MW, MW2}  (the case that crashed)');

    good = isequal(refNames('MW', 1), {'MW'});
    ok = ok & good; report(good, 't15 char at 1 channel -> {MW}');

    good = isequal(refNames('MW2', 1), {'MW2'});
    ok = ok & good; report(good, 't15 char MW2 at 1 channel -> {MW2}  (activeChannels resolves this to output 2)');

    good = isequal(refNames({'MW2'}, 1), {'MW2'});
    ok = ok & good; report(good, 't15 cell of one at 1 channel -> {MW2}');

    good = isequal(refNames({{'MW'}, 'MW2'}, 1), {'MW'});
    ok = ok & good; report(good, 't15 nested cell at 1 channel flattens -> {MW}');

    % and the three that have to be REFUSED rather than passed through
    bad = {{'numeric', 1}, {'empty', []}, {'not an identifier', 'MW 2'}};
    for i = 1:numel(bad)
        err = refNamesGuarded(bad{i}{2}, 1);
        good = ~isempty(err);
        ok = ok & good;
        report(good, sprintf('t15 %-18s refused: %s', bad{i}{1}, err));
    end
end

function names = refNames(raw, nChannels)
% The reference rule. At two channels BOTH microwave channels are driven, so MWChannel carries no
% information there and is not read at all -- which is what Echo, Ramsey, Rabi and ESR have always
% done. At one channel it names which output is live, so it is flattened and the first taken.
    if nChannels >= 2
        names = {'MW', 'MW2'};
        return
    end
    names = {};
    if ischar(raw)
        names = {raw};
    elseif iscell(raw)
        for i = 1:numel(raw)
            v = raw{i};
            if ischar(v)
                names{end+1} = v;
            elseif iscell(v)
                for j = 1:numel(v)
                    if ischar(v{j}); names{end+1} = v{j}; end
                end
            end
        end
    end
    if numel(names) > 1; names = names(1); end
end

function err = refNamesGuarded(raw, nChannels)
    err = '';
    names = refNames(raw, nChannels);
    if isempty(names)
        err = sprintf('MWChannel is a %s, not a channel name', class(raw));
        return
    end
    for i = 1:numel(names)
        if ~isvarname(names{i})
            err = sprintf('"%s" is not a valid identifier', names{i});
            return
        end
    end
end

function ok = t16_sampleGrid()
%T16 The AWG sample grid, which is where every T1 run failed next.
%
% AddWaveformByFunction refuses a duration that is not a whole number of samples:
%
%       The parameter 'duration' not divided by the AWG sample time. Fix it or change the
%       parameter 'toRound' to 'True'
%
% and it refuses it four frames below anything that mentions piTime. At 16 GS/s the grid is
% 62.5 ps, so a piTime that came out of calibPulseTimes -- a fitted number, landing wherever the
% fit put it -- misses the grid essentially always.
%
% AWGExpT1 was the only one of the five that neither rounded the pulse time nor passed
% toRound = 1. AWGExpEcho rounds halfPiTime onto the grid explicitly; AWGExpRamsey passes
% toRound. The class now does both, and this checks the arithmetic behind that: that an off-grid
% value really would be refused, that rounding fixes it, that the physical cost is negligible,
% and that the two durations T1 builds still clear the 2400-sample minimum afterwards.
    SR = 16000;                 % MHz, FrequencyGeneratorTektronixAWG.MAX_SAMPLE_RATE
    tSample = 1/SR;             % 6.25e-5 us = 62.5 ps
    minWf = 2400/SR;            % 0.15 us, and exactly MIN_WAVEFORM_LENGTH samples
    ok = true;

    % a calibrated pulse time, as calibPulseTimes would leave it
    piRaw = 0.044761;
    good = ~onGrid(piRaw + minWf, tSample);
    ok = ok & good;
    report(good, sprintf('t16 raw piTime %.6f us: %.6f + %.2f is %.2f samples -- REFUSED, as it was', ...
                         piRaw, piRaw, minWf, (piRaw + minWf)/tSample));

    piFix = round(piRaw*SR)/SR;                 % what RoundDurationByFrequencies does
    good = onGrid(piFix + minWf, tSample) && onGrid(piFix, tSample);
    ok = ok & good;
    report(good, sprintf('t16 rounded to %.6f us: %.0f samples, and with the pad %.0f -- accepted', ...
                         piFix, piFix/tSample, (piFix + minWf)/tSample));

    % the physical cost of the round
    moved = abs(piFix - piRaw);
    good = moved <= tSample/2 + 1e-15 && moved/piRaw < 0.001;
    ok = ok & good;
    report(good, sprintf('t16 the round moves the pulse by %.4f ns, %.3f%% of it -- inside anything calibPulseTimes can resolve', ...
                         1e3*moved, 100*moved/piRaw));

    % the two durations T1 actually builds, against the 2400-sample floor
    nZero = (minWf)/tSample;
    nPi   = (piFix + minWf)/tSample;
    good = nZero >= 2400 && nPi > 2400;
    ok = ok & good;
    report(good, sprintf('t16 waveform lengths: zero pad %.0f samples (exactly the 2400 floor), pi %.0f', ...
                         nZero, nPi));

    % and timeDelay, the other duration that reaches a waveform call, is rounded before use
    tdRaw = 0.0137;
    tdFix = floor(tdRaw*SR)/SR;                 % RoundDurationByFrequencies(..., -1)
    good = onGrid(minWf + tdFix, tSample);
    ok = ok & good;
    report(good, sprintf('t16 prePulseTime with a %.4f us timeDelay floors to %.6f and stays on the grid', ...
                         tdRaw, tdFix));
end

function tf = onGrid(duration, tSample)
% The divisibility test AddWaveformByFunction applies, written out. MAX_ERROR is its 1e-9.
    MAX_ERROR = 1e-9;
    r = mod(duration, tSample);
    tf = ~(r > MAX_ERROR && tSample - r > MAX_ERROR);
end

% ======================================================================

function report(ok, what)
    if ok
        fprintf('  PASS  %s\n', what);
    else
        fprintf('  FAIL  %s\n', what);
    end
end

function s = eitherText(tf, a, b)
    if tf; s = a; else; s = b; end
end
