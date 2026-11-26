function out = calibMwPair(opts)
%CALIBMWPAIR amplitudeDiff and channelPhase for the two microwave channels, from Rabi.
%
%   out = calibMwPair
%   out = calibMwPair(struct('amplitude', -16))
%   out = calibMwPair(struct('phases', 0:30:330, 'nBalanceReps', 3))
%
% LEVEL 4. Run after calibB1Rabi, with the coils ON and at the field it chose, so the NVs
% sit at the antenna's response peak and the drive is as strong as this setup gets.
% Replaces calibMwBalance and calibMwPhase, which used CW dip depth; this uses the Rabi
% frequency, which is a direct and absolutely calibrated measure of B1.
%
% Produces two numbers for expProperties.json:
%
%       amplitudeDiff = [0 delta]        dB to add to channel 2 so the two match
%       channelPhase  = [0 phiMax]       degrees on channel 2 for maximum combined drive
%
% and requires timeDelay = 0, since timeDelay is the other, much coarser phase mechanism
% and the two are additive.
%
% ===================================================================
% WHY ONE CHANNEL AT A TIME IS THE RIGHT WAY TO DO THE BALANCE
% ===================================================================
%
% The worry is reasonable: if each channel is characterised alone, does the answer still
% hold when both are driven into the same resonator? It does, and the reason is that
% everything between the AWG outputs and the diamond is a LINEAR network. The transfer
% function from AWG output n to B1 at the NV -- call it H1 and H2 -- is a fixed property of
% the cables, circulators, combiner and resonator. It cannot depend on what the other port
% is doing. So the field is the phasor sum a1*H1 + a2*H2, each term exactly what was
% measured with the other channel off, and equalising |a1*H1| with |a2*H2| separately makes
% them equal when both are on. The resonator has no way to know the second channel arrived.
%
% On this setup each channel has its OWN amplifier and its OWN circulator, which makes the
% argument stronger rather than weaker. A circulator is precisely the component that stops
% channel 2's reflected power reaching channel 1's amplifier output, so the one path that
% could have coupled the two chains is terminated in a load. And because each amplifier
% only ever sees its own channel, turning the other one on does not move its operating
% point at all.
%
% What per-channel amplifiers DO mean is that the balance is LEVEL dependent: two
% amplifiers compress at different points, so a delta measured at -30 dBm need not hold at
% -16. So this calibrates at the amplitude you will run at, and complains if you ask for one
% and then use another. That is the real caveat, and it has nothing to do with both channels
% being on.
%
% None of the above has to be believed, because stage 3 measures it. In a linear network at
% the optimum phase,
%
%       Omega_both  =  Omega1 + Omega2      exactly
%
% so R = Omega_both/(Omega1+Omega2) is a direct test of superposition. R = 1 means the
% one-at-a-time calibration is valid. R < 1 means either the phase is off the optimum or
% something both channels share is compressing, and stage 2 has already found the phase, so
% a shortfall at the fitted maximum is compression.
%
% ===================================================================
% EACH CHANNEL IS MEASURED AT THE LEVEL IT WILL PLAY  (2026-11-23)
% ===================================================================
%
% Stage 1 used to measure both channels at the SAME dBm and then play channel 2 delta dB
% harder, assuming its Rabi frequency would grow by exactly 10^(delta/20). Below compression
% that holds; above it, it does not, and the 11-10 data show it: at -16 dBm the closure
% came back R = 0.910 +/- 0.003, and the numbers fit channel 2 gaining only ~0.4 dB of
% field for its extra 2.1 dB of drive -- the pair ran ~1.7 dB unbalanced while the json said
% balanced. At -32 dBm, R = 0.986 +/- 0.007: linear, and the old way was fine there.
%
% So stage 1 now ITERATES (iterateBalance, default true). Round 0 is the old measurement,
% both channels at the operating dBm, which gives a first delta. Every later round measures
% channel 1 at the operating dBm and channel 2 at operating + delta -- the level it will
% actually play -- and corrects delta by the remaining mismatch divided by channel 2's
% local slope (dB of Rabi frequency per dB of drive, from its last two levels). It stops
% when the mismatch is below max(2 sigma, balanceTolDb), usually after one round in the
% linear regime and two or three in compression. The delta it returns is one at which
% channel 2 was MEASURED, not extrapolated to, and the closure test (stage 3) compares the
% both-on Rabi frequency with those measured values -- a pure superposition test.
%
% calibMwQuadrature (level4.m cell 4.2b) then checks the result with BOTH channels on,
% which no single-channel measurement can.
%
% ===================================================================
% WHY THE ANTIPHASE NULL IS A CHECK HERE AND NOT THE MEASUREMENT
% ===================================================================
%
% At antiphase Omega_total = |Omega1 - Omega2|, which is exactly zero at balance, and a
% null measured with both channels driven is the most direct answer to the question above.
% It is stage 4. But it is a CHECK, not the primary measurement, and it is worth being clear
% why -- the opposite conclusion would have been reasonable a calibration ago.
%
% When the observable was CW dip depth, one-at-a-time was hopeless: the depth is some
% unknown monotone function of B1, so comparing two channels meant comparing two points on
% a curve nobody had measured. A null needs no such curve, which is why calibMwBalance was
% built as one.
%
% The Rabi frequency changes that completely. It IS B1, up to one fixed constant, and it
% comes out of the fit to better than half a percent. So comparing Omega1 with Omega2
% directly gives the ratio to about 0.7 %, i.e. 0.06 dB, in four runs.
%
% The null cannot match that, for a reason that is physics and not technique. Near balance
% the residual drive is small, so the residual Rabi frequency is slow, and a slow
% oscillation cannot be measured on a sweep shorter than its period -- nor beyond where the
% Rabi envelope has decayed. A 1 % imbalance on a 16 MHz combined drive leaves 0.08 MHz, a
% 12 us period. Over the standard 1 us sweep nothing is resolvable below about 0.5 MHz,
% which bounds the imbalance only to roughly 0.5 dB. Stage 4 therefore uses a LONGER sweep
% (antiphaseTauMax) to push that bound down, reports it honestly as a bound when the null is
% deeper than the sweep can see, and never overrides the direct measurement.
%
% So: stage 1 measures the balance, stage 4 confirms it is not wildly wrong and confirms the
% null is where the model says. If they disagree, that is a real finding and it is printed.
%
% ===================================================================
% WHY THE PHASE IS NOT CIRCULAR WITH THE BALANCE
% ===================================================================
%
% Omega_total^2 = Omega1^2 + Omega2^2 + 2*Omega1*Omega2*cos(phi - phi0) is extremal at phi0
% and phi0+180 for ANY Omega1 and Omega2. The phase optimum simply does not depend on the
% balance, so there is no loop to iterate: a rough balance is enough to find the phase
% exactly, and the phase does not have to be revisited after the balance is refined.
%
% The phase is fitted from the WHOLE curve rather than from a peak or a minimum -- see
% phaseScanFit, where the model is linear in three coefficients and therefore has no start
% point and no local optima. It is worth knowing that the peak alone would be a poor target:
% the derivative vanishes at both extrema, and a 20 degree phase error costs only 1.5 % of
% the drive at the peak, so peak-finding is insensitive. The null is where the phase
% information lives, and fitting the whole curve uses it.
%
% ===================================================================
% OPTIONS
% ===================================================================
%
%   amplitude      operating dBm for channel 1, the level everything is calibrated at
%                  (default -16, the Rabi standard). The balance is level dependent because
%                  the amplifiers are per channel, so this should be what you run at.
%   frequency      MHz. Default [] means use whatever loadJsonParams supplies, which after
%                  calibB1Rabi is the resonator peak. Printed either way.
%   deltaStart     skip stage 1 and use this dB offset instead (default [] = measure it)
%   phases         degrees to scan on channel 2 (default 0:45:315, 8 points). The model has
%                  three parameters and is exact, so 8 is generous; 12 buys a factor 1.4 on
%                  the phase error for 50 % more time.
%   nBalanceReps   how many channel-1/channel-2 pairs in round 0 of stage 1 (default 2). The
%                  order alternates within each pair so a steady drift cancels.
%   iterateBalance measure channel 2 at the level it will play and iterate (default true).
%                  false = the old single round, with channel 2's Rabi frequency at
%                  operating + delta EXTRAPOLATED linearly from the operating level.
%   balanceTolDb   stop when |mismatch| <= max(2 sigma, this) (default 0.05 dB of field)
%   maxBalanceRounds  refinement rounds after round 0 (default 4)
%   roundReps      channel-1/channel-2 pairs per refinement round (default 1)
%   singleChannel  how stage 1 measures one channel (default 'nChannels1'):
%                  'nChannels1'  a ONE-channel Rabi, as the run script's single-channel cell:
%                                nChannels = 1, MWChannel = 'MW' (output 1) or 'MW2'
%                                (output 2), amplitude = the operating level. Only that
%                                output plays and only its MW switch opens. With one channel
%                                the waveform is scaled by amplitude(1) whichever output
%                                plays it, so both channels get exactly the same dBm.
%                  'offAmplitude' the old way: a two-channel Rabi with the other channel
%                                held at offAmplitude. Kept for comparison; it did not give
%                                usable stage-1 Rabis on this setup.
%   singleDetection  us, detection window for the stage-1 single-channel Rabis. [] (default)
%                  = the json's detectionDuration, as the single-channel cell uses. The Rabi
%                  FREQUENCY does not depend on the window, only the contrast and noise do,
%                  so this does not bias the balance.
%   offAmplitude   dBm for the held-off channel, 'offAmplitude' mode only (default -70).
%   antiphaseTauMax  us, the sweep length for stage 4 (default 4.0). Longer resolves a
%                  smaller residual and so a tighter bound on the imbalance, until the Rabi
%                  decay makes it pointless -- the measured decay is used to say where that is.
%   doCoarse, doPhase, doClosure, doAntiphase   run each stage (default all true)
%
%   RABI STANDARD -- same as calibB1Rabi. Changing these changes what is being compared.
%   tau, repeats, averages, useSequenceWalk, passesPerRun, walkWarmupTriggers,
%   resyncEveryPass, rabiDetectionMHz, nFreq, minContrast, maxFErrFrac
%   newOrderEveryAverage  walk-mode Rabi runs: a new random tau order before every average
%                  (true/false). Default [] = AWGExpRabi's own setting (false). Costs one
%                  sequence upload + arm per average.
%
%   plotIt, isTracking, jsonLocation
%
% Changes no calibration file itself. level4.m cell 4.3 writes amplitudeDiff, channelPhase and
% timeDelay = 0 into expProperties.json with expPropsWrite (only when out.ok).
% Leaves the field alone entirely -- it never touches the Helmholtz.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct( ...
        'amplitude', -16, 'frequency', [], 'deltaStart', [], ...
        'phases', 0:45:315, 'nBalanceReps', 2, 'offAmplitude', -70, ...
        'iterateBalance', true, 'balanceTolDb', 0.05, 'maxBalanceRounds', 4, 'roundReps', 1, ...
        'singleChannel', 'nChannels1', 'singleDetection', [], ...
        'antiphaseTauMax', 4.0, ...
        'doCoarse', true, 'doPhase', true, 'doClosure', true, 'doAntiphase', true, ...
        'tau', [], 'repeats', 300, 'averages', 2, 'useSequenceWalk', true, ...
        'passesPerRun', [], 'walkWarmupTriggers', 1, 'resyncEveryPass', true, ...
        'newOrderEveryAverage', [], ...
        'rabiDetectionMHz', 0.3, 'nFreq', 1, 'minContrast', 0.005, 'maxFErrFrac', 0.05, ...
        'plotIt', true, 'isTracking', 0, ...
        'jsonLocation', ''), ...
        'calibMwPair');

    % ---- where the calibration lives ------------------------------------
    % Resolved by setup7Paths rather than written down here. Twelve files used to carry
    % their own copy of the path, so moving the results folder meant twelve edits, and
    % missing one silently loaded the OLD calibration into this one experiment while
    % everything else used the new one -- two calibrations live at once, with no symptom
    % except data that will not reproduce.
    %
    % This function CANNOT run without a calibration, so a failure to resolve is fatal here
    % and setup7Paths says what is wrong with the folder.
    %
    % Pass jsonLocation explicitly ONLY to pin a run to an older calibration on purpose:
    %     setup7Paths('2026-08-14').json
    if isempty(opts.jsonLocation)
        opts.jsonLocation = setup7Paths().json;
    end

    if isempty(opts.tau); opts.tau = round(0.02 : 0.01 : 1.0, 3); end
    out = struct('ok', false, 'opts', opts);

    fprintf('\n=== calibMwPair ===\n');
    nRefine = opts.iterateBalance * opts.maxBalanceRounds * 2 * opts.roundReps;
    nRuns = opts.doCoarse*2*opts.nBalanceReps + opts.doPhase*numel(opts.phases) + ...
            opts.doClosure + opts.doAntiphase;
    fprintf(['Channel 1 at %g dBm throughout. %d Rabi runs in all, plus up to %d for the\n' ...
             'balance rounds -- every one reloads the AWG, because both amplitude and\n' ...
             'channelPhase change the waveform samples, so budget roughly %.0f-%.0f minutes.\n'], ...
            opts.amplitude, nRuns, nRefine, 1.5*nRuns, 1.5*(nRuns + nRefine));
    fprintf(['The field is NOT touched. Run calibB1Rabi first and leave the coils where it\n' ...
             'put them, so this measures the pair at the operating point.\n']);

    % ---------------- STAGE 1: the balance, one channel at a time -------------
    delta = opts.deltaStart; deltaErr = NaN;
    O1 = NaN; O1e = NaN; O2 = NaN; O2e = NaN;
    O2play = NaN; O2playE = NaN;            % channel 2 MEASURED at operating + delta
    out.coarse = [];
    out.balanceRounds = [];
    out.balanceConverged = false;
    out.slope2 = NaN;
    if opts.doCoarse && isempty(opts.deltaStart)
        if ~any(strcmp(opts.singleChannel, {'nChannels1', 'offAmplitude'}))
            error('calibMwPair: singleChannel must be ''nChannels1'' or ''offAmplitude''.');
        end
        if strcmp(opts.singleChannel, 'nChannels1')
            fprintf(['\n--- stage 1: each channel alone (one-channel Rabi: nChannels = 1, ' ...
                     'MWChannel MW / MW2) ---\n']);
        else
            fprintf('\n--- stage 1: each channel alone (other channel held at %g dBm) ---\n', ...
                    opts.offAmplitude);
        end
        fprintf('round 0: both channels at %g dBm\n', opts.amplitude);
        fprintf('%6s %8s %12s %12s %9s %8s  %s\n', 'rep', 'channel', 'amps/dBm', ...
                'fRabi/MHz', 'pi/ns', 'contr%', 'note');
        [r1, r2, ms] = balancePairs(opts.nBalanceReps, 0, opts);
        out.coarse = [out.coarse ms];
        [O1, O1e] = combine(r1);
        [O2, O2e] = combine(r2);
        if ~isfinite(O1) || ~isfinite(O2)
            error(['calibMwPair: one of the channels gave no usable Rabi. Check the notes ' ...
                   'above -- a channel with no contrast at all usually means its switch or ' ...
                   'amplifier is off, not that the calibration failed.']);
        end
        % Omega is proportional to amplitude in VOLTS, and dB between two amplitudes is
        % 20*log10 of their ratio. Channel 2 needs delta MORE dB when it is the weaker one.
        delta = 20*log10(O1/O2);
        deltaErr = (20/log(10)) * hypot(O1e/O1, O2e/O2);
        fprintf('\n  channel 1  %.4f +/- %.4f MHz   (%d runs)\n', O1, O1e, numel(r1));
        fprintf('  channel 2  %.4f +/- %.4f MHz   (%d runs)   both at %g dBm\n', O2, O2e, ...
                numel(r2), opts.amplitude);
        fprintf('  round 0: delta = %+.3f +/- %.3f dB\n', delta, deltaErr);
        spreadWarn(r1, 'channel 1'); spreadWarn(r2, 'channel 2');
        rounds = struct('round', 0, 'delta', 0, 'omega1', O1, 'omega1Err', O1e, ...
                        'omega2', O2, 'omega2Err', O2e, 'mismatchDb', delta, ...
                        'mismatchErr', deltaErr, 'slope', NaN);
        r1all = r1;

        if ~opts.iterateBalance
            % The old way: channel 2's Rabi frequency at operating + delta EXTRAPOLATED, which
            % is exact only below compression. Kept for comparison.
            O2play = O2 * 10^(delta/20);  O2playE = O2e * 10^(delta/20);
            fprintf(['  iterateBalance = false: channel 2 at %+.3f dB was NOT measured; its Rabi\n' ...
                     '  frequency there is extrapolated linearly (%.4f MHz). Exact only below\n' ...
                     '  compression.\n'], delta, O2play);
        else
            % Channel 2 at the level it will play, until the mismatch is gone. Every round
            % gives a point (delta_k, mismatch_k) -- the round's own channel 1 against its own
            % channel 2, taken back to back, so drift between rounds cancels. The next delta
            % is the zero of a curve through ALL the points so far (balanceRoot: a line for
            % two, a weighted quadratic from three). A secant through the last two points
            % undershoots in compression, from the same side every round.
            dl = 0;  mm = delta;  mE = deltaErr;
            slope = 1;
            for k = 1:opts.maxBalanceRounds
                fprintf('\nround %d: channel 1 at %g dBm, channel 2 at %g %+.3f = %.3f dBm\n', ...
                        k, opts.amplitude, opts.amplitude, delta, opts.amplitude + delta);
                fprintf('%6s %8s %12s %12s %9s %8s  %s\n', 'rep', 'channel', 'amps/dBm', ...
                        'fRabi/MHz', 'pi/ns', 'contr%', 'note');
                [q1, q2, ms] = balancePairs(opts.roundReps, delta, opts);
                out.coarse = [out.coarse ms];
                r1all = [r1all q1];                                     %#ok<AGROW>
                [a1, a1e] = combine(q1);
                [a2, a2e] = combine(q2);
                if ~isfinite(a1) || ~isfinite(a2)
                    fprintf('  round %d gave no usable Rabi on one channel -- stopping here.\n', k);
                    break
                end
                m  = 20*log10(a1/a2);
                me = (20/log(10)) * hypot(a1e/a1, a2e/a2);
                dl(end+1) = delta;  mm(end+1) = m;  mE(end+1) = me;     %#ok<AGROW>
                [root, rootErr, sl] = balanceRoot(dl, mm, mE);
                if isfinite(sl); slope = sl; end
                rounds(end+1) = struct('round', k, 'delta', delta, 'omega1', a1, ...
                    'omega1Err', a1e, 'omega2', a2, 'omega2Err', a2e, 'mismatchDb', m, ...
                    'mismatchErr', me, 'slope', slope);                 %#ok<AGROW>
                fprintf(['  channel 2 at %+.3f dB: %.4f MHz against channel 1 %.4f MHz -> ' ...
                         'mismatch %+.3f +/- %.3f dB of field (slope %.2f dB/dB)\n'], ...
                        delta, a2, a1, m, me, slope);
                if ~isfinite(root)
                    root = delta + m / slope;
                    rootErr = me / slope;
                end
                if abs(m) <= max(2*me, opts.balanceTolDb)
                    % Converged: the mismatch at the level just played is inside the noise
                    % (or the tolerance). The curve's zero is still the best estimate of the
                    % balance -- it uses every round, not one -- and it is within about
                    % 2 sigma of the measured level, so channel 2's Rabi frequency there is
                    % the measured one moved along the measured slope.
                    out.balanceConverged = true;
                    step = root - delta;
                    O2play = a2 * 10^(slope * step / 20);  O2playE = a2e * 10^(slope * step / 20);
                    fprintf(['  converged: mismatch inside the noise; zero of the fit %+.3f ' ...
                             '(%+.3f from the level played)\n'], root, step);
                    delta = root;
                    % Floor: one round's own precision, divided by the slope. The fit's
                    % covariance alone came out up to 2x optimistic in heavy compression
                    % (simulated), where channel 2 moves only ~0.3 dB of field per dB of drive.
                    deltaErr = max(rootErr, me / slope);
                    if ~isfinite(deltaErr); deltaErr = me / slope; end
                    break
                end
                step = root - delta;
                step = sign(step) * min(abs(step), 6);
                delta = delta + step;
            end
            out.slope2 = slope;
            if ~out.balanceConverged
                fprintf(['\n  *** NOT CONVERGED after %d rounds. The last measured mismatch was ' ...
                         '%+.3f dB at delta %+.3f.\n  Usually drift between runs larger than ' ...
                         'the tolerance; raise roundReps or maxBalanceRounds.\n'], ...
                        opts.maxBalanceRounds, rounds(end).mismatchDb, rounds(end).delta);
                % The best measured point, so the stages below play a level that was measured.
                [~, ib] = min(abs([rounds(2:end).mismatchDb]));
                if ~isempty(ib)
                    rb = rounds(ib + 1);
                    delta = rb.delta;  O2play = rb.omega2;  O2playE = rb.omega2Err;
                    deltaErr = hypot(rb.mismatchErr, rb.mismatchDb) / slope;
                end
            end
            [O1, O1e] = combine(r1all);       % channel 1 was at the operating dBm every time
        end
        out.balanceRounds = rounds;
        fprintf('\n  amplitudeDiff = [0 %+.3f] +/- %.3f dB\n', delta, deltaErr);
        if delta > 0
            fprintf('  i.e. channel 2 needs %+.3f dB MORE than channel 1 to match it.\n', delta);
        else
            fprintf('  i.e. channel 2 needs %.3f dB LESS than channel 1 to match it.\n', -delta);
        end
    end
    % ---- still no balance? fall back rather than refusing ----------------
    % Skipping stage 1 to run a phase scan on its own is a normal thing to want, and it used
    % to error here. It should not: the phase optimum does NOT depend on the balance --
    % Omega^2 = O1^2 + O2^2 + 2*O1*O2*cos(phi - phi0) is extremal at phi0 for ANY O1 and O2 --
    % so an approximate delta, or none at all, still finds the right phase. All the balance
    % changes is how deep the null goes, and the fit uses the whole curve rather than the
    % null alone.
    %
    % So: take it from the measured amplitude curve if there is one, otherwise from the json,
    % otherwise run unbalanced and say so.
    if isempty(delta)
        try
            dv = amplitudeDiffAt(opts.amplitude, struct('verbose', false));
            delta = dv(2);
            fprintf(['\nNo stage 1 this run, so amplitudeDiff came from the measured curve: ' ...
                     '%+.3f dB at %g dBm.\n'], delta, opts.amplitude);
        catch
            try
                j = jsondecode(fileread(opts.jsonLocation));
                if isfield(j, 'amplitudeDiff') && numel(j.amplitudeDiff) >= 2
                    delta = j.amplitudeDiff(2);
                    fprintf(['\nNo stage 1 this run and no amplitude curve, so amplitudeDiff ' ...
                             'came from the json: %+.3f dB.\n'], delta);
                end
            catch
            end
        end
    end
    if isempty(delta)
        delta = 0;
        fprintf(['\nNo balance available from anywhere, so running with amplitudeDiff = 0.\n' ...
                 'That is fine for a PHASE scan -- the optimum phase does not depend on the\n' ...
                 'balance -- but stage 1 or a deltaStart is needed for anything else.\n']);
    end
    out.delta = delta; out.deltaErr = deltaErr;
    out.omega1 = O1; out.omega1Err = O1e; out.omega2 = O2; out.omega2Err = O2e;
    out.omega2AtPlay = O2play; out.omega2AtPlayErr = O2playE;

    ampsBal = [opts.amplitude, opts.amplitude + delta];
    % Stages 2-4 play channel 2 at operating + delta, so the sum they are compared with is
    % Omega1 + Omega2 AT THAT LEVEL. Stage 1 measured it there (iterateBalance); without
    % stage 1 it is not known, and the sum stays NaN.
    O2bal = O2play; O2balE = O2playE;
    out.sumBalanced = O1 + O2bal; out.sumBalancedErr = hypot(O1e, O2balE);

    % ---------------- STAGE 2: the phase, both channels on -------------------
    out.phaseScan = []; pf = [];
    if opts.doPhase
        fprintf('\n--- stage 2: phase scan, both channels driven at [%g %+.3f] dBm ---\n', ...
                ampsBal(1), ampsBal(2));
        fprintf('%9s %12s %9s %8s %8s  %s\n', 'phase/deg', 'fRabi/MHz', 'pi/ns', ...
                'contr%', 'R^2', 'note');
        for i = 1:numel(opts.phases)
            m = runRabi(ampsBal, opts.phases(i), opts.tau, opts);
            m.phaseDeg = opts.phases(i);
            fprintf('%9.1f %6.3f+-%-5.3f %9.1f %8.2f %8.4f  %s\n', ...
                    m.phaseDeg, m.f, m.fErr, m.piNs, 100*m.contrast, m.r2, m.note);
            out.phaseScan = [out.phaseScan m];                   %#ok<AGROW>
        end
        okp = logical([out.phaseScan.ok]);
        pf = phaseScanFit([out.phaseScan(okp).phaseDeg], [out.phaseScan(okp).f], ...
                          [out.phaseScan(okp).fErr], struct('verbose', true));
        out.fit = pf;
        if pf.ok && opts.doCoarse
            crossCheck(pf, O1, O2bal, delta);
        end
    end

    phiMax = 0; phiNull = 180;
    if isstruct(pf) && isfield(pf, 'ok') && pf.ok
        phiMax = pf.phiMax; phiNull = pf.phiNull;
    elseif opts.doPhase
        fprintf(['\n  The phase fit failed, so the phase is left at 0 for the stages below and\n' ...
                 '  no channelPhase is recommended. Fix stage 2 before using either number.\n']);
    end
    out.phiMax = phiMax; out.phiNull = phiNull;

    % ---------------- STAGE 3: the closure test ------------------------------
    out.closure = [];
    if opts.doClosure
        fprintf('\n--- stage 3: both channels at the fitted maximum (%.1f deg) ---\n', phiMax);
        m = runRabi(ampsBal, phiMax, opts.tau, opts);
        out.closure = m;
        fprintf('  fRabi both channels  %.4f +/- %.4f MHz   (pi = %.1f ns, contrast %.2f %%)\n', ...
                m.f, m.fErr, m.piNs, 100*m.contrast);
        if opts.doCoarse && m.ok
            S = O1 + O2bal; Se = hypot(O1e, O2balE);   % at the stage-3 amplitudes
            R = m.f/S; Re = R*hypot(m.fErr/m.f, Se/S);
            out.closureRatio = R; out.closureRatioErr = Re;
            reportClosure(R, Re, m.f, S, Se);
        end
    end

    % ---------------- STAGE 4: the antiphase check --------------------------
    out.antiphase = [];
    if opts.doAntiphase
        dtau = max(round(opts.antiphaseTauMax/numel(opts.tau), 3), 0.005);
        tauLong = round(dtau : dtau : opts.antiphaseTauMax, 3);
        fprintf('\n--- stage 4: antiphase (%.1f deg), %g us sweep ---\n', ...
                phiNull, max(tauLong));
        fprintf(['A longer sweep than the standard, because the point is to resolve a SLOW\n' ...
                 'residual: the smallest frequency a sweep can establish is about one half\n' ...
                 'cycle over its length, so %g us reaches %.3f MHz.\n'], ...
                max(tauLong), 1/(2*(max(tauLong) - min(tauLong))));
        m = runRabi(ampsBal, phiNull, tauLong, opts);
        out.antiphase = m;
        S4 = O1 + O2bal;                            % at the stage-4 amplitudes
        if ~isfinite(S4) && isstruct(pf) && isfield(pf, 'sumOmega'); S4 = pf.sumOmega; end
        if isfinite(S4)
            reportAntiphase(m, tauLong, S4, delta, deltaErr, out);
        else
            fprintf(['  No value for Omega1+Omega2 (stage 1 and stage 2 were both skipped), so\n' ...
                     '  the residual cannot be turned into an imbalance. Measured residual: '...
                     '%.4f MHz, note "%s".\n'], m.f, m.note);
        end
    end

    % ---------------- what to write ----------------------------------------
    out.amplitudeDiff = [0 delta];
    out.channelPhase  = [0 phiMax];
    out.ok = isfinite(delta) && (~opts.doPhase || (isstruct(pf) && isfield(pf,'ok') && pf.ok));
    % A balance that was iterated but never closed is not a calibration: the level written
    % would be one channel 2 was never shown to match at.
    ranIteration = opts.doCoarse && isempty(opts.deltaStart) && opts.iterateBalance;
    if ranIteration && ~out.balanceConverged
        out.ok = false;
    end
    printHandover(out, opts);

    if opts.plotIt && ~isempty(out.phaseScan); drawIt(out, opts); end
end

% ======================================================================
% ONE RABI
% ======================================================================

function m = runRabi(amps, phaseDeg, tauList, opts)
% One Rabi at explicit per-channel amplitudes and an explicit relative phase.
%
% amplitudeDiff is forced to [0 0] and the amplitudes handed over as a full vector. The
% amplitude setter ADDS amplitudeDiff to whatever it is given, so a non-zero amplitudeDiff
% left in place would be applied on top of the offset being tested, and the run would
% measure a different delta from the one in the table. Explicit vectors cannot do that.
    m = struct('ok', false, 'note', '', 'f', NaN, 'fErr', NaN, 'contrast', NaN, ...
               'r2', NaN, 'r2Ceiling', NaN, 'chi2', NaN, 'piNs', NaN, 'piErrNs', NaN, ...
               'T', NaN, 'residAutocorr', NaN, 'amps', amps, 'phaseDeg', phaseDeg, ...
               'seconds', NaN, 'tauMax', max(tauList), 'frequency', NaN);
    % 'frequency' is initialised HERE and not only on the success path. These structs get
    % concatenated into an array, and MATLAB requires identical field sets to do that -- a
    % field added after an early return would make one failed run break the whole scan.

    b = AWGExpRabi;
    b.loadJsonParams('all', opts.jsonLocation);
    b.nChannels = 2;                    % before amplitude, always
    b.amplitudeDiff = [0 0];
    b.amplitude = amps;
    if ~isempty(opts.frequency); b.frequency = opts.frequency; end
    b.timeDelay = 0;                    % channelPhase is the phase knob; they are additive
    b.channelPhase = [0 phaseDeg];
    b.repeats   = opts.repeats;
    b.averages  = opts.averages;
    b.isTracking = opts.isTracking;
    b.balancedSequence = 1;
    b.detectionDuration = opts.rabiDetectionMHz;
    b.referenceDetectionDuration = b.detectionDuration;
    b.shouldAutosave = false;
    b.tau = tauList;
    b.useSequenceWalk    = opts.useSequenceWalk;
    b.passesPerRun       = opts.passesPerRun;
    b.walkWarmupTriggers = opts.walkWarmupTriggers;
    b.resyncEveryPass    = opts.resyncEveryPass;
    if ~isempty(opts.newOrderEveryAverage)
        b.newOrderEveryAverage = logical(opts.newOrderEveryAverage);
    end
    m = runAndFit(b, m, opts);
end

function [root, rootErr, slope] = balanceRoot(dl, mm, mE)
% Where the mismatch (dB of field) crosses zero as a function of the channel-2 offset
% played (dB of drive), from EVERY balance round so far. Two points: the straight line
% through them. Three or more: a weighted quadratic through all of them -- in compression
% the curve is concave, and any straight line through points left of the zero is too steep
% there, so it falls short from the same side every round (simulated on this setup's
% compression: line bias -0.1 to -0.2 dB, quadratic bias 0.00 dB, rms 0.05 dB linear and
% 0.17 dB at -16 dBm). slope = the channel's local dB-per-dB at the zero, clamped to
% [0.15 1.5].
    root = NaN;  rootErr = NaN;  slope = NaN;
    n = numel(dl);
    if n < 2; return; end
    x = dl(:);  y = mm(:);  e = mE(:);
    good = e > 0 & isfinite(e);
    if ~any(good); return; end
    e(~good) = max(e(good));
    if max(x) - min(x) < 0.02; return; end
    dLast = x(end);
    if n >= 3
        X = [ones(n, 1), x, x.^2];
    else
        X = [ones(n, 1), x];
    end
    W = diag(1 ./ e.^2);
    Cv = inv(X' * W * X);
    c = Cv * (X' * W * y);
    if n > size(X, 2)
        chi2nu = sum(((y - X * c) ./ e).^2) / (n - size(X, 2));
        Cv = Cv * max(1, chi2nu);
    end
    if n >= 3 && c(3) ~= 0
        disc = c(2)^2 - 4 * c(3) * c(1);
        if disc >= 0
            r = (-c(2) + [-1 1] * sqrt(disc)) / (2 * c(3));
            [~, i] = min(abs(r - dLast));
            r = r(i);
        else
            r = NaN;
        end
        der = c(2) + 2 * c(3) * r;           % d(mismatch)/d(offset) at the zero
        J = [1, r, r^2];
    else
        r = -c(1) / c(2);
        der = c(2);
        J = [1, r];
    end
    if isfinite(r) && der < 0 && -der >= 0.15
        root = r;
        slope = min(-der, 1.5);
        rootErr = sqrt(max(J * Cv * J', 0)) / abs(der);
    else
        % No sensible zero (noise between close levels, or a flat stretch): step from the
        % last point with a clamped local slope instead.
        slope = 0.5;
        if n >= 2 && x(end) ~= x(end-1)
            slope = min(max(-(y(end) - y(end-1)) / (x(end) - x(end-1)), 0.15), 1.5);
        end
        root = x(end) + y(end) / slope;
        rootErr = e(end) / slope;
    end
end

function [r1, r2, allRuns] = balancePairs(nReps, delta, opts)
% nReps channel-1/channel-2 pairs: channel 1 at the operating dBm, channel 2 at operating +
% delta. Which channel goes first alternates from pair to pair: if the drive or the
% fluorescence drifts steadily, a fixed order puts all of the drift into the ratio, and
% alternating puts half into each channel, where it cancels.
    r1 = []; r2 = []; allRuns = [];
    amp = [opts.amplitude, opts.amplitude + delta];
    for rep = 1:nReps
        order = [1 2];
        if mod(rep, 2) == 0; order = [2 1]; end
        for k = order
            if strcmp(opts.singleChannel, 'nChannels1')
                m = runRabiSingle(k, opts.tau, opts, amp(k));
                ampLabel = sprintf('out%d %.2f', k, amp(k));
            else
                amps = opts.offAmplitude * ones(1, 2);
                amps(k) = amp(k);
                m = runRabi(amps, 0, opts.tau, opts);
                ampLabel = mat2str(amps, 4);
            end
            fprintf('%6d %8d %12s %6.3f+-%-5.3f %9.1f %8.2f  %s\n', rep, k, ...
                    ampLabel, m.f, m.fErr, m.piNs, 100*m.contrast, m.note);
            if k == 1; r1 = [r1 m]; else; r2 = [r2 m]; end   %#ok<AGROW>
            allRuns = [allRuns m];                           %#ok<AGROW>
        end
    end
end

function m = runRabiSingle(k, tauList, opts, ampDbm)
% One Rabi on ONE AWG output (k = 1 or 2), built like the run script's single-channel cell:
%   nChannels = 1 BEFORE amplitude, MWChannel = 'MW' / 'MW2', amplitude = ampDbm (default:
%   the operating level).
% The other output is not played and its MW switch stays shut, so nothing of it reaches the
% resonator. amplitudeDiff and channelPhase are zeroed: with one channel the waveform uses
% amplitude(1) and phase(1) whichever output plays it, and a leftover amplitudeDiff would
% only be added to an element nothing reads -- zeroing makes that explicit.
    if ~exist('ampDbm', 'var') || isempty(ampDbm); ampDbm = opts.amplitude; end
    amps = nan(1, 2); amps(k) = ampDbm;
    m = struct('ok', false, 'note', '', 'f', NaN, 'fErr', NaN, 'contrast', NaN, ...
               'r2', NaN, 'r2Ceiling', NaN, 'chi2', NaN, 'piNs', NaN, 'piErrNs', NaN, ...
               'T', NaN, 'residAutocorr', NaN, 'amps', amps, 'phaseDeg', 0, ...
               'seconds', NaN, 'tauMax', max(tauList), 'frequency', NaN);
    chans = {'MW', 'MW2'};

    b = AWGExpRabi;
    b.loadJsonParams('all', opts.jsonLocation);
    b.nChannels = 1;                    % BEFORE amplitude (set.nChannels re-enters set.amplitude)
    b.MWChannel = chans{k};             % 'MW' = AWG output 1, 'MW2' = output 2
    b.amplitudeDiff = [0 0];
    b.amplitude = ampDbm;
    b.channelPhase = [0 0];
    if ~isempty(opts.frequency); b.frequency = opts.frequency; end
    b.timeDelay = 0;
    b.repeats   = opts.repeats;
    b.averages  = opts.averages;
    b.isTracking = opts.isTracking;
    b.balancedSequence = 1;
    if ~isempty(opts.singleDetection)
        b.detectionDuration = opts.singleDetection;
        b.referenceDetectionDuration = b.detectionDuration;
    end                                 % else: the json's window, as the single-channel cell
    b.shouldAutosave = false;
    b.tau = tauList;
    b.useSequenceWalk    = opts.useSequenceWalk;
    b.passesPerRun       = opts.passesPerRun;
    b.walkWarmupTriggers = opts.walkWarmupTriggers;
    b.resyncEveryPass    = opts.resyncEveryPass;
    if ~isempty(opts.newOrderEveryAverage)
        b.newOrderEveryAverage = logical(opts.newOrderEveryAverage);
    end
    if ~isequal(b.activeChannels, k)
        m.note = sprintf('AWG output %s would play, not %d -- not run', mat2str(b.activeChannels), k);
        return
    end
    m = runAndFit(b, m, opts);
end

function m = runAndFit(b, m, opts)
% Run a configured AWGExpRabi (its output captured) and fit it with rabiFitRobust.
    txt = ''; tId = tic;
    try
        txt = evalc('b.run');
    catch err
        m.seconds = toc(tId); m.note = ['run failed: ' err.message];
        printIndented(txt); return
    end
    m.seconds = toc(tId);
    m.frequency = b.frequency;

    try
        o = rabiFitRobust(struct('exp', b, 'plotIt', false, 'verbose', false, ...
                                 'nFreq', opts.nFreq));
    catch err
        m.note = ['fit threw: ' err.message]; return
    end
    if ~isstruct(o) || ~isfield(o, 'ok') || ~o.ok
        m.note = 'fit failed';
        if isstruct(o) && isfield(o, 'why'); m.note = ['fit failed: ' o.why]; end
        return
    end
    m.f = o.f; m.fErr = o.fErr; m.contrast = o.contrast; m.r2 = o.r2;
    m.chi2 = o.chi2Reduced; m.residAutocorr = o.residAutocorr; m.T = o.T;
    m.piNs = 1e3*o.pulses(2); m.piErrNs = 1e3*o.pulsesErr(2);
    try
        m.r2Ceiling = 1 - sum(o.e.^2)/max(sum((o.y - mean(o.y)).^2), eps);
    catch
    end

    if m.contrast < opts.minContrast
        m.note = sprintf('contrast %.2f %% below %.2f %%', 100*m.contrast, 100*opts.minContrast);
        return
    end
    if m.fErr/max(m.f, eps) > opts.maxFErrFrac
        m.note = sprintf('frequency only %.0f %% determined', 100*(1 - m.fErr/max(m.f, eps)));
        return
    end
    m.ok = true; m.note = 'ok';
end

function [v, e] = combine(ms)
% Inverse-variance mean over the usable repeats of one channel.
    v = NaN; e = NaN;
    if isempty(ms); return; end
    ok = logical([ms.ok]);
    if ~any(ok); return; end
    f = [ms(ok).f]; s = [ms(ok).fErr];
    s(~isfinite(s) | s <= 0) = median(s(isfinite(s) & s > 0));
    if all(~isfinite(s)); v = mean(f); e = std(f)/sqrt(numel(f)); return; end
    w = 1./s.^2;
    v = sum(w.*f)/sum(w);
    e = 1/sqrt(sum(w));
    % If the repeats scatter by more than their error bars, the error bars are not the whole
    % story -- drift is -- and the scatter is the better estimate. Take the larger.
    if numel(f) > 1
        e = max(e, std(f)/sqrt(numel(f)));
    end
end

function spreadWarn(ms, name)
    if isempty(ms); return; end
    ok = logical([ms.ok]);
    if sum(ok) < 2; return; end
    f = [ms(ok).f]; s = [ms(ok).fErr];
    sep = abs(diff(f([1 end])));
    comb = hypot(s(1), s(end));
    if sep > 3*comb
        fprintf(['  *** %s repeated at %.4f then %.4f MHz, %.1f sigma apart. The drive or the\n' ...
                 '  fluorescence is drifting across the run, so the balance is being measured\n' ...
                 '  against a moving reference. Raise nBalanceReps -- the alternating order\n' ...
                 '  cancels a LINEAR drift, but only when there are enough pairs to average.\n'], ...
                name, f(1), f(end), sep/max(comb, eps));
    end
end

% ======================================================================
% REPORTING
% ======================================================================

function crossCheck(pf, O1, O2, delta)
% The phase scan measures Omega1+Omega2 independently of stage 1. They should agree.
% O2 here is channel 2 at the stage-2 amplitude, i.e. already scaled by 10^(delta/20).
    S = O1 + O2;
    d = pf.sumOmega - S;
    s = abs(d)/max(hypot(pf.sumOmegaErr, 0.01*S), eps);
    fprintf('\n  cross-check on the sum:\n');
    fprintf('    stage 1, one at a time   %.4f MHz\n', S);
    fprintf('    stage 2, from the curve  %.4f +/- %.4f MHz   (%+.4f, %.1f sigma)\n', ...
            pf.sumOmega, pf.sumOmegaErr, d, s);
    if s > 4
        fprintf(['    *** These should be the same number measured two ways. That they are not\n' ...
                 '    means one of the assumptions has broken: either the two channels are not\n' ...
                 '    adding linearly, or a channel measured alone is not the same channel when\n' ...
                 '    both play, or something drifted between the stages. This is\n' ...
                 '    the most informative disagreement on the page -- do not ignore it.\n']);
    end
    if pf.balanceResolved
        fprintf('    imbalance from the curve %.2f dB against %.3f dB one at a time\n', ...
                pf.imbalanceDb, abs(delta));
    else
        fprintf(['    the curve bounds the imbalance at under %.2f dB, consistent with the\n' ...
                 '    %.3f dB measured directly (it cannot do better near balance, by\n' ...
                 '    construction -- see phaseScanFit)\n'], pf.imbalanceDbBound, abs(delta));
    end
end

function reportClosure(R, Re, both, S, Se)
    fprintf('\n  CLOSURE TEST -- does superposition hold?\n');
    fprintf('    Omega(both)              %.4f MHz\n', both);
    fprintf('    Omega1 + Omega2 (ch2 +delta) %.4f +/- %.4f MHz\n', S, Se);
    fprintf('    R = ratio                %.4f +/- %.4f\n', R, Re);
    if abs(R - 1) < max(3*Re, 0.02)
        fprintf(['    CONSISTENT WITH 1. The two channels add exactly as a linear network must,\n' ...
                 '    so characterising them one at a time is valid and the amplitudeDiff above\n' ...
                 '    is the offset that really equalises them when both are driven. This is the\n' ...
                 '    measurement that settles the question rather than arguing it.\n']);
    elseif R < 1
        fprintf(['    *** R IS %.1f %% BELOW 1. The pair delivers less than the sum of its parts,\n' ...
                 '    and the phase is already at its fitted maximum, so a phase error is not the\n' ...
                 '    explanation. With iterateBalance (the default) each channel was MEASURED at\n' ...
                 '    the level it plays here, so single-channel compression is not it either:\n' ...
                 '    the two channels interact when both are on -- something they share\n' ...
                 '    compresses or heats, or one leaks back into the other''s chain. That CAN\n' ...
                 '    unbalance them, so run calibMwQuadrature (level4.m 4.2b): it measures the\n' ...
                 '    balance with both channels on. (With iterateBalance = false, channel 2 at\n' ...
                 '    operating + delta was extrapolated, and its own compression is the first\n' ...
                 '    suspect -- rerun with iterateBalance = true.)\n'], 100*(1 - R));
    else
        fprintf(['    *** R IS %.1f %% ABOVE 1, which no linear network can do -- the sum is an\n' ...
                 '    upper bound. Something is inconsistent between the stages rather than\n' ...
                 '    interesting: most likely the drive drifted upward between the stages, or\n' ...
                 '    (offAmplitude mode only) the held-off channel was still contributing.\n'], ...
                100*(R - 1));
    end
end

function reportAntiphase(m, tauList, S, delta, deltaErr, out)                %#ok<INUSD>
    span = max(tauList) - min(tauList);
    % The smallest frequency this sweep can establish: about a half cycle over the sweep.
    % But the Rabi envelope also decays, and a sweep longer than the coherence is not
    % measuring anything, so the effective length is the shorter of the two.
    Tdec = NaN;
    if isstruct(out.closure) && isfield(out.closure, 'T'); Tdec = out.closure.T; end
    effSpan = span;
    if isfinite(Tdec) && 2*Tdec < span; effSpan = 2*Tdec; end
    fMin = 1/(2*max(effSpan, eps));

    fprintf('\n');
    if isfinite(Tdec)
        fprintf('  Rabi decay measured in stage 3: %.2f us.\n', Tdec);
        if 2*Tdec < span
            fprintf(['  That is shorter than the %.1f us sweep, so the useful length is about\n' ...
                     '  2T = %.1f us and the sweep cannot resolve below %.3f MHz. Lengthening\n' ...
                     '  antiphaseTauMax past this buys nothing -- there is no signal out there.\n'], ...
                    span, effSpan, fMin);
        end
    end

    if m.ok
        implied = 20*log10((S + m.f)/max(S - m.f, eps));
        fprintf('  residual Omega at antiphase  %.4f +/- %.4f MHz\n', m.f, m.fErr);
        fprintf('  implies an imbalance of      %.3f dB\n', implied);
        fprintf('  measured directly in stage 1 %.3f dB +/- %.3f\n', abs(delta), deltaErr);
        phiErr = NaN;
        if isstruct(out) && isfield(out, 'fit') && isstruct(out.fit) && isfield(out.fit, 'phiMaxErr')
            phiErr = out.fit.phiMaxErr;
        end
        if abs(implied - abs(delta)) > max(4*deltaErr, 0.15) && isfinite(phiErr)
            fprintf(['    *** The two disagree. The null is the more direct measurement in\n' ...
                     '    principle -- both channels driven, nothing extrapolated -- but it is\n' ...
                     '    also the one that a residual PHASE error inflates: if the phase is off\n' ...
                     '    the null by dphi, the residual picks up 2*Omega*sin(dphi/2) whether or\n' ...
                     '    not the amplitudes match. With the phase determined to %.1f deg that\n' ...
                     '    floor is %.3f MHz. Compare it against the residual above before\n' ...
                     '    believing either number.\n'], ...
                    phiErr, S*abs(sind(phiErr/2)));
        else
            fprintf(['    The null and the one-at-a-time measurement agree. Two different\n' ...
                     '    measurements, one with a single channel driven and one with both, give\n' ...
                     '    the same imbalance -- which is the same conclusion as the closure test\n' ...
                     '    reached from the other direction.\n']);
        end
    else
        bound = 20*log10((S + fMin)/max(S - fMin, eps));
        fprintf(['  NO OSCILLATION RESOLVED AT ANTIPHASE, and that is the good outcome. The\n' ...
                 '  residual drive is below what a %.1f us sweep can see (%.3f MHz), so\n' ...
                 '        |Omega1 - Omega2| < %.3f MHz   ->   imbalance under %.2f dB\n' ...
                 '  which is a bound, not a value. Stage 1 measured %.3f dB directly, and the\n' ...
                 '  two are %s.\n'], ...
                effSpan, fMin, fMin, bound, abs(delta), ...
                BooleanHelper.ifTrueElse(abs(delta) <= bound, 'consistent', ...
                    'NOT consistent -- stage 1 claims more imbalance than the null allows'));
        fprintf('  (fit note: %s)\n', m.note);
    end
end

function printHandover(out, opts)
    fprintf('\n=== what level4.m cell 4.3 writes to expProperties.json (expPropsWrite) ===\n');
    if ~out.ok
        fprintf(['  NOT SETTLED -- read the stages above. Cell 4.3 writes nothing while this is\n' ...
                 '  the case; fix the failing stage and re-run.\n']);
        return
    end
    fprintf('  amplitudeDiff  [0, %+.3f]   dB, measured at %g dBm on channel 1\n', ...
            out.delta, opts.amplitude);
    if out.balanceConverged
        fprintf(['                 channel 2 measured AT %.3f dBm, the level it plays (%d ' ...
                 'balance round(s))\n'], opts.amplitude + out.delta, numel(out.balanceRounds) - 1);
    end
    fprintf('  channelPhase   [0, %.2f]    degrees\n', out.phiMax);
    fprintf('  timeDelay      0            REQUIRED -- see below\n');
    fprintf(['\n  timeDelay MUST go to zero. It is the other phase mechanism and the two are\n' ...
             '  ADDITIVE inside the waveform builder, so leaving both non-zero means the real\n' ...
             '  phase is their sum and nothing in the json records that. The stored 327 ps is\n' ...
             '  within 6 %% of a full period at 2870 MHz anyway, so dropping it costs nothing.\n']);
    fprintf(['\n  The balance was measured AT %g dBm and the amplifiers are per channel, so it\n' ...
             '  is a level-dependent number. If you later run at a materially different\n' ...
             '  amplitude, re-run stage 1 alone -- level4.m cell 4.7 does that and writes the\n' ...
             '  new amplitudeDiff. The PHASE does not need redoing -- it is set by path\n' ...
             '  lengths, not by level.\n'], opts.amplitude);
    if isstruct(out.closure) && isfield(out.closure, 'piNs') && isfinite(out.closure.piNs)
        fprintf('\n  pi pulse with both channels: %.1f +/- %.1f ns.\n', ...
                out.closure.piNs, out.closure.piErrNs);
    end
end

% ======================================================================
% PLOT
% ======================================================================

function drawIt(out, opts)
    C = struct('surface', [252 252 251]/255, 'ink', [11 11 11]/255, ...
               'ink2', [82 81 78]/255, 'data', [42 120 214]/255, ...
               'fit', [235 104 52]/255, 'drop', [140 138 133]/255);
    P = out.phaseScan; ok = logical([P.ok]);
    ph = [P(ok).phaseDeg]; w = [P(ok).f]; we = [P(ok).fErr];

    fig = figure('Name', 'calibMwPair', 'NumberTitle', 'off', 'Color', C.surface);
    ax = axes(fig); hold(ax, 'on');
    hE = errorbar(ax, ph, w, we);
    set(hE, 'LineStyle', 'none', 'Color', C.data, 'LineWidth', 0.75);
    try; set(hE, 'CapSize', 4); catch; end
    try; hE.Color(4) = 0.35; catch; end
    plot(ax, ph, w, 'o', 'MarkerSize', 6, 'LineWidth', 1, 'Color', C.data, ...
         'MarkerFaceColor', C.surface);
    if any(~ok)
        B = P(~ok); bp = [B.phaseDeg]; bw = [B.f]; k = isfinite(bw);
        if any(k)
            plot(ax, bp(k), bw(k), 'x', 'MarkerSize', 9, 'LineWidth', 1.2, 'Color', C.drop);
        end
    end
    if isstruct(out.fit) && isfield(out.fit, 'model') && isa(out.fit.model, 'function_handle')
        pp = linspace(0, 360, 720);
        % The fit is in Omega^2; the measurement is Omega, so plot the square root of the
        % model rather than fitting the square root.
        plot(ax, pp, sqrt(max(out.fit.model(pp), 0)), '-', 'LineWidth', 2, 'Color', C.fit);
        yl = get(ax, 'YLim');
        plot(ax, [out.phiMax out.phiMax], yl, '-',  'LineWidth', 1, 'Color', C.fit);
        plot(ax, [out.phiNull out.phiNull], yl, '--', 'LineWidth', 1, 'Color', C.fit);
        set(ax, 'YLim', yl);
    end
    % The single-channel measurements, as a horizontal reference: the curve must peak at
    % their sum and null at their difference, and seeing all three on one axes is the whole
    % argument of this calibration in one picture.
    if isfinite(out.omega1) && isfinite(out.omega2)
        S = out.sumBalanced;                        % at the phase-scan amplitudes
        plot(ax, [0 360], [S S], ':', 'LineWidth', 1.5, 'Color', C.ink2);
        text(ax, 5, S, sprintf('  \\Omega_1 + \\Omega_2 = %.3f', S), 'Color', C.ink2, ...
             'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
    grid(ax, 'on');
    set(ax, 'Color', C.surface, 'XColor', C.ink2, 'YColor', C.ink2, 'Box', 'off', ...
            'GridLineStyle', '-', 'Layer', 'top', 'TickDir', 'out', 'FontSize', 10);
    try; set(ax, 'GridAlpha', 0.18); catch; end
    xlim(ax, [0 360]); set(ax, 'XTick', 0:45:360);
    xlabel(ax, 'channel 2 phase  (deg)', 'Color', C.ink2);
    ylabel(ax, '\Omega_{Rabi}  (MHz)', 'Color', C.ink2);
    ttl = sprintf('\\phi_{max} = %.2f \\pm %.2f deg    amplitudeDiff = [0 %+.3f] dB', ...
                  out.phiMax, out.fit.phiMaxErr, out.delta);
    title(ax, ttl, 'Color', C.ink, 'FontWeight', 'normal');
    if isfield(out, 'closureRatio') && isfinite(out.closureRatio)
        sub = sprintf('R = \\Omega_{both}/(\\Omega_1+\\Omega_2) = %.4f \\pm %.4f    at %g dBm', ...
                      out.closureRatio, out.closureRatioErr, opts.amplitude);
        try
            subtitle(ax, sub, 'Color', C.ink2);
        catch
            title(ax, {ttl, sub}, 'Color', C.ink, 'FontWeight', 'normal');
        end
    end
end

% ======================================================================

function printIndented(txt)
    if isempty(txt); return; end
    L = strsplit(txt, newline);
    for i = 1:numel(L)
        if ~isempty(strtrim(L{i})); fprintf('      | %s\n', L{i}); end
    end
end

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d); given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('%s: unknown option "%s". Valid: %s.', who, bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || isempty(o.(known{i}))
            o.(known{i}) = d.(known{i});
        end
    end
end
