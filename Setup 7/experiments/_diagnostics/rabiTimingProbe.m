function out = rabiTimingProbe(opts)
%RABITIMINGPROBE Split a Rabi run into fixed cost, per-pass cost and photon time.
%
%   out = rabiTimingProbe
%   out = rabiTimingProbe(struct('nTau', 20, 'repeats', 3000))
%
% Four short runs. What comes out: the cost that is paid once per run no matter what,
% the cost of one pass over all the taus, the share of that pass which is actually
% photons, and what the autosave costs.
%
% ===================================================================
% WHY THIS DIFFERENCES OVER AVERAGES, AND NOT OVER POINTS
% ===================================================================
%
% The first version of this probe compared runs with N and 2N taus and with R and R/2
% repeats, and solved for a per-point overhead. On real data it returned a NEGATIVE
% overhead, which is not scatter -- it is the model being wrong. Measured:
%
%       20 taus x 3000 repeats : 7.9 s
%       40 taus x 3000 repeats : 7.8 s
%       20 taus x 1500 repeats : 7.5 s
%
% Doubling the number of points changed the wall clock by 0.1 s and halving the repeats
% changed it by 0.4 s. Nothing in that experiment scaled with the amount of measuring
% being done, because a FIXED cost of roughly seven and a half seconds swamped all of
% it. Differencing over N or R cannot separate a fixed cost from anything -- it appears
% identically in both runs and cancels out of the wrong side of the subtraction.
%
% Averages are the one axis that multiplies ONLY the measurement. One run at A = 1 and
% one at A = 3, everything else identical:
%
%       t(A) = fixed + A * onePass
%
% so onePass = (t3 - t1)/2 and fixed = t1 - onePass, with no modelling assumptions at
% all. Then a third run at half the repeats splits onePass into the part that scales
% with repeats (photons) and the part that does not (per-point instrument traffic).
%
% The lesson generalises: never difference across an axis that a suspected fixed cost
% also lives on.
%
% ===================================================================
% WHAT THE PRINTED "each average taking X seconds" IS, AND IS NOT
% ===================================================================
%
% It is a PREDICTION, not a measurement:
%
%       averageTime = obj.repeats * pg.sequenceDuration() * nRepeatsRuns
%
% -- the pulse generator's nominal sequence duration times the number of repeats. On
% your run it read 3.3 s for 20 taus and 6.6 s for 40, i.e. a sequence duration of
% about 55 us per repeat. Compare that with the 7.9 s and 7.8 s actually taken and the
% two are not consistent with each other: the wall clock did not double when the
% prediction did. So one of the two is wrong, and this probe measures which by timing
% the same run at two different average counts.
%
% If onePass comes out well BELOW the prediction, the sequence is not really taking 55
% us per repeat and the counter is streaming the sweep rather than stopping at each
% point. If it comes out ABOVE, the overhead is real and sits between the points.
%
% OPTIONS
%   nTau, repeats, averagesHigh   the two average counts are 1 and averagesHigh
%                                 (defaults 20, 3000, 4)
%   repeatFactor  the second repeat count is this times the first (default 8). It has
%                 to be a LARGE ratio: the first run of this probe used R and R/2 and
%                 the difference was buried, because the measurement turned out to be
%                 a tenth of a second while the autosave was six seconds.
%   amplitude, frequency, tauMax, nChannels, detectionMHz, jsonLocation, which
%   testAutosave   also run once with shouldAutosave = false, to price it (default true)
%
% Changes no calibration and saves nothing of its own.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('nTau', 20, 'repeats', 3000, 'averagesHigh', 4, ...
        'repeatFactor', 8, ...
        'amplitude', -16, 'frequency', [], 'tauMax', 0.5, 'nChannels', 2, ...
        'detectionMHz', 0.3, 'which', 'AWGExpRabi', 'testAutosave', true, ...
        'jsonLocation', ''));

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

    N = opts.nTau; R = opts.repeats; Ah = opts.averagesHigh;
    Rbig = round(R*opts.repeatFactor);
    taus = round(linspace(0.02, opts.tauMax, N)*1000)/1000;

    fprintf('\n=== rabiTimingProbe (%s) ===\n', opts.which);
    fprintf('%d taus. Timing runs have AUTOSAVE OFF -- measured at ~6 s per run, it is\n', N);
    fprintf('both the largest single cost and the largest source of scatter, so it is\n');
    fprintf('priced once at the end instead of being carried through every subtraction.\n');

    % 1 and 2 differ ONLY in averages, so every fixed cost cancels in the difference.
    t1 = timeRun(opts, taus, R, 1,  false);
    fprintf('  A = 1, R = %6d, autosave off : %6.2f s\n', R, t1);
    t2 = timeRun(opts, taus, R, Ah, false);
    fprintf('  A = %d, R = %6d, autosave off : %6.2f s\n', Ah, R, t2);

    onePass   = (t2 - t1)/max(Ah - 1, 1);
    fixedCost = t1 - onePass;

    % 3 differs from 2 only in repeats, and by a LARGE factor so the photon term is
    % bigger than the run-to-run scatter.
    t3 = timeRun(opts, taus, Rbig, Ah, false);
    fprintf('  A = %d, R = %6d, autosave off : %6.2f s\n', Ah, Rbig, t3);
    onePassBig = (t3 - fixedCost)/Ah;

    p = (onePassBig - onePass) / max(N*(Rbig - R), 1);     % s per repeat
    c = onePass/N - R*p;                                   % s per point

    % 4 prices the autosave.
    t4 = NaN;
    if opts.testAutosave
        t4 = timeRun(opts, taus, R, 1, true);
        fprintf('  A = 1, R = %6d, autosave ON  : %6.2f s\n', R, t4);
    end

    out = struct('t1', t1, 't2', t2, 't3', t3, 't4', t4, 'onePass', onePass, ...
                 'onePassBig', onePassBig, 'fixedCost', fixedCost, 'perRepeat', p, ...
                 'perPoint', c, 'nTau', N, 'repeats', R, 'repeatsBig', Rbig);

    fprintf('\n--- breakdown ---\n');
    fprintf('  fixed, once per run        %6.2f s   (AWG load, prepare, pulse generator)\n', fixedCost);
    fprintf('  one pass over %3d taus     %6.2f s at R = %d\n', N, onePass, R);
    fprintf('                             %6.2f s at R = %d\n', onePassBig, Rbig);
    if isfinite(t4)
        out.autosaveCost = t4 - t1;
        fprintf('  autosave                   %6.2f s per run\n', out.autosaveCost);
    end

    % A pass smaller than the scatter is not a measurement of a pass. Say so rather
    % than dividing by it and printing a ratio with fourteen digits.
    scatter = 0.3;
    if onePass < scatter
        fprintf(['\n  ONE PASS IS BELOW THE RUN-TO-RUN SCATTER (%.2f s < %.1f s). Adding\n' ...
                 '  averages costs nothing measurable, so the measurement is NOT what this\n' ...
                 '  run spends its time on -- the fixed cost is. Do not optimise the sweep.\n'], ...
                onePass, scatter);
        out.passMeasurable = false;
    else
        out.passMeasurable = true;
        fprintf('    of which photons         %6.2f s   (%.1f us per repeat)\n', N*R*p, 1e6*p);
        fprintf('    of which traffic         %6.2f s   (%.1f ms per tau)\n', N*c, 1e3*c);
        if c > R*p && c > 0
            fprintf(['\n  INSTRUMENT TRAFFIC DOMINATES a pass by %.1fx -- the single-sequence\n' ...
                     '  jumpToStep scheme is worth writing.\n'], c/max(R*p, eps));
        else
            fprintf(['\n  PHOTONS DOMINATE a pass by %.1fx. Shorten the sequence rather than the\n' ...
                     '  AWG traffic: laserInitializationDuration is paid twice per repeat under\n' ...
                     '  balancedSequence, and calibDetectionWindow measured how short it can be.\n'], ...
                    (R*p)/max(c, eps));
        end
    end

    % ---- does `repeats` do what it says? ---------------------------------
    % The experiment predicts repeats * pg.sequenceDuration() * nParams per average. If
    % raising the repeat count by repeatFactor does not raise the pass time by about the
    % same factor, then either that prediction is wrong or the counter is not
    % integrating for as many repeats as asked -- and that is an SNR question, not a
    % speed one, so it matters more than anything else on this page.
    fprintf('\n--- does `repeats` cost what it should? ---\n');
    ratio = onePassBig / max(onePass, eps);
    fprintf('  repeats x%.0f changed a pass by x%.2f (should be about x%.0f)\n', ...
            opts.repeatFactor, ratio, opts.repeatFactor);
    if onePass >= scatter && (ratio < 0.5*opts.repeatFactor)
        fprintf(['  *** IT DOES NOT SCALE. You are probably not integrating the repeats you\n' ...
                 '      asked for. Check pg.sequenceDuration() against a scope trace, and\n' ...
                 '      check the SNR against sqrt(repeats) -- if the time does not scale but\n' ...
                 '      the noise does, the time is being spent elsewhere; if neither scales,\n' ...
                 '      the repeats are not happening.\n']);
    end

    fprintf('\n--- for a 97-tau, 3-average run, on these numbers ---\n');
    nF = 97; aF = 3;
    loadScale = fixedCost;      % the sequence load grows with the tau count
    now_ = loadScale + aF*max(onePass, 0)*nF/N;
    fprintf('  as it stands, autosave off       %6.0f s\n', now_);
    if isfinite(t4)
        fprintf('  with autosave on                 %6.0f s\n', now_ + out.autosaveCost);
    end
    fprintf(['\n  Note the fixed cost is mostly the AWG SEQUENCE LOAD, which grows with the\n' ...
             '  number of taus and is paid again whenever frequency, amplitude or the tau\n' ...
             '  list changes. AWGExpRabiCont holds one waveform per channel and no per-tau\n' ...
             '  sequences at all, so for a scan that steps FREQUENCY it avoids this entirely.\n']);
end

function t = timeRun(opts, taus, reps, avgs, autosave)
    b = makeExp(opts, taus, reps, avgs);
    b.shouldAutosave = autosave;
    % ID FORM, tId = tic ... toc(tId), NOT bare tic/toc.
    %
    % A bare tic sets a single global stopwatch, so ANY tic called inside b.run resets
    % it and the following toc measures from that internal tic instead of from here.
    % That is not hypothetical: this function reported 1.9 s for a run whose own log
    % showed 18.5 s of sequence settle and 16.4 s to arm, and the bench experience was
    % two minutes. Every number rabiTimingProbe produced was wrong for the same reason.
    % rabiStageTimer happened to use the id form and its 264 s is the one measurement of
    % the three that can be trusted.
    tId = tic; b.run; t = toc(tId);
end

function b = makeExp(opts, taus, reps, avgs)
    switch lower(opts.which)
        case 'awgexprabicont'; b = AWGExpRabiCont;
        otherwise;             b = AWGExpRabi;
    end
    b.loadJsonParams('all', opts.jsonLocation);
    if isprop(b, 'useSequenceWalk'); b.useSequenceWalk = false; end   % times the classic path
    b.nChannels = opts.nChannels;                     % BEFORE amplitude, always
    if opts.nChannels == 2
        b.amplitude = opts.amplitude + [0 0];
    else
        b.amplitude = opts.amplitude;
    end
    if ~isempty(opts.frequency); b.frequency = opts.frequency; end
    b.repeats  = reps;
    b.averages = avgs;
    b.isTracking = 0;
    b.balancedSequence = 1;
    b.detectionDuration = opts.detectionMHz;
    b.referenceDetectionDuration = b.detectionDuration;
    b.tau = taus;
end

function o = withDefaults(o, d)
    if isempty(o); o = struct(); end
    known = fieldnames(d);
    given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('rabiTimingProbe: unknown option "%s". Valid: %s.', bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || isempty(o.(known{i}))
            o.(known{i}) = d.(known{i});
        end
    end
end
