function out = calibMwQuadrature(opts)
%CALIBMWQUADRATURE Is the pair balanced while BOTH channels are playing?  (level 4, cell 4.2b)
%
%   q = calibMwQuadrature(struct('amplitude', -32, 'delta', out.delta, 'phiMax', out.phiMax))
%   q = calibMwQuadrature                         % everything from the json
%   q = calibMwQuadrature(struct('replay', q))    % re-analyse, no hardware
%   q3 = calibMwQuadrature(struct(..., 'quadAngle', 1.5));  qc = mwQuadratureCombine(q, q3)
%
% calibMwPair measures each channel ALONE, at the level it will play, and checks with the
% closure test that the two add as a linear network must. What no single-channel
% measurement can show is the balance with both channels on -- and that is the only
% balance an experiment ever uses. This measures it directly, signed, in dB.
%
% HOW (see AWGExpMwQuadrature for the physics): a pulse with the two channels 90 degrees
% apart has its rotation axis at atan(Omega2/Omega1), i.e. at 45 degrees exactly when the
% two deliver the same field at the NVs. An in-phase pi/2 about 135 / -45 degrees reads the
% axis tilt as population. Played with channel 2 leading and with channel 1 leading, the
% difference Q is first order in the imbalance and blind, to first order, to a
% channelPhase error. Every configuration is played at three channel-2 levels around
% delta (offsetsDb), so the slope and its sign are measured; the zero crossing of the
% fitted line is the both-on balance:
%
%       deltaBoth = delta + d*
%
% Simulated (14N triplet, 1 MHz inhomogeneous, +/-10% opposite field gradients from the two
% PCBs): slope ~0.056 of the full pi contrast per dB at 1.7 MHz per channel, ~0.11 at
% 8 MHz; zero within 0.05 dB; unmoved by a 3 degree phase error, the gradient, or +/-10%
% pulse-time errors.
%
% CHANNEL SKEW (2026-11-26). What the simulation above left out: a delay between the two
% channels' pulse envelopes at the NVs. For that long one channel plays alone at each pulse
% edge, and Q reads it as an imbalance -- about 0.08 dB per ns per MHz of per-channel drive
% (0.55 dB per ns at 6.7 MHz), the same for a gap between the pulses or none. Nothing
% aligns the envelopes (timeDelay = 0; channelPhase only fixes the carrier, modulo 0.35 ns),
% and 1-2 ns is 20-40 cm of cable. With the quadrature pulse at 3pi/2 (quadAngle 1.5) a
% real imbalance reads the same and the skew bias changes sign, so the two runs agree
% (within 0.07 dB, simulated) only when there is no skew. mwQuadratureCombine(q, q3) uses
% that as the skew test, and gives a balance only when it passes: the mean of the two is
% NOT a reliable skew correction (how far the bias flips depends on the pulse length).
% level4.m 4.2b runs both. A single run's verdict below cannot tell skew from imbalance.
%
% VERDICT (out.verdict):
%   'agree'         |d*| <= 3 sigma + tolDb, sigma = both errors combined: the
%                   single-channel balance holds when both
%                   play. level4.m 4.3 writes.
%   'disagree'      it does not. The pair interacts when both are on. Re-run the phase scan
%                   at deltaBoth (the order rule: amplitude first, then phase -- one round
%                   is enough, the cross-talk between the two is second order), then this
%                   again. level4.m prints the two calls.
%   'inconclusive'  the slope or the zero crossing is not measured well enough (sigma >
%                   maxErrDb). More averages, or a stronger drive.
%
% OPTIONS (all optional)
%   amplitude      dBm on channel 1 ([] = the json's amplitude)
%   delta          dB on channel 2 relative to channel 1 ([] = json amplitudeDiff(2))
%   deltaErr       the error on delta from calibMwPair (out.deltaErr). The comparison uses both
%                  errors: in compression stage 1's error grows (it is divided by channel 2's
%                  slope), and ignoring it calls noise a disagreement. Default 0.
%   phiMax         deg, channel 2 phase for maximum drive ([] = json channelPhase(2))
%   piTime         us, the PAIR's pi time at this amplitude ([] = the json's piTime, with a
%                  warning if the json amplitude differs). quadT2 = piTime/2 (in-phase pi/2),
%                  quadT1 = sqrt(2)*quadT2 (pi/2 for |Omega1 + i Omega2|). The result is
%                  insensitive to +/-10% errors here.
%   quadAngle      the quadrature pulse's rotation in units of pi: 0.5 (default) or 1.5,
%                  i.e. quadT1 = sqrt(2)*quadT2 or 3*sqrt(2)*quadT2. Same levels either way.
%                  1.5 has the opposite channel-skew bias (see CHANNEL SKEW). It needs a
%                  drive well above the 2.16 MHz hyperfine: simulated at 1.7 MHz per channel
%                  with pulses 10% long it lost the slope entirely; at 3.5 and 6.7 MHz fine.
%   offsetsDb      channel-2 levels around delta (default [-1 0 1])
%   repeats        passes per average (default 10000; walk: repeats = passes)
%   averages       default 20
%   detectionDuration  us ([] = the json's, the pulsed window level 2 chose)
%   maxDigitizerSegments  per pg.run (default 20000)
%   useSequenceWalk, newOrderEveryAverage, walkWarmupTriggers, resyncEveryPass
%                  (defaults true, true, 1, true)
%   tolDb          added to 3 sigma in the agreement test (default 0.1 dB). It covers this
%                  test's own small biases, measured in simulation: +0.03-0.05 dB from the 14N
%                  hyperfine at weak drive, up to +0.09 dB in heavy compression (a straight
%                  line through a curved response). Without it, 8% of balanced pairs in
%                  compression were called 'disagree'.
%   maxErrDb       above this error on d* the answer is 'inconclusive' (default 0.3)
%   frequency, jsonLocation, isTracking, plotIt, replay

    if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
    opts = withDefaults(opts, struct( ...
        'amplitude', [], 'delta', [], 'deltaErr', 0, 'phiMax', [], 'piTime', [], 'quadAngle', 0.5, ...
        'offsetsDb', [-1 0 1], ...
        'repeats', 10000, 'averages', 20, 'detectionDuration', [], ...
        'maxDigitizerSegments', 20000, 'useSequenceWalk', true, 'newOrderEveryAverage', true, ...
        'walkWarmupTriggers', 1, 'resyncEveryPass', true, 'tolDb', 0.1, 'maxErrDb', 0.3, ...
        'frequency', [], 'jsonLocation', '', 'isTracking', 0, 'plotIt', true, 'replay', []), ...
        'calibMwQuadrature');

    if ~isempty(opts.replay)
        R0 = opts.replay;
        if ~isstruct(R0) || ~isfield(R0, 'raw')
            error('calibMwQuadrature: replay must be an earlier output (with a raw field).');
        end
        fprintf('\n=== calibMwQuadrature: REPLAY, no hardware ===\n');
        out = R0;
        out.opts.tolDb = opts.tolDb;  out.opts.maxErrDb = opts.maxErrDb;
        out = analyse(out);
        report(out);
        if opts.plotIt; drawIt(out); end
        return
    end

    if ~(isscalar(opts.quadAngle) && any(abs(opts.quadAngle - [0.5 1.5]) < 1e-9))
        error('calibMwQuadrature: quadAngle must be 0.5 (pi/2) or 1.5 (3pi/2).');
    end
    if isempty(opts.jsonLocation); opts.jsonLocation = setup7Paths().json; end
    J = jsondecode(fileread(opts.jsonLocation));
    amp   = pickVal(opts.amplitude, J, 'amplitude', [], 1);
    delta = pickVal(opts.delta, J, 'amplitudeDiff', 0, 2);
    phi   = pickVal(opts.phiMax, J, 'channelPhase', 0, 2);
    if isempty(amp)
        error('calibMwQuadrature: no amplitude given and none in the json.');
    end
    piT = opts.piTime;
    if isempty(piT)
        piT = pickVal([], J, 'piTime', [], 1);
        if isempty(piT)
            error('calibMwQuadrature: no piTime given and none in the json.');
        end
        if isfield(J, 'amplitude') && abs(J.amplitude(1) - amp) > 0.01
            fprintf(['NOTE the json piTime (%.4f us) belongs to %g dBm, not %g. The test is ' ...
                     'insensitive to +/-10%% pulse errors, but pass piTime\n     from the ' ...
                     'stage-3 Rabi (calibMwPair out.closure.piNs/1e3) if they differ more.\n'], ...
                    piT, J.amplitude(1), amp);
        end
    end
    T2 = round(piT / 2 * 1e4) / 1e4;
    T1 = round(sqrt(2) * T2 * (opts.quadAngle / 0.5) * 1e4) / 1e4;
    angleName = ternary(opts.quadAngle > 1, '3pi/2', 'pi/2');

    fprintf('\n=== calibMwQuadrature: the balance with BOTH channels playing ===\n');
    fprintf(['channel 1 %g dBm, channel 2 %g %+.3f dB = %.3f dBm, channel 2 phase %.2f deg\n' ...
             'quadrature pulse %.1f ns (%s), in-phase pi/2 %.1f ns (pair pi %.1f ns); channel-2 ' ...
             'offsets %s dB\n'], amp, amp, delta, amp + delta, phi, 1e3*T1, angleName, 1e3*T2, ...
            1e3*piT, mat2str(opts.offsetsDb));

    b = AWGExpMwQuadrature;
    b.loadJsonParams('all', opts.jsonLocation);
    b.nChannels = 2;                    % before amplitude, always
    b.amplitudeDiff = [0 0];            % explicit vector below: no double-adding
    b.amplitude = [amp, amp + delta];
    if ~isempty(opts.frequency); b.frequency = opts.frequency; end
    b.timeDelay = 0;
    b.channelPhase = [0 phi];
    b.quadT2 = T2;
    b.quadT1 = T1;
    b.quadOffsetsDb = opts.offsetsDb;
    b.repeats = opts.repeats;
    b.averages = opts.averages;
    b.isTracking = opts.isTracking;
    b.balancedSequence = 1;
    if ~isempty(opts.detectionDuration)
        b.detectionDuration = opts.detectionDuration;
        b.referenceDetectionDuration = opts.detectionDuration;
    end
    b.shouldAutosave = false;
    b.useSequenceWalk = opts.useSequenceWalk;
    b.maxDigitizerSegments = opts.maxDigitizerSegments;
    b.walkWarmupTriggers = opts.walkWarmupTriggers;
    b.resyncEveryPass = opts.resyncEveryPass;
    b.newOrderEveryAverage = logical(opts.newOrderEveryAverage);

    nCfg = numel(b.configs);
    fprintf(['%d configurations, %d passes x %d averages = %d repeats each. Expect roughly %.0f ' ...
             'minutes.\n'], nCfg, opts.repeats, opts.averages, opts.repeats * opts.averages, ...
            estimateMinutes(b, opts, nCfg));

    tId = tic;
    txt = '';
    try
        txt = evalc('b.run');
    catch err
        printIndented(txt);
        error('calibMwQuadrature: the run failed: %s', err.message);
    end
    fprintf('run took %.1f min\n', toc(tId) / 60);

    raw = pullRaw(b);
    out = struct('opts', opts, 'raw', raw, 'amplitude', amp, 'delta', delta, ...
                 'deltaErr', nanToZero(opts.deltaErr), 'phiMax', phi, ...
                 'piTime', piT, 'quadAngle', opts.quadAngle, 'quadT1', T1, 'quadT2', T2, ...
                 'configs', b.configs, ...
                 'frequency', b.frequency, 'when', datestr(now, 31)); %#ok<TNOW1,DATST>
    out = analyse(out);
    report(out);
    if opts.plotIt; drawIt(out); end
end

% ======================================================================

function raw = pullRaw(b)
    % Per-average sig/ref ratio per configuration, as rabiFitRobust rebuilds it.
    sig = b.signal;  ste = b.sterr;
    filled = squeeze(all(all(sig ~= 0, 1), 2));
    nA = find(~filled(:)', 1) - 1;
    if isempty(nA); nA = size(sig, 3); end
    nA = max(nA, 1);
    s1 = reshape(sig(1, :, 1:nA), [], nA);
    s2 = reshape(sig(2, :, 1:nA), [], nA);
    R = s1 ./ s2;
    E = [];
    try
        E = b.signalParam.sterr(:);
    catch
    end
    raw = struct('R', R, 'E', E, 's1', s1, 's2', s2, ...
                 'd1', reshape(ste(1, :, 1:nA), [], nA), 'd2', reshape(ste(2, :, 1:nA), [], nA));
end

function out = analyse(out)
    f = mwQuadratureFit(out.configs, out.raw.R, out.raw.E);
    out.fit = f;
    out.verdict = 'inconclusive';
    out.extrapolated = false;
    out.dStar = NaN;  out.dStarErr = NaN;  out.deltaBoth = NaN;  out.deltaBothErr = NaN;
    if f.ok
        out.dStar = f.dStar;  out.dStarErr = f.dStarErr;
        out.deltaBoth = out.delta + f.dStar;  out.deltaBothErr = f.dStarErr;
        out.extrapolated = ~f.inside;
        % Both numbers measure the same thing when they agree; the inverse-variance mean is
        % reported (not written -- 4.3 writes calibMwPair's).
        e1 = getDeltaErr(out);
        if e1 > 0
            w1 = 1 / e1^2;  w2 = 1 / f.dStarErr^2;
            out.deltaCombined = (w1 * out.delta + w2 * out.deltaBoth) / (w1 + w2);
            out.deltaCombinedErr = 1 / sqrt(w1 + w2);
        else
            out.deltaCombined = NaN;  out.deltaCombinedErr = NaN;
        end
        if f.dStarErr > out.opts.maxErrDb
            out.verdict = 'inconclusive';
        elseif abs(f.dStar) <= 3 * hypot(f.dStarErr, getDeltaErr(out)) + out.opts.tolDb
            out.verdict = 'agree';
        else
            % Outside the offsets played it is still a clear disagreement -- only the value
            % is an extrapolation (and in compression a straight line undershoots it).
            out.verdict = 'disagree';
        end
    end
end

function report(out)
    f = out.fit;
    fprintf('\n  offset/dB      Q (x1e-4)           Q / contrast        P / contrast   (P: phase-sensitive, info)\n');
    for k = 1:numel(f.offsets)
        fprintf('  %+6.2f    %+8.2f +/- %-6.2f   %+8.4f +/- %-7.4f   %+8.4f\n', f.offsets(k), ...
                1e4 * f.Q(k), 1e4 * f.Qerr(k), f.Qnorm(k), f.QnormErr(k), f.Pnorm(k));
    end
    fprintf('  contrast (bright - dark) %.4f +/- %.4f\n', f.contrast, f.contrastErr);
    if isfield(f, 'slope')
        fprintf('  line: Q = %+.3g + (%+.3g +/- %.2g) * offset;  slope %.4f of the contrast per dB', ...
                f.intercept, f.slope, f.slopeErr, f.slopeNorm);
        if isfinite(f.chi2nu); fprintf(', chi2/dof %.2f', f.chi2nu); end
        fprintf('\n');
    end
    fprintf('\n==================== BOTH CHANNELS ON ====================\n');
    if f.ok
        fprintf('  zero crossing  d* = %+.3f +/- %.3f dB  (relative to the delta played)\n', ...
                out.dStar, out.dStarErr);
        fprintf('  balance with both on:  delta = %+.3f +/- %.3f dB   (single-channel: %+.3f +/- %.3f)\n', ...
                out.deltaBoth, out.deltaBothErr, out.delta, getDeltaErr(out));
        if isfield(out, 'deltaCombined') && isfinite(out.deltaCombined) && strcmp(out.verdict, 'agree')
            fprintf('  the two combined:       delta = %+.3f +/- %.3f dB   (for information)\n', ...
                    out.deltaCombined, out.deltaCombinedErr);
        end
    else
        fprintf('  no zero crossing: %s\n', f.why);
    end
    switch out.verdict
        case 'agree'
            fprintf(['  AGREE: the single-channel balance holds with both channels playing. ' ...
                     '4.3 may write it.\n']);
        case 'disagree'
            fprintf(['  DISAGREE by %.1f sigma: the balance is not the same with both channels ' ...
                     'on. Re-run the phase\n  scan at delta = %+.3f (stage 1 skipped), then this ' ...
                     'again; level4.m prints the calls.\n'], ...
                    abs(out.dStar) / max(hypot(out.dStarErr, getDeltaErr(out)), eps), out.deltaBoth);
            if out.extrapolated
                fprintf(['  The crossing is OUTSIDE the offsets played (%s dB), so %+.3f is an ' ...
                         'extrapolation -- in\n  compression a straight line falls short. The ' ...
                         'next 4.2b, centred on it, measures it properly.\n'], ...
                        mat2str(out.fit.offsets), out.deltaBoth);
            end
        otherwise
            fprintf(['  INCONCLUSIVE: %s\n  More averages, or a stronger drive (the slope ' ...
                     'doubles once Omega >> the 2.16 MHz hyperfine).\n'], ...
                    ternary(isempty(f.why), sprintf('d* error %.2f dB > maxErrDb %.2f', ...
                    out.dStarErr, out.opts.maxErrDb), f.why));
    end
    if isfield(out, 'quadAngle') && out.quadAngle > 1
        fprintf(['  NOTE 3pi/2: a channel skew biases this the OPPOSITE way to pi/2. Compare the ' ...
                 'two:\n  qc = mwQuadratureCombine(q, q3).\n']);
    else
        fprintf(['  NOTE a timing skew between the channels'' pulses reads as an imbalance here ' ...
                 '(~0.08 dB\n  per ns per MHz of drive). A 3pi/2 run (quadAngle 1.5) and ' ...
                 'mwQuadratureCombine tell the two apart.\n']);
    end
    fprintf('==========================================================\n');
end

function drawIt(out)
    f = out.fit;
    figure('Name', 'calibMwQuadrature', 'Color', 'w');
    errorbar(f.offsets, f.Qnorm, f.QnormErr, 'o', 'LineWidth', 1); hold on;
    if isfield(f, 'line')
        xx = linspace(min(f.offsets) - 0.3, max(f.offsets) + 0.3, 50);
        plot(xx, f.line(xx) / f.contrast, '-', 'LineWidth', 1.5);
    end
    plot(xlim, [0 0], 'k:');
    if isfinite(out.dStar)
        plot(out.dStar * [1 1], ylim, '--');
    end
    grid on;
    xlabel(sprintf('channel 2 offset from %+.3f dB', out.delta));
    ylabel('Q / contrast');
    ang = 'pi/2';
    if isfield(out, 'quadAngle') && out.quadAngle > 1; ang = '3pi/2'; end
    title(sprintf('%s: delta = %+.3f +/- %.3f dB  (%s)', ang, out.deltaBoth, out.deltaBothErr, ...
                  out.verdict));
end

function m = estimateMinutes(b, opts, nCfg)
    % Sequence time per pass ~ 2 halves x (laser init + pulse gate + overhead), plus one
    % re-arm (~0.6 s) per pg.run.
    tHalf = b.laserInitializationDuration + max(b.tau) + 1;            % us
    perPass = nCfg + opts.walkWarmupTriggers;
    passesPerRun = max(1, floor(opts.maxDigitizerSegments / (2 * perPass)));
    nRuns = ceil(opts.repeats / passesPerRun) * opts.averages;
    m = (2 * tHalf * 1e-6 * perPass * opts.repeats * opts.averages + 0.6 * nRuns) / 60;
end

function e = getDeltaErr(out)
    e = 0;
    if isfield(out, 'deltaErr') && isfinite(out.deltaErr); e = out.deltaErr; end
end

function x = nanToZero(x)
    if isempty(x) || ~isfinite(x); x = 0; end
end

function v = pickVal(given, J, field, dflt, idx)
    v = given;
    if ~isempty(v); return; end
    v = dflt;
    if isfield(J, field)
        x = J.(field);
        if numel(x) >= idx; v = x(idx); end
    end
end

function printIndented(txt)
    if isempty(txt); return; end
    L = strsplit(txt, newline);
    for i = 1:numel(L)
        if ~isempty(strtrim(L{i})); fprintf('      | %s\n', L{i}); end
    end
end

function out = ternary(c, a, b)
    if c; out = a; else; out = b; end
end

function o = withDefaults(o, d, who)
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d); given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('%s: unknown option "%s". Valid: %s.', who, bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || (isempty(o.(known{i})) && ~isempty(d.(known{i})))
            o.(known{i}) = d.(known{i});
        end
    end
end
