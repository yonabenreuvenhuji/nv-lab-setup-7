function out = echoGateTest(opts)
%ECHOGATETEST Does Echo's fixed MW gate cost contrast? Run this BEFORE trusting a walked Echo.
%
%   out = echoGateTest
%   out = echoGateTest(struct('repeats', 3000))
%   out = echoGateTest(struct('selfTest', true))     % 1 s, no instruments
%
% ===================================================================
% WHY THIS IS NOT OPTIONAL FOR ECHO, THE WAY IT NEARLY WAS FOR RAMSEY
% ===================================================================
%
% Walk mode cannot rewrite the microwave gate per point, so it stays at the LONGEST tau. The
% excess -- gate open, AWG idling at zero -- is what has to be harmless.
%
%       Ramsey, measured: 4.05 us of excess cost nothing. 106% +/- 5% of classic contrast.
%       Echo, a 1-50 us list at cycles = 1: 98 us of excess. TWENTY-FIVE TIMES more.
%       Echo, a 1-200 us list at cycles = 8: 3.2 ms. Eight hundred times more.
%
% Nothing measured on this setup speaks to that. Leakage through an open switch is a rate, so
% whatever it does over 4 us it does twenty-five times as much over 98, and the effect is a
% loss of contrast at SHORT tau -- which on an Echo is indistinguishable from a shorter T2.
% That is the failure to fear: not a crash, a plausible number.
%
% The good news is that Echo is otherwise the easiest of the four to walk. tau lives in the
% AWG sequence as the repeat count of a zero step, so the waveform set is fixed and tiny and
% walking costs no extra upload at all. The gate is the only thing standing in the way, and
% this measures it.
%
% ===================================================================
% THE TEST
% ===================================================================
%
% Vary the excess gate and NOTHING else. Four runs on the SAME tau list:
%
%   A  classic                    reference; the gate follows tau, as always
%   B  walk, gate = max(tau)      the excess you would actually run with
%   C  walk, gate forced LONGER   by appending one long tau that is then dropped
%   D  classic again              the drift control
%
% AWGExpEcho builds the gate from max(obj.tau), so appending one long point lengthens the gate
% while every other point stays exactly where it was. B and C then differ in the excess gate
% and in nothing else -- which is what makes the answer conclusive rather than suggestive.
%
%   B keeps its contrast, C loses it   -> the gate leaks, and the leak scales with the excess.
%                                         Band the tau list: walk a few narrow ranges, each
%                                         with its own gate. Or drop walk for Echo.
%   both lose it                       -> the gate leaks even at B's excess. Do not walk Echo
%                                         on this tau range at all.
%   both fine                          -> the gate is harmless out to C's excess. Walk freely
%                                         up to that much, and re-test if you go further.
%
% Everything else -- the drift control, the error-bar factor, the feature gate, the power
% check -- is carried over from ramseyGateTest for the reasons documented there. The one
% Echo-specific difference is what counts as a feature: an Echo has a DECAY, not a fringe, so
% the check asks whether the trace actually falls rather than whether it oscillates.
%
% OPTIONS
%   tauMax (30)       us. The top of the list, and B's gate.
%   tauLongMax (120)  us. C's gate, set by the appended point. Push it to the longest tau you
%                     would realistically walk.
%   nTau (16), cycles (1), pulsePhase ('Y')
%   repeats (3000), averages (2)   averages must be >= 2 or sterr means nothing
%   passesPerRun ([]), maxRunSeconds (30)
%   initUs (16), amplitude (-16), nChannels (2), jsonLocation
%   plotIt (true), selfTest (false), synthSeed (1)
%
% Writes nothing, saves nothing, leaves AWGExpEcho on the classic path.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('jsonLocation', '', 'amplitude', -16, 'nChannels', 2, ...
        'initUs', 16, 'tauMax', 30, 'tauLongMax', 120, 'nTau', 16, 'cycles', 1, ...
        'pulsePhase', 'Y', 'repeats', 3000, 'averages', 2, 'passesPerRun', [], ...
        'maxRunSeconds', 30, 'plotIt', true, 'selfTest', false, 'synthSeed', 1), ...
        'echoGateTest');

    if isempty(opts.jsonLocation) && ~opts.selfTest
        opts.jsonLocation = setup7Paths().json;
    end
    if opts.averages < 2
        error(['echoGateTest: averages must be at least 2. At one average sterr is the ' ...
               'within-average scatter with nothing to compare it against, the drift control ' ...
               'cannot measure the error-bar factor, and every sigma printed is a guess.']);
    end
    if opts.tauLongMax <= opts.tauMax
        error(['echoGateTest: tauLongMax (%g) must exceed tauMax (%g), or B and C have the ' ...
               'same gate and the test compares nothing.\nTry: echoGateTest(struct(' ...
               '''tauMax'', %g, ''tauLongMax'', %g))'], ...
              opts.tauLongMax, opts.tauMax, opts.tauMax, 4*opts.tauMax);
    end

    % Log spacing, because an Echo decay is what is being measured and the interesting part is
    % the knee. A linear grid spends most of its points on the flat tail.
    tau = unique(roundTo(logspace(log10(1), log10(opts.tauMax), opts.nTau), 3));
    if numel(tau) < 8
        error('echoGateTest: only %d distinct taus. Raise nTau or tauMax.', numel(tau));
    end

    banner(tau, opts);

    A = runOne(tau, [], opts, false, 'A  classic');
    pw = powerCheck(A, opts);
    reportPower(pw, opts);
    if ~pw.usable
        out = struct('opts', opts, 'tau', tau, 'A', A, 'ok', false, 'power', pw);
        fprintf(['\nSTOPPED after run A rather than spend three more runs on a test that\n' ...
                 'cannot tell a 50%% contrast loss from none. Re-run with repeats = %d.\n\n'], ...
                pw.needRepeats);
        return
    end

    B = runOne(tau, [],                opts, true,  'B  walk, gate = max(tau)');
    C = runOne(tau, opts.tauLongMax,   opts, true,  'C  walk, gate forced long');
    D = runOne(tau, [],                opts, false, 'D  classic again (drift control)');

    out = struct('opts', opts, 'tau', tau, 'A', A, 'B', B, 'C', C, 'D', D, ...
                 'ok', false, 'power', pw);

    n = numel(tau); runs = {A, B, C, D};
    for i = 1:numel(runs)
        if numel(runs{i}.value) ~= n
            fprintf(['\n*** the %s run returned %d points, expected %d. Look above for the ' ...
                     'underlying error -- Experiment.run catches per-average failures and only ' ...
                     'warns. Nothing below would mean anything, so stopping.\n'], ...
                    runs{i}.tag, numel(runs{i}.value), n);
            verdictGate(out); return
        end
    end

    out.ref   = combineTwo(A, D);
    out.drift = driftBetween(A, D);
    reportDrift(out.drift);

    dAD = A.value(:) - D.value(:);
    out.refNoiseVar = var(dAD)/4;
    fprintf(['\n  MEASURED NOISE (from A - D, no reference to sterr): per-run rms %.6f,\n' ...
             '  reference rms %.6f. Used for the dilution correction.\n'], ...
            sqrt(var(dAD)/2), sqrt(out.refNoiseVar));

    out.feature = decayStrength(tau, out.ref.value, out.ref.sterr, A, D);
    reportFeature(out.feature);

    out.B = addCompare(B, out.ref, out.refNoiseVar);
    out.C = addCompare(C, out.ref, out.refNoiseVar);
    reportCompare(out.B, out.drift, out.feature);
    reportCompare(out.C, out.drift, out.feature);

    out.ok = true;
    if opts.plotIt
        try; drawGate(out); catch err
            fprintf('\n(plot failed: %s -- the numbers above are unaffected)\n', err.message);
        end
    end
    verdictGate(out);
end

% ======================================================================
% RUNNING
% ======================================================================

function g = runOne(tauList, padTo, opts, walk, tag)
% padTo: [] for none, or a tau APPENDED purely to lengthen the gate. AWGExpEcho builds the
% gate from max(obj.tau), so one long point forces a long gate while every other point keeps
% its own tau. Its data is dropped afterwards, so B and C are compared on identical taus and
% the gate is the only difference between them.
    fprintf('\n--- %s ---\n', tag);
    tau = tauList;
    if ~isempty(padTo); tau = [tauList, padTo]; end

    if opts.selfTest
        g = synthRun(tauList, padTo, walk, tag, 10*opts.synthSeed + double(tag(1)), opts);
        return
    end

    e = AWGExpEcho;
    e.loadJsonParams('all', opts.jsonLocation);
    e.nChannels = opts.nChannels;                   % BEFORE amplitude, always
    e.amplitude = opts.amplitude + zeros(1, opts.nChannels);
    e.timeDelay = 0;                                % channelPhase is the phase knob
    if ~isempty(opts.initUs); e.laserInitializationDuration = opts.initUs; end
    e.repeats          = opts.repeats;
    e.averages         = opts.averages;
    e.cycles           = opts.cycles;
    e.pulsePhase       = opts.pulsePhase;
    e.constantTime     = 1;
    e.balancedSequence = 1;
    e.isTracking       = 0;
    e.shouldAutosave   = false;
    e.tau              = tau;
    e.useSequenceWalk  = walk;
    if walk
        e.passesPerRun       = opts.passesPerRun;
        e.walkWarmupTriggers = 1;
        e.resyncEveryPass    = true;
        e.maxRunSeconds      = opts.maxRunSeconds;
    end

    styleWas = '';
    try; styleWas = e.AWG.progressStyle; e.AWG.progressStyle = 'milestone'; catch; end
    tId = tic;                     % ID FORM: a bare tic is reset by any tic inside e.run
    try
        e.run;
    catch err
        fprintf('  *** the run threw: %s\n', err.message);
    end
    secs = toc(tId);
    try; e.AWG.progressStyle = styleWas; catch; end
    try; e.useSequenceWalk = false; catch; end

    g = expGrab(e, tag);
    g.walk = walk; g.seconds = secs; g.padded = ~isempty(padTo);
    g.gateUs = 2*max(tau)*opts.cycles + 2*halfPiOf(g);
    if isempty(g.x); g.x = 2*tauList*opts.cycles; end

    if ~isempty(padTo)
        k = numel(tauList);
        g.x = firstN(g.x, k); g.value = firstN(g.value, k); g.sterr = firstN(g.sterr, k);
    end
    fprintf('  %.0f s, gate %.2f us, span %.5f, median sterr %.6f\n', ...
            secs, g.gateUs, spanOf(g.value), median(g.sterr));
end

function h = halfPiOf(g)
    h = 0.027;
    try
        if isfield(g.settings, 'halfPiTime') && ~isempty(g.settings.halfPiTime)
            h = g.settings.halfPiTime;
        end
    catch
    end
end

function v = firstN(v, k)
    if numel(v) > k; v = v(1:k); end
end

function s = spanOf(v)
    if isempty(v); s = NaN; else; s = max(v) - min(v); end
end

% ======================================================================
% SYNTHETIC, FOR selfTest
% ======================================================================

function g = synthRun(tauList, padTo, walk, tag, seed, opts)
% Injects a KNOWN gate effect: the long-gate run loses 55% of its contrast, the short-gate one
% keeps all of it. A passing selfTest must print "THE GATE LEAKS". If it prints anything else,
% the verdict is not reading the decomposition and the real run would be four minutes wasted.
    n = numel(tauList);
    x = 2*tauList*opts.cycles;
    decay = 0.012 * exp(-(x/40).^1.5);          % a plausible Echo: 1.2% at t=0, T2 ~ 40 us
    noise = 2.163e-4 * sqrt(300 / max(opts.repeats*opts.averages/2, 1));
    longGate = ~isempty(padTo);
    amp = 1; if longGate; amp = 0.45; end
    v = 0.995 + amp*decay + noise*synthNoise(n, seed);
    gate = 2*max(tauList)*opts.cycles + 0.054;
    if longGate; gate = 2*padTo*opts.cycles + 0.054; end
    g = struct('tag', tag, 'x', x, 'value', v, 'sterr', noise*ones(1,n), ...
               'value2', [], 'sterr2', [], 'signal', [], 'rawSterr', [], ...
               'walk', walk, 'seconds', 40, 'padded', longGate, 'gateUs', gate, ...
               'settings', struct('halfPiTime', 0.027), 'name', 'synthetic');
    fprintf('  selfTest, gate %.2f us, injected contrast %.0f%%\n', gate, 100*amp);
end

function r = synthNoise(n, seed)
% Written-out LCG plus Box-Muller: rng() does not exist in Octave, where this is
% parse-checked, and randn('seed',k) is deprecated in MATLAB.
    m = 2^31 - 1;
    x = mod(1103515245*(seed + 12345) + 12345, m);
    u = zeros(1, 2*n);
    for i = 1:2*n
        x = mod(1103515245*x + 12345, m);
        u(i) = (x + 0.5)/m;
    end
    r = sqrt(-2*log(u(1:n))) .* cos(2*pi*u(n+1:2*n));
end

% ======================================================================
% POWER, DRIFT, THE FEATURE GATE
% ======================================================================

function pw = powerCheck(A, opts)
% Will the contrast scale come back with a small enough error to be useful? Answered from run
% A alone, before three more runs are spent. An underpowered version of this test does not
% fail -- it returns "both kept their contrast" with error bars wide enough to hide the whole
% effect, which reads as good news and is worse than not running.
    pw = struct('scaleErr', NaN, 'sigRms', NaN, 'usable', false, 'good', false, ...
                'needRepeats', opts.repeats, 'n', numel(A.value));
    if isempty(A.value) || numel(A.value) < 6; return; end
    v = A.value(:); pw.sigRms = std(v);
    e = A.sterr(isfinite(A.sterr) & A.sterr > 0);
    if isempty(e) || pw.sigRms <= 0; return; end
    pw.scaleErr = sqrt(2)*median(e) / (pw.sigRms * sqrt(pw.n));
    pw.good   = pw.scaleErr <= 0.06;
    pw.usable = pw.scaleErr <= 0.11;
    pw.needRepeats = max(opts.repeats, 500*ceil(opts.repeats*(pw.scaleErr/0.06)^2/500));
end

function reportPower(pw, opts)
    if ~isfinite(pw.scaleErr)
        fprintf('\n  POWER CHECK: could not be made (run A returned nothing usable).\n'); return
    end
    fprintf('\n  POWER CHECK from run A: decay rms %.6f over %d points -> the contrast scale\n', ...
            pw.sigRms, pw.n);
    fprintf('  will come back to about +/- %.3f.\n', pw.scaleErr);
    if pw.good
        fprintf('  -> good. The 3-sigma band is %.0f%%, well inside a 50%% effect.\n', ...
                100*3*pw.scaleErr);
    elseif pw.usable
        fprintf(['  -> MARGINAL. The 3-sigma band is %.0f%%. B versus C will separate if the\n' ...
                 '     gate leaks, but "B kept its contrast" will be a weak statement.\n' ...
                 '     repeats = %d would make it clean.\n'], 100*3*pw.scaleErr, pw.needRepeats);
    else
        fprintf(['  -> *** NOT ENOUGH. The 3-sigma band is %.0f%%, wider than the effect, so\n' ...
                 '     this configuration cannot answer its own question whatever it measures.\n'], ...
                100*3*pw.scaleErr);
    end
end

function C = combineTwo(a, b)
    ea = a.sterr(:)'; eb = b.sterr(:)';
    ea(~isfinite(ea) | ea <= 0) = Inf; eb(~isfinite(eb) | eb <= 0) = Inf;
    wa = 1./ea.^2; wb = 1./eb.^2; sw = wa + wb;
    C.value = (wa.*a.value(:)' + wb.*b.value(:)') ./ max(sw, eps);
    C.sterr = 1 ./ sqrt(max(sw, eps));
    C.x = a.x;
end

function d = driftBetween(a, b)
    d = struct('measured', true, 'scale', NaN, 'scaleErr', NaN, 'sigmaDiff', NaN, 'errScale', NaN);
    p = walkAuditDecompose(a.value, a.sterr, b.value, b.sterr, [], ...
                           var(a.value(:) - b.value(:))/2);
    if p.ok; d.scale = p.scale; d.scaleErr = p.scaleErr; end
    sig = sqrt(a.sterr.^2 + b.sterr.^2); sig(sig <= 0) = eps;
    d.sigmaDiff = median(abs(a.value - b.value) ./ sig);
    d.errScale = d.sigmaDiff / 0.674;
end

function F = decayStrength(tau, y, e, a, b)
%DECAYSTRENGTH Is there an Echo DECAY here, or only a baseline?
%
% The Echo analogue of ramseyGateTest's fringe check, and the one place the two differ: an
% Echo does not oscillate, it falls. So the question is whether the trace departs from FLAT
% by enough to measure a contrast on -- a straight line, not a quadratic, because a quadratic
% would absorb most of a decay and leave nothing.
    F = struct('sigma', NaN, 'repro', NaN, 'resolved', false, 'medErr', NaN);
    y = y(:); e = e(:);
    if numel(y) < 6; return; end
    med = median(e(isfinite(e) & e > 0));
    if isempty(med) || ~isfinite(med) || med <= 0; return; end
    r = y - mean(y);
    F.medErr = med;
    F.sigma = (max(y) - min(y)) / med;
    r1 = a.value(:) - mean(a.value); r2 = b.value(:) - mean(b.value);
    if numel(r1) == numel(r2) && numel(r1) > 4; F.repro = corrCol(r1, r2); end
    F.resolved = (F.sigma >= 6) && (~isfinite(F.repro) || F.repro >= 0.5);
end

function r = corrCol(a, b)
    a = a(:) - mean(a(:)); b = b(:) - mean(b(:));
    d = sqrt(sum(a.^2)*sum(b.^2));
    if d == 0; r = 0; else; r = sum(a.*b)/d; end
end

function g = addCompare(g, ref, noiseVar)
    g.decomp  = walkAuditDecompose(ref.value, ref.sterr, g.value, g.sterr, [], noiseVar);
    g.metrics = walkCompareMetrics(ref.value, ref.sterr, g.value, g.sterr);
end

% ======================================================================
% REPORTING
% ======================================================================

function banner(tau, opts)
    gB = 2*max(tau)*opts.cycles + 0.054;
    gC = 2*opts.tauLongMax*opts.cycles + 0.054;
    fprintf('\n============================================================\n');
    fprintf('echoGateTest  --  does the fixed MW gate leak?\n');
    fprintf('============================================================\n');
    fprintf('%d taus, %.1f to %.1f us, log spaced. cycles = %d, phases %s.\n', ...
            numel(tau), min(tau), max(tau), opts.cycles, opts.pulsePhase);
    fprintf('%d repeats x %d averages.\n', opts.repeats, opts.averages);
    fprintf('\n  A  classic                  the gate follows tau, as always\n');
    fprintf('  B  walk, gate %7.1f us    excess up to %.1f us at the shortest tau\n', ...
            gB, gB - (2*min(tau)*opts.cycles + 0.054));
    fprintf('  C  walk, gate %7.1f us    same taus as B, %.1fx the excess\n', ...
            gC, (gC - 2*min(tau)*opts.cycles)/max(gB - 2*min(tau)*opts.cycles, eps));
    fprintf('  D  classic again            the drift control\n');
    fprintf(['\nFor scale: ramseyGateTest cleared 4.05 us of excess gate. B is at %.0f us and\n' ...
             'C at %.0f us, so this is measuring something genuinely untested.\n'], ...
            gB - (2*min(tau)*opts.cycles + 0.054), gC - (2*min(tau)*opts.cycles + 0.054));
end

function reportDrift(d)
    fprintf('\n  DRIFT CONTROL (A vs D): contrast scale %+.3f +/- %.3f, median %.2f sigma\n', ...
            d.scale, d.scaleErr, d.sigmaDiff);
    fprintf('  per point. Any gate effect has to beat this.\n');
    fprintf('  ERROR BARS: two identical runs should differ by 0.67 sigma per point, so\n');
    if d.errScale > 2
        fprintf(['  sterr here is UNDERSTATED BY %.0fx. The contrast scales below survive it --\n' ...
                 '  walkAuditDecompose Birge-scales its covariance -- but do not read any OTHER\n' ...
                 '  sigma in this output at face value.\n'], d.errScale);
    elseif d.errScale < 0.5
        fprintf('  sterr is CONSERVATIVE by about %.1fx. Safe, if generous.\n', 1/max(d.errScale, eps));
    else
        fprintf('  sterr is honest to within %.1fx.\n', max(d.errScale, 1/max(d.errScale, eps)));
    end
end

function reportFeature(F)
    fprintf('\n  DECAY CHECK: the trace spans %.1f error bars', F.sigma);
    if isfinite(F.repro)
        fprintf(',\n  reproducing between the two classic runs at r = %+.2f.\n', F.repro);
    else
        fprintf('.\n');
    end
    if F.resolved
        fprintf('  -> there is a decay, so a contrast ratio means something.\n');
    else
        fprintf(['  -> *** NO RESOLVED DECAY. B and C would be a ratio between two baselines.\n' ...
                 '     Raise repeats, shorten tauMax so the decay is inside the window, and\n' ...
                 '     check piTime -- a wrong pi costs Echo contrast on every path at once.\n']);
    end
end

function reportCompare(g, d, F)
    p = g.decomp; m = g.metrics;
    fprintf('\n  --- %s  (gate %.2f us) ---\n', g.tag, g.gateUs);
    fprintf('  correlation %+.4f against a ceiling of %+.4f (ratio %.2f), best shift %+d\n', ...
            m.correlation, m.rMax, m.rRatio, m.bestShift);
    if ~p.ok
        fprintf('  decomposition failed: %s\n', p.why); return
    end
    if ~F.resolved
        fprintf('  CONTRAST SCALE withheld -- no decay in the reference. It would have printed %+.3f.\n', ...
                p.scale);
        return
    end
    fprintf('  CONTRAST SCALE %+.3f +/- %.3f   (%.1f sigma from zero)\n', ...
            p.scale, p.scaleErr, abs(p.scale)/max(p.scaleErr, eps));
    fprintf('     contrast is %.0f%% of classic''s.\n', 100*(1 + p.scale));
    fprintf('  reduced chi2 %.2f', p.chi2);
    if isfinite(p.attenuation)
        fprintf(', dilution: uncorrected %+.3f -> corrected %+.3f (attenuation %.2f, noise %s)', ...
                p.scaleRaw, p.scale, p.attenuation, p.noiseFrom);
    end
    fprintf('\n');
    if isfinite(d.scale)
        fprintf('  effect is %.1fx the classic-to-classic drift.\n', ...
                abs(p.scale)/max(abs(d.scale), eps));
    end
end

% ======================================================================
% VERDICT
% ======================================================================

function verdictGate(out)
    fprintf('\n============================================================\n');
    fprintf('VERDICT\n');
    fprintf('============================================================\n');
    if ~out.ok
        fprintf('\nThe runs did not all complete. Nothing about the gate follows.\n\n'); return
    end
    if ~out.feature.resolved
        fprintf(['\nINCONCLUSIVE -- no measurable decay in the classic reference, so there is no\n' ...
                 'contrast for B and C to be a fraction of. See the decay check above.\n\n']);
        return
    end
    pB = out.B.decomp; pC = out.C.decomp;
    keptB = abs(pB.scale) < max(3*pB.scaleErr, 0.10);
    keptC = abs(pC.scale) < max(3*pC.scaleErr, 0.10);
    sep = abs(pB.scale - pC.scale) / max(sqrt(pB.scaleErr^2 + pC.scaleErr^2), eps);

    fprintf('\n  B  gate %7.2f us   contrast %3.0f%% of classic   -> %s\n', ...
            out.B.gateUs, 100*(1 + pB.scale), eitherText(keptB, 'KEPT', 'LOST'));
    fprintf('  C  gate %7.2f us   contrast %3.0f%% of classic   -> %s\n', ...
            out.C.gateUs, 100*(1 + pC.scale), eitherText(keptC, 'KEPT', 'LOST'));
    fprintf('  B and C differ by %.1f sigma.\n', sep);
    fprintf('  SENSITIVITY: a contrast change above %.0f%% would have shown at 3 sigma.\n', ...
            100*3*max(pB.scaleErr, pC.scaleErr));

    if keptB && keptC
        fprintf(['\nTHE GATE IS HARMLESS OUT TO %.0f us OF EXCESS. Both walks keep their contrast\n' ...
                 '(%.0f%% and %.0f%%) and they agree to %.1f sigma, so leakage through the open\n' ...
                 'switch is below what this run can see even at C''s gate.\n\n' ...
                 'SO WALK ECHO FREELY up to about that much excess. Echo is the cheapest of the\n' ...
                 'three to walk -- tau lives in the sequence repeat counts, so the waveform set\n' ...
                 'never grows -- and the only thing that was in the way is now measured.\n\n' ...
                 'Re-test if you go substantially longer: leakage is a rate, so the excess is\n' ...
                 'what sets the dose, and a 3 ms XY8 gate is a different regime again.\n'], ...
                out.C.gateUs, 100*(1+pB.scale), 100*(1+pC.scale), sep);
    elseif keptB && ~keptC && sep > 3
        fprintf(['\nTHE GATE LEAKS, AND THE LEAK SCALES WITH THE EXCESS. B and C differ in nothing\n' ...
                 'except how long the switch stays open past the end of the pulse train: B keeps\n' ...
                 'its contrast (%.0f%%) and C loses it (%.0f%%).\n\n' ...
                 'WHAT TO DO, best first:\n' ...
                 '  1. BAND THE TAU LIST. Walk two or three narrow ranges, each with its own\n' ...
                 '     gate, so no block sees much excess. Costs one extra arm per band and\n' ...
                 '     keeps most of the saving. B''s gate is known to be safe, so bands no\n' ...
                 '     wider than that are known to be safe too.\n' ...
                 '  2. Keep walk only for short-tau Echoes, where max(tau) is small enough that\n' ...
                 '     the excess stays inside B''s.\n' ...
                 '  3. Drop walk for Echo. Classic costs one arm per tau, which on a 16-point\n' ...
                 '     list is about 20 s per average -- annoying, not prohibitive.\n'], ...
                100*(1+pB.scale), 100*(1+pC.scale));
    elseif ~keptB && ~keptC
        fprintf(['\nTHE GATE LEAKS ALREADY AT B (%.0f%% of classic), which is the excess you would\n' ...
                 'actually run with. DO NOT WALK ECHO ON THIS TAU RANGE.\n\n' ...
                 'Band the list into ranges narrow enough that the excess is small -- the Ramsey\n' ...
                 'measurement says 4 us is certainly safe -- or use classic. And note this is a\n' ...
                 'useful physics result in its own right: it says the microwave switch is not\n' ...
                 'clean over tens of microseconds, which is worth knowing for any long sequence,\n' ...
                 'walked or not.\n'], 100*(1+pB.scale));
    else
        fprintf(['\nBACKWARDS: the SHORT gate lost contrast (%.0f%%) and the long one kept it\n' ...
                 '(%.0f%%). No gate mechanism does that, so suspect the run -- check the drift\n' ...
                 'line above, and whether B landed on a bad stretch. Repeat it.\n'], ...
                100*(1+pB.scale), 100*(1+pC.scale));
    end
    fprintf('\n============================================================\n\n');
end

% ======================================================================
% PLOT
% ======================================================================

function drawGate(out)
    S = fitStyle();
    fig = S.newFigure('echoGateTest', true);
    axM = S.axes(subplot(3, 1, 1:2));
    H = []; L = {};
    H(end+1) = S.data(axM, out.D.x, out.D.value, out.D.sterr, S.col.s3); L{end+1} = 'D classic again';
    H(end+1) = S.data(axM, out.ref.x, out.ref.value, out.ref.sterr, S.col.fit); L{end+1} = 'A+D classic ref';
    H(end+1) = S.data(axM, out.B.x, out.B.value, out.B.sterr, S.col.s1);
    L{end+1} = sprintf('B walk, gate %.0f \\mus', out.B.gateUs);
    H(end+1) = S.data(axM, out.C.x, out.C.value, out.C.sterr, S.col.s2);
    L{end+1} = sprintf('C walk, gate %.0f \\mus', out.C.gateUs);
    xlabel(axM, 'total evolution time  [\mus]', 'Color', S.col.ink2);
    ylabel(axM, S.yLabel, 'Color', S.col.ink2);
    S.tight(axM, out.ref.x); S.headroom(axM, 0.30); S.legend(axM, H, L);
    S.title(axM, gateTitle(out));

    axR = S.axes(subplot(3, 1, 3));
    S.zeroLine(axR, [min(out.ref.x) max(out.ref.x)]);
    S.resid(axR, out.B.x, out.B.value - out.ref.value, S.col.s1);
    S.resid(axR, out.C.x, out.C.value - out.ref.value, S.col.s2);
    xlabel(axR, 'total evolution time  [\mus]', 'Color', S.col.ink2);
    ylabel(axR, 'walk - classic', 'Color', S.col.ink2);
    S.tight(axR, out.ref.x);
end

function T = gateTitle(out)
    if ~out.feature.resolved
        T = {sprintf('NO RESOLVED DECAY -- %.1f error bars of span', out.feature.sigma), ...
             'contrast scales withheld'};
        return
    end
    pB = out.B.decomp; pC = out.C.decomp;
    T = {sprintf(['B  gate %.0f \\mus: %.0f%% of classic (%+.3f \\pm %.3f)      ' ...
                  'C  gate %.0f \\mus: %.0f%% of classic (%+.3f \\pm %.3f)'], ...
                 out.B.gateUs, 100*(1+pB.scale), pB.scale, pB.scaleErr, ...
                 out.C.gateUs, 100*(1+pC.scale), pC.scale, pC.scaleErr), ...
         sprintf('drift %+.3f   |   B vs C differ by %.1f\\sigma', out.drift.scale, ...
                 abs(pB.scale - pC.scale)/max(sqrt(pB.scaleErr^2 + pC.scaleErr^2), eps))};
end

% ======================================================================

function y = roundTo(x, n)
% round(x, n) is R2014b and later in MATLAB and absent in Octave, where this is
% parse-checked. Written out so the file behaves the same in both.
    f = 10^n;
    y = round(x .* f) ./ f;
end

function s = eitherText(tf, a, b)
    if tf; s = a; else; s = b; end
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
        if ~isfield(o, known{i}); o.(known{i}) = d.(known{i}); end
    end
end
