function out = ramseyGateTest(opts)
%RAMSEYGATETEST Is the fixed MW gate what costs the Ramsey walk its contrast? ~3.5 min.
%
%   out = ramseyGateTest
%   out = ramseyGateTest(struct('repeats', 3000))     % if the fringe comes out marginal
%   out = ramseyGateTest(struct('selfTest', true))    % 1 s, no instruments
%
% ===================================================================
% WHERE THIS SITS
% ===================================================================
%
% ramseyWalkDiag, at 300 passes and 34 taus:
%
%       run           contrast     mean sterr
%       classic       0.001107     0.0002163
%       walk P=1      0.0004962    0.0002356     <-- 0.45 of classic
%       walk P=228    0.0003453    0.0002177     <-- 0.31
%
% THE ERROR BARS MATCH AND THE CONTRAST DOES NOT. That rules out the block-to-tau
% bookkeeping: averaging, smearing and mis-mapping all move the NOISE, because walkReduce
% derives the error bar from the scatter across passes. Noise right and signal halved means
% different PULSES, not different labels.
%
% walkAudit then removed most of the suspects. Rabi walk came back at 98% of classic
% contrast with a gate-tail tilt of +0.0003 +/- 0.0002 -- flat -- on a 1.0 us excess gate,
% and ESR walk at 99%. So on Rabi, at a tail of the same order as Ramsey's:
%
%   * the gate left open over an idling AWG does NOT measurably leak
%   * the sequencer DOES advance one block per trigger
%   * the digitizer segments ARE aligned with the blocks
%
% Which is why this test is still worth running rather than assumed: Ramsey's tail is 2.5x
% longer than the one Rabi cleared, and -- the part Rabi cannot speak to -- during Ramsey's
% tail the state sits in a SUPERPOSITION rather than driven or along z. A tail that leaves a
% z-aligned state untouched still dephases a superposition, and that is a contrast-only loss
% with unchanged shot noise, which is exactly the observation.
%
% ===================================================================
% THE TEST
% ===================================================================
%
% Vary the excess gate and NOTHING else. Four runs on the SAME short tau list:
%
%   A  classic                          reference; the gate follows tau, as always
%   B  walk, gate = tauShort            little excess
%   C  walk, gate forced to tauLong     by appending ONE long tau that is then dropped
%   D  classic again                    the drift control
%
% The gate comes from tau(END) -- AWGExpRamsey builds
% S.addEvent(obj.tau(end) + 2*obj.halfPiTime, MW, 'tau') -- so appending one long point
% lengthens the gate while every other point stays exactly where it was. B and C then differ
% in the excess gate and in literally nothing else, which is what makes this conclusive
% rather than suggestive.
%
%   B keeps its contrast, C loses it   -> THE GATE, conclusively
%   both B and C lose it               -> the gate is innocent, and the fault is something
%                                         Ramsey has that Rabi does not. The four detection
%                                         windows are then the place to look.
%   both fine                          -> the loss needs the WIDE tau list; raise tauShort
%                                         and find where it starts
%
% ===================================================================
% WHAT IT BORROWED FROM walkAudit, AND WHY
% ===================================================================
%
% Every one of these was learned the hard way in the last two runs on the setup:
%
%   THE DRIFT CONTROL. Two runs minutes apart cannot separate "walk is worse" from "the
%   second one is worse". A and D bracket B and C in time, so their difference is a measured
%   drift that any gate effect has to beat.
%
%   THE ERROR-BAR FACTOR. Two identical runs should differ by 0.674 sigma per point -- the
%   median of a folded standard normal. ESR came back at 12.3, i.e. sterr understated by
%   18x. Nothing quoted in sigma means anything until that factor is known, so it is
%   measured here from A against D and printed.
%
%   THE FRINGE GATE. On the first walkAudit run the arithmetic returned a confident
%   "-0.789 +/- 0.099, 8.0 sigma" on a trace with no resonance in it at all. Ramsey is far
%   more exposed to that than ESR was: its real contrast is about 0.1%, a twentieth of an
%   ESR dip. If A has no visible fringe then B and C would be a ratio between two baselines,
%   and this refuses to print one.
%
%   CONTRAST AS A FITTED SCALE, not a MAD statistic. walkAuditDecompose projects
%   (walk - classic) onto (classic - mean classic) and returns the scale WITH an error bar,
%   Birge-scaled so an understated sterr cannot manufacture significance. The earlier
%   moving-median MAD metric gave a bare number, which is how 0.45 and 0.31 came to be
%   compared as though the difference between them meant something.
%
% ===================================================================
% OPTIONS
% ===================================================================
%
%   tauShort (0.93)  us. The short list's span, and B's gate. Two full hyperfine beats
%                    (0.463 us each) so there is a real fringe to measure.
%   tauLong (2.5)    us. C's gate, set by the appended point. Match it to the tau span of
%                    the Ramsey runs you actually care about.
%   perBeat (15)     tau points per hyperfine beat.
%   detuneMHz (0)    on resonance, where the 2.16 MHz beat IS the carrier and the grid is
%                    coarsest. Detuning multiplies the point count for no extra T2*.
%   repeats (3000), averages (2)   averages must be >= 2 or sterr means nothing.
%                    3000 is not caution, it is the number this test NEEDS. The error on the
%                    contrast scale goes as 1/sqrt(repeats * nTau), and at 1000 repeats with
%                    26 points it comes out near 0.14 -- so "B kept its contrast" would mean
%                    "B is somewhere between 58% and 142%", which cannot distinguish the
%                    0.45 that started this investigation from 1.00. selfTest caught that
%                    before any beam time was spent on it. The audit prints the sensitivity
%                    it actually achieved, and there is a power check after run A that stops
%                    the test rather than let it produce an unusable answer.
%   passesPerRun ([]) [] lets AWGExpRamsey choose from its digitizer ceiling
%   initUs (16)      laserInitializationDuration. 81% of a Ramsey shot is laser
%                    initialisation and level 2 put the SNR-per-unit-time optimum at 15-20;
%                    the json has 29.4.
%   amplitude (-16), nChannels (2), jsonLocation, plotIt (true)
%   synthSeed (1)    which synthetic noise realisation selfTest uses. Vary it to check the
%                    verdict is not riding on one lucky draw.
%   selfTest (false) drive the whole pipeline on synthetic traces, no instruments. The
%                    long-gate run has 55% of its contrast removed on purpose, so a passing
%                    selfTest MUST print "THE GATE IS THE CAUSE".
%
% Writes nothing, saves nothing, leaves AWGExpRamsey on the classic path.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('jsonLocation', '', 'detuneMHz', 0, ...
        'amplitude', -16, 'nChannels', 2, 'initUs', 16, 'tauShort', 0.93, ...
        'tauLong', 2.5, 'perBeat', 15, 'repeats', 3000, 'averages', 2, ...
        'passesPerRun', [], 'plotIt', true, 'selfTest', false, ...
        'synthSeed', 1), 'ramseyGateTest');

    if isempty(opts.jsonLocation) && ~opts.selfTest
        opts.jsonLocation = setup7Paths().json;
    end
    if opts.averages < 2
        error(['ramseyGateTest: averages must be at least 2. With one average sterr is the ' ...
               'within-average scatter with nothing to compare it against, the drift control ' ...
               'cannot measure the error-bar factor, and every sigma printed would be a guess.']);
    end
    if opts.tauLong <= opts.tauShort
        % Name the fix rather than just the rule. Raising tauShort without raising tauLong is
        % the obvious thing to try and it lands here, so say what to type.
        error(['ramseyGateTest: tauLong (%g) must exceed tauShort (%g), or B and C have the ' ...
               'same gate and the test compares nothing.\n' ...
               'Try:  ramseyGateTest(struct(''tauShort'', %g, ''tauLong'', %g))'], ...
              opts.tauLong, opts.tauShort, opts.tauShort, opts.tauShort + 1.5);
    end

    % The grid has to resolve the FASTEST hyperfine component, |detuning| + 2.16 MHz, not
    % the detuning. On resonance that is 2.16 MHz, a 0.463 us beat.
    beat = 1 / (abs(opts.detuneMHz) + 2.16);
    step = beat / opts.perBeat;
    tauS = 0.005 : step : opts.tauShort;
    if numel(tauS) < 12
        error(['ramseyGateTest: tauShort = %g gives only %d points. The decomposition needs ' ...
               'enough to separate a contrast scale from a baseline -- raise tauShort to at ' ...
               'least %.2f us (two beats) or perBeat above %d.'], ...
              opts.tauShort, numel(tauS), 2*beat, opts.perBeat);
    end

    banner(tauS, beat, opts);

    % ---- the four runs. expGrab after each, BEFORE the next constructor -----
    % AWGExpRamsey is a singleton keyed on class name: the second constructor DELETES the
    % first and the old variable becomes a dead handle. There is no way to hold two.
    A = runOne(tauS, [],           opts, false, 'A  classic');

    % ---- CAN THIS RUN ANSWER THE QUESTION? ------------------------------
    % Checked after A and before spending three more runs, because the answer is knowable
    % from A alone: the error on the contrast scale is set by A's noise against A's fringe,
    % and neither B nor C will change it much. An underpowered version of this test does not
    % fail -- it returns "both kept their contrast" with error bars wide enough to hide the
    % effect, which reads as good news and is worse than no run at all.
    pw = powerCheck(A, opts);
    reportPower(pw, opts);
    if ~pw.usable
        out = struct('opts', opts, 'tau', tauS, 'A', A, 'ok', false, 'power', pw, ...
                     'when', nowText());
        fprintf(['\nSTOPPED after run A rather than spend three more runs on a test that\n' ...
                 'cannot distinguish 45%% from 100%%. Re-run with repeats = %d.\n\n'], pw.needRepeats);
        return
    end

    B = runOne(tauS, [],           opts, true,  'B  walk, short gate');
    C = runOne(tauS, opts.tauLong, opts, true,  'C  walk, gate forced long');
    D = runOne(tauS, [],           opts, false, 'D  classic again (drift control)');

    out = struct('opts', opts, 'tau', tauS, 'A', A, 'B', B, 'C', C, 'D', D, ...
                 'ok', false, 'power', pw, 'when', nowText());

    n = numel(tauS);
    runs = {A, B, C, D};
    for i = 1:numel(runs)
        if numel(runs{i}.value) ~= n
            fprintf(['\n*** the %s run returned %d points, expected %d. Look further up for ' ...
                     'the underlying error -- Experiment.run catches per-average failures and ' ...
                     'only warns. Nothing below would mean anything, so stopping.\n'], ...
                    runs{i}.tag, numel(runs{i}.value), n);
            verdictGate(out);
            return
        end
    end

    % ---- the classic reference, straddling B and C in time ------------------
    out.ref   = combineTwo(A, D);
    out.drift = driftBetween(A, D);
    reportDrift(out.drift);

    % ---- is there a fringe at all? -----------------------------------------
    out.feature = fringeStrength(tauS, out.ref.value, out.ref.sterr, A, D);
    reportFringe(out.feature);

    % ---- B and C against it -------------------------------------------------
    % A - D is pure noise, so it measures what sterr is supposed to say. The reference is
    % their mean, whose noise variance is a quarter of var(A - D).
    dAD = A.value(:) - D.value(:);
    out.refNoiseVar = var(dAD) / 4;
    out.runNoiseVar = var(dAD) / 2;
    fprintf(['\n  MEASURED NOISE (from A - D, no reference to sterr): per-run rms %.6f,\n' ...
             '  reference rms %.6f. Used for the dilution correction below.\n'], ...
            sqrt(out.runNoiseVar), sqrt(out.refNoiseVar));

    out.B = addCompare(B, out.ref, out.refNoiseVar);
    out.C = addCompare(C, out.ref, out.refNoiseVar);
    reportCompare(out.B, out.drift, out.feature);
    reportCompare(out.C, out.drift, out.feature);

    out.ok = true;
    if opts.plotIt
        try
            drawGateTest(out);
        catch err
            fprintf('\n(plot failed: %s -- the numbers above are unaffected)\n', err.message);
        end
    end
    verdictGate(out);
end

% ======================================================================
% RUNNING
% ======================================================================

function g = runOne(tauList, padTo, opts, walk, tag)
% padTo: [] for none, or a tau APPENDED purely to lengthen the MW gate. AWGExpRamsey builds
% the gate from tau(END), so one long point at the end forces a long gate while every other
% point keeps its own tau. Its data is dropped afterwards, so B and C are compared on
% identical taus and the gate is the only difference between them.
    fprintf('\n--- %s ---\n', tag);
    tau = tauList;
    if ~isempty(padTo); tau = [tauList, padTo]; end

    if opts.selfTest
        % Noise scaled from the real measurement: ramseyWalkDiag saw sterr 2.163e-4 at 300
        % repeats, and it goes as 1/sqrt(repeats * averages).
        noise = 2.163e-4 * sqrt(300 / max(opts.repeats * opts.averages / 2, 1));
        % Seed from the run's LETTER, so A, B, C and D are four independent realisations.
        % Deriving it from (walk, padded) gave A and D the same seed -- identical traces, a
        % drift of exactly zero, and a printed "effect is 2.2e15x the drift".
        g = synthRun(tauList, padTo, walk, tag, 10*opts.synthSeed + double(tag(1)), noise);
        return
    end

    r = AWGExpRamsey;
    r.loadJsonParams('all', opts.jsonLocation);
    r.nChannels = opts.nChannels;                       % BEFORE amplitude, always
    r.amplitude = opts.amplitude + zeros(1, opts.nChannels);
    r.frequency = r.frequency + opts.detuneMHz;
    if ~isempty(opts.initUs); r.laserInitializationDuration = opts.initUs; end
    r.repeats          = opts.repeats;
    r.averages         = opts.averages;
    r.constantTime     = 1;      % or the laser duty cycle changes across the scan and
                                 % biases T2* short, differently on the two paths
    r.balancedSequence = 1;
    r.isTracking       = 0;
    r.shouldAutosave   = false;
    r.tau              = tau;    % tau has NO set method, so it cannot be changed on a
                                 % reused object -- see the note in expGrab
    r.useSequenceWalk  = walk;
    if walk
        r.passesPerRun       = opts.passesPerRun;
        r.walkWarmupTriggers = 1;
        r.resyncEveryPass    = true;
    end

    styleWas = '';
    try; styleWas = r.AWG.progressStyle; r.AWG.progressStyle = 'milestone'; catch; end
    % ID FORM. A bare tic sets one global stopwatch that anything inside r.run resets, and
    % then this toc measures from there instead of from here.
    tId = tic;
    try
        r.run;
    catch err
        fprintf('  *** the run threw: %s\n', err.message);
    end
    secs = toc(tId);
    try; r.AWG.progressStyle = styleWas; catch; end
    try; r.useSequenceWalk = false; catch; end

    g = expGrab(r, tag);
    g.walk = walk; g.seconds = secs; g.padded = ~isempty(padTo);
    g.gateUs = tau(end) + 2*halfPiOf(g);
    if isempty(g.x); g.x = tau; end

    % Drop the padding point so all four runs are compared on identical taus.
    if ~isempty(padTo)
        k = numel(tauList);
        g.x     = firstN(g.x, k);
        g.value = firstN(g.value, k);
        g.sterr = firstN(g.sterr, k);
    end

    fprintf('  %.0f s, gate %.3f us, span %.5f, median sterr %.6f\n', ...
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
% SYNTHETIC RUNS, FOR selfTest
% ======================================================================

function g = synthRun(tauList, padTo, walk, tag, seedBase, noise)
% Injects a KNOWN gate effect: the long-gate run loses 55% of its contrast, the short-gate
% one keeps all of it. A passing selfTest must therefore print "THE GATE IS THE CAUSE"; if it
% prints anything else, the verdict is not reading the decomposition and the real run would
% be six minutes wasted.
    n = numel(tauList);
    t = tauList;
    fringe = 5e-4 * exp(-t/1.8) .* cos(2*pi*2.16*t);
    longGate = ~isempty(padTo);
    % A DISTINCT seed per run. The first version derived the seed from (walk, padded) alone,
    % which gave runs A and D the same one -- identical traces, a drift of exactly zero, and
    % a reported "effect is 2.2e15x the drift". Two runs that cannot differ are not a drift
    % control.
    seed = seedBase;
    amp = 1; if longGate; amp = 0.45; end
    v = 1 + amp*fringe + noise*synthNoise(n, seed);
    gate = 0.054 + max(t);
    if longGate; gate = 0.054 + padTo; end
    g = struct('tag', tag, 'x', t, 'value', v, 'sterr', noise*ones(1, n), ...
               'value2', [], 'sterr2', [], 'signal', [], 'rawSterr', [], ...
               'walk', walk, 'seconds', 30, 'padded', longGate, 'gateUs', gate, ...
               'settings', struct('halfPiTime', 0.027), 'name', 'synthetic');
    fprintf('  selfTest, gate %.3f us, injected contrast %.0f%%\n', gate, 100*amp);
end

function r = synthNoise(n, seed)
% Written-out LCG plus Box-Muller. rng() does not exist in Octave, where this file is
% parse-checked, and randn('seed',k) is deprecated in MATLAB -- so neither is used and
% selfTest is bit-reproducible on any release.
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
% STATISTICAL POWER
% ======================================================================

function pw = powerCheck(A, opts)
%POWERCHECK Will the contrast scale come back with a small enough error to be useful?
%
% The scale is a weighted projection of (walk - classic) onto (classic - mean classic), so
% its error is
%
%       scaleErr  ~  sigma_diff / ( rms(fringe) * sqrt(N) )
%
% with sigma_diff = sqrt(2) * the per-point error (walk and classic each contribute) and
% rms(fringe) the size of the thing being scaled. Both come straight out of run A.
%
% The bar: to call "B kept its contrast" the 3-sigma band has to be narrower than the
% effect being looked for, and the effect is 0.55. A band of +/- 0.18 does that comfortably,
% so the target is scaleErr <= 0.06. Anything worse than 0.11 cannot separate 0.45 from 1.00
% at 3 sigma even in the B-versus-C comparison, and is not worth the beam time.
    pw = struct('scaleErr', NaN, 'fringeRms', NaN, 'sigmaDiff', NaN, 'usable', false, ...
                'good', false, 'needRepeats', opts.repeats, 'n', numel(A.value));
    if isempty(A.value) || numel(A.value) < 8; return; end
    r = baselineResid(A.x(:), A.value(:));
    pw.fringeRms = std(r);
    e = A.sterr(isfinite(A.sterr) & A.sterr > 0);
    if isempty(e) || pw.fringeRms <= 0; return; end
    pw.sigmaDiff = sqrt(2) * median(e);
    pw.scaleErr = pw.sigmaDiff / (pw.fringeRms * sqrt(pw.n));

    pw.good   = pw.scaleErr <= 0.06;
    pw.usable = pw.scaleErr <= 0.11;
    % scaleErr scales as 1/sqrt(repeats), so the repeats needed for the 0.06 target is a
    % square. Rounded up to something round, because 4137 is false precision.
    need = opts.repeats * (pw.scaleErr / 0.06)^2;
    pw.needRepeats = max(opts.repeats, 500*ceil(need/500));
end

function reportPower(pw, opts)
    if ~isfinite(pw.scaleErr)
        fprintf('\n  POWER CHECK: could not be made (run A returned nothing usable).\n');
        return
    end
    fprintf('\n  POWER CHECK from run A: fringe rms %.6f against %.6f per-point error over\n', ...
            pw.fringeRms, pw.sigmaDiff/sqrt(2));
    fprintf('  %d points -> the contrast scale will come back to about +/- %.3f.\n', ...
            pw.n, pw.scaleErr);
    if pw.good
        fprintf(['  -> good. The 3-sigma band is %.0f%%, comfortably inside the 55%% effect this\n' ...
                 '     test is looking for.\n'], 100*3*pw.scaleErr);
    elseif pw.usable
        fprintf(['  -> MARGINAL. The 3-sigma band is %.0f%%, against a 55%% effect. B versus C\n' ...
                 '     will separate if the gate is guilty, but "B kept its contrast" will be a\n' ...
                 '     weak statement. Continuing; repeats = %d would make it clean.\n'], ...
                100*3*pw.scaleErr, pw.needRepeats);
    else
        fprintf(['  -> *** NOT ENOUGH. The 3-sigma band is %.0f%%, wider than the 55%% effect,\n' ...
                 '     so this configuration cannot tell 45%% from 100%% whatever it measures.\n' ...
                 '     It would return "both kept their contrast" with error bars big enough to\n' ...
                 '     hide the whole effect -- which reads as good news and is worse than not\n' ...
                 '     running at all.\n'], 100*3*pw.scaleErr);
        if pw.fringeRms < 2e-4
            fprintf(['     Note the fringe itself is small (rms %.6f). Before raising repeats,\n' ...
                     '     check halfPiTime in the json -- a wrong pi/2 costs contrast on every\n' ...
                     '     path at once and no number of repeats fixes it.\n'], pw.fringeRms);
        end
    end
end

% ======================================================================
% COMBINING, DRIFT, THE FRINGE GATE
% ======================================================================

function C = combineTwo(a, b)
% Inverse-variance mean of the two classic runs. The right reference precisely BECAUSE it
% straddles B and C in time, so it is contemporaneous with them to first order in the drift.
    ea = a.sterr(:)'; eb = b.sterr(:)';
    ea(~isfinite(ea) | ea <= 0) = Inf;
    eb(~isfinite(eb) | eb <= 0) = Inf;
    wa = 1./ea.^2; wb = 1./eb.^2; sw = wa + wb;
    C.value = (wa.*a.value(:)' + wb.*b.value(:)') ./ max(sw, eps);
    C.sterr = 1 ./ sqrt(max(sw, eps));
    C.x = a.x;
end

function d = driftBetween(a, b)
    d = struct('measured', true, 'scale', NaN, 'scaleErr', NaN, 'sigmaDiff', NaN, ...
               'errScale', NaN);
    % The basis here is run a alone, whose noise variance is var(a-b)/2.
    p = walkAuditDecompose(a.value, a.sterr, b.value, b.sterr, [], var(a.value(:)-b.value(:))/2);
    if p.ok; d.scale = p.scale; d.scaleErr = p.scaleErr; end
    sig = sqrt(a.sterr.^2 + b.sterr.^2); sig(sig <= 0) = eps;
    d.sigmaDiff = median(abs(a.value - b.value) ./ sig);
    d.errScale = d.sigmaDiff / 0.674;
end

function F = fringeStrength(t, y, e, a, b)
%FRINGESTRENGTH Is there a Ramsey fringe here, or only a baseline?
%
% Ramsey needs this more than ESR did. Its real contrast is about 0.1% -- a twentieth of an
% ESR dip -- so "the run produced numbers" and "the run produced a fringe" are much further
% apart, and the arithmetic downstream will happily report a contrast ratio between two
% things that are both baseline.
%
% Two criteria, because at 0.1% neither alone is enough:
%   sigma   peak departure from a quadratic baseline, in error bars. A quadratic absorbs
%           laser drift and leaves an oscillation nearly untouched.
%   repro   correlation of that residual between the two classic runs. A real fringe
%           reproduces; a wander does not. This is the criterion that carries the weight
%           here -- at this contrast one run cannot tell a fringe from drift, two can.
    F = struct('sigma', NaN, 'repro', NaN, 'resolved', false, 'ptp', NaN, 'medErr', NaN);
    t = t(:); y = y(:); e = e(:);
    if numel(t) < 8; return; end
    good = isfinite(e) & e > 0;
    med = median(e(good));
    if isempty(med) || ~isfinite(med) || med <= 0; return; end
    r = baselineResid(t, y);
    F.medErr = med;
    F.ptp    = max(r) - min(r);
    F.sigma  = max(abs(r)) / med;
    r1 = baselineResid(a.x(:), a.value(:));
    r2 = baselineResid(b.x(:), b.value(:));
    if numel(r1) == numel(r2) && numel(r1) > 4; F.repro = corrCol(r1, r2); end
    F.resolved = (F.sigma >= 4) && (~isfinite(F.repro) || F.repro >= 0.4);
end

function r = baselineResid(x, y)
    x = x(:); y = y(:);
    u = (x - mean(x)) / max(max(x) - min(x), eps);
    V = [ones(numel(u),1), u, u.^2];
    r = y - V*(V \ y);
end

function r = corrCol(a, b)
    a = a(:) - mean(a(:)); b = b(:) - mean(b(:));
    d = sqrt(sum(a.^2)*sum(b.^2));
    if d == 0; r = 0; else; r = sum(a.*b)/d; end
end

function g = addCompare(g, ref, noiseVar)
% noiseVar is the MEASURED noise variance of the reference, from the two classic runs, so the
% dilution correction does not depend on sterr -- which on this setup is not to be trusted.
    g.decomp  = walkAuditDecompose(ref.value, ref.sterr, g.value, g.sterr, [], noiseVar);
    g.metrics = walkCompareMetrics(ref.value, ref.sterr, g.value, g.sterr);
end

% ======================================================================
% REPORTING
% ======================================================================

function banner(tauS, beat, opts)
    fprintf('\n============================================================\n');
    fprintf('ramseyGateTest  --  is the fixed MW gate the cause?\n');
    fprintf('============================================================\n');
    fprintf('%d taus, %.4f us step, spanning %.3f us (%.1f hyperfine beats of %.3f us).\n', ...
            numel(tauS), median(diff(tauS)), max(tauS), max(tauS)/beat, beat);
    fprintf('%d repeats x %d averages.\n', opts.repeats, opts.averages);
    fprintf('\n  A  classic              the gate follows tau, as always\n');
    fprintf('  B  walk, short gate     gate ~ %.3f us   (little excess)\n', opts.tauShort + 0.054);
    fprintf('  C  walk, long gate      gate ~ %.3f us   (same taus as B)\n', opts.tauLong + 0.054);
    fprintf('  D  classic again        the drift control\n');
    fprintf(['\nB and C differ in the excess gate and in NOTHING else -- same taus, same\n' ...
             'repeats, same everything. That is what makes the answer conclusive.\n']);
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
        fprintf('  sterr is CONSERVATIVE by about %.1fx here. Safe, if generous.\n', ...
                1/max(d.errScale, eps));
    else
        fprintf('  sterr is honest to within %.1fx. Sigmas can be read at face value.\n', ...
                max(d.errScale, 1/max(d.errScale, eps)));
    end
end

function reportFringe(F)
    fprintf('\n  FRINGE CHECK: peak departure from a smooth baseline %.1f error bars', F.sigma);
    if isfinite(F.repro)
        fprintf(',\n  reproducing between the two classic runs at r = %+.2f.\n', F.repro);
    else
        fprintf('.\n');
    end
    if F.resolved
        fprintf('  -> there is a fringe, so a contrast ratio means something.\n');
    else
        fprintf(['  -> *** NO RESOLVED FRINGE in the classic reference. Ramsey contrast is\n' ...
                 '     about 0.1%%, so this is the ordinary failure of a short run rather than\n' ...
                 '     an exotic one. No contrast scale will be reported: with no fringe in the\n' ...
                 '     reference, B and C would be a ratio between two baselines.\n' ...
                 '     RAISE repeats (try 3000), AND check that halfPiTime in the json is still\n' ...
                 '     right -- a wrong pi/2 costs contrast on every path at once and would\n' ...
                 '     look exactly like this.\n']);
    end
end

function reportCompare(g, d, F)
    p = g.decomp; m = g.metrics;
    fprintf('\n  --- %s  (gate %.3f us) ---\n', g.tag, g.gateUs);
    fprintf('  correlation %+.4f against a noise-limited ceiling of %+.4f (ratio %.2f)\n', ...
            m.correlation, m.rMax, m.rRatio);
    fprintf('  best shift %+d\n', m.bestShift);
    if ~p.ok
        fprintf('  decomposition failed: %s\n', p.why);
        return
    end
    if ~F.resolved
        fprintf(['  CONTRAST SCALE withheld -- no fringe in the reference. It would have\n' ...
                 '  printed %+.3f.\n'], p.scale);
        return
    end
    fprintf('  CONTRAST SCALE %+.3f +/- %.3f   (%.1f sigma from zero)\n', ...
            p.scale, p.scaleErr, abs(p.scale)/max(p.scaleErr, eps));
    fprintf('     contrast is %.0f%% of classic''s.\n', 100*(1 + p.scale));
    fprintf('  reduced chi2 of the decomposition %.2f\n', p.chi2);
    if isfinite(p.attenuation)
        fprintf('  dilution: uncorrected %+.3f -> corrected %+.3f (attenuation %.2f, noise %s)\n', ...
                p.scaleRaw, p.scale, p.attenuation, p.noiseFrom);
        if ~isempty(p.why); fprintf('     %s\n', p.why); end
    end
    if isfinite(d.scale)
        fprintf('  effect is %.1fx the classic-to-classic drift.\n', ...
                abs(p.scale)/max(abs(d.scale), eps));
    end
end

% ======================================================================
% THE VERDICT
% ======================================================================

function verdictGate(out)
    fprintf('\n============================================================\n');
    fprintf('VERDICT\n');
    fprintf('============================================================\n');

    if ~out.ok
        if isfield(out, 'power') && isstruct(out.power) && ~out.power.usable
            fprintf(['\nNOT RUN -- the configuration could not have answered the question. See\n' ...
                     'the power check above: repeats = %d.\n\n'], out.power.needRepeats);
        else
            fprintf('\nThe runs did not all complete. Nothing about the gate follows.\n\n');
        end
        return
    end
    if ~out.feature.resolved
        fprintf(['\nINCONCLUSIVE -- the classic reference has no measurable fringe, so there is\n' ...
                 'no contrast for B and C to be a fraction of. See the fringe check above for\n' ...
                 'what to change. Note this is the correct outcome and not a failure of the\n' ...
                 'test: the arithmetic would happily have printed confident ratios between two\n' ...
                 'baselines, which is exactly what walkAudit did on its first run.\n\n']);
        return
    end

    pB = out.B.decomp; pC = out.C.decomp;
    % "Kept" and "lost" are decided against each coefficient's OWN error and a 10% floor,
    % not a fixed threshold. Birge scaling inside the decomposition means these errors
    % already absorb an understated sterr.
    keptB = abs(pB.scale) < max(3*pB.scaleErr, 0.10);
    keptC = abs(pC.scale) < max(3*pC.scaleErr, 0.10);
    sep = abs(pB.scale - pC.scale) / max(sqrt(pB.scaleErr^2 + pC.scaleErr^2), eps);

    fprintf('\n  B  short gate (%.2f us)   contrast %3.0f%% of classic   -> %s\n', ...
            out.B.gateUs, 100*(1 + pB.scale), eitherText(keptB, 'KEPT', 'LOST'));
    fprintf('  C  long gate  (%.2f us)   contrast %3.0f%% of classic   -> %s\n', ...
            out.C.gateUs, 100*(1 + pC.scale), eitherText(keptC, 'KEPT', 'LOST'));
    fprintf('  B and C differ by %.1f sigma.\n', sep);
    fprintf('  SENSITIVITY: a contrast change above %.0f%% would have shown at 3 sigma.\n', ...
            100*3*max(pB.scaleErr, pC.scaleErr));

    if keptB && ~keptC && sep > 3
        fprintf(['\nTHE GATE IS THE CAUSE, conclusively. B and C differ in NOTHING except how\n' ...
                 'long the microwave gate stays open past the end of the waveform: B keeps its\n' ...
                 'contrast (%.0f%%) and C loses it (%.0f%%).\n\n' ...
                 'And this is CONSISTENT with Rabi being clean, which is worth saying because it\n' ...
                 'looked like a contradiction: Rabi''s state is along z during its tail while\n' ...
                 'Ramsey''s is in a superposition. The tail dephases a superposition and leaves\n' ...
                 'a z-aligned state alone, so the same open gate costs Ramsey contrast and costs\n' ...
                 'Rabi nothing.\n\n' ...
                 'THE OPTIONS, best first:\n' ...
                 '  1. RIGHT-ALIGN the pulses in a fixed-length waveform, so the closing pi/2\n' ...
                 '     always ends at the gate close. Removes the excess gate entirely and keeps\n' ...
                 '     walk. Costs upload -- every waveform becomes max length.\n' ...
                 '  2. BAND the sweep: a few narrow tau ranges, each walked with its own gate.\n' ...
                 '     Small upload, one extra arm per band, no change to the waveform code.\n' ...
                 '  3. DROP walk for Ramsey. Classic with the on-resonance grid (167 -> 34\n' ...
                 '     points) and initialisation at 16 us is already about 5x faster than what\n' ...
                 '     you were running, with no assumptions at all.\n\n' ...
                 'Option 3 costs nothing and is available today. Option 2 is the cheapest way to\n' ...
                 'keep walk if the speed matters later.\n'], 100*(1+pB.scale), 100*(1+pC.scale));
    elseif ~keptB && ~keptC
        fprintf(['\nTHE GATE IS INNOCENT. Contrast is down to %.0f%% even with the gate barely\n' ...
                 'longer than the waveform, and B and C agree to %.1f sigma -- so the excess gate\n' ...
                 'is not what costs it.\n\n' ...
                 'Combined with walkAudit that is a tight corner. Rabi walk is clean at 98%%, so\n' ...
                 'the sequencer advances correctly, the segments line up with the blocks, and the\n' ...
                 'open gate does not leak. The fault has to be in something Ramsey has that Rabi\n' ...
                 'does not, and one structural difference is left:\n\n' ...
                 '  RAMSEY HAS FOUR DETECTION WINDOWS PER REPEAT. RABI HAS TWO.\n\n' ...
                 'signalParam is formed from rows 1-2 alone (publishRamsey), so a TWO-window\n' ...
                 'offset in the walk stream swaps the |1> and |0> projections. In Rabi a\n' ...
                 'two-window offset is a whole block and cannot happen, which is exactly why\n' ...
                 'Rabi''s clean result cannot rule this out. It is testable OFFLINE against a\n' ...
                 'synthetic four-window stream -- far cheaper than another night on the setup,\n' ...
                 'and walkReduce already has a test harness to extend.\n\n' ...
                 'Meanwhile: drop walk for Ramsey and take the 5x that classic already gives.\n'], ...
                100*(1+pB.scale), sep);
    elseif keptB && keptC
        fprintf(['\nBOTH WALKS ARE FINE HERE (%.0f%% and %.0f%% of classic), so the loss needs the\n' ...
                 'WIDE tau list to appear -- this test ran to %.2f us, ramseyWalkDiag ran to\n' ...
                 '2.5 us with 34 points.\n\n' ...
                 'That points away from the gate and toward something that scales with the tau\n' ...
                 'SPAN or the POINT COUNT: waveform length, sequence size, or blocks per pass.\n' ...
                 'Note C already had a %.2f us gate and was fine, which rules out the gate\n' ...
                 'length by itself -- so bisect on the number of taus, not on the gate. Re-run\n' ...
                 'with tauShort raised toward 2.5 and watch where it starts.\n'], ...
                100*(1+pB.scale), 100*(1+pC.scale), max(out.tau), out.C.gateUs);
    else
        fprintf(['\nBACKWARDS: the SHORT gate lost contrast (%.0f%%) and the long one kept it\n' ...
                 '(%.0f%%). No gate mechanism does that, so suspect the run rather than the\n' ...
                 'physics -- check the drift line above, and whether B landed on a bad stretch.\n' ...
                 'Repeat it. If it reproduces, the ORDER of the runs is doing something and that\n' ...
                 'is worth knowing on its own.\n'], 100*(1+pB.scale), 100*(1+pC.scale));
    end
    fprintf('\n============================================================\n\n');
end

% ======================================================================
% THE PLOT
% ======================================================================

function drawGateTest(out)
% Top: the traces. Bottom: both differences against the classic reference, which is where
% the answer is legible -- a contrast loss is a difference shaped like the FRINGE, and on the
% overlay panel a 55% loss and a perfect match both look like "close enough" at 0.1%.
    COL = struct('classic', [0.11 0.35 0.60], 'B', [0.85 0.37 0.10], ...
                 'C', [0.48 0.24 0.58], 'second', [0.55 0.62 0.70], ...
                 'band', [0.86 0.88 0.91]);

    figure('Color', 'w', 'Name', 'ramseyGateTest', 'Position', [60 60 1220 780]);

    subplot(2, 1, 1); hold on; grid on; box off; axStyle();
    H = []; L = {};
    H(end+1) = errbar(out.D.x, out.D.value, out.D.sterr, COL.second, 1.0, 4);
    L{end+1} = 'D classic again';
    H(end+1) = errbar(out.ref.x, out.ref.value, out.ref.sterr, COL.classic, 2.0, 6);
    L{end+1} = 'A+D classic reference';
    H(end+1) = errbar(out.B.x, out.B.value, out.B.sterr, COL.B, 1.8, 6);
    L{end+1} = sprintf('B walk, gate %.2f \\mus', out.B.gateUs);
    H(end+1) = errbar(out.C.x, out.C.value, out.C.sterr, COL.C, 1.8, 6);
    L{end+1} = sprintf('C walk, gate %.2f \\mus', out.C.gateUs);
    xlabel('\tau (\mus)'); ylabel('ratio');
    headroom(0.28); legendFrom(H, L); xlim(padRange(out.ref.x));
    title('Ramsey: does the excess MW gate cost contrast?', 'FontWeight', 'bold');

    subplot(2, 1, 2); hold on; grid on; box off; axStyle();
    H = []; L = {};
    sg = out.B.decomp.sigma;
    if ~isempty(sg) && numel(sg) == numel(out.ref.x)
        H(end+1) = fill([out.ref.x, fliplr(out.ref.x)], [sg, -fliplr(sg)], COL.band, ...
                        'EdgeColor', 'none');
        L{end+1} = '\pm1\sigma';
    end
    plot(padRange(out.ref.x), [0 0], '-', 'Color', [0.65 0.65 0.65]);
    H(end+1) = errbar(out.B.x, out.B.value - out.ref.value, out.B.decomp.sigma, COL.B, 1.6, 6);
    L{end+1} = 'B - classic  (short gate)';
    H(end+1) = errbar(out.C.x, out.C.value - out.ref.value, out.C.decomp.sigma, COL.C, 1.6, 6);
    L{end+1} = 'C - classic  (long gate)';
    xlabel('\tau (\mus)'); ylabel('walk - classic');
    headroom(0.28); legendFrom(H, L); xlim(padRange(out.ref.x));
    title(gateTitle(out), 'FontWeight', 'normal', 'FontSize', 9);
end

function T = gateTitle(out)
    if ~out.feature.resolved
        T = {sprintf('NO RESOLVED FRINGE -- %.1f error bars, r = %+.2f between the two classics', ...
                     out.feature.sigma, out.feature.repro), ...
             'contrast scales withheld: no fringe for them to be a fraction of'};
        return
    end
    pB = out.B.decomp; pC = out.C.decomp;
    T = {sprintf(['B  gate %.2f \\mus: %.0f%% of classic (%+.3f \\pm %.3f)      ' ...
                  'C  gate %.2f \\mus: %.0f%% of classic (%+.3f \\pm %.3f)'], ...
                 out.B.gateUs, 100*(1+pB.scale), pB.scale, pB.scaleErr, ...
                 out.C.gateUs, 100*(1+pC.scale), pC.scale, pC.scaleErr), ...
         sprintf(['classic-to-classic drift %+.3f   |   %s   |   B vs C differ by %.1f\\sigma'], ...
                 out.drift.scale, errBarPhrase(out.drift.errScale), ...
                 abs(pB.scale - pC.scale)/max(sqrt(pB.scaleErr^2 + pC.scaleErr^2), eps))};
end

function s = errBarPhrase(k)
% "error bars understated 1x" is not a statement, it is a rounding artefact. Say which of
% the three cases it is instead.
    if ~isfinite(k)
        s = 'error bars not checked';
    elseif k > 2
        s = sprintf('error bars understated %.0f\\times', k);
    elseif k < 0.5
        s = sprintf('error bars conservative %.1f\\times', 1/max(k, eps));
    else
        s = 'error bars honest';
    end
end

function h = errbar(x, y, e, col, lw, ms)
    if isempty(e) || numel(e) ~= numel(y); e = zeros(size(y)); end
    try
        h = errorbar(x, y, e, 'o-', 'Color', col, 'LineWidth', lw, 'MarkerSize', ms, ...
                     'MarkerFaceColor', col, 'MarkerEdgeColor', 'w');
    catch
        % Octave's errorbar takes no property/value pairs after the format string.
        h = errorbar(x, y, e, 'o-');
        try
            set(h, 'Color', col, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', col);
        catch
        end
    end
    try; set(h, 'CapSize', 0); catch; end
end

function axStyle()
    try
        set(gca, 'GridAlpha', 0.12, 'Layer', 'top', 'TickDir', 'out', 'FontSize', 9);
    catch
        set(gca, 'FontSize', 9);           % GridAlpha is R2014b and later
    end
end

function legendFrom(H, L)
% Fixed location with headroom reserved for it, not 'best': 'best' is only as good as the
% renderer's collision search and it put the legend through the data in the first version of
% walkAudit's figure.
    try
        lg = legend(H, L, 'Location', 'northeast');
        try; set(lg, 'Box', 'off'); catch; end
    catch
    end
end

function headroom(frac)
    yl = ylim; span = max(yl(2) - yl(1), eps);
    ylim([yl(1), yl(2) + frac*span]);
end

function r = padRange(x)
    lo = min(x); hi = max(x); w = max(hi - lo, eps);
    r = [lo - 0.02*w, hi + 0.02*w];
end

% ======================================================================

function s = eitherText(tf, a, b)
    if tf; s = a; else; s = b; end
end

function s = nowText()
    try
        s = datestr(now, 'yyyy-mm-dd HH:MM:SS');   %#ok<TNOW1,DATST>
    catch
        s = '';
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
        if ~isfield(o, known{i}); o.(known{i}) = d.(known{i}); end
    end
end
