function out = ramseyStageTimer(opts)
%RAMSEYSTAGETIMER Measure the REAL per-tau cost of a Ramsey, by running three short ones.
%
%   out = ramseyStageTimer                     % ~60 s of real runs
%
% Unlike ramseySpeedProbe, this RUNS THE EXPERIMENT. It has to: the per-tau cost is not the
% AWG re-arm alone, and there is no way to measure the rest from outside.
%
% ===================================================================
% WHY THIS EXISTS -- MY PREVIOUS NUMBER WAS WRONG
% ===================================================================
%
% ramseySpeedProbe timed awg.Run and called it "the per-tau overhead". It is not. Every tau in
% AWGExpRamsey.perform pays for ALL of this:
%
%       pg.changeSequence('tau', ...)            reprogram the pulse generator
%       pg.changeSequence('lastDelay', ...)      again, when constantTime is on
%       assignSequence  x nChannels              switch the AWG to this tau's sequence
%       awg.Run                                  <-- the only part I measured
%       getRawData -> pg.run + digitizer read    arm the counter, run, read back
%
% Timing one of five and projecting from it is how I ended up telling you the walk trick was
% worth 1.3x. The 0.6 s per point I had been quoting for Rabi came from a real Rabi timing
% probe, which included the whole chain -- so I compared a partial measurement against a
% complete one. That is the whole error.
%
% AND IT IS EXACTLY WHY WALK WON FOR RABI. Walk does not just re-arm less often; it removes
% the per-tau overhead ENTIRELY -- one pg.changeSequence never happens, one assignSequence
% covers the whole sweep, and one pg.run of (nTau x passes) triggers replaces nTau separate
% arm-run-read cycles. If the true per-tau cost is around a second, walk is worth 5-10x on
% Ramsey and not 1.3x.
%
% ===================================================================
% HOW IT SEPARATES THE TERMS WITHOUT INSTRUMENTING perform
% ===================================================================
%
% perform is protected, so the loop cannot be timed from here. But the total obeys
%
%       T(nTau, repeats) = A(nTau) + nTau * averages * (B + repeats * C)
%
%   A   upload and prepare, which the measured upload rate already predicts
%   B   the per-tau FIXED cost -- everything in the list above
%   C   the per-shot cost, i.e. the sequence duration
%
% Three short runs pin B and C:
%
%   M1  nTau = n0, repeats = rLow      M2 - M1 isolates C:  n0*(rHigh-rLow)*C
%   M2  nTau = n0, repeats = rHigh     M3 - M1 isolates B:  dn*(B + rLow*C)
%   M3  nTau = n0 + dn, repeats = rLow
%
% Everything is a FRESH object, because tau has no setter and a reused one would not reload.
%
% OPTIONS
%   jsonLocation, nTau0 (8), dnTau (8), repeatsLow (500), repeatsHigh (4000)
%   tauMax (2.5)   the span, so waveform lengths are representative of a real run
%   initUs ([])    override laserInitializationDuration; [] = the json's
%   detuneMHz (0)
%   amplitude (-16)
%
% OUTPUT
%   out.perTauSeconds   B, the number that decides everything
%   out.perShotUs       C in microseconds, to compare against the sequence duration
%   out.measured        the three (nTau, repeats, seconds) triples
%   out.projection      full-run estimates, classic and walk, from the REAL B

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('jsonLocation', '', 'nTau0', 8, 'dnTau', 8, ...
        'repeatsLow', 500, 'repeatsHigh', 4000, 'tauMax', 2.5, 'initUs', [], ...
        'detuneMHz', 0, 'amplitude', -16, 'averages', 1, ...
        'armSeconds', 0.222), 'ramseyStageTimer');   % what ramseySpeedProbe measured for
                                                     % awg.Run alone, for the comparison below

    if isempty(opts.jsonLocation); opts.jsonLocation = setup7Paths().json; end

    fprintf('\n============================================================\n');
    fprintf('ramseyStageTimer   --  THREE SHORT REAL RUNS, about 60 s\n');
    fprintf('============================================================\n');
    fprintf(['This runs the experiment. It has to: the per-tau cost is pg.changeSequence x2 +\n' ...
             'assignSequence + awg.Run + pg.run + the digitizer read, and only awg.Run can be\n' ...
             'measured from outside. ramseySpeedProbe timed that one and I wrongly called it\n' ...
             'the whole cost.\n']);

    n0 = opts.nTau0; dn = opts.dnTau;
    rL = opts.repeatsLow; rH = opts.repeatsHigh;

    fprintf('\n  M1  %2d taus, %5d repeats\n', n0,    rL);
    fprintf('  M2  %2d taus, %5d repeats      (M2-M1 gives the per-SHOT cost)\n', n0,    rH);
    fprintf('  M3  %2d taus, %5d repeats      (M3-M1 gives the per-TAU cost)\n', n0+dn, rL);

    t1 = runOne(n0,    rL, opts, 'M1');
    t2 = runOne(n0,    rH, opts, 'M2');
    t3 = runOne(n0+dn, rL, opts, 'M3');

    out = struct();
    out.measured = [n0 rL t1; n0 rH t2; n0+dn rL t3];

    avg = opts.averages;
    % C from M2 - M1: the only thing that changed is repeats.
    C = (t2 - t1) / max(n0*avg*(rH - rL), eps);              % seconds per shot
    % B from M3 - M1: dn extra taus, each costing B + rL*C. The extra upload for dn more
    % waveforms is in here too; it is small at these sizes and is reported separately below.
    B = (t3 - t1) / max(dn*avg, eps) - rL*C;

    out.perTauSeconds = B;
    out.perShotUs     = 1e6*C;

    fprintf('\n--- the decomposition ---\n');
    fprintf('  M1 %6.1f s     M2 %6.1f s     M3 %6.1f s\n', t1, t2, t3);
    fprintf('  per SHOT   %8.1f us   (the sequence duration, as the run actually pays it)\n', 1e6*C);
    fprintf('  per TAU    %8.3f s    <-- EVERYTHING walk removes\n', B);
    if B > 0
        fprintf('  ratio      %8.0f    shots'' worth of time in one tau''s overhead\n', B/max(C,eps));
    end
    if B <= 0
        fprintf(['\n  The per-tau cost came out non-positive, which means the two runs were too\n' ...
                 '  close together to separate. Raise dnTau and repeat.\n']);
        return
    end

    % ---- what this means for a full run ----
    J = readJsonSafe(opts.jsonLocation);
    step  = min(1/(6*(abs(opts.detuneMHz) + 2.16)), 0.075);
    nFull = numel(0.005:step:opts.tauMax);
    nNow  = numel(0.005:0.015:opts.tauMax);
    avgF  = 2; repF = 3000;

    classicNow  = nNow  * avgF * (B + repF*C);
    classicFast = nFull * avgF * (B + repF*C);

    % Walk: one pg.run covers (nTau + warmup) x P blocks, so the per-tau overhead is paid
    % once per RUN rather than once per tau.
    maxScans = 32000;      % measured by ramseySpeedProbe
    P = max(1, min(floor(maxScans/(4*(nFull+1))), repF));
    nRuns = ceil(repF/P) * avgF;
    walkFast = nRuns * (B + P*nFull*C);

    out.projection = struct('classicNow', classicNow, 'classicFast', classicFast, ...
                            'walkFast', walkFast, 'nFull', nFull, 'passesPerRun', P);

    fprintf('\n--- projected full runs, from the MEASURED per-tau cost ---\n');
    fprintf('  %-46s %7.0f s  (%.1f min)\n', ...
            sprintf('classic, %d taus (what you run now)', nNow), classicNow, classicNow/60);
    fprintf('  %-46s %7.0f s  (%.1f min)   %.1fx\n', ...
            sprintf('classic, %d taus (on resonance)', nFull), classicFast, classicFast/60, ...
            classicNow/max(classicFast,eps));
    fprintf('  %-46s %7.0f s  (%.1f min)   %.1fx\n', ...
            sprintf('WALK, %d taus, passesPerRun %d', nFull, P), walkFast, walkFast/60, ...
            classicNow/max(walkFast,eps));

    fprintf('\n============================================================\n');
    fprintf('VERDICT\n');
    fprintf('============================================================\n');
    if B > 20*C*repF/nFull || classicFast/max(walkFast,eps) > 2
        fprintf(['\nBUILD THE WALK. The per-tau cost is %.3f s, which is %.0fx the time it takes\n' ...
                 'to collect one shot -- so the run is dominated by the overhead of STARTING a\n' ...
                 'tau point, not by measuring it. That is precisely the term walk deletes, and\n' ...
                 'it is why walk won for Rabi.\n'], B, B/max(C,eps));
        fprintf('\nOn resonance with walk: %.0f s against the %.0f s you run now, %.1fx.\n', ...
                walkFast, classicNow, classicNow/max(walkFast,eps));
        fprintf('passesPerRun comes out %d from the 32000-scan ceiling and 4 windows per repeat.\n', P);
    else
        fprintf(['\nThe per-tau cost is %.3f s, only %.0fx one shot, so the run really is\n' ...
                 'photon-dominated and walk would buy little. Cut tau points and shorten the\n' ...
                 'initialisation instead.\n'], B, B/max(C,eps));
    end

    fprintf(['\nAND COMPARE WITH ramseySpeedProbe: it measured %.3f s for awg.Run alone. This\n' ...
             'says the real per-tau cost is %.3f s, so awg.Run is %.0f%% of it and the other\n' ...
             '%.0f%% is the pulse generator, the sequence assign and the digitizer -- none of\n' ...
             'which that probe could see.\n'], opts.armSeconds, B, ...
            100*opts.armSeconds/max(B,eps), 100*(1-opts.armSeconds/max(B,eps)));
    fprintf('\n============================================================\n\n');
end

% ======================================================================

function seconds = runOne(nTau, repeats, opts, tag)
% One short Ramsey, timed end to end. FRESH object every time: tau has no set method, so
% assigning it does not raise changeFlag and a reused object would silently keep the old
% waveforms.
    fprintf('\n--- %s: %d taus, %d repeats ---\n', tag, nTau, repeats);
    r = AWGExpRamsey;
    r.loadJsonParams('all', opts.jsonLocation);
    r.nChannels = 2;
    r.amplitude = opts.amplitude + [0 0];
    r.frequency = r.frequency + opts.detuneMHz;
    if ~isempty(opts.initUs); r.laserInitializationDuration = opts.initUs; end
    r.repeats          = repeats;
    r.averages         = opts.averages;
    r.constantTime     = 1;
    r.balancedSequence = 1;
    r.smallDelay       = 0.5;
    r.isTracking       = 0;
    r.shouldAutosave   = false;
    % Span the real range so the waveform LENGTHS are representative -- the upload is part of
    % what is being separated, and short taus would understate it.
    r.tau = linspace(0.005, opts.tauMax, nTau);

    % Quiet the per-waveform progress: three runs of it is noise here.
    styleWas = r.AWG.progressStyle;
    r.AWG.progressStyle = 'quiet';
    t0 = tic;
    r.run;
    seconds = toc(t0);
    r.AWG.progressStyle = styleWas;
    fprintf('  %.1f s  (%.3f s per tau, all costs included)\n', seconds, seconds/nTau);
end

function J = readJsonSafe(f)
    J = struct();
    if isempty(f) || ~exist(f, 'file'); return; end
    try; J = jsondecode(fileread(f)); catch; end
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
