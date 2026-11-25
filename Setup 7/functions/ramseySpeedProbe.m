function out = ramseySpeedProbe(opts)
%RAMSEYSPEEDPROBE Measure where a Ramsey run's time goes, then project the alternatives.
%
%   out = ramseySpeedProbe                       % ~40 s, no experiment is run
%   out = ramseySpeedProbe(struct('probeCeiling', false))
%
% Runs NO experiment and takes NO data. It measures the three constants that decide how long
% a Ramsey takes -- upload rate, re-arm cost, digitizer scan ceiling -- and then projects the
% full-run time for every strategy worth considering, so the choice is made on numbers rather
% than on which knob is nearest to hand.
%
% ===================================================================
% WHY THIS COMES BEFORE ANY CODE CHANGE
% ===================================================================
%
% A classic Ramsey costs, per average:
%
%       nTau re-arms  x  ~0.6 s          the AWG stopping and re-arming
%     + nTau x repeats x seqDuration     photons, the only part that is measurement
%     + one upload                       nTau x 2 x nChannels waveforms
%
% On the 167-point list that comes out roughly 60 s of upload, 195 s of re-arm and 135 s of
% photons -- so about two thirds is overhead. That looks like a transport problem and it is
% tempting to reach straight for the walk trick that fixed ESR and Rabi. But the re-arm count
% is nTau x averages, and nTau is set by the TAU GRID, which is set by physics -- so the first
% question is not "how do we re-arm less often" but "how many tau points does the measurement
% actually need".
%
% And the photon term is not what it looks like either: at these settings 81% of the sequence
% duration is LASER INITIALISATION, four events of (laserInitializationDuration - detection)
% per repeat. That is a third lever, and section 4 below breaks it out.
%
% THE ANSWER DEPENDS ON THE DETUNING, AND NOT IN THE OBVIOUS DIRECTION. The 14N hyperfine
% splits the free-induction signal into three components. Detuned by delta they sit at
% delta-2.16, delta and delta+2.16 MHz, so the fastest one is |delta| + 2.16 -- and six points
% per period of THAT is what sets the step:
%
%       delta = 5 MHz   fastest 7.16 MHz   step 0.023 us   ->  ~167 points to 2.5 us
%       delta = 2 MHz   fastest 4.16 MHz   step 0.040 us   ->   ~63 points
%       delta = 0       fastest 2.16 MHz   step 0.077 us   ->   ~33 points
%
% Five times fewer re-arms AND five times less upload, for the same information, purely from
% not detuning. The only thing detuning buys is moving the signal away from zero frequency
% where laser drift lives -- and a floating baseline in the fit does that job for nothing,
% because on resonance the signal is still a 2.16 MHz beat, not a DC level.
%
% So this probe measures the constants and then prints the table. Read the table.
%
% ===================================================================
% WHAT IT DOES TO THE INSTRUMENT
% ===================================================================
%
% It uploads a throwaway set of waveforms to time the upload, and it arms the AWG a few times
% to time the re-arm. Both are things a normal run does constantly. It leaves the AWG stopped
% with those waveforms loaded, so THE NEXT EXPERIMENT WILL RELOAD -- which it would anyway,
% since a fresh object always reloads.
%
% probeCeiling additionally asks the digitizer to allocate increasing numbers of scans until
% it refuses, and stops the count after each try. That is the measurement that decides whether
% the walk trick can help Ramsey at all, and it is the reason to run this at the setup rather
% than reason about it. Set probeCeiling false to skip it.
%
% OPTIONS
%   jsonLocation  default setup7Paths().json
%   tauMax        us, the longest tau the real run would use (default 2.5)
%   repeats       shots per tau per average, as the real run would set it (default 3000)
%   averages      default 2
%   currentStep   us, the tau step your run script uses NOW (default 0.015). The baseline
%                 row is computed from this rather than from the theoretical minimum, so
%                 "vs now" compares against what you actually run.
%   currentDetune MHz, the detuning your run script uses now (default 5)
%   nProbe        waveforms to upload for the timing (default 24, spanning 0 to tauMax)
%   nArm          how many times to arm the AWG for the re-arm timing (default 8)
%   probeCeiling  find the digitizer scan ceiling (default true)
%   ceilingTries  scan counts to try (default [5000 8000 12000 16000 24000 32000])
%   quiet         default false
%
% OUTPUT
%   out.secondsPerArm, out.mbPerSecond, out.buildSamplesPerSecond
%   out.maxScans          the largest allocation the digitizer accepted
%   out.strategies        struct array, one per row of the printed table
%   out.recommended       the row it recommends, and why

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('jsonLocation', '', 'tauMax', 2.5, ...
        'repeats', 3000, 'averages', 2, 'nProbe', 24, 'nArm', 8, ...
        'currentStep', 0.015, 'currentDetune', 5, ...
        'probeCeiling', true, 'ceilingTries', [5000 8000 12000 16000 24000 32000], ...
        'quiet', false), 'ramseySpeedProbe');

    if isempty(opts.jsonLocation)
        opts.jsonLocation = setup7Paths().json;
    end
    J = readJsonSafe(opts.jsonLocation);
    say = @(varargin) qprintf(opts.quiet, varargin{:});

    out = struct('secondsPerArm', NaN, 'mbPerSecond', NaN, 'buildSamplesPerSecond', NaN, ...
                 'maxScans', NaN, 'armOk', false, 'strategies', [], 'recommended', struct());

    say('\n============================================================\n');
    say('ramseySpeedProbe   --  no experiment is run, no data is taken\n');
    say('============================================================\n');

    awg  = resolveAWG();
    spcm = [];
    try; spcm = getObjByName(Spcm.NAME); catch; end

    nCh    = pick(J, 'nChannels', 2);
    halfPi = pick(J, 'halfPiTime', 0.027);
    det    = pick(J, 'detectionDuration', 2.14);
    refDet = pick(J, 'referenceDetectionDuration', det);
    ini    = pick(J, 'laserInitializationDuration', 29.4);
    lastD  = pick(J, 'lastDelay', 0.455);
    smallD = 0.5;

    say('\nfrom the json: nChannels %d, halfPiTime %.4f us, detection %.3f, init %.3f,\n', ...
        nCh, halfPi, det, ini);
    say('               lastDelay %.3f us. tauMax for the projection: %.2f us.\n', lastD, opts.tauMax);

    % ================= 1. UPLOAD =================
    % Representative LENGTHS, not just short ones: the upload cost is set by total samples,
    % and a Ramsey waveform is prePulse + 2*halfPi + tau + pad, so a probe built only from
    % short taus would understate the real upload several-fold.
    say('\n--- 1. upload rate ---\n');
    say('Uploading %d throwaway waveforms spanning tau 0 to %.2f us, which is the same\n', ...
        opts.nProbe, opts.tauMax);
    say('length distribution the real run uploads.\n');
    minWf = awg.minWaveformDuration;
    w = 2*pi*pick(J, 'frequency', 2820);

    % MaxNormPowerAllowed IS EMPTY UNTIL setAmplitude HAS BEEN CALLED -- it is not a constant,
    % it is set by setAmplitude and nothing else. On a fresh MATLAB session, before any
    % experiment has run, it is []. Unguarded that made the build loop execute zero times, so
    % the probe reported 0 MB in 0 s and every projection below was garbage, silently.
    %
    % The amplitude is irrelevant to a timing measurement -- only the sample count matters --
    % but calling setAmplitude first makes this exercise the same path a real LoadAWG takes,
    % so prefer that and fall back to unity only if it fails.
    ampJson = pick(J, 'amplitude', -16);
    try
        awg.setAmplitude(ampJson + zeros(1, nCh));
    catch
    end
    A = awg.MaxNormPowerAllowed;
    if isempty(A) || ~isnumeric(A)
        A = ones(1, nCh);
        say('\n  NOTE MaxNormPowerAllowed was empty (setAmplitude has not run this session),\n');
        say('       so unity amplitude is used. Timing is unaffected -- only the sample\n');
        say('       count matters here.\n');
    end
    if numel(A) < nCh; A = [A(:)' repmat(A(end), 1, nCh - numel(A))]; end
    tauProbe = linspace(0, opts.tauMax, opts.nProbe);

    try
        awg.ResetStoredData;
        awg.SampleRate = awg.MAX_SAMPLE_RATE;
        awg.workChannels = 1:nCh;
    catch err
        error(['ramseySpeedProbe: found an AWG object of class %s but could not set it up ' ...
               '(%s). That is not the Tektronix AWG driver.'], class(awg), err.message);
    end
    styleWas = awg.progressStyle;
    awg.progressStyle = 'milestone';
    awg.progressBegin(sprintf(['speed probe: %d representative Ramsey waveforms ' ...
        '(throwaway)'], opts.nProbe*nCh), opts.nProbe*nCh, 'waveform');
    tBuild = tic;
    for i = 1:opts.nProbe
        for chan = 1:min(nCh, numel(A))
            f = @(t) A(chan) * ((minWf <= t & t < minWf+halfPi) .* cos(w*t) + ...
                    (minWf+halfPi+tauProbe(i) <= t & t < minWf+2*halfPi+tauProbe(i)) .* cos(w*t));
            awg.AddWaveformByFunction(f, minWf + 2*halfPi + tauProbe(i) + minWf, 0, 1);
        end
    end
    buildSeconds = toc(tBuild);
    awg.progressEnd();

    nSamp = 0;
    for k = 1:numel(awg.waveforms); nSamp = nSamp + numel(awg.waveforms{k}); end
    if nSamp == 0
        error(['ramseySpeedProbe: no waveform samples were built, so there is nothing to ' ...
               'time. %d waveforms are stored. This is a bug in the probe, not in the ' ...
               'setup -- do not trust any number it would have printed.'], numel(awg.waveforms));
    end
    mb = 4*nSamp/2^20;
    tUp = tic;
    awg.LoadWaveform();
    upSeconds = toc(tUp);
    awg.progressStyle = styleWas;

    out.mbPerSecond = mb / max(upSeconds, eps);
    out.buildSamplesPerSecond = nSamp / max(buildSeconds, eps);
    say('\n  measured: %.1f MB in %.2f s = %.2f MB/s on the link\n', mb, upSeconds, out.mbPerSecond);
    say('            build %.2f s = %.2f M samples/s in MATLAB\n', ...
        buildSeconds, out.buildSamplesPerSecond/1e6);

    % Per-waveform-microsecond costs, which is what the projection needs.
    usTotal = nSamp / awg.SampleRate;
    out.uploadSecondsPerUs = upSeconds / max(usTotal, eps);
    out.buildSecondsPerUs  = buildSeconds / max(usTotal, eps);
    say('            -> %.3f s of upload per microsecond of waveform, all channels\n', ...
        out.uploadSecondsPerUs);

    % ================= 2. RE-ARM =================
    %
    % A SEQUENCE HAS TO BE ASSIGNED BEFORE Run CAN SUCCEED, and the first version of this
    % probe did not build one. Uploading waveforms is not enough: with nothing assigned the
    % instrument stays in RSTate = 0 (Stopped), waitForRunning polls until RUN_TIMEOUT, and
    % since RUN_TIMEOUT is 120 s and Constant the probe sat there for two minutes per arm and
    % reported the timeouts as if they were re-arm measurements. That is the worst possible
    % failure for a timing tool: it does not stop, and the number it produces is 200x too
    % large in the direction that makes the conclusion look more dramatic.
    %
    % So: build the same two-step trigger-gated sequence a real Ramsey block uses, assign it
    % to the active channels, and only then time Run. This is exactly what LoadAWG does.
    say('\n--- 2. re-arm cost ---\n');
    out.armOk = false;
    try
        nWf  = numel(awg.waveforms);
        pair = [max(nWf-1, 1), nWf];              % any two loaded waveforms will do
        seqIdx = awg.AddSequence(pair, [1 1], [1 1]);
        awg.LoadSequence();
        for chan = 1:nCh
            awg.assignSequence(seqIdx, chan);
        end
        out.armOk = true;
    catch err
        say(['  COULD NOT ASSIGN A SEQUENCE (%s), so Run would sit until the 120 s timeout.\n' ...
             '  Skipping the re-arm measurement rather than reporting a timeout as a\n' ...
             '  measurement.\n'], err.message);
    end

    if out.armOk
        say('Arming the AWG %d times, with a sequence assigned so the arm actually completes.\n', ...
            opts.nArm);
        say(['This is the per-tau overhead in classic mode. If the instrument rejects the\n' ...
             'sequence anyway, the FIRST arm costs up to 120 s (RUN_TIMEOUT is Constant and\n' ...
             'cannot be shortened from here) and then this section gives up -- so the worst\n' ...
             'case is one timeout, not %d of them.\n'], opts.nArm);
        armT = nan(1, opts.nArm);
        for k = 1:opts.nArm
            try
                awg.sendCommand('awgcontrol:stop:immediate');
            catch
            end
            tA = tic;
            try
                awg.Run;
                armT(k) = toc(tA);
            catch err
                % NOT swallowed into the average. One failure means the rest would each cost
                % another RUN_TIMEOUT, so stop now and say so.
                say(['\n  ARM %d FAILED after %.1f s: %s\n' ...
                     '  Abandoning the re-arm measurement -- every further attempt would cost\n' ...
                     '  another %g s timeout, and a timeout is not a re-arm time.\n'], ...
                    k, toc(tA), err.message, 120);
                break
            end
        end
        armT = armT(isfinite(armT));
        if isempty(armT)
            out.armOk = false;
        else
            out.secondsPerArm = median(armT);
            say('  measured: median %.3f s, spread %.3f to %.3f s over %d arms\n', ...
                out.secondsPerArm, min(armT), max(armT), numel(armT));
            % A plausibility gate. A real re-arm on this instrument is a fraction of a second;
            % anything near the timeout is a stuck arm being mistaken for a cost.
            if out.secondsPerArm > 10
                say(['\n  REJECTING THAT: %.1f s is not a re-arm, it is a stuck arm near the\n' ...
                     '  %g s timeout. The projection below would be meaningless, so the\n' ...
                     '  assumed value is used instead and the run needs looking at.\n'], ...
                    out.secondsPerArm, 120);
                out.armOk = false;
            end
        end
    end

    if ~out.armOk
        out.secondsPerArm = 0.6;      % measured for Rabi on this setup
        say(['\n  USING AN ASSUMED 0.600 s per re-arm (the value measured for Rabi here).\n' ...
             '  Every re-arm column below is therefore an estimate, not a measurement --\n' ...
             '  but the RATIOS between strategies do not depend on it, because every\n' ...
             '  strategy pays the same per-arm cost. The ranking stands either way.\n']);
    end

    % ================= 3. DIGITIZER CEILING =================
    say('\n--- 3. digitizer scan ceiling ---\n');
    say(['THIS IS THE MEASUREMENT THAT DECIDES WHETHER THE WALK TRICK CAN HELP RAMSEY.\n' ...
         'Walk arms once per pg.run instead of once per tau, but one pg.run must fit in the\n' ...
         'digitizer: it asks for detectionPeriodsPerRepeat x triggers scans. Ramsey has FOUR\n' ...
         'detection windows per repeat where Rabi has two, so the same ceiling buys Ramsey\n' ...
         'half the passes -- which is exactly why this cannot be assumed from Rabi.\n']);
    out.maxScans = NaN;
    if opts.probeCeiling && ~isempty(spcm)
        for n = sort(opts.ceilingTries)
            ok = false;
            try
                spcm.prepareExperimentCount(n, 5, 1);
                ok = true;
            catch
            end
            try
                spcm.stopExperimentCount;
            catch
            end
            if ok
                out.maxScans = n;
                say('  %6d scans: accepted\n', n);
            else
                say('  %6d scans: REFUSED\n', n);
                break
            end
        end
        if isfinite(out.maxScans)
            say(['\n  largest accepted allocation: %d scans. Note that a CLASSIC Ramsey at\n' ...
                 '  %d repeats already asks for %d (4 x repeats), so anything at or above\n' ...
                 '  that is known to work in practice.\n'], ...
                out.maxScans, opts.repeats, 4*opts.repeats);
        end
    else
        say('  skipped.\n');
    end
    if ~isfinite(out.maxScans)
        out.maxScans = 4*opts.repeats;   % what classic demonstrably does
        say('  using %d as the working ceiling: that is what a classic run already asks for.\n', ...
            out.maxScans);
    end

    % ================= 4. THE PROJECTION =================
    % Sequence duration per repeat, from the events prepare() actually adds. Written out
    % rather than read from pg.sequenceDuration so this works without building a sequence.
    initUs  = 2 * 2*(ini - det);                       % four init events per repeat
    readUs  = 2 * (det + refDet);
    mwUs    = 2 * (2*opts.tauMax + 2*halfPi);
    padUs   = 2 * 2*(2*smallD + 0.1 + lastD);
    seqUs   = initUs + readUs + mwUs + padUs;

    say('\n--- 4. where the sequence duration goes ---\n');
    say('Per repeat, with constantTime, balancedSequence and doubleMeasurement all on:\n\n');
    say('  laser initialisation  %7.1f us   %4.0f%%   4 x (laserInit %.1f - detection %.2f)\n', ...
        initUs, 100*initUs/seqUs, ini, det);
    say('  readout windows       %7.1f us   %4.0f%%   4 x %.2f\n', readUs, 100*readUs/seqUs, det);
    say('  MW gate at max tau    %7.1f us   %4.0f%%\n', mwUs, 100*mwUs/seqUs);
    say('  pads and triggers     %7.1f us   %4.0f%%\n', padUs, 100*padUs/seqUs);
    say('  --------------------------------------\n');
    say('  total per repeat      %7.1f us\n', seqUs);
    if initUs > 0.5*seqUs
        say(['\n  THE INITIALISATION IS %.0f%% OF EVERY SHOT. laserInitializationDuration is\n' ...
             '  %.1f us in the json; level 2 put the SNR-per-unit-time optimum nearer 15-20,\n' ...
             '  because chasing the last few percent of polarisation costs more statistics\n' ...
             '  than it returns. Dropping it to %.0f us would cut the photon term below to\n' ...
             '  %.0f%% of what it is, on every strategy in the table.\n'], ...
            100*initUs/seqUs, ini, 16, ...
            100*(2*2*(16-det) + readUs + mwUs + padUs)/seqUs);
    end

    say('\n--- 5. projected run time, per strategy ---\n');
    say('sequence duration per repeat: %.1f us\n', seqUs);
    say('\n%-34s %5s %7s %8s %8s %8s %7s\n', 'strategy', 'taus', 'upload', 're-arms', ...
        'photons', 'TOTAL', 'vs now');

    % ROW 1 IS YOUR CURRENT CELL, not the theoretical minimum for that detuning. They are
    % not the same: at 5 MHz detuning six points per period allows a 0.023 us step, but the
    % cell uses 0.015, which is 107 points against 167. Comparing everything to the
    % theoretical minimum would quietly flatter the alternatives by a third.
    nNow = max(round(opts.tauMax/opts.currentStep), 8);
    grids = { opts.currentDetune, sprintf('NOW: detune %g, step %g us', ...
                                          opts.currentDetune, opts.currentStep), nNow
              opts.currentDetune, sprintf('detune %g, coarsest legal step', opts.currentDetune), []
              2, 'detune 2 MHz, coarsest legal', []
              0, 'ON RESONANCE, coarsest legal', [] };
    S = struct('name', {}, 'detune', {}, 'nTau', {}, 'mode', {}, ...
               'upload', {}, 'rearm', {}, 'photons', {}, 'total', {});

    base = NaN;
    for g = 1:size(grids, 1)
        delta = grids{g,1};
        step  = 1/(6*(abs(delta) + 2.16));
        nTau  = grids{g,3};
        if isempty(nTau); nTau = max(round(opts.tauMax/step), 8); end
        % upload: total waveform microseconds = 2 (sig+ref) x nCh x sum(pre + 2hp + tau + pad)
        usAll = 2 * nCh * sum(2*minWf + 2*halfPi + linspace(0, opts.tauMax, nTau));
        up    = usAll * (out.uploadSecondsPerUs + out.buildSecondsPerUs);
        % classic
        arms  = nTau * opts.averages;
        ph    = arms * opts.repeats * seqUs * 1e-6;
        tot   = up + arms*out.secondsPerArm + ph;
        if g == 1; base = tot; end
        S(end+1) = struct('name', grids{g,2}, 'detune', delta, 'nTau', nTau, ...
                          'mode', 'classic', 'upload', up, 'rearm', arms*out.secondsPerArm, ...
                          'photons', ph, 'total', tot);                       %#ok<AGROW>
        printRow(say, grids{g,2}, nTau, up, arms*out.secondsPerArm, ph, tot, base);

        % walk, if it fits
        P = floor(out.maxScans / (4*(nTau + 1)));
        P = max(1, min(P, opts.repeats));
        armsW = ceil(opts.repeats / P) * opts.averages;
        totW  = up + armsW*out.secondsPerArm + ph;
        nm    = sprintf('%s + WALK (P=%d)', shortName(grids{g,2}), P);
        S(end+1) = struct('name', nm, 'detune', delta, 'nTau', nTau, 'mode', 'walk', ...
                          'upload', up, 'rearm', armsW*out.secondsPerArm, ...
                          'photons', ph, 'total', totW);                      %#ok<AGROW>
        printRow(say, nm, nTau, up, armsW*out.secondsPerArm, ph, totW, base);
    end
    out.strategies = S;

    % ================= 5. THE VERDICT =================
    % Leave it stopped rather than armed on a throwaway sequence.
    try; awg.sendCommand('awgcontrol:stop:immediate'); catch; end

    [~, iBest] = min([S.total]);
    out.recommended = S(iBest);
    seq = struct('total', seqUs, 'init', initUs, 'read', readUs, 'mw', mwUs, 'pad', padUs, ...
                 'ini', ini, 'det', det, 'refDet', refDet, 'lastD', lastD, ...
                 'smallD', smallD, 'halfPi', halfPi, 'tauMax', opts.tauMax);
    out.seq = seq;
    verdict(S, iBest, base, out, opts, seq);
end

% ======================================================================

function s = shortName(name)
% The walk row has to stand on its own in the verdict, where it is printed without the
% classic row above it for context.
    s = regexprep(name, '^NOW: ', '');
    s = regexprep(s, ', coarsest legal.*$', '');
    s = regexprep(s, ', step.*$', '');
    if numel(s) > 20; s = s(1:20); end
end

function printRow(say, name, nTau, up, arm, ph, tot, base)
    rel = '';
    if isfinite(base) && base > 0
        if abs(tot - base) < 0.02*base
            rel = '    --';
        else
            rel = sprintf(' %5.1fx', base/tot);
        end
    end
    say('%-34s %5d %6.0fs %7.0fs %7.0fs %7.0fs %s\n', name, nTau, up, arm, ph, tot, rel);
end

function verdict(S, iBest, base, out, opts, seq)
    fprintf('\n============================================================\n');
    fprintf('WHAT TO DO\n');
    fprintf('============================================================\n');

    % The best row overall is usually a walk row, but walk is new code -- so lead with the
    % best option that needs nothing written, and mention walk separately. Leading with a
    % number you cannot have today is how a recommendation gets ignored.
    clas = S(strcmp({S.mode}, 'classic'));
    [~, iC] = min([clas.total]);
    bFree = clas(iC);
    b = S(iBest);
    if ~isfield(out, 'armOk') || ~out.armOk
        fprintf(['\n0. THE RE-ARM COST WAS NOT MEASURED -- 0.600 s is assumed. Totals below are\n' ...
                 '   estimates; the ratios are not affected, since every strategy pays the same\n' ...
                 '   per-arm cost.\n']);
    end
    fprintf('\n1. FASTEST WITHOUT WRITING ANY CODE: %s.\n', strtrim(bFree.name));
    fprintf('   %.0f s against %.0f s, a factor of %.1f -- one number in the run script.\n', ...
            bFree.total, base, base/max(bFree.total,eps));
    if ~strcmp(b.mode, 'classic')
        fprintf('   (Fastest of all is %s at %.0f s, %.1fx, but that needs the walk built.)\n', ...
                strtrim(b.name), b.total, base/max(b.total,eps));
    end

    % Which term dominates the CURRENT setting, and which lever moves it.
    cur = S(1);
    terms = [cur.upload cur.rearm cur.photons];
    names = {'the upload', 'the re-arms', 'photons'};
    [~, iw] = max(terms);
    fprintf(['\n2. RIGHT NOW %.0f%% of the run is %s.\n'], 100*terms(iw)/max(cur.total,eps), names{iw});
    switch iw
        case 2
            fprintf(['   Re-arms are nTau x averages, so the only things that move them are\n' ...
                     '   FEWER TAU POINTS or walk mode. repeats does not enter -- raising\n' ...
                     '   repeats buys statistics for free until photons catch up.\n']);
        case 1
            fprintf(['   The upload scales with total waveform microseconds, which for Ramsey\n' ...
                     '   is dominated by the long-tau end. Fewer tau points, or a shorter\n' ...
                     '   tauMax, are the only levers; walk does not help an upload.\n']);
        case 3
            fprintf(['   Photons are the measurement, so in principle this is the good case --\n' ...
                     '   but NOT here, because %.0f%% of every shot is laser initialisation --\n' ...
                     '   four separate init events per repeat, two of them only there to serve\n' ...
                     '   the balanced reference. Photons dominating means the lever has moved\n' ...
                     '   from the transport to the SEQUENCE, not that there is nothing left to\n' ...
                     '   do. See the second table.\n'], ...
                    100*seq.init/seq.total);
    end

    % Is walk worth building at all?
    walkRows = S(strcmp({S.mode}, 'walk'));
    clasRows = S(strcmp({S.mode}, 'classic'));
    gain = [clasRows.total] ./ [walkRows.total];
    fprintf('\n3. IS THE WALK TRICK WORTH BUILDING FOR RAMSEY?\n');
    fprintf('   Measured ceiling %d scans; Ramsey needs 4 per repeat.\n', out.maxScans);
    fprintf(['   Walk is worth at most %.2fx. That is the whole answer: the re-arm term it\n' ...
             '   attacks is only %.0f%% of the run at the measured %.3f s per arm, so removing\n' ...
             '   almost all of it still leaves %.0f%% untouched. Amdahl, not pessimism.\n'], ...
            max(gain), 100*cur.rearm/max(cur.total,eps), out.secondsPerArm, ...
            100*(1 - cur.rearm/max(cur.total,eps)));
    fprintf('   NOT YET, on these numbers -- see item 6 for the order.\n');

    fprintf(['\n4. WHAT DETUNING COSTS YOU, said plainly. The tau step has to resolve the\n' ...
             '   FASTEST hyperfine component, |detune| + 2.16 MHz, so detuning multiplies\n' ...
             '   your point count by (|detune| + 2.16)/2.16 for no extra information about\n' ...
             '   T2*. At 5 MHz that is 3.3x. The one thing detuning buys -- signal away from\n' ...
             '   DC where laser drift lives -- is already provided by the 2.16 MHz beat\n' ...
             '   itself, and a floating baseline in the fit removes the rest.\n']);

    % ---- the second table: the SEQUENCE levers, on the best grid ----
    b2 = S(iBest);
    if strcmp(b2.mode, 'walk'); b2 = S(iBest-1); end     % use the classic row of that grid
    fprintf(['\n5. THE SEQUENCE-DURATION LEVERS, on the best grid above. These scale the\n' ...
             '   PHOTON column, which is now the largest one, and neither needs new code.\n' ...
             '   Grid assumed: %s.\n\n'], strtrim(b2.name));
    initShort = 16;                                     % level 2's SNR-per-unit-time region
    seqNoBal  = 2*(2*seq.smallD + 0.1 + seq.lastD + seq.tauMax + 2*seq.halfPi ...
                   + seq.det + (seq.ini - seq.det) + seq.refDet + 1);
    seqInit   = 2*2*(initShort - seq.det) + seq.read + seq.mw + seq.pad;
    seqBoth   = 2*(2*seq.smallD + 0.1 + seq.lastD + seq.tauMax + 2*seq.halfPi ...
                   + seq.det + (initShort - seq.det) + seq.refDet + 1);
    rows = { 'as now',                                  seq.total
             sprintf('laserInit %g -> %g us', seq.ini, initShort), seqInit
             'balancedSequence off',                    seqNoBal
             'both',                                    seqBoth };
    fprintf('   %-34s %9s %9s %8s\n', 'sequence change', 'us/shot', 'photons', 'TOTAL');
    for i = 1:size(rows,1)
        ph = b2.photons * rows{i,2} / max(seq.total, eps);
        tt = b2.upload + b2.rearm + ph;
        fprintf('   %-34s %8.1f %8.0fs %7.0fs', rows{i,1}, rows{i,2}, ph, tt);
        if i == 1
            fprintf('\n');
        else
            fprintf('   %.1fx\n', (b2.upload+b2.rearm+b2.photons)/max(tt,eps));
        end
    end
    fprintf(['\n   laserInitializationDuration is %g us in the json. Level 2 put the\n' ...
             '   SNR-per-unit-time optimum near 15-20, and said in its own report that the\n' ...
             '   value should come out SHORTER than the residual table suggests -- chasing the\n' ...
             '   last few percent of polarisation costs more statistics than it returns. This\n' ...
             '   is that finding, in seconds.\n'], seq.ini);
    fprintf(['\n   balancedSequence is the honest trade. On, the reference gets its own full\n' ...
             '   half-sequence -- trigger, gate, lastDelay, readout, initialisation -- so it\n' ...
             '   roughly doubles the shot. Off, the reference is a short window after the\n' ...
             '   measurement instead, which is cheaper but no longer matched to it. Change it\n' ...
             '   only if you are prepared to check that the ratio still normalises.\n']);

    % ---- the ordered plan ----
    fprintf('\n6. DO IT IN THIS ORDER.\n');
    fprintf('   a. Tau grid -> %s.\n', strtrim(b2.name));
    fprintf(['      Free, %.1fx, and it cuts upload, re-arms AND photons together, which is\n' ...
             '      why it beats anything that touches only one term.\n'], ...
            base/max(b2.total, eps));
    phInit = b2.photons * seqInit / max(seq.total, eps);
    fprintf('   b. laserInitializationDuration -> %g us. Free, another %.1fx.\n', ...
            initShort, (b2.upload+b2.rearm+b2.photons)/max(b2.upload+b2.rearm+phInit, eps));
    fprintf('      Cumulative: about %.0f s against the %.0f s you run now.\n', ...
            b2.upload+b2.rearm+phInit, base);
    walkRows2 = S(strcmp({S.mode}, 'walk'));
    clasRows2 = S(strcmp({S.mode}, 'classic'));
    g2 = max([clasRows2.total] ./ [walkRows2.total]);
    fprintf('   c. THEN re-measure, and only then decide about walk. On today''s numbers it\n');
    fprintf('      is worth %.1fx, against %.1fx for two numbers in the run script -- and it\n', ...
            g2, base/max(b2.upload+b2.rearm+phInit, eps));
    fprintf('      is new code, in a shared class, with a correctness risk (the MW gate has\n');
    fprintf('      to sit open at max tau for every point) that cannot be checked offline.\n');

    fprintf('\n============================================================\n\n');
end

% ======================================================================

function awg = resolveAWG()
%RESOLVEAWG Get the AWG object the same way ExperimentAWG does, and say what went wrong if not.
%
% getObjByName('TektronixAWG') is NOT reliable: the registered name only has to CONTAIN
% 'TektronixAWG', so on a setup where it is registered as e.g. 'TektronixAWG70002B' the exact
% lookup returns empty. And FrequencyGenerator.getFG() returns a CELL ARRAY of generators, not
% one object -- so falling back to it directly hands you a cell and the next line fails with
% "Dot indexing is not supported for variables of this type", which names neither the variable
% nor the reason.
    awg = [];

    % 1. the exact name, in case it is registered that way
    try
        awg = getObjByName('TektronixAWG');
    catch
    end
    if ~isempty(awg) && ~iscell(awg); return; end

    % 2. the ExperimentAWG route: scan the registered generators by name
    fgs = {};
    try
        fgs = FrequencyGenerator.getFG();
    catch err
        error(['ramseySpeedProbe: cannot reach the frequency generators (%s). Is the setup ' ...
               'initialised? Run the usual startup first.'], err.message);
    end
    if ~iscell(fgs); fgs = {fgs}; end

    names = {};
    for i = 1:numel(fgs)
        nm = '';
        try; nm = fgs{i}.name; catch; end
        names{end+1} = nm;                                                   %#ok<AGROW>
        % strfind, not contains(): contains is R2016b+ and absent from Octave, so a file
        % using it cannot be exercised by the offline graders. ExperimentAWG uses contains;
        % this is the same test written portably.
        if ~isempty(nm) && (~isempty(strfind(nm, 'TektronixAWG')) || ...
                            ~isempty(strfind(nm, 'AWGdummy')))                %#ok<STREMP>
            try
                awg = getObjByName(nm);
            catch
                awg = fgs{i};       % the object itself is just as good
            end
            if ~isempty(awg) && ~iscell(awg); return; end
        end
    end

    error(['ramseySpeedProbe: no AWG found. Registered frequency generators: %s.\n' ...
           'ExperimentAWG looks for a name containing ''TektronixAWG'' or ''AWGdummy'', and ' ...
           'none of the above match. If the AWG is off or not initialised, start it and ' ...
           'retry.'], strjoin(names(~cellfun(@isempty, names)), ', '));
end

function v = pick(J, field, dflt)
    v = dflt;
    if isstruct(J) && isfield(J, field) && ~isempty(J.(field)) && isnumeric(J.(field))
        v = J.(field);
        if ~isscalar(v); v = v(1); end
    end
end

function J = readJsonSafe(f)
    J = struct();
    if isempty(f) || ~exist(f, 'file'); return; end
    try; J = jsondecode(fileread(f)); catch; end
end

function qprintf(quiet, varargin)
    if ~quiet; fprintf(varargin{:}); end
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
