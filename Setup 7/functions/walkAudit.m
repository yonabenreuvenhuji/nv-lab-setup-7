function out = walkAudit(opts)
%WALKAUDIT Has walk mode been costing ESR and Rabi contrast all along? ~6 min.
%
%   out = walkAudit
%   out = walkAudit(struct('esr', false))          % Rabi only, ~3.5 min
%   out = walkAudit(struct('rabi', false))         % ESR only,  ~2.5 min
%
% ===================================================================
% WHY THIS EXISTS
% ===================================================================
%
% ramseyWalkDiag returned something that is not explainable by bookkeeping. At
% passesPerRun = 1 the walk's error bars MATCHED classic (0.000236 vs 0.000216) and its
% contrast was HALVED (0.45x). Any mis-mapping, smearing or averaging of tau would have
% raised the scatter, because walkReduce derives the error bar from exactly that scatter
% across passes. Noise right and signal halved means walk is applying DIFFERENT PULSES,
% not mislabelling the same ones.
%
% Which raises the obvious question about every Rabi and every ESR already taken on the
% walk path: is the same thing quietly in that data?
%
% ===================================================================
% WHAT THE CODE PREDICTS, BEFORE ANY MEASUREMENT
% ===================================================================
%
% This is worth reading first, because the two experiments are NOT equally exposed and
% the prediction is checkable.
%
% ESR IS STRUCTURALLY IMMUNE. Grep AWGExpESR for changeSequence: there are no hits. The
% ESR sequence is
%
%       S.addEvent(TRIG_DURATION,                    [greenLaser, triggers, MW]);
%       S.addEvent(detectionDuration - TRIG_DURATION,[greenLaser, MW, detector]);
%       ... reference half, no MW ...
%
% and NOTHING in it depends on the swept variable. ESR sweeps FREQUENCY at fixed drive
% duration, so the pulse-generator sequence handed to the hardware is byte-identical on
% the classic and walk paths. The only thing walk changes is which waveform the sequencer
% plays and how the digitizer stream is parsed. There is no gate to get wrong.
%
% RABI IS EXPOSED, and its own source says where. AWGExpRabi.prepare:
%
%       S.addEvent(obj.tau(end), MW, 'MW')      % the gate, built at the LONGEST tau
%
% Classic then rewrites it per point -- pg.changeSequence('MW','duration',obj.tau(t)) at
% line ~1010 -- and pads lastDelay to hold the total duration constant. Walk cannot: the
% pulse generator fires one trigger per block and never changes between them, so the gate
% stays at max(tau) for every point. The class documents the assumption it is making and
% asks for exactly this test:
%
%       "What you are trusting instead is that the AWG really does idle at zero while the
%        switch is still open... it is a real assumption and it is worth one check: run a
%        known pi time in both modes and compare the contrast."
%
% So: if the fixed gate is the mechanism, ESR agrees and Rabi does not. If BOTH disagree,
% the mechanism is in the concatenated sequence itself and predates the gate question. If
% both agree, walk is sound on these two and Ramsey's loss is specific to Ramsey.
% Every one of those three outcomes is informative, which is what makes this worth 6 min.
%
% ===================================================================
% THE ONE CONTROL THAT MAKES THE NUMBERS MEAN ANYTHING
% ===================================================================
%
% A contrast comparison between two runs taken at different times is confoundable by
% drift, and on this setup drift is not small -- it is the reason esrFitRobust takes the
% LARGER of the propagated error and the across-average scatter. Whichever run happens to
% go second looks worse, and a single classic-then-walk pair cannot tell that apart from
% walk being at fault.
%
% So each experiment runs THREE times: classic, walk, classic again. That gives a
% measured classic-to-classic difference over the same interval, and the walk deficit is
% judged against it rather than against zero. If the two classics differ by as much as
% walk does, the answer is "drift, not walk" and no amount of staring at a single pair
% would have told you. Set repeatClassic = false to skip it, and then do not quote the
% result.
%
% ===================================================================
% HOW THE DIFFERENCE IS READ -- NOT BY COMPARING TWO FITS
% ===================================================================
%
% Comparing fitted contrasts throws away where the difference sits. Instead the residual
% (walk - classic) is projected onto the two shapes that the two candidate mechanisms
% produce, in one weighted least squares:
%
%       walk - classic  =  offset  +  k * (classic - mean classic)  +  m * gateTail
%
%   k   a pure CONTRAST SCALE. walk = classic scaled by (1 + k) about its own mean.
%       k = -0.55 is the halved-contrast signature Ramsey showed. This is what a
%       uniform loss of coherence or a diluted reference looks like, and it is
%       INDEPENDENT of tau.
%
%   m   a TILT along the excess gate time, gateTail = (max(tau) - tau)/max(tau),
%       mean-removed. This is the signature of the fixed gate specifically: the tail
%       with the switch open and the AWG idling is longest at SHORT tau and zero at
%       max(tau), so a leak during it must show up as a tau-dependent tilt and cannot
%       show up as a flat scale. For ESR the column is omitted -- there is no varying
%       gate to tilt along, which is the same statement as ESR being immune.
%
% k and m are orthogonalised against each other and against the constant, so each comes
% back with an honest standard error and one cannot absorb the other. That is the whole
% point: a flat contrast deficit and a gate leak are different physics and they are
% separable in a single 50-point residual. Comparing two fitted amplitudes conflates them.
%
% The offset is fitted and then ignored -- a common-mode shift in the ratio is the laser
% or the reference level, not the pulses.
%
% The arithmetic lives in walkAuditDecompose.m, on its own, because testWalkAudit.m grades
% it: a known scale, a known tilt, both together, a pure baseline shift that must NOT leak
% into either, and 400 null realisations to show the reported error bar matches the actual
% spread of the coefficient. Run testWalkAudit before this -- it takes a second and needs
% no hardware.
%
% ===================================================================
% THE GATE: IS THERE ANYTHING TO COMPARE?
% ===================================================================
%
% A contrast scale is a scale OF something, and the arithmetic does not know that. The
% first real run of this file returned "-0.789 +/- 0.099, 8.0 sigma, walk contrast is 21%
% of classic's" on an ESR sweep with NO RESONANCE IN IT -- a smooth monotonic ramp, the
% microwave chain's transmission envelope across 60 MHz. Every digit of that was correct
% arithmetic on a baseline, and it read as a catastrophic finding.
%
% So the contrast scale is now WITHHELD, in the console and on the figure, unless both:
%
%   * the peak departure from a quadratic baseline is at least 5 error bars. A quadratic
%     absorbs a smooth envelope almost entirely and leaves a resonance or an oscillation
%     nearly untouched, so this separates a line from an envelope rather than signal from
%     noise.
%   * the line fit on the classic reference CONVERGED. esrFitRobust and rabiFitRobust
%     carry their own depth and width checks, and on that run esrFitRobust said exactly the
%     right thing -- "the only line left does not survive its checks" -- while being
%     treated as a footnote. It is a veto now.
%
% A fitter that cannot be CALLED (missing from the path) is reported loudly and does not
% veto, because that is a statement about the installation and not about the data.
%
% Also reported: whether the structure REPRODUCES between the two classic runs. That
% distinguishes reproducible wiggles no line fit will take -- something is there and the
% readout is diluting it -- from non-reproducing wiggles, which are drift.
%
% HOW SENSITIVE IS THIS? At the default settings the error on k comes out around 0.01,
% so a CLEAN verdict here rules out a contrast change of more than about 3%, and the
% Ramsey-sized 55% deficit would appear at over 40 sigma. The exact number for the run you
% actually did is printed in the verdict, because it depends on tonight's photons.
%
% ===================================================================
% OPTIONS
% ===================================================================
%
%   esr, rabi          run each (default both true)
%   repeatClassic      the drift control described above (default true). Turning it off
%                      makes the run 1/3 shorter and the conclusion unquotable.
%   passesPerRun       the walk setting to AUDIT (default 25), which is what the Rabi and
%                      ESR walk cells actually use. The question is about data already
%                      taken, so the default is what you have been running -- not a
%                      setting chosen to make walk look good or bad.
%   alsoP1             add a second Rabi walk run at passesPerRun = 1 (default false).
%                      Ramsey's loss was already there at P = 1, so this separates
%                      "walk is wrong" from "many passes per arm is wrong".
%   nTau, tauMax, rabiRepeats, rabiAverages     Rabi run (51, 1.0, 500, 2)
%   nFreq, freqSpanMHz, esrRepeats, esrAverages ESR run  (41, 60, 500, 2)
%   frequency          ESR sweep centre in MHz. Default: the json's frequency.
%   esrDetectionUs     ESR detection window in MICROSECONDS (default 100).
%                      NOT the json's value, and not a short pulsed-experiment window.
%                      ESR here is CW -- the MW is on for the whole detection window and
%                      the sequence has no changeSequence in it -- so a short window
%                      collects almost no spin signal and what is left is the microwave
%                      chain's transmission envelope across the span. That is not a
%                      hypothesis: the first version of this file used 0.3 us and returned
%                      a smooth monotonic ramp across 60 MHz with no line in it, on a setup
%                      whose resonance had been calibrated to the middle of that span. 10 us
%                      was already established to come out "extremely weird" here and 100 us
%                      to work, which is why 100 is the default.
%   rabiDetectionUs    Rabi detection window in microseconds. Default [] = LEAVE THE JSON
%                      ALONE, because level 2 calibrated it and a pulsed readout window has
%                      an SNR optimum that this file has no business overriding.
%   esrCoilsOff        centre the ESR sweep on the zero-field 2870 MHz with a 120 MHz span
%                      (default false). Moves the SWEEP only -- turn the coils off
%                      yourself. Use it when the calibrated centre no longer has a line in
%                      it, which is what happens if the field is not what it was at
%                      calibration.
%   amplitude, nChannels, jsonLocation
%   plotIt             (default true)
%   selfTest           run the whole pipeline on SYNTHETIC traces with an injected -0.45
%                      contrast scale, touching no instruments (default false). Use it to
%                      check the reporting, the fits and the plot before spending six
%                      minutes of setup time: everything in this file except the four run
%                      calls is arithmetic and printing, and that is historically where
%                      these scripts die -- after the measurement, with the data in a local
%                      that goes out of scope when the error propagates. Under selfTest the
%                      verdict must read AFFECTED; if it reads CLEAN the decomposition is
%                      not reaching the verdict.
%
% Runs nothing but experiments, writes no json, saves no .mat, and leaves both classes
% on the CLASSIC path when it finishes -- including when it throws.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct( ...
        'esr', true, 'rabi', true, 'repeatClassic', true, 'passesPerRun', 25, ...
        'alsoP1', false, ...
        'nTau', 51, 'tauMax', 1.0, 'rabiRepeats', 500, 'rabiAverages', 2, ...
        'nFreq', 41, 'freqSpanMHz', 60, 'esrRepeats', 500, 'esrAverages', 2, ...
        'frequency', [], 'amplitude', -16, 'nChannels', 2, ...
        'esrDetectionUs', 100, 'rabiDetectionUs', [], 'esrCoilsOff', false, ...
        'jsonLocation', '', 'plotIt', true, 'selfTest', false), 'walkAudit');

    if isempty(opts.jsonLocation) && ~isOn(opts.selfTest)
        opts.jsonLocation = setup7Paths().json;
    end
    if ~opts.esr && ~opts.rabi
        error('walkAudit: both esr and rabi are false, so there is nothing to run.');
    end
    if opts.rabiAverages < 2 || opts.esrAverages < 2
        error(['walkAudit: averages must be at least 2. With a single average sterr is ' ...
               'not a standard error -- it is the within-average scatter with nothing to ' ...
               'compare it against -- and every sigma printed below would be meaningless.']);
    end

    banner(opts);

    out = struct('opts', opts, 'esr', [], 'rabi', [], 'when', nowText());

    % Both experiments get the same treatment. The order is ESR first because it is the
    % one predicted to be clean: if it comes back dirty, the Rabi runs are answering a
    % different question and you will want to know that before waiting for them.
    if opts.esr
        out.esr = auditOne('esr', opts);
    end
    if opts.rabi
        out.rabi = auditOne('rabi', opts);
    end

    if opts.plotIt
        try
            drawAudit(out, opts);
        catch err
            fprintf('\n(plot failed: %s -- the numbers above are unaffected)\n', err.message);
        end
    end

    verdict(out, opts);
end

% ======================================================================
% ONE EXPERIMENT, THREE OR FOUR RUNS
% ======================================================================

function A = auditOne(which, opts)
    A = struct('which', which, 'runs', [], 'ok', false, 'status', 'did-not-finish');

    switch which
        case 'esr'
            xLabelText = 'MHz';
            [x, label] = esrAxis(opts);
        case 'rabi'
            xLabelText = '\tau (\mus)';
            [x, label] = rabiAxis(opts);
    end
    A.x = x; A.xLabel = xLabelText;

    fprintf('\n============================================================\n');
    fprintf('%s\n', upper(label));
    fprintf('============================================================\n');

    % ---- the three runs, each grabbed BEFORE the next object is built ----
    % expGrab is not optional here. Experiment objects are singletons keyed by class
    % name, so constructing the second AWGExpRabi DELETES the first -- the variable
    % survives as a dead handle and the next field access is "Invalid or deleted
    % object". The data has to come out before the next constructor runs.
    if isOn(opts.selfTest)
        R = synthRuns(which, x, opts);
    else
        R = {};
        R{end+1} = runOne(which, x, false, NaN, opts, 'classic #1');
        R{end+1} = runOne(which, x, true,  opts.passesPerRun, opts, ...
                          sprintf('walk, P = %d', opts.passesPerRun));
        if opts.repeatClassic
            R{end+1} = runOne(which, x, false, NaN, opts, 'classic #2 (drift control)');
        end
        if strcmp(which, 'rabi') && opts.alsoP1
            R{end+1} = runOne(which, x, true, 1, opts, 'walk, P = 1');
        end
    end
    A.runs = R;

    % ---- align every run onto the same axis -------------------------------
    n = numel(x);
    for i = 1:numel(R)
        if isempty(R{i}.x); R{i}.x = x; end
        if numel(R{i}.value) ~= n
            fprintf(['\n*** the %s run returned %d points, expected %d. Look further up ' ...
                     'for the underlying error -- Experiment.run catches per-average ' ...
                     'failures and only warns. Aborting this experiment.\n'], ...
                    R{i}.tag, numel(R{i}.value), n);
            return
        end
    end

    % ---- the classic reference, and the drift it measures ----------------
    iC = find(cellfun(@(r) ~r.walk, R));
    iW = find(cellfun(@(r)  r.walk, R));
    C = combineClassics(R(iC));
    A.classic = C;
    A.drift = driftBetweenClassics(R(iC));
    reportDrift(A.drift, opts.repeatClassic);

    % ---- IS THERE ANYTHING TO COMPARE? -----------------------------------
    % This gate exists because the first real run of this file produced a confident
    % "-0.789 +/- 0.099, 8.0 sigma, walk contrast is 21% of classic's" on an ESR sweep that
    % contained NO RESONANCE AT ALL -- just a smooth monotonic ramp, the microwave chain's
    % transmission envelope across 60 MHz. A contrast scale is a scale OF something, and
    % with no line there is nothing for it to be a scale of: the number was real
    % arithmetic on a featureless trace, and it was 8 sigma because the ramp is large and
    % the error bars are small.
    %
    % The discriminator is a LINE versus an ENVELOPE, not signal versus noise. A quadratic
    % baseline absorbs a smooth ramp almost entirely and leaves a resonance or an
    % oscillation nearly untouched, so the peak residual after removing one, measured in
    % error bars, separates the two cleanly. walkAuditDecompose's rank check cannot: the
    % column is not degenerate, it is just meaningless.
    % The two classic runs differ only by noise, so their difference measures it -- and the
    % reference is their mean, whose noise variance is a quarter of that difference's.
    A.refNoiseVar = [];
    if numel(iC) >= 2
        A.refNoiseVar = var(R{iC(1)}.value(:) - R{iC(2)}.value(:)) / 4;
        fprintf('\n  MEASURED NOISE (from the two classic runs, no reference to sterr):\n');
        fprintf('  reference rms %.6f, against a median sterr of %.6f -- a factor of %.1f.\n', ...
                sqrt(A.refNoiseVar), median(C.sterr), sqrt(A.refNoiseVar)/max(median(C.sterr), eps));
    end

    A.feature = featureStrength(C.x, C.value, C.sterr, R(iC));

    % The authoritative half of the gate is THE FIT, not a heuristic. esrFitRobust and
    % rabiFitRobust each carry their own depth and width checks and each drop a component
    % that fits nothing -- and on the run that motivated this gate, esrFitRobust said
    % exactly the right thing ("the only line left does not survive its checks") while the
    % decomposition went ahead and printed a confident 8 sigma contrast deficit anyway. So
    % the fitter's refusal is promoted from a footnote to a veto. The heuristic stays as a
    % second opinion: it catches a pure envelope in one line and it works on Rabi too.
    A.fits = struct('classic', fitTrace(which, C.x, C.value, C.sterr, 'classic'), 'walk', []);
    A.feature.fitOk    = A.fits.classic.ok;
    A.feature.fitCalled = A.fits.classic.available;
    if A.fits.classic.available
        A.feature.resolved = A.feature.resolved && A.fits.classic.ok;
    end            % else: fall back to the heuristic alone, having said so loudly
    reportFeature(A.feature, which);

    % ---- gate tail column: Rabi only, and that is the point --------------
    if strcmp(which, 'rabi')
        g = (max(x) - x) / max(max(x), eps);
    else
        g = [];
    end

    % ---- each walk run against that reference ----------------------------
    A.walk = {};
    for j = 1:numel(iW)
        W = R{iW(j)};
        cmp = compareTo(C, W, g, A.refNoiseVar);
        cmp.tag = W.tag;
        cmp.seconds = W.seconds;
        cmp.feature = A.feature;
        A.walk{end+1} = cmp;
        reportCompare(cmp, A.drift, A.feature);
    end

    % ---- the fits, which are the physics in a form you can quote ---------
    % classic was already fitted above, because its success is part of the gate.
    A.fits.walk = fitTrace(which, R{iW(1)}.x, R{iW(1)}.value, R{iW(1)}.sterr, 'walk');
    reportFits(which, A.fits);

    A.ok = true;
    if A.feature.resolved
        A.status = 'measured';
    else
        A.status = 'no-feature';
    end
end

% ======================================================================
% IS THERE A LINE, OR ONLY AN ENVELOPE?
% ======================================================================

function F = featureStrength(x, y, e, Rc)
%FEATURESTRENGTH Peak departure from a smooth baseline, in error bars.
%
% .sigma      peak residual after removing a quadratic, over the median error bar
% .repro      correlation of those residuals between the two classic runs, when there are
%             two. A real line reproduces; a drift excursion does not. NaN with one run.
% .resolved   sigma >= 5, ANDed by the caller with the line fit having converged
%
% Five error bars is chosen so a real ESR dip clears it by a wide margin -- a 2% dip on
% 2.5e-4 error bars is 80 -- while a transmission envelope that a quadratic almost
% perfectly absorbs does not. .sigma alone is NOT sufficient and is not used alone: on the
% run that motivated this gate the trace was a ramp with reproducible 6-sigma wiggles in
% it, which passes a size test and is still not a line you can measure a contrast on. That
% is why the caller also requires the fitter to have converged.
    F = struct('sigma', NaN, 'resolved', false, 'ptp', NaN, 'medErr', NaN, ...
               'ptpRaw', NaN, 'absorbed', NaN, 'repro', NaN, 'fitOk', true, ...
               'fitCalled', true);
    x = x(:); y = y(:); e = e(:);
    ok = isfinite(x) & isfinite(y);
    if sum(ok) < 8; return; end
    x = x(ok); y = y(ok);
    if numel(e) == numel(ok); e = e(ok); else; e = []; end

    med = median(e(isfinite(e) & e > 0));
    if isempty(med) || ~isfinite(med) || med <= 0; med = std(diff(y))/1.41; end
    if ~isfinite(med) || med <= 0; return; end

    % Quadratic baseline in a centred, scaled variable, so the normal equations stay
    % conditioned whether x is 2790-2850 or 0.02-1.0.
    u = (x - mean(x)) / max(max(x) - min(x), eps);
    V = [ones(numel(u),1), u, u.^2];
    th = V \ y;
    resid = y - V*th;

    F.medErr = med;
    F.ptpRaw = max(y) - min(y);
    F.ptp    = max(resid) - min(resid);
    F.sigma  = max(abs(resid)) / med;
    F.absorbed = 1 - F.ptp / max(F.ptpRaw, eps);
    F.resolved = F.sigma >= 5;

    % Does the structure REPRODUCE between the two classic runs? Reported rather than
    % judged, because it answers a different and useful question: reproducible wiggles that
    % no line fit will take mean there IS something on the setup and the readout is
    % diluting it, whereas non-reproducing wiggles are drift and there is nothing to find.
    if nargin >= 4 && iscell(Rc) && numel(Rc) >= 2
        r1 = baselineResid(Rc{1}.x, Rc{1}.value);
        r2 = baselineResid(Rc{2}.x, Rc{2}.value);
        if numel(r1) == numel(r2) && numel(r1) > 4
            F.repro = corrCol(r1, r2);
        end
    end
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

function reportFeature(F, which)
    if ~isfinite(F.sigma); return; end
    fprintf('\n  FEATURE CHECK: peak departure from a smooth baseline is %.1f error bars\n', ...
            F.sigma);
    fprintf('  (a quadratic absorbs %.0f%% of the raw span, leaving %.4f of %.4f)\n', ...
            100*F.absorbed, F.ptp, F.ptpRaw);
    if isfinite(F.repro)
        fprintf('  that structure reproduces between the two classic runs at r = %+.2f\n', ...
                F.repro);
    end
    if isfield(F, 'fitCalled') && ~F.fitCalled
        fprintf('  line fit on the classic reference: NOT AVAILABLE, heuristic only\n');
    else
        fprintf('  line fit on the classic reference: %s\n', ...
                eitherText(F.fitOk, 'converged', '*** REFUSED (see its reason above)'));
    end
    if F.resolved
        fprintf('  -> there is a real feature here, so a contrast comparison means something.\n');
    else
        fprintf(['  -> *** NO RESOLVED FEATURE. A quadratic accounts for essentially all of\n' ...
                 '     this trace, which is what a transmission envelope looks like and NOT\n' ...
                 '     what a %s looks like. No contrast scale will be reported: a scale has\n' ...
                 '     to be a scale OF something.\n'], ...
                eitherText(strcmp(which, 'esr'), 'resonance', 'Rabi oscillation'));
    end
end

% ======================================================================
% RUNNING
% ======================================================================

function g = runOne(which, x, walk, P, opts, tag)
    fprintf('\n--- %s ---\n', tag);
    switch which
        case 'esr';  e = buildESR(x, walk, P, opts);
        case 'rabi'; e = buildRabi(x, walk, P, opts);
    end

    % Quiet the per-waveform progress. The audit prints its own timings and the load
    % chatter would bury them; progressStyle is restored immediately afterwards so this
    % does not silently change how the next real run reports itself.
    styleWas = '';
    try; styleWas = e.AWG.progressStyle; e.AWG.progressStyle = 'milestone'; catch; end
    try; e.AWG.resetCommandCount; catch; end

    % ID FORM. A bare tic sets one global stopwatch, so any tic inside e.run resets it
    % and the following toc measures from there instead of from here -- that is how an
    % earlier probe reported 1.9 s for a run that took two minutes.
    tId = tic;
    ranOk = true; why = '';
    try
        e.run;
    catch err
        ranOk = false; why = err.message;
        fprintf('  *** the run threw: %s\n', why);
    end
    secs = toc(tId);

    nCmd = NaN;
    try; nCmd = e.AWG.commandCount; catch; end
    try; e.AWG.progressStyle = styleWas; catch; end

    g = expGrab(e, tag);
    g.walk = walk; g.passesPerRun = P; g.seconds = secs; g.commands = nCmd;
    g.ranOk = ranOk; g.why = why;
    if isempty(g.x); g.x = x; end

    % Leave nothing on the walk path. If this function throws or the caller aborts, the
    % object is about to be deleted by the next constructor anyway, but a run that ends
    % here must not leave a configured experiment armed on an unproven path.
    try; e.useSequenceWalk = false; catch; end

    fprintf('  %.0f s', secs);
    if isfinite(nCmd); fprintf(', %d SCPI commands', nCmd); end
    if ~isempty(g.value)
        fprintf(', span %.4f, median sterr %.5f', ...
                max(g.value) - min(g.value), median(g.sterr));
    end
    fprintf('\n');
end

function R = synthRuns(which, x, opts)
%SYNTHRUNS Fabricated runs for selfTest. No instruments, no experiment objects.
%
% This exists because every part of this file EXCEPT the four run() calls is arithmetic and
% printing, and that part has historically been where these scripts die -- after the
% measurement, with the data in a local variable that goes out of scope when the error
% propagates. Six minutes of setup time then buys nothing. selfTest drives the whole
% pipeline -- combine, drift, decompose, both fitters, the report, the plot, the verdict --
% on synthetic traces with a KNOWN injected deficit, so a formatting bug or a field-name
% typo surfaces in one second at the desk instead of at the end of the real run.
%
% The walk trace has a deliberate -0.45 contrast scale and a small gate tilt for Rabi, so
% the verdict printed under selfTest should read AFFECTED for both. A selfTest that
% printed CLEAN would mean the decomposition is not wired into the verdict at all.
    n = numel(x);
    flat = ischar(opts.selfTest) && strcmp(opts.selfTest, 'nofeature');
    switch which
        case 'esr'
            f0 = median(x);
            if flat
                % A smooth monotonic ramp and nothing else -- the microwave chain's
                % transmission envelope, which is exactly what the first real run of this
                % file measured when the ESR detection window was 0.3 us. The feature check
                % must refuse to report a contrast scale on this, and selfTest('nofeature')
                % is how that refusal stays tested.
                base = 1 - 0.012*(x - min(x))/max(max(x) - min(x), eps);
            else
                % TWO nearly-equal dips, 33 MHz apart, deliberately. That is the real
                % spectrum this setup produced (2855 and 2888, 1.915% and 1.908%), and it is
                % the case that broke the fit comparison: ranking lines by depth made the two
                % runs report different lines and print a 45 sigma centre disagreement on
                % data that agreed perfectly. Making it the DEFAULT synthetic spectrum means
                % the frequency pairing is exercised on every selfTest rather than only when
                % someone remembers to ask for it.
                base = 1 - 0.0300 ./ (1 + ((x - (f0 - 16.5))/6).^2) ...
                         - 0.0299 ./ (1 + ((x - (f0 + 16.5))/6).^2);
            end
            noise = 8e-4;
        case 'rabi'
            base = 1 - 0.06 * exp(-x/1.2) .* (1 - cos(2*pi*5.5*x)) / 2;
            noise = 8e-4;
    end
    e = noise*ones(1, n);
    c1 = base + noise*synthNoise(n, 1);
    c2 = base + noise*synthNoise(n, 2) + 0.0004;      % a little drift, as there always is
    kInject = -0.45;
    w = mean(base) + (1 + kInject)*(base - mean(base)) + noise*synthNoise(n, 3);
    if strcmp(which, 'rabi')
        g = (max(x) - x)/max(x);
        w = w - 0.002*(g - mean(g));                  % an injected gate tilt
    end

    R = {synthRun(x, c1, e, false, NaN, 'classic #1'), ...
         synthRun(x, w, e, true, opts.passesPerRun, ...
                  sprintf('walk, P = %d', opts.passesPerRun))};
    if opts.repeatClassic
        R{end+1} = synthRun(x, c2, e, false, NaN, 'classic #2 (drift control)');
    end
    fprintf('\n  selfTest: synthetic traces with an INJECTED contrast scale of %+.2f.\n', kInject);
    fprintf('  The verdict below should say AFFECTED. If it says CLEAN, the decomposition\n');
    fprintf('  is not reaching the verdict and the real run would be worthless.\n');
end

function g = synthRun(x, v, e, walkFlag, P, tag)
% The same field set runOne returns, so nothing downstream can tell the difference. If a
% field is added to runOne's output it has to be added here too, or selfTest silently
% stops covering whatever uses it.
    g = struct('tag', tag, 'x', x, 'value', v, 'sterr', e, ...
               'value2', [], 'sterr2', [], 'signal', [], 'rawSterr', [], ...
               'walk', walkFlag, 'passesPerRun', P, 'seconds', 42, 'commands', 1234, ...
               'ranOk', true, 'why', '', 'settings', struct(), 'name', 'synthetic');
end

function r = synthNoise(n, seed)
% A written-out LCG plus Box-Muller. rng() does not exist in Octave, where this file is
% parse-checked, and randn('seed',k) is deprecated in MATLAB -- so neither is used, and
% selfTest is bit-reproducible on any release.
    m = 2^31 - 1;
    xx = mod(1103515245*(seed + 12345) + 12345, m);
    u = zeros(1, 2*n);
    for i = 1:2*n
        xx = mod(1103515245*xx + 12345, m);
        u(i) = (xx + 0.5)/m;
    end
    r = sqrt(-2*log(u(1:n))) .* cos(2*pi*u(n+1:2*n));
end

function e = buildESR(freqs, walk, P, opts)
    e = AWGExpESR;
    e.loadJsonParams('all', opts.jsonLocation);
    e.nChannels = opts.nChannels;                 % BEFORE amplitude, always
    e.amplitude = opts.amplitude + zeros(1, opts.nChannels);
    e.repeats  = opts.esrRepeats;
    e.averages = opts.esrAverages;
    e.isTracking = 0;
    e.shouldAutosave = false;
    % A LONG window, on purpose. See esrDetectionUs in the header: CW ESR with a short
    % detection window measures the microwave chain, not the spins.
    e.detectionDuration = opts.esrDetectionUs;
    e.referenceDetectionDuration = e.detectionDuration;
    e.frequency = freqs;
    e.useSequenceJump = false;                    % never both
    e.useSequenceWalk = walk;
    if walk
        e.passesPerRun       = P;
        e.walkWarmupTriggers = 1;
        e.resyncEveryPass    = true;
    end
end

function e = buildRabi(taus, walk, P, opts)
    e = AWGExpRabi;
    e.loadJsonParams('all', opts.jsonLocation);
    e.nChannels = opts.nChannels;                 % BEFORE amplitude, always
    e.amplitude = opts.amplitude + zeros(1, opts.nChannels);
    e.repeats  = opts.rabiRepeats;
    e.averages = opts.rabiAverages;
    e.isTracking = 0;
    e.balancedSequence = 1;
    e.shouldAutosave = false;
    % Only if explicitly asked. The json's detectionDuration is the calibrated one and
    % overriding it here would compare walk against classic at a readout window neither
    % was calibrated for.
    if ~isempty(opts.rabiDetectionUs)
        e.detectionDuration = opts.rabiDetectionUs;
        e.referenceDetectionDuration = e.detectionDuration;
    end

    % constantTime ON for BOTH paths, deliberately. Walk holds the sequence duration
    % fixed for free -- the gate does not move -- so a classic run with constantTime off
    % would have a laser duty cycle that changes across the sweep while walk's does not,
    % and that difference alone would show up as a tilt and be blamed on the gate.
    try; e.constantTime = 1; catch; end

    e.tau = taus;
    e.useSequenceWalk = walk;
    if walk
        e.passesPerRun       = P;
        e.walkWarmupTriggers = 1;
        e.resyncEveryPass    = true;
    end
end

function [x, label] = esrAxis(opts)
    f0 = opts.frequency;
    span = opts.freqSpanMHz;
    if opts.esrCoilsOff
        % Zero field: both ms = +/-1 branches sit on D = 2870 MHz, split only by strain,
        % so one wide sweep centred there finds the line without needing to know the field
        % at all. A wider span than the default because at zero field there is nothing to
        % have calibrated the centre against.
        if isempty(f0); f0 = 2870; end
        if isempty(opts.frequency) && span < 120; span = 120; end
        fprintf(['\nesrCoilsOff: centring on %.0f MHz with a %.0f MHz span. This moves the\n' ...
                 'SWEEP only -- turn the coils off yourself, the option cannot do it.\n'], ...
                f0, span);
    end
    if isempty(f0) && isOn(opts.selfTest); f0 = 2760; end
    if isempty(f0)
        % Read the centre the same way the experiment will, so the sweep is centred on
        % the CALIBRATION rather than on a number written down in this file. The throwaway
        % object is safe: it is an AWGExpESR, so it deletes any earlier AWGExpESR -- and
        % the first real run below constructs a fresh one anyway.
        try
            tmp = AWGExpESR;
            tmp.loadJsonParams('all', opts.jsonLocation);
            f0 = median(tmp.frequency(:));
        catch
            f0 = [];
        end
        if isempty(f0) || ~isfinite(f0) || f0 <= 0
            f0 = 2760;
            fprintf(['\nNOTE could not read a sweep centre from the json, so the ESR sweep is\n' ...
                     '     centred on %g MHz. If the resonance is not in the span the ESR half\n' ...
                     '     of this audit is meaningless -- pass frequency explicitly.\n'], f0);
        end
    end
    half = span/2;
    x = linspace(f0 - half, f0 + half, opts.nFreq);
    label = sprintf(['ESR: %d points, %.0f MHz span about %.1f MHz, repeats %d, ' ...
                     'averages %d, detection %g us'], ...
                    opts.nFreq, span, f0, opts.esrRepeats, opts.esrAverages, ...
                    opts.esrDetectionUs);
end

function [x, label] = rabiAxis(opts)
    x = roundTo(linspace(0.02, opts.tauMax, opts.nTau), 3);
    label = sprintf('RABI: %d taus 0.02 to %.2f us, repeats %d, averages %d', ...
                    opts.nTau, opts.tauMax, opts.rabiRepeats, opts.rabiAverages);
end

% ======================================================================
% COMBINING, DRIFT, COMPARING
% ======================================================================

function C = combineClassics(Rc)
% Inverse-variance average of the classic runs, which is the right reference for the walk
% run precisely BECAUSE it straddles it in time: the walk run sits between classic #1 and
% classic #2, so their mean is contemporaneous with it to first order in the drift.
    n = numel(Rc{1}.value);
    W = zeros(numel(Rc), n); V = zeros(numel(Rc), n);
    for i = 1:numel(Rc)
        e = Rc{i}.sterr(:)';
        e(~isfinite(e) | e <= 0) = Inf;
        W(i,:) = 1./e.^2;
        V(i,:) = Rc{i}.value(:)';
    end
    sw = sum(W, 1);
    C.value = sum(W .* V, 1) ./ max(sw, eps);
    C.sterr = 1 ./ sqrt(max(sw, eps));
    C.x = Rc{1}.x;
    C.n = numel(Rc);
    C.seconds = sum(cellfun(@(r) r.seconds, Rc));
end

function d = driftBetweenClassics(Rc)
    d = struct('measured', false, 'scale', NaN, 'scaleErr', NaN, 'sigmaDiff', NaN, ...
               'errScale', NaN);
    if numel(Rc) < 2; return; end
    a = Rc{1}; b = Rc{2};
    % The same decomposition used for walk, so the two numbers are directly comparable:
    % how much did the CONTRAST change between two identical classic runs?
    p = walkAuditDecompose(a.value, a.sterr, b.value, b.sterr, [], ...
                           var(a.value(:) - b.value(:))/2);
    d.measured  = true;
    d.scale     = p.scale;
    d.scaleErr  = p.scaleErr;
    sig = sqrt(a.sterr.^2 + b.sterr.^2); sig(sig <= 0) = eps;
    d.sigmaDiff = median(abs(a.value - b.value) ./ sig);
    % HOW WRONG ARE THE ERROR BARS? For two correct measurements of the same thing,
    % median|difference|/sigma is 0.674 -- the median of a folded standard normal. Two
    % IDENTICAL classic runs therefore measure the understatement factor directly, with no
    % model and no assumption about what the spectrum looks like. This is worth having
    % because on this setup the answer is not 1: ESR came back at 12.3, i.e. sterr
    % understated by a factor of 18, while Rabi came back at 0.49, i.e. honest.
    d.errScale = d.sigmaDiff / 0.674;
    d.secondsApart = b.seconds;
end

function cmp = compareTo(C, W, g, noiseVar)
    cmp = struct();
    cmp.metrics = walkCompareMetrics(C.value, C.sterr, W.value, W.sterr);
    % noiseVar is MEASURED from the two classic runs where there are two, so the dilution
    % correction does not lean on sterr -- which on ESR here runs 18x low.
    if ~exist('noiseVar', 'var'); noiseVar = []; end
    cmp.decomp  = walkAuditDecompose(C.value, C.sterr, W.value, W.sterr, g, noiseVar);
    cmp.isEsr = isempty(g);          % no gate-tail column <=> ESR
    cmp.x = C.x; cmp.classic = C.value; cmp.classicErr = C.sterr;
    cmp.walk = W.value; cmp.walkErr = W.sterr;
end

% ======================================================================
% REPORTING
% ======================================================================

function banner(opts)
    fprintf('\n');
    fprintf('############################################################\n');
    fprintf('# walkAudit -- has walk mode been costing ESR and Rabi\n');
    fprintf('#              contrast, the way it costs Ramsey contrast?\n');
    fprintf('############################################################\n');
    nRuns = 0;
    if opts.esr;  nRuns = nRuns + 2 + double(opts.repeatClassic); end
    if opts.rabi; nRuns = nRuns + 2 + double(opts.repeatClassic) + double(opts.alsoP1); end
    fprintf('%d runs. Nothing is written, nothing is saved, both classes are left on\n', nRuns);
    fprintf('the classic path at the end.\n\n');
    fprintf('What the code predicts before any measurement:\n');
    fprintf('  ESR  cannot be affected. Its pulse-generator sequence does not depend on\n');
    fprintf('       the swept variable at all -- no changeSequence anywhere in the class --\n');
    fprintf('       so classic and walk hand the hardware the same sequence.\n');
    fprintf('  RABI can be. Its gate is built at max(tau) and classic rewrites it per\n');
    fprintf('       point; walk cannot, so the switch stays open while the AWG idles.\n');
    if ~opts.repeatClassic
        fprintf('\n*** repeatClassic is OFF. Without a second classic run there is no drift\n');
        fprintf('    control, and a contrast difference cannot be attributed to walk.\n');
    end
end

function reportDrift(d, wanted)
    if ~d.measured
        if wanted
            fprintf('\n  (no drift control -- the second classic run did not complete)\n');
        end
        return
    end
    fprintf('\n  DRIFT CONTROL: two identical classic runs differ by a contrast scale of\n');
    fprintf('  %+.3f +/- %.3f, median %.2f sigma per point.\n', d.scale, d.scaleErr, d.sigmaDiff);
    fprintf('  Any walk deficit has to beat this to mean anything.\n');
    if isfinite(d.errScale)
        fprintf('  ERROR BARS: two identical runs should differ by 0.67 sigma per point, so\n');
        if d.errScale > 2
            fprintf(['  sterr on this experiment is UNDERSTATED BY A FACTOR OF %.0f. Every\n' ...
                     '  "N sigma" and every chi2 computed from it is wrong by that factor --\n' ...
                     '  including any you have quoted from earlier runs of this experiment.\n' ...
                     '  The contrast scale above survives it: walkAuditDecompose inflates its\n' ...
                     '  covariance by the reduced chi2, which is what that Birge scaling is for.\n'], ...
                    d.errScale);
        elseif d.errScale < 0.5
            fprintf('  sterr is CONSERVATIVE by about %.1fx here. Errors are safe, if generous.\n', ...
                    1/max(d.errScale, eps));
        else
            fprintf('  sterr is honest to within %.1fx. Sigmas below can be read at face value.\n', ...
                    max(d.errScale, 1/max(d.errScale, eps)));
        end
    end
end

function reportCompare(cmp, drift, F)
    m = cmp.metrics; p = cmp.decomp;
    fprintf('\n  --- %s vs classic ---\n', cmp.tag);
    fprintf('  correlation %+.4f against a noise-limited ceiling of %+.4f (ratio %.2f)\n', ...
            m.correlation, m.rMax, m.rRatio);
    fprintf('  best shift  %+d          median |diff|/sigma %.2f\n', m.bestShift, m.sigmaDiff);
    % scatterRatio is reported for ESR only. It compares the point-to-point step to the
    % error bar, and on a Rabi that oscillates several times across the sweep a large value
    % is the SIGNAL, not scrambling -- this run printed 13.8 on a Rabi that agrees with
    % classic to 0.9 sigma per point. Printing it there invites exactly the wrong reading.
    if isfield(cmp, 'isEsr') && cmp.isEsr
        fprintf('  scatter ratio %.2f', m.scatterRatio);
        if drift.measured && isfinite(drift.errScale) && drift.errScale > 2
            % Divided by the error-bar understatement the drift control measured. It is an
            % independent check on that factor: the first real multi-line run gave 27.8 raw
            % and 1.5 corrected, against the 1.4 a real spectrum gives -- so the two ways of
            % finding the same understatement agree.
            fprintf('   -> %.2f once the %.0fx error-bar understatement is taken out', ...
                    m.scatterRatio/drift.errScale, drift.errScale);
        end
        fprintf('   (~1.4 for a real spectrum)\n');
    end
    if ~p.ok
        fprintf('  decomposition failed: %s\n', p.why);
        return
    end
    if nargin >= 3 && isstruct(F) && isfinite(F.sigma) && ~F.resolved
        fprintf(['  CONTRAST SCALE  withheld. The feature check above found nothing for a\n' ...
                 '  contrast to be the contrast of, so the number this would print (%+.3f)\n' ...
                 '  is arithmetic on a baseline and would be believed. Fix the run first.\n'], ...
                p.scale);
        return
    end
    fprintf('  CONTRAST SCALE  %+.3f +/- %.3f', p.scale, p.scaleErr);
    nSig = abs(p.scale)/max(p.scaleErr, eps);
    fprintf('   (%.1f sigma from zero)\n', nSig);
    fprintf('      walk contrast is %.0f%% of classic''s.\n', 100*(1 + p.scale));
    if isfinite(p.tilt)
        fprintf('  GATE-TAIL TILT  %+.4f +/- %.4f   (%.1f sigma)\n', ...
                p.tilt, p.tiltErr, abs(p.tilt)/max(p.tiltErr, eps));
        fprintf('      a tilt along (max(tau) - tau)/max(tau), which is the ONLY shape the\n');
        fprintf('      fixed gate can produce. Flat means the gate is not leaking.\n');
    end
    fprintf('  reduced chi2 of the decomposition %.2f\n', p.chi2);
    if isfinite(p.attenuation)
        fprintf('  dilution: uncorrected %+.3f -> corrected %+.3f (attenuation %.2f, noise %s)\n', ...
                p.scaleRaw, p.scale, p.attenuation, p.noiseFrom);
        if ~isempty(p.why); fprintf('     %s\n', p.why); end
    end
    if p.chi2 > 5
        fprintf(['      *** chi2 that large means the difference is NOT a scale plus a tilt.\n' ...
                 '      Something else is going on and the two coefficients above are a bad\n' ...
                 '      summary of it -- look at the difference panel before quoting them.\n']);
    end
    if drift.measured && isfinite(drift.scale)
        rel = abs(p.scale) / max(abs(drift.scale), eps);
        fprintf('  walk deficit is %.1fx the classic-to-classic drift.\n', rel);
    end
end

function f = fitTrace(which, x, y, e, tag)
    % 'available' distinguishes "the fitter RAN and refused this data" from "the fitter
    % could not be called at all". Only the first is evidence about the data. Without the
    % distinction, a missing file or a path problem makes the audit announce that the run
    % measured nothing -- which is the same words for a completely different situation, and
    % it would send you to the setup to fix a resonance that was never the problem.
    f = struct('ok', false, 'available', true, 'tag', tag, 'why', '');
    try
        switch which
            case 'esr'
                o = esrFitRobust(struct('f', x, 'y', y, 'e', e, ...
                                        'plotIt', false, 'verbose', false));
                if o.ok
                    % KEEP EVERY LINE and let reportFits pair them BY FREQUENCY.
                    %
                    % This used to take the deepest line, which is a trap on a real
                    % spectrum: with two dips of 1.915% and 1.908% the ranking is decided by
                    % noise, and the first run of this on a genuine multi-line spectrum
                    % reported "centre 2855.094 vs 2887.897, 45.4 sigma" -- two DIFFERENT
                    % lines, 33 MHz apart, while the depths agreed to 0.007% and the widths
                    % to 0.7 sigma. Nothing was wrong except which line each fit was asked
                    % about. Depth-rank pairing is the same "compare two fits" fragility the
                    % decomposition exists to avoid, and it had been reintroduced here.
                    L = o.lines;
                    f.lines = L;
                    f.nLines = numel(L);
                    f.r2 = o.r2; f.ok = true;
                else
                    f.why = o.why;
                end
            case 'rabi'
                o = rabiFitRobust(struct('t', x, 'y', y, 'e', e, ...
                                         'plotIt', false, 'verbose', false));
                if o.ok
                    f.freq = o.f;               f.freqErr = o.fErr;
                    f.piTime = 1e3*o.pulses(2); f.piErr = 1e3*o.pulsesErr(2);
                    f.depth = 100*o.contrast;   f.depthErr = 100*o.contrastErr;
                    f.r2 = o.r2; f.chi2 = o.chi2Reduced; f.ok = true;
                else
                    f.why = o.why;
                end
        end
    catch err
        f.why = err.message;
        f.available = false;
    end
    if ~f.ok && f.available
        fprintf('  (%s fit did not converge: %s)\n', tag, f.why);
    elseif ~f.ok
        fprintf(['  (*** the %s fitter could not be CALLED at all: %s\n' ...
                 '   That is an environment problem, not a statement about the data, so it\n' ...
                 '   does not veto anything below -- but fix it, because the fit is half of\n' ...
                 '   the check on whether this run measured a line.)\n'], tag, f.why);
    end
end

function reportFits(which, F)
    if ~F.classic.ok || ~F.walk.ok; return; end
    if strcmp(which, 'esr')
        reportEsrLines(F);
        return
    end
    fprintf('\n  --- the fits, for quoting ---\n');
    fprintf('  %-16s %12s %12s %14s\n', 'quantity', 'classic', 'walk', 'difference');
    rows = {'Rabi f / MHz', 'freq',   'freqErr'; ...
            'pi time / ns', 'piTime', 'piErr'; ...
            'contrast / %', 'depth',  'depthErr'};
    for k = 1:size(rows,1)
        vC = F.classic.(rows{k,2}); eC = F.classic.(rows{k,3});
        vW = F.walk.(rows{k,2});    eW = F.walk.(rows{k,3});
        vC = vC(1); eC = eC(1); vW = vW(1); eW = eW(1);
        sd = sqrt(eC^2 + eW^2);
        fprintf('  %-16s %12.3f %12.3f   %+8.3f (%.1f sigma)\n', ...
                rows{k,1}, vC, vW, vW - vC, abs(vW-vC)/max(sd, eps));
    end
    fprintf('  R^2              %12.4f %12.4f\n', F.classic.r2, F.walk.r2);
end

function reportEsrLines(F)
% Pair classic lines to walk lines by NEAREST CENTRE, greedily, closest pair first. Then
% every row compares the same physical line in the two runs, which is the only comparison
% that means anything on a multi-line spectrum.
    A = F.classic.lines; B = F.walk.lines;
    if isempty(A) || isempty(B); return; end
    ca = [A.centre]; cb = [B.centre];
    tol = 8;                                  % MHz. Wider than any real run-to-run shift,
                                              % narrower than the smallest line spacing
                                              % this setup shows (33 MHz here).
    [ia, ib] = pairByCentre(ca, cb, tol);

    fprintf('\n  --- the lines, paired by frequency (%d classic, %d walk) ---\n', ...
            numel(A), numel(B));
    fprintf('  %10s %10s   %8s %8s %7s   %7s %7s\n', ...
            'classic', 'walk', 'depth C', 'depth W', 'sigma', 'FWHM C', 'FWHM W');
    for k = 1:numel(ia)
        a = A(ia(k)); b = B(ib(k));
        sd = sqrt(a.depthErr^2 + b.depthErr^2);
        fprintf('  %10.2f %10.2f   %7.3f%% %7.3f%% %7.1f   %7.2f %7.2f\n', ...
                a.centre, b.centre, 100*a.depth, 100*b.depth, ...
                abs(b.depth - a.depth)/max(sd, eps), a.fwhm, b.fwhm);
    end
    nUnA = numel(A) - numel(ia); nUnB = numel(B) - numel(ib);
    if nUnA > 0 || nUnB > 0
        fprintf(['  %d classic line(s) and %d walk line(s) had no partner within %g MHz.\n' ...
                 '  On a spectrum this dense that is usually the fitter splitting or merging\n' ...
                 '  a shoulder, not a line appearing or vanishing -- check the overlay panel.\n'], ...
                nUnA, nUnB, tol);
    end
    fprintf('  R^2              %12.4f %12.4f\n', F.classic.r2, F.walk.r2);
end

function [ia, ib] = pairByCentre(ca, cb, tol)
    ia = []; ib = [];
    D = abs(ca(:) - cb(:)');
    D(D > tol) = Inf;
    while true
        [m, idx] = min(D(:));
        if ~isfinite(m); break; end
        [r, c] = ind2sub(size(D), idx);
        ia(end+1) = r; ib(end+1) = c;          %#ok<AGROW>
        D(r, :) = Inf; D(:, c) = Inf;
    end
    [~, ord] = sort(ca(ia));
    ia = ia(ord); ib = ib(ord);
end

% ======================================================================
% THE VERDICT
% ======================================================================

function verdict(out, opts)
    fprintf('\n############################################################\n');
    fprintf('# VERDICT\n');
    fprintf('############################################################\n');

    [okE, sE, eE, dE] = pick(out.esr);
    [okR, sR, eR, dR] = pick(out.rabi);

    guiltyE = okE && abs(sE) > 3*eE && abs(sE) > 0.10 && bigger(sE, dE);
    guiltyR = okR && abs(sR) > 3*eR && abs(sR) > 0.10 && bigger(sR, dR);

    if ranButNoFeature(out.esr)
        noFeatureAdvice('ESR', out.esr);
    end
    if ranButNoFeature(out.rabi)
        noFeatureAdvice('RABI', out.rabi);
    end

    if okE
        fprintf('\nESR   contrast scale %+.3f +/- %.3f  -> %s\n', sE, eE, tick(~guiltyE));
    end
    if okR
        fprintf('RABI  contrast scale %+.3f +/- %.3f  -> %s\n', sR, eR, tick(~guiltyR));
        t = out.rabi.walk{1}.decomp;
        if isfinite(t.tilt)
            fprintf('      gate-tail tilt %+.4f +/- %.4f -> %s\n', t.tilt, t.tiltErr, ...
                    tick(abs(t.tilt) <= 3*t.tiltErr));
        end
    end

    % A "clean" verdict is only worth anything with a stated sensitivity: clean means
    % "no deficit larger than this", and without the number it means nothing at all.
    % These are the actual measured error bars, so they already include whatever the
    % photon budget and the drift did tonight.
    fprintf('\nSENSITIVITY -- what a CLEAN result here rules out:\n');
    if okE
        fprintf('  ESR   a contrast change larger than %.1f%% would have shown at 3 sigma.\n', ...
                100*3*eE);
    end
    if okR
        fprintf('  RABI  a contrast change larger than %.1f%% would have shown at 3 sigma.\n', ...
                100*3*eR);
    end
    fprintf('  For scale, the Ramsey deficit was 55%%.\n');

    fprintf('\n');
    if ~okE && okR
        % The commonest real outcome of a first run, and it still answers the question that
        % prompted the audit -- because Rabi is the EXPOSED one. ESR being unmeasured costs
        % the cross-check, not the headline.
        if ~guiltyR
            fprintf(['RABI WALK IS CLEAN, and Rabi is the half that was at risk.\n\n' ...
                     'Rabi is the experiment whose MW gate walk cannot rewrite, so it is the\n' ...
                     'one that could have shown the Ramsey failure -- and it does not. Walk\n' ...
                     'reproduces classic to %.0f%% of contrast (%.1f sigma from identical), the\n' ...
                     'correlation sits at the noise-limited ceiling, and the gate-tail tilt is\n' ...
                     'flat. THE RABI DATA YOU HAVE ALREADY TAKEN ON THE WALK PATH IS FINE, and\n' ...
                     'that includes the pi times in the json.\n\n' ...
                     'The flat tilt is the more interesting half of that. The tilt is the only\n' ...
                     'shape a gate left open over an idling AWG can produce, and at %.4f +/- %.4f\n' ...
                     'it is not there. So the leading explanation for Ramsey -- leakage during\n' ...
                     'the excess gate -- just got much weaker, on a %.2f us tail of the same\n' ...
                     'order as Ramsey''s 2.5 us. What Ramsey has and Rabi does not:\n' ...
                     '  * the state sits in a SUPERPOSITION during the tail, not along z. A\n' ...
                     '    tail that leaves a z-aligned state untouched still dephases a\n' ...
                     '    superposition, and that is a contrast-only loss.\n' ...
                     '  * four detection windows per repeat instead of two\n' ...
                     '  * a closing pi/2 whose PHASE matters, where Rabi has no closing pulse\n\n' ...
                     'The first of those is now the leading candidate, and ramseyGateTest --\n' ...
                     'which varies the excess gate and nothing else -- is the test that\n' ...
                     'separates it from the remaining gate story.\n\n' ...
                     'The ESR half measured nothing; see above for why and what to change.\n'], ...
                    100*(1 + sR), abs(sR)/max(eR, eps), ...
                    out.rabi.walk{1}.decomp.tilt, out.rabi.walk{1}.decomp.tiltErr, ...
                    max(out.rabi.x));
        else
            fprintf(['*** RABI WALK DATA IS AFFECTED (contrast %.0f%% of classic). The ESR half\n' ...
                     'measured nothing, so there is no cross-check tonight -- but Rabi is the\n' ...
                     'exposed half and it is the one that matters. Stop using walk for Rabi,\n' ...
                     'and re-run level 4.4 on the classic path.\n'], 100*(1 + sR));
        end
    elseif okE && ~okR
        fprintf(['Only the ESR half is judgeable tonight. Rabi is the half that was at risk --\n' ...
                 'it is the one with the gate walk cannot rewrite -- so re-run with\n' ...
                 'walkAudit(struct(''esr'', false)) before concluding anything about walk.\n']);
    elseif ~okE && ~okR
        fprintf(['NEITHER HALF MEASURED ANYTHING. Nothing about walk follows from this run.\n' ...
                 'Fix what the notes above say and run it again -- and note that the audit\n' ...
                 'declining to give you a number is the correct outcome here, not a failure:\n' ...
                 'the arithmetic would happily have printed a confident 8 sigma contrast\n' ...
                 'deficit measured on a baseline.\n']);
    elseif ~guiltyE && ~guiltyR
        fprintf(['WALK IS CLEAN ON ESR AND RABI. Neither shows a contrast deficit bigger\n' ...
                 'than three times its own error and bigger than the drift between two\n' ...
                 'identical classic runs. The ESR and Rabi data you have already taken on\n' ...
                 'the walk path does not need re-taking.\n\n' ...
                 'That also localises the Ramsey problem rather than leaving it open. Rabi\n' ...
                 'shares the fixed gate with Ramsey and is fine, so an open switch over an\n' ...
                 'idling AWG is NOT by itself enough to destroy contrast -- which means the\n' ...
                 'thing Ramsey has and Rabi does not is what to look at:\n' ...
                 '  * two free-evolution windows instead of one driven window, so the state\n' ...
                 '    is in a superposition during the tail rather than along z\n' ...
                 '  * four detection windows per repeat instead of two\n' ...
                 '  * a second pi/2 whose PHASE matters, where Rabi has no closing pulse\n' ...
                 'A superposition dephases during a tail that leaves a z-aligned state\n' ...
                 'untouched, so the first of those explains a contrast-only loss and the\n' ...
                 'Rabi result here is consistent with it.\n\n' ...
                 'NEXT: ramseyGateTest, which forces the excess gate to differ while holding\n' ...
                 'everything else fixed, is now the test worth running.\n']);
    elseif guiltyR && ~guiltyE
        fprintf(['*** RABI WALK DATA IS AFFECTED. ESR IS NOT.\n\n' ...
                 'That is exactly the split the code predicts: ESR hands the hardware the\n' ...
                 'same sequence on both paths, Rabi does not -- its gate is stuck at\n' ...
                 'max(tau) in walk mode.\n\n' ...
                 'WHAT TO DO:\n' ...
                 '  1. ESR walk data is fine. Keep it.\n' ...
                 '  2. Any Rabi taken with useSequenceWalk = true reads a contrast %.0f%% of\n' ...
                 '     the truth. A pi TIME from such a fit is much less affected than the\n' ...
                 '     contrast -- the fitted period barely moves under a uniform amplitude\n' ...
                 '     scale -- so calibrations that only used the pi time are probably\n' ...
                 '     usable. Check the tilt line above: if the tilt is significant the\n' ...
                 '     period IS affected and the pi times need re-measuring.\n' ...
                 '  3. Set useSequenceWalk = false in the Rabi cells of\n' ...
                 '     runScriptBasicExpAWG and in runRabiBlock, today.\n' ...
                 '  4. Re-run level 4.4 (calibPulseTimes) on the classic path and compare\n' ...
                 '     the pi time with what is in the json now. That is the one number\n' ...
                 '     everything downstream depends on.\n'], 100*(1+sR));
    elseif guiltyE && guiltyR
        fprintf(['*** BOTH ARE AFFECTED, AND THAT RULES OUT THE GATE.\n\n' ...
                 'ESR hands the hardware an identical pulse-generator sequence on both\n' ...
                 'paths -- there is no changeSequence in the class and nothing in its\n' ...
                 'sequence depends on frequency. So if ESR loses contrast in walk mode the\n' ...
                 'cause cannot be the MW gate, and it cannot be tau-dependent either. What\n' ...
                 'is left is common to all three classes: the concatenated sequence itself.\n\n' ...
                 'The candidates, in the order worth checking:\n' ...
                 '  * the AWG does not settle between trigger-gated steps, so a block plays\n' ...
                 '    a waveform that has not finished loading its start point\n' ...
                 '  * the digitizer segments are not aligned with the blocks, so every point\n' ...
                 '    is a fixed blend of two neighbours -- which dilutes contrast without\n' ...
                 '    adding scatter, and is the one bookkeeping error that CAN do that\n' ...
                 '  * the warm-up block is not enough after the idle gap\n\n' ...
                 'WHAT TO DO: stop using walk on all three classes until this is settled,\n' ...
                 'and re-take anything whose CONTRAST you have quoted. The blend hypothesis\n' ...
                 'is testable offline against a synthetic stream and that is cheaper than\n' ...
                 'another night on the setup.\n']);
    else
        fprintf(['*** ESR IS AFFECTED AND RABI IS NOT, WHICH IS BACKWARDS.\n\n' ...
                 'Rabi is the exposed one -- it has the fixed gate -- and ESR is the one that\n' ...
                 'hands the hardware an identical sequence on both paths. Getting this pair\n' ...
                 'of answers means the ESR result is more likely a fault in the run than in\n' ...
                 'walk: check the printed correlation and best shift for the ESR comparison,\n' ...
                 'and whether the sweep straddled the resonance at all. Re-run the ESR half\n' ...
                 'with a wider span before drawing any conclusion.\n']);
    end

    if okR && ~opts.alsoP1 && guiltyR
        fprintf(['\nAlso worth one extra run: walkAudit(struct(''esr'', false, ''alsoP1'', true))\n' ...
                 'adds a Rabi walk at passesPerRun = 1. If the deficit survives at P = 1 it is\n' ...
                 'the walk sequence; if it only appears at P = %d it is the many-passes-per-arm\n' ...
                 'path, which is a different and easier fix.\n'], opts.passesPerRun);
    end
    fprintf('\n############################################################\n\n');
end

function [ok, s, e, d] = pick(A)
% ok is false when there is nothing to judge -- either the runs did not complete, or the
% feature check found no line. Those are different from "clean", and collapsing them into
% clean is how a failed run gets reported as a passing audit.
    ok = false; s = NaN; e = NaN; d = NaN;
    if isempty(A) || ~isstruct(A) || ~A.ok || isempty(A.walk); return; end
    if isfield(A, 'feature') && isstruct(A.feature) && isfinite(A.feature.sigma) ...
            && ~A.feature.resolved
        return
    end
    p = A.walk{1}.decomp;
    if ~p.ok; return; end
    ok = true; s = p.scale; e = p.scaleErr;
    if A.drift.measured; d = A.drift.scale; end
end

function tf = ranButNoFeature(A)
    tf = ~isempty(A) && isstruct(A) && isfield(A, 'status') && ...
         strcmp(A.status, 'no-feature');
end

function tf = isOn(v)
% selfTest accepts true, or a string naming a variant ('nofeature'). Anything non-empty and
% not false means on.
    if ischar(v); tf = ~isempty(v); return; end
    tf = ~isempty(v) && any(logical(v(:)));
end

function s = eitherText(tf, a, b)
% Local, not BooleanHelper.ifTrueElse, so selfTest runs outside the lab codebase.
    if tf; s = a; else; s = b; end
end

function tf = bigger(s, d)
% Is the walk deficit bigger than the drift between two identical classic runs? With no
% drift control this returns true, so the guilty verdict can still be reached -- but the
% banner has already said the conclusion is unquotable in that case.
    if ~isfinite(d); tf = true; return; end
    tf = abs(s) > 2*abs(d);
end

function noFeatureAdvice(name, A)
    F = A.feature;
    fprintf('\n%s   NO RESOLVED FEATURE -- this half measured nothing.\n', name);
    fprintf('      Peak departure from a smooth baseline: %.1f error bars (needs 5).\n', F.sigma);
    fprintf('      A quadratic absorbed %.0f%% of the span, so what was measured is a\n', ...
            100*F.absorbed);
    fprintf('      transmission envelope, not a line.\n');
    if strcmp(name, 'ESR')
        fprintf(['\n      TWO THINGS TO CHANGE, IN THIS ORDER:\n' ...
                 '      1. detectionDuration. ESR here is CW -- the MW is on for the whole\n' ...
                 '         detection window -- so a short window collects almost no spin\n' ...
                 '         signal. 10 us was already known to come out weird on this setup\n' ...
                 '         and 100 us to work. walkAudit now defaults esrDetectionUs to 100;\n' ...
                 '         if you are reading this from a run that used less, that alone is\n' ...
                 '         enough to explain a featureless sweep.\n' ...
                 '      2. WHERE the resonance actually is. The sweep centre came from the\n' ...
                 '         json, i.e. from the last calibration. If the coils are off now, or\n' ...
                 '         set differently, the line is not there any more. The shape of the\n' ...
                 '         trace tells you which way: still falling at the top of the span\n' ...
                 '         means the line is ABOVE it, and with the coils off it goes to the\n' ...
                 '         zero-field 2870 MHz.\n' ...
                 '\n      COILS OFF, WIDE SWEEP AT ZERO FIELD:\n' ...
                 '         >> walkAudit(struct(''rabi'', false, ''esrCoilsOff'', true))\n' ...
                 '      which centres on 2870 MHz with a 120 MHz span. Turn the coils off\n' ...
                 '      first -- the option only moves the sweep, it cannot move the field.\n']);
    end
end

function s = tick(good)
    if good; s = 'CLEAN'; else; s = 'AFFECTED'; end
end

% ======================================================================
% THE PLOT
% ======================================================================

function drawAudit(out, opts)
% Two rows, one per experiment: the traces, then the difference with the decomposition
% drawn on it. The DIFFERENCE panel is the one that answers the question -- on the overlay
% panel a halved contrast and a perfect match both just look like "close enough", whereas
% in the difference a contrast deficit looks like the SIGNAL and a gate leak looks like a
% RAMP, and those are visually unmistakable.
    have = {};
    if ~isempty(out.esr)  && out.esr.ok;  have{end+1} = 'esr';  end
    if ~isempty(out.rabi) && out.rabi.ok; have{end+1} = 'rabi'; end
    if isempty(have); return; end

    % Deep blue against orange: the one categorical pair that stays separable under every
    % common colour vision deficiency, and in greyscale print the lightness differs too.
    COL = struct('classic', [0.11 0.35 0.60], ...
                 'walk',    [0.85 0.37 0.10], ...
                 'second',  [0.55 0.62 0.70], ...
                 'model',   [0.15 0.15 0.15], ...
                 'band',    [0.86 0.88 0.91]);

    figure('Color', 'w', 'Name', 'walkAudit', ...
           'Position', [60 60 1220 400*numel(have)]);
    for r = 1:numel(have)
        A = out.(have{r});
        cmp = A.walk{1};
        p = cmp.decomp;

        % ---- left: the two traces ------------------------------------
        subplot(numel(have), 2, 2*r-1); hold on; grid on; box off
        axStyle();
        H = []; L = {};
        iC = find(cellfun(@(q) ~q.walk, A.runs));
        if numel(iC) >= 2
            H(end+1) = errbar(A.runs{iC(2)}.x, A.runs{iC(2)}.value, ...
                              A.runs{iC(2)}.sterr, COL.second, 1.0, 4, 'classic #2');
            L{end+1} = 'classic #2';
        end
        H(end+1) = errbar(cmp.x, cmp.classic, cmp.classicErr, COL.classic, 2.0, 6, 'classic');
        L{end+1} = 'classic';
        H(end+1) = errbar(cmp.x, cmp.walk, cmp.walkErr, COL.walk, 2.0, 6, cmp.tag);
        L{end+1} = cmp.tag;
        xlabel(A.xLabel); ylabel('ratio');
        % The legend is built from EXPLICIT handles. Letting it read DisplayName off
        % whatever is in the axes picks up the extra line objects errorbar creates and
        % produces a legend with "data1", "data2", "data3" in it -- six entries for three
        % series, and the reader cannot tell which colour is which.
        headroom(0.24);          % make room for the legend BEFORE placing it
        legendFrom(H, L);
        title(sprintf('%s   classic vs walk', upper(have{r})), 'FontWeight', 'bold');
        xlim(padRange(cmp.x));

        % ---- right: the difference, which is the answer ---------------
        subplot(numel(have), 2, 2*r); hold on; grid on; box off
        axStyle();
        d = cmp.walk - cmp.classic;
        H = []; L = {};
        if isfield(p, 'sigma') && ~isempty(p.sigma) && numel(p.sigma) == numel(cmp.x)
            % the band the difference is allowed to wander in, drawn first so the points
            % sit on top of it
            hb = fill([cmp.x, fliplr(cmp.x)], [p.sigma, -fliplr(p.sigma)], COL.band, ...
                      'EdgeColor', 'none');
            H(end+1) = hb; L{end+1} = '\pm1\sigma';
        end
        plot(padRange(cmp.x), [0 0], '-', 'Color', [0.65 0.65 0.65]);
        H(end+1) = errbar(cmp.x, d, p.sigma, COL.walk, 1.6, 6, 'walk - classic');
        L{end+1} = 'walk - classic';
        if p.ok
            H(end+1) = plot(cmp.x, p.model, '-', 'Color', COL.model, 'LineWidth', 2.0);
            L{end+1} = 'fitted scale + tilt';
        end
        xlabel(A.xLabel); ylabel('walk - classic');
        headroom(0.26);
        legendFrom(H, L);
        xlim(padRange(cmp.x));

        % Values ON the graph. Two short lines rather than one long one: a single-line
        % title of this length runs off the panel and collides with its neighbour, and a
        % number the reader cannot see is the same as a number that was never printed.
        Ffeat = struct('sigma', NaN, 'resolved', true);
        if isfield(A, 'feature'); Ffeat = A.feature; end
        title(titleLines(p, A.drift, cmp.metrics, Ffeat), 'FontWeight', 'normal', ...
              'FontSize', 8);
    end
end

function T = titleLines(p, drift, m, F)
    if nargin >= 4 && isstruct(F) && isfinite(F.sigma) && ~F.resolved
        % Say it on the figure too. A withheld number in the console and a confident one in
        % the title is the worst of both, and the figure is what gets pasted into a message.
        T = {sprintf('NO RESOLVED FEATURE -- %.1f error bars from a smooth baseline, needs 5', ...
                     F.sigma), ...
             'contrast scale withheld: nothing here for a contrast to be the contrast of'};
        return
    end
    a = sprintf('contrast scale %+.3f \\pm %.3f (%.1f\\sigma)   |   walk contrast = %.0f%% of classic', ...
                p.scale, p.scaleErr, abs(p.scale)/max(p.scaleErr, eps), 100*(1 + p.scale));
    bits = {};
    if isfinite(p.tilt)
        bits{end+1} = sprintf('gate-tail tilt %+.4f \\pm %.4f (%.1f\\sigma)', ...
                              p.tilt, p.tiltErr, abs(p.tilt)/max(p.tiltErr, eps));
    end
    if drift.measured
        bits{end+1} = sprintf('drift %+.3f', drift.scale);
    end
    bits{end+1} = sprintf('r %+.3f of %+.3f, shift %+d', m.correlation, m.rMax, m.bestShift);
    T = {a, strjoin(bits, '   |   ')};
end

function axStyle()
    try
        set(gca, 'GridAlpha', 0.12, 'Layer', 'top', 'TickDir', 'out', 'FontSize', 9);
    catch
        set(gca, 'FontSize', 9);          % GridAlpha is R2014b and later
    end
end

function legendFrom(H, L)
% 'northeast' with headroom made for it, NOT 'best'. 'best' is only as good as the
% renderer's collision search, and in the first version of this figure it put the legend
% straight through the walk trace on both left panels -- the series label was struck
% through by the data it labels. Reserving the space and then placing the legend in it is
% deterministic and looks the same on every release.
    try
        lg = legend(H, L, 'Location', 'northeast');
        try; set(lg, 'Box', 'off'); catch; end
    catch
    end
end

function headroom(frac)
% Expand the y axis upward by a fraction of its span, so the legend has empty space to sit
% in instead of covering the data.
    yl = ylim;
    span = max(yl(2) - yl(1), eps);
    ylim([yl(1), yl(2) + frac*span]);
end

function h = errbar(x, y, e, col, lw, ms, name)
% CapSize is R2016b and later. Set separately and in a try, so an older release draws the
% same plot with caps rather than failing and losing the whole figure -- a plot that does
% not appear is a worse outcome than a plot with serifs on the error bars.
    if isempty(e) || numel(e) ~= numel(y); e = zeros(size(y)); end
    try
        h = errorbar(x, y, e, 'o-', 'Color', col, 'LineWidth', lw, 'MarkerSize', ms, ...
                     'MarkerFaceColor', col, 'MarkerEdgeColor', 'w');
    catch
        % Octave's errorbar accepts no property/value pairs after the format string, so
        % the properties go on afterwards. This path exists so the figure still draws when
        % the file is run outside MATLAB, which is where it gets parse-checked.
        h = errorbar(x, y, e, 'o-');
        try
            set(h, 'Color', col, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', col);
        catch
        end
    end
    try; set(h, 'CapSize', 0); catch; end
    try; set(h, 'DisplayName', name); catch; end
end

function r = padRange(x)
    lo = min(x); hi = max(x); w = max(hi - lo, eps);
    r = [lo - 0.02*w, hi + 0.02*w];
end

% ======================================================================

function y = roundTo(x, n)
% round(x, n) is R2014b and absent in Octave, and this file is parse-checked there.
    f = 10^n;
    y = round(x .* f) ./ f;
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
