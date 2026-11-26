function out = esrPulsedVsCwCompare(opts)
%ESRPULSEDVSCWCOMPARE Does pulsed ESR work here, and is it better than CW? Measure both, decide.
%
%   out = esrPulsedVsCwCompare                          % everything from expProperties.json
%   out = esrPulsedVsCwCompare(struct('quick', true))   % fewer settings, about half the time
%   out = esrPulsedVsCwCompare(struct('data', out))     % re-analyse a saved run, no hardware
%
% DIAGNOSTIC. Writes nothing to the json. Saves `out` to the results folder as
% esrPulsedVsCwCompare.mat, prints a table and a verdict, draws two figures.
%
% ======================================================================
% WHAT IT TAKES FROM expProperties.json  (nothing to type)
% ======================================================================
%   frequency           the line every sweep is centred on: the line piTime was calibrated
%                       on, so it needs the calibrated field (field = 'auto', the default)
%   amplitude, piTime   the calibrated pair (-16 dBm, 37 ns). The amplitude is the top of
%                       the pulsed amplitude scan; piTime is used ONLY to place the first
%                       Rabi's tau range. Every pi time used is MEASURED (stage 2).
%   amplitudeDiff, channelPhase, timeDelay, nChannels, MWChannel    via loadJsonParams
%   detectionDuration, laserInitializationDuration, lastDelay       the pulsed readout
%   tauPol, laserInitializationDuration   -> the CW window, via esrWindowPick, exactly as
%                       the run script's "ESR window" cell does it
%
% ======================================================================
% THE PHYSICS: WHAT EACH METHOD TRADES, AND THEREFORE WHAT IS MEASURED
% ======================================================================
% CW. Laser and MW are on together. The laser keeps re-polarising while the MW mixes, so
%   the line is broadened by the MW (power broadening) and by the optical pumping, and
%   contrast and width grow TOGETHER with MW power:
%       FWHM^2 = FWHM0^2 + a*P          C = Cmax * P/(P + Psat)          (P in mW)
%   Narrow line or full contrast -- not both.
% PULSED. Polarise, laser off, ONE pi pulse of length tau, read out. The line is the
%   pulse's own spectral response, so it narrows as 1/tau at FULL contrast:
%       FWHM^2 = (0.80/tau)^2 + FWHMinh^2      (0.80 is exact for a square pi pulse)
%   until tau reaches T2* ~ 1/(pi*FWHMinh): then the width stops falling and the
%   contrast drops. No laser during the MW, so no optical broadening at all.
%
% So one setting against one setting proves nothing. What is compared is each method's
% CONTRAST-vs-LINEWIDTH curve (CW: MW power scan; pulsed: tau scan), and the sensitivity
% each reaches at its best point:
%       eta = 0.77 * FWHM * noise / C      MHz/sqrt(Hz)   (0.77: Lorentzian steepest slope)
% with the noise MEASURED (per-average scatter, so drift counts), per second of sequence
% time (eta_data: the physics) and per second of wall clock (eta_wall: also every re-arm,
% upload and transfer). Lower is better. 1 MHz/sqrt(Hz) = 36 uT/sqrt(Hz).
%
% BOTH METHODS RUN IN WALK MODE (useSequenceWalk = true, newOrderEveryAverage = true): one
% sequence holds the whole sweep, armed once per pg.run, a new random order every average.
% That is how CW ESR is run here, and it is the fast path the pulsed class now has too, so
% the comparison is pulsed-walk against CW-walk, wall clock included. Because pulsed walk is
% new, stage 3 first checks it against the pulsed CLASSIC path (the one already seen
% working) on the same sweep, at piTime and on the narrowest line.
%
% Also compared: line centres (must agree -- same transition), PL level, the off-resonance
% ratio (should be 1), and the dip in the REFERENCE trace (the slow-NV artefact found by
% esrRefDipScan; the pulsed reference sits in the same repeat with the same optical
% history, so it should be ~0 there).
%
% ======================================================================
% STAGES
% ======================================================================
%   0 survey    CW, +-40 MHz: finds the line at the json frequency (STOPS if it is not
%               there -- usually the field is off) and how close its neighbours are, which
%               caps every later sweep so only this one line is fitted
%   1 CW        MW power scan                 -> width, contrast, sensitivity per power
%   2 pulsed    MW amplitude scan. At each amplitude: a Rabi ON the line (two channels,
%               as calibrated), tau range covering ~4 oscillations, fitted by rabiFitRobust
%               -> the MEASURED pi time; then the pulsed spectrum with that pi pulse ->
%               width, contrast, sensitivity. Nothing is scaled: the range of each Rabi is
%               placed from the previous measured one, and re-placed once if the fit shows
%               fewer than 2 or more than 10 oscillations in it.
%   3 walk      pulsed classic vs pulsed walk, same sweep, at piTime and at the narrowest
%               line: correlation, contrast ratio, centre shift, wall-clock speed-up
%   4 knobs     at the best tau, a FULL spectrum (the same sweep as stage 2) for each
%               readout window, dark time before the pulse (singlet emptying, ~0.2-0.3 us)
%               and laser (re)initialisation time (compare tauPol): contrast, and its SNR
%               per sqrt(second)
%   5 pi check  OFF by default (piCheck): -3/0/+3 dB around a measured pi. The stage-2
%               Rabi already measures pi at every amplitude, so this is only a cross-check.
%   Stage 4/5 contrast: each spectrum is fitted twice. Free fit (centre, width, depth) for
%   the table, and a fit with the centre and width FIXED to the stage-2 line (only depth and
%   a linear baseline free) for the comparison between settings: these knobs change the
%   depth, not the line, and a fixed shape gives the same estimator at every setting
%   instead of a free fit that may wander on a weak line.
%
% EVERY MEASUREMENT IS A NEW EXPERIMENT OBJECT (constructor + loadJsonParams + settings), so
% the framework never finds changed parameters on an object that holds data, and never
% opens its "Parameters changed -- restart?" dialog. Side effect: the 'AWG ESR' and 'AWG ESR
% Pulsed' objects in your workspace are replaced (make new ones after). The AWG load cache
% lives in the AWG, not in the object, so re-measuring the same sweep still skips the upload.
%
% Every sweep re-centres and sizes itself: span ~3 FWHM, capped by the neighbours, and one
% automatic re-run when the first sweep was too narrow or too coarse for the line it found.
%
% ======================================================================
% THE HYPERFINE: THE SHARPEST TEST OF WHICH METHOD IS BETTER
% ======================================================================
% 14N splits every line into three, 2.16 MHz apart. Whether a method RESOLVES that triplet,
% and at what contrast per line, is a cleaner verdict than any single width or sensitivity
% number: it needs narrow lines AND contrast at the same time, which is exactly the trade
% the two methods make differently. CW can only narrow its lines by lowering the power,
% which costs contrast, and its optical broadening never goes away; pulsed narrows them by
% lengthening the pulse at full contrast, until T2*.
%
% So EVERY CW and stage-2 pulsed spectrum is put through esrFitRobust's multiplet model:
% 0, 1, 2 or 3 dips, shared width and shape, the spacing FITTED (1-5 MHz) and only then
% compared with 2.16, chosen by BIC. Pulsed spectra are also offered the square-pulse
% response (the pulse length is known), which separates the pulse's own width from the
% inhomogeneous one. Per spectrum: the dip count and confidence, dBIC(3 vs 1 dip) -- a
% continuous "how detectable" number, > 0 favours the triplet, > 10 is decisive -- the
% splitting, the per-line FWHM and contrast, and eta_hf, the sensitivity on ONE hyperfine
% line. The summary says at which CW powers and which pulse lengths / amplitudes the triplet
% is resolved, and which method resolves it better.
%
% A sweep must sample 2.16 MHz to resolve it, so every stage 1/2 sweep gets at least one point
% per hfStepMHz (0.5 MHz: ~4 points per spacing), up to hfMaxPoints. Same rule for both
% methods, same sequence time per spectrum (fewer passes per point), so the comparison stays
% fair; set hfStepMHz = 0 for the old fixed nPoints.
%
% ======================================================================
% OPTIONS (struct fields, all optional)
% ======================================================================
%   jsonLocation   '' = setup7Paths().json
%   field          'auto' (default) h.control('AUTO'), the calibrated field, as the run
%                  script's first cell | 'off' (zero field) | 'asIs' (do not touch it)
%   center         MHz. [] = the json frequency ('off': 2870)
%   cwAmplitudes   dBm. [] = [-44 -38 -32 -26 -20]
%   pulseAmplitudes  dBm, stage 2. [] = json amplitude + [0 -6 -12 -18 -24 -28]
%   rabiPeriods    oscillations each stage-2 Rabi covers (default 4)
%   rabiPoints     taus per Rabi (default 100)
%   rabiRepeats, rabiAverages   per Rabi (default 300, 2: the run script's Rabi cell)
%   windows        us, stage 4 readout windows. [] = [0.3 0.6 1 1.5 2.5 4] (+ the json's)
%   darks          us, stage 4 darkBeforeMW.    [] = [0.1 0.3 1 3]
%   inits          us, stage 4 laser init.      [] = [5 10 20] + the json's
%   nPoints        points per spectrum, at least (default 41)
%   hfStepMHz      stage 1/2 sweeps get at least one point per this many MHz, so the 2.16
%                  MHz hyperfine can be resolved (default 0.5; 0 = off)
%   hfMaxPoints    cap on that (default 121)
%   dataSeconds    sequence time per spectrum, the SAME for both methods (default 6)
%   scanSeconds    sequence time per stage 4/5 spectrum (default 4)
%   averages       default 3 (>= 2, so the error bars include drift)
%   validateWalk, knobs, piCheck   stages 3/4/5 on or off (default true, true, false)
%   walkCheckTaus  us, stage 3 pulse lengths. [] = piTime and the narrowest pulsed line
%   maxDigitizerSamples  per pg.run (default 4.8e6 = 32000 segments of 1.5 us at 100 MS/s,
%                  the AWGExpRamsey measurement); sets CW passesPerRun
%   quick          fewer settings and averages
%   plot, save     default true
%   data           a previous out: re-analyse and re-plot only
%
% OUTPUT
%   out.spectra    one entry per spectrum: method, setting, raw per-average traces, and
%                  after analysis f0, fwhm, contrast, PL, refDip, eta_data, eta_wall, ...
%                  and the hyperfine: hfNDips, hfIsotope, hfResolved, hfDbic31, hfA,
%                  hfFwhm (per line), hfC (contrast per line), etaHf, ...
%   out.model.hyperfine   where each method resolves the triplet, and its best spectrum
%   out.knobs      stage 4/5 points
%   out.model      the CW and pulsed physics fits
%   out.best       best settings found for each method
%   out.verdict    the printed conclusions

    DEFAULTS = struct('jsonLocation', '', 'field', 'auto', 'center', [], ...
        'cwAmplitudes', [], 'pulseAmplitudes', [], 'rabiPeriods', 4, 'rabiPoints', 100, ...
        'rabiRepeats', 300, 'rabiAverages', 2, ...
        'windows', [], 'darks', [], 'inits', [], ...
        'nPoints', 41, 'hfStepMHz', 0.5, 'hfMaxPoints', 121, ...
        'dataSeconds', 6, 'scanSeconds', 4, 'averages', 3, ...
        'validateWalk', true, 'walkCheckTaus', [], 'knobs', true, 'piCheck', false, ...
        'maxDigitizerSamples', 4.8e6, 'quick', false, 'plot', true, 'save', true, 'data', []);
    if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
    opts = withDefaults(opts, DEFAULTS);

    if ~isempty(opts.data)
        out = opts.data;
        fprintf('\n=== esrPulsedVsCwCompare (REPLAY of saved data, no hardware) ===\n');
    else
        % The RAW data is saved as it is taken (after every spectrum) and again before any
        % analysis, in its own time-stamped file, so neither a hardware error halfway nor a
        % bug in the analysis can lose a run. Re-analyse with struct('data', out).
        out = acquire(opts);
        checkpoint(out);
    end
    try
        out = analyse(out);
        report(out);
    catch err
        fprintf(2, ['\nThe ANALYSIS failed -- the measurement is safe%s.\n%s\n' ...
                    'Once fixed, re-analyse without measuring again:\n' ...
                    '    load(''%s''); out = esrPulsedVsCwCompare(struct(''data'', out));\n'], ...
                BooleanHelper.ifTrueElse(isfield(out, 'file') && ~isempty(out.file), ...
                    [' in ' getFieldOr(out, 'file', '')], ''), err.getReport(), getFieldOr(out, 'file', '<file>'));
        return
    end
    if isempty(opts.data) && opts.save
        checkpoint(out);
        if isfield(out, 'file') && ~isempty(out.file)
            fprintf('saved to %s   (re-analyse: esrPulsedVsCwCompare(struct(''data'', out)))\n\n', out.file);
        end
    end
    if opts.plot
        try; plotAll(out); catch err; fprintf('plot failed: %s\n', err.message); end
    end
end

function checkpoint(out)
% Save `out` to its time-stamped file (set at the start of acquire). Never throws.
    if ~isfield(out, 'file') || isempty(out.file); return; end
    try
        save(out.file, 'out');
    catch err
        fprintf('could not save to %s (%s)\n', out.file, err.message);
    end
end

function v = getFieldOr(s, f, d)
    if isstruct(s) && isfield(s, f); v = s.(f); else; v = d; end
end

% ======================================================================
%  ACQUISITION
% ======================================================================

function out = acquire(opts)
    S = settingsFromJson(opts);
    opts = S.opts;
    setField(opts.field);
    out = struct('opts', rmfieldSafe(opts, 'data'), 'S', S, 'spectra', [], ...
                 'discarded', [], 'knobs', [], 'walk', [], 'started', datestr(now), 'file', '');
    if opts.save
        try
            out.file = fullfile(setup7Paths().folder, ...
                                sprintf('esrPulsedVsCwCompare_%s.mat', datestr(now, 'yyyymmdd_HHMM')));
            fprintf('Data will be saved as it is taken to %s\n', out.file);
        catch err
            fprintf('No results folder (%s): this run will NOT be saved.\n', err.message);
        end
    end
    printPlan(S, opts);

    % No shared experiment objects: runCw / runPulsed build a new one for every measurement.

    % ---- 0. survey ---------------------------------------------------------------------
    fprintf('\n--- [0] survey: CW at %g dBm, %g +- 40 MHz ---\n', S.surveyAmp, S.center);
    sv = runCw(S, opts, S.center + linspace(-40, 40, 81), S.surveyAmp, 'survey', opts.dataSeconds/2);
    out.survey = sv; checkpoint(out);
    [S, why] = locateLine(S, sv);
    if ~isempty(why)
        out.S = S; out.abort = why;
        fprintf(2, '\nSTOPPED: %s\n', why);
        return
    end
    fprintf(['line at %.2f MHz, FWHM %.1f MHz at %g dBm; sweeps stay within +-%.1f MHz ' ...
             '(neighbouring lines)\n'], S.f0, S.fwhmSurvey, S.surveyAmp, S.halfFree);

    % ---- 1. CW power scan ----------------------------------------------------------------
    for i = 1:numel(opts.cwAmplitudes)
        amp = opts.cwAmplitudes(i);
        fprintf('\n--- [1.%d/%d] CW %g dBm ---\n', i, numel(opts.cwAmplitudes), amp);
        h0 = clamp(3*S.fwhmSurvey, 4, S.halfFree);
        meas = @(h) runCw(S, opts, S.f0 + linspace(-h, h, sweepPoints(h, opts)), amp, ...
                          sprintf('CW %g dBm', amp), opts.dataSeconds);
        [sp, out] = measureAdaptive(meas, h0, S, out);
        sp.setting = amp;
        out.spectra = appendS(out.spectra, sp);
        checkpoint(out);
    end
    out.S = S;

    % ---- 2. pulsed amplitude scan: Rabi -> measured pi -> spectrum ------------------------
    gInh = guessInhomogeneous(out.spectra);
    piPrev = S.piTime; aPrev = S.amp0;      % where to LOOK for the first Rabi, nothing more
    out.rabi = [];
    amps = sort(opts.pulseAmplitudes, 'descend');   % strong first: each range from the last
    for i = 1:numel(amps)
        A = amps(i);
        fprintf('\n--- [2.%d/%d] pulsed at %.1f dBm: Rabi on the line, then the spectrum ---\n', ...
                i, numel(amps), A);
        R = measureRabi(S, opts, A, piPrev * 10^((aPrev - A)/20));
        out.rabi = appendS(out.rabi, R); checkpoint(out);
        if ~R.ok
            fprintf('   no usable Rabi at %.1f dBm (%s), so no pulsed spectrum there.\n', A, R.note);
            continue
        end
        piPrev = R.piUs; aPrev = A;
        tau = R.piUs;
        fprintf('   measured pi = %.1f +- %.1f ns at %.1f dBm -> pulsed spectrum with that pulse\n', ...
                1e3*tau, 1e3*R.piErrUs, A);
        h0 = clamp(3*sqrt((0.8/tau)^2 + gInh^2), 4, S.halfFree);
        meas = @(h) runPulsed(S, opts, S.f0 + linspace(-h, h, sweepPoints(h, opts)), ...
                    struct('tau', tau, 'amplitude', A), ...
                    sprintf('pulsed %.0f dBm pi %.0f ns', A, 1e3*tau), opts.dataSeconds);
        [sp, out] = measureAdaptive(meas, h0, S, out);
        sp.setting = tau;
        out.spectra = appendS(out.spectra, sp);
        checkpoint(out);
    end

    % ---- 3. walk vs classic --------------------------------------------------------------
    if opts.validateWalk
        isP = find(strcmp({out.spectra.method}, 'pulsed') & [out.spectra.ok]);
        taus = [out.spectra(isP).setting];
        pick = [];
        if ~isempty(opts.walkCheckTaus)
            for t = opts.walkCheckTaus
                [~, k] = min(abs(taus - t)); pick(end+1) = k; %#ok<AGROW>
            end
        elseif ~isempty(isP)
            [~, k1] = min(abs(taus - S.piTime));              % the setting seen working classic
            g = inf(size(isP));
            for j = 1:numel(isP)
                q = quickFit(out.spectra(isP(j)), S.f0);
                if q.ok && q.C > 3*q.CErr; g(j) = q.fwhm; end
            end
            [~, k2] = min(g);                                  % the narrowest line
            pick = unique([k1, k2]);
        end
        W = struct('walkIndex', {}, 'classic', {});
        for k = pick
            ref = out.spectra(isP(k));
            fprintf(['\n--- [3] walk check: pulsed tau %.3g us, the SAME sweep on the CLASSIC path ' ...
                     '(walk took %.0f s) ---\n'], ref.setting, ref.tWall);
            cl = runPulsed(S, opts, ref.f, struct('tau', ref.setting, 'amplitude', ref.cfg.amplitude, ...
                           'walk', false), [ref.label ' CLASSIC'], opts.dataSeconds);
            W(end+1) = struct('walkIndex', isP(k), 'classic', cl); %#ok<AGROW>
            out.walk = W; checkpoint(out);
        end
        out.walk = W;
    end

    % ---- 4/5. knob scans at the best tau --------------------------------------------------
    if opts.knobs || opts.piCheck
        [tauB, fB, gB, ampB] = bestPulsed(out.spectra);
        if isnan(tauB)
            fprintf('\nNo usable pulsed line, so stages 4 and 5 are skipped.\n');
        else
            K = struct('label', {}, 'kind', {}, 'value', {}, 'tau', {}, 'refIndex', {}, 'raw', {});
            base = struct('tau', tauB, 'amplitude', ampB, 'window', S.det, 'dark', 1, ...
                          'init', S.init, 'ampShift', 0);
            if opts.knobs
                [fK, iRef] = sweepFor(out, S, opts, tauB, gB, fB);
                fprintf(['\n--- [4] pulsed knobs at %.1f dBm, measured pi %.1f ns: full spectra, %d ' ...
                         'points over %.2f +- %.1f MHz ---\n'], ampB, 1e3*tauB, numel(fK), mean(fK), ...
                        (max(fK) - min(fK))/2);
                lists = {'window', opts.windows; 'dark', opts.darks; 'init', opts.inits};
                for L = 1:size(lists, 1)
                    for v = lists{L, 2}
                        c = base; c.(lists{L, 1}) = v;
                        lab = sprintf('%s %g us', lists{L, 1}, v);
                        fprintf('\n   [4] %s\n', lab);
                        r = runPulsed(S, opts, fK, c, lab, opts.scanSeconds);
                        K(end+1) = struct('label', lab, 'kind', lists{L, 1}, 'value', v, ...
                                          'tau', tauB, 'refIndex', iRef, 'raw', r); %#ok<AGROW>
                        out.knobs = K; checkpoint(out);
                    end
                end
            end
            if opts.piCheck
                taus = unique([S.piTime, max([out.spectra(strcmp({out.spectra.method}, 'pulsed')).setting])]);
                for t = taus
                    [fT, iRef] = sweepFor(out, S, opts, t, gB, fB);
                    fprintf('\n--- [5] pi check at tau %.3g us: -3 / 0 / +3 dB, full spectra ---\n', t);
                    for d = [-3 0 3]
                        c = base; c.tau = t; c.ampShift = d;
                        c.amplitude = out.spectra(iRef).cfg.amplitude;
                        lab = sprintf('pi check %.3g us %+d dB', t, d);
                        r = runPulsed(S, opts, fT, c, lab, opts.scanSeconds);
                        K(end+1) = struct('label', lab, 'kind', 'pi', 'value', d, 'tau', t, ...
                                          'refIndex', iRef, 'raw', r); %#ok<AGROW>
                    end
                end
            end
            out.knobs = K;
        end
    end
    out.finished = datestr(now);
end

function S = settingsFromJson(opts)
    if isempty(opts.jsonLocation)
        opts.jsonLocation = setup7Paths().json;
    end
    J = jsondecode(fileread(opts.jsonLocation));
    need = {'frequency', 'piTime', 'amplitude', 'detectionDuration', 'laserInitializationDuration'};
    for i = 1:numel(need)
        if ~isfield(J, need{i}) || isempty(J.(need{i}))
            error('esrPulsedVsCwCompare: expProperties.json has no %s (%s).', need{i}, opts.jsonLocation);
        end
    end
    S = struct();
    S.json = opts.jsonLocation;
    S.piTime = J.piTime;
    S.amp0 = J.amplitude(1);
    S.det = J.detectionDuration;
    S.init = J.laserInitializationDuration;
    S.lastDelay = 0.05; if isfield(J, 'lastDelay'); S.lastDelay = J.lastDelay; end
    S.tauPol = NaN; if isfield(J, 'tauPol'); S.tauPol = J.tauPol; end
    S.jsonFrequency = J.frequency(1);
    if ~isempty(opts.center)
        S.center = opts.center;
    elseif strcmpi(opts.field, 'off')
        S.center = 2870;
    else
        S.center = S.jsonFrequency;
    end

    % CW window: the run script's "ESR window" cell, verbatim in effect.
    cw = esrWindowPick(struct('jsonLocation', opts.jsonLocation, 'quiet', true));
    S.cwWin  = cw.detectionDuration;
    S.cwInit = min(cw.laserInit, max(10, floor(S.cwWin + 2*cw.gap)));
    S.cwGap  = cw.gap;

    if opts.quick
        if isempty(opts.cwAmplitudes); opts.cwAmplitudes = [-42 -32 -22]; end
        if isempty(opts.pulseAmplitudes); opts.pulseAmplitudes = S.amp0 + [0 -12 -24]; end
        if isempty(opts.windows);      opts.windows = [0.5 1.5 4]; end
        if isempty(opts.darks);        opts.darks = [0.2 1]; end
        if isempty(opts.inits);        opts.inits = [10 S.init]; end
        opts.averages = min(opts.averages, 2);
        opts.dataSeconds = min(opts.dataSeconds, 4);
    end
    if isempty(opts.cwAmplitudes); opts.cwAmplitudes = [-44 -38 -32 -26 -20]; end
    if isempty(opts.pulseAmplitudes); opts.pulseAmplitudes = S.amp0 + [0 -6 -12 -18 -24 -28]; end
    if isempty(opts.windows);      opts.windows = [0.3 0.6 1 1.5 2.5 4 S.det]; end
    if isempty(opts.darks);        opts.darks = [0.1 0.3 1 3]; end
    if isempty(opts.inits);        opts.inits = [5 10 20 S.init]; end
    opts.pulseAmplitudes = sort(unique(roundTo(opts.pulseAmplitudes(:)', 2)), 'descend');
    opts.windows = unique(opts.windows(opts.windows <= S.init - 1));
    opts.darks = unique(opts.darks(opts.darks >= 0.05));
    opts.inits = unique(opts.inits(opts.inits >= max(opts.windows) + 1));
    opts.averages = max(2, round(opts.averages));
    S.surveyAmp = median(opts.cwAmplitudes);
    S.opts = opts;
end

function setField(mode)
    switch lower(mode)
        case 'asis'
            fprintf('Field left as it is (field = ''asIs'').\n');
        case {'auto', 'off'}
            h = [];
            try; h = getObjByName('Helmholtz'); catch; end
            if isempty(h)
                fprintf(['No Helmholtz object in this session: the field is left as it is. ' ...
                         'Run the run script''s first cell (h.control(''auto'')) if it is off.\n']);
                return
            end
            h.control(upper(mode));     % upper case: control() branches case-sensitively
            fprintf('Helmholtz control set to %s.\n', upper(mode));
        otherwise
            error('esrPulsedVsCwCompare: field must be ''auto'', ''off'' or ''asIs''.');
    end
end

function printPlan(S, opts)
    nCw = numel(opts.cwAmplitudes); nP = numel(opts.pulseAmplitudes);
    cwRep = 2*(S.cwWin + S.cwGap);
    pRep = 2*(1 + S.piTime + S.lastDelay + S.init);
    nRunsCw = ceil(passesFor(opts.dataSeconds, opts.averages, opts.nPoints, cwRep) / ...
                   cwPasses(1e9, opts.nPoints, S.cwWin, opts.maxDigitizerSamples));
    tCw = nCw*1.3*(opts.dataSeconds + opts.averages*nRunsCw*0.25 + 8);
    tR  = nP*1.3*(opts.rabiPoints*opts.rabiRepeats*opts.rabiAverages*2*(S.init + 2)*1e-6 + 15);
    tP  = nP*1.3*(opts.dataSeconds + 6) + tR;
    tW  = opts.validateWalk * 2 * opts.nPoints*opts.averages*0.7;
    nK  = opts.knobs*(numel(opts.windows) + numel(opts.darks) + numel(opts.inits)) + opts.piCheck*6;
    tK  = nK*(opts.scanSeconds + 6);       % full spectra now: data + object + (cached) load
    fprintf('\n=== esrPulsedVsCwCompare ===\n');
    fprintf('json: %s\n', S.json);
    fprintf('centre %.3f MHz (%s), piTime %.4f us at %.1f dBm\n', S.center, ...
            BooleanHelper.ifTrueElse(isempty(opts.center), 'json frequency', 'given'), S.piTime, S.amp0);
    fprintf('pulsed readout: window %.2f us, laser %.1f us, lastDelay %.3f us; repeat ~%.1f us\n', ...
            S.det, S.init, S.lastDelay, pRep);
    fprintf('CW (esrWindowPick): window %g us, laserInit %g us, gap %g us; repeat %.0f us\n', ...
            S.cwWin, S.cwInit, S.cwGap, cwRep);
    fprintf(['CW powers %s dBm | pulsed amplitudes %s dBm, each: a Rabi (%d taus over ~%g ' ...
             'oscillations) -> measured pi -> spectrum\n'], mat2str(opts.cwAmplitudes), ...
            mat2str(opts.pulseAmplitudes), opts.rabiPoints, opts.rabiPeriods);
    fprintf('%d points/spectrum, %g s of sequence each, %d averages. Rough time: %.0f min.\n', ...
            opts.nPoints, opts.dataSeconds, opts.averages, (30 + tCw + tP + tW + tK)/60);
    if opts.hfStepMHz > 0
        fprintf(['Stage 1/2 sweeps: at least one point per %g MHz (up to %d points) so the 2.16 MHz ' ...
                 'hyperfine can be resolved; every one is tested for it.\n'], ...
                opts.hfStepMHz, opts.hfMaxPoints);
    end
end

function n = sweepPoints(h, opts)
% Points for a +-h MHz sweep: nPoints, or enough for one per hfStepMHz -- whichever is more,
% capped at hfMaxPoints. Without it a +-15 MHz CW sweep stepped at 0.75 MHz, a 2.16 MHz
% triplet got under three points per spacing, and "not resolved" would have been the
% sampling's verdict rather than the method's.
    n = opts.nPoints;
    if isfield(opts, 'hfStepMHz') && opts.hfStepMHz > 0
        n = max(n, min(opts.hfMaxPoints, ceil(2*h / opts.hfStepMHz) + 1));
    end
end

function [S, why] = locateLine(S, sv)
    why = '';
    S.f0 = NaN; S.fwhmSurvey = NaN; S.halfFree = 40;
    [y, e] = ratioOf(sv);
    fit = esrFitRobust(struct('f', sv.f, 'y', y, 'e', e, 'baselineOrder', 2, ...
                              'plotIt', false, 'verbose', false));
    if ~fit.ok || isempty(fit.lines)
        why = sprintf(['no ESR line within +-40 MHz of %.2f MHz. With field = ''auto'' this ' ...
                       'is the json frequency, so check the field is really on (and that the ' ...
                       'json is from this session''s calibration).'], S.center);
        return
    end
    c = [fit.lines.centre]; w = [fit.lines.fwhm];
    [d, k] = min(abs(c - S.center));
    if d > max(8, 0.5*w(k))
        why = sprintf(['the nearest line is %.1f MHz away from %.2f MHz (at %.2f). The field ' ...
                       'is not where the json says -- is it on, and is the json current?'], d, S.center, c(k));
        return
    end
    S.f0 = c(k); S.fwhmSurvey = w(k);
    others = setdiff(1:numel(c), k);
    if ~isempty(others)
        S.halfFree = clamp(min(abs(c(others) - S.f0) - w(others)), 4, 40);
    end
end

function [sp, out] = measureAdaptive(meas, h0, S, out)
% One sweep, and at most one re-run when the line it found does not suit the sweep: too
% wide for the span (no baseline left) or too narrow for the step (under-sampled).
    sp = meas(h0);
    m = quickFit(sp, S.f0);
    if ~m.ok; return; end
    step = 2*h0/(numel(sp.f) - 1);
    h1 = clamp(3*m.fwhm, 4, S.halfFree);
    tooWide   = m.fwhm > 0.6*h0 && h1 > 1.2*h0;
    tooCoarse = m.fwhm < 4*step && h1 < 0.8*h0;
    if tooWide || tooCoarse
        fprintf('   FWHM %.2f MHz does not suit +-%.1f MHz: re-measuring at +-%.1f MHz.\n', m.fwhm, h0, h1);
        out.discarded = appendS(out.discarded, sp);
        sp = meas(h1);
    end
end

function sp = runCw(S, opts, freqs, amp, label, seconds)
    % A NEW object every time: on a fresh object (no data yet) the framework restarts
    % without asking, so no "Parameters changed" dialog.
    a = AWGExpESR;
    a.loadJsonParams('all', S.json);
    N = numel(freqs); K = opts.averages;
    repUs = 2*(S.cwWin + a.delayBetweenDetections);
    passes = passesFor(seconds, K, N, repUs);
    a.frequency = freqs;
    a.amplitude = amp;
    a.detectionDuration = S.cwWin;
    a.referenceDetectionDuration = S.cwWin;
    a.laserInitializationDuration = S.cwInit;
    a.repeats = passes;
    a.averages = K;
    a.isTracking = 0;
    a.useSequenceWalk = true;
    a.useSequenceJump = false;
    a.useFgen = false;
    a.driveMWSwitches = true;
    a.randomizeFrequencyOrder = true;
    a.newOrderEveryAverage = true;
    a.passesPerRun = cwPasses(passes, N, S.cwWin, opts.maxDigitizerSamples);
    sp = runAndCollect(a, 'cw', label, freqs, repUs, passes, ...
                       struct('amplitude', amp, 'window', S.cwWin));
end

function sp = runPulsed(S, opts, freqs, c, label, seconds)
    % A NEW object every time (see runCw). loadJsonParams also records piTimeAmplitude.
    p = AWGExpESRPulsed;
    p.loadJsonParams('all', S.json);
    c = withDefaults(c, struct('tau', S.piTime, 'amplitude', NaN, 'window', S.det, 'dark', 1, ...
                               'init', S.init, 'ampShift', 0, 'walk', true));
    if ~isfinite(c.amplitude)
        error('esrPulsedVsCwCompare: runPulsed needs the amplitude its pi time was MEASURED at.');
    end
    N = numel(freqs); K = opts.averages;
    repUs = 2*(c.dark + c.tau + S.lastDelay + c.init);
    passes = passesFor(seconds, K, N, repUs);
    if ~c.walk; passes = min(passes, 2000); end       % classic: repeats per point
    p.frequency = freqs;
    p.mwPulseDuration = c.tau;
    p.amplitude = c.amplitude + c.ampShift;        % the Rabi's amplitude: pi was measured here
    p.detectionDuration = c.window;
    p.darkBeforeMW = c.dark;
    p.laserInitializationDuration = c.init;
    p.referenceMode = 'noMW';
    p.repeats = passes;
    p.averages = K;
    p.isTracking = 0;
    p.useSequenceWalk = c.walk;
    p.randomizeFrequencyOrder = true;
    p.newOrderEveryAverage = c.walk;
    p.passesPerRun = [];
    sp = runAndCollect(p, 'pulsed', label, freqs, repUs, passes, ...
                       struct('tau', c.tau, 'amplitude', p.amplitude(1), 'window', c.window, ...
                              'dark', c.dark, 'init', c.init, 'ampShift', c.ampShift, 'walk', c.walk));
end

function R = measureRabi(S, opts, A, piLook)
% The pi time at amplitude A, MEASURED: a Rabi on the line, fitted by rabiFitRobust (which
% scans every frequency the sampling resolves, so it needs no frequency guess). piLook only
% decides WHERE to look: the tau range covers opts.rabiPeriods periods of 2*piLook. If the
% fit then shows fewer than 2 or more than 10 oscillations in the range, the range is
% re-placed from the FITTED frequency and the Rabi is run again (at most twice).
    R = struct('amplitude', A, 'ok', false, 'note', '', 'f', NaN, 'fErr', NaN, 'piUs', NaN, ...
               'piErrUs', NaN, 'contrast', NaN, 'T', NaN, 'nOsc', NaN, 'r2', NaN, 'tau', [], ...
               'S1', [], 'S2', [], 'E1', [], 'E2', [], 'tries', 0, 'seconds', 0, ...
               'scaledPiUs', S.piTime * 10^((S.amp0 - A)/20));
    period = 2 * piLook;
    for attempt = 1:3
        tauMax = opts.rabiPeriods * period;
        step = max(0.001, roundTo(tauMax / opts.rabiPoints, 3));
        tau = unique(roundTo(step : step : tauMax, 3));
        fprintf('   Rabi %d at %.1f dBm on %.2f MHz: tau %.3f-%.3f us, %d points\n', attempt, A, ...
                S.f0, min(tau), max(tau), numel(tau));
        [o, b, secs, err] = runRabiAt(S, opts, A, tau);
        R.tries = attempt; R.tau = tau; R.seconds = R.seconds + secs;
        if ~isempty(b)
            try
                K = b.currIter;
                R.S1 = reshape(b.signal(1, :, 1:K), [], K); R.S2 = reshape(b.signal(2, :, 1:K), [], K);
                R.E1 = reshape(b.sterr(1, :, 1:K), [], K);  R.E2 = reshape(b.sterr(2, :, 1:K), [], K);
            catch
            end
        end
        if ~isempty(err)
            R.note = err; return
        end
        if ~isstruct(o) || ~isfield(o, 'ok') || ~o.ok
            R.note = 'Rabi fit failed';
            if isstruct(o) && isfield(o, 'why'); R.note = ['Rabi fit failed: ' o.why]; end
            if attempt == 1
                period = 2 * period;        % slower than looked for, most likely: look wider
                fprintf('   %s -- trying a range twice as long.\n', R.note);
                continue
            end
            return
        end
        nOsc = max(tau) * o.f;
        R.f = o.f; R.fErr = o.fErr; R.contrast = o.contrast; R.T = o.T; R.nOsc = nOsc;
        R.piUs = o.pulses(2); R.piErrUs = o.pulsesErr(2);
        try; R.r2 = o.r2; catch; end
        if (nOsc < 2 || nOsc > 10) && attempt < 3
            fprintf('   %.1f oscillations in the range (f = %.3f MHz): re-placing it for %g.\n', ...
                    nOsc, o.f, opts.rabiPeriods);
            period = 1 / o.f;
            continue
        end
        if o.fErr / max(o.f, eps) > 0.05
            R.note = sprintf('Rabi frequency only %.0f %% determined', 100*(1 - o.fErr/o.f));
            return
        end
        R.ok = true; R.note = 'ok';
        fprintf('   Rabi: %.3f +- %.3f MHz, pi = %.1f +- %.1f ns, contrast %.2f %%, %.1f oscillations\n', ...
                o.f, o.fErr, 1e3*R.piUs, 1e3*R.piErrUs, 100*o.contrast, nOsc);
        return
    end
end

function [o, b, secs, err] = runRabiAt(S, opts, A, tau)
% One Rabi as the run script's Rabi cell does it, with the drive the pulsed ESR will use:
% two channels from the json, amplitude A (the setter adds amplitudeDiff, exactly as it does
% for the pulsed ESR), channelPhase and the detection window from the json, ON the line.
% A new object, so no "parameters changed" dialog.
    o = []; err = ''; t0 = tic;
    b = AWGExpRabi;
    b.loadJsonParams('all', S.json);
    b.nChannels = 2;                    % BEFORE amplitude (set.nChannels re-enters set.amplitude)
    b.amplitude = A;
    b.frequency = S.f0;
    b.repeats   = opts.rabiRepeats;
    b.averages  = opts.rabiAverages;
    b.isTracking = 0;
    b.balancedSequence = 1;
    b.shouldAutosave = false;
    b.tau = tau;
    b.useSequenceWalk    = true;
    b.passesPerRun       = [];
    b.walkWarmupTriggers = 1;
    b.resyncEveryPass    = true;
    b.newOrderEveryAverage = true;
    try
        b.run;
        o = rabiFitRobust(struct('exp', b, 'plotIt', false, 'verbose', false));
    catch e
        err = ['Rabi run failed: ' e.message];
    end
    secs = toc(t0);
end

function sp = runAndCollect(e, method, label, freqs, repUs, passes, cfg)
    sp = struct('method', method, 'label', label, 'setting', NaN, 'cfg', cfg, 'f', freqs(:), ...
                'S1', [], 'S2', [], 'E1', [], 'E2', [], 'tWall', NaN, 'repUs', repUs, ...
                'passes', passes, 'averages', NaN, 'ok', false, 'err', '');
    t0 = tic;
    try
        e.run;
        if isprop(e, 'timing') && isstruct(e.timing) && isfield(e.timing, 'repeatDuration')
            sp.repUs = e.timing.repeatDuration;          % the repeat as built, not the estimate
        end
        K = e.currIter;
        if K < 1; error('no average completed'); end
        sp.S1 = reshape(e.signal(1, :, 1:K), [], K);
        sp.S2 = reshape(e.signal(2, :, 1:K), [], K);
        sp.E1 = reshape(e.sterr(1, :, 1:K), [], K);
        sp.E2 = reshape(e.sterr(2, :, 1:K), [], K);
        sp.averages = K;
        sp.ok = true;
    catch err
        sp.err = err.message;
        fprintf(2, '   %s FAILED: %s\n', label, err.message);
    end
    sp.tWall = toc(t0);
end

% ======================================================================
%  ANALYSIS
% ======================================================================

function out = analyse(out)
    if isfield(out, 'abort') && ~isempty(out.abort)
        out.verdict = {['Stopped before measuring: ' out.abort]};
        return
    end
    S = out.S;
    if isempty(out.spectra)
        out.verdict = {'Nothing was measured.'}; out.model = []; out.best = [];
        return
    end
    for i = 1:numel(out.spectra)
        m = spectrumMetrics(out.spectra(i), S);
        out.spectra = setFields(out.spectra, i, m);
    end
    sp = out.spectra;
    isCw = strcmp({sp.method}, 'cw') & [sp.ok];
    isP  = strcmp({sp.method}, 'pulsed') & [sp.ok];

    % ---- physics fits ---------------------------------------------------------------------
    M = struct();
    M.cw = cwModel(sp(isCw));
    M.pulsed = pulsedModel(sp(isP));
    M.hyperfine = hyperfineSummary(sp, isCw, isP);
    out.model = M;

    % ---- walk check -----------------------------------------------------------------------
    for j = 1:numel(out.walk)
        cl = out.walk(j).classic;
        m = spectrumMetrics(cl, S);
        f = fieldnames(m);
        for q = 1:numel(f); cl.(f{q}) = m.(f{q}); end
        out.walk(j).classicMetrics = cl;
        out.walk(j).metrics = walkCompare(sp(out.walk(j).walkIndex), cl);
    end

    % ---- knobs ----------------------------------------------------------------------------
    if ~isempty(out.knobs)
        for i = 1:numel(out.knobs)
            ref = out.spectra(out.knobs(i).refIndex);
            out.knobs(i).m = knobMetrics(out.knobs(i).raw, S, ref);
        end
    end

    out.best = bestPoints(sp, isCw, isP);
    out.verdict = verdict(out);
end

function m = spectrumMetrics(sp, S, doHf)
% doHf (default true): also test the spectrum for the nitrogen hyperfine. The stage 4/5 knob
% spectra skip it -- they are judged on contrast at a fixed line shape -- but every spectrum
% that gets it carries the same fields, so the tables can be concatenated.
    if ~exist('doHf', 'var'); doHf = true; end
    m = struct('fitOk', false, 'f0', NaN, 'f0Err', NaN, 'fwhm', NaN, 'fwhmErr', NaN, ...
               'C', NaN, 'CErr', NaN, 'B0', NaN, 'chi2', NaN, 'y', [], 'e', [], 'PL', NaN, ...
               'PLerr', NaN, 'refDip', NaN, 'refDipErr', NaN, 'sigDip', NaN, 'sigDipErr', NaN, ...
               'noiseDensity', NaN, 'etaData', NaN, 'etaWall', NaN, 'tData', NaN);
    m = setFields(m, 1, hfBlankMetrics());
    if ~sp.ok; return; end
    [y, e] = ratioOf(sp);
    m.y = y; m.e = e;
    % PL: the reference trace is the unperturbed level, in the detector's units
    r = mean(sp.S2, 1);
    m.PL = mean(r); m.PLerr = std(r) / sqrt(numel(r));
    q = quickFit(sp, S.f0);
    if ~q.ok; return; end
    m.fitOk = true;
    m.f0 = q.f0; m.f0Err = q.f0Err; m.fwhm = q.fwhm; m.fwhmErr = q.fwhmErr;
    m.B0 = q.B0; m.C = q.C; m.CErr = q.CErr; m.chi2 = q.chi2;

    % Dips in each trace, per average against its own linear wing baseline, projected on
    % the fitted line shape (the same projection esrRefDipScan uses).
    T = 1 ./ (1 + ((sp.f - m.f0) / (m.fwhm/2)).^2);
    wing = abs(sp.f - m.f0) > 1.5*m.fwhm;
    if nnz(wing) < 4
        [~, o] = sort(abs(sp.f - m.f0), 'descend'); wing = false(size(sp.f)); wing(o(1:min(6, end))) = true;
    end
    [m.sigDip, m.sigDipErr] = projectDip(sp.S1, sp.f, wing, T);
    [m.refDip, m.refDipErr] = projectDip(sp.S2, sp.f, wing, T);

    % Sensitivity. tData = sequence time actually spent on ONE frequency point; the measured
    % per-point error times its square root is the noise density in ratio units per sqrt(Hz).
    K = sp.averages;
    tPoint = K * sp.passes * sp.repUs * 1e-6;
    m.tData = tPoint * numel(sp.f);
    m.noiseDensity = median(e) * sqrt(tPoint);
    m.etaData = 0.770 * m.fwhm * m.noiseDensity / m.C;
    m.etaWall = m.etaData * sqrt(sp.tWall / m.tData);

    if doHf
        m = setFields(m, 1, hyperfineMetrics(sp, y, e, m.B0, m.noiseDensity));
    end
end

function h = hfBlankMetrics()
    h = struct('hfTested', false, 'hfWhy', '', 'hfStep', NaN, 'hfNDips', NaN, ...
               'hfIsotope', '', 'hfConfidence', '', 'hfResolved', false, 'hfVerdict', '', ...
               'hfShape', '', 'hfA', NaN, 'hfAErr', NaN, 'hfFwhm', NaN, 'hfFwhmErr', NaN, ...
               'hfInh', NaN, 'hfC', NaN, 'hfCErr', NaN, 'hfP3', NaN, 'hfDbic31', NaN, ...
               'hfDbicVs1', NaN, 'hfRes', NaN, 'etaHf', NaN, 'hfCurveF', [], 'hfCurveY', []);
end

function h = hyperfineMetrics(sp, y, e, B0, noiseDensity)
% Is the nitrogen hyperfine resolved in this spectrum, and how well? esrFitRobust's multiplet
% model decides: 0-3 dips, shared width and shape, spacing FITTED in 1-5 MHz and only then
% compared with 14N (2.16) / 15N (3.03), shapes lorentz and gauss -- plus the square-pulse
% response for a pulsed spectrum, whose pulse length is known -- chosen by BIC.
%
%   hfDbic31   BIC(1 dip) - BIC(3 dips): the continuous "how detectable" number. Positive
%              favours the triplet; above 10 the evidence is decisive. It keeps meaning
%              something where the yes/no verdict flips, which is what places the onset.
%   hfFwhm     the width of ONE hyperfine line (pulse and inhomogeneous together for the
%              pulse shape; hfInh is the inhomogeneous part alone)
%   hfC        contrast of one line, relative to the baseline
%   hfRes      A / hfFwhm: above ~0.6 a Lorentzian triplet shows separate dips
%   etaHf      0.77 * hfFwhm * noise / hfC: the sensitivity on ONE hyperfine line, when the
%              triplet is resolved -- the number that says which method you would sense with
    h = hfBlankMetrics();
    h.hfStep = median(diff(sp.f));
    tau = [];
    if strcmp(sp.method, 'pulsed') && isfield(sp.cfg, 'tau') && isfinite(sp.cfg.tau)
        tau = sp.cfg.tau;
    end
    try
        o = esrFitRobust(struct('f', sp.f, 'y', y, 'e', e, 'nLines', 1, 'baselineOrder', 1, ...
            'tau', tau, 'hyperfine', true, 'resonance', false, 'field', false, ...
            'plotIt', false, 'verbose', false));
    catch err
        h.hfWhy = ['hyperfine fit error: ' err.message]; return
    end
    if ~o.ok; h.hfWhy = o.why; return; end
    hf = o.hyperfine;
    if ~hf.tested || ~hf.ok; h.hfWhy = hf.why; return; end
    h.hfTested = true;
    h.hfNDips = hf.nDips; h.hfIsotope = hf.isotope; h.hfConfidence = hf.confidence;
    h.hfResolved = hf.resolved; h.hfVerdict = hf.verdict; h.hfShape = hf.shape;
    h.hfP3 = hf.prob(4);
    c = hf.candidates;
    nd = [c.nDips]; b = [c.bic];
    if any(nd == 1) && any(nd == 3)
        h.hfDbic31 = b(nd == 1) - b(nd == 3);
    end
    if any(nd == 1) && any(nd == hf.nDips)
        h.hfDbicVs1 = b(nd == 1) - b(nd == hf.nDips);   % the chosen count against one line
    end
    h.hfFwhm = hf.fwhm; h.hfFwhmErr = hf.fwhmErr;
    if strcmp(hf.shape, 'pulse'); h.hfInh = hf.widthParam; end
    if hf.nDips >= 2
        h.hfA = hf.splitting; h.hfAErr = hf.splittingErr;
        h.hfRes = hf.splitting / max(hf.fwhm, eps);
    end
    if isfinite(B0) && B0 > 0
        h.hfC = hf.depth / B0; h.hfCErr = hf.depthErr / B0;
    end
    if h.hfResolved && h.hfC > 0
        h.etaHf = 0.770 * h.hfFwhm * noiseDensity / h.hfC;
    end
    ff = linspace(min(hf.f), max(hf.f), 600)';
    h.hfCurveF = ff; h.hfCurveY = hf.modelAt(ff);
end

function [y, e] = ratioOf(sp)
% Ratio per average first (level drift cancels), then the mean; the error is the LARGER of
% the scatter across averages and the propagated error -- esrPullData's rule.
    R = sp.S1 ./ sp.S2;
    K = size(R, 2);
    y = mean(R, 2);
    prop = sqrt(sum((R .* sqrt((sp.E1 ./ sp.S1).^2 + (sp.E2 ./ sp.S2).^2)).^2, 2)) / K;
    if K >= 2
        e = max(std(R, 0, 2) / sqrt(K), prop);
    else
        e = prop;
    end
    e(~isfinite(e) | e <= 0) = median(e(isfinite(e) & e > 0));
end

function q = quickFit(sp, fGuess)
    q = struct('ok', false, 'fit', [], 'f0', NaN, 'f0Err', NaN, 'fwhm', NaN, 'fwhmErr', NaN, ...
               'B0', NaN, 'C', NaN, 'CErr', NaN, 'chi2', NaN);
    if ~sp.ok; return; end
    [y, e] = ratioOf(sp);
    span = max(sp.f) - min(sp.f);
    try
        % One line, the plain model: no resonance, field or hyperfine stage, so every
        % spectrum is measured the same way (hyperfineMetrics asks that question separately).
        fit = esrFitRobust(struct('f', sp.f, 'y', y, 'e', e, 'nLines', 1, ...
            'seedCentres', fGuess, 'baselineOrder', 1, 'maxFwhmMHz', max(40, 1.5*span), ...
            'minDepthSigma', 2, 'resonance', false, 'field', false, 'hyperfine', false, ...
            'plotIt', false, 'verbose', false, 'name', sp.label));
    catch err
        fprintf('   fit error: %s\n', err.message);
        return
    end
    q.fit = fit;
    if ~fit.ok || isempty(fit.lines); return; end
    [~, k] = min(abs([fit.lines.centre] - fGuess));
    L = fit.lines(k);
    b = fit.baseline(:);
    B0 = sum(b .* (L.centre - fit.xref).^((0:numel(b)-1)'));
    q.ok = true; q.f0 = L.centre; q.f0Err = L.centreErr; q.fwhm = L.fwhm; q.fwhmErr = L.fwhmErr;
    q.B0 = B0; q.C = L.depth / B0; q.CErr = L.depthErr / B0; q.chi2 = fit.chi2Reduced;
end

function [a, ea] = projectDip(Y, f, wing, T)
    K = size(Y, 2); D = zeros(size(Y));
    for k = 1:K
        pp = polyfit(f(wing) - mean(f), Y(wing, k), 1);
        D(:, k) = 1 - Y(:, k) ./ polyval(pp, f - mean(f));
    end
    d = mean(D, 2);
    if K >= 2; s = std(D, 0, 2) / sqrt(K); else; s = std(d(wing)) * ones(size(d)); end
    s(~isfinite(s) | s <= 0) = median(s(isfinite(s) & s > 0));
    w = 1 ./ s.^2;
    a = sum(w .* T .* d) / sum(w .* T.^2);
    ea = 1 / sqrt(sum(w .* T.^2));
    chi = sum(w .* (d - a*T).^2) / max(numel(d) - 1, 1);
    ea = ea * sqrt(max(chi, 1));
end

function M = cwModel(sp)
% FWHM^2 = FWHM0^2 + a*P and C = Cmax*P/(P+Psat), P in mW.
    M = struct('ok', false, 'fwhm0', NaN, 'a', NaN, 'Cmax', NaN, 'PsatDbm', NaN, 'n', 0);
    ok = [sp.fitOk]; sp = sp(ok);
    M.n = numel(sp);
    if numel(sp) < 3; return; end
    P = 10.^([sp.setting] / 10);
    G = [sp.fwhm]; Ge = max([sp.fwhmErr], 0.02*G);
    pf = wlinfit(P, G.^2, 1 ./ (2*G.*Ge).^2);
    M.a = pf(1); M.fwhm0 = sqrt(max(pf(2), 0));
    C = [sp.C]; Ce = max([sp.CErr], 0.02*C);
    best = inf;
    for Ps = logspace(log10(min(P)/30), log10(max(P)*30), 400)
        g = P ./ (P + Ps);
        w = 1 ./ Ce.^2;
        Cm = sum(w .* g .* C) / sum(w .* g.^2);
        chi = sum(w .* (C - Cm*g).^2);
        if chi < best; best = chi; M.Cmax = Cm; M.PsatDbm = 10*log10(Ps); end
    end
    M.ok = true;
end

function M = pulsedModel(sp)
% FWHM^2 = (k/tau)^2 + FWHMinh^2. k = 0.80 for a square pi pulse.
    M = struct('ok', false, 'k', NaN, 'fwhmInh', NaN, 'T2star', NaN, 'Cmax', NaN, ...
               'CdropLong', NaN, 'n', 0);
    ok = [sp.fitOk]; sp = sp(ok);
    M.n = numel(sp);
    if numel(sp) < 2; return; end
    tau = [sp.setting]; G = [sp.fwhm]; Ge = max([sp.fwhmErr], 0.02*G);
    pf = wlinfit(1 ./ tau.^2, G.^2, 1 ./ (2*G.*Ge).^2);
    M.k = sqrt(max(pf(1), 0));
    M.fwhmInh = sqrt(max(pf(2), 0));
    M.T2star = 1 / (pi * max(M.fwhmInh, eps));         % us, Lorentzian estimate
    C = [sp.C];
    M.Cmax = max(C);
    [~, iL] = max(tau);
    M.CdropLong = C(iL) / M.Cmax;
    M.ok = true;
end

function H = hyperfineSummary(sp, isCw, isP)
% Where each method resolves the triplet, its best resolved spectrum, and the splitting.
% Replaces a single 3-line fit at the longest pulse, seeded AT 2.16 MHz and judged by a
% chi2 ratio: that asked one spectrum one question, and it could not say where the triplet
% appears or whether CW resolves it at all.
    H = struct('cw', [], 'pulsed', [], 'A', NaN, 'AErr', NaN, 'isotope', '', ...
               'inconsistent', [], 'nCoarse', 0);
    tested = (isCw | isP) & [sp.hfTested];
    res0 = tested & [sp.hfResolved];

    % ONE isotope per sample. Each spectrum's call is made on its own, and a coarse or noisy
    % one can land on the other nitrogen: on a saved run a 1.35 MHz-stepped pulsed spectrum
    % came out "15N, A = 3.58" beside four clean 2.15-2.19 MHz triplets. So the sample's
    % isotope is the one with the STRONGEST evidence against a single line, and a spectrum
    % that disagrees is reported as inconsistent rather than counted as resolved.
    okRes = res0;
    if any(res0)
        iso = {sp.hfIsotope};
        u = unique(iso(res0));
        score = -Inf(size(u));
        for i = 1:numel(u)
            score(i) = max([sp(res0 & strcmp(iso, u{i})).hfDbicVs1]);
        end
        [~, k] = max(score);
        H.isotope = u{k};
        okRes = res0 & strcmp(iso, H.isotope);
        H.inconsistent = find(res0 & ~okRes);
        r = find(okRes);
        a = [sp(r).hfA]; ae = [sp(r).hfAErr];
        ok = isfinite(a) & isfinite(ae) & ae > 0;
        if any(ok)
            w = 1 ./ ae(ok).^2;
            H.A = sum(w .* a(ok)) / sum(w); H.AErr = 1 / sqrt(sum(w));
        end
    end
    % DETECTED but not cleanly resolved: a multiplet beats one line decisively (dBIC > 10),
    % but the count or the spacing is not settled. That is the onset region, and the answer
    % to "when does it become detectable" lives there as much as in the clean calls.
    detected = tested & ~okRes & [sp.hfDbic31] > 10;
    H.cw = hfMethodSummary(sp, find(isCw), 'cw', okRes, detected);
    H.pulsed = hfMethodSummary(sp, find(isP), 'pulsed', okRes, detected);
    H.nCoarse = nnz(tested & [sp.hfStep] > 0.7);
end

function Q = hfMethodSummary(sp, idx, kind, okRes, detected)
    Q = struct('kind', kind, 'nTested', 0, 'resolved', [], 'detected', [], ...
               'notResolved', [], 'best', [], 'narrowest', [], 'onset', NaN, 'onsetAmp', NaN);
    if isempty(idx); return; end
    t = idx([sp(idx).hfTested]);
    Q.nTested = numel(t);
    if isempty(t); return; end
    [~, k] = min([sp(t).hfFwhm]);
    Q.narrowest = t(k);
    r = t(okRes(t));
    Q.resolved = r;
    Q.detected = t(detected(t));
    Q.notResolved = setdiff(t, [r, Q.detected]);
    if isempty(r); return; end
    eta = [sp(r).etaHf]; eta(~isfinite(eta)) = Inf;
    [~, k] = min(eta);
    Q.best = r(k);
    % The onset: for CW the STRONGEST power that still resolves it (more power, wider lines);
    % for pulsed the SHORTEST pulse that does (shorter pulse, wider line), and its amplitude.
    if strcmp(kind, 'cw')
        Q.onset = max([sp(r).setting]);
    else
        [Q.onset, k] = min([sp(r).setting]);
        Q.onsetAmp = sp(r(k)).cfg.amplitude;
    end
end

function s = settingText(sp)
    if strcmp(sp.method, 'cw')
        s = sprintf('%g dBm', sp.setting);
    else
        s = sprintf('%.0f ns (%.1f dBm)', 1e3*sp.setting, sp.cfg.amplitude);
    end
end

function s = settingList(sp, idx)
    if isempty(idx); s = 'none'; return; end
    parts = cell(1, numel(idx));
    for i = 1:numel(idx); parts{i} = settingText(sp(idx(i))); end
    s = strjoin(parts, ', ');
end

function v = hyperfineVerdict(sp, H, Mp)
    v = {'3. THE 14N HYPERFINE: WHICH METHOD RESOLVES IT? (every stage 1/2 spectrum tested)'};
    names = struct('cw', 'CW', 'pulsed', 'pulsed');
    for kind = {'cw', 'pulsed'}
        Q = H.(kind{1}); nm = names.(kind{1});
        if Q.nTested == 0
            v{end+1} = sprintf('   %s: no spectrum could be tested.', nm); %#ok<AGROW>
            continue
        end
        seen = '';
        if ~isempty(Q.detected)
            seen = sprintf('; a multiplet is seen (dBIC 3-vs-1 > 10) but not cleanly resolved at %s', ...
                           settingList(sp, Q.detected));
        end
        if isempty(Q.resolved)
            s = sp(Q.narrowest);
            v{end+1} = sprintf(['   %s: NOT resolved in any of %d spectra%s. Narrowest: %.2f MHz per ' ...
                'line at %s, dBIC 3-vs-1 %+.1f.'], nm, Q.nTested, seen, s.hfFwhm, ...
                settingText(s), s.hfDbic31); %#ok<AGROW>
            continue
        end
        v{end+1} = sprintf('   %s: resolved at %s%s; not at %s.', nm, settingList(sp, Q.resolved), ...
                           seen, settingList(sp, Q.notResolved)); %#ok<AGROW>
        if strcmp(kind{1}, 'cw')
            v{end+1} = sprintf('      -> CW resolves it up to %g dBm.', Q.onset); %#ok<AGROW>
        else
            v{end+1} = sprintf(['      -> pulsed resolves it from %.0f ns pulses (%.1f dBm) and ' ...
                                'longer.'], 1e3*Q.onset, Q.onsetAmp); %#ok<AGROW>
        end
        b = sp(Q.best);
        inh = '';
        if isfinite(b.hfInh); inh = sprintf(' (inhomogeneous %.2f)', b.hfInh); end
        v{end+1} = sprintf(['      best: %s -- %.2f +- %.2f MHz per line%s, %.2f%% contrast per ' ...
            'line, A = %.3f +- %.3f MHz, dBIC %+.0f, eta_hf %.3g MHz/rtHz.'], settingText(b), ...
            b.hfFwhm, b.hfFwhmErr, inh, 100*b.hfC, b.hfA, b.hfAErr, b.hfDbic31, b.etaHf); %#ok<AGROW>
    end
    bc = H.cw.best; bp = H.pulsed.best;
    if ~isempty(bc) && ~isempty(bp)
        c = sp(bc); p = sp(bp);
        r = p.etaHf / c.etaHf;
        v{end+1} = sprintf(['   => on one hyperfine line pulsed is %s by %.1fx (eta_hf %.3g vs %.3g ' ...
            'MHz/rtHz), with %.1fx the CW contrast per line at %.2fx its width.'], ...
            okWord(r < 1, 'BETTER', 'worse'), max(r, 1/r), p.etaHf, c.etaHf, ...
            p.hfC / c.hfC, p.hfFwhm / c.hfFwhm);
    elseif ~isempty(bp)
        v{end+1} = ['   => only PULSED resolves the triplet here: the clearest single sign that it ' ...
                    'beats CW on this sample.'];
    elseif ~isempty(bc)
        v{end+1} = ['   => only CW resolves it here: the pulsed lines stay too wide or too shallow ' ...
                    'at these pulse lengths -- longer pulses (weaker amplitude) are the next thing ' ...
                    'to try.'];
    elseif H.cw.nTested + H.pulsed.nTested > 0
        v{end+1} = '   => neither resolves it at these settings.';
    end
    if isfinite(H.A)
        v{end+1} = sprintf(['   splitting over the resolved spectra: A = %.3f +- %.3f MHz (%s; 14N ' ...
                            '2.16, 15N 3.03).'], H.A, H.AErr, H.isotope);
    end
    for k = H.inconsistent(:)'
        s = sp(k);
        v{end+1} = sprintf(['   %s %s came out %s (A = %.2f MHz, step %.2f MHz) on a %s sample: not ' ...
                            'counted -- too coarse or too noisy to call.'], upper(s.method), ...
                           settingText(s), s.hfIsotope, s.hfA, s.hfStep, H.isotope); %#ok<AGROW>
    end
    % The pulsed model fits ONE-Lorentzian widths, so once the triplet is there its
    % "inhomogeneous" floor is the multiplet's envelope (~2A + w), not a line width.
    if ~isempty(Mp) && isfield(Mp, 'ok') && Mp.ok && ~isempty(H.pulsed.resolved)
        b = sp(H.pulsed.best);
        inh = '';
        if isfinite(b.hfInh); inh = sprintf(', inhomogeneous %.2f MHz from the pulse-shape fit', b.hfInh); end
        v{end+1} = sprintf(['   NOTE the pulsed model''s %.2f MHz "inhomogeneous" width (and T2* from it) ' ...
            'is fitted to one-Lorentzian widths and includes the hyperfine envelope; one line is ' ...
            '%.2f MHz at %s%s.'], Mp.fwhmInh, b.hfFwhm, settingText(b), inh);
    end
    if H.nCoarse > 0
        v{end+1} = sprintf(['   NOTE %d spectra are stepped coarser than 0.7 MHz: "not resolved" there ' ...
                            'is partly the sampling (hfStepMHz).'], H.nCoarse);
    end
end

function W = walkCompare(w, cl)
% Same sweep, two paths. A correlation alone cannot judge this: on a noisy spectrum two
% CORRECT measurements correlate at 0.8, and a scrambled one can still correlate if the
% line is wide. So the test is whether the DIFFERENCE is noise: chi2 of (walk - classic)
% against the combined error bars, after removing a common offset (drift between the two
% runs), plus contrast and centre agreeing within 3 sigma.
    W = struct('r', NaN, 'Cratio', NaN, 'Cratio_err', NaN, 'df0Sigma', NaN, 'chi2diff', NaN, ...
               'ok', false, 'speedup', NaN);
    if isempty(w.y) || isempty(cl.y); return; end
    a = 1 - w.y(:); b = 1 - cl.y(:);
    cc = corrcoef(a, b); W.r = cc(1, 2);
    d = w.y(:) - cl.y(:); sd = sqrt(w.e(:).^2 + cl.e(:).^2);
    d = d - sum(d ./ sd.^2) / sum(1 ./ sd.^2);
    W.chi2diff = sum((d ./ sd).^2) / max(numel(d) - 1, 1);
    W.Cratio = cl.C / w.C;
    W.Cratio_err = W.Cratio * sqrt((cl.CErr/cl.C)^2 + (w.CErr/w.C)^2);
    W.df0Sigma = abs(cl.f0 - w.f0) / sqrt(cl.f0Err^2 + w.f0Err^2);
    % chi2 of n-1 = 30 points has a 1% tail at ~1.7; 2 leaves room for imperfect error bars
    W.ok = W.chi2diff < 2 && abs(W.Cratio - 1) < 3*W.Cratio_err + 0.05 && W.df0Sigma < 3;
    W.speedup = cl.tWall / w.tWall;     % same sweep, same repeats and averages
end

function m = knobMetrics(sp, S, ref)
% A stage 4/5 spectrum. C/CErr: depth with the line shape FIXED to the stage-2 line `ref`
% (only depth and a linear baseline free). Free-fit values kept alongside for the table.
% fom = C / CErr / sqrt(sequence seconds): the SNR per sqrt(second) to maximise.
    m = struct('C', NaN, 'CErr', NaN, 'fom', NaN, 'PL', NaN, 'freeOk', false, ...
               'freeC', NaN, 'freeCErr', NaN, 'freeFwhm', NaN, 'freeF0', NaN, 'etaData', NaN);
    if ~sp.ok; return; end
    m.PL = mean(sp.S2(:));
    tData = sp.averages * sp.passes * numel(sp.f) * sp.repUs * 1e-6;
    if ref.fitOk
        [m.C, m.CErr] = fixedShape(sp, ref.f0, ref.fwhm);
        m.fom = m.C / m.CErr / sqrt(tData);
    end
    f = spectrumMetrics(sp, S, false);
    if f.fitOk
        m.freeOk = true; m.freeC = f.C; m.freeCErr = f.CErr; m.freeFwhm = f.fwhm;
        m.freeF0 = f.f0; m.etaData = f.etaData;
    end
end

function [C, CErr] = fixedShape(sp, f0, g)
% y = b0 + b1*(f - f0) - A*L(f), L = Lorentzian (f0, FWHM g) normalised to 1 at f0.
% Linear least squares; C = A/b0, error scaled by sqrt(chi2/dof) when that is above 1.
    [y, e] = ratioOf(sp);
    x = sp.f(:) - f0;
    L = 1 ./ (1 + (x / (g/2)).^2);
    X = [ones(size(x)), x, -L];
    w = 1 ./ e(:);
    Xw = X .* w; yw = y(:) .* w;
    b = Xw \ yw;
    cov = inv(Xw' * Xw);
    chi = sum((yw - Xw*b).^2) / max(numel(y) - 3, 1);
    cov = cov * max(chi, 1);
    C = b(3) / b(1);
    CErr = abs(C) * sqrt(cov(3,3)/b(3)^2 + cov(1,1)/b(1)^2);
end

function [f, iRef] = sweepFor(out, S, opts, tau, gB, fB)
% The stage-2 sweep of this tau (already sized to its line), and its index as the shape
% reference. Without one (a custom walkCheckTaus), a sweep sized like stage 2 would.
    isP = find(strcmp({out.spectra.method}, 'pulsed') & [out.spectra.ok]);
    iRef = [];
    if ~isempty(isP)
        [d, k] = min(abs([out.spectra(isP).setting] - tau));
        if d < 1e-6; iRef = isP(k); end
    end
    if ~isempty(iRef)
        f = out.spectra(iRef).f(:)';
    else
        h = clamp(3*sqrt((0.8/tau)^2 + gB^2), 4, S.halfFree);
        f = fB + linspace(-h, h, opts.nPoints);
        [~, k] = min(abs([out.spectra(isP).setting] - tau));   % nearest line as the shape
        iRef = isP(k);
    end
end

function B = bestPoints(sp, isCw, isP)
    B = struct('cw', [], 'pulsed', []);
    for kind = {'cw', 'pulsed'}
        if strcmp(kind{1}, 'cw'); sel = find(isCw); else; sel = find(isP); end
        good = sel([sp(sel).fitOk] & [sp(sel).C] > 3*[sp(sel).CErr]);
        if isempty(good); continue; end
        [~, k] = min([sp(good).etaData]);
        [~, kw] = min([sp(good).etaWall]);
        [~, kn] = min([sp(good).fwhm]);
        B.(kind{1}) = struct('etaData', sp(good(k)).etaData, 'atData', sp(good(k)).setting, ...
                             'etaWall', sp(good(kw)).etaWall, 'atWall', sp(good(kw)).setting, ...
                             'fwhmMin', sp(good(kn)).fwhm, 'atFwhmMin', sp(good(kn)).setting, ...
                             'CatFwhmMin', sp(good(kn)).C, 'idxData', good(k));
    end
end

function [tauB, fB, gB, ampB] = bestPulsed(sp)
    tauB = NaN; fB = NaN; gB = NaN; ampB = NaN;
    if isempty(sp); return; end
    isP = find(strcmp({sp.method}, 'pulsed') & [sp.ok]);
    eta = inf(size(isP)); f0 = nan(size(isP)); g = nan(size(isP));
    for j = 1:numel(isP)
        q = quickFit(sp(isP(j)), mean(sp(isP(j)).f));
        if q.ok && q.C > 3*q.CErr
            [y, e] = ratioOf(sp(isP(j)));  %#ok<ASGLU>
            tPoint = sp(isP(j)).averages * sp(isP(j)).passes * sp(isP(j)).repUs * 1e-6;
            eta(j) = q.fwhm * median(e) * sqrt(tPoint) / q.C;
            f0(j) = q.f0; g(j) = q.fwhm;
        end
    end
    if all(isinf(eta)); return; end
    [~, j] = min(eta);
    tauB = sp(isP(j)).setting; fB = f0(j); gB = g(j); ampB = sp(isP(j)).cfg.amplitude;
end

function g = guessInhomogeneous(sp)
    g = 1;
    isCw = find(strcmp({sp.method}, 'cw') & [sp.ok]);
    w = [];
    for j = isCw
        q = quickFit(sp(j), mean(sp(j).f));
        if q.ok; w(end+1) = q.fwhm; end %#ok<AGROW>
    end
    if ~isempty(w); g = max(1, 0.6*min(w)); end
end

% ======================================================================
%  VERDICT
% ======================================================================

function v = verdict(out)
    v = {};
    sp = out.spectra; S = out.S; M = out.model; B = out.best;
    isCw = strcmp({sp.method}, 'cw') & [sp.fitOk];
    isP  = strcmp({sp.method}, 'pulsed') & [sp.fitOk];
    nPall = nnz(strcmp({sp.method}, 'pulsed'));

    % ---------------- 1. does pulsed work? ----------------
    v{end+1} = '1. DOES PULSED ESR WORK?';
    works = true;
    if ~any(isP)
        v{end+1} = '   NO: no pulsed spectrum produced a line. Check the MW switches, the trigger delay and piTime.';
        works = false;
    else
        sig = [sp(isP).C] ./ [sp(isP).CErr];
        v{end+1} = sprintf('   lines found in %d of %d pulsed spectra; contrast %.2f-%.2f%% at %.0f-%.0f sigma.', ...
            nnz(isP), nPall, 100*min([sp(isP).C]), 100*max([sp(isP).C]), min(sig), max(sig));
        if any(isCw)
            iCw = find(isCw); iPp = find(isP);
            [~, iC] = min([sp(iCw).f0Err]); c = sp(iCw(iC));
            [~, iP] = min([sp(iPp).f0Err]); q = sp(iPp(iP));
            d = q.f0 - c.f0; sd = sqrt(q.f0Err^2 + c.f0Err^2);
            agree = abs(d) < max(3*sd, 0.1*min(q.fwhm, c.fwhm));
            works = works && agree;
            v{end+1} = sprintf('   centre: pulsed %.3f +- %.3f vs CW %.3f +- %.3f MHz -> %+.3f MHz (%.1f sigma) %s', ...
                q.f0, q.f0Err, c.f0, c.f0Err, d, abs(d)/sd, okWord(agree, 'same line', 'DISAGREE'));
        end
        if M.pulsed.ok
            kOk = M.pulsed.k > 0.55 && M.pulsed.k < 1.15;
            works = works && kOk;
            v{end+1} = sprintf(['   width vs pulse length: FWHM^2 = (%.2f/tau)^2 + %.2f^2 MHz^2 (theory 0.80/tau) %s; ' ...
                                'inhomogeneous FWHM %.2f MHz -> T2* ~ %.2f us.'], M.pulsed.k, M.pulsed.fwhmInh, ...
                                okWord(kOk, '-> the line IS the pulse''s spectrum', '-> NOT the pi-pulse scaling'), ...
                                M.pulsed.fwhmInh, M.pulsed.T2star);
            if M.pulsed.CdropLong < 0.7
                v{end+1} = sprintf(['   contrast falls to %.0f%% of its best at the longest pulse: pulses past T2* ' ...
                                    '(or a pi amplitude that no longer holds -- see the pi check).'], 100*M.pulsed.CdropLong);
            end
        end
        sg = isP & [sp.sigDip] > 5*[sp.sigDipErr];
        if any(sg)
            rd = [sp(sg).refDip] ./ [sp(sg).sigDip];
            v{end+1} = sprintf('   reference trace dip: %.1f%% of the signal dip (median of %d clear pulsed lines).', ...
                               100*median(rd), nnz(sg));
        end
    end
    for j = 1:numel(out.walk)
        if ~isfield(out.walk(j), 'metrics') || isempty(out.walk(j).metrics); continue; end
        W = out.walk(j).metrics;
        w = sp(out.walk(j).walkIndex);
        v{end+1} = sprintf(['   pulsed WALK vs CLASSIC at tau %.3g us: difference chi2/dof %.2f (r = %.2f), contrast ' ...
            'classic/walk %.2f +- %.2f, centres %.1f sigma apart; walk %.1fx faster -> %s'], w.setting, ...
            W.chi2diff, W.r, W.Cratio, W.Cratio_err, W.df0Sigma, W.speedup, ...
            okWord(W.ok, 'same spectrum, walk is safe', 'WALK DISAGREES: use classic until understood'));
        works = works && W.ok;
    end
    if isfield(out, 'rabi') && ~isempty(out.rabi)
        okR = [out.rabi.ok];
        if any(okR)
            dev = 100 * ([out.rabi(okR).piUs] ./ [out.rabi(okR).scaledPiUs] - 1);
            v{end+1} = sprintf(['   pi measured by Rabi at %d of %d amplitudes; 20*log10 scaling from the ' ...
                                'json pair would have been off by %+.0f to %+.0f %% (%s).'], nnz(okR), ...
                                numel(okR), min(dev), max(dev), ...
                                BooleanHelper.ifTrueElse(max(abs(dev)) < 10, 'MW chain close to linear', ...
                                'MW chain NOT linear over this range -- measuring was needed'));
        end
        if any(~okR)
            v{end+1} = sprintf('   no usable Rabi at %s dBm, so no pulsed spectrum there.', ...
                               mat2str([out.rabi(~okR).amplitude]));
        end
    end
    v{end+1} = sprintf('   => %s', okWord(works && any(isP), 'YES, pulsed ESR works.', 'NOT YET -- see the lines above.'));

    % ---------------- 2. is it better? ----------------
    v{end+1} = '2. IS IT BETTER THAN CW?';
    if ~any(isCw) || ~any(isP)
        v{end+1} = '   cannot say: one of the methods has no usable spectrum.';
    else
        % contrast at equal linewidth
        gC = [sp(isCw).fwhm]; cC = [sp(isCw).C];
        [gC, o] = sort(gC); cC = cC(o);
        gP = [sp(isP).fwhm]; cP = [sp(isP).C];
        inR = gP >= min(gC) & gP <= max(gC);
        if any(inR) && numel(gC) >= 2
            cAt = exp(interp1(log(gC), log(cC), log(gP(inR))));
            ratio = median(cP(inR) ./ cAt);
            v{end+1} = sprintf('   at the SAME linewidth pulsed has %.1fx the CW contrast (%d widths compared, %.1f-%.1f MHz).', ...
                ratio, nnz(inR), min(gP(inR)), max(gP(inR)));
        end
        if M.cw.ok
            okM = gP > 1.05*M.cw.fwhm0;
            if any(okM)
                Pm = (gP(okM).^2 - M.cw.fwhm0^2) / M.cw.a;          % CW power giving that width
                cM = M.cw.Cmax * Pm ./ (Pm + 10^(M.cw.PsatDbm/10));
                v{end+1} = sprintf(['   same comparison from the CW MODEL (any power): pulsed/CW contrast %s at ' ...
                                    'FWHM %s MHz.'], mat2str(roundTo(cP(okM) ./ cM, 2)), mat2str(roundTo(gP(okM), 1)));
            end
            if any(~okM)
                v{end+1} = sprintf(['   %d pulsed line(s) are narrower than CW can ever be (%.2f MHz at zero power: ' ...
                                    'optical + inhomogeneous broadening).'], nnz(~okM), M.cw.fwhm0);
            end
        end
        if min(gP) < min(gC)
            [~, iN] = min(gP); cN = cP(iN);
            v{end+1} = sprintf(['   pulsed reaches %.2f MHz at %.2f%% contrast; the narrowest CW line here is %.2f MHz ' ...
                                '(at %.2f%%)%s.'], min(gP), 100*cN, min(gC), 100*cC(1), ...
                                cwFloorNote(M.cw));
        else
            v{end+1} = sprintf('   narrowest lines: pulsed %.2f MHz, CW %.2f MHz.', min(gP), min(gC));
        end
        if ~isempty(B.cw) && ~isempty(B.pulsed)
            rD = B.pulsed.etaData / B.cw.etaData;
            rW = B.pulsed.etaWall / B.cw.etaWall;
            v{end+1} = sprintf(['   best sensitivity per second of SEQUENCE: pulsed %.3g (tau %.3g us) vs CW %.3g MHz/rtHz ' ...
                                '(%g dBm) -> pulsed %s by %.1fx.'], B.pulsed.etaData, B.pulsed.atData, ...
                                B.cw.etaData, B.cw.atData, okWord(rD < 1, 'better', 'worse'), max(rD, 1/rD));
            v{end+1} = sprintf(['   best sensitivity per second of WALL CLOCK: pulsed %.3g vs CW %.3g MHz/rtHz -> pulsed %s ' ...
                                'by %.1fx (%.3g / %.3g uT/rtHz).'], B.pulsed.etaWall, B.cw.etaWall, ...
                                okWord(rW < 1, 'better', 'worse'), max(rW, 1/rW), ...
                                1e2*B.pulsed.etaWall/2.8025, 1e2*B.cw.etaWall/2.8025);
        end
        ovP = median([sp(isP).tWall] ./ [sp(isP).tData]); ovC = median([sp(isCw).tWall] ./ [sp(isCw).tData]);
        v{end+1} = sprintf(['   walk overhead (wall clock / sequence time, uploads included): pulsed %.1fx, CW %.1fx.'], ...
                           ovP, ovC);
        plP = median([sp(isP).PL]); plC = median([sp(isCw).PL]);
        v{end+1} = sprintf('   PL level (reference trace): pulsed %.4g, CW %.4g (ratio %.2f, detector units per us).', ...
                           plP, plC, plP/plC);
        clear_ = [sp.sigDip] > 5*[sp.sigDipErr];
        selP = isP & clear_; selC = isCw & clear_;
        if any(selP) && any(selC)
            v{end+1} = sprintf(['   reference-trace dip / signal dip: pulsed %.1f%%, CW %.1f%% (clear lines; 0 = a clean ' ...
                                'reference, CW''s is the slow-NV artefact).'], ...
                               100*median([sp(selP).refDip] ./ [sp(selP).sigDip]), ...
                               100*median([sp(selC).refDip] ./ [sp(selC).sigDip]));
        end
    end

    % ---------------- 3. the hyperfine ----------------
    v = [v, hyperfineVerdict(sp, M.hyperfine, M.pulsed)];

    % ---------------- 4. knobs and pi check ----------------
    if ~isempty(out.knobs)
        v{end+1} = '4. PULSED SETTINGS (stage 4/5, at the best tau)';
        K = out.knobs;
        for kind = {'window', 'dark', 'init'}
            sel = find(strcmp({K.kind}, kind{1}));
            if isempty(sel); continue; end
            ms = [K(sel).m]; val = [K(sel).value];
            fom = [ms.fom]; C = [ms.C];
            if all(~isfinite(fom)); continue; end
            [~, j] = max(fom);
            switch kind{1}
                case 'window'
                    v{end+1} = sprintf(['   readout window: best SNR/rt(s) at %g us (json %g us); contrast %.2f%% at ' ...
                        '%g us -> %.2f%% at %g us.'], val(j), S.det, 100*C(1), val(1), 100*C(end), val(end));
                case 'dark'
                    k90 = find(C >= 0.9*max(C), 1);
                    v{end+1} = sprintf(['   dark before the pulse: 90%% of the contrast by %g us (singlet ~0.2-0.3 us); ' ...
                        'best SNR/rt(s) at %g us.'], val(k90), val(j));
                case 'init'
                    tp = fitSaturation(val, C);
                    v{end+1} = sprintf(['   laser init: best SNR/rt(s) at %g us (json %g us); contrast saturates with ' ...
                        'tau ~ %.1f us (json tauPol %.1f us).'], val(j), S.init, tp, S.tauPol);
            end
        end
        selPi = find(strcmp({K.kind}, 'pi'));
        tl = unique([K(selPi).tau]);
        for t = tl(:)'                          % (:)': unique of [] is 0x1, and MATLAB loops ONCE over that
            s = selPi([K(selPi).tau] == t);
            if isempty(s); continue; end
            ms = [K(s).m];
            if max([ms.C]) < 4*median([ms.CErr])
                v{end+1} = sprintf(['   pi check at tau %.3g us: contrast too small to judge (%.2f%% +- %.2f%% at best).'], ...
                    t, 100*max([ms.C]), 100*median([ms.CErr]));
                continue
            end
            [k, ~] = fitRotation([K(s).value], [ms.C], [ms.CErr]);
            v{end+1} = sprintf(['   pi check at tau %.3g us: rotation %.2f x pi at its Rabi-measured amplitude -> %s'], ...
                t, k, rotationAdvice(k));
        end
    end

    % ---------------- 5. recommendation ----------------
    if ~isempty(B.pulsed)
        b = sp(B.pulsed.idxData);
        v{end+1} = '5. TO USE PULSED ESR AT ITS BEST HERE';
        v{end+1} = sprintf(['   p.amplitude = %.1f;  p.mwPulseDuration = %.4f;   %% pi MEASURED by Rabi at this ' ...
                            'amplitude; FWHM %.2f MHz'], b.cfg.amplitude, b.setting, b.fwhm);
        v{end+1} = sprintf(['   p.frequency = %.2f + (%.2f:%.2f:%.2f);   %% ~3 FWHM each side, 5 points per FWHM'], ...
                           b.f0, -3*b.fwhm, b.fwhm/5, 3*b.fwhm);
        v{end+1} = '   p.useSequenceWalk = true;  p.newOrderEveryAverage = true;';
    end
end

function s = cwFloorNote(Mcw)
    if Mcw.ok
        s = sprintf('; CW cannot go below ~%.2f MHz at any power (its zero-power width)', Mcw.fwhm0);
    else
        s = '';
    end
end

function t = fitSaturation(x, C)
% C = Cinf*(1 - exp(-x/t)), grid over t.
    t = NaN; best = inf;
    if numel(x) < 3 || any(~isfinite(C)); return; end
    for tt = logspace(log10(min(x)/10), log10(max(x)*3), 300)
        g = 1 - exp(-x/tt);
        Ci = sum(g .* C) / sum(g.^2);
        r = sum((C - Ci*g).^2);
        if r < best; best = r; t = tt; end
    end
end

function [k, Cm] = fitRotation(dB, C, Ce)
% C = Cm * sin^2(pi/2 * k * 10^(dB/20)): k = rotation in units of pi at 0 dB.
    k = NaN; Cm = NaN; best = inf;
    if numel(dB) < 3 || any(~isfinite(C)); return; end
    w = 1 ./ max(Ce, eps).^2;
    for kk = 0.3:0.005:1.8
        g = sin(pi/2 * kk * 10.^(dB/20)).^2;
        c = sum(w .* g .* C) / sum(w .* g.^2);
        r = sum(w .* (C - c*g).^2);
        if r < best; best = r; k = kk; Cm = c; end
    end
end

function s = rotationAdvice(k)
    if ~isfinite(k)
        s = 'no fit';
    elseif abs(k - 1) <= 0.12
        s = 'still a pi pulse (MW chain linear over this range).';
    else
        s = sprintf(['%s: correct the amplitude by %+.1f dB at this length (or recalibrate piTime ' ...
                     'with a Rabi there).'], BooleanHelper.ifTrueElse(k > 1, 'OVER-rotated', 'UNDER-rotated'), ...
                     -20*log10(k));
    end
end

function s = okWord(tf, a, b)
    if tf; s = a; else; s = b; end
end

% ======================================================================
%  REPORT AND PLOTS
% ======================================================================

function report(out)
    if isfield(out, 'abort') && ~isempty(out.abort)
        fprintf('\n%s\n', out.verdict{1}); return
    end
    sp = out.spectra;
    rows = sp;
    for j = 1:numel(out.walk)
        if isfield(out.walk(j), 'classicMetrics') && ~isempty(out.walk(j).classicMetrics)
            c = out.walk(j).classicMetrics; c.setting = sp(out.walk(j).walkIndex).setting;
            rows = appendS(rows, orderfields(c, rows(1)));
        end
    end
    fprintf('\nAll spectra in WALK mode except the CLASSIC rows (stage 3).\n');
    fprintf('%-26s %18s %14s %14s %10s %9s %11s %11s %7s %7s\n', 'spectrum', 'f0 (MHz)', 'FWHM (MHz)', ...
            'contrast (%)', 'PL', 'refDip %', 'eta_data', 'eta_wall', 'wall s', 'wall/');
    fprintf('%-26s %18s %14s %14s %10s %9s %11s %11s %7s %7s\n', '', '', '', '', '', 'of sig', 'MHz/rtHz', ...
            'MHz/rtHz', '', 'seq');
    for i = 1:numel(rows)
        s = rows(i);
        if ~s.ok
            fprintf('%-26s  FAILED: %s\n', s.label, s.err); continue
        end
        if ~s.fitOk
            fprintf('%-26s  no line fitted (PL %.4g)\n', s.label, s.PL); continue
        end
        fprintf('%-26s %10.3f+-%-6.3f %6.2f+-%-5.2f %6.3f+-%-5.3f %10.4g %9.1f %11.3g %11.3g %7.0f %7.1f\n', ...
                s.label, s.f0, s.f0Err, s.fwhm, s.fwhmErr, 100*s.C, 100*s.CErr, s.PL, ...
                100*s.refDip/max(s.sigDip, eps), s.etaData, s.etaWall, s.tWall, s.tWall/s.tData);
    end
    reportHyperfine(sp, out.model);
    if isfield(out, 'rabi') && ~isempty(out.rabi)
        fprintf('\nStage 2 Rabis (two channels as calibrated, on the line). pi is MEASURED; the\n');
        fprintf('scaled column is what 20*log10 scaling from the json pair would have used instead.\n');
        fprintf('%10s %16s %16s %14s %10s %9s %5s  %s\n', 'dBm', 'tau range (us)', 'fRabi (MHz)', ...
                'pi meas (ns)', 'scaled', 'contr %', 'osc', 'note');
        for i = 1:numel(out.rabi)
            r = out.rabi(i);
            rng = '-';
            if ~isempty(r.tau); rng = sprintf('%.3f-%.3f', min(r.tau), max(r.tau)); end
            fprintf('%10.1f %16s %8.3f+-%-6.3f %7.1f+-%-5.1f %10.1f %9.2f %5.1f  %s\n', r.amplitude, rng, ...
                    r.f, r.fErr, 1e3*r.piUs, 1e3*r.piErrUs, 1e3*r.scaledPiUs, 100*r.contrast, r.nOsc, r.note);
        end
    end
    M = out.model;
    if M.cw.ok
        fprintf('\nCW model:     FWHM^2 = %.2f^2 + %.3g*P(mW) MHz^2;  contrast = %.2f%% * P/(P + Psat), Psat = %.1f dBm\n', ...
                M.cw.fwhm0, M.cw.a, 100*M.cw.Cmax, M.cw.PsatDbm);
    end
    if M.pulsed.ok
        fprintf('pulsed model: FWHM^2 = (%.3f/tau)^2 + %.2f^2 MHz^2  (T2* ~ %.2f us), best contrast %.2f%%\n', ...
                M.pulsed.k, M.pulsed.fwhmInh, M.pulsed.T2star, 100*M.pulsed.Cmax);
    end
    if ~isempty(out.knobs)
        fprintf(['\nStage 4/5: full spectra. Contrast with the line shape fixed to the stage-2 line; ' ...
                 'free fit alongside.\n']);
        fprintf('%-26s %16s %10s   %-28s\n', 'spectrum', 'contrast (%)', 'SNR/rt(s)', 'free fit: contrast (%), FWHM');
        for i = 1:numel(out.knobs)
            m = out.knobs(i).m;
            if m.freeOk
                free = sprintf('%.3f+-%.3f, %.2f MHz', 100*m.freeC, 100*m.freeCErr, m.freeFwhm);
            else
                free = 'no line';
            end
            fprintf('%-26s %7.3f+-%-7.3f %10.1f   %-28s\n', out.knobs(i).label, 100*m.C, 100*m.CErr, m.fom, free);
        end
    end
    fprintf('\nWHAT IT SAYS\n');
    for i = 1:numel(out.verdict); fprintf('%s\n', out.verdict{i}); end
    fprintf('\n');
end

function reportHyperfine(sp, M)
    if ~isfield(sp, 'hfTested') || ~any([sp.hfTested]); return; end
    iso = '';
    if isstruct(M) && isfield(M, 'hyperfine'); iso = M.hyperfine.isotope; end
    fprintf(['\nHyperfine per spectrum (esrFitRobust multiplet model, BIC over 0-3 dips). dBIC 3v1 > 0 ' ...
             'favours the triplet, > 10 is decisive;\nFWHM and contrast are per HYPERFINE line; ' ...
             'eta_hf is the sensitivity on one of them.\n']);
    fprintf('%-26s %6s %-11s %-10s %9s %16s %15s %14s %6s %10s\n', 'spectrum', 'step', 'verdict', ...
            'confidence', 'dBIC 3v1', 'A (MHz)', 'line FWHM', 'per line (%)', 'A/FW', 'eta_hf');
    for i = 1:numel(sp)
        s = sp(i);
        if ~s.ok || ~s.fitOk; continue; end
        if ~s.hfTested
            fprintf('%-26s  not tested: %s\n', s.label, s.hfWhy); continue
        end
        if s.hfResolved && (isempty(iso) || strcmp(s.hfIsotope, iso))
            vd = s.hfIsotope;
        elseif s.hfResolved
            vd = [s.hfIsotope '?!'];          % disagrees with the sample: see the verdict
        elseif s.hfDbic31 > 10
            vd = 'seen';
        elseif s.hfNDips >= 2
            vd = sprintf('%d dips ?', s.hfNDips);
        else
            vd = 'no';
        end
        aTxt = '-';
        if isfinite(s.hfA); aTxt = sprintf('%.3f+-%.3f', s.hfA, s.hfAErr); end
        fprintf('%-26s %6.2f %-11s %-10s %+9.1f %16s %7.2f+-%-5.2f %6.3f+-%-5.3f %6.2f %10.3g\n', ...
                s.label, s.hfStep, vd, s.hfConfidence, s.hfDbic31, aTxt, s.hfFwhm, s.hfFwhmErr, ...
                100*s.hfC, 100*s.hfCErr, s.hfRes, s.etaHf);
    end
end

function plotHyperfine(out, cCw, cP)
% Three panels: is the triplet there (dBIC 3 vs 1 against the line width), what each method
% gets per hyperfine line (contrast against width, the hyperfine analogue of the main
% figure's contrast-vs-linewidth), and the best spectrum of each with its multiplet fit.
    sp = out.spectra;
    if ~isfield(sp, 'hfTested'); return; end
    sel = {find(strcmp({sp.method}, 'cw') & [sp.fitOk] & [sp.hfTested]), ...
           find(strcmp({sp.method}, 'pulsed') & [sp.fitOk] & [sp.hfTested])};
    if isempty(sel{1}) && isempty(sel{2}); return; end
    cols = {cCw, cP}; mk = {'o', 's'}; nm = {'CW', 'pulsed'};
    resAll = find([sp.hfResolved]);         % the summary's list when there is one (isotope-consistent)
    if isfield(out, 'model') && isstruct(out.model) && isfield(out.model, 'hyperfine')
        resAll = [out.model.hyperfine.cw.resolved, out.model.hyperfine.pulsed.resolved];
    end
    figure('Color', 'w', 'Position', [80 80 1300 400], 'Name', 'esrPulsedVsCwCompare hyperfine');

    ax = subplot(1, 3, 1); hold(ax, 'on'); styleAx(ax);
    hs = []; ls = {};
    for m = 1:2
        j = sel{m}; if isempty(j); continue; end
        [~, o] = sort([sp(j).setting]); j = j(o);     % joined in power / pulse-length order
        hs(end+1) = plot(ax, [sp(j).fwhm], [sp(j).hfDbic31], [mk{m} '-'], 'Color', cols{m}); %#ok<AGROW>
        ls{end+1} = nm{m}; %#ok<AGROW>
        r = j(ismember(j, resAll));
        plot(ax, [sp(r).fwhm], [sp(r).hfDbic31], mk{m}, 'Color', cols{m}, 'MarkerFaceColor', cols{m});
        for k = j
            text(ax, sp(k).fwhm, sp(k).hfDbic31, ['  ' shortSetting(sp(k))], 'FontSize', 7, ...
                 'Color', cols{m});
        end
    end
    xl = xlim(ax);
    plot(ax, xl, [0 0], ':', 'Color', [0.4 0.4 0.4]);
    plot(ax, xl, [10 10], '--', 'Color', [0.4 0.4 0.4]);
    xlabel(ax, 'line FWHM, one-Lorentzian fit (MHz)'); ylabel(ax, '\DeltaBIC, 3 dips vs 1');
    title(ax, 'is the triplet there? (filled: resolved; dashed: decisive)', 'FontWeight', 'bold');
    if ~isempty(hs); legend(ax, hs, ls, 'Location', 'northeast', 'Box', 'off'); end

    ax = subplot(1, 3, 2); hold(ax, 'on'); styleAx(ax); set(ax, 'XScale', 'log');
    for m = 1:2
        j = sel{m}; if isempty(j); continue; end
        ok = isfinite([sp(j).hfC]) & isfinite([sp(j).hfFwhm]); j = j(ok);
        if isempty(j); continue; end
        [~, o] = sort([sp(j).hfFwhm]); j = j(o);
        h = errorbar(ax, [sp(j).hfFwhm], 100*[sp(j).hfC], 100*[sp(j).hfCErr], [mk{m} '-']);
        set(h, 'Color', cols{m});                    % open markers; resolved ones filled below
        r = j(ismember(j, resAll));
        plot(ax, [sp(r).hfFwhm], 100*[sp(r).hfC], mk{m}, 'Color', cols{m}, 'MarkerFaceColor', cols{m});
    end
    yl = ylim(ax);
    plot(ax, [2.16 2.16]/0.577, yl, ':', 'Color', [0.4 0.4 0.4]);   % Lorentzian resolution limit
    ylim(ax, yl);
    xlabel(ax, 'one hyperfine line FWHM (MHz)'); ylabel(ax, 'contrast per line (%)');
    title(ax, 'per hyperfine line (higher-left wins; dotted: A/0.58)', 'FontWeight', 'bold');

    ax = subplot(1, 3, 3); hold(ax, 'on'); styleAx(ax);
    H = []; if isfield(out, 'model') && isfield(out.model, 'hyperfine'); H = out.model.hyperfine; end
    hs = []; ls = {}; off = 0; keys = {'cw', 'pulsed'};
    for m = 1:2
        k = [];
        if ~isempty(H)
            Q = H.(keys{m});
            k = Q.best; if isempty(k); k = Q.narrowest; end
        end
        if isempty(k); continue; end
        s = sp(k);
        C = max(s.C, eps);
        yN = off + 1 - (1 - s.y/s.B0) / C;          % each normalised to its own depth
        plot(ax, s.f - s.f0, yN, '.', 'Color', cols{m}, 'MarkerSize', 9);
        if ~isempty(s.hfCurveF)
            hs(end+1) = plot(ax, s.hfCurveF - s.f0, off + 1 - (1 - s.hfCurveY/s.B0) / C, '-', ...
                             'Color', cols{m}, 'LineWidth', 1.6); %#ok<AGROW>
            ls{end+1} = sprintf('%s %s: %s', nm{m}, settingText(s), ...
                                okWord(s.hfResolved, s.hfIsotope, 'not resolved')); %#ok<AGROW>
        end
        off = off - 1.4;
    end
    yl = ylim(ax);
    for a = [-2.16 0 2.16]; plot(ax, [a a], yl, ':', 'Color', [0.6 0.6 0.6]); end
    ylim(ax, yl);
    xlabel(ax, 'f - f_0 (MHz)'); ylabel(ax, 'normalised dip (offset per method)');
    title(ax, 'best of each method, with the multiplet fit', 'FontWeight', 'bold');
    if ~isempty(hs); legend(ax, hs, ls, 'Location', 'southoutside', 'Box', 'off', 'FontSize', 7); end
end

function s = shortSetting(sp)
    if strcmp(sp.method, 'cw'); s = sprintf('%g dBm', sp.setting);
    else; s = sprintf('%.0f ns', 1e3*sp.setting); end
end

function plotAll(out)
    if isfield(out, 'abort') && ~isempty(out.abort); return; end
    sp = out.spectra;
    isCw = find(strcmp({sp.method}, 'cw') & [sp.fitOk]);
    isP  = find(strcmp({sp.method}, 'pulsed') & [sp.fitOk]);
    cCw = [0.85 0.33 0.10]; cP = [0.15 0.42 0.69];
    figure('Color', 'w', 'Position', [60 60 1300 760], 'Name', 'esrPulsedVsCwCompare');

    for pn = 1:2
        ax = subplot(2, 3, pn); hold(ax, 'on'); styleAx(ax);
        if pn == 1; sel = isCw; col = cCw; ttl = 'CW (by MW power)'; else; sel = isP; col = cP; ttl = 'pulsed (by pulse length)'; end
        n = numel(sel); hs = []; ls = {};
        for j = 1:n
            s = sp(sel(j)); sh = 0.35 + 0.65*(j-1)/max(n-1, 1);
            cc = 1 - sh*(1 - col);
            yN = 1 - (1 - s.y/s.B0) / max(s.C, eps);      % each normalised to its own depth
            hs(end+1) = plot(ax, s.f - s.f0, yN, '.-', 'Color', cc); %#ok<AGROW>
            ls{end+1} = s.label; %#ok<AGROW>
        end
        xlabel(ax, 'f - f_0 (MHz)'); ylabel(ax, 'normalised dip');
        title(ax, ttl, 'FontWeight', 'bold');
        if ~isempty(hs); legend(ax, hs, ls, 'Location', 'southeast', 'Box', 'off', 'FontSize', 7); end
    end

    ax = subplot(2, 3, 3); hold(ax, 'on'); styleAx(ax); set(ax, 'XScale', 'log');
    errbar(ax, [sp(isCw).fwhm], 100*[sp(isCw).C], 100*[sp(isCw).CErr], cCw, 'o-');
    errbar(ax, [sp(isP).fwhm], 100*[sp(isP).C], 100*[sp(isP).CErr], cP, 's-');
    xlabel(ax, 'FWHM (MHz)'); ylabel(ax, 'contrast (%)');
    title(ax, 'contrast vs linewidth (higher-left wins)', 'FontWeight', 'bold');
    legend(ax, {'CW', 'pulsed'}, 'Location', 'best', 'Box', 'off');

    ax = subplot(2, 3, 4); hold(ax, 'on'); styleAx(ax);
    if ~isempty(isP)
        tau = [sp(isP).setting];
        errbar(ax, 1./tau, [sp(isP).fwhm], [sp(isP).fwhmErr], cP, 's');
        M = out.model.pulsed;
        if M.ok
            u = linspace(0, max(1./tau)*1.05, 200);
            plot(ax, u, sqrt((M.k*u).^2 + M.fwhmInh^2), '-', 'Color', cP);
            plot(ax, u, 0.8*u, ':', 'Color', [0.4 0.4 0.4]);
        end
    end
    xlabel(ax, '1/\tau (1/\mus)'); ylabel(ax, 'FWHM (MHz)');
    title(ax, 'pulsed width vs 1/\tau (dotted: 0.8/\tau)', 'FontWeight', 'bold');

    ax = subplot(2, 3, 5); hold(ax, 'on'); styleAx(ax); set(ax, 'XScale', 'log', 'YScale', 'log');
    plot(ax, [sp(isCw).fwhm], [sp(isCw).etaData], 'o-', 'Color', cCw, 'MarkerFaceColor', cCw);
    plot(ax, [sp(isP).fwhm], [sp(isP).etaData], 's-', 'Color', cP, 'MarkerFaceColor', cP);
    plot(ax, [sp(isCw).fwhm], [sp(isCw).etaWall], 'o--', 'Color', cCw);
    plot(ax, [sp(isP).fwhm], [sp(isP).etaWall], 's--', 'Color', cP);
    xlabel(ax, 'FWHM (MHz)'); ylabel(ax, '\eta (MHz/\surdHz)');
    title(ax, 'sensitivity (solid: sequence time, dashed: wall clock)', 'FontWeight', 'bold');
    legend(ax, {'CW', 'pulsed'}, 'Location', 'best', 'Box', 'off');

    ax = subplot(2, 3, 6); hold(ax, 'on'); styleAx(ax);
    lab = {}; vals = []; errs = []; cols = [];
    for j = [isCw, isP]
        lab{end+1} = sp(j).label; vals(end+1) = 100*sp(j).refDip/max(sp(j).sigDip, eps); %#ok<AGROW>
        errs(end+1) = 100*sp(j).refDipErr/max(sp(j).sigDip, eps); %#ok<AGROW>
        cols(end+1, :) = BooleanHelper.ifTrueElse(strcmp(sp(j).method, 'cw'), cCw, cP); %#ok<AGROW>
    end
    for j = 1:numel(vals)
        barh(ax, j, vals(j), 0.6, 'FaceColor', cols(j, :), 'EdgeColor', 'none');
        plot(ax, vals(j) + [-1 1]*errs(j), [j j], 'k-');
    end
    set(ax, 'YTick', 1:numel(lab), 'YTickLabel', lab, 'FontSize', 7);
    xlabel(ax, 'reference-trace dip (% of signal dip)');
    title(ax, 'reference dip (0 = clean reference)', 'FontWeight', 'bold');

    try
        plotHyperfine(out, cCw, cP);
    catch err
        fprintf('hyperfine plot failed: %s\n', err.message);
    end

    if isfield(out, 'rabi') && ~isempty(out.rabi) && any([out.rabi.ok])
        r = out.rabi([out.rabi.ok]);
        figure('Color', 'w', 'Position', [120 120 620 400], 'Name', 'esrPulsedVsCwCompare Rabi');
        ax = axes; hold(ax, 'on'); styleAx(ax);
        x = 10.^(([r.amplitude] - out.S.amp0)/20);          % field relative to the json amplitude
        xx = linspace(0, 1.05*max([x 1]), 50);
        plot(ax, xx, xx/(2*out.S.piTime), ':', 'Color', [0.4 0.4 0.4]);
        errbar(ax, x, [r.f], [r.fErr], cP, 's');
        for i = 1:numel(r)
            text(ax, x(i), r(i).f, sprintf('  %.0f dBm: \\pi %.0f ns', r(i).amplitude, 1e3*r(i).piUs), 'FontSize', 8);
        end
        xlabel(ax, sprintf('MW field relative to %.0f dBm  (10^{\\Delta dB/20})', out.S.amp0));
        ylabel(ax, 'f_{Rabi} (MHz)');
        title(ax, 'measured Rabi frequency vs MW amplitude (dotted: linear from the json pi time)', ...
              'FontWeight', 'bold', 'FontSize', 9);
    end
    if ~isempty(out.knobs)
        K = out.knobs;
        figure('Color', 'w', 'Position', [90 90 1300 330], 'Name', 'esrPulsedVsCwCompare knobs');
        kinds = {'window', 'readout window (\mus)'; 'dark', 'dark before pulse (\mus)'; ...
                 'init', 'laser init (\mus)'; 'pi', 'amplitude offset (dB)'};
        for k = 1:4
            ax = subplot(1, 4, k); hold(ax, 'on'); styleAx(ax);
            sel = find(strcmp({K.kind}, kinds{k, 1}));
            if isempty(sel); continue; end
            ms = [K(sel).m];
            if k < 4
                yyaxis(ax, 'left');
                errbar(ax, [K(sel).value], 100*[ms.C], 100*[ms.CErr], cP, 's-'); ylabel(ax, 'contrast (%)');
                yyaxis(ax, 'right');
                plot(ax, [K(sel).value], [ms.fom], 'o--', 'Color', [0.3 0.3 0.3]); ylabel(ax, 'SNR/\surds');
                set(ax, 'XScale', 'log');
            else
                tl = unique([K(sel).tau]);
                for t = tl(:)'
                    s = sel([K(sel).tau] == t); mm = [K(s).m];
                    errbar(ax, [K(s).value], 100*[mm.C], 100*[mm.CErr], [], 's-');
                end
                ylabel(ax, 'contrast (%)');
                legend(ax, arrayfun(@(t) sprintf('\\tau %.3g \\mus', t), unique([K(sel).tau]), 'UniformOutput', false), ...
                       'Location', 'best', 'Box', 'off');
            end
            xlabel(ax, kinds{k, 2});
        end
    end
end

function styleAx(ax)
    set(ax, 'Box', 'off', 'TickDir', 'out', 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 9);
    try; set(ax, 'GridAlpha', 0.12); catch; end
end

function errbar(ax, x, y, e, col, spec)
    if isempty(x); return; end
    h = errorbar(ax, x, y, e, spec);
    if ~isempty(col)
        set(h, 'Color', col);
        try; set(h, 'MarkerFaceColor', col); catch; end
    end
end

% ======================================================================
%  SMALL HELPERS
% ======================================================================

function n = passesFor(seconds, K, N, repUs)
% Passes (= repeats in walk mode) so that K averages over N points take `seconds` of
% sequence time. At least 20, so walkReduce's drop-the-first-pass rule (> 10) applies.
    n = max(20, ceil(seconds / (K * N * repUs * 1e-6)));
end

function P = cwPasses(passes, N, win, maxSamples)
% CW passesPerRun: as many passes per pg.run as the digitizer holds. One segment is the
% detection window at 100 MS/s, so long CW windows fit far fewer passes than pulsed ones.
    P = floor(maxSamples / (2*(N + 1) * ceil(win*100)));
    P = max(1, min(P, passes));
end

function y = roundTo(x, n)
    y = round(x * 10^n) / 10^n;
end

function x = clamp(x, lo, hi)
    x = min(max(x, lo), hi);
end

function pf = wlinfit(x, y, w)
    x = x(:); y = y(:); w = w(:);
    A = [x, ones(size(x))] .* sqrt(w);
    pf = (A \ (y .* sqrt(w)))';
end

function s = appendS(s, x)
    if isempty(s); s = x; else; s(end+1) = x; end
end

function arr = setFields(arr, i, m)
    f = fieldnames(m);
    for k = 1:numel(f)
        arr(i).(f{k}) = m.(f{k});
    end
end

function o = withDefaults(o, d)
    if isempty(o); o = struct(); end
    known = fieldnames(d); given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('esrPulsedVsCwCompare: unknown option "%s". Valid: %s.', bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || isempty(o.(known{i}))
            o.(known{i}) = d.(known{i});
        end
    end
end

function s = rmfieldSafe(s, f)
    if isfield(s, f); s = rmfield(s, f); end
end
