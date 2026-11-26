function out = calibCoilMatrix(opts)
%CALIBCOILMATRIX Calibrate the Helmholtz coils against the NVs: the 3x3 current-to-field matrix.
%
%   cc = calibCoilMatrix(struct('dryRun', true))     % the plan and the time, nothing moves
%   cc = calibCoilMatrix                              % measure (~20 ESR sweeps) and fit
%   cc = calibCoilMatrix(struct('replay', cc))        % re-fit the same spectra, no hardware
%   cc = calibCoilMatrix(struct('replay', 'G:\...\level 3\coilMatrix_20261127_101500.mat'))
%
% LEVEL 3. Replaces the angle scan (calibFieldAlign): once the coils are calibrated against
% the NVs, a field along an NV axis is a CALCULATION -- theta = 54.7356, phi = 45 in the NV
% frame -- not a search.
%
% ===================================================================
% WHAT IT MEASURES, AND WHY THE NVs ARE THE RULER
% ===================================================================
%
% At each of ~20 coil-current settings, spread over the supplies' range, one wide CW ESR.
% Each spectrum's eight lines fix the field there in the NV crystal frame (up to the
% crystal's symmetries), through the exact spin Hamiltonian with gamma = 2.8025 MHz/G -- a
% physical constant, so the field comes out in absolute Gauss with nothing assumed about the
% coils. coilMatrixFit then finds the one linear model
%
%       B = M * I + b0          (B in the NV frame, M 3x3 in G/A, b0 in G)
%
% that predicts all ~150 line positions at once. That answers both open questions about the
% coils at the same time:
%
%   - the current-to-field RATIOS were never measured (the 5.382 / 8.034 / 13.649 G/A in
%     setupInfo.json are nominal), and
%   - a current in one coil makes field along the OTHER crystal axes too: the diamond is
%     not square to the coils and the coils are not exactly square to each other. The
%     off-diagonal terms of M are exactly that.
%
% The NVs, not the coils, decide what x, y and z are: the cubic axes of the crystal. The
% coils are just three vectors in that frame (columns of M), whatever their physical tilt.
%
% ===================================================================
% CHOICES
% ===================================================================
%
% CURRENTS spread over the whole box (coilCurrentDesign), not small steps around one point:
% the fit is linear in I, so leverage is spread. Kept where the model in use predicts
% Brange (default 12-40 G): below that the lines crowd inside two linewidths of D, above it
% the outer pair leaves the sweep. The model in use is the nominal one before the first
% calibration, so the sweep is padded 20% for its scale error.
%
% ONE FIXED FREQUENCY LIST for every point: an unchanged list skips the AWG re-upload (and
% the "parameters have changed" dialog). Wide enough for the outer pair at 1.2*Bmax.
%
% LOW MW AMPLITUDE (-30 dBm). The line positions are all that is used, and power broadening
% is what blends lines -- every line below 2850 came out at 18-23 MHz at higher power here.
%
% RAW DATA SAVED AFTER EVERY POINT to <results>\level 3\coilMatrix_<stamp>.mat, so a crash at
% point 17 of 20 does not cost the first 16, and so the fit can be redone with 'replay'
% after any change to the analysis.
%
% Restores the coil currents (and the control state) at the end and on any error.
%
% ===================================================================
% OPTIONS (struct fields, all optional)
% ===================================================================
%   currents       N x 3 list of currents to measure at, A (default: coilCurrentDesign)
%   nPoints        (default 20)
%   Brange         [min max] |B| for the design, G (default [12 40])
%   frac           [lo hi] fraction of each supply's range to use (default [0.05 1])
%   freqs          MHz (default: D0 +/- (1.2*gamma*Brange(2) + 10), step freqStep)
%   freqStep       MHz (default 1.5)
%   amplitude      dBm per channel (default -30)
%   repeats, averages   (default 100, 3)
%   settle         s after each current change (default 2)
%   jsonLocation   expProperties.json loaded into each ESR (default setup7Paths().json)
%   nominalM, nominalB0, gaugeRef    passed to coilMatrixFit. Defaults from the Helmholtz
%                  object: its nominal model, and as gaugeRef its present calibration if it
%                  has one (so the NV axes keep their names across recalibrations)
%   fitOpts        struct of extra coilMatrixFit options
%   esrOpts        struct of extra esrFitRobust options
%   replay         a previous output, a struct with .spectra (I, freqs, signal, sterr), or
%                  the path of a saved coilMatrix_*.mat: fit again with no hardware
%   dryRun         print the design, the sweep and an estimate, touch nothing (default false)
%   saveRaw        write the raw spectra as they come in (default true)
%   restore        put the currents back at the end (default true)
%   isTracking     (default 0)
%   plotIt, verbose  (default true)
%
% OUTPUT
%   ok, why
%   fit            coilMatrixFit's full output (M, b0, D, errors, per-spectrum table, ...)
%   calibration    fit.calibration: what h.setCoilCalibration takes and the json keeps
%   spectra(j)     I, freqs, signal, sterr, centres, centreErr, nLines, fitOk, why
%   design, freqs, rawFile, nvReach (reach along each NV direction, new model)
%
% THEN (level3.m 3.2b):
%   h.setCoilCalibration(cc.calibration);  h.save(setup7Paths().json);

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('currents', [], 'nPoints', 20, 'Brange', [12 40], ...
        'frac', [0.05 1], 'freqs', [], 'freqStep', 1.5, 'amplitude', -30, 'repeats', 100, ...
        'averages', 3, 'settle', 2, 'jsonLocation', '', 'nominalM', [], 'nominalB0', [], ...
        'gaugeRef', [], 'fitOpts', struct(), 'esrOpts', struct(), 'replay', [], ...
        'dryRun', false, 'saveRaw', true, 'restore', true, 'isTracking', 0, ...
        'plotIt', true, 'verbose', true, 'D0', 2870));

    GAMMA = 2.8025;
    out = struct('ok', false, 'why', '', 'opts', opts);
    if opts.verbose
        fprintf('\n=== calibCoilMatrix: B = M*I + b0 against the NVs ===\n');
    end

    % ---------------------------------------------------------------------
    %   REPLAY: no hardware at all
    % ---------------------------------------------------------------------
    if ~isempty(opts.replay)
        [spectra, saved] = unpackReplay(opts.replay);
        if opts.verbose
            fprintf('REPLAY of %d spectra. No experiment constructed, no current changed.\n', ...
                    numel(spectra));
        end
        % The models the live run used, so a replay names the NV axes the same way: the
        % run's own gaugeRef and nominal model if it recorded them, else the Helmholtz
        % object's, else setup 7's setupInfo.json values.
        h = [];
        try; h = getObjByName('Helmholtz'); catch; end
        for fn = {'nominalM', 'nominalB0', 'gaugeRef'}
            if isempty(opts.(fn{1})) && isfield(saved, fn{1}) && ~isempty(saved.(fn{1}))
                opts.(fn{1}) = saved.(fn{1});
            end
        end
        if isempty(opts.nominalM)
            if ~isempty(h); opts.nominalM = h.nominalMatrix; else; opts.nominalM = diag([5.382 8.034 13.649]); end
        end
        if isempty(opts.nominalB0)
            if ~isempty(h); opts.nominalB0 = h.nominalOffset; else; opts.nominalB0 = [-0.373 0.179 2.121]; end
        end
        spectra = fitAllLines(spectra, opts);
        out = analyse(out, spectra, opts);
        return
    end

    % ---------------------------------------------------------------------
    %   LIVE
    % ---------------------------------------------------------------------
    h = getObjByName('Helmholtz');
    if isempty(h)
        error('calibCoilMatrix: no Helmholtz object. Start one (setup7Session) first.');
    end
    if isempty(opts.jsonLocation)
        try; opts.jsonLocation = setup7Paths().json; catch; opts.jsonLocation = ''; end
    end
    if isempty(opts.nominalM); opts.nominalM = h.nominalMatrix; end
    if isempty(opts.nominalB0); opts.nominalB0 = h.nominalOffset; end
    if isempty(opts.gaugeRef)
        opts.gaugeRef = h.coilMatrix;       % the present calibration, or the nominal model
    end

    % ---- design ----------------------------------------------------------
    if isempty(opts.currents)
        d = coilCurrentDesign(struct('M', h.coilMatrix, 'b0', h.coilOffset, ...
            'Imin', h.minCurrent, 'Imax', h.maxCurrent, 'frac', opts.frac, ...
            'nPoints', opts.nPoints, 'Brange', opts.Brange));
        Ilist = d.I;
    else
        Ilist = reshape(opts.currents, [], 3);
    end
    if any(any(Ilist < h.minCurrent - 1e-9 | Ilist > h.maxCurrent + 1e-9))
        error('calibCoilMatrix: some requested currents are outside the supplies'' range.');
    end
    Bpred = h.fieldFor(Ilist);
    BmaxPred = max(sqrt(sum(Bpred.^2, 2)));
    if isempty(opts.freqs)
        span = 1.2 * GAMMA * max(opts.Brange(2), BmaxPred) + 10;
        opts.freqs = round(opts.D0 - span) : opts.freqStep : round(opts.D0 + span);
    end
    out.design = Ilist; out.freqs = opts.freqs;

    if opts.verbose
        fprintf('model in use for the design: %s\n', h.calibrationSummary());
        fprintf('%d points, sweep %.0f-%.0f MHz in %.1f MHz steps (%d frequencies), %d dBm\n', ...
                size(Ilist, 1), min(opts.freqs), max(opts.freqs), median(diff(opts.freqs)), ...
                numel(opts.freqs), opts.amplitude);
        fprintf('\n  %3s  %21s  %s\n', '#', 'I (A)', '|B| predicted by the model in use (G)');
        for j = 1:size(Ilist, 1)
            fprintf('  %3d  [%5.3f %5.3f %5.3f]  %6.2f\n', j, Ilist(j, :), norm(Bpred(j, :)));
        end
    end
    if opts.dryRun
        fprintf(['\nDRY RUN: nothing measured, no current changed. Each point is one CW ESR of\n' ...
                 '%d frequencies x %d repeats x %d averages; the first real point prints the\n' ...
                 'time it took and an estimate for the rest.\n'], numel(opts.freqs), ...
                 opts.repeats, opts.averages);
        out.ok = true;
        return
    end

    % ---- where the raw data goes ------------------------------------------
    out.rawFile = '';
    if opts.saveRaw
        try
            P = setup7Paths();
            folder = fullfile(P.folder, 'level 3');
            if ~isfolder(folder); mkdir(folder); end
            out.rawFile = fullfile(folder, ['coilMatrix_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']); %#ok<TNOW1,DATST>
        catch
            out.rawFile = fullfile(tempdir, ['coilMatrix_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']); %#ok<TNOW1,DATST>
        end
        if opts.verbose; fprintf('\nraw spectra saved as they come in: %s\n', out.rawFile); end
    end

    I0 = h.I; ctrl0 = h.helmControl;
    if strcmpi(ctrl0, 'OFF')
        h.control('AUTO');
        if opts.verbose; fprintf('Helmholtz control was OFF; set to AUTO for the run.\n'); end
    end
    spectra = struct('I', {}, 'freqs', {}, 'signal', {}, 'sterr', {}, 'centres', {}, ...
                     'centreErr', {}, 'nLines', {}, 'fitOk', {}, 'why', {}, 'Bsolve', {});
    tStart = tic;
    try
        for j = 1:size(Ilist, 1)
            h.setI(Ilist(j, :));
            pause(opts.settle);
            t0 = tic;
            a = makeExp(opts);
            a.run;
            fit = esrFitRobust(esrArgs(struct('exp', a), opts));
            spectra(j) = spectrumFrom(Ilist(j, :), fit);                       %#ok<AGROW>
            if opts.saveRaw && ~isempty(out.rawFile)
                coilRaw = struct('spectra', spectra, 'opts', opts, 'Ilist', Ilist); %#ok<NASGU>
                save(out.rawFile, 'coilRaw');
            end
            if opts.verbose
                dt = toc(t0);
                left = (size(Ilist, 1) - j) * toc(tStart) / j;
                fprintf('  %2d/%d  I = [%5.3f %5.3f %5.3f] A   %d lines, |B| ~ %5.2f G   (%.0f s; ~%.0f min left)\n', ...
                        j, size(Ilist, 1), Ilist(j, :), spectra(j).nLines, ...
                        spectra(j).Bsolve, dt, left/60);
            end
        end
    catch err
        if opts.restore; restoreCoils(h, I0, ctrl0); end
        fprintf(2, ['\ncalibCoilMatrix stopped at point %d: %s\nThe spectra taken so far are ' ...
                    'in %s -- fit them with calibCoilMatrix(struct(''replay'', thatFile)).\n'], ...
                numel(spectra) + 1, err.message, out.rawFile);
        rethrow(err);
    end
    if opts.restore; restoreCoils(h, I0, ctrl0); end

    out = analyse(out, spectra, opts);
end

% ======================================================================

function out = analyse(out, spectra, opts)
    out.opts = opts;            % as finally used, so a replay of this output can reuse them
    out.spectra = spectra;
    use = [spectra.fitOk] & [spectra.nLines] >= 3;
    if opts.verbose
        fprintf('\n%d of %d spectra have 3+ lines and go into the fit.\n', nnz(use), numel(spectra));
    end
    data = struct('I', {}, 'centres', {}, 'centreErr', {}, 'sweep', {});
    for j = find(use)
        s = spectra(j);
        data(end+1) = struct('I', s.I, 'centres', s.centres, 'centreErr', s.centreErr, ...
                             'sweep', [min(s.freqs) max(s.freqs)]);           %#ok<AGROW>
    end
    fo = opts.fitOpts;
    fo.nominalM = opts.nominalM; fo.nominalB0 = opts.nominalB0;
    if ~isempty(opts.gaugeRef); fo.gaugeRef = opts.gaugeRef; end
    fo.verbose = opts.verbose;
    fit = coilMatrixFit(data, fo);
    out.fit = fit;
    out.ok = fit.ok;
    out.why = fit.why;
    if ~fit.ok
        if opts.verbose; fprintf('\nNO CALIBRATION: %s\n', fit.why); end
        return
    end
    out.calibration = fit.calibration;

    % ---- what the new model means for alignment ----------------------------
    Imin = [0 0 0]; Imax = [5.2 5.2 5.2];
    h = [];
    try; h = getObjByName('Helmholtz'); catch; end
    if ~isempty(h); Imin = h.minCurrent; Imax = h.maxCurrent; end
    dirs = Helmholtz.nvDirections();
    reach = zeros(8, 2);
    for k = 1:8
        [lo, hi] = Helmholtz.rayInterval(fit.M, fit.b0, Imin, Imax, dirs(k, :));
        reach(k, :) = [lo hi];
    end
    out.nvReach = reach;
    if opts.verbose
        fprintf('\nWith this calibration, along each NV direction (crystal frame):\n');
        fprintf('  %2s  %12s  %7s  %8s   %s\n', 'k', 'direction', 'theta', 'phi', '|B| reachable (G)');
        for k = 1:8
            [~, t, p] = Helmholtz.cart2sphPhys(dirs(k, :));
            if isnan(reach(k, 1)); rs = 'not reachable';
            else; rs = sprintf('%.1f to %.1f', reach(k, 1), reach(k, 2)); end
            fprintf('  %2d  [%+d %+d %+d]  %7.3f  %8.2f   %s\n', k, round(dirs(k, :)*sqrt(3)), t, p, rs);
        end
        fprintf(['\nTo use it:\n' ...
                 '   h = getObjByName(''Helmholtz'');\n' ...
                 '   h.setCoilCalibration(cc.calibration);      %% the currents do not move\n' ...
                 '   h.save(setup7Paths().json);                  %% into expProperties.json, "Helmholtz"\n' ...
                 '   h.alignToNV(k, Bmag);                        %% k from the table above\n']);
    end
    if opts.plotIt; plotFit(out); end
end

function out = fitAllLines(spectra, opts)
% Replay: refit the lines of every spectrum from its raw data (so a fitter change applies).
% A spectrum given as line centres only (no freqs) is kept as it is.
    out = struct('I', {}, 'freqs', {}, 'signal', {}, 'sterr', {}, 'centres', {}, ...
                 'centreErr', {}, 'nLines', {}, 'fitOk', {}, 'why', {}, 'Bsolve', {});
    for j = 1:numel(spectra)
        s = spectra(j);
        if ~isfield(s, 'freqs') || isempty(s.freqs)
            t = struct('I', s.I, 'freqs', [], 'signal', [], 'sterr', [], ...
                       'centres', s.centres, 'centreErr', [], 'nLines', numel(s.centres), ...
                       'fitOk', true, 'why', '', 'Bsolve', NaN);
            if isfield(s, 'centreErr'); t.centreErr = s.centreErr; end
            if isfield(s, 'sweep'); t.freqs = s.sweep; end
            out(j) = t;
            continue
        end
        e = []; if isfield(s, 'sterr'); e = s.sterr(:); end
        fit = esrFitRobust(esrArgs(struct('f', s.freqs(:), 'y', s.signal(:), 'e', e), opts));
        out(j) = spectrumFrom(s.I, fit);
        if opts.verbose
            fprintf('  %2d  I = [%5.3f %5.3f %5.3f] A   %d lines\n', j, s.I, out(j).nLines);
        end
    end
end

function args = esrArgs(args, opts)
    args.plotIt = false;
    args.verbose = false;
    args.refine = true;
    args.maxLines = 10;
    args.D0 = opts.D0;
    args.name = 'coil calibration';
    fn = fieldnames(opts.esrOpts);
    for i = 1:numel(fn); args.(fn{i}) = opts.esrOpts.(fn{i}); end
end

function s = spectrumFrom(I, fit)
    s = struct('I', I, 'freqs', [], 'signal', [], 'sterr', [], 'centres', [], ...
               'centreErr', [], 'nLines', 0, 'fitOk', false, 'why', '', 'Bsolve', NaN);
    if isfield(fit, 'f'); s.freqs = fit.f(:)'; end
    if isfield(fit, 'y'); s.signal = fit.y(:)'; end
    if isfield(fit, 'e'); s.sterr = fit.e(:)'; end
    s.fitOk = isfield(fit, 'ok') && fit.ok;
    if isfield(fit, 'why'); s.why = fit.why; end
    if s.fitOk && isfield(fit, 'lines') && ~isempty(fit.lines)
        s.centres = [fit.lines.centre];
        s.centreErr = [fit.lines.centreErr];
        s.nLines = numel(s.centres);
    end
    if isfield(fit, 'B') && isfinite(fit.B); s.Bsolve = fit.B; end
end

function [spectra, saved] = unpackReplay(r)
% The spectra, and the options of the run that took them (empty struct if unknown).
    saved = struct();
    if ischar(r) || isstring(r)
        S = load(char(r));
        if isfield(S, 'coilRaw'); r = S.coilRaw; else; error('calibCoilMatrix: no coilRaw in %s', r); end
    end
    if isstruct(r) && isfield(r, 'opts') && isstruct(r.opts)
        saved = r.opts;
    end
    if isstruct(r) && isfield(r, 'spectra')
        spectra = r.spectra;
    elseif isstruct(r) && isfield(r, 'I')
        spectra = r;
    else
        error('calibCoilMatrix: replay must be a previous output, a struct with .spectra, or a coilMatrix_*.mat');
    end
end

function restoreCoils(h, I0, ctrl0)
    try; h.setI(I0); catch; end
    try
        if ~strcmpi(h.helmControl, ctrl0); h.control(ctrl0); end
    catch
    end
end

function a = makeExp(opts)
% A FRESH AWGExpESR every point, every mode flag set explicitly -- the same reasons as in
% calibFieldSurvey: no "parameters have changed" dialog (a new object has no previous state
% to compare with), and nothing inherited from whatever ran last. The AWG skips the
% re-upload anyway, because the frequency list never changes.
    a = AWGExpESR;
    if ~isempty(opts.jsonLocation)
        a.loadJsonParams('all', opts.jsonLocation);
    else
        a.setCwWindow();
    end
    a.frequency  = opts.freqs(:)';
    a.averages   = opts.averages;
    a.repeats    = opts.repeats;
    if ~isempty(opts.amplitude); a.amplitude = opts.amplitude; end
    a.isTracking = opts.isTracking;
end

function plotFit(out)
    fit = out.fit;
    S = fit.spectra;
    n = numel(S);
    fig = figure('Name', 'calibCoilMatrix', 'Color', 'w');
    ax1 = subplot(2, 2, [1 2], 'Parent', fig);
    hold(ax1, 'on');
    for j = 1:n
        plot(ax1, j + 0.25*[-1 1]' * ones(1, 8), [1; 1] * S(j).allLines, '-', ...
             'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        ok = ~S(j).outlier;
        plot(ax1, j*ones(1, nnz(ok)), S(j).centres(ok), 'o', 'Color', [0 0.35 0.75], ...
             'MarkerFaceColor', [0 0.35 0.75], 'MarkerSize', 4);
        if any(~ok)
            plot(ax1, j*ones(1, nnz(~ok)), S(j).centres(~ok), 'x', 'Color', [0.85 0.2 0.1], ...
                 'MarkerSize', 7, 'LineWidth', 1.5);
        end
    end
    hold(ax1, 'off');
    xlim(ax1, [0.3 n + 0.7]);
    xlabel(ax1, 'spectrum'); ylabel(ax1, 'line centre (MHz)');
    title(ax1, sprintf(['measured lines (dots), model (grey ticks), outliers (x): ' ...
                        'rms %.3f MHz'], fit.rmsMHz));
    ax2 = subplot(2, 2, 3, 'Parent', fig);
    hold(ax2, 'on');
    for j = 1:n
        ok = ~S(j).outlier;
        plot(ax2, j*ones(1, nnz(ok)), S(j).resid(ok), 'o', 'Color', [0 0.35 0.75], 'MarkerSize', 4);
    end
    plot(ax2, [0.3 n + 0.7], [0 0], 'k-');
    hold(ax2, 'off');
    xlim(ax2, [0.3 n + 0.7]);
    xlabel(ax2, 'spectrum'); ylabel(ax2, 'measured - model (MHz)');
    ax3 = subplot(2, 2, 4, 'Parent', fig);
    own = [S.BmagSolve]; mdl = [S.Bmag];
    plot(ax3, own, mdl - own, 'o', 'Color', [0 0.35 0.75], 'MarkerSize', 5);
    hold(ax3, 'on'); plot(ax3, xlim(ax3), [0 0], 'k-'); hold(ax3, 'off');
    xlabel(ax3, '|B| from that spectrum alone (G)'); ylabel(ax3, 'model - own (G)');
    title(ax3, 'one point far off = a spectrum the model overrules');
end

function o = withDefaults(o, d)
    if isempty(o); o = struct(); end
    known = fieldnames(d); given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('calibCoilMatrix: unknown option "%s". Valid: %s.', bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || isempty(o.(known{i})); o.(known{i}) = d.(known{i}); end
    end
end
