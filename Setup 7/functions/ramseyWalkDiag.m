function out = ramseyWalkDiag(opts)
%RAMSEYWALKDIAG Find out WHY walk is smoother than classic, by measuring contrast vs passesPerRun.
%
%   out = ramseyWalkDiag                      % ~4 min of real runs
%
% ===================================================================
% THE TWO CLUES, AND WHAT THEY RULE OUT
% ===================================================================
%
% In the walk-vs-classic comparison, walk came back SMOOTHER (less fringe contrast) and with
% error bars about THREE TIMES SMALLER. Less signal and less noise together.
%
% That rules out the failure everybody reaches for first. A scrambled block-to-tau mapping
% would make each slot a MIXTURE of different taus, and mixing raises the scatter across
% passes -- walkReduce computes the error from exactly that scatter, so a shuffled or
% mis-mapped walk has LARGER error bars, not smaller. testRamseyWalk also checks the mapping
% offline and passes.
%
% What produces less signal AND less noise is a MOVING AVERAGE: each tau slot collecting a
% little of its neighbours. Averaging neighbours flattens the fringe and reduces the scatter
% at the same time, which is precisely the signature in the plot.
%
% AND THERE IS AN OBVIOUS MECHANISM FOR IT. resyncEveryPass resyncs once per pg.run, but one
% pg.run holds passesPerRun = 228 passes chained through the goto wrap. If the block alignment
% slips by even a fraction of a block per pass -- a dropped trigger, a trigger arriving while
% the sequencer is still playing, an off-by-one at the wrap -- it accumulates over 228 passes
% and walks each slot into its neighbours. The smear width grows with passesPerRun.
%
% THAT IS A TESTABLE PREDICTION, and it is what this does: measure the fringe contrast at
% several values of passesPerRun and see whether it falls as P rises.
%
%   contrast FLAT in P, and equal to classic   -> walk is sound; the disagreement was
%                                                 something else (check the settings match)
%   contrast FALLS as P rises                  -> alignment drift within a pg.run. The fix is
%                                                 a cap on P, and this measures where the cap
%                                                 belongs.
%   contrast LOW even at P = 1                 -> not drift. The fixed gate or the sequence
%                                                 layout is wrong, and no P will save it.
%
% P = 1 resyncs before EVERY pass, so it is the control: at P = 1 walk cannot drift, and it is
% also no faster than classic. If contrast at P = 1 matches classic, everything downstream of
% the mapping is right.
%
% OPTIONS
%   jsonLocation, tauMax (2.5), detuneMHz (0), amplitude (-16), initUs (16)
%   repeats (300)      passes per point, kept the SAME for every P so contrast is comparable
%   passList ([1 10 50 228])
%   doClassic (true)   run the classic reference too
%   plotIt (true)

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('jsonLocation', '', 'tauMax', 2.5, 'detuneMHz', 0, ...
        'amplitude', -16, 'initUs', 16, 'repeats', 300, 'passList', [1 10 50 228], ...
        'doClassic', true, 'plotIt', true), 'ramseyWalkDiag');
    if isempty(opts.jsonLocation); opts.jsonLocation = setup7Paths().json; end

    step = min(1/(6*(abs(opts.detuneMHz) + 2.16)), 0.075);
    tau  = 0.005 : step : opts.tauMax;

    fprintf('\n============================================================\n');
    fprintf('ramseyWalkDiag  --  contrast vs passesPerRun\n');
    fprintf('============================================================\n');
    fprintf('%d taus, %d passes each, at passesPerRun = %s\n', ...
            numel(tau), opts.repeats, mat2str(opts.passList));
    fprintf(['P = 1 is the CONTROL: it resyncs before every pass, so it cannot drift. If its\n' ...
             'contrast matches classic, everything except the multi-pass alignment is right.\n']);
    nArms = sum(ceil(opts.repeats ./ opts.passList)) + opts.doClassic*numel(tau);
    fprintf('about %d arms in total, so roughly %.0f min plus photons.\n', ...
            nArms, nArms*0.588/60);

    R = struct('label', {}, 'P', {}, 'contrast', {}, 'meanSterr', {}, ...
               'x', {}, 'value', {}, 'sterr', {});

    % ---- classic reference FIRST, so the fringe is known before anything is judged -------
    if opts.doClassic
        g = runRamsey(tau, opts, false, [], 'classic');
        [c, ms] = contrastOf(g);
        R(end+1) = struct('label', 'classic', 'P', NaN, 'contrast', c, 'meanSterr', ms, ...
                          'x', g.x, 'value', g.value, 'sterr', g.sterr);       %#ok<AGROW>
        fprintf('\n  classic          contrast %.3g, mean sterr %.3g\n', c, ms);
    end

    for P = opts.passList(:)'
        g = runRamsey(tau, opts, true, P, sprintf('walk P=%d', P));
        [c, ms] = contrastOf(g);
        R(end+1) = struct('label', sprintf('walk P=%d', P), 'P', P, 'contrast', c, ...
                          'meanSterr', ms, 'x', g.x, 'value', g.value, 'sterr', g.sterr); %#ok<AGROW>
        fprintf('  walk P=%-4d      contrast %.3g, mean sterr %.3g\n', P, c, ms);
    end

    out.results = R;
    verdictWalk(R, opts);
    if opts.plotIt; drawIt(R); end
end

% ======================================================================

function g = runRamsey(tau, opts, walk, P, tag)
% One run, then GRAB THE DATA OUT. Experiment objects are singletons keyed by class name, so
% the next constructor deletes this one -- see expGrab.
    fprintf('\n--- %s ---\n', tag);
    r = AWGExpRamsey;
    r.loadJsonParams('all', opts.jsonLocation);
    r.nChannels = 2;
    r.amplitude = opts.amplitude + [0 0];
    r.frequency = r.frequency + opts.detuneMHz;
    if ~isempty(opts.initUs); r.laserInitializationDuration = opts.initUs; end
    r.repeats          = opts.repeats;
    r.averages         = 1;
    r.balancedSequence = 1;
    r.smallDelay       = 0.5;
    r.isTracking       = 0;
    r.shouldAutosave   = false;
    r.constantTime     = 1;          % ignored in walk mode, needed in classic
    r.tau              = tau;
    if walk
        r.useSequenceWalk    = true;
        r.passesPerRun       = P;
        r.walkWarmupTriggers = 1;
        r.resyncEveryPass    = true;
    else
        r.useSequenceWalk    = false;
    end
    styleWas = r.AWG.progressStyle;
    r.AWG.progressStyle = 'quiet';
    t0 = tic;
    r.run;
    secs = toc(t0);
    r.AWG.progressStyle = styleWas;
    g = expGrab(r, tag);
    g.seconds = secs;
    fprintf('  %.0f s\n', secs);
end

function [c, ms] = contrastOf(g)
% Fringe amplitude, robustly, with the smooth rise removed first.
%
% Peak-to-peak of the raw trace would be dominated by the monotonic rise at small tau, which
% is not the fringe. A 5-point moving median is subtracted to leave the oscillation, and the
% contrast is the robust spread of what remains -- 1.4826*MAD, scaled to a peak-to-peak
% equivalent, so it is not hostage to one outlier.
    c = NaN; ms = NaN;
    if isempty(g.value) || numel(g.value) < 7; return; end
    v = g.value(:)';
    trend = movMedian(v, 5);
    resid = v - trend;
    c = 2 * 1.4826 * median(abs(resid - median(resid)));
    ms = median(g.sterr(:));
end

function m = movMedian(v, w)
    n = numel(v); m = zeros(1, n); h = floor(w/2);
    for i = 1:n
        lo = max(1, i-h); hi = min(n, i+h);
        m(i) = median(v(lo:hi));
    end
end

function verdictWalk(R, opts)
    fprintf('\n============================================================\n');
    fprintf('VERDICT\n');
    fprintf('============================================================\n');

    isC = strcmp({R.label}, 'classic');
    W = R(~isC);
    if isempty(W); fprintf('\nNo walk runs to judge.\n'); return; end

    fprintf('\n  %-14s %10s %12s %10s\n', 'run', 'contrast', 'mean sterr', 'vs classic');
    cRef = NaN;
    if any(isC)
        cRef = R(isC).contrast;
        fprintf('  %-14s %10.4g %12.4g %10s\n', 'classic', cRef, R(isC).meanSterr, '--');
    end
    for i = 1:numel(W)
        rel = '';
        if isfinite(cRef) && cRef > 0; rel = sprintf('%.2f', W(i).contrast/cRef); end
        fprintf('  %-14s %10.4g %12.4g %10s\n', W(i).label, W(i).contrast, W(i).meanSterr, rel);
    end

    cs = [W.contrast]; Ps = [W.P];
    [~, iLo] = min(Ps); [~, iHi] = max(Ps);
    fallsWithP = isfinite(cs(iLo)) && isfinite(cs(iHi)) && cs(iHi) < 0.7*cs(iLo);
    lowAtOne   = isfinite(cRef) && isfinite(cs(iLo)) && cs(iLo) < 0.7*cRef;

    if lowAtOne
        fprintf(['\nNOT DRIFT. Contrast is already down to %.2f of classic at passesPerRun = %d,\n' ...
                 'where walk resyncs before every pass and cannot drift at all. So the problem\n' ...
                 'is upstream of the multi-pass alignment: either the FIXED GATE is costing more\n' ...
                 'than the T1 estimate said, or the concatenated sequence is not playing what we\n' ...
                 'think per trigger. No value of passesPerRun will fix it -- do not use walk\n' ...
                 'until this is understood.\n'], cs(iLo)/cRef, Ps(iLo));
        fprintf(['\nNext test, and it is cheap: run walk with TWO taus only, one at a fringe\n' ...
                 'maximum and one at a minimum taken from the classic trace. Two slots cannot\n' ...
                 'smear into a smooth curve, so if their difference still comes out small the\n' ...
                 'pulses themselves are wrong, not the bookkeeping.\n']);
    elseif fallsWithP
        safe = Ps(iLo);
        for i = 1:numel(W)
            if isfinite(cs(i)) && isfinite(cRef) && cs(i) > 0.9*cRef
                safe = max(safe, Ps(i));
            end
        end
        fprintf(['\nALIGNMENT DRIFT WITHIN A pg.run, confirmed. Contrast falls from %.4g at\n' ...
                 'passesPerRun = %d to %.4g at %d -- a factor of %.2f -- which is what a slow\n' ...
                 'slip in block alignment does: each slot averages a little of its neighbours,\n' ...
                 'flattening the fringe and shrinking the scatter at the same time.\n'], ...
                cs(iLo), Ps(iLo), cs(iHi), Ps(iHi), cs(iLo)/max(cs(iHi), eps));
        fprintf(['\nLARGEST passesPerRun THAT KEEPS FULL CONTRAST: %d. Set it explicitly:\n' ...
                 '    r.passesPerRun = %d;\n' ...
                 'and leave maxDigitizerSegments alone -- the digitizer was never the problem.\n'], ...
                safe, safe);
        fprintf(['\nNote what this costs. Walk pays only when passesPerRun > repeats/nTau, so at\n' ...
                 '%d passes and %d taus the threshold is %.0f. If the safe P is below that, walk\n' ...
                 'cannot beat classic on this tau list and the honest answer is to keep classic\n' ...
                 'and take the 5x from the tau grid instead.\n'], ...
                opts.repeats, numel(W(1).x), opts.repeats/max(numel(W(1).x),1));
    else
        fprintf(['\nCONTRAST IS FLAT IN passesPerRun and matches classic, so walk is measuring\n' ...
                 'the same thing after all. The earlier disagreement was not the walk mechanism:\n' ...
                 'check that the two runs really had identical settings -- laserInitialization,\n' ...
                 'constantTime, repeats, averages and the tau list all have to match, and\n' ...
                 'constantTime changes the classic sequence duration but not the walk one.\n']);
    end
    fprintf('\n============================================================\n\n');
end

function drawIt(R)
    try
        figure('Color', 'w');
        subplot(2,1,1); hold on; grid on; box off
        for i = 1:numel(R)
            errorbar(R(i).x, R(i).value, R(i).sterr, 'o-', 'DisplayName', R(i).label, ...
                     'MarkerSize', 3, 'LineWidth', 1);
        end
        xlabel('\tau (\mus)'); ylabel('ratio');
        legend('Location', 'southeast', 'Box', 'off');
        title('Ramsey: classic against walk at several passesPerRun');

        subplot(2,1,2); hold on; grid on; box off
        isC = strcmp({R.label}, 'classic');
        W = R(~isC);
        if ~isempty(W)
            plot([W.P], [W.contrast], 'o-', 'LineWidth', 1.4, 'DisplayName', 'walk');
            if any(isC)
                yline_(R(isC).contrast, 'classic');
            end
            set(gca, 'XScale', 'log');
            xlabel('passesPerRun'); ylabel('fringe contrast');
            legend('Location', 'southwest', 'Box', 'off');
            title('contrast falling with passesPerRun = alignment drift within a pg.run');
        end
    catch err
        fprintf('(plot skipped: %s)\n', err.message);
    end
end

function yline_(y, name)
    xl = xlim;
    plot(xl, [y y], '--', 'LineWidth', 1.2, 'DisplayName', name);
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
