function out = rabiStageTimer(opts)
%RABISTAGETIMER Time every stage of a complete Rabi cell, not just b.run.
%
%   out = rabiStageTimer
%   out = rabiStageTimer(struct('nTau', 97, 'repeats', 3000, 'averages', 3))
%
% Prints a stage-by-stage table and a GRAND TOTAL that should match a stopwatch held
% against the whole cell. Compare the two: anything the stopwatch sees that this does
% not is time being spent somewhere neither of us has thought to look yet.
%
% ===================================================================
% WHY THIS EXISTS: THE LAST PROBE MEASURED THE WRONG INTERVAL
% ===================================================================
%
% rabiTimingProbe wrapped tic/toc around b.run and nothing else. It reported 2.8 to
% 8.8 s per scan while the bench experience is over 30 s, and the discrepancy is not
% mysterious -- everything OUTSIDE b.run was excluded by construction:
%
%   b = AWGExpRabi;                  object construction, device lookups, the
%                                    "Creating a new AWG Rabi Experiment" message
%   b.loadJsonParams('all', ...)     a file read plus dozens of property setters
%   b.nChannels = 2; b.amplitude = ...; b.tau = ...
%                                    every setter raises changeFlag, and some talk to
%                                    the instrument
%   rabiFit(...)                     the FIRST call into the Curve Fitting Toolbox in
%                                    a MATLAB session can cost several seconds on its
%                                    own, before it fits anything
%   figure / errorbar / annotation   plotting is not free either
%
% Timing a fraction of a cell and then reasoning about the whole cell is how I ended up
% recommending a change (the jumpToStep port) that would have saved nothing. So this
% times the whole thing, from before the object exists to after the plot is drawn, and
% attributes every second to a named stage.
%
% ===================================================================
% WHAT IT SEPARATES, AND WHY EACH ONE IS ITS OWN LINE
% ===================================================================
%
%   construct        cheap in principle, and worth confirming
%   loadJsonParams   a file read and many setters. If this is slow, every experiment in
%                    the lab pays it, not just Rabi
%   set parameters   timed individually. set.nChannels re-enters set.amplitude, and
%                    set.tau runs checkTimeVector over the whole vector
%   run, COLD        includes prepare and LoadAWG: the waveform and sequence upload
%   run, WARM        same object, nothing changed, so LoadAWG should be skipped. The
%                    difference between cold and warm IS the load
%   run, no autosave same again with shouldAutosave = false, which prices the .mat write
%   fit              the analysis, separately for the old and the new fitter, because
%                    the toolbox load belongs to the old one and should not be blamed
%                    on the measurement
%   plot             drawing
%
% OPTIONS
%   nTau, repeats, averages   the run to time (defaults 97, 3000, 3 -- a real fine Rabi)
%   tauMax, amplitude, frequency, nChannels, detectionMHz, jsonLocation
%   which          'AWGExpRabi' or 'AWGExpRabiCont' (default 'AWGExpRabi')
%   compareOldFit  also time rabiFit for comparison (default true)
%   doPlot         time the plotting too (default true)
%   checkRepeats   two extra FRESH cold runs at nTauSmall taus with R and 10R repeats
%                  (default true). The load is identical between them because the tau
%                  count is, so the difference is the measurement and nothing else --
%                  which is the only clean way to ask whether `repeats` costs what it
%                  claims. A WARM run cannot answer it: perform() skips taus that
%                  already hold data for the current average, so a second run on the
%                  same object measures nothing at all. That is why the warm run above
%                  comes back in well under a second and must NOT be read as the cost
%                  of the measurement.
%   nTauSmall      taus for those two runs (default 12, to keep the load small)
%
% Leaves the field alone and saves nothing of its own.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('nTau', 97, 'repeats', 3000, 'averages', 3, ...
        'tauMax', 0.5, 'amplitude', -16, 'frequency', [], 'nChannels', 2, ...
        'detectionMHz', 0.3, 'which', 'AWGExpRabi', 'compareOldFit', true, ...
        'doPlot', true, 'checkRepeats', true, 'nTauSmall', 12, ...
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

    S = struct('name', {}, 'seconds', {});
    tAll = tic;

    fprintf('\n=== rabiStageTimer: %d taus, %d repeats, %d averages ===\n', ...
            opts.nTau, opts.repeats, opts.averages);
    fprintf('START a stopwatch now, and compare it with the grand total at the end.\n\n');

    taus = round(linspace(0.02, opts.tauMax, opts.nTau)*1000)/1000;

    % ---- the transport knobs, FIRST --------------------------------------
    % These were printed at the END of the previous version, so the one run where they
    % mattered died before reaching them. Anything that decides how to read everything
    % below it belongs at the top.
    reportTransport();

    % ---- construction ----------------------------------------------------
    t = tic;
    switch lower(opts.which)
        case 'awgexprabicont'; b = AWGExpRabiCont;
        otherwise;             b = AWGExpRabi;
    end
    S = add(S, 'construct object', toc(t));

    % ---- json ------------------------------------------------------------
    t = tic;
    b.loadJsonParams('all', opts.jsonLocation);
    S = add(S, 'loadJsonParams', toc(t));
    if isprop(b, 'useSequenceWalk'); b.useSequenceWalk = false; end   % times the classic path

    % ---- setters, individually -------------------------------------------
    t = tic; b.nChannels = opts.nChannels;              S = add(S, 'set nChannels', toc(t));
    t = tic;
    if opts.nChannels == 2
        b.amplitude = opts.amplitude + [0 0];
    else
        b.amplitude = opts.amplitude;
    end
    S = add(S, 'set amplitude', toc(t));
    if ~isempty(opts.frequency)
        t = tic; b.frequency = opts.frequency;          S = add(S, 'set frequency', toc(t));
    end
    t = tic; b.tau = taus;                              S = add(S, 'set tau', toc(t));
    t = tic;
    b.repeats  = opts.repeats;
    b.averages = opts.averages;
    b.isTracking = 0;
    b.balancedSequence = 1;
    b.detectionDuration = opts.detectionMHz;
    b.referenceDetectionDuration = b.detectionDuration;
    S = add(S, 'set the rest', toc(t));

    % ---- the runs --------------------------------------------------------
    try; b.AWG.resetCommandCount; catch; end
    t = tic; b.run;                                     S = add(S, 'run, COLD (with load)', toc(t));
    cmdCold = safeCount(b);

    % A second run on the SAME object with nothing changed. LoadAWG should be skipped.
    % NOTE this may also skip re-measuring points that already hold data -- perform()
    % checks sum(signal(:,t,currIter)) ~= 0 and continues past filled points -- so read
    % this line as "what a repeat run costs", not as a second independent measurement.
    try; b.AWG.resetCommandCount; catch; end
    t = tic; b.run;                                     S = add(S, 'run, WARM (load skipped)', toc(t));
    cmdWarm = safeCount(b);

    % Same again with the autosave off.
    b.shouldAutosave = false;
    try; b.AWG.resetCommandCount; catch; end
    t = tic; b.run;                                     S = add(S, 'run, WARM, autosave off', toc(t));
    b.shouldAutosave = true;

    % ---- the analysis ----------------------------------------------------
    t = tic;
    try
        % The diagnostics struct is the SIXTH output, not the first -- the first is
        % pulses. Binding one output and then reading .ok is what produced the "Dot
        % indexing is not supported" error attributed to the fitter.
        oNew = rabiFitRobust(struct('exp', b, 'plotIt', false, 'verbose', false));
        newOk = oNew.ok;
    catch err
        newOk = false; fprintf('  (rabiFitRobust failed: %s)\n', err.message);
    end
    S = add(S, 'fit, rabiFitRobust', toc(t));

    if opts.compareOldFit
        t = tic;
        try
            rabiFit(b, 0, 0, 0);
        catch err
            fprintf('  (rabiFit failed: %s)\n', err.message);
        end
        S = add(S, 'fit, rabiFit (old)', toc(t));
        % The Curve Fitting Toolbox is loaded by now, so a second call shows the fit
        % cost without the one-off toolbox load. The difference is worth seeing: it is
        % charged to the first analysis of every MATLAB session.
        t = tic;
        try; rabiFit(b, 0, 0, 0); catch; end
        S = add(S, 'fit, rabiFit (2nd call)', toc(t));
    end

    if opts.doPlot
        t = tic;
        try; rabiFitRobust(struct('exp', b, 'plotIt', true, 'verbose', false)); catch; end
        S = add(S, 'fit + plot', toc(t));
    end

    % ---- does `repeats` cost what it claims? ------------------------------
    rc = struct('ok', false);
    if opts.checkRepeats
        nS = opts.nTauSmall;
        tausS = round(linspace(0.02, opts.tauMax, nS)*1000)/1000;
        fprintf('\n  repeats check: two fresh runs at %d taus, %d and %d repeats\n', ...
                nS, opts.repeats, 10*opts.repeats);
        tR1 = timeFreshCold(opts, tausS, opts.repeats);
        S = add(S, sprintf('cold, %d taus, R', nS), tR1);
        tR2 = timeFreshCold(opts, tausS, 10*opts.repeats);
        S = add(S, sprintf('cold, %d taus, 10R', nS), tR2);
        rc = struct('ok', true, 't1', tR1, 't2', tR2, 'delta', tR2 - tR1);
    end

    total = toc(tAll);

    % ---- report ----------------------------------------------------------
    fprintf('\n%-28s %9s %8s\n', 'stage', 'seconds', 'share');
    fprintf('%s\n', repmat('-', 1, 48));
    for k = 1:numel(S)
        fprintf('%-28s %9.2f %7.1f%%\n', S(k).name, S(k).seconds, 100*S(k).seconds/total);
    end
    fprintf('%s\n', repmat('-', 1, 48));
    fprintf('%-28s %9.2f\n', 'GRAND TOTAL', total);
    fprintf(['\nCompare that with your stopwatch. If the stopwatch is materially larger,\n' ...
             'the missing time is outside every stage above -- the MATLAB event loop, the\n' ...
             'GUI redrawing, or the tracker -- and I have still not found it.\n']);

    cold = seconds(S, 'run, COLD (with load)');
    warm = seconds(S, 'run, WARM (load skipped)');
    noSv = seconds(S, 'run, WARM, autosave off');
    fprintf('\n--- what that says ---\n');
    fprintf('  AWG load           %6.2f s   (cold run minus warm run)\n', cold - warm);
    fprintf('  autosave           %6.2f s   (warm run minus warm-no-autosave)\n', warm - noSv);
    fprintf('  the measurement    %6.2f s   (what is left of the warm run)\n', noSv);
    fprintf('  setup before run   %6.2f s   (construct + json + setters)\n', ...
            sumUpTo(S, 'run, COLD (with load)'));
    fprintf('  analysis + plot    %6.2f s\n', sumFrom(S, 'fit, rabiFitRobust'));
    if opts.compareOldFit
        f1 = seconds(S, 'fit, rabiFit (old)');
        f2 = seconds(S, 'fit, rabiFit (2nd call)');
        if isfinite(f1) && isfinite(f2) && f1 > f2 + 0.5
            fprintf(['  of which %.1f s was the Curve Fitting Toolbox loading on first use --\n' ...
                     '  paid once per MATLAB session, and rabiFitRobust does not pay it at all.\n'], ...
                    f1 - f2);
        end
    end
    if isfinite(cmdCold) && isfinite(cmdWarm)
        fprintf('  SCPI commands      %6d cold, %6d warm\n', cmdCold, cmdWarm);
    end
    if rc.ok
        fprintf('\n--- does `repeats` cost what it claims? ---\n');
        fprintf('  10x the repeats at the same tau count changed the run by %+.2f s\n', rc.delta);
        fprintf('  (the load is identical between those two, so this IS the measurement)\n');
        if abs(rc.delta) < 1
            fprintf(['  *** IT DOES NOT SCALE. 10x the integration cost nothing, so you are\n' ...
                     '      almost certainly not integrating the repeats you ask for. Every\n' ...
                     '      "repeats is free, raise it" argument we have used depends on this,\n' ...
                     '      so it is worth settling: check the SNR against sqrt(repeats), and\n' ...
                     '      check pg.sequenceDuration() against a scope trace of the sequence.\n']);
        else
            fprintf('      Scales, so repeats are real: %.1f us per repeat per tau.\n', ...
                    1e6*rc.delta/(9*opts.repeats*opts.nTauSmall));
        end
    end
    fprintf('\n--- the AWG transport, which is where the load time goes ---\n');
    if isfinite(cmdCold)
        fprintf('  %d commands in the cold run, %.0f ms each on average.\n', ...
                cmdCold, 1e3*(cold - warm)/max(cmdCold, 1));
        if 1e3*(cold-warm)/max(cmdCold,1) > 20
            fprintf(['  *** %.0f ms per SCPI command is a transport problem, not a waveform\n' ...
                     '      problem. Nothing about pushing a few thousand samples justifies it.\n'], ...
                    1e3*(cold-warm)/max(cmdCold,1));
        end
    end
    try; b.AWG.printInstrumentation; catch err
        fprintf('  (printInstrumentation unavailable: %s)\n', err.message);
    end

    out = struct('stages', S, 'total', total, 'repeatsCheck', rc, 'loadSeconds', cold - warm, ...
                 'autosaveSeconds', warm - noSv, 'measureSeconds', noSv, ...
                 'setupSeconds', sumUpTo(S, 'run, COLD (with load)'), ...
                 'fitOk', newOk);
end

function reportTransport()
% Read the three knobs that decide how fast every AWG command is, and say plainly when
% one is in its slow state. main() rebuilds the device objects, so these are lost every
% time it is re-run -- which is exactly what makes one session's timings incomparable
% with another's, and why this is printed FIRST rather than last.
    fprintf('--- AWG transport, before anything else ---\n');
    awg = [];
    try
        awg = getObjByName('TektronixAWG');
        if isempty(awg)
            fgs = FrequencyGenerator.getFG();
            for i = 1:numel(fgs)
                if contains(fgs{i}.name, 'TektronixAWG'); awg = fgs{i}; break; end
            end
        end
    catch
    end
    if isempty(awg)
        fprintf('  (no TektronixAWG found by name)\n\n'); return
    end
    vd = NaN; cdel = NaN; srd = NaN;
    try; vd   = awg.useVisadev;           catch; end
    try; cdel = awg.commDelayTime;        catch; end
    try; srd  = awg.skipRedundantDeletes; catch; end
    fprintf('  useVisadev %s | commDelayTime %s s | skipRedundantDeletes %s\n', ...
            num2str(vd), num2str(cdel), num2str(srd));
    if isequal(vd, false) || isequal(vd, 0)
        fprintf(['  *** useVisadev is FALSE. Every command is on the legacy visa()\n' ...
                 '      interface, which was measured as the single biggest speed\n' ...
                 '      difference on this instrument. Set awg.useVisadev = true.\n']);
    end
    if isnumeric(cdel) && isfinite(cdel) && cdel > 0
        fprintf(['  *** commDelayTime = %g s. At the ~1400 commands a 97-tau load takes,\n' ...
                 '      that is %.1f s of pure busy-waiting. Set awg.commDelayTime = 0.\n'], ...
                cdel, cdel*1400);
    end
    if isequal(srd, false) || isequal(srd, 0)
        fprintf('  *** skipRedundantDeletes is FALSE -- every waveform is deleted twice.\n');
    end
    fprintf('\n');
end

function t = timeFreshCold(opts, taus, reps)
% A brand new object every time, so LoadAWG definitely runs and nothing is skipped
% because it already holds data. perform() continues past any tau that already has a
% value for the current average, so a warm rerun measures nothing -- which is why the
% repeats test has to build a fresh object for each of its two runs.
    t = timeRun(opts, taus, reps, 1, false);
end

function S = add(S, name, sec)
    S(end+1) = struct('name', name, 'seconds', sec);
    fprintf('  %-28s %6.2f s\n', name, sec);
end

function v = seconds(S, name)
    i = find(strcmp({S.name}, name), 1);
    if isempty(i); v = NaN; else; v = S(i).seconds; end
end

function v = sumUpTo(S, name)
    i = find(strcmp({S.name}, name), 1);
    if isempty(i); v = NaN; else; v = sum([S(1:i-1).seconds]); end
end

function v = sumFrom(S, name)
    i = find(strcmp({S.name}, name), 1);
    if isempty(i); v = NaN; else; v = sum([S(i:end).seconds]); end
end

function n = safeCount(b)
    n = NaN;
    try; n = b.AWG.commandCount; catch; end
end

function o = withDefaults(o, d)
    if isempty(o); o = struct(); end
    known = fieldnames(d);
    given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('rabiStageTimer: unknown option "%s". Valid: %s.', bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || isempty(o.(known{i}))
            o.(known{i}) = d.(known{i});
        end
    end
end
