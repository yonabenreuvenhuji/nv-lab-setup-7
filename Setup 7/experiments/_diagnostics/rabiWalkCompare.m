function out = rabiWalkCompare(opts)
%RABIWALKCOMPARE Classic against walk: same settings, both timed, both fitted.
%
%   out = rabiWalkCompare
%   out = rabiWalkCompare(struct('nTau', 97, 'repeats', 300, 'averages', 1))
%
% Runs the SAME Rabi twice, once on each path, and reports two things that have to be
% checked together: how long each took, and whether they agree on the answer. A faster
% path that measures something slightly different is not faster, it is broken.
%
% ===================================================================
% THE STATISTICS ARE IDENTICAL, WHICH IS WHY THIS IS A FAIR TEST
% ===================================================================
%
% `repeats` means different things on the two paths and yet the photon count comes out
% the same, so no scaling is needed:
%
%   classic  pg.repeats = repeats, and pg.run happens once per tau
%              -> repeats triggers per tau per average
%   walk     pg.repeats = (nTau + warmup) * passesPerRun, one trigger per tau per pass,
%            and nRuns = ceil(repeats / passesPerRun)
%              -> repeats triggers per tau per average
%
% Same triggers, same integration, same sequence. What differs is only the ORDER they
% arrive in: classic takes a tau's repeats back to back, walk spreads them across
% passes. Spreading is better for drift, which is the one respect in which the two are
% not equivalent -- and it favours walk, so it cannot flatter the comparison.
%
% ===================================================================
% WHAT TO LOOK AT, IN ORDER
% ===================================================================
%
% 1. THE PHYSICS, FIRST. pi time and contrast, with the difference quoted in sigma.
%    Walk holds the MW gate fixed at max(tau) and lets the AWG waveform define the
%    pulse, so the switch stays open while the AWG idles after the burst. If that tail
%    leaks, walk reads a LOWER contrast at the same drive. Two sigma or less on both
%    numbers and the paths are measuring the same thing.
%
% 2. THE CLOCK. Wall time and SCPI command count. The command count is the cleaner
%    number because it is immune to whatever else the PC is doing.
%
% 3. THE WAVEFORM COUNT, which walk mode does NOT reduce and which may now dominate.
%    See the note printed at the end: it is set by the tau STEP, not by the number of
%    taus, and that is not obvious.
%
% OPTIONS
%   nTau, repeats, averages, tauMax    the run (defaults 99, 300, 2, 1.0 -- the same as
%                                      the walk cell in runRabiBlock)
%   tau        pass the tau vector explicitly to compare on your own grid
%   amplitude, frequency, nChannels, detectionMHz, jsonLocation
%   passesPerRun, walkWarmupTriggers   walk settings (defaults 25, 1). passesPerRun = 1
%                                      makes walk re-arm once per PASS and it will lose
%                                      badly -- see the re-arm arithmetic printed at the end
%   autosave    leave the .mat autosave on (default false -- it measured ~6 s per run
%               and would sit identically in both numbers, adding noise and nothing else)
%
% Saves nothing and changes no calibration.

    if ~exist('opts', 'var'); opts = struct(); end
    % DEFAULTS MATCHED TO THE WALK CELL IN runRabiBlock, because a comparison run at
    % different settings from the thing being compared is not a comparison. In
    % particular averages was 1 here and 2 there, and that alone explains why both
    % fits looked wrong in the comparison and fine in the block: with a single average
    % sterr is not a standard error at all -- it is the within-average scatter with
    % nothing to compare it against -- so every point's weight is wrong and the fit is
    % pulled by whichever points happened to get the smallest bar.
    opts = withDefaults(opts, struct('nTau', 99, 'repeats', 300, 'averages', 2, ...
        'tauMax', 1.0, 'amplitude', -16, 'frequency', [], 'nChannels', 2, ...
        'detectionMHz', 0.3, 'passesPerRun', 25, 'walkWarmupTriggers', 1, ...
        'autosave', false, 'tau', [], ...
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

    % Pass tau explicitly to compare on exactly the grid you run day to day. The
    % default reproduces the block's 0.02:0.01:tauMax, which also keeps the residue
    % waveform count at 0.2/step rather than one per tau.
    if isempty(opts.tau)
        taus = round((0.02 : (opts.tauMax-0.02)/(opts.nTau-1) : opts.tauMax)*1000)/1000;
    else
        taus = round(opts.tau(:)'*1000)/1000;
        opts.nTau = numel(taus);
    end
    fprintf('\n=== rabiWalkCompare: %d taus, repeats %d, averages %d, passesPerRun %d ===\n', ...
            opts.nTau, opts.repeats, opts.averages, opts.passesPerRun);
    if opts.averages < 2
        fprintf(['WARNING averages = 1. sterr is then not a standard error, so both fits\n' ...
                 '        get wrong weights and both will look off. Use 2 or more.\n']);
    end
    fprintf('Both paths get the same repeats, which is the same number of triggers per\n');
    fprintf('tau on each -- see the header. Nothing is scaled between them.\n\n');

    % ---- classic ---------------------------------------------------------
    fprintf('--- classic path ---\n');
    [tC, bC, nC] = runOne(opts, taus, false);
    fC = fitOne(bC);
    fprintf('  %.1f s, %d SCPI commands\n', tC, nC);

    % ---- walk ------------------------------------------------------------
    fprintf('\n--- walk path ---\n');
    [tW, bW, nW] = runOne(opts, taus, true);
    fW = fitOne(bW);
    fprintf('  %.1f s, %d SCPI commands\n', tW, nW);

    % ---- the physics -----------------------------------------------------
    fprintf('\n--- do they agree? ---\n');
    fprintf('%-18s %14s %14s %12s\n', 'quantity', 'classic', 'walk', 'difference');
    rows = {'f',        'MHz',  'f',        'fErr'; ...
            'pi time',  'ns',   'piTime',   'piErr'; ...
            'contrast', '%',     'contrast', 'contrastErr'};
    worst = 0; agree = true;
    for k = 1:size(rows, 1)
        vC = fC.(rows{k,3}); eC = fC.(rows{k,4});
        vW = fW.(rows{k,3}); eW = fW.(rows{k,4});
        d  = vW - vC;
        sd = sqrt(eC^2 + eW^2);
        nSig = abs(d)/max(sd, eps);
        worst = max(worst, nSig);
        fprintf('%-18s %8.3f%-6s %8.3f%-6s %+8.3f (%.1f sigma)\n', ...
                [rows{k,1} ' /' rows{k,2}], vC, '', vW, '', d, nSig);
        if nSig > 2; agree = false; end
    end
    fprintf('  R^2                %8.4f       %8.4f\n', fC.r2, fW.r2);
    fprintf('  reduced chi2       %8.2f       %8.2f\n', fC.chi2, fW.chi2);

    if agree
        fprintf(['\n  The two paths agree to %.1f sigma on everything, so walk is measuring\n' ...
                 '  the same Rabi. In particular the contrast holds up, which is the test\n' ...
                 '  that the fixed MW gate is not leaking through the switch during the\n' ...
                 '  tail after the burst.\n'], worst);
    else
        fprintf('\n  *** THEY DISAGREE at %.1f sigma.\n', worst);
        if fW.contrast > fC.contrast + 2*sqrt(fC.contrastErr^2 + fW.contrastErr^2)
            % Walk HIGHER, which is the interesting direction and not the one to fear.
            fprintf(['  Walk reads a HIGHER contrast than classic (%.2f%% vs %.2f%%). That is not\n' ...
                     '  walk inventing signal -- it is more likely classic LOSING some. Classic\n' ...
                     '  re-arms at every tau (assignSequence stops playback), and any repeats that\n' ...
                     '  fire before the AWG resumes come back at ratio ~= 1 and dilute the\n' ...
                     '  contrast. That dilution is a fixed time per point, so it hurts more at low\n' ...
                     '  repeats. Test it: run classic at 300, 1000 and 3000 repeats -- if the\n' ...
                     '  contrast climbs toward walk''s, the dilution is real and walk is the more\n' ...
                     '  accurate path, not just the faster one.\n'], fW.contrast, fC.contrast);
        elseif fC.contrast > fW.contrast
            fprintf(['  Walk reads a LOWER contrast. Suspect the fixed MW gate: the switch stays\n' ...
                     '  open while the AWG idles after the burst, and anything leaking through\n' ...
                     '  that tail reaches both windows.\n']);
        else
            fprintf(['  The frequency moved but the contrast did not. Suspect the tau mapping:\n' ...
                     '  walk plays blocks in obj.walkOrder and performWalk maps them back, so an\n' ...
                     '  aborted pass with resyncEveryPass off would shift it.\n']);
        end
    end

    % ---- the re-arm arithmetic, which is what actually sets the time ------
    % Measured from two independent runs: a classic point-visit costs 0.57 to 0.70 s,
    % almost all of it the assignSequence-stops-playback re-arm the AWG class documents
    % as ~509 ms. So the run time is set by HOW MANY re-arms happen, and that is
    % arithmetic, not mystery:
    %
    %   classic  nTau * averages                        re-arms
    %   walk     ceil(repeats/passesPerRun) * averages   re-arms
    %
    % Walk only wins when passesPerRun is large enough that the second number is the
    % smaller one. At passesPerRun = 1 it is repeats, which is usually far WORSE.
    reC = opts.nTau * opts.averages;
    reW = ceil(opts.repeats/opts.passesPerRun) * opts.averages;
    fprintf('\n--- re-arms, at ~0.6 s each ---\n');
    fprintf('  classic  %d taus x %d averages          = %4d re-arms (~%.0f s)\n', ...
            opts.nTau, opts.averages, reC, 0.6*reC);
    fprintf('  walk     ceil(%d/%d) x %d averages      = %4d re-arms (~%.0f s)\n', ...
            opts.repeats, opts.passesPerRun, opts.averages, reW, 0.6*reW);
    fprintf('  walk digitizer load: 2*(%d+%d)*%d = %d segments per pg.run\n', ...
            opts.nTau, opts.walkWarmupTriggers, opts.passesPerRun, ...
            2*(opts.nTau+opts.walkWarmupTriggers)*opts.passesPerRun);
    if reW >= reC
        fprintf(['  *** walk is configured to re-arm MORE than classic. Raise passesPerRun\n' ...
                 '      to at least %d before drawing any conclusion about walk.\n'], ...
                max(ceil(opts.repeats/max(reC-1,1)), 2));
    end

    % ---- the clock -------------------------------------------------------
    fprintf('\n--- the clock ---\n');
    fprintf('  classic %7.1f s, %6d commands\n', tC, nC);
    fprintf('  walk    %7.1f s, %6d commands\n', tW, nW);
    if tW < tC
        fprintf('  walk is %.1fx faster (%.0f s saved), and sends %.1fx fewer commands\n', ...
                tC/max(tW, eps), tC - tW, nC/max(nW, 1));
    else
        fprintf('  walk is NOT faster here (%.1fx). Read the waveform note below.\n', tW/max(tC,eps));
    end

    % ---- what walk does NOT fix -----------------------------------------
    % The number of RESIDUE waveforms is set by the tau STEP against the 0.2 us base
    % waveform, not by the number of taus, and that is not obvious from the code. Each
    % tau needs a waveform for mod(tau, 0.2), so a uniform grid of step s produces
    % 0.2/s distinct residues however many taus there are -- and one per channel.
    step = median(diff(taus));
    nRes = numel(unique(round(mod(taus, 0.2)*1e6)/1e6));
    fprintf('\n--- what walk does NOT fix: the waveform count ---\n');
    fprintf('  tau step %.4f us -> %d distinct residues -> %d waveforms (x%d channels)\n', ...
            step, nRes, nRes*opts.nChannels + 2*opts.nChannels, opts.nChannels);
    fprintf(['  Walk replaces %d AddSequence calls with %d, but the waveforms are unchanged.\n' ...
             '  Their count follows the STEP, not the number of taus: each tau needs a\n' ...
             '  waveform for mod(tau, 0.2 us), so a uniform grid of step s needs 0.2/s of\n' ...
             '  them whatever nTau is. Doubling the step to %.4f us would need %d residues\n' ...
             '  instead of %d, at %d taus over the same range -- half the upload for half\n' ...
             '  the points, and 12 points per period is already plenty.\n'], ...
            opts.nChannels*opts.nTau, opts.nChannels, 2*step, ...
            max(round(0.2/(2*step)), 1), nRes, ceil(opts.nTau/2));

    out = struct('classicSeconds', tC, 'walkSeconds', tW, ...
                 'classicCommands', cmdCount(bC), 'walkCommands', cmdCount(bW), ...
                 'classic', fC, 'walk', fW, 'agree', agree, 'worstSigma', worst, ...
                 'residues', nRes, 'tauStep', step);
end

% ======================================================================

function [t, b, nCmd] = runOne(opts, taus, useWalk)
% The command count is read IMMEDIATELY after the run, not later. Both experiment
% objects share the same AWG object, so resetCommandCount in the second run resets the
% counter the first run's number would have been read from -- which is why the summary
% previously reported the classic run's count as NaN, or worse, as the walk run's.
    b = AWGExpRabi;
    b.loadJsonParams('all', opts.jsonLocation);
    b.nChannels = opts.nChannels;                 % BEFORE amplitude, always
    if opts.nChannels == 2
        b.amplitude = opts.amplitude + [0 0];
    else
        b.amplitude = opts.amplitude;
    end
    if ~isempty(opts.frequency); b.frequency = opts.frequency; end
    b.repeats  = opts.repeats;
    b.averages = opts.averages;
    b.isTracking = 0;
    b.balancedSequence = 1;
    b.detectionDuration = opts.detectionMHz;
    b.referenceDetectionDuration = b.detectionDuration;
    b.tau = taus;
    b.shouldAutosave = opts.autosave;
    if useWalk
        b.useSequenceWalk     = true;
        b.passesPerRun        = opts.passesPerRun;
        b.walkWarmupTriggers  = opts.walkWarmupTriggers;
        b.resyncEveryPass     = true;
    end
    try; b.AWG.resetCommandCount; catch; end
    % ID FORM, tId = tic ... toc(tId), NOT bare tic/toc.
    %
    % A bare tic sets a single global stopwatch, so ANY tic called inside b.run resets
    % it and the following toc measures from that internal tic instead of from here.
    % That is not hypothetical: this function reported 1.9 s for a run whose own log
    % showed 18.5 s of sequence settle and 16.4 s to arm, and the bench experience was
    % two minutes. Every number rabiTimingProbe produced was wrong for the same reason.
    % rabiStageTimer happened to use the id form and its 264 s is the one measurement of
    % the three that can be trusted.
    tId = tic; b.run; t = toc(tId);
    nCmd = cmdCount(b);
end

function f = fitOne(b)
% One fit, reduced to the four numbers the comparison needs. rabiFitRobust rather than
% rabiFit, because a comparison is worthless if the fitter can land on a bound: the old
% one returned exactly 1.5x its start point on this data.
    f = struct('f', NaN, 'fErr', NaN, 'piTime', NaN, 'piErr', NaN, ...
               'contrast', NaN, 'contrastErr', NaN, 'chi2', NaN, 'r2', NaN);
    try
        % SIX outputs, and the struct is the SIXTH. Calling it with one output binds
        % `pulses` -- a 1x3 double -- and then o.ok raises "Dot indexing is not supported
        % for variables of this type". That error was blamed on the fitter's data
        % extraction for three rounds; it was this line the whole time.
        o = rabiFitRobust(struct('exp', b, 'plotIt', false, 'verbose', false));
        if o.ok
            f.f = o.f;                  f.fErr = o.fErr;
            f.piTime = 1e3*o.pulses(2); f.piErr = 1e3*o.pulsesErr(2);
            f.contrast = 100*o.contrast; f.contrastErr = 100*o.contrastErr;
            f.chi2 = o.chi2Reduced; f.r2 = o.r2;
        else
            fprintf('  fit failed: %s\n', o.why);
        end
    catch err
        fprintf('  fit threw: %s\n', err.message);
    end
end

function n = cmdCount(b)
    n = NaN;
    try; n = b.AWG.commandCount; catch; end
end

function o = withDefaults(o, d)
    if isempty(o); o = struct(); end
    known = fieldnames(d);
    given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('rabiWalkCompare: unknown option "%s". Valid: %s.', bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || isempty(o.(known{i}))
            o.(known{i}) = d.(known{i});
        end
    end
end
