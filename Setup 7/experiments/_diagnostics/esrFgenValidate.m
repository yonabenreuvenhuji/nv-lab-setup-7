function results = esrFgenValidate(a)
%ESRFGENVALIDATE Does FGEN mode give the same spectrum as uploaded waveforms?
%
%   results = esrFgenValidate(a)
%
% Runs the same sweep twice -- once with uploaded waveforms, once with the AWG's
% built-in generator -- and compares them. Always leaves useFgen = false.
%
% CONTRAST IS THE CRITERION HERE, unlike esrWalkValidate where it deliberately is
% not. The reason is what could go wrong in each case. A sequencer that misassigns
% frequencies preserves contrast exactly, which is why the walk and jump checks
% lean on correlation and shape. FGEN cannot misassign anything -- there is no
% sequencer -- but it does leave the microwave on continuously at the switch input,
% so the one plausible failure is leakage through the closed switch putting MW into
% the reference window. That shows up as REDUCED CONTRAST with the shape intact.
% So both are checked, and contrast is the one to read first.
%
% Set a.averages = 1 and give it enough points to resolve the resonance.

    if ~isa(a, 'AWGExpESR')
        error('esrFgenValidate expects an AWGExpESR, got a %s', class(a));
    end
    if length(a.frequency) < 20
        error('Use at least 20 frequency points so the dip is resolved.');
    end

    fprintf('\n=== esrFgenValidate: %d points, repeats = %d ===\n', ...
            length(a.frequency), a.repeats);
    fprintf('timeDelay = %g us', a.timeDelay);
    if a.timeDelay == 0
        fprintf('  (zero, so phase plays no part in this comparison)\n');
    else
        fprintf('  -> %.2f deg at %g MHz in FGEN mode\n', ...
                360*max(a.frequency)*abs(a.timeDelay), max(a.frequency));
    end

    results = struct();

    % Remember the sequencer flags and clear them for the FGEN half. The first
    % version did not, so on a setup where useSequenceWalk was left on -- which is
    % the recommended default now -- the FGEN run hit the mutual-exclusion error
    % instead of running. The AWG-mode reference deliberately KEEPS whatever mode
    % was set, walk or classic, because that is the path being compared against.
    savedWalk = a.useSequenceWalk;
    savedJump = a.useSequenceJump;
    results.referencePath = BooleanHelper.ifTrueElse(savedWalk, 'walk', ...
                            BooleanHelper.ifTrueElse(savedJump, 'jump', 'classic'));
    fprintf('reference path: %s\n', results.referencePath);

    try
        % ---- reference: uploaded waveforms -----------------------------------
        fprintf('\nreference run (uploaded waveforms, %s)...\n', results.referencePath);
        a.useFgen = false;
        t0 = tic;
        a.run;
        results.awgSeconds = toc(t0);
        ref  = a.signalParam.value(:);
        refE = a.signalParam.sterr(:);
        refContrast = 100*(1 - min(ref));
        % Mean level in the REFERENCE window (microwave off). This is the number
        % that tells the two zero-contrast failures apart -- see the note below.
        results.awgRefLevel = mean(a.signal(2, :, 1));
        fprintf('  contrast %.3f%%, dip %.1f sigma, %.1f s, reference level %.1f\n', ...
                refContrast, (1-min(ref))/median(refE), results.awgSeconds, results.awgRefLevel);

        % ---- FGEN ------------------------------------------------------------
        fprintf('\nFGEN run...\n');
        a.useSequenceWalk = false;      % no sequencer exists in FGEN mode
        a.useSequenceJump = false;
        a.useFgen = true;
        t0 = tic;
        a.run;
        results.fgenSeconds = toc(t0);
        fg  = a.signalParam.value(:);
        fgE = a.signalParam.sterr(:);
        a.useFgen = false;
        a.useSequenceWalk = savedWalk;
        a.useSequenceJump = savedJump;

        if isempty(fg) || numel(fg) ~= numel(ref)
            error(['The FGEN run did not produce a full spectrum: got %d values, ' ...
                   'expected %d. Look above for the underlying error. Left ' ...
                   'useFgen = false.'], numel(fg), numel(ref));
        end
        fgContrast = 100*(1 - min(fg));
        results.fgenRefLevel = mean(a.signal(2, :, 1));
        fprintf('  contrast %.3f%%, dip %.1f sigma, %.1f s, reference level %.1f\n', ...
                fgContrast, (1-min(fg))/median(fgE), results.fgenSeconds, results.fgenRefLevel);

        % ---- compare ---------------------------------------------------------
        mm = walkCompareMetrics(ref, refE, fg, fgE);
        results.metrics     = mm;
        results.reference   = ref;
        results.fgen        = fg;
        results.refContrast = refContrast;
        results.fgenContrast = fgContrast;

        % Contrast, in units of its own uncertainty. The dip depth's error is
        % roughly the error on one point, so comparing two independent dips uses
        % the two combined -- otherwise a 5% contrast change at 50 sigma looks
        % alarming when it is ordinary noise.
        dContrast = fgContrast - refContrast;
        sContrast = 100*sqrt(median(refE)^2 + median(fgE)^2);
        results.contrastRatio = fgContrast / refContrast;
        results.contrastSigma = dContrast / max(sContrast, eps);

        fprintf('\n  contrast   AWG %.3f%%  ->  FGEN %.3f%%   (%.0f%% of it, %+.1f sigma)\n', ...
                refContrast, fgContrast, 100*results.contrastRatio, results.contrastSigma);
        fprintf('  correlation %+.4f against a noise ceiling of %+.4f (ratio %.2f)\n', ...
                mm.correlation, mm.rMax, mm.rRatio);
        fprintf('  |diff| / sigma %.2f, best shift %+d\n', mm.sigmaDiff, mm.bestShift);
        fprintf('  time       AWG %.1f s  ->  FGEN %.1f s\n', ...
                results.awgSeconds, results.fgenSeconds);
        fprintf(['             (if the AWG run printed a "Cold load breakdown" line it paid\n' ...
                 '              for the upload too, so compare against a warm run before\n' ...
                 '              drawing a speed conclusion)\n']);

        % Pass needs BOTH: the spectra agree, and the contrast has not been eaten.
        % 0.9 on the ratio is a judgement call rather than a measured threshold --
        % it is there to catch leakage worth worrying about, not to be precise. If
        % it lands just under, look at the sigma figure before believing it.
        results.shapeOk    = mm.rRatio > 0.7 && mm.sigmaDiff < 2 && mm.bestShift == 0;
        results.contrastOk = results.contrastRatio > 0.9 || results.contrastSigma > -2;
        results.pass       = results.shapeOk && results.contrastOk;

        fprintf('\n  %s\n', BooleanHelper.ifTrueElse(results.pass, ...
                'PASS - FGEN mode reproduces the spectrum.', ...
                'FAIL - see below.'));
        if ~results.shapeOk
            fprintf(['  The SHAPE disagrees. FGEN has no sequencer, so this is not a\n' ...
                     '  frequency mapping error. Raise a.fgenSettleTime and re-run: if that\n' ...
                     '  changes nothing, settling is not the cause either.\n']);
        end
        if ~results.contrastOk
            % Near-zero contrast has TWO possible causes and they call for opposite
            % responses, so do not guess between them -- the reference-window level
            % separates them. No microwave at all leaves both windows at full
            % brightness, so the reference level matches the AWG-mode run. Microwave
            % in BOTH windows dims both, so the reference level drops.
            levelRatio = results.fgenRefLevel / max(results.awgRefLevel, eps);
            fprintf('\n  CONTRAST is %.0f%% of the AWG-mode value. Reference level ratio %.3f.\n', ...
                    100*results.contrastRatio, levelRatio);
            if levelRatio > 0.97
                fprintf(['  The reference window is just as bright as in AWG mode, so the NVs\n' ...
                         '  are seeing NO microwave in either window -- this is a plumbing or\n' ...
                         '  output-state problem, not switch leakage. Check the FGEN readback\n' ...
                         '  line above: output state, run state, frequency and amplitude.\n']);
            elseif levelRatio < 0.9
                fprintf(['  The reference window is %.0f%% dimmer than in AWG mode, so microwave\n' ...
                         '  is present in BOTH windows. The switch is not gating the AWG output,\n' ...
                         '  which means the burst itself was the only gate and FGEN cannot work\n' ...
                         '  on this setup without a switch in that path.\n'], 100*(1-levelRatio));
            else
                fprintf(['  The reference level moved a little but not decisively. Repeat with\n' ...
                         '  more repeats before reading anything into it.\n']);
            end
        end

        figure; hold on
        errorbar(a.frequency, ref, refE, '.-');
        errorbar(a.frequency, fg,  fgE,  '.-');
        legend('uploaded waveforms', 'FGEN'); xlabel('MHz'); ylabel('normalized FL');
        title(sprintf('FGEN vs AWG: contrast %.2f%% vs %.2f%%, r = %+.3f', ...
                      fgContrast, refContrast, mm.correlation));

    catch err
        a.useFgen = false;
        a.useSequenceWalk = savedWalk;
        a.useSequenceJump = savedJump;
        fprintf('\n  esrFgenValidate aborted: %s\n', err.message);
        rethrow(err);
    end

    a.useFgen = false;
    a.useSequenceWalk = savedWalk;
    a.useSequenceJump = savedJump;
    fprintf('  Left useFgen = false, useSequenceWalk = %d.\n\n', a.useSequenceWalk);
end
