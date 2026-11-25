function results = esrWalkValidate(a)
%ESRWALKVALIDATE Check that walk mode reproduces the classic spectrum.
%
%   results = esrWalkValidate(a)
%
% 'a' must be a configured AWGExpESR (frequency, repeats, amplitude, detection
% durations already set) with enough points to resolve the resonance. Leaves the
% experiment on the classic path when it finishes, whatever happened.
%
% WHAT IT CHECKS, and why each one is here:
%
% All the metrics come from walkCompareMetrics -- see that file for what each one
% means and why. The important design point: agreement is judged against the
% NOISE-LIMITED correlation ceiling, not against 1.
%
% Two correct but noisy measurements of the same spectrum cannot correlate
% arbitrarily well; their own error bars cap it. At repeats = 20 the ceiling is
% around 0.8, so an earlier flat "r > 0.9" test failed walk mode three times over
% while its plotted spectrum visibly matched the reference. rRatio (observed r
% divided by that ceiling) and sigmaDiff (agreement in units of sigma) are the
% honest criteria; raw correlation is reported but not judged in isolation.
%
% Contrast is deliberately NOT a criterion: a frequency misassignment preserves it
% exactly. The jump path produced the right 2.16% contrast at r = +0.04.

    if ~isa(a, 'AWGExpESR')
        error('esrWalkValidate expects an AWGExpESR, got a %s', class(a));
    end
    if length(a.frequency) < 20
        error(['Use at least 20 frequency points -- a short sweep cannot ' ...
               'distinguish a shifted or scrambled spectrum from a correct one.']);
    end

    awg = a.AWG;
    fprintf('\n=== esrWalkValidate: %d points, repeats = %d ===\n', ...
            length(a.frequency), a.repeats);

    results = struct();

    try
        % ---- reference: the classic path -------------------------------------
        fprintf('\nreference run (classic path)...\n');
        a.useSequenceWalk = false;
        a.useSequenceJump = false;
        awg.resetCommandCount;
        a.run;
        ref  = a.signalParam.value(:);
        refE = a.signalParam.sterr(:);
        tRef = awg.timing;
        fprintf('  contrast %.2f%%, dip %.1f sigma, %.1f s\n', ...
                100*(1-min(ref)), (1-min(ref))/median(refE), toc(tRef.wallStart));

        % ---- walk path -------------------------------------------------------
        fprintf('\nwalk run...\n');
        a.useSequenceWalk = true;
        awg.resetCommandCount;
        a.run;
        walk  = a.signalParam.value(:);
        walkE = a.signalParam.sterr(:);
        tWalk = awg.timing;
        walkSeconds = toc(tWalk.wallStart);

        % Check the walk run actually produced a spectrum before computing any
        % metrics. If it errored mid-sweep, signalParam is empty or short, and
        % the comparison below would fail with an unhelpful size error instead of
        % telling you the run itself broke.
        if isempty(walk) || numel(walk) ~= numel(ref)
            a.useSequenceWalk = false;
            error(['The walk run did not produce a full spectrum: got %d values, ' ...
                   'expected %d. Look further up for the underlying error -- ' ...
                   'Experiment.run catches per-average failures and only warns. ' ...
                   'Left useSequenceWalk = false.'], numel(walk), numel(ref));
        end

        fprintf('  contrast %.2f%%, %.1f s, %d SCPI commands\n', ...
                100*(1-min(walk)), walkSeconds, awg.commandCount);

        % ---- metrics ---------------------------------------------------------
        % Judged against the NOISE-LIMITED correlation ceiling, not against 1.
        % See walkCompareMetrics: two correct but noisy measurements of the same
        % spectrum cannot reach r = 0.9 when their own error bars cap them lower,
        % and a flat 0.9 threshold produced repeated false failures here.
        mm = walkCompareMetrics(ref, refE, walk, walkE);
        results = mm;
        results.reference   = ref;
        results.walk        = walk;
        results.walkSeconds = walkSeconds;

        fprintf('\n  correlation            : %+.4f\n', mm.correlation);
        fprintf('  noise-limited ceiling  : %+.4f   (the best two noisy runs could do)\n', mm.rMax);
        fprintf('  ratio to ceiling       : %.2f     (>0.7 is agreement)\n', mm.rRatio);
        fprintf('  best shift             : %+d      (r = %+.4f)\n', mm.bestShift, mm.bestShiftCorr);
        fprintf('  |diff| / sigma         : %.2f     (~0.8 for a correct pair)\n', mm.sigmaDiff);
        fprintf('  scatter ratio          : %.2f     (~1.4 for a real spectrum)\n', mm.scatterRatio);
        fprintf('  first point            : %.1f sigma\n', mm.firstPointSigma);

        fprintf('\n  %s\n', BooleanHelper.ifTrueElse(mm.pass, ...
                'PASS - walk mode reproduces the classic spectrum.', ...
                'FAIL - do not use walk mode yet.'));

        if ~mm.pass
            if mm.bestShift ~= 0 && mm.bestShiftCorr > 0.7*mm.rMax
                fprintf('  The spectrum is CORRECT but SHIFTED by %+d points -- the sweep is\n', mm.bestShift);
                fprintf('  gaining or losing %d trigger(s) per pass. An exact, fixable off-by-%d.\n', ...
                        abs(mm.bestShift), abs(mm.bestShift));
            elseif mm.scatterRatio < 0.8
                fprintf('  The walk spectrum is SMOOTHER than its own error bars, which means it\n');
                fprintf('  is stuck on one frequency rather than scrambled across many.\n');
            elseif mm.rRatio < 0.3
                fprintf('  No shift recovers it, so the sequencer is not advancing one block per\n');
                fprintf('  trigger.\n');
            end
        end
        if mm.firstPointSigma > 4
            fprintf('  NOTE the first frequency disagrees by %.1f sigma. The classic path drops\n', mm.firstPointSigma);
            fprintf('  the first repeat AT EVERY frequency; walk mode drops the first PASS, so\n');
            fprintf('  frequency 1 is the one point that can keep an un-discarded warm-up.\n');
        end

        figure; hold on
        errorbar(a.frequency, ref,  refE,  '.-');
        errorbar(a.frequency, walk, walkE, '.-');
        legend('classic', 'walk'); xlabel('MHz'); ylabel('normalized FL');
        title(sprintf('walk vs classic: r = %+.3f (ceiling %+.3f), shift %+d', ...
                      mm.correlation, mm.rMax, mm.bestShift));

    catch err
        a.useSequenceWalk = false;
        fprintf('\n  esrWalkValidate aborted: %s\n', err.message);
        rethrow(err);
    end

    a.useSequenceWalk = false;      % never leave the setup on an unproven path
    fprintf('  Left useSequenceWalk = false.\n\n');
end

function r = corr2col(x, y)
    x = x(:) - mean(x(:));
    y = y(:) - mean(y(:));
    d = sqrt(sum(x.^2) * sum(y.^2));
    if d == 0
        r = 0;
    else
        r = sum(x .* y) / d;
    end
end
