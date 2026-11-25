function results = esrJumpUnscramble(a, opts)
%ESRJUMPUNSCRAMBLE Is jump mode's scrambling a fixed, invertible permutation?
%
%   results = esrJumpUnscramble(a)
%   results = esrJumpUnscramble(a, struct('lagRange', -5:5))
%
% Four tests on one question: jump mode collects the right contrast and files it under
% the wrong frequencies. If that misfiling is a DETERMINISTIC permutation, it can be
% inverted in software and jump mode becomes usable with no new hardware at all. If it
% differs run to run, it is a race and no amount of reordering will ever fix it.
%
% ===================================================================
% WHY NOT JUST MATCH THE VALUES AND REORDER
% ===================================================================
%
% The obvious test -- find which jump point matches which reference point, reorder, see
% if it agrees -- CANNOT FAIL, and that is its problem. Fitting a permutation to data
% and then checking that the fitted permutation fits is circular. With 96 points there
% are 96 assignment choices to play with, and an optimal matching will reproduce almost
% any target curve from almost any set of values, including pure noise. A test that
% passes on noise tells you nothing when it passes on your data.
%
% There is a second, physical degeneracy on top of the statistical one: an ESR dip is
% symmetric, so f0 - d and f0 + d have the SAME fluorescence. Value matching cannot
% distinguish the two flanks even in principle. Half the ordering information simply is
% not in the values.
%
% So the tests below are built to be falsifiable instead:
%
%   1. DISTRIBUTION (permutation-invariant, no fitting). A permutation preserves the
%      multiset of values, so sorted(jump) must equal sorted(reference). This one is
%      already known to pass -- 2.16% contrast against 2.10%, min 0.9784, median 0.9921
%      -- and it is repeated here because it is the precondition for everything else.
%      If it FAILS, the values are not a permutation of the truth, no reordering can
%      help, and the other three tests are meaningless.
%
%   2. REPRODUCIBILITY (no fitting at all, and the decisive one). Run jump twice and
%      correlate the two runs POINT BY POINT. No permutation is inferred; the two runs
%      are simply compared. r near 1 means the scramble is the same every time, so it
%      is deterministic and invertible. r near 0 means it is a race, and jump mode is
%      dead however good the other numbers look.
%
%   3. LAG (one integer, tested over a small range -- a hypothesis, not a fit). The
%      diagnosed suspect is that the value recorded at visit k belongs to the frequency
%      visited at k-1. The visit order is now recorded in a.lastVisitOrder, so that
%      remap is fully determined and has NO free parameters beyond the integer lag
%      itself. If lag = -1 lifts the correlation from 0.04 to ~0.99, the diagnosis is
%      proved and the fix is an index shift.
%
%   4. HELD-OUT MATCHING (the original idea, made honest). Infer the assignment from
%      jump run A, then apply that same assignment to jump run B and compare B with the
%      reference. The permutation is fitted on one dataset and tested on another, so
%      agreement means something. Noise cannot pass this: a permutation fitted to
%      noise in A is unrelated to the noise in B.
%
% ===================================================================
% BEFORE RUNNING THIS
% ===================================================================
%
% Set a.randomizeFrequencyOrder = false if you want test 3 to also show up in the
% ordinary best-shift metric: in random visit order a one-step lag scatters the spectrum
% and no frequency-space shift describes it, which is why "best shift" always reported 0.
% Sequential order turns the same lag into a plain shift. The tests here work either way,
% because they use the recorded order rather than assuming one.
%
% Field OFF, so the dip sits at 2870.
%
% ===================================================================
% OPTIONS (struct fields, all optional)
%   lagRange     integer lags to test in visit order (default -4:4)
%   nJump        jump runs to take (default 2; 3 sharpens test 2)
%   referencePath 'walk' or 'classic' (default: whichever a has set, else walk)
%   plotIt       draw the comparison (default true)
%
% Leaves useSequenceJump = false and restores every flag it touched, including on error.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('lagRange', -4:4, 'nJump', 2, ...
                                     'referencePath', '', 'plotIt', true));
    if opts.nJump < 2
        error('nJump must be at least 2 -- test 2 needs two independent jump runs.');
    end

    savedJump = a.useSequenceJump;
    savedWalk = a.useSequenceWalk;
    savedFgen = a.useFgen;
    f = a.frequency(:);
    n = numel(f);

    fprintf('\n=== esrJumpUnscramble: %d points ===\n', n);
    fprintf('visit order: %s\n', BooleanHelper.ifTrueElse(a.randomizeFrequencyOrder, ...
            'RANDOM (a lag will not appear as a frequency shift)', ...
            'SEQUENTIAL (a lag appears as a plain frequency shift)'));

    results = struct();
    try
        % ---- reference ------------------------------------------------------
        a.useFgen = false;
        a.useSequenceJump = false;
        if isempty(opts.referencePath)
            a.useSequenceWalk = savedWalk;
        else
            a.useSequenceWalk = strcmpi(opts.referencePath, 'walk');
        end
        fprintf('\nreference run (%s)...\n', ...
                BooleanHelper.ifTrueElse(a.useSequenceWalk, 'walk', 'classic'));
        a.run;
        ref  = a.signalParam.value(:);
        refE = a.signalParam.sterr(:);
        refOrder = a.lastVisitOrder;
        fprintf('  contrast %.3f%%\n', 100*(1 - min(ref)));

        % ---- jump runs ------------------------------------------------------
        a.useSequenceWalk = false;
        a.useSequenceJump = true;
        J = nan(n, opts.nJump); JE = nan(n, opts.nJump); JO = nan(n, opts.nJump);
        for r = 1:opts.nJump
            fprintf('jump run %d of %d...\n', r, opts.nJump);
            a.run;
            J(:, r)  = a.signalParam.value(:);
            JE(:, r) = a.signalParam.sterr(:);
            ord = a.lastVisitOrder;
            if isempty(ord)
                error(['a.lastVisitOrder is empty. This test needs the patched ' ...
                       'AWGExpESR that records the visit order.']);
            end
            JO(:, r) = ord(:);
            fprintf('  contrast %.3f%%, correlation with reference %+.4f\n', ...
                    100*(1 - min(J(:, r))), corrSafe(J(:, r), ref));
        end
        a.useSequenceJump = false;

        results.frequency = f;
        results.reference = ref;  results.referenceErr = refE;
        results.jump = J;         results.jumpErr = JE;
        results.visitOrder = JO;  results.referenceVisitOrder = refOrder;

        % =================================================================
        % TEST 1 -- distribution. Permutation-invariant, nothing fitted.
        % =================================================================
        % Sorted values, compared against the noise you would expect anyway. A
        % permutation cannot change the multiset, so disagreement here is fatal to the
        % whole idea and everything below would be meaningless.
        fprintf('\n--- 1. distribution (is it even a permutation?) ---\n');
        sr = sort(ref); sj = sort(J(:, 1));
        dSorted = sj - sr;
        sigExp  = sqrt(median(refE)^2 + median(JE(:,1))^2);
        results.sortedRms    = sqrt(mean(dSorted.^2));
        results.sortedRmsSig = results.sortedRms / max(sigExp, eps);
        fprintf('  sorted values differ by %.5f rms  (%.1f x the expected noise %.5f)\n', ...
                results.sortedRms, results.sortedRmsSig, sigExp);
        fprintf('  contrast   reference %.3f%%  jump %.3f%%\n', ...
                100*(1-min(ref)), 100*(1-min(J(:,1))));
        fprintf('  median     reference %.5f  jump %.5f\n', median(ref), median(J(:,1)));
        results.isPermutation = results.sortedRmsSig < 3;
        fprintf('  => %s\n', BooleanHelper.ifTrueElse(results.isPermutation, ...
                'the values ARE a permutation of the reference. Reordering can work.', ...
                ['the values are NOT a permutation. Jump is losing or adding signal, ' ...
                 'not just misfiling it -- reordering cannot fix this.']));

        % =================================================================
        % TEST 2 -- reproducibility. The decisive test, and nothing is fitted.
        % =================================================================
        fprintf('\n--- 2. reproducibility (is the scramble the same every time?) ---\n');
        rJJ = nan(opts.nJump, opts.nJump);
        for p = 1:opts.nJump
            for q = 1:opts.nJump
                rJJ(p, q) = corrSafe(J(:, p), J(:, q));
            end
        end
        offDiag = rJJ(~eye(opts.nJump));
        results.jumpToJumpR = median(offDiag);
        % The ceiling: two identical measurements still disagree by their noise, so a
        % perfect correlation is not 1. Estimated the same way walkCompareMetrics does.
        results.jumpToJumpCeiling = corrCeiling(ref, median(JE, 2));
        fprintf('  jump-to-jump correlation %+.4f against a ceiling of %+.4f\n', ...
                results.jumpToJumpR, results.jumpToJumpCeiling);
        sameOrder = results.jumpToJumpR > 0.7 * results.jumpToJumpCeiling;
        results.deterministic = sameOrder;
        if sameOrder
            fprintf(['  => DETERMINISTIC. The same scramble happens every run, so it is a\n' ...
                     '     fixed permutation and can be inverted in software. No flag\n' ...
                     '     wiring needed.\n']);
        else
            fprintf(['  => NOT REPRODUCIBLE. The scramble differs between runs, so it is a\n' ...
                     '     race, not a permutation. Jump mode cannot be rescued by\n' ...
                     '     reordering, and a flag output would only tell you how badly it\n' ...
                     '     failed each time. This is where I would stop.\n']);
        end

        % =================================================================
        % TEST 3 -- lag in visit order. One integer, over a small range.
        % =================================================================
        % If the value filed under frequency visitOrder(k) actually belongs to
        % visitOrder(k+lag), undo it. The visit order is measured, so this has no free
        % continuous parameters -- just a small integer, scanned so the shape of the
        % result is visible rather than a single number to believe.
        fprintf('\n--- 3. lag in visit order (the diagnosed suspect) ---\n');
        lags = opts.lagRange(:)';
        rLag = nan(size(lags));
        for il = 1:numel(lags)
            unl = applyVisitLag(J(:, 1), JO(:, 1), lags(il));
            rLag(il) = corrSafe(unl, ref);
        end
        results.lags = lags; results.lagR = rLag;
        [rBest, ib] = max(rLag);
        results.bestLag = lags(ib); results.bestLagR = rBest;
        for il = 1:numel(lags)
            fprintf('    lag %+d : r = %+.4f%s\n', lags(il), rLag(il), ...
                    ternary(il == ib, '   <-- best', ''));
        end
        if rBest > 0.9
            fprintf(['  => EXPLAINED. A lag of %+d in visit order lifts the correlation\n' ...
                     '     from %+.4f to %+.4f. That is a one-line index fix, and it needs\n' ...
                     '     no settle time and no hardware.\n'], ...
                    results.bestLag, rLag(lags == 0), rBest);
        else
            fprintf(['  => not a simple lag. Best is %+d at r = %+.4f. The misfiling is\n' ...
                     '     something more structured than an off-by-k in visit order.\n'], ...
                    results.bestLag, rBest);
        end

        % =================================================================
        % TEST 4 -- held-out matching. The original idea, made falsifiable.
        % =================================================================
        % Infer the assignment from jump run 1, APPLY IT TO RUN 2, and compare run 2
        % with the reference. Fitted on A, tested on B. Noise cannot pass: a permutation
        % fitted to run 1's noise has nothing to do with run 2's noise.
        fprintf('\n--- 4. held-out matching (fit the order on run 1, test it on run 2) ---\n');
        perm = greedyAssign(J(:, 1), ref);
        inSample  = corrSafe(J(perm, 1), ref);
        outSample = corrSafe(J(perm, 2), ref);
        results.matchPerm = perm;
        results.matchInSample = inSample;
        results.matchOutSample = outSample;
        fprintf('  in-sample  (fitted and tested on run 1) r = %+.4f  <-- always high, ignore\n', ...
                inSample);
        fprintf('  OUT-OF-SAMPLE (fitted on 1, tested on 2) r = %+.4f  <-- this is the number\n', ...
                outSample);
        % A control: the same procedure against SHUFFLED reference values. Whatever this
        % scores is what the method returns from no real structure at all, so the
        % out-of-sample number has to beat it to mean anything.
        rc = nan(1, 20);
        for t = 1:20
            shuf = ref(randperm(n));
            pc = greedyAssign(J(:, 1), shuf);
            rc(t) = corrSafe(J(pc, 2), shuf);
        end
        results.matchControl = median(rc);
        fprintf('  control (same method on shuffled reference) r = %+.4f\n', results.matchControl);
        if outSample > 0.9 && outSample > results.matchControl + 0.3
            fprintf('  => the recovered order is REAL and transfers between runs.\n');
        else
            fprintf(['  => does not transfer. The matching was fitting noise, which is\n' ...
                     '     exactly what the control shows it does anyway.\n']);
        end

        % ---- overall --------------------------------------------------------
        fprintf('\n--- verdict ---\n');
        if ~results.isPermutation
            fprintf(['  Jump mode is not misfiling data, it is losing it. Reordering is the\n' ...
                     '  wrong tool. Stop here.\n']);
        elseif ~results.deterministic
            fprintf(['  Jump mode collects the right values in an order that CHANGES every\n' ...
                     '  run. Nothing in software fixes that, and a flag output would only\n' ...
                     '  document the failure rather than prevent it. Stop here.\n']);
        elseif rBest > 0.9
            fprintf(['  Jump mode is a fixed lag of %+d in visit order. Correct the index and\n' ...
                     '  it works -- no flags, no settle time, no new hardware.\n'], results.bestLag);
        elseif outSample > 0.9
            fprintf(['  Jump mode is a fixed permutation that is not a simple lag, and it\n' ...
                     '  transfers between runs. Worth finding its structure; a flag output\n' ...
                     '  would identify it directly.\n']);
        else
            fprintf(['  Reproducible but not a lag and not recoverable by matching. This is\n' ...
                     '  the case where a flag output would actually earn its wiring.\n']);
        end

        if opts.plotIt; drawIt(results); end

    catch err
        a.useSequenceJump = savedJump;
        a.useSequenceWalk = savedWalk;
        a.useFgen         = savedFgen;
        rethrow(err);
    end
    a.useSequenceJump = false;
    a.useSequenceWalk = savedWalk;
    a.useFgen         = savedFgen;
    fprintf('\n(left useSequenceJump = false)\n');
end

% ======================================================================

function y = applyVisitLag(v, order, lag)
% v is indexed by FREQUENCY. order(k) is the frequency index visited kth. If the value
% filed under frequency order(k) actually belongs to frequency order(k+lag), rebuild the
% spectrum accordingly. Points whose partner falls off either end are left as NaN rather
% than wrapped, because wrapping would invent a data point.
    n = numel(v);
    y = nan(n, 1);
    for k = 1:n
        kk = k + lag;
        if kk >= 1 && kk <= n
            y(order(kk)) = v(order(k));
        end
    end
end

function perm = greedyAssign(v, target)
% Assign each target point the unused v point closest to it, taking the most
% distinctive targets first. A greedy stand-in for a full Hungarian assignment -- good
% enough here, because the point of test 4 is that even a GENEROUS matching fails
% out-of-sample when there is no real order to recover. A better matcher would only
% raise the in-sample number, which is already known to be meaningless.
    n = numel(v);
    perm = nan(n, 1);
    used = false(n, 1);
    % Most extreme targets first: they are the least ambiguous, so they get first pick.
    [~, tOrder] = sort(abs(target - median(target)), 'descend');
    for t = tOrder(:)'
        d = abs(v - target(t));
        d(used) = Inf;
        [~, j] = min(d);
        perm(t) = j;
        used(j) = true;
    end
end

function r = corrSafe(x, y)
    ok = isfinite(x) & isfinite(y);
    if sum(ok) < 3; r = NaN; return; end
    c = corrcoef(x(ok), y(ok));
    r = c(1, 2);
end

function rMax = corrCeiling(ref, err)
% What correlation two perfect repeats of this measurement would show, given the noise.
% Without it, "r = 0.97" cannot be read as good or bad.
    s = std(ref, 'omitnan');
    e = median(err, 'omitnan');
    if ~isfinite(s) || s <= 0 || ~isfinite(e); rMax = 1; return; end
    rMax = s^2 / (s^2 + e^2);
end

function v = ternary(c, a, b)
    if c; v = a; else; v = b; end
end

function drawIt(R)
    figure('Name', 'esrJumpUnscramble', 'Color', 'w');

    subplot(2, 2, 1); hold on; box on;
    plot(R.frequency, R.reference, '.-');
    plot(R.frequency, R.jump(:, 1), '.-');
    xlabel('MHz'); ylabel('normalized FL');
    title(sprintf('as measured (r = %+.3f)', corrSafe(R.jump(:,1), R.reference)));
    legend('reference', 'jump', 'Location', 'best');

    subplot(2, 2, 2); hold on; box on;
    plot(sort(R.reference), '.-');
    plot(sort(R.jump(:, 1)), '.-');
    xlabel('rank'); ylabel('normalized FL');
    title(sprintf('sorted: %.1fx noise apart', R.sortedRmsSig));
    legend('reference', 'jump', 'Location', 'best');

    subplot(2, 2, 3); hold on; box on;
    plot(R.lags, R.lagR, 'o-');
    yline(0, 'k:');
    xlabel('lag in visit order'); ylabel('r with reference');
    title(sprintf('best lag %+d at r = %+.3f', R.bestLag, R.bestLagR));

    subplot(2, 2, 4); hold on; box on;
    plot(R.frequency, R.reference, '.-');
    best = applyVisitLag(R.jump(:, 1), R.visitOrder(:, 1), R.bestLag);
    plot(R.frequency, best, '.-');
    xlabel('MHz'); ylabel('normalized FL');
    title(sprintf('jump un-lagged by %+d', R.bestLag));
    legend('reference', 'un-lagged', 'Location', 'best');
end

function o = withDefaults(o, d)
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(o, fn{i}) || isempty(o.(fn{i})); o.(fn{i}) = d.(fn{i}); end
    end
end
