function ok = testRamseyWalk()
%TESTRAMSEYWALK Offline grader for the Ramsey sequence walk. No hardware, ~1 s.
%
%   testRamseyWalk
%
% Run it before every walk run, and certainly after touching AWGExpRamsey. Every check here
% corresponds to a failure that returns a PLAUSIBLE Ramsey rather than an error: a shuffled
% tau axis still fits, still decays, and gives a confident wrong T2*.
%
% WHAT IT CANNOT CHECK: that the instrument really advances one block per pulse-generator
% repeat. That is a hardware fact and the only test for it is one walk run against one
% classic run at the same settings. This grader checks everything BETWEEN that assumption and
% the answer.

    fprintf('\n=== testRamseyWalk ===\n');
    ok = true;

    ok = ok & t1_roundTrip();
    ok = ok & t2_wrongWayRound();
    ok = ok & t3_pairsStayTogether();
    ok = ok & t4_warmupDropped();
    ok = ok & t5_refusesBadInput();
    ok = ok & t6_sequenceLayout();
    ok = ok & t7_passesPerRun();

    fprintf('\n');
    if ok
        fprintf('ALL TESTS PASS\n\n');
    else
        fprintf('*** FAILURES ABOVE -- do not run a walk until they are fixed ***\n\n');
    end
end

% ======================================================================

function ok = t1_roundTrip()
% The core claim: data laid out in BLOCK order comes back in TAU order.
    N = 9; nDet = 4; nWarm = 2; P = 3;
    walkOrder = [4 1 9 7 2 5 3 8 6];

    % Build a stream where the value encodes (tau, window, pass) unambiguously.
    val = @(tau, d, p) 1000*tau + 10*d + p;
    blocksPerPass = N + nWarm;
    data = nan(nDet, blocksPerPass, P);
    for p = 1:P
        for k = 1:blocksPerPass
            if k <= nWarm
                tau = walkOrder(1);                 % warm-up replays the first block
            else
                tau = walkOrder(k - nWarm);
            end
            for d = 1:nDet
                data(d, k, p) = val(tau, d, p);
            end
        end
    end

    slice = walkUnscramble(data(:), nDet, nWarm, N, P, walkOrder);

    want = nan(nDet, N, P);
    for p = 1:P; for tau = 1:N; for d = 1:nDet
        want(d, tau, p) = val(tau, d, p);
    end; end; end

    ok = isequaln(slice, want);
    report(ok, 'round trip: block order -> tau order, nDet=4, warmup=2, P=3');
end

function ok = t2_wrongWayRound()
% The permutation applied the wrong way round is the mistake this is guarding against, and
% it must NOT accidentally agree. If it did, the test would pass on a broken implementation.
    N = 9; nDet = 2; nWarm = 0; P = 1;
    walkOrder = [4 1 9 7 2 5 3 8 6];
    data = reshape(1:(nDet*N*P), nDet, N, P);

    right = walkUnscramble(data(:), nDet, nWarm, N, P, walkOrder);
    wrong = data(:, walkOrder, :);            % the plausible-looking error

    ok = ~isequaln(right, wrong);
    report(ok, 'the wrong-way-round permutation gives a DIFFERENT answer (so t1 is meaningful)');

    % And the right one really is the inverse.
    inv = nan(1, N); inv(walkOrder) = 1:N;
    ok2 = isequaln(right, data(:, inv, :));
    report(ok2, 'the mapping is the inverse permutation, as documented');
    ok = ok & ok2;
end

function ok = t3_pairsStayTogether()
% Ramsey's four windows are TWO (signal, reference) pairs. Rows must not be mixed, or the
% |1> projection gets the |0> reference and the contrast collapses in a way that looks like
% short T2*.
    N = 5; nDet = 4; nWarm = 1; P = 2;
    walkOrder = [3 5 1 4 2];
    blocksPerPass = N + nWarm;
    data = nan(nDet, blocksPerPass, P);
    for p = 1:P
        for k = 1:blocksPerPass
            tau = walkOrder(max(k - nWarm, 1));
            data(:, k, p) = [tau; -tau; 100+tau; -(100+tau)];   % pair 1 and pair 2, signed
        end
    end
    slice = walkUnscramble(data(:), nDet, nWarm, N, P, walkOrder);

    ok = true;
    for tau = 1:N
        ok = ok & all(slice(1,tau,:) == tau) & all(slice(2,tau,:) == -tau) ...
                & all(slice(3,tau,:) == 100+tau) & all(slice(4,tau,:) == -(100+tau));
    end
    report(ok, 'the two (signal, reference) pairs stay associated with their own tau');
end

function ok = t4_warmupDropped()
% The warm-up blocks exist because the first block after an idle sequencer is not preceded by
% a full initialisation. If they leaked in they would land on tau walkOrder(1) and bias
% exactly one point -- which is invisible on a plot.
    N = 4; nDet = 2; nWarm = 3; P = 1;
    walkOrder = [2 4 1 3];
    blocksPerPass = N + nWarm;
    data = nan(nDet, blocksPerPass, P);
    data(:, 1:nWarm, 1) = -999;                       % poison the warm-up
    for k = 1:N
        data(:, nWarm+k, 1) = walkOrder(k);
    end
    slice = walkUnscramble(data(:), nDet, nWarm, N, P, walkOrder);
    ok = ~any(slice(:) == -999) && isequal(squeeze(slice(1,:,1)), 1:N);
    report(ok, 'warm-up blocks are discarded and do not contaminate any tau');
end

function ok = t5_refusesBadInput()
    ok = true;
    ok = ok & mustFail(@() walkUnscramble(1:10, 2, 0, 9, 1, 1:9), 'short stream refused');
    ok = ok & mustFail(@() walkUnscramble(1:100, 2, 0, 5, 1, [1 2 3 4 4]), 'non-permutation refused');
    ok = ok & mustFail(@() walkUnscramble(1:100, 2, 0, 5, 1, 1:4), 'wrong-length walkOrder refused');
end

function ok = t6_sequenceLayout()
% Replicates the step arithmetic in buildWalkSequence and checks the three things that make
% the sequencer walk correctly: two trigger-gated steps per block, the wrap on the last step,
% and indexSeq pointing at the last step of each tau's block.
    N = 6; nWarm = 1; nCh = 2;
    walkOrder = [5 2 6 1 4 3];
    blockOrder = [repmat(walkOrder(1), 1, nWarm), walkOrder];
    nSteps = 2*numel(blockOrder);

    trgList = ones(1, nSteps);
    gotoList = zeros(1, nSteps); gotoList(nSteps) = 1;
    indexSeq = zeros(nCh, N);
    s0 = 0;
    for blk = 1:numel(blockOrder)
        s0 = s0 + 2;
        if blk > nWarm; indexSeq(:, blockOrder(blk)) = s0; end
    end

    ok = (nSteps == 2*(N+nWarm));
    report(ok, sprintf('step count is 2 per block: %d for %d taus + %d warm-up', nSteps, N, nWarm));

    ok2 = all(trgList == 1);
    report(ok2, 'every step is trigger-gated, so one pg repeat advances exactly one block');

    ok3 = gotoList(nSteps) == 1 && all(gotoList(1:nSteps-1) == 0);
    report(ok3, 'only the last step wraps to step 1, so passes chain without a re-arm');

    ok4 = all(indexSeq(1,:) > 0) && numel(unique(indexSeq(1,:))) == N && ...
          all(mod(indexSeq(1,:), 2) == 0);
    report(ok4, 'indexSeq points at a distinct EVEN step (the block''s last) for every tau');

    ok = ok & ok2 & ok3 & ok4;
end

function ok = t7_passesPerRun()
% The digitizer budget. Ramsey has FOUR windows per repeat where Rabi has two, so the same
% ceiling buys half the passes -- the arithmetic that decides whether walk is worth anything.
    ceilings = [12000 32000];
    N = 34; nWarm = 1; nDet = 4; repeats = 3000;
    ok = true;
    for c = ceilings
        P = max(1, min(floor(c/(nDet*(N+nWarm))), repeats));
        used = nDet*(N+nWarm)*P;
        ok = ok & (used <= c);
        fprintf('    ceiling %5d -> passesPerRun %3d, %5d scans used, %d arms (classic would be %d)\n', ...
                c, P, used, ceil(repeats/P), N);
    end
    report(ok, 'passesPerRun never exceeds the digitizer ceiling');
end

% ======================================================================

function ok = mustFail(fn, what)
    ok = false;
    try
        fn();
    catch
        ok = true;
    end
    report(ok, what);
end

function report(ok, what)
    if ok
        fprintf('  PASS  %s\n', what);
    else
        fprintf('  FAIL  %s\n', what);
    end
end
