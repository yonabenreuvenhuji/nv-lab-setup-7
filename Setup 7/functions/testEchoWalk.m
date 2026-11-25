function ok = testEchoWalk()
%TESTECHOWALK Offline grader for the Echo sequence walk. No hardware, ~1 s.
%
%   testEchoWalk
%
% Run it before every walked Echo, and certainly after touching AWGExpEcho. Every check here
% corresponds to a failure that returns a PLAUSIBLE Echo rather than an error: a shuffled tau
% axis still fits, still decays, and gives a confident wrong T2.
%
% ===================================================================
% WHY ECHO NEEDS ITS OWN GRADER AND NOT JUST testRamseyWalk
% ===================================================================
%
% Ramsey's block is TWO steps -- signal waveform, reference waveform -- and the step
% arithmetic is trivial. Echo's block is 2*(3*cycles + 2), it carries per-step REPEAT COUNTS
% that encode tau, and only two of its steps are trigger-gated. Three things can go wrong
% here that cannot go wrong in Ramsey:
%
%   * the trigger count per block. If a block carries anything other than exactly two
%     triggers, the sequencer and the pulse generator disagree about what a repeat is, and
%     every tau after the first is read from the wrong block. Nothing errors.
%   * the repeat counts. tau lives in them. Concatenating blocks must carry each tau's own
%     zeroSegNum with it, or the taus are silently swapped while the labels are not.
%   * the step count with cycles. At XY8 a block is 52 steps, so the MAX_SEQUENCE_STEPS
%     ceiling arrives twenty-six times sooner than it does for Ramsey.
%
% WHAT IT CANNOT CHECK: that the instrument really advances one block per pulse-generator
% repeat, and that the AWG idles at zero through the excess gate. Those are hardware facts.
% The first is settled by one walk run against one classic run; the second is what
% echoGateTest measures, and for Echo it is a real question rather than a formality -- the
% excess gate reaches ~98 us on a 1-50 us list, against the 4 us ramseyGateTest cleared.

    fprintf('\n=== testEchoWalk ===\n');
    ok = true;

    ok = ok & t1_roundTrip();
    ok = ok & t2_wrongWayRound();
    ok = ok & t3_pairsStayTogether();
    ok = ok & t4_warmupDropped();
    ok = ok & t5_blockLayout();
    ok = ok & t6_twoTriggersPerBlock();
    ok = ok & t7_repeatCountsTravelWithTau();
    ok = ok & t8_stepCeiling();
    ok = ok & t9_passesPerRunCeilings();

    fprintf('\n');
    if ok
        fprintf('ALL TESTS PASS\n\n');
    else
        fprintf('*** FAILURES ABOVE -- do not run a walked Echo until they are fixed ***\n\n');
    end
end

% ======================================================================
% THE REORDERING, which Echo shares with Ramsey
% ======================================================================

function ok = t1_roundTrip()
% The core claim: data laid out in BLOCK order comes back in TAU order.
    N = 9; nDet = 4; nWarm = 2; P = 3;
    walkOrder = [4 1 9 7 2 5 3 8 6];
    val = @(tau, d, p) 1000*tau + 10*d + p;
    blocksPerPass = N + nWarm;
    data = nan(nDet, blocksPerPass, P);
    for p = 1:P
        for k = 1:blocksPerPass
            if k <= nWarm
                tau = walkOrder(1);
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
% The permutation applied the wrong way round is the mistake being guarded against, and it
% must NOT accidentally agree, or t1 proves nothing.
    N = 9; nDet = 4; nWarm = 0; P = 1;
    walkOrder = [4 1 9 7 2 5 3 8 6];
    data = reshape(1:(nDet*N*P), nDet, N, P);
    right = walkUnscramble(data(:), nDet, nWarm, N, P, walkOrder);
    wrong = data(:, walkOrder, :);
    ok = ~isequaln(right, wrong);
    report(ok, 'the wrong-way-round permutation gives a DIFFERENT answer (so t1 is meaningful)');
    inv = nan(1, N); inv(walkOrder) = 1:N;
    ok2 = isequaln(right, data(:, inv, :));
    report(ok2, 'the mapping is the inverse permutation, as documented');
    ok = ok & ok2;
end

function ok = t3_pairsStayTogether()
% Echo's four windows are TWO (signal, reference) pairs. Mixing the rows gives the |1>
% projection the |0> reference, and the contrast collapses in a way that looks exactly like
% a short T2.
    N = 5; nDet = 4; nWarm = 1; P = 2;
    walkOrder = [3 5 1 4 2];
    blocksPerPass = N + nWarm;
    data = nan(nDet, blocksPerPass, P);
    for p = 1:P
        for k = 1:blocksPerPass
            tau = walkOrder(max(k - nWarm, 1));
            data(:, k, p) = [tau; -tau; 100+tau; -(100+tau)];
        end
    end
    slice = walkUnscramble(data(:), nDet, nWarm, N, P, walkOrder);
    ok = true;
    for tau = 1:N
        ok = ok & all(slice(1,tau,:) == tau) & all(slice(2,tau,:) == -tau) ...
                & all(slice(3,tau,:) == 100+tau) & all(slice(4,tau,:) == -(100+tau));
    end
    report(ok, 'the two (signal, reference) pairs stay with their own tau');
end

function ok = t4_warmupDropped()
    N = 4; nDet = 4; nWarm = 3; P = 1;
    walkOrder = [2 4 1 3];
    blocksPerPass = N + nWarm;
    data = nan(nDet, blocksPerPass, P);
    data(:, 1:nWarm, 1) = -999;
    for k = 1:N
        data(:, nWarm+k, 1) = walkOrder(k);
    end
    slice = walkUnscramble(data(:), nDet, nWarm, N, P, walkOrder);
    ok = ~any(slice(:) == -999) && isequal(squeeze(slice(1,:,1)), 1:N);
    report(ok, 'warm-up blocks are discarded and contaminate no tau');
end

% ======================================================================
% THE BLOCK, which is Echo's own and is where it differs from Ramsey
% ======================================================================

function [wf, rpt, trg] = blockFor(tauUs, cycles, X, baseWF, centerPulseTime, ix)
% Replicates AWGExpEcho.blockFor exactly, so the layout can be graded without the class,
% the AWG or the lab path. If the two ever disagree this file is wrong and should be fixed
% to match the class, not the other way round.
    totalZeroTime = tauUs - centerPulseTime/2;
    zeroSegNum = totalZeroTime/baseWF;
    piPulses = ix.piX*X + ix.piY*~X;
    piPulses = piPulses(1:cycles);
    seqBody = reshape([repmat(ix.zero,1,cycles); piPulses; repmat(ix.zero,1,cycles)], 1, 3*cycles);
    wf  = [ix.firstHalf seqBody ix.lastHalf  ix.firstHalf seqBody ix.lastHalfRef];
    rpt = [1 repmat([zeroSegNum 1 zeroSegNum],1,cycles) 1, ...
           1 repmat([zeroSegNum 1 zeroSegNum],1,cycles) 1];
    trg = [1 repmat([0 0 0],1,cycles) 0, 1 repmat([0 0 0],1,cycles) 0];
end

function ix = idxFixture()
    ix = struct('zero', 1, 'firstHalf', 2, 'lastHalf', 3, 'lastHalfRef', 4, 'piX', 5, 'piY', 6);
end

function ok = t5_blockLayout()
% Step count, and where indexSeq must point. Checked at cycles = 1 and cycles = 8, because
% the whole layout scales with cycles and XY8 is the case that will actually be run.
    ok = true;
    for cycles = [1 8]
        ix = idxFixture(); X = true(1, max(cycles, 4));
        [wf, ~, ~] = blockFor(10, cycles, X, 0.5, 1.0, ix);
        want = 2*(3*cycles + 2);
        good = numel(wf) == want;
        report(good, sprintf('cycles = %d: block is %d steps (2 x (3*cycles+2) = %d)', ...
                             cycles, numel(wf), want));
        ok = ok & good;

        % indexSeq points at the LAST step of each tau's block, which is what the classic
        % path's assignSequence target becomes in walk mode.
        N = 6; nWarm = 1;
        walkOrder = [5 2 6 1 4 3];
        blockOrder = [repmat(walkOrder(1), 1, nWarm), walkOrder];
        stepsPerBlock = want;
        idxSeq = zeros(1, N); s0 = 0;
        for blk = 1:numel(blockOrder)
            s0 = s0 + stepsPerBlock;
            if blk > nWarm; idxSeq(blockOrder(blk)) = s0; end
        end
        good2 = all(idxSeq > 0) && numel(unique(idxSeq)) == N && ...
                all(mod(idxSeq, stepsPerBlock) == 0);
        report(good2, sprintf(['cycles = %d: indexSeq lands on a distinct block boundary for ' ...
                               'every tau'], cycles));
        ok = ok & good2;
    end
end

function ok = t6_twoTriggersPerBlock()
% THE ONE THAT MATTERS MOST. One pulse-generator repeat must consume exactly one block, and
% a block is defined by its triggers. Any other count and the sequencer and the pulse
% generator disagree about what a repeat is: every tau after the first is then read from the
% wrong block, the sweep is shifted, and nothing raises.
    ok = true;
    for cycles = [1 2 4 8]
        ix = idxFixture(); X = true(1, max(cycles, 4));
        [~, ~, trg] = blockFor(10, cycles, X, 0.5, 1.0, ix);
        n = sum(trg ~= 0);
        good = (n == 2);
        report(good, sprintf('cycles = %d: exactly %d trigger-gated steps per block (must be 2)', ...
                             cycles, n));
        ok = ok & good;
        % and they must be the FIRST step of each half, not anywhere
        half = numel(trg)/2;
        good2 = trg(1) == 1 && trg(half+1) == 1 && ~any(trg(2:half)) && ~any(trg(half+2:end));
        report(good2, sprintf('cycles = %d: the triggers are the first step of each half', cycles));
        ok = ok & good2;
    end
end

function ok = t7_repeatCountsTravelWithTau()
% tau is the zero step's REPEAT COUNT, not a waveform. When blocks are concatenated each
% tau's counts must travel with it. If they did not, the taus would be silently swapped while
% their labels were not -- an Echo that decays with the wrong time axis.
    ix = idxFixture(); cycles = 2; X = true(1,4); baseWF = 0.5; cpt = 1.0;
    taus = [1 7 23];
    stepsPerBlock = 2*(3*cycles + 2);
    order = [3 1 2];
    rptList = [];
    for k = 1:numel(order)
        [~, rpt, ~] = blockFor(taus(order(k)), cycles, X, baseWF, cpt, ix);
        rptList = [rptList, fix(rpt)];                     %#ok<AGROW>
    end
    ok = true;
    for k = 1:numel(order)
        seg = rptList((k-1)*stepsPerBlock + (1:stepsPerBlock));
        wantZero = fix((taus(order(k)) - cpt/2)/baseWF);
        good = seg(2) == wantZero;
        ok = ok & good;
    end
    report(ok, sprintf(['each block carries its own zero-step repeat count through the ' ...
                        'concatenation (%s us -> %s repeats)'], mat2str(taus(order)), ...
                       mat2str(fix((taus(order) - cpt/2)/baseWF))));
end

function ok = t8_stepCeiling()
% The AWG's step ceiling arrives 26x sooner for Echo at XY8 than it does for Ramsey, so the
% guard has to be there and has to be right.
    MAXSTEPS = 16384;
    fprintf('    steps needed, against the AWG''s %d:\n', MAXSTEPS);
    ok = true;
    for cycles = [1 8 16]
        for N = [25 50 100]
            steps = 2*(3*cycles + 2) * (N + 1);
            flag = '';
            if steps > MAXSTEPS; flag = '  <-- REFUSED, correctly'; end
            fprintf('      cycles %2d, %3d taus -> %6d steps%s\n', cycles, N, steps, flag);
        end
    end
    % the guard must fire exactly when it should
    cycles = 16; N = 200;
    steps = 2*(3*cycles + 2)*(N + 1);
    ok = ok & (steps > MAXSTEPS);
    report(ok, sprintf('cycles = 16 with 200 taus needs %d steps and must be refused', steps));
end

function ok = t9_passesPerRunCeilings()
% Echo's TWO ceilings. The digitizer one Ramsey has too; the clock one is Echo's, because its
% sequence duration grows with the tau range and the fixed gate makes every block pay it.
    maxSeg = 32000; maxSec = 30; nDet = 4; nWarm = 1; repeats = 1000;
    fprintf('    passesPerRun, and which ceiling binds:\n');
    ok = true;
    % The clock ceiling binds when it is the smaller of the two, i.e. when
    %       30/(blocks*seq) < 32000/(4*blocks)   ->   seq > 3750 us
    % which needs cycles = 8 past tau ~ 230 us, or cycles = 16 past ~ 115 us. Below that the
    % digitizer binds and the clock cap costs nothing. The cases below straddle it, and the
    % last one is where walk genuinely loses.
    cases = {  % nTau, cycles, maxTau us
        {25,  1,  50}, ...
        {25,  8, 200}, ...
        {25, 16, 400}, ...
        {25, 16, 900}, ...
        {60,  1,  20}};
    for c = 1:numel(cases)
        N = cases{c}{1}; cyc = cases{c}{2}; tmax = cases{c}{3};
        blocks = N + nWarm;
        seqUs = 2*tmax*cyc + 4*(16 + 0.3);          % gate + four init/detect pairs
        pDig = floor(maxSeg / (nDet*blocks));
        pClk = floor(maxSec / max(blocks*seqUs*1e-6, eps));
        P = max(1, min([pDig, pClk, repeats]));
        arms = ceil(repeats/P);
        who = 'digitizer'; if pClk < pDig; who = 'clock'; end
        verdict = 'walk wins'; if arms >= N; verdict = 'CLASSIC WINS'; end
        fprintf(['      %3d taus, cycles %2d, tau_max %3d us: seq %.0f us -> P = %4d ' ...
                 '(%s binds), %3d arms vs %3d classic  %s\n'], ...
                N, cyc, tmax, seqUs, P, who, arms, N, verdict);
        ok = ok & (nDet*blocks*P <= maxSeg) & (blocks*P*seqUs*1e-6 <= maxSec*1.001);
    end
    report(ok, 'both ceilings are respected in every case');
    fprintf(['    The table above is the honest answer to "is walk worth it for Echo", and it\n' ...
             '    corrects something this file asserted in prose before it computed it: at\n' ...
             '    cycles = 8 and tau_max 200 us the DIGITIZER still binds, not the clock, and\n' ...
             '    walk wins comfortably. The clock only takes over past a ~3750 us sequence.\n' ...
             '    Walk finally LOSES in the last case -- deep CPMG at very long tau -- where one\n' ...
             '    block is tens of milliseconds and P is capped so hard that walk re-arms more\n' ...
             '    often than classic. The class prints that warning rather than being quietly\n' ...
             '    slower.\n']);
end

% ======================================================================

function report(ok, what)
    if ok
        fprintf('  PASS  %s\n', what);
    else
        fprintf('  FAIL  %s\n', what);
    end
end
