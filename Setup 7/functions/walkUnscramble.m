function slice = walkUnscramble(data, nDet, nWarm, N, P, walkOrder)
%WALKUNSCRAMBLE Turn one walk pg.run's digitizer stream back into tau order.
%
%   slice = walkUnscramble(data, nDet, nWarm, N, P, walkOrder)
%
% This is the single place a sequence walk can corrupt data SILENTLY, which is why it is a
% function with a test rather than four lines inside performWalk.
%
% In walk mode the AWG sequencer advances one BLOCK per pulse-generator repeat, and the
% blocks were built in walkOrder (a randperm), not in tau order. So the digitizer stream
% arrives as
%
%       [ warm-up blocks ][ block 1 ][ block 2 ] ... [ block N ]   x P passes
%
% with block k holding the tau walkOrder(k), and each block contributing nDet detection
% windows. Getting the reordering wrong does not raise: it produces a Ramsey whose points
% are shuffled, which still fits, still decays, and returns a T2* that is simply wrong.
%
% INPUTS
%   data       the raw stream, at least nDet*(N+nWarm)*P long. Longer is fine; the head is
%              used, which is what the digitizer returns.
%   nDet       detection windows per pg repeat (2 for Rabi/ESR, 4 for Ramsey/Echo/T1 with
%              doubleMeasurement on)
%   nWarm      warm-up blocks at the START of each pass, discarded
%   N          number of tau points
%   P          passes in this pg.run
%   walkOrder  1xN permutation: walkOrder(k) is the tau index that block k played
%
% OUTPUT
%   slice      [nDet x N x P], indexed by TAU, ready to drop into the accumulator
%
% The inverse mapping is the subtle part. walkOrder(k) is "which tau block k played", so the
% assignment is  slice(:, walkOrder, :) = blocks  -- NOT slice = blocks(:, walkOrder, :),
% which applies the permutation the wrong way round and is wrong for every walkOrder except
% an involution. Both look equally plausible in review; only one survives the test.

    if ~isvector(walkOrder) || numel(walkOrder) ~= N
        error('walkUnscramble: walkOrder must have N = %d entries, got %d', N, numel(walkOrder));
    end
    if ~isequal(sort(walkOrder(:)'), 1:N)
        error('walkUnscramble: walkOrder must be a permutation of 1:%d', N);
    end
    blocksPerPass = N + nWarm;
    needed = nDet * blocksPerPass * P;
    if numel(data) < needed
        error(['walkUnscramble: the pass returned %d reads, expected at least %d ' ...
               '(%d windows x %d blocks x %d passes)'], numel(data), needed, nDet, ...
              blocksPerPass, P);
    end

    blocks = reshape(data(1:needed), nDet, blocksPerPass, P);
    blocks = blocks(:, nWarm+1:end, :);         % drop the warm-up
    slice = nan(nDet, N, P);
    slice(:, walkOrder, :) = blocks;
end
