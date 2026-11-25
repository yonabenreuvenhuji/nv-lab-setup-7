function m = walkCompareMetrics(ref, refE, walk, walkE)
%WALKCOMPAREMETRICS Compare two spectra, accounting for their own noise.
%
%   m = walkCompareMetrics(ref, refE, walk, walkE)
%
% Decides whether two measurements of the same spectrum agree. Written because a
% flat "correlation > 0.9" test gave false failures: correlation between two NOISY
% measurements of the same thing is bounded well below 1, and the bound depends on
% the SNR of the sweep. At repeats = 20 a correct pair legitimately correlates at
% only ~0.83, which is exactly what walk mode produced while its graph visibly
% matched the reference.
%
% For ref = s + n1 and walk = s + n2 with independent noise,
%
%     var(ref)      = var(s) + e1^2
%     var(walk)     = var(s) + e2^2
%     cov(ref,walk) = var(s)
%     =>  r = var(s) / sqrt(var(ref) * var(walk))
%
% so the achievable correlation follows from the error bars. var(s) is estimated by
% subtracting the noise variance from each measured variance and averaging. The
% test then asks whether the observed r reaches that ceiling, not whether it
% reaches 1.
%
% FIELDS
%   correlation   Pearson r at zero shift
%   rMax          the noise-limited ceiling described above
%   rRatio        correlation / rMax. ~1 means "as good as the noise allows"
%   bestShift     circular shift APPLIED TO walk that maximises r, reported as
%                 non-zero ONLY when it beats zero shift by more than the sampling
%                 noise on r (see shiftMargin below). Non-zero means the sweep is
%                 gaining or losing triggers -- an exact, fixable off-by-N. Sign
%                 convention: a walk lagging ref by one point is realigned by
%                 shifting it -1, so it reports bestShift = -1
%   bestShiftCorr r at that shift
%   rawBestShift, rawBestShiftCorr  the bare argmax, before the significance test
%   shiftMargin   how much a shift had to win by to count
%   sigmaDiff     median |ref - walk| / combined sigma. ~0.8 for a correct pair
%   scatterRatio  point-to-point roughness over what the error bars predict.
%                 ~1.4 for a real spectrum; >>1 means scrambled, <<1 means stuck
%   firstPointSigma  |ref(1) - walk(1)| / sigma. Called out separately because the
%                 walk path discards the first PASS, whereas the classic path
%                 discards the first REPEAT AT EVERY FREQUENCY -- so frequency 1
%                 is the one point that could carry an un-discarded warm-up
%   pass          all of: bestShift == 0, rRatio > 0.7, sigmaDiff < 2

    ref = ref(:); refE = refE(:); walk = walk(:); walkE = walkE(:);
    if numel(ref) ~= numel(walk)
        error('walkCompareMetrics: %d reference values vs %d walk values', numel(ref), numel(walk));
    end
    n = numel(ref);

    m.correlation = corr2col(walk, ref);

    % --- noise-limited ceiling --------------------------------------------
    e1sq = mean(refE.^2);
    e2sq = mean(walkE.^2);
    vRef  = var(ref);
    vWalk = var(walk);
    vSig  = mean([max(vRef - e1sq, 0), max(vWalk - e2sq, 0)]);
    if vRef > 0 && vWalk > 0
        m.rMax = min(1, vSig / sqrt(vRef * vWalk));
    else
        m.rMax = 0;
    end
    if m.rMax > 0
        m.rRatio = m.correlation / m.rMax;
    else
        m.rRatio = 0;
    end

    % --- circular shift: catches an off-by-N trigger count -----------------
    % A plain argmax over 11 shifts is NOT a shift test. On a 96-point sweep it
    % reported bestShift = -1 at r = +0.7559 against r = +0.7556 at zero shift --
    % a difference of 0.0003, which is nothing, and it turned a run whose graph
    % matched the reference into a FAIL. With 11 candidates one of them beats zero
    % shift by a hair essentially always.
    %
    % So a shift now has to WIN by more than the sampling noise on r. The standard
    % error of a correlation is about (1-r^2)/sqrt(n), and two shifted versions of
    % the same data are highly dependent, so 2 of those is a generous bar; 0.05 is
    % a floor for short sweeps. A real off-by-one on a resolved resonance clears
    % this easily -- it moves r by tenths, not thousandths.
    maxShift = min(5, n - 1);
    shifts = -maxShift:maxShift;
    rs = zeros(size(shifts));
    for i = 1:length(shifts)
        rs(i) = corr2col(circshift(walk, shifts(i)), ref);
    end
    [rBest, iBest] = max(rs);
    m.rawBestShift     = shifts(iBest);        % argmax, kept for diagnostics
    m.rawBestShiftCorr = rBest;
    m.shiftMargin = max(0.05, 2 * (1 - m.correlation^2) / sqrt(max(n, 4)));
    if m.rawBestShift ~= 0 && rBest > m.correlation + m.shiftMargin
        m.bestShift     = m.rawBestShift;
        m.bestShiftCorr = rBest;
    else
        m.bestShift     = 0;
        m.bestShiftCorr = m.correlation;
    end

    % --- agreement within the error bars ----------------------------------
    sigma = sqrt(refE.^2 + walkE.^2);
    sigma(sigma <= 0) = eps;
    m.sigmaDiff = median(abs(ref - walk) ./ sigma);

    % --- shape sanity ------------------------------------------------------
    if n > 1
        m.scatterRatio = median(abs(diff(walk))) / max(1.4*median(walkE), eps);
    else
        m.scatterRatio = NaN;
    end

    m.firstPointSigma = abs(ref(1) - walk(1)) / sigma(1);

    m.pass = (m.bestShift == 0) && (m.rRatio > 0.7) && (m.sigmaDiff < 2);
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
