function [value, sterr] = ratioAcrossAverages(S1, S2, S1sterr, S2sterr, nParams)
%RATIOACROSSAVERAGES Normalized signal S1/S2, combined across averages.
%
%   [value, sterr] = ratioAcrossAverages(S1, S2, S1sterr, S2sterr, nParams)
%
% Replaces Experiment.getRatioDistributionValues for experiments whose
% observable is the RATIO of two measured channels (ESR, Rabi, ... : signal
% with MW over reference without MW).
%
% WHY: the original combined the two channels across averages separately,
% via getCombinedSterr, and only formed the ratio at the end. getCombinedSterr
% is a pooled-variance formula, so it charges the difference in MEAN LEVEL
% between averages as though it were statistical scatter:
%
%       + n*m*(meanSoFar - means(:,i)).^2 ./ ((n+m)^2*(n+m-1))
%
% But the whole point of dividing by the reference is that the overall
% fluorescence level is ALLOWED to drift (NV brightness, laser power, focus).
% So the error bar grew with drift even though the ratio was rock steady --
% up to ~10x inflation at 20% drift, while the plotted value stayed correct.
% Worse on the photodiode/digitizer path, where processData deliberately sets
% the reference column to a constant (so its sterr is exactly 0, all the noise
% having been folded into channel 1). getCombinedSterr then receives all-zero
% errors and returns a number driven ENTIRELY by drift, with no statistical
% content at all.
%
% The fix is to form the ratio FIRST, per average -- which is where the level
% drift cancels -- and only then combine across averages.
%
% INPUTS
%   S1, S2            - signal and reference. Either nParams x nAverages, or a
%                       1 x nParams row (see the shape note below).
%   S1sterr, S2sterr  - matching standard errors.
%   nParams           - number of swept parameter points, i.e.
%                       obj.getTotalNumberOfParams(). Required, because the
%                       incoming shape is ambiguous without it.
%
% OUTPUTS
%   value, sterr      - nParams x 1 columns.
%
% SHAPE NOTE -- this is the trap that makes nParams mandatory. perform() passes
% squeeze(obj.signal(k, :, 1:currIter)), and squeeze does NOT touch a 2-D
% array, so the orientation FLIPS with the number of averages:
%
%   currIter = 1 : signal(1,:,1)   is 1 x nParams        -> squeeze -> 1 x nParams   (ROW)
%   currIter = 2 : signal(1,:,1:2) is 1 x nParams x 2    -> squeeze -> nParams x 2   (MATRIX)
%
% So you cannot infer "how many averages" from size(...,2): at one average
% that is nParams, not 1. Reshaping against nParams pins it down for both.

    % --- pin the orientation down, whatever we were handed -----------------
    S1      = reshape(S1,      nParams, []);
    S2      = reshape(S2,      nParams, []);
    S1sterr = reshape(S1sterr, nParams, []);
    S2sterr = reshape(S2sterr, nParams, []);
    nAvg    = size(S1, 2);

    % --- per-average ratio and its error ----------------------------------
    % Standard propagation for a quotient. On the photodiode path S2sterr is
    % identically zero and this reduces to ste(refSignal), which is exactly
    % the error of the normalized signal -- the original single-average
    % behaviour, preserved.
    r  = S1 ./ S2;
    rs = abs(r) .* sqrt((S1sterr./S1).^2 + (S2sterr./S2).^2);

    if nAvg == 1
        value = r;
        sterr = rs;
    else
        % --- combine the RATIOS across averages ---------------------------
        % Inverse-variance weighted mean. Level drift has already cancelled
        % in r, so it cannot leak into the error bar here.
        usableElem = isfinite(r) & isfinite(rs) & (rs > 0);
        w = zeros(size(rs));
        w(usableElem) = 1 ./ rs(usableElem).^2;
        sumW = sum(w, 2);
        haveW = sumW > 0;

        value     = zeros(nParams, 1);
        sterrStat = zeros(nParams, 1);

        value(haveW)     = sum(r(haveW,:) .* w(haveW,:), 2) ./ sumW(haveW);
        sterrStat(haveW) = 1 ./ sqrt(sumW(haveW));

        if any(~haveW)
            % No usable error anywhere on this row (e.g. every sterr came back
            % as exactly 0). Fall back to an unweighted mean of the finite
            % entries and let the scatter term below carry the uncertainty.
            for k = find(~haveW).'
                finiteR = r(k, isfinite(r(k,:)));
                if isempty(finiteR)
                    value(k) = NaN;
                else
                    value(k) = mean(finiteR);
                end
            end
        end

        % --- stay honest about genuine disagreement ------------------------
        % Birge ratio (the PDG scale factor): if the per-average ratios scatter
        % by more than their own error bars predict, chi2/dof exceeds 1 and the
        % error is scaled by sqrt(chi2/dof). Reduces to the plain
        % inverse-variance error when the averages agree, so it costs nothing
        % in the normal case.
        %
        % Crucially, level drift has already cancelled in r, so drift does NOT
        % drive chi2. Only a real shift in the NORMALIZED signal does -- a
        % resonance moving because the field drifted, say. That is exactly the
        % discrimination the old code failed to make.
        %
        % Needs dof >= 2 to be worth anything: with two averages there is one
        % degree of freedom and chi2/dof fluctuates so wildly (it is a chi2_1
        % variate, so values of 4+ are common by chance) that scaling by it
        % would inject noise into the error bar rather than information -- and
        % would reproduce the original complaint in miniature, with the bar
        % getting BIGGER on the second average. So below three averages the
        % statistical error is reported as-is.
        sterr = sterrStat;
        if nAvg >= 3
            resid = r - value;              % implicit expansion, nParams x nAvg
            resid(~usableElem) = 0;
            chi2  = sum(w .* resid.^2, 2);
            dof   = max(sum(usableElem, 2) - 1, 1);
            sterr = sterrStat .* sqrt(max(1, chi2 ./ dof));
        end
    end

    % --- shape guard -------------------------------------------------------
    % Fail loudly rather than handing a wrong-sized array to ExpParameter,
    % whose error message ("Trying to define value for ... failed!") does not
    % say what the size actually was.
    if numel(value) ~= nParams || numel(sterr) ~= nParams
        error('ratioAcrossAverages: produced %d values and %d errors for %d parameter points', ...
              numel(value), numel(sterr), nParams);
    end
    value = reshape(value, nParams, 1);
    sterr = reshape(sterr, nParams, 1);
end
