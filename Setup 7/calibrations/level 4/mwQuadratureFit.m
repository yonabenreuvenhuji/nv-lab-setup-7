function f = mwQuadratureFit(cfg, R, E)
%MWQUADRATUREFIT Where the both-on quadrature signal crosses zero, in dB of channel 2.
%
%   f = mwQuadratureFit(cfg, R)        R: nConfig x nAverages, the sig/ref ratio per average
%   f = mwQuadratureFit(cfg, R, E)     E: nConfig x 1 errors, used when there are < 3 averages
%
% cfg is AWGExpMwQuadrature's configuration list (fields kind, offsetDb, lead, refPhaseDeg).
% No hardware; calibMwQuadrature calls it, and testMwQuadrature grades it.
%
% For every channel-2 offset d, per average:
%
%   Q(d) = [A(135) - A(-45)] - [B(135) - B(-45)]     A: ch2 leads by 90, B: ch1 leads
%   P(d) = [A(135) - A(-45)] + [B(135) - B(-45)]     the phase-sensitive combination
%
% Q is first order in the imbalance and cancels a channelPhase error to first order; P is
% the reverse (it carries what both orders share: detuning and any channelPhase error), and
% is reported, not used. Q is formed WITHIN each average, so slow drift between averages
% cancels before the averages are combined.
%
% A weighted straight line Q = a + b*d through the offsets gives the zero crossing
% d* = -a/b, with its error from the full covariance (Birge-scaled when chi2/dof > 1).
% The balance with both channels on is then delta + d*. The sign of b is whatever the
% phase convention makes it -- it is measured, not assumed, which is the point of
% playing several offsets.
%
% Contrast = bright - dark (no pulse minus an in-phase pi), to put Q in units of the full
% spin contrast. The zero crossing does not depend on it.

    if ~exist('E', 'var'); E = []; end
    kinds = {cfg.kind};
    quad = strcmp(kinds, 'quad');
    offs = unique([cfg(quad).offsetDb]);
    nA = size(R, 2);
    f = struct('ok', false, 'why', '', 'offsets', offs, 'nAverages', nA);

    iBright = find(strcmp(kinds, 'bright'), 1);
    iDark   = find(strcmp(kinds, 'dark'), 1);
    C = nan(1, nA);
    if ~isempty(iBright) && ~isempty(iDark)
        C = R(iBright, :) - R(iDark, :);
    end
    f.contrast = nmean(C);
    f.contrastErr = nsem(C);

    nO = numel(offs);
    Qm = nan(1, nO);  Qe = nan(1, nO);  Pm = nan(1, nO);  Pe = nan(1, nO);
    nUsed = zeros(1, nO);
    for k = 1:nO
        d = offs(k);
        iAp = pick(cfg, d, 2, 135);  iAm = pick(cfg, d, 2, -45);
        iBp = pick(cfg, d, 1, 135);  iBm = pick(cfg, d, 1, -45);
        if any([iAp iAm iBp iBm] == 0)
            f.why = sprintf('configuration list incomplete at offset %+g dB', d);
            return
        end
        dA = R(iAp, :) - R(iAm, :);
        dB = R(iBp, :) - R(iBm, :);
        q = dA - dB;  p = dA + dB;
        good = isfinite(q);
        if nnz(good) >= 5
            % One bad average (a glitched pass) should not move the line: > 5 robust sigma.
            med = median(q(good));
            mad = 1.4826 * median(abs(q(good) - med));
            if mad > 0
                good = good & abs(q - med) <= 5 * mad;
            end
        end
        nUsed(k) = nnz(good);
        Qm(k) = mean(q(good));  Pm(k) = mean(p(good));
        if nnz(good) >= 3
            Qe(k) = std(q(good)) / sqrt(nnz(good));
            Pe(k) = std(p(good)) / sqrt(nnz(good));
        elseif ~isempty(E)
            Qe(k) = sqrt(E(iAp)^2 + E(iAm)^2 + E(iBp)^2 + E(iBm)^2);
            Pe(k) = Qe(k);
        end
    end
    f.Q = Qm;  f.Qerr = Qe;  f.P = Pm;  f.Perr = Pe;  f.nUsed = nUsed;
    f.Qnorm = Qm / f.contrast;  f.QnormErr = Qe / abs(f.contrast);
    f.Pnorm = Pm / f.contrast;

    ok = isfinite(Qm) & isfinite(Qe) & Qe > 0;
    if nnz(ok) < 2
        f.why = 'fewer than two offsets with a usable Q and error (need >= 3 averages)';
        return
    end
    x = offs(ok)';  y = Qm(ok)';  e = Qe(ok)';
    X = [ones(size(x)), x];
    W = diag(1 ./ e.^2);
    Cv = inv(X' * W * X);
    c = Cv * (X' * W * y);
    r = (y - X * c) ./ e;
    dof = numel(y) - 2;
    chi2nu = NaN;
    if dof > 0
        chi2nu = sum(r.^2) / dof;
        Cv = Cv * max(1, chi2nu);
    end
    a = c(1);  b = c(2);
    f.intercept = a;  f.slope = b;
    f.interceptErr = sqrt(Cv(1, 1));  f.slopeErr = sqrt(Cv(2, 2));
    f.chi2nu = chi2nu;
    f.slopeNorm = b / f.contrast;
    f.line = @(d) a + b * d;

    if ~(abs(b) > 4 * f.slopeErr)
        f.why = sprintf(['no measurable slope (%.3g +/- %.3g per dB): the test cannot see the ' ...
                         'balance -- too little contrast, too few averages, or wrong pulse times'], ...
                        b, f.slopeErr);
        f.dStar = NaN;  f.dStarErr = NaN;
        return
    end
    dS = -a / b;
    v = (Cv(1, 1) + dS^2 * Cv(2, 2) + 2 * dS * Cv(1, 2)) / b^2;
    f.dStar = dS;
    f.dStarErr = sqrt(max(v, 0));

    % Curvature: sin(tilt) is not a straight line far from balance. Reported, not used.
    f.curvature = NaN;  f.curvatureErr = NaN;
    if nnz(ok) >= 3
        X2 = [ones(size(x)), x, x.^2];
        Cv2 = inv(X2' * W * X2);
        c2 = Cv2 * (X2' * W * y);
        f.curvature = c2(3);  f.curvatureErr = sqrt(Cv2(3, 3));
    end
    inside = dS >= min(offs) - 0.25 * range(offs) && dS <= max(offs) + 0.25 * range(offs);
    f.inside = inside;
    f.ok = true;
    if ~inside
        f.why = sprintf(['the zero crossing (%+.2f dB) is outside the offsets played; the line ' ...
                         'is extrapolated -- rerun centred on it'], dS);
    end
end

function i = pick(cfg, d, lead, ref)
    i = 0;
    for k = 1:numel(cfg)
        if strcmp(cfg(k).kind, 'quad') && abs(cfg(k).offsetDb - d) < 1e-9 && ...
                cfg(k).lead == lead && abs(cfg(k).refPhaseDeg - ref) < 1e-9
            i = k;
            return
        end
    end
end

function m = nmean(x)
    x = x(isfinite(x));
    if isempty(x); m = NaN; else; m = mean(x); end
end

function s = nsem(x)
    x = x(isfinite(x));
    if numel(x) < 2; s = NaN; else; s = std(x) / sqrt(numel(x)); end
end

function r = range(x)
    r = max(x) - min(x);
end
