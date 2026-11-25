function rel = coilLineRelation(rho, fLine, fLineErr, opts)
%COILLINERELATION Fit (coil setting -> NV line frequency) from a scan, and invert it.
%
%   rel = coilLineRelation(rho, fLine, fLineErr)
%   rel = coilLineRelation(rho, fLine, fLineErr, struct('tag', 'coarse'))
%   [rhoWanted, err] = rel.invert(2823.0, 0.4)
%
% This is the coil calibration in the only form a frequency target needs: which Brho puts
% the line at a given frequency. It is fitted from the scan that is asking the question,
% not imported from an earlier one, and that is deliberate.
%
% ===================================================================
% WHY IT IS REFITTED EVERY TIME RATHER THAN SAVED
% ===================================================================
%
% The slope and the offset of a straight line through a narrow x range are strongly
% collinear -- change the slope and you can absorb almost all of it by moving the offset.
% Over 15 to 20.6 G the correlation between them is 0.995. Two honest runs of calibB1Peak
% produced
%
%       |B| = 0.974 * Brho - 0.813          and          |B| = 1.000 * Brho - 1.176
%
% Both fit their own data well. They agree to better than 0.1 G near 18 G, where the data
% were, and diverge outside it. Neither slope is "the" coil ratio and neither offset is
% "the" offset; the PAIR is an interpolator over the range measured, and that is all.
%
% So: fit it here, use it to interpolate inside the range just scanned, and let the caller
% check whether the answer landed inside that range. Never save the coefficients and reuse
% them, and never quote the slope on its own as a calibration constant.
%
% ===================================================================
% THE SLOPE IS PREDICTED, SO IT IS A CHECK AND NOT JUST A FIT
% ===================================================================
%
% f_line = D - gamma*|B| for the lower transition of the most-aligned orientation, and
% |B| = k*Brho + m for any linear coil, so
%
%       slope = -gamma * k  =  about -2.72 MHz/G for k near 0.97
%
% A slope far from that does not mean a poor fit -- it means the points are not all
% following the same line. That happens when a field point's "lowest line" belongs to a
% different NV orientation, which is exactly the failure that once turned a coil fit into
% 0.13*Brho + 9.7 and dragged a reported peak 18 MHz. It is checked here rather than left
% to be noticed downstream.
%
% OPTIONS
%   gamma        MHz/G (default 2.8025)
%   coilRatio    expected |B|/Brho, only for the slope sanity check (default 0.97)
%   ratioLo, ratioHi   how far the implied ratio may stray before slopeSane goes false
%                      (defaults 0.7, 1.15)
%   verbose      print the report (default true)
%   tag          a label for the report
%
% OUTPUT
%   ok, why      did it fit
%   a, b         f_line = a*rho + b
%   cov          2x2 covariance of [a b], Birge-scaled
%   n, rms, chi2Reduced
%   coilRatio    -a/gamma, the implied |B|/Brho
%   slopeSane    is that ratio physically plausible
%   invert       function handle: [rho, rhoErr] = rel.invert(fTarget, fTargetErr)
%   range        [min max] of the rho values fitted, for the caller's extrapolation check
%
% SELF-CONTAINED: base MATLAB only.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('gamma', 2.8025, 'coilRatio', 0.97, ...
                                     'ratioLo', 0.7, 'ratioHi', 1.15, ...
                                     'verbose', true, 'tag', ''), 'coilLineRelation');

    rho = rho(:); fLine = fLine(:);
    if isempty(fLineErr); fLineErr = ones(size(rho)); end
    fLineErr = fLineErr(:);

    rel = struct('ok', false, 'why', '', 'a', NaN, 'b', NaN, 'cov', nan(2), ...
                 'n', 0, 'rms', NaN, 'chi2Reduced', NaN, 'coilRatio', NaN, ...
                 'slopePredicted', -opts.gamma*opts.coilRatio, 'slopeSane', false, ...
                 'range', [NaN NaN], 'opts', opts, 'invert', @(f, e) deal(NaN, NaN));

    keep = isfinite(rho) & isfinite(fLine);
    rho = rho(keep); fLine = fLine(keep); fLineErr = fLineErr(keep);
    if numel(rho) < 3
        rel.why = sprintf('%d usable line positions, need at least 3', numel(rho));
        if opts.verbose; report(rel); end
        return
    end
    % A missing or zero error bar becomes the median rather than dropping the point. A
    % fitted centre with no error is still a measured centre; a ZERO error would give it
    % infinite weight, which is the real hazard.
    okE = isfinite(fLineErr) & fLineErr > 0;
    if ~any(okE)
        fLineErr = ones(size(fLineErr));
    else
        fLineErr(~okE) = median(fLineErr(okE));
    end
    if numel(unique(rho)) < 2
        rel.why = 'all the points are at the same field';
        if opts.verbose; report(rel); end
        return
    end

    X  = [rho, ones(size(rho))];
    Xw = X ./ fLineErr;  yw = fLine ./ fLineErr;
    c  = Xw \ yw;
    if any(~isfinite(c))
        rel.why = 'the weighted solve failed';
        if opts.verbose; report(rel); end
        return
    end

    resid = fLine - X*c;
    dof   = max(numel(rho) - 2, 1);
    chi2r = sum((resid./fLineErr).^2)/dof;
    birge = max(sqrt(chi2r), 0.05);         % two-sided and floored, as everywhere else
    C = nan(2);
    try
        JtJ = Xw'*Xw;
        if rcond(JtJ) > 1e-14; C = inv(JtJ) * birge^2; end
    catch
    end

    rel.a = c(1); rel.b = c(2); rel.cov = C;
    rel.n = numel(rho); rel.rms = sqrt(mean(resid.^2)); rel.chi2Reduced = chi2r;
    rel.range = [min(rho) max(rho)];
    rel.coilRatio = -rel.a/opts.gamma;
    rel.slopeSane = rel.a < 0 && rel.coilRatio > opts.ratioLo && rel.coilRatio < opts.ratioHi;
    rel.ok = true;
    rel.resid = resid;
    % The x values the residuals belong to. Without this a caller plotting resid has to
    % reconstruct which points were used, and its own "good" mask is not necessarily the
    % same one -- calibB1Rabi fits this on (finite fLine AND aligned) but plots on P.good.
    rel.rho = rho(:)';

    % Closure over the fitted values, so the caller cannot invert a relation that has not
    % been fitted, and cannot accidentally invert it with someone else's coefficients.
    a = rel.a; b = rel.b; Cc = C;
    rel.invert = @(fT, fTErr) invertLine(a, b, Cc, fT, fTErr);

    if opts.verbose; report(rel); end
end

% ======================================================================

function [rho, rhoErr] = invertLine(a, b, C, fT, fTErr)
% rho = (fT - b)/a, with the a-b covariance carried through.
%
% The cross term is not a refinement, it is most of the answer. In any straight-line fit a
% steeper line and a shifted line are strongly anticorrelated, and adding the variances of
% a and b in quadrature while ignoring cov(a,b) overstates the error on the inverse by a
% factor of several -- which then makes a perfectly determined field setting look like a
% guess.
    rho = NaN; rhoErr = NaN;
    if ~isfinite(a) || a == 0 || ~isfinite(fT); return; end
    rho = (fT - b)/a;
    if nargin < 5 || isempty(fTErr) || ~isfinite(fTErr); fTErr = 0; end
    if all(isfinite(C(:)))
        % d(rho)/d(fT) = 1/a,  d/db = -1/a,  d/da = -rho/a
        v = (fTErr^2 + rho^2*C(1,1) + C(2,2) + 2*rho*C(1,2)) / a^2;
        rhoErr = sqrt(max(v, 0));
    end
end

function report(rel)
    tag = rel.opts.tag;
    if isempty(tag); tag = 'line position vs field'; else; tag = ['line position vs field, ' tag]; end
    fprintf('\n--- %s ---\n', tag);
    if ~rel.ok
        fprintf('  no relation: %s\n', rel.why);
        return
    end
    fprintf('  f_line = %.4f * Brho %+.2f MHz   (%d points over %.2f to %.2f G)\n', ...
            rel.a, rel.b, rel.n, rel.range(1), rel.range(2));
    fprintf('  rms residual %.2f MHz, reduced chi2 %.1f\n', rel.rms, rel.chi2Reduced);
    fprintf('  implied |B|/Brho = %.4f   (slope expected near %.3f MHz/G)\n', ...
            rel.coilRatio, rel.slopePredicted);
    if ~rel.slopeSane
        fprintf(['  *** THE SLOPE IS NOT PHYSICAL. f_line = D - gamma*|B| fixes it near %.2f\n' ...
                 '  MHz/G for any sane coil, and %.3f is not that. The scan is not following one\n' ...
                 '  orientation''s line across the fields -- most likely some points picked up a\n' ...
                 '  different orientation''s line. Nothing downstream of this relation is\n' ...
                 '  trustworthy until that is fixed.\n'], rel.slopePredicted, rel.a);
    end
    if rel.chi2Reduced > 10
        fprintf(['  Reduced chi2 %.0f: the positions do not lie on a straight line within their\n' ...
                 '  error bars. Either the coils are nonlinear over this range, or a point is\n' ...
                 '  tracking the wrong line. The inversion assumes straight, so check the\n' ...
                 '  residuals (rel.resid) before using it far from the middle of the range.\n'], ...
                rel.chi2Reduced);
    end
end

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d);
    given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('%s: unknown option "%s". Valid: %s.', who, bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || isempty(o.(known{i}))
            o.(known{i}) = d.(known{i});
        end
    end
end
