function p = walkAuditDecompose(c, cE, w, wE, g, noiseVar)
%WALKAUDITDECOMPOSE Split (walk - classic) into a contrast scale and a gate-tail tilt.
%
%   p = walkAuditDecompose(classic, classicErr, walk, walkErr, gateTail)
%   p = walkAuditDecompose(classic, classicErr, walk, walkErr, [])   % no tilt column
%   p = walkAuditDecompose(c, cE, w, wE, g, noiseVar)   % MEASURED reference noise variance
%
% ===================================================================
% WHY A DECOMPOSITION RATHER THAN TWO FITS
% ===================================================================
%
% The usual way to ask "does walk lose contrast" is to fit both traces and subtract the
% two fitted amplitudes. That throws away WHERE the difference sits, and where it sits is
% the only thing that distinguishes the two candidate mechanisms.
%
% So the residual is projected onto the two shapes those mechanisms produce:
%
%       walk - classic  =  offset  +  scale * (classic - mean classic)  +  tilt * gateTail
%
%   scale   a pure CONTRAST SCALE. walk = classic scaled by (1 + scale) about its own
%           mean. scale = -0.55 is the halved-contrast signature Ramsey showed. This is
%           what a uniform loss of coherence, or a reference window that has picked up
%           some signal, looks like -- and it is INDEPENDENT of the swept variable.
%
%   tilt    a ramp along gateTail = (max(tau) - tau)/max(tau), mean-removed. This is the
%           signature of the FIXED MW GATE specifically. In walk mode the gate is built at
%           max(tau) and never rewritten, so the interval with the switch open and the AWG
%           idling is longest at short tau and zero at max(tau). Anything leaking during
%           that interval MUST be tau-dependent, so it must appear here and it cannot
%           appear as a flat scale.
%
%           Pass g = [] to drop the column. ESR gets [], because ESR's pulse-generator
%           sequence does not depend on the swept variable at all -- there is no gate to
%           tilt along. tilt then returns NaN, deliberately, not 0: a printed 0 would read
%           as evidence about a gate ESR does not have.
%
%   offset  fitted and then ignored. A common-mode shift in the ratio is the laser or the
%           reference level, not the pulses. It is in the model so that it cannot leak
%           into scale, which is the one thing that would make this whole audit cry wolf
%           on every pair of runs taken ten minutes apart.
%
% The columns are not orthogonal on real data -- a damped oscillation has a downward trend
% of its own, which is a partial tilt -- so the fit is a genuine joint solve rather than
% three separate projections, and each coefficient's error accounts for its correlation
% with the others. testWalkAudit injects a known scale, a known tilt, and both together,
% and checks that neither column absorbs the other.
%
% ===================================================================
% REGRESSION DILUTION, AND THE CORRECTION FOR IT
% ===================================================================
%
% The scale column is (classic - mean classic), and the classic trace is itself MEASURED --
% it carries noise. Regressing on a noisy predictor biases the slope toward zero, and here
% that bias is not toward zero, it is NEGATIVE: it invents a contrast loss.
%
% Write s for the true structure, so classic = s + e1 and walk = m + (1+k)s + e2. Then
% y = walk - classic = m + k*s + (e2 - e1), and regressing y on x = classic gives
%
%       slope = ( k*var(s) - var(e1) ) / ( var(s) + var(e1) )
%
% At k = 0 -- walk and classic identical -- that is NOT zero, it is -var(e1)/var(x). A trace
% whose noise variance is a quarter of its structure variance therefore reports a 20%
% contrast loss that does not exist. testContrastMetric caught this: on a deliberately
% noise-heavy 34-point grid the uncorrected scale read -0.446 for a pair with IDENTICAL
% contrast.
%
% The correction is exact and needs nothing that is not already known. With V = var(x)
% measured, Ve = mean(cE.^2) the known noise variance, and Vs = V - Ve,
%
%       k = ( slope*V + Ve ) / Vs
%
% and the error scales by V/Vs. Both are reported: .scale is corrected, .scaleRaw is the bare
% slope, and .attenuation = Vs/V says how much correcting mattered. When attenuation drops
% below 0.8 the reference trace is too noisy for its own structure and .why says so -- the
% right response is more repeats on the reference, not a bigger correction.
%
% DO NOT TRUST sterr FOR THIS. Ve defaults to mean(cE.^2), and sterr on this setup is not
% reliable -- ESR's runs 18x low. Both directions hurt: too small under-corrects and leaves
% a negative bias, too large over-corrects and manufactures "walk is BETTER than classic".
%
% So pass noiseVar when you can measure it, and both callers can. With two independent
% classic runs A and D of the same thing, A - D contains only noise, so
%
%       per-run noise variance      v = var(A - D) / 2
%       noise variance of their mean  = v / 2
%
% with no reference to sterr at all. That is what ramseyGateTest and walkAudit pass, and it is
% why their contrast scales survive an error-bar factor of 18.
%
% THE ONE SAFETY PROPERTY WORTH REMEMBERING: the dilution bias is always NEGATIVE. So
% .scaleRaw is a lower bound on the true scale whatever the noise estimate does, and a CLEAN
% verdict read off .scaleRaw is trustworthy even if the correction is wrong.
%
% WEIGHTS AND ERRORS
%   Weighted by the combined per-point sigma, sqrt(cE^2 + wE^2). The covariance is Birge
%   scaled by the reduced chi2 when that exceeds 1, which is the conservative choice: if
%   the residual scatter is larger than the error bars then the error bars are understated,
%   and without the scaling a 3x understated sterr turns a 1 sigma effect into a 3 sigma
%   "discovery". That has happened on this setup before -- AWGExpESR's sterr read 42x the
%   point-to-point scatter -- so the correction is not hypothetical.
%
% FIELDS
%   scale, scaleErr    the CORRECTED contrast scale and its error
%   scaleRaw           the uncorrected slope, for comparison
%   attenuation        Vs/V. 1 means the correction did nothing; below 0.8 is a warning
%   tilt, tiltErr, offset, chi2, resid, model, sigma, names, ok, why

    c = c(:); cE = cE(:); w = w(:); wE = wE(:);
    p = struct('scale', NaN, 'scaleErr', NaN, 'scaleRaw', NaN, 'attenuation', NaN, ...
               'tilt', NaN, 'tiltErr', NaN, ...
               'offset', NaN, 'chi2', NaN, 'resid', [], 'model', [], 'sigma', [], ...
               'names', {{}}, 'ok', false, 'why', '', 'noiseFrom', 'sterr');

    if numel(c) ~= numel(w)
        p.why = sprintf('%d classic points against %d walk points', numel(c), numel(w));
        return
    end
    if numel(c) < 5
        p.why = sprintf('only %d points -- too few to separate a scale from a tilt', numel(c));
        return
    end

    y = w - c;

    % ---- the weights -----------------------------------------------------
    if numel(cE) ~= numel(c); cE = zeros(size(c)); end
    if numel(wE) ~= numel(w); wE = zeros(size(w)); end
    sig = sqrt(cE.^2 + wE.^2);
    good = isfinite(sig) & sig > 0;
    if ~any(good)
        sig = ones(size(y));                     % unweighted rather than refusing: the
    else                                         % shape information is still there
        sig(~good) = median(sig(good));
    end

    % ---- the design matrix ------------------------------------------------
    cols = [ones(numel(y),1), c - mean(c)];
    names = {'offset', 'scale'};
    if ~isempty(g)
        gg = g(:);
        if numel(gg) ~= numel(y)
            p.why = sprintf('gateTail has %d entries, the traces have %d', numel(gg), numel(y));
            return
        end
        cols = [cols, gg - mean(gg)];
        names{end+1} = 'tilt';
    end
    p.names = names;

    % ---- weighted normal equations ----------------------------------------
    % Formed explicitly and the conditioning CHECKED. A rank-deficient solve through
    % backslash returns a confident coefficient with no warning worth noticing, and the
    % degenerate case here is real: a classic trace with no structure (all points equal,
    % e.g. a run where the resonance was outside the span) makes column 2 identically
    % zero, and the reported scale would then be arbitrary.
    iv = 1 ./ sig.^2;
    A = cols' * (cols .* iv);
    b = cols' * (y .* iv);
    if rcond(A) < 1e-12
        p.why = ['the design matrix is rank deficient. The usual cause is a classic trace ' ...
                 'with no structure in it -- an ESR sweep that missed the resonance, or a ' ...
                 'Rabi with no visible oscillation -- and then there is no contrast for a ' ...
                 'scale to be a scale OF.'];
        return
    end
    th = A \ b;

    resid = y - cols*th;
    dof = max(numel(y) - numel(th), 1);
    chi2 = sum((resid ./ sig).^2) / dof;
    Cov = inv(A) * max(chi2, 1);
    se = sqrt(abs(diag(Cov)));

    p.offset   = th(1);
    p.scaleRaw = th(2);

    % ---- undo the regression dilution ------------------------------------
    x  = cols(:, 2);
    V  = mean(x.^2);                 % x is already mean-removed
    if exist('noiseVar', 'var') && ~isempty(noiseVar) && isfinite(noiseVar) && noiseVar >= 0
        Ve = noiseVar;               % MEASURED, from two independent runs of the reference
        p.noiseFrom = 'measured';
    else
        Ve = mean(cE.^2);            % fallback: believe sterr, and see the warning above
        p.noiseFrom = 'sterr';
    end
    Vs = V - Ve;
    if V <= 0 || Vs <= 0.05*V
        % The reference has no more structure than its own noise. Correcting here would
        % divide by something indistinguishable from zero and return a spectacular number.
        p.attenuation = max(Vs, 0) / max(V, eps);
        p.scale    = th(2);
        p.scaleErr = se(2);
        p.why = sprintf(['the reference trace is dominated by its own noise (structure ' ...
                         'variance is %.0f%% of the total), so the contrast scale cannot be ' ...
                         'corrected for regression dilution and is biased LOW. Raise repeats ' ...
                         'on the reference run.'], 100*max(Vs,0)/max(V,eps));
    else
        p.attenuation = Vs / V;
        p.scale    = (th(2)*V + Ve) / Vs;
        p.scaleErr = se(2) * V / Vs;
        if p.attenuation < 0.8
            p.why = sprintf(['regression dilution was significant (attenuation %.2f): the ' ...
                             'uncorrected scale was %+.3f and the corrected one is %+.3f. ' ...
                             'The correction is exact but it inflates the error bar by %.2fx, ' ...
                             'so more repeats on the reference would be worth more than the ' ...
                             'correction.'], p.attenuation, th(2), p.scale, V/Vs);
        end
    end
    if numel(th) >= 3
        p.tilt    = th(3);
        p.tiltErr = se(3);
    end
    p.chi2  = chi2;
    p.resid = resid(:)';
    p.model = (cols*th)';
    p.sigma = sig(:)';
    p.ok    = true;
end
