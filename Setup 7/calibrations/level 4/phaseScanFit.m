function out = phaseScanFit(phaseDeg, omega, omegaErr, opts)
%PHASESCANFIT Fit the two-channel interference curve. Exact, linear, no start point.
%
%   out = phaseScanFit(phaseDeg, fRabi, fRabiErr)
%
% Two coherent drives add as phasors, so the total Rabi frequency obeys
%
%       Omega_total^2(phi) = Omega1^2 + Omega2^2 + 2*Omega1*Omega2*cos(phi - phi0)
%
% exactly -- not to first order, exactly, for any pair of amplitudes. Expanding the
% cosine,
%
%       Omega^2 = c0 + c1*cos(phi) + c2*sin(phi)
%
% which is LINEAR in (c0, c1, c2). So this is one weighted least-squares solve. There is
% no start point to get wrong, no iteration to fail to converge, and no local optimum,
% and the fit uses every point rather than the three nearest an extremum.
%
% Fitting in Omega SQUARED rather than Omega is what makes that true. Omega itself is
% |Omega1 + Omega2*exp(i*phi)|, a square root of the above, which is not a sinusoid and
% would need a nonlinear fit for no gain.
%
% ===================================================================
% WHAT COMES OUT, AND WHICH PARTS TO TRUST
% ===================================================================
%
% From (c0, c1, c2), writing P = c0 and Q = sqrt(c1^2 + c2^2):
%
%       phiMax  = atan2(c2, c1)        the phase of maximum drive -- the operating point
%       phiNull = phiMax + 180         the antiphase null
%       P = Omega1^2 + Omega2^2        so  Omega1 + Omega2 = sqrt(P + Q)
%       Q = 2*Omega1*Omega2            and |Omega1 - Omega2| = sqrt(P - Q)
%
% phiMax and the SUM are well determined. The DIFFERENCE is not, near balance, and that
% is arithmetic rather than bad luck. Writing Omega1 = Omega(1+e), Omega2 = Omega(1-e),
%
%       (P - Q) / (P + Q)  =  e^2
%
% so a 1 % amplitude imbalance (e = 0.005) makes P - Q equal 2.5e-5 of P + Q. Extracting
% that from two fitted numbers each known to a tenth of a percent is impossible, and
% P - Q will often come out slightly NEGATIVE, which is not an error -- it means the
% channels are balanced better than this scan can resolve. The imbalance is then reported
% as an upper bound rather than a value, and the honest measurement of a small imbalance
% is a direct one at the null, not this.
%
% The phase optimum, by contrast, does NOT depend on the balance: the expression above is
% extremal at phi0 and phi0+180 for ANY Omega1 and Omega2. So this fit gives a trustworthy
% phase even on a badly balanced pair, which is what makes the phase safe to calibrate
% first.
%
% CALIBRATE ON THE NULL, NOT THE PEAK. Both are found here at once because the whole curve
% is fitted, but note why the peak alone would be a poor target: d(Omega^2)/d(phi) vanishes
% at both extrema, while the FRACTIONAL change is huge at the null and tiny at the peak. A
% 20 degree phase error costs 1.5 % of the drive at the peak and would be invisible; the
% same error at the null raises Omega from near zero to 0.35 of a single channel.
%
% OPTIONS
%   secondHarmonic  also fit cos(2phi), sin(2phi) and report whether they are significant
%                   (default true). A pure two-drive interference has NO second harmonic,
%                   so a significant one means the model is incomplete -- a nonlinearity
%                   both channels share, or a phase that is not doing what we think.
%   verbose         print the report (default true)
%
% OUTPUT
%   ok, why, phiMax, phiMaxErr, phiNull, sumOmega, sumOmegaErr, diffOmega, diffOmegaErr,
%   diffBound, imbalanceDb, imbalanceDbBound, balanceResolved, visibility,
%   chi2Reduced, r2, harmonic2, harmonic2Sigma, model
%
% SELF-CONTAINED: base MATLAB only.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('secondHarmonic', true, 'verbose', true), 'phaseScanFit');

    phaseDeg = phaseDeg(:); omega = omega(:);
    if isempty(omegaErr); omegaErr = ones(size(omega)); else; omegaErr = omegaErr(:); end

    out = struct('ok', false, 'why', '', 'opts', opts, ...
                 'phiMax', NaN, 'phiMaxErr', NaN, 'phiNull', NaN, ...
                 'sumOmega', NaN, 'sumOmegaErr', NaN, ...
                 'diffOmega', NaN, 'diffOmegaErr', NaN, 'diffBound', NaN, ...
                 'imbalanceDb', NaN, 'imbalanceDbBound', NaN, ...
                 'balanceResolved', false, 'visibility', NaN, ...
                 'chi2Reduced', NaN, 'r2', NaN, ...
                 'harmonic2', NaN, 'harmonic2Sigma', NaN);

    keep = isfinite(phaseDeg) & isfinite(omega) & isfinite(omegaErr) & omega > 0 & omegaErr > 0;
    phaseDeg = phaseDeg(keep); omega = omega(keep); omegaErr = omegaErr(keep);
    if numel(omega) < 4
        out.why = sprintf('%d usable phase points, need at least 4 for a 3-parameter fit', ...
                          numel(omega));
        if opts.verbose && ~isempty(out.why); fprintf('\nphaseScanFit: %s\n', out.why); end
        return
    end
    if max(phaseDeg) - min(phaseDeg) < 90
        out.why = sprintf(['the phases span only %.0f degrees. Below about 180 the three ' ...
            'parameters are not separable -- a shifted curve and a taller one fit the same ' ...
            'arc.'], max(phaseDeg) - min(phaseDeg));
        if opts.verbose; fprintf('\nphaseScanFit: %s\n', out.why); end
        return
    end

    ph = deg2rad(phaseDeg);
    y  = omega.^2;
    sy = 2*omega.*omegaErr;             % exact to first order, and the Rabi frequency is
                                        % known to well under a percent, so that is enough
    X  = [ones(size(ph)), cos(ph), sin(ph)];
    [c, C, chi2r] = wls(X, y, sy);
    if isempty(c); out.why = 'the weighted solve failed'; return; end

    P = c(1); Q = hypot(c(2), c(3));
    out.chi2Reduced = chi2r;
    resid = y - X*c;
    out.r2 = 1 - sum(resid.^2)/max(sum((y - mean(y)).^2), eps);
    out.resid = resid;

    if Q <= 0 || ~isfinite(Q)
        out.why = 'no phase dependence at all -- the two channels are not interfering';
        if opts.verbose; fprintf('\nphaseScanFit: %s\n', out.why); end
        return
    end

    out.phiMax  = mod(rad2deg(atan2(c(3), c(2))), 360);
    out.phiNull = mod(out.phiMax + 180, 360);
    % d(phiMax)/dc = [0, -c2/Q^2, c1/Q^2], in radians.
    g = [0, -c(3)/Q^2, c(2)/Q^2];
    out.phiMaxErr = rad2deg(sqrt(max(g*C*g', 0)));

    % ---- sum and difference of the two channel amplitudes ----------------
    sSq = P + Q; dSq = P - Q;
    gS = [1, c(2)/Q, c(3)/Q];               % d(P+Q)/dc
    gD = [1, -c(2)/Q, -c(3)/Q];             % d(P-Q)/dc
    vS = max(gS*C*gS', 0); vD = max(gD*C*gD', 0);

    out.sumOmega = sqrt(max(sSq, 0));
    out.sumOmegaErr = sqrt(vS)/max(2*out.sumOmega, eps);
    out.visibility = sqrt(max(dSq, 0))/max(out.sumOmega, eps);

    % P - Q is (Omega1 - Omega2)^2, a difference of two nearly equal fitted numbers. It is
    % reported as a VALUE only when it is significantly positive, and as an upper BOUND
    % otherwise -- because a negative or marginal P - Q means "better balanced than this can
    % see", which is a result, and printing sqrt of a negative number as 0.000 would hide it.
    sigD = sqrt(vD);
    out.balanceResolved = isfinite(sigD) && dSq > 2*sigD && dSq > 0;
    if out.balanceResolved
        out.diffOmega = sqrt(dSq);
        out.diffOmegaErr = sigD/max(2*out.diffOmega, eps);
        out.diffBound = out.diffOmega + 2*out.diffOmegaErr;
    else
        out.diffOmega = NaN; out.diffOmegaErr = NaN;
        out.diffBound = sqrt(max(dSq, 0) + 2*max(sigD, 0));   % 2-sigma bound on |O1 - O2|
    end
    out.imbalanceDb      = dbFromSumDiff(out.sumOmega, out.diffOmega);
    out.imbalanceDbBound = dbFromSumDiff(out.sumOmega, out.diffBound);

    % ---- is a second harmonic there? -------------------------------------
    % Two coherent drives produce a pure first harmonic in Omega^2. Anything at 2*phi is
    % not interference: it is a nonlinearity the pair shares, or the phase knob not doing
    % what it says. Cheap to test and it is the only internal check on the model.
    if opts.secondHarmonic && numel(y) >= 6
        X2 = [X, cos(2*ph), sin(2*ph)];
        [c2f, C2, ~] = wls(X2, y, sy);
        if ~isempty(c2f)
            amp2 = hypot(c2f(4), c2f(5));
            v2 = C2(4,4)*c2f(4)^2 + C2(5,5)*c2f(5)^2 + 2*C2(4,5)*c2f(4)*c2f(5);
            s2 = sqrt(max(v2, 0))/max(amp2, eps);
            out.harmonic2 = amp2/max(Q, eps);          % as a fraction of the real harmonic
            out.harmonic2Sigma = amp2/max(s2, eps);
        end
    end

    out.model = @(pd) c(1) + c(2)*cos(deg2rad(pd(:))) + c(3)*sin(deg2rad(pd(:)));
    out.coeffs = c; out.cov = C; out.P = P; out.Q = Q;
    out.ok = true;
    if opts.verbose; report(out); end
end

% ======================================================================

function db = dbFromSumDiff(s, d)
% Convert (Omega1+Omega2, |Omega1-Omega2|) into the dB offset between the channels.
% Omega1 = (s+d)/2, Omega2 = (s-d)/2, and dB is 20*log10 of an amplitude ratio.
    db = NaN;
    if ~isfinite(s) || ~isfinite(d) || s <= 0 || d < 0 || d >= s; return; end
    db = 20*log10((s + d)/(s - d));
end

function [c, C, chi2r] = wls(X, y, sy)
% Weighted least squares with a Birge-scaled covariance. Column-scaled before inverting,
% because the columns here are O(1) while y is O(100) and an unscaled normal matrix loses
% significant figures.
    c = []; C = []; chi2r = NaN;
    Xw = X ./ sy; yw = y ./ sy;
    c = Xw \ yw;
    if any(~isfinite(c)); c = []; return; end
    dof = max(numel(y) - size(X,2), 1);
    r = Xw*c - yw;
    chi2r = sum(r.^2)/dof;
    birge = max(sqrt(chi2r), 0.05);
    C = nan(size(X,2));
    try
        M = Xw'*Xw;
        sc = sqrt(diag(M))'; sc(sc == 0 | ~isfinite(sc)) = 1;
        Ms = M ./ (sc'*sc);
        if rcond(Ms) > 1e-14
            C = (inv(Ms) ./ (sc'*sc)) * birge^2;
        end
    catch
    end
end

function report(o)
    fprintf('\n--- two-channel interference fit ---\n');
    fprintf('  maximum drive at phase   %.2f +/- %.2f deg   <- the operating point\n', ...
            o.phiMax, o.phiMaxErr);
    fprintf('  null at                  %.2f deg\n', o.phiNull);
    fprintf('  Omega1 + Omega2          %.3f +/- %.3f MHz\n', o.sumOmega, o.sumOmegaErr);
    fprintf('  null visibility          %.4f   (0 = perfectly balanced and cancelling)\n', ...
            o.visibility);
    if o.balanceResolved
        fprintf('  |Omega1 - Omega2|        %.3f +/- %.3f MHz  ->  imbalance %.2f dB\n', ...
                o.diffOmega, o.diffOmegaErr, o.imbalanceDb);
    else
        fprintf(['  |Omega1 - Omega2|        below what this scan resolves; 2-sigma bound\n' ...
                 '                           %.3f MHz, i.e. imbalance under %.2f dB\n'], ...
                o.diffBound, o.imbalanceDbBound);
        fprintf(['    That is expected near balance and is not a failure. (P-Q)/(P+Q) equals\n' ...
                 '    the squared fractional imbalance, so 1 %% imbalance hides at 2.5e-5 of\n' ...
                 '    the fitted scale. A small imbalance has to be measured directly, one\n' ...
                 '    channel at a time, not squeezed out of this fit.\n']);
    end
    fprintf('  R^2 %.5f, reduced chi2 %.2f\n', o.r2, o.chi2Reduced);

    if isfinite(o.harmonic2Sigma)
        fprintf('  second harmonic          %.1f %% of the first, %.1f sigma\n', ...
                100*o.harmonic2, o.harmonic2Sigma);
        if o.harmonic2Sigma > 3 && o.harmonic2 > 0.05
            fprintf(['    *** A SIGNIFICANT SECOND HARMONIC SHOULD NOT EXIST. Two coherent\n' ...
                     '    drives give Omega^2 = P + Q*cos(phi - phi0) exactly, with nothing at\n' ...
                     '    2*phi. Something is off the model: a shared nonlinearity, a phase\n' ...
                     '    knob that also changes amplitude, or drift across the scan that\n' ...
                     '    happened to correlate with phase. Check the point order.\n']);
        end
    end
    if o.chi2Reduced > 5
        fprintf(['\n  Reduced chi2 %.1f: the points scatter well beyond their error bars. On a\n' ...
                 '  model that is exact, that points at drift between points rather than at the\n' ...
                 '  fit -- the amplitude or the field moving during the scan.\n'], o.chi2Reduced);
    end
    if o.phiMaxErr > 10
        fprintf(['\n  The phase is only determined to %.0f degrees. At the peak that costs little\n' ...
                 '  drive (a 20 degree error is 1.5 %%), but if you need the NULL to be deep --\n' ...
                 '  for a balance measurement -- take more phase points or more averages.\n'], ...
                o.phiMaxErr);
    end
end

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d); given = fieldnames(o);
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
