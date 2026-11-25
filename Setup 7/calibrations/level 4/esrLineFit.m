function out = esrLineFit(freqs, signal, sterr, opts)
%ESRLINEFIT Fit one CW ESR line and return centre, width, depth and area with errors.
%
%   out = esrLineFit(freqs, signal)
%   out = esrLineFit(freqs, signal, sterr)
%   out = esrLineFit(freqs, signal, sterr, opts)
%
% LEVEL 3 INFRASTRUCTURE. Not a calibration itself -- every Level-3 and Level-4
% calibration calls this, so it is written once, here, with error bars.
%
% ===================================================================
% WHY NOT JUST USE myESRFit
% ===================================================================
%
% yyMagneticFieldAlignment has a local myESRFit that fits
%
%   c - a1*g1^2/((x-f1)^2 + g1^2)
%
% which is the right lineshape. Three things are added here.
%
% 1. NO HARDCODED START POINT. myESRFit uses StartPoint = [0.0065 1 2805 4.5], so
%    the centre guess is fixed at 2805 MHz. It works for that script's sweep and
%    silently misbehaves anywhere else -- at zero field the line is at 2870, 65 MHz
%    and many linewidths away from the guess. Here every start point is derived
%    from the data.
%
% 2. A SLOPED BASELINE. myESRFit's baseline is the constant c. But the microwave
%    amplitude is a single scalar for the whole sweep while the resonator response
%    varies across it, so the off-resonance level really does tilt -- their own
%    resonator fit gives a ~20 MHz-wide Lorentzian, and yyMagneticFieldAlignment
%    sweeps 2780-2950, several linewidths of it. Forcing a flat baseline under a
%    tilted one pulls the fitted CENTRE, which is exactly the quantity the field
%    alignment depends on. Fitting a slope absorbs the tilt instead. The fitted
%    slope is returned, because a large one is itself a warning.
%
% 3. EVERYTHING RETURNED WITH AN UNCERTAINTY, including the width and the area,
%    not just the centre. Width is what the Level-3 calibrations optimise against,
%    so a width without an error bar is not usable.
%
% ===================================================================
% WHY WIDTH RATHER THAN DEPTH IS THE USEFUL OBSERVABLE
% ===================================================================
%
% For tuning the microwave chain you want a monotone, single-valued function of the
% drive field B1. Depth is not one: contrast saturates once the drive is strong
% enough to equalise the populations, so above modest power it stops responding and
% becomes a useless optimiser target. Width keeps growing -- power broadening gives
% roughly Gamma^2 = Gamma0^2 + c*Omega^2 -- so it stays monotone in B1 over a wide
% range. Area (depth x width) is intermediate and slightly more noise-tolerant.
%
% All three are returned. Which one to optimise is the caller's decision, but the
% default advice is width.
%
% ===================================================================
% ERROR BARS AND HOW MUCH TO BELIEVE THEM
% ===================================================================
%
% Parameter errors come from confint, then rescaled by sqrt(reduced chi2) -- in BOTH
% directions. This is the third time on this setup that an error bar has been
% correct for the quantity it describes and wrong for the question being asked, so
% it is worth stating once properly.
%
% The sterr AWGExpESR reports is dominated by drift BETWEEN averages. That drift is
% common to every frequency in a sweep, and calibLastDelay measured it at 42 times
% the actual point-to-point scatter. confint scales linearly with the input errors,
% so it inherits the same factor, and every parameter comes back with an error bar
% tens of times too large.
%
% An earlier version of this file scaled up on chi2 > 1 and refused to scale down on
% chi2 < 1, on the grounds that shrinking would manufacture precision. That is
% backwards, and it threw away an entire channel of calibMwDrive -- seven dips
% rising cleanly from 0.75% to 3.71%, every one rejected as "only 2 sigma". Refusing
% to scale down is not conservative when the input scale is known to be inflated; it
% just makes the test unpassable.
%
% So: trust the SHAPE of the input errors, not their SCALE, and let the residuals
% set the scale. The scaling is floored at 0.05 -- below that the input errors are
% broken rather than mis-scaled, and it should not be papered over.
%
% A KNOWN LIMIT ON ABSOLUTE WIDTHS. The AWG emits a whole number of microwave
% cycles in a burst whose sample count is rounded, so the emitted frequency can
% differ from the requested one by up to ~0.5 MHz, varying point to point. That
% jitter adds in quadrature to the true linewidth. It is negligible against a
% 20 MHz line and it is not negligible against a 3 MHz one, so treat fitted widths
% below about 5 MHz as upper bounds until that rounding is fixed.
%
% ===================================================================
% INPUTS
%   freqs    MHz
%   signal   normalized FL (the dip goes DOWN)
%   sterr    optional; used as fit weights 1/sterr^2
%   opts     .minSpanInWidths  refuse to trust a fit whose sweep is narrower than
%                              this many FWHM (default 4)
%            .plot             draw the fit (default false)
%            .name             label for the plot
%
% OUTPUT (struct)
%   centre, centreErr      MHz
%   fwhm, fwhmErr          MHz
%   depth, depthErr        in units of the signal (the dip is positive)
%   area, areaErr          depth*fwhm*pi/2
%   baseline, slope        off-resonance level and its tilt per MHz
%   chi2Reduced, birge
%   ok                     false if there is no usable line
%   why                    why not, when ok is false
%   fitobject, gof

    if ~exist('sterr', 'var'); sterr = []; end
    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('minSpanInWidths', 4, 'plot', false, 'name', ''));

    x = freqs(:); y = signal(:);
    good = isfinite(x) & isfinite(y);
    if ~isempty(sterr)
        e = sterr(:); e(~isfinite(e) | e <= 0) = NaN;
    else
        e = nan(size(x));
    end
    x = x(good); y = y(good); e = e(good);

    out = struct('ok', false, 'why', '', 'centre', NaN, 'centreErr', NaN, ...
                 'fwhm', NaN, 'fwhmErr', NaN, 'depth', NaN, 'depthErr', NaN, ...
                 'area', NaN, 'areaErr', NaN, 'baseline', NaN, 'slope', NaN, ...
                 'chi2Reduced', NaN, 'birge', NaN, 'fitobject', [], 'gof', []);

    if numel(x) < 8
        out.why = sprintf('only %d usable points; a 5-parameter fit needs more', numel(x));
        return
    end
    [x, ix] = sort(x); y = y(ix); e = e(ix);

    % ---- start points, all from the data ----------------------------------
    % Baseline and slope from the outer quarter of the sweep at each end, which is
    % where the line is not. Using ALL the points would drag the baseline down into
    % the dip and then the depth start point comes out too small.
    nOut = max(2, round(0.25*numel(x)));
    outer = [1:nOut, numel(x)-nOut+1:numel(x)];
    pB = polyfit(x(outer) - mean(x), y(outer), 1);
    slope0 = pB(1); base0 = polyval(pB, 0);

    resid0 = (base0 + slope0*(x - mean(x))) - y;      % positive inside the dip
    noise0 = 1.4826 * median(abs(resid0(outer) - median(resid0(outer))));
    noise0 = max(noise0, eps);

    depth0 = max(resid0);
    if depth0 < 3*noise0
        out.why = sprintf('no dip: deepest excursion is %.1f sigma of the baseline scatter', ...
                          depth0/noise0);
        return
    end

    % Centre from the DEFICIT-WEIGHTED CENTROID, not from argmin. argmin is one
    % sample and one glitch sets it -- the same failure mode as the min()-based
    % contrast metrics that produced a confident wrong conclusion in the AWG work.
    wCen = max(resid0, 0); wCen(wCen < 0.2*depth0) = 0;
    centre0 = sum(wCen .* x) / sum(wCen);

    % Width from where the deficit crosses half its peak, measured directly.
    aboveHalf = x(resid0 >= 0.5*depth0);
    if numel(aboveHalf) >= 2
        fwhm0 = max(aboveHalf) - min(aboveHalf);
    else
        fwhm0 = 4 * median(diff(x));
    end
    fwhm0 = max(fwhm0, 2*median(diff(x)));
    g0 = fwhm0/2;

    % ---- the fit -----------------------------------------------------------
    xc = mean(x);
    ft = fittype('c + m*(x - xc0) - a*g^2/((x-f)^2 + g^2)', ...
                 'independent', 'x', 'dependent', 'y', ...
                 'problem', 'xc0', 'coefficients', {'a','c','f','g','m'});
    fo = fitoptions('Method', 'NonlinearLeastSquares');
    fo.Display = 'Off';
    %            a        c       f              g          m
    fo.StartPoint = [depth0, base0, centre0,      g0,        slope0];
    fo.Lower      = [0,     -Inf,   min(x),       0,        -Inf];
    fo.Upper      = [Inf,    Inf,   max(x),       range(x), Inf];
    if any(isfinite(e))
        w = 1 ./ e.^2;
        w(~isfinite(w)) = median(w(isfinite(w)));
        fo.Weights = w;
    end

    try
        [fr, gof] = fit(x, y, ft, fo, 'problem', xc);
    catch err
        out.why = sprintf('fit failed: %s', err.message);
        return
    end
    out.fitobject = fr; out.gof = gof;

    cv = coeffvalues(fr);           % a c f g m, alphabetical
    ci = confint(fr);
    ce = (ci(2,:) - ci(1,:)) / (2*1.96);        % confint is 95%, want 1 sigma

    a = cv(1); c = cv(2); f = cv(3); g = cv(4); m = cv(5);

    % ---- chi2, and the one-sided Birge inflation ---------------------------
    yhat = c + m*(x - xc) - a*g^2./((x-f).^2 + g^2);
    r = y - yhat;
    dof = max(numel(x) - 5, 1);
    if any(isfinite(e))
        ee = e; ee(~isfinite(ee)) = median(e(isfinite(e)));
        chi2 = sum((r./ee).^2) / dof;
    else
        chi2 = NaN;
    end
    % TWO-SIDED. The first version inflated on chi2 > 1 and did nothing on chi2 < 1,
    % with the reasoning that common-mode drift divides out of a lineshape fit so a
    % low chi2 should not be rewarded. That reasoning was wrong, and it cost a whole
    % channel of calibMwDrive: every point came back "REJECTED: depth is only 2 sigma"
    % on dips that were plainly there and rising monotonically, 0.75% to 3.71%.
    %
    % Here is the actual situation. confint scales linearly with the input error
    % bars. The sterr from AWGExpESR is dominated by drift BETWEEN averages, which
    % is common to every frequency in the sweep -- measured at 42x the real
    % point-to-point scatter in calibLastDelay. So the input errors are ~40x too
    % large for the residuals, chi2 comes out near zero, and confint inherits the
    % same factor of 40. Refusing to scale down does not make the answer
    % conservative, it makes it useless: a 3-sigma test that can never be passed
    % rejects real data.
    %
    % The standard treatment when the SHAPE of the error bars is trusted but their
    % SCALE is not is to scale parameter errors by sqrt(reduced chi2) in both
    % directions. The residuals are the honest measure of the noise; chi2 just
    % converts the claimed scale into the real one. The constant part of the drift
    % is absorbed by the baseline parameter c, so it does not leak into the width or
    % the depth.
    %
    % Floored at 0.05, i.e. the errors can shrink by at most 20x. Below that
    % something is wrong with the input errors rather than merely mis-scaled, and
    % silently claiming 100x precision would be worse than the problem being fixed.
    birge = 1;
    if isfinite(chi2) && chi2 > 0
        birge = max(sqrt(chi2), 0.05);
    end
    ceRaw = ce;
    ce = ce * birge;
    out.chi2Reduced = chi2; out.birge = birge;
    out.errScaledBy = birge;
    out.centreErrRaw = ceRaw(3);        % before scaling, for diagnosis

    out.centre = f;      out.centreErr = ce(3);
    out.fwhm   = 2*g;    out.fwhmErr   = 2*ce(4);
    out.depth  = a;      out.depthErr  = ce(1);
    out.baseline = c;    out.slope     = m;
    % Integrated area of a Lorentzian of peak a and HWHM g is pi*a*g. Errors on a
    % and g are correlated, and confint does not give the covariance, so this adds
    % them in quadrature -- an approximation, and usually a slight overestimate
    % since depth and width trade off against each other in the fit.
    out.area = pi * a * g;
    out.areaErr = out.area * sqrt((ce(1)/max(a,eps))^2 + (ce(4)/max(g,eps))^2);

    % ---- is the answer trustworthy -----------------------------------------
    span = max(x) - min(x);
    out.ok = true;
    if out.fwhm <= 0 || ~isfinite(out.fwhm)
        out.ok = false; out.why = 'width did not converge';
    elseif span < opts.minSpanInWidths * out.fwhm
        out.ok = false;
        out.why = sprintf(['sweep spans %.1f FWHM, under the %g required. The baseline and ' ...
                           'the Lorentzian wings are degenerate at that span, so the width ' ...
                           'and depth are not separable. Widen the sweep.'], ...
                          span/out.fwhm, opts.minSpanInWidths);
    elseif out.fwhm < 3 * median(diff(x))
        out.ok = false;
        out.why = sprintf(['FWHM %.2f MHz is under 3 frequency steps (%.2f MHz). The line is ' ...
                           'not resolved. Use a finer step.'], out.fwhm, median(diff(x)));
    elseif out.depth < 3 * out.depthErr
        out.ok = false;
        out.why = sprintf('depth is only %.1f sigma', out.depth/max(out.depthErr,eps));
    elseif abs(m) * span > 0.5 * a
        % Not a failure, but worth flagging: the baseline tilt across the sweep is
        % comparable to the dip itself, which is the situation that biases centres.
        out.why = sprintf(['WARNING baseline tilts by %.1f%% of the dip depth across the sweep. ' ...
                           'The slope is fitted so the centre should be safe, but this is the ' ...
                           'resonator response and it argues for the per-frequency amplitude ' ...
                           'table.'], 100*abs(m)*span/a);
    end

    if opts.plot
        figure; hold on
        if any(isfinite(e)); errorbar(x, y, e, '.'); else; plot(x, y, '.'); end
        xf = linspace(min(x), max(x), 500);
        plot(xf, c + m*(xf-xc) - a*g^2./((xf-f).^2 + g^2), 'r-', 'LineWidth', 1.5);
        plot(xf, c + m*(xf-xc), 'k--');
        xlabel('frequency (MHz)'); ylabel('FL');
        title(sprintf('%s  centre %.3f +/- %.3f, FWHM %.2f +/- %.2f MHz', ...
              opts.name, out.centre, out.centreErr, out.fwhm, out.fwhmErr));
    end
end

function o = withDefaults(o, d)
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(o, f{i}) || isempty(o.(f{i})); o.(f{i}) = d.(f{i}); end
    end
end
