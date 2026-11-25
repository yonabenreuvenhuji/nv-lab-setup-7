function out = esrFitRobust(args)
%ESRFITROBUST Robust ESR fit. ONE STRUCT IN, ONE STRUCT OUT.
%
%   out = esrFitRobust(struct('exp', a))
%   out = esrFitRobust(struct('exp', a, 'nLines', 1))
%   out = esrFitRobust(struct('f', f, 'y', y, 'e', e))        % plain vectors
%
% The data goes in the same struct as the options and everything comes back in `out`, with
% out.ok saying whether to believe it. The old positional outputs are kept as NAMED fields --
% out.f0, out.width, out.contrast, out.SNR, out.coeffError -- so nothing that read them has
% to change, and out.lines is unchanged.
%
% Note out.f is the FREQUENCY AXIS, as it always was, not the first line centre. The centres
% are out.f0 and out.lines(k).centre.
%
%ESRFITROBUST Robust ESR fit. A front end to esrMultiLineFit, not a new fitter.
%
%   (see the usage block at the top -- one struct in, one struct out)
%
% out.f0, out.width and out.contrast carry what ESRfit returned positionally, and with more
% than one line they are vectors ordered by frequency. out.lines is the per-line detail.
%
% ===================================================================
% THE FITTING IS ALREADY WRITTEN. THIS IS THE PART THAT WAS MISSING.
% ===================================================================
%
% esrMultiLineFit has been in level 3 since the field survey was built, and it is already
% the fit this setup needs: N Lorentzians on a polynomial baseline, seeded from a robust
% envelope rather than a guess, Gauss-Newton with an explicit covariance and two-sided
% Birge scaling, and per-line depth and width gates so a component that fits nothing gets
% dropped and the rest refitted without it. It has been graded against synthetic spectra
% with known fields by testFieldPipeline since the day two real bugs in it were found that
% way.
%
% What it did not have was a way to call it on an experiment. It takes (freqs, signal,
% sterr) as vectors, so the only route from an AWGExpESR to it was to write the extraction
% out by hand at every call site. That is what this file is: extraction, then the fit,
% then a plot. The fitting is unchanged.
%
% So the answer to "can we just use esrMultiLineFit" is yes, and it belongs in
% functions\ rather than level 3, because it now has two callers -- calibFieldSurvey and
% this.
%
% ===================================================================
% WHY NOT extract1Ddata, AND NOT extractESRdata
% ===================================================================
%
% extractESRdata is a four-line wrapper whose whole body is a call to extract1Ddata (or
% extract2Ddata when mirrorSweepAround is set), so the objection is really to that.
%
% extract1Ddata applies three successive cuts to the ratio across averages, at 5, then 3,
% then ONE standard deviation:
%
%       s(abs(s - mean(s,2)) > 1*mean(std(s,1,2))) = NaN;
%
% A one-sigma cut discards about a third of perfectly good Gaussian data. Worse, it runs
% after two harsher passes have already shrunk the spread it measures against, and it
% biases the survivors toward the mean -- so the scatter comes out smaller than it is and
% every error bar downstream is too small. rabiFitRobust refused it for the same reason
% and this does too.
%
% What happens here instead: the reduced signal and its sterr are read from signalParam,
% and where the per-average ratios are reachable the error is taken as the LARGER of the
% propagated error and the scatter across averages. Those answer different questions --
% photon statistics within an average, versus everything that drifts between them -- and
% on this setup the second is routinely an order of magnitude bigger. That is the same
% drift that made AWGExpESR's sterr read 42x the point-to-point scatter. Taking the larger
% is the only choice that does not understate the uncertainty. Nothing is discarded.
%
% ===================================================================
% OPTIONS
% ===================================================================
%
% Everything esrMultiLineFit accepts is passed straight through: nLines, maxLines,
% baselineOrder, snr, minSepMHz, minDepthSigma, maxFwhmMHz, seedCentres. Plus:
%
%   baselineOrder   polynomial order of the shared background (default 2). The fitter
%                   documents 1 for a narrow window and 2 for sweeps past ~100 MHz, where
%                   the resonator envelope curves; 2 is the safe default for a run-script
%                   sweep. Higher orders work but start competing with the lines.
%   plotIt      draw it (default true)
%   verbose     print the per-line table (default true)
%   avgOutlierMad   drop an AVERAGE whose ratio is this many MADs from the median, and
%                   only when there are at least 5 averages to judge with (default 4,
%                   set Inf to keep everything). Note the contrast with the one-sigma
%                   cut above: a MAD rule on 5+ samples, applied once, to whole averages
%                   rather than to individual points.
%
% SELF-CONTAINED apart from esrMultiLineFit. No Curve Fitting Toolbox, no extract1Ddata,
% no findpeaks.

    % ---- ONE STRUCT IN ---------------------------------------------------
    DEFAULTS = struct('nLines', [], 'maxLines', 8, 'baselineOrder', 2, 'snr', 4, ...
                      'minSepMHz', 6, 'minDepthSigma', 3, 'maxFwhmMHz', 40, ...
                      'seedCentres', [], 'plotIt', true, 'verbose', true, ...
                      'avgOutlierMad', 4, 'exp', [], 'f', [], 'y', [], 'e', [], ...
                      'showResid', false, 'yLabel', '');

    if ~exist('args', 'var') || isempty(args); args = struct(); end
    if ~isstruct(args)
        error(['esrFitRobust takes ONE struct: esrFitRobust(struct(''exp'', a)) or ' ...
               'esrFitRobust(struct(''f'', f, ''y'', y, ''e'', e)). The positional form ' ...
               '(Exp, opts) and the six outputs are gone -- which also retires ESRfit''s ' ...
               'unreadable signature for good.']);
    end
    opts = withDefaults(args, DEFAULTS, 'esrFitRobust');

    out = struct('ok', false, 'why', '', 'opts', opts);
    f0 = NaN; width = NaN; contrast = NaN; SNR = NaN; coeffError = nan(1,4);
    Exp = opts.exp;
    fromVectors = isempty(Exp);

    % ---- the data --------------------------------------------------------
    if fromVectors
        if isempty(opts.f) || isempty(opts.y)
            out.why = ['nothing to fit: give exp, or f and y. ' ...
                       'esrFitRobust(struct(''exp'', a)).'];
            fprintf('\nesrFitRobust FAILED: %s\n', out.why); return
        end
        f = opts.f(:); y = opts.y(:); e = opts.e(:);
        reps = 1; nAvg = 1;
    else
        [f, y, e, reps, nAvg, pullInfo] = esrPullData(Exp, ...
            struct('avgOutlierMad', opts.avgOutlierMad, 'verbose', opts.verbose, ...
                   'name', 'esrFitRobust'));
        dropped = pullInfo.dropped;
    end

    keep = isfinite(f) & isfinite(y);
    f = f(keep); y = y(keep);
    if isempty(e) || all(~isfinite(e)); e = []; else; e = e(keep); end

    % ---- the fit ---------------------------------------------------------
    fitOpts = struct('nLines', opts.nLines, 'maxLines', opts.maxLines, ...
                     'baselineOrder', opts.baselineOrder, 'snr', opts.snr, ...
                     'minSepMHz', opts.minSepMHz, 'minDepthSigma', opts.minDepthSigma, ...
                     'maxFwhmMHz', opts.maxFwhmMHz, 'seedCentres', opts.seedCentres, ...
                     'plot', false, 'name', 'esrFitRobust');
    out = esrMultiLineFit(f, y, e, fitOpts);
    out.f = f; out.y = y; out.e = e; out.repeats = reps; out.averages = nAvg;
    if ~fromVectors; out.dropped = dropped; end
    out.opts = opts;

    if ~out.ok
        fprintf('\nesrFitRobust FAILED: %s\n', out.why);
        out.f0 = f0; out.width = width; out.contrast = contrast; out.SNR = SNR;
        out.coeffError = coeffError; out.opts = opts; return
    end

    L = out.lines;
    [~, ord] = sort([L.centre]);
    L = L(ord); out.lines = L;
    f0       = [L.centre];
    width    = [L.fwhm];
    contrast = [L.depth];
    coeffError = [ [L.depthErr]; [L.centreErr]; [L.fwhmErr] ]';

    % Same definition as ESRfit's printed SNR, so the two are comparable: depth over the
    % mean error scaled by the total shots. It is a figure of merit, not a fitted quantity.
    if ~isempty(e)
        SNR = contrast ./ (mean(e)*sqrt(max(nAvg,1)*max(reps,1)));
    end
    % The old positional outputs, as named fields. Nothing that read them changes.
    out.f0 = f0; out.width = width; out.contrast = contrast; out.coeffError = coeffError;
    out.SNR = SNR;
    out.r2 = 1 - sum((y - out.model).^2)/max(sum((y - mean(y)).^2), eps);
    % The ceiling the noise imposes: the best R^2 any correct model could reach on this
    % data. R^2 = 0.82 against a ceiling of 0.84 is as good as the photons allow; the same
    % 0.82 against 0.99 is a real failure, and one number cannot tell those apart.
    if ~isempty(e)
        out.r2Ceiling = 1 - sum(e.^2)/max(sum((y - mean(y)).^2), eps);
    else
        out.r2Ceiling = NaN;
    end

    if opts.verbose; report(out); end
    if opts.plotIt;  drawIt(out);  end
end

% ======================================================================

function report(o)
    fprintf('\n--- esrFitRobust: %d line(s) ---\n', numel(o.lines));
    fprintf('%10s %14s %12s %12s %10s\n', 'line', 'centre/MHz', 'FWHM/MHz', 'depth/%', 'SNR');
    for i = 1:numel(o.lines)
        L = o.lines(i);
        sn = NaN; if numel(o.SNR) >= i; sn = o.SNR(i); end
        fprintf('%10d %7.3f+-%-5.3f %6.2f+-%-4.2f %6.3f+-%-4.3f %10.1f\n', ...
                i, L.centre, L.centreErr, L.fwhm, L.fwhmErr, ...
                100*L.depth, 100*L.depthErr, sn);
    end
    fprintf('R^2 %.5f (noise ceiling %.5f), reduced chi2 %.2f\n', ...
            o.r2, o.r2Ceiling, o.chi2Reduced);
    if isfinite(o.r2Ceiling) && o.r2 < o.r2Ceiling - 0.05
        fprintf(['  R^2 sits %.3f below what the noise allows, so this is a MODEL problem ' ...
                 'and not a\n  quiet spectrum. Try more lines (nLines), or a higher ' ...
                 'baselineOrder if the\n  background is curved.\n'], o.r2Ceiling - o.r2);
    end
    if numel(o.lines) > 1
        c = sort([o.lines.centre]);
        fprintf('splittings: %s MHz\n', numList(diff(c), '%.3f', ', '));
    end
end

function drawIt(o)
% ONE axes. The look comes from fitStyle, shared with rabiFitRobust and ramseyFitRobust --
% these three had drifted into three different blue/orange assignments, two grid styles and
% two y labels, and one file is the only way that stays fixed.
    S = fitStyle();
    fig = S.newFigure('ESR', o.opts.showResid);
    if o.opts.showResid
        ax = S.axes(subplot(4, 1, 1:3));
    else
        ax = S.axes(fig);
    end

    H = []; L = {};
    H(end+1) = S.data(ax, o.f, o.y, o.e, S.col.s1);  L{end+1} = 'measured';
    ff = linspace(min(o.f), max(o.f), 4000)';
    H(end+1) = S.model(ax, ff, o.modelAt(ff));       L{end+1} = 'fit';

    % A tick UNDER each fitted centre, not a line through the data. The line positions are
    % the result, and drawing them over the points they came from makes both harder to read.
    S.tight(ax, o.f); S.headroom(ax, 0.26);
    yl = ylim(ax);
    tickTop = yl(1) + 0.05*(yl(2) - yl(1));
    for i = 1:numel(o.lines)
        plot(ax, [o.lines(i).centre o.lines(i).centre], [yl(1) tickTop], '-', ...
             'LineWidth', 1.6, 'Color', S.col.flag);
    end
    ylim(ax, yl);

    xlabel(ax, 'frequency  [MHz]', 'Color', S.col.ink2);
    ylabel(ax, esrYLabel(o, S), 'Color', S.col.ink2);
    S.legend(ax, H, L);
    S.title(ax, esrTitle(o));

    if o.opts.showResid
        axR = S.axes(subplot(4, 1, 4));
        S.zeroLine(axR, [min(o.f) max(o.f)]);
        r = o.y - o.model;
        sg = o.e; if isempty(sg); sg = ones(size(r)); end
        S.resid(axR, o.f, r ./ max(sg, eps), S.col.s1);
        xlabel(axR, 'frequency  [MHz]', 'Color', S.col.ink2);
        ylabel(axR, 'residual / \sigma', 'Color', S.col.ink2);
        S.tight(axR, o.f);
        S.title(axR, sprintf('\\chi^2_\\nu = %.2f', o.chi2Reduced));
    end
end

function s = esrYLabel(o, S)
    s = S.yLabel;
    if isfield(o.opts, 'yLabel') && ~isempty(o.opts.yLabel); s = o.opts.yLabel; end
end

function T = esrTitle(o)
    if numel(o.lines) == 1
        L = o.lines(1);
        a = sprintf(['f_0 = %.3f \\pm %.3f MHz     FWHM = %.2f \\pm %.2f MHz     ' ...
                     'depth = %.3f \\pm %.3f %%'], ...
                    L.centre, L.centreErr, L.fwhm, L.fwhmErr, ...
                    100*L.depth, 100*L.depthErr);
    else
        a = sprintf('%d lines at %s MHz', numel(o.lines), ...
                    numList([o.lines.centre], '%.2f', ',  '));
    end
    b = sprintf('R^2 = %.4f of %.4f     |     \\chi^2_\\nu = %.2f     |     %d points', ...
                o.r2, o.r2Ceiling, o.chi2Reduced, numel(o.f));
    T = {a, b};
end

function s = numList(v, fmt, sep)
% strjoin(compose(...)) needs R2016b and compose does not exist in Octave, so the offline
% tests could not run the report. Same output, no version floor.
    parts = cell(1, numel(v));
    for i = 1:numel(v); parts{i} = sprintf(fmt, v(i)); end
    s = strjoin(parts, sep);
end

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
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
