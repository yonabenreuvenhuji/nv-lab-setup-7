classdef ExpPlotStyle
    %EXPPLOTSTYLE The look of the live plot in an experiment window (ESR, Rabi, Ramsey, ...).
    %
    %   ExpPlotStyle.apply(gAxes, nDim, xVector)
    %
    % It is the robust fits' look (Setup 7/functions/fitStyle.m) without the fit, so the plot
    % that fills in while a run is going and the fit figure drawn at the end read as the same
    % measurement: same palette, same markers, same grid, same labels.
    %
    % THE PALETTE COMES FROM fitStyle WHEN IT IS ON THE PATH (setup 7 adds it with
    % setup7AddPath). That keeps one source for the colours. Other setups do not have
    % fitStyle, so a copy of the same values is kept below as a fallback -- if you change the
    % palette in fitStyle, change it here too.
    %
    % It restyles the axes after Experiment.plotResults has drawn them, rather than drawing
    % anything itself. The experiment axes keep NextPlot = 'replacechildren', so every live
    % update throws the old lines away and draws plain new ones; apply() is called after each
    % update and puts the look back. AxesHelper is untouched, so the scan, counter and tracker
    % plots, which share it, keep their own look.

    properties (Constant)
        % Labels in the robust fits' wording. The y label is DIMENSIONLESS (signal photons over
        % reference photons, 1 = no contrast) -- see the y-label note in fitStyle for why it is
        % neither 'ratio' nor 'a.u.'.
        Y_NORMALISED = 'normalised fluorescence   (signal / reference)'
        Y_CONTRAST   = 'contrast   |0\rangle - |1\rangle   (normalised)'
        X_FREQUENCY  = 'frequency  [MHz]'
        X_RABI       = 'pulse duration \tau  [\mus]'
        X_RAMSEY     = 'free evolution \tau  [\mus]'
        X_ECHO       = 'total evolution time  2\tau\timescycles  [\mus]'
        X_T1         = '\tau  [\mus]'
    end

    methods (Static)
        function C = col()
            % The palette, as a struct: surface, ink, ink2, s1, s2, s3 (+ fitStyle's others).
            persistent cached
            if isempty(cached)
                cached = ExpPlotStyle.fallbackCol();
                try
                    if exist('fitStyle', 'file') == 2
                        S = fitStyle();
                        cached = S.col;
                    end
                catch
                    % fitStyle exists but failed: the copy below is the same palette.
                end
            end
            C = cached;
        end

        function apply(gAxes, nDim, xVector)
            % Restyle the experiment axes and everything drawn on them.
            if isempty(gAxes) || ~isgraphics(gAxes)
                return
            end
            if ~exist('xVector', 'var'); xVector = []; end
            C = ExpPlotStyle.col();

            % --- axes: secondary ink, ticks out, no box, solid 12% grid behind the data
            set(gAxes, 'FontSize', 9, 'XColor', C.ink2, 'YColor', C.ink2, ...
                'TickDir', 'out', 'Color', C.surface, 'Box', 'off');
            gAxes.XLabel.Color = C.ink2;
            gAxes.YLabel.Color = C.ink2;
            if nDim == 1
                grid(gAxes, 'on');
                set(gAxes, 'GridLineStyle', '-', 'GridAlpha', 0.12, 'Layer', 'bottom');
            else
                % An image fills the axes: a grid would only draw over it.
                grid(gAxes, 'off');
                set(gAxes, 'Layer', 'top');
                cb = ExpPlotStyle.findColorbar(gAxes);
                if ~isempty(cb)
                    set(cb, 'Color', C.ink2, 'FontSize', 9, 'TickDirection', 'out', 'Box', 'off');
                end
                return
            end

            % --- the series. Children are newest-first, so flip to drawing order: the first
            % signal gets orange (s1), the second purple (s2), the third slate (s3).
            kids = flipud(gAxes.Children(:));
            k = 0;
            for i = 1:numel(kids)
                h = kids(i);
                if any(strcmp(h.Type, {'errorbar', 'line'}))
                    k = k + 1;
                    ExpPlotStyle.styleSeries(h, k, C);
                end
            end

            % --- x held to the sweep, with a 2% pad. Auto-limits round outward to a tick and
            % leave an empty strip that reads as "the sweep stopped early".
            x = xVector(isfinite(xVector));
            if numel(x) >= 2 && max(x) > min(x)
                pad = 0.02*(max(x) - min(x));
                xlim(gAxes, [min(x) - pad, max(x) + pad]);
            end

            % --- legend: fixed corner, no box, with headroom reserved so it never sits on
            % the data ('best' is only as good as MATLAB's collision search).
            lg = gAxes.Legend;
            ylim(gAxes, 'auto');
            if ~isempty(lg) && isvalid(lg) && strcmp(lg.Visible, 'on')
                set(lg, 'Location', 'northeast', 'Box', 'off', 'TextColor', C.ink2, 'FontSize', 8);
                lg.String = ExpPlotStyle.ketify(lg.String);
                yl = ylim(gAxes);
                span = yl(2) - yl(1);
                if all(isfinite(yl)) && span > 0
                    ylim(gAxes, [yl(1), yl(2) + 0.28*span]);
                end
            end
        end

        function styleFigure(fig)
            % For the copies: 'Open In Figure' and the saved .png/.fig.
            if ~isempty(fig) && isgraphics(fig)
                set(fig, 'Color', ExpPlotStyle.col().surface);
            end
        end

        function styleTitle(gAxes, str)
            % fitStyle's title: normal weight, 9 pt, ink.
            title(gAxes, str, 'FontWeight', 'normal', 'FontSize', 9, ...
                'Color', ExpPlotStyle.col().ink, 'Interpreter', 'none');
        end
    end

    methods (Static, Access = private)
        function styleSeries(h, k, C)
            % fitStyle.data: filled markers with a white edge, a thin connecting line, error
            % bars at full opacity in the marker's own hue and no caps.
            slots = {C.s1, C.s2, C.s3};
            c = slots{mod(k - 1, numel(slots)) + 1};
            set(h, 'Color', c, 'LineWidth', 1.1, 'LineStyle', '-', 'Marker', 'o', ...
                'MarkerSize', 6, 'MarkerFaceColor', c, 'MarkerEdgeColor', C.surface);
            if k > numel(slots)
                % Only ESR's raw view with a mirrored sweep has four series. The 4th repeats a
                % hue, so it changes mark too: dashed, hollow.
                set(h, 'LineStyle', '--', 'MarkerFaceColor', C.surface, 'MarkerEdgeColor', c);
            end
            if isprop(h, 'CapSize')
                h.CapSize = 0;
            end
        end

        function cb = findColorbar(gAxes)
            cb = [];
            try
                cb = gAxes.Colorbar;    % R2018b+
            catch
            end
        end

        function s = ketify(s)
            % '|1>' -> '|1\rangle' so the legend shows a ket, not a "greater than".
            s = regexprep(s, '\|(\w+)>', '|$1\\rangle');
        end

        function C = fallbackCol()
            % A copy of fitStyle's palette, for setups without Setup 7/functions on the path.
            C = struct( ...
                'surface', [252 252 251]/255, ...
                'ink',     [ 20  20  20]/255, ...
                'ink2',    [ 82  81  78]/255, ...
                's1',      [235 104  52]/255, ...   % orange, default
                's2',      [122  61 148]/255, ...   % purple
                's3',      [107 120 133]/255, ...   % slate
                'fit',     [ 20  20  20]/255, ...
                'env',     [107 120 133]/255, ...
                'flag',    [ 20  20  20]/255, ...
                'band',    [219 224 232]/255);
        end
    end
end
