classdef ExpFit
    %EXPFIT The robust fits, in the experiment window.
    %
    % An experiment opts in by returning a spec from fitSpec:
    %
    %     function spec = fitSpec(obj)
    %         spec = ExpFit.esr(obj);
    %     end
    %
    % and then a.fit, the window's Fit button and Auto-fit all come here. Each is ONE call to
    % the experiment's robust fitter (esrFitRobust, rabiFitRobust, ...), made exactly as a
    % script makes it, and the window draws that call's curve and title. Open In Figure makes
    % the same call with plotIt = true, so the figure is the fitter's own.
    %
    % WHAT IS FITTED follows the view the window shows:
    %   normal      the fitter's own data path, struct('exp', a): error bars rebuilt from the
    %               averages, both projections where there are two. The same numbers as
    %               o = esrFitRobust(struct('exp', a)) in a script.
    %   alternate,  every displayed trace on its own, as vectors (f or t, y, e): ESR's With MW
    %   raw         and Without MW, the raw detection windows, Ramsey's |0> - |1>. The fitters
    %               expect a ratio near 1 (their contrast is amplitude over baseline), so raw
    %               counts are divided by their median, and a difference centred on 0 is
    %               fitted as 1 + difference; the curve is mapped back for the window. A
    %               reference window that should be flat comes back "no line" or with a
    %               contrast consistent with zero -- which is the point of fitting it.
    %
    % THE SPEC (what fitSpec returns):
    %   fitter      the fitter's name
    %   xField      what the fitter calls the x vector: 'f' (ESR) or 't'
    %   controls    the window's option widgets, from ExpFit.popup / check / edit
    %
    % THE RECORD (a.fits.normal / .alternate / .raw; saved with the experiment):
    %   view, fitter, options, time, averages (how many the fit saw),
    %   traces(k): name, ok, why, x, y (the curve; one row per projection), scale,
    %              titleLines (the fitter's figure title), out (the fitter's out, without
    %              function handles or the experiment object, so it saves cleanly)

    properties (Constant)
        TAG = 'ExpFitOverlay'       % the fit curves on the window's axes
        VIEWS = {'normal', 'alternate', 'raw'}
    end

    %% Specs: one per kind of experiment. Options are the fitter's own (see its help).
    methods (Static)
        function s = esr(~)
            s = ExpFit.newSpec('esrFitRobust', 'f');
            s.controls = [ ...
                ExpFit.popup('nLines', 'NV lines', {'auto', '1', '2', '3', '4', '5', '6', '7', '8'}, ...
                    {[], 1, 2, 3, 4, 5, 6, 7, 8}, ...
                    'Force the number of NV lines. The resonator''s dip is not counted.'), ...
                ExpFit.popup('hyperfine', 'Hyperfine', {'auto', 'test', 'skip'}, {[], true, false}, ...
                    ['14N / 15N multiplet. auto: tested on one line with a step of 0.6 MHz or less ' ...
                     '(pulsed sweeps).']), ...
                ExpFit.popup('resonance', 'Resonator dip', {'auto', 'look for it', 'none'}, {[], true, false}, ...
                    'The resonator''s broad dip, fitted as its own component and kept out of the field.'), ...
                ExpFit.popup('field', 'Field solve', {'auto', 'off'}, {[], false}, ...
                    '|B|, the tilt from the aligned axis and D, from the lines.'), ...
                ExpFit.popup('baselineOrder', 'Baseline', {'auto', 'order 0', 'order 1', 'order 2', 'order 3'}, ...
                    {[], 0, 1, 2, 3}, 'Polynomial baseline. auto: 2 for spans of 60 MHz or more, else 1.'), ...
                ExpFit.edit('range', 'Range (MHz)', 'fmin fmax: fit only this part of the sweep. Empty: all of it.'), ...
                ExpFit.check('showResid', 'Residuals panel', false, 'In the fit figure (Open In Figure).', 'figure')];
        end

        function s = rabi(~)
            s = ExpFit.newSpec('rabiFitRobust', 't');
            s.controls = [ ...
                ExpFit.popup('nFreq', 'Frequencies', {'1', '2'}, {1, 2}, ...
                    ['2 when the residual autocorrelation is over 0.4: a second orientation driven ' ...
                     'off resonance, or the hyperfine.']), ...
                ExpFit.popup('q', 'Decay', {'search', 'exp (q = 1)', 'gauss (q = 2)'}, {[], 1, 2}, ...
                    'Decay exponent q of exp(-(t/T)^q). search: tries 1 and 2.'), ...
                ExpFit.check('linearDrift', 'Linear drift', true, ...
                    'A sloped baseline, so slow laser drift is not taken for a decay.'), ...
                ExpFit.check('robust', 'Outlier pass', true, 'One 4-sigma outlier pass after the fit, then refit.'), ...
                ExpFit.check('showResid', 'Residuals panel', false, 'In the fit figure (Open In Figure).', 'figure'), ...
                ExpFit.check('showScan', 'Frequency scan', false, ...
                    'The periodogram, as a second figure (Open In Figure).', 'figure')];
        end

        function s = ramsey(~)
            s = ExpFit.newSpec('ramseyFitRobust', 't');
            s.controls = [ ...
                ExpFit.popup('nHyperfine', 'Hyperfine lines', {'3 (14N)', '2 (15N)', '1 (none)'}, {3, 2, 1}, ...
                    '3 for 14N, 2 for 15N, 1 for a plain damped cosine.'), ...
                ExpFit.projection(), ...
                ExpFit.check('autoHyperfine', 'Free A if it wins', true, ...
                    'Also fit the hyperfine splitting A, and keep it only if it wins by 3 sigma.'), ...
                ExpFit.check('fitHyperfine', 'Always fit A', false, 'Free the hyperfine splitting A unconditionally.'), ...
                ExpFit.check('robust', 'Robust (Tukey)', true, 'Tukey bisquare re-weighting after the first pass.'), ...
                ExpFit.edit('tMin', 't min (us)', 'Drop taus below this. Empty: 2 x halfPiTime from the experiment.'), ...
                ExpFit.check('showResid', 'Residuals panel', false, 'In the fit figure (Open In Figure).', 'figure'), ...
                ExpFit.check('showSpectrum', 'Fringe spectrum', false, 'In the fit figure (Open In Figure).', 'figure')];
        end

        function s = echo(~)
            s = ExpFit.newSpec('echoFitRobust', 't');
            s.controls = [ ...
                ExpFit.popup('revivals', '13C revivals', {'auto', 'on', 'off'}, {[], 'on', 'off'}, ...
                    'The 13C revival comb. auto: kept only when the data asks for it.'), ...
                ExpFit.popup('pFixed', 'Stretch p', {'fit', '1', '2', '3'}, {[], 1, 2, 3}, ...
                    'Fix the stretch exponent (3: a dipolar 13C bath). fit: free, 0.5 to 4.'), ...
                ExpFit.projection(), ...
                ExpFit.check('fitDrift', 'Linear drift', false, ...
                    'Free the linear drift. It trades off against T2: read echoFitRobust''s note first.'), ...
                ExpFit.check('robust', 'Robust (Tukey)', true, 'Tukey bisquare re-weighting after the first pass.'), ...
                ExpFit.edit('tMin', 't min (us)', 'Drop points below this t (the window''s x axis). Empty: from the pulse times.'), ...
                ExpFit.edit('tMax', 't max (us)', 'Drop points above this t. Empty: keep them all.'), ...
                ExpFit.check('showResid', 'Residuals panel', false, 'In the fit figure (Open In Figure).', 'figure')];
        end

        function s = t1(~)
            s = ExpFit.newSpec('t1FitRobust', 't');
            s.controls = [ ...
                ExpFit.popup('model', 'Model', {'stretched', 'double', 'auto'}, {'stretched', 'double', 'auto'}, ...
                    'One stretched exponential, two exponentials, or let the data choose.'), ...
                ExpFit.popup('pFixed', 'Stretch p', {'fit', '1 (exp)'}, {[], 1}, ...
                    'Fix the stretch exponent. 1: a plain exponential.'), ...
                ExpFit.projection(), ...
                ExpFit.check('robust', 'Robust (Tukey)', true, 'Tukey bisquare re-weighting after the first pass.'), ...
                ExpFit.edit('tMin', 'tau min (us)', 'Drop taus below this. Empty: keep them all.'), ...
                ExpFit.edit('tMax', 'tau max (us)', 'Drop taus above this. Empty: keep them all.'), ...
                ExpFit.check('logX', 'Log tau axis', true, 'In the fit figure (Open In Figure).', 'figure'), ...
                ExpFit.check('showResid', 'Residuals panel', false, 'In the fit figure (Open In Figure).', 'figure')];
        end
    end

    %% Building specs
    methods (Static)
        function s = newSpec(fitter, xField)
            s = struct('fitter', fitter, 'xField', xField, 'controls', ExpFit.emptyControls());
        end

        function c = popup(key, label, items, values, tip, scope)
            % A drop-down. values{k} is what items{k} passes to the fitter; [] passes nothing,
            % so the fitter's own default applies.
            if ~exist('scope', 'var'); scope = 'all'; end
            c = ExpFit.control(key, label, 'popup', items, values, [], tip, scope);
        end

        function c = check(key, label, default, tip, scope)
            % A check box. default is the fitter's own default, shown until it is changed.
            if ~exist('scope', 'var'); scope = 'all'; end
            c = ExpFit.control(key, label, 'check', {}, {}, logical(default), tip, scope);
        end

        function c = edit(key, label, tip, scope)
            % Numbers, separated by spaces or commas. Empty passes nothing.
            if ~exist('scope', 'var'); scope = 'all'; end
            c = ExpFit.control(key, label, 'edit', {}, {}, [], tip, scope);
        end

        function c = projection()
            % useProjection: only meaningful on the experiment's own data (the normal view),
            % where the fitter reads both projections.
            c = ExpFit.popup('useProjection', 'Projections', {'both', '|1> only', '|0> only'}, ...
                {'both', 'one', 'zero'}, 'Which projections the fit uses (the normal view only).', 'normal');
        end
    end

    methods (Static, Access = private)
        function c = control(key, label, style, items, values, default, tip, scope)
            c = struct('key', key, 'label', label, 'style', style, 'items', {items}, ...
                'values', {values}, 'default', default, 'tip', tip, 'scope', scope);
        end

        function c = emptyControls()
            c = ExpFit.control('', '', '', {}, {}, [], '', '');
            c = c([]);
        end
    end

    %% Running a fit
    methods (Static)
        function tf = isAvailable(spec)
            % The fitters live in Setup 7\functions (setup7AddPath puts them on the path).
            tf = ~isempty(spec) && exist(spec.fitter, 'file') == 2;
        end

        function [rec, outs] = run(exp, spec, view, userOpts, plotIt, verbose)
            % One fit of one view. outs{k} is the fitter's out for trace k, untouched.
            if ~ExpFit.isAvailable(spec)
                error('ExpFit:noFitter', '%s is not on the path. Run setup7AddPath first.', spec.fitter);
            end
            o = userOpts;
            if isempty(o); o = struct(); end
            o.plotIt = logical(plotIt);
            o.verbose = logical(verbose);
            fn = spec.fitter;

            if strcmp(view, 'normal')
                o.exp = exp;
                outs = {feval(fn, o)};
                names = {''};
                scales = 1;
                passed = rmfield(o, 'exp');
            else
                o = ExpFit.dropScope(o, spec, 'normal');
                passed = o;
                [x, tr] = exp.viewTraces(view);
                if isempty(tr)
                    error('ExpFit:noData', 'Nothing to fit in this view yet.');
                end
                n = numel(tr);
                outs = cell(1, n); names = {tr.name}; scales = ones(1, n); shifts = zeros(1, n);
                for k = 1:n
                    [ok, scales(k), shifts(k)] = ExpFit.vectorArgs(o, spec, x, tr(k));
                    try
                        outs{k} = feval(fn, ok);
                    catch err
                        outs{k} = struct('ok', false, 'why', err.message);
                    end
                end
            end
            if ~exist('shifts', 'var'); shifts = 0; end
            rec = ExpFit.record(exp, spec, view, passed, outs, names, scales, shifts);
        end

        function rec = record(exp, spec, view, passed, outs, names, scales, shifts)
            % What the window draws and the file keeps.
            passed = rmfield(passed, intersect(fieldnames(passed), {'plotIt', 'verbose'}));
            rec = struct('view', view, 'fitter', spec.fitter, ...
                'options', ExpFit.sanitize(passed), ...
                'time', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...    %#ok<TNOW1,DATST>
                'averages', exp.currIter, 'traces', []);
            tr = struct('name', {}, 'ok', {}, 'why', {}, 'x', {}, 'y', {}, 'scale', {}, 'shift', {}, ...
                'titleLines', {}, 'out', {});
            for k = 1:numel(outs)
                out = outs{k};
                t = struct('name', names{k}, 'ok', isstruct(out) && isfield(out, 'ok') && isequal(out.ok, true), ...
                    'why', '', 'x', [], 'y', [], 'scale', scales(k), 'shift', shifts(k), ...
                    'titleLines', {{}}, 'out', []);
                if isfield(out, 'why') && ischar(out.why); t.why = out.why; end
                if t.ok
                    try
                        [t.x, t.y] = ExpFit.curveOf(out, spec);
                        t.y = (t.y - shifts(k)) * scales(k);    % back to what the window plots
                    catch err
                        t.why = sprintf('fitted, but the curve could not be drawn: %s', err.message);
                    end
                    if isfield(out, 'titleLines') && iscell(out.titleLines)
                        t.titleLines = out.titleLines;
                    end
                end
                t.out = ExpFit.sanitize(out);
                tr(k) = t;
            end
            rec.traces = tr;
        end

        function [x, y] = curveOf(out, spec)
            % The fitted curve over the span the fit used, dense enough to look smooth.
            xd = out.(spec.xField);
            xd = xd(isfinite(xd));
            lo = min(xd); hi = max(xd);
            if lo > 0 && hi/lo > 20
                x = logspace(log10(lo), log10(hi), 1500);   % log-spaced sweeps (T1)
            else
                x = linspace(lo, hi, 1500);
            end
            if isfield(out, 'curveAt') && isa(out.curveAt, 'function_handle')
                y = out.curveAt(x);
            elseif isfield(out, 'modelAt') && isa(out.modelAt, 'function_handle')
                y = reshape(out.modelAt(x(:)), 1, []);
            else
                % A fitter without a curve function: its model at the data points.
                x = reshape(out.(spec.xField), 1, []);
                y = out.model;
                if size(y, 2) ~= numel(x); y = y.'; end
            end
        end
    end

    %% Drawing
    methods (Static)
        function drawOverlay(gAxes, rec, nowAverages)
            % The fit curves and the fit's title on the experiment axes. A fit that saw fewer
            % averages than the plot now shows is drawn dashed and says so.
            delete(findobj(gAxes, 'Tag', ExpFit.TAG));
            if isempty(rec)
                title(gAxes, '');
                return
            end
            C = ExpPlotStyle.col();
            stale = rec.averages ~= nowAverages;
            for k = 1:numel(rec.traces)
                t = rec.traces(k);
                if ~t.ok || isempty(t.y); continue; end
                for r = 1:size(t.y, 1)
                    h = line('Parent', gAxes, 'XData', t.x, 'YData', t.y(r, :), ...
                        'Color', C.fit, 'LineWidth', 2.2, 'LineStyle', '-', ...
                        'Tag', ExpFit.TAG, 'PickableParts', 'none', 'HitTest', 'off');
                    if stale
                        set(h, 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 0.55*C.fit + 0.45*C.surface);
                    end
                    try
                        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    catch
                    end
                end
            end
            lines = ExpFit.titleOf(rec);
            if stale
                lines{end+1} = sprintf('{\\color[rgb]{0.55 0.35 0.05}fit of %d average%s, the plot now shows %d -- refit}', ...
                    rec.averages, ExpFit.plural(rec.averages), nowAverages);
            end
            title(gAxes, lines, 'FontWeight', 'normal', 'FontSize', 8.5, 'Color', C.ink, ...
                'Interpreter', 'tex');
        end

        function lines = titleOf(rec)
            % The fitter's own title. With several traces (alternate and raw views), the first
            % line of each, named.
            lines = {};
            tr = rec.traces;
            if numel(tr) == 1
                t = tr(1);
                if t.ok && ~isempty(t.titleLines)
                    lines = reshape(t.titleLines, 1, []);
                elseif ~t.ok
                    lines = {['fit failed: ' ExpFit.texEscape(t.why)]};
                end
                return
            end
            for k = 1:numel(tr)
                name = ExpFit.texEscape(tr(k).name);
                if tr(k).ok && ~isempty(tr(k).titleLines)
                    lines{end+1} = sprintf('%s:   %s', name, tr(k).titleLines{1}); %#ok<AGROW>
                elseif tr(k).ok
                    lines{end+1} = sprintf('%s:   fitted', name); %#ok<AGROW>
                else
                    lines{end+1} = sprintf('%s:   no fit (%s)', name, ExpFit.texEscape(ExpFit.shorten(tr(k).why, 70))); %#ok<AGROW>
                end
            end
        end

        function txt = summary(rec, nowAverages)
            % The fit as plain text, for the window's readout and the clipboard.
            if isempty(rec); txt = {}; return; end
            txt = {sprintf('%s, %s view, %d average%s, %s', rec.fitter, rec.view, rec.averages, ...
                ExpFit.plural(rec.averages), rec.time(12:end))};
            if exist('nowAverages', 'var') && rec.averages ~= nowAverages
                txt{end+1} = sprintf('STALE: the plot now shows %d averages. Refit.', nowAverages);
            end
            for k = 1:numel(rec.traces)
                t = rec.traces(k);
                if ~isempty(t.name); txt{end+1} = ''; txt{end+1} = [t.name ':']; end %#ok<AGROW>
                if t.ok
                    for j = 1:numel(t.titleLines)
                        txt = [txt, ExpFit.plainLines(t.titleLines{j})]; %#ok<AGROW>
                    end
                    if isempty(t.titleLines); txt{end+1} = 'fitted'; end %#ok<AGROW>
                else
                    txt{end+1} = ['no fit: ' t.why]; %#ok<AGROW>
                end
            end
        end
    end

    %% Helpers
    methods (Static)
        function v = sanitize(v, depth)
            % Strip what must not go into a .mat: function handles (they would carry their
            % workspace, the experiment included) and handle objects (the experiment itself,
            % in opts.exp). Values and structs are kept.
            if ~exist('depth', 'var'); depth = 0; end
            if depth > 20; v = []; return; end
            if isa(v, 'function_handle')
                v = [];
            elseif isstruct(v)
                f = fieldnames(v);
                for i = 1:numel(v)
                    for j = 1:numel(f)
                        v(i).(f{j}) = ExpFit.sanitize(v(i).(f{j}), depth + 1);
                    end
                end
            elseif iscell(v)
                for i = 1:numel(v)
                    v{i} = ExpFit.sanitize(v{i}, depth + 1);
                end
            elseif isobject(v) && ~(isstring(v) || isdatetime(v) || isduration(v) || iscategorical(v))
                v = sprintf('<%s>', class(v));
            end
        end

        function o = dropScope(o, spec, scope)
            % Remove the options a control marks as belonging to another view.
            for c = spec.controls
                if strcmp(c.scope, scope) && isfield(o, c.key)
                    o = rmfield(o, c.key);
                end
            end
        end

        function s = plainLines(texStr)
            % A tex title line as plain text, one quantity per line: the titles separate
            % quantities by '|' or by a run of spaces.
            s = ExpFit.texToPlain(texStr);
            s = strtrim(regexp(s, '(\s*\|\s*)|(\s{3,})', 'split'));
            s = s(~cellfun(@isempty, s));
            s = regexprep(s, '\s+', ' ');
        end

        function s = texToPlain(s)
            pairs = { ...
                '\pm', char(177);  '\mus', [char(181) 's'];  '\mu', char(181);  '^\circ', char(176); ...
                '\chi^2_\nu', 'chi2_nu';  '\tau', 'tau';  '\pi', 'pi';  '\Delta', 'Delta'; ...
                '\delta', 'delta';  '\gamma', 'gamma';  '\Omega', 'Omega';  '\langle', '<'; ...
                '\rangle', '>';  '_{||}', '_par';  'T_2^*', 'T2*';  'T_2', 'T2';  'T_1', 'T1'; ...
                'R^2', 'R2';  '\times', 'x';  '\color', ''};
            for i = 1:size(pairs, 1)
                s = strrep(s, pairs{i, 1}, pairs{i, 2});
            end
            s = regexprep(s, '\[rgb\]\{[^}]*\}', '');
            s = regexprep(s, '\\([a-zA-Z]+)', '$1');
            s = strrep(strrep(s, '{', ''), '}', '');
        end

        function s = texEscape(s)
            % Text from outside (trace names, failure messages) shown under the tex interpreter.
            s = regexprep(s, '([_^\\{}])', '\\$1');
        end

        function s = shorten(s, n)
            if numel(s) > n; s = [s(1:n-3) '...']; end
        end

        function s = plural(n)
            if n == 1; s = ''; else; s = 's'; end
        end
    end

    methods (Static, Access = private)
        function [o, scale, shift] = vectorArgs(o, spec, x, tr)
            % One trace as the fitter's vector input, which the fitters expect to be a ratio
            % near 1 (their contrast is the amplitude over the baseline). Raw counts are
            % divided by their median. A difference trace centred on 0 (|0> - |1>) is fitted
            % as 1 + trace: on 0 the contrast is amplitude/~0 and comes out as nonsense.
            x = reshape(x, 1, []); y = reshape(tr.y, 1, []); e = reshape(tr.e, 1, []);
            if numel(x) ~= numel(y)
                error('ExpFit:size', '%s: %d x values for %d points.', tr.name, numel(x), numel(y));
            end
            keep = isfinite(x) & isfinite(y);
            x = x(keep); y = y(keep);
            if numel(e) == numel(keep); e = e(keep); else; e = []; end
            scale = 1;
            shift = 0;
            if tr.raw
                s = median(abs(y));
                if isfinite(s) && s > 0; scale = s; end
            elseif median(abs(y)) < 0.5
                shift = 1;
            end
            y = y / scale + shift;
            e = e / scale;
            % The y label of the fitter's figure, if it draws one (Open In Figure)
            if scale ~= 1
                o.yLabel = sprintf('%s / %.4g', tr.label, scale);
            elseif shift ~= 0
                o.yLabel = sprintf('1 + %s', tr.label);
            elseif ~isempty(tr.label)
                o.yLabel = tr.label;
            end
            if isempty(e) || ~all(isfinite(e)) || any(e <= 0)
                % No usable error bars (one average): the point-to-point scatter, the same for
                % every point. The fitters Birge-scale their errors anyway.
                e = max(std(diff(y))/sqrt(2), eps) * ones(size(y));
            end
            o.(spec.xField) = x;
            o.y = y;
            o.e = e;
        end
    end
end
