classdef MainGuiStyle
    %MAINGUISTYLE The look of the main setup window, and the factory for its controls.
    %
    %   S = MainGuiStyle;                       % (all members are static)
    %   lbl = MainGuiStyle.label(parent, 'From (µm)');
    %   edt = MainGuiStyle.edit(parent, '12.5');
    %   btn = MainGuiStyle.button(parent, 'Scan', 'primary');
    %
    % One palette for the whole window, taken from the fit figures (fitStyle / ExpPlotStyle),
    % so the live photodiode trace, the scan image and the experiment plots read as one
    % system: blue for the first series and for the main action, amber for the second series
    % and for anything "pending" (typed but not applied yet), brick red for stop/halt.
    %
    % Only the main window uses this. Experiment and Helmholtz windows keep GuiComponent's
    % PROP_* look.
    %
    % SIZES ARE EXPLICIT. Every row of every section has a known height (ROW_H), and a
    % section reserves TITLE_H for its title. The old panels computed their heights from
    % guesses that were too small once R2025a enlarged uicontrol fonts, so the GUI Layout
    % Toolbox squeezed the content under the panel titles (the hidden "From / To / # Pts"
    % header row). Here a section is always at least as tall as its rows.

    properties (Constant)
        % --- sizes (pixels)
        ROW_H = 26          % one row of controls
        HEAD_H = 18         % a column-header row
        TITLE_H = 22        % a section's title row, with the space under it
        PAD = 6             % padding inside a section
        GAP = 6             % spacing between rows / sections
        FONT = 9            % labels, buttons
        FONT_EDIT = 10      % numbers being typed
        FONT_BIG = 18       % live readout
    end

    methods (Static)
        function C = col()
            % Palette. Data colours from fitStyle (via ExpPlotStyle, which keeps a copy for
            % setups without Setup 7 on the path); interface colours below.
            persistent cached
            if isempty(cached)
                D = ExpPlotStyle.col();
                cached = struct( ...
                    'window',   [0.937 0.937 0.930], ...    % figure background
                    'panel',    [0.980 0.980 0.975], ...    % section background
                    'field',    [1 1 1], ...                % edit fields
                    'line',     [0.80 0.80 0.78], ...       % section borders
                    'ink',      D.ink, ...
                    'ink2',     D.ink2, ...
                    'ink3',     [0.50 0.50 0.48], ...       % hints, units, disabled text
                    'readonly', [0.925 0.925 0.915], ...    % read-only values
                    'off',      [0.90 0.90 0.89], ...       % disabled fields
                    'accent',   D.s1, ...                   % main action, series 1
                    'accent2',  D.s2, ...                   % series 2
                    'pending',  [1.00 0.93 0.76], ...       % typed, not yet applied (pale amber)
                    'pendingInk', [0.55 0.33 0.00], ...
                    'stop',     D.fit, ...                  % stop / halt (brick red)
                    'ok',       [0.20 0.55 0.30], ...       % saved / done
                    'warn',     [0.72 0.45 0.00], ...
                    'err',      [0.70 0.12 0.10], ...
                    'surface',  D.surface, ...
                    'band',     [0.86 0.88 0.91]);
            end
            C = cached;
        end

        %% ---- containers
        function [panel, inner, titleText] = section(parent, title, varargin)
            % A titled section: a bordered panel, its title as the first row, and an inner
            % VBox for the content. Full height: MainGuiStyle.sectionHeight(contentHeight).
            %
            % The title is a text row, NOT the panel's Title: a uix.Panel with a Title lays
            % its content out for the panel's inner height measured before R2025a+ knows the
            % title's height, so the content ends up under the title -- the "labels cover
            % the content" bug. An untitled panel has no title height to get wrong.
            C = MainGuiStyle.col();
            panel = uix.Panel('Parent', parent, ...
                'Padding', MainGuiStyle.PAD, 'BorderType', 'line', ...
                'BackgroundColor', C.panel, varargin{:});
            MainGuiStyle.borderColor(panel, C.line);
            outer = uix.VBox('Parent', panel, 'Spacing', MainGuiStyle.TITLE_H - 16, ...
                'BackgroundColor', C.panel);
            titleText = uicontrol('Parent', outer, 'Style', 'text', 'String', title, ...
                'FontSize', MainGuiStyle.FONT + 1, 'FontWeight', 'bold', ...
                'ForegroundColor', C.ink, 'BackgroundColor', C.panel, ...
                'HorizontalAlignment', 'left', 'Tag', 'sectionTitle');
            inner = uix.VBox('Parent', outer, 'Spacing', MainGuiStyle.GAP, ...
                'BackgroundColor', C.panel);
            outer.Heights = [16 -1];
        end

        function borderColor(panel, c)
            % BorderColor (R2023a+); HighlightColor before. (ShadowColor is deprecated.)
            if isprop(panel, 'BorderColor')
                set(panel, 'BorderColor', c);
            else
                set(panel, 'HighlightColor', c);
            end
        end

        function makeInteractive(ax)
            % Zoom, pan and data tips on a plot: scroll = zoom, drag = pan, hover / click =
            % data tip, double-click = back to the full view, and a toolbar (top right, on
            % hover) with the same plus export.
            %
            % Needs NO figure-wide mode to be on: the old "Marker" cursor switched the
            % window's data-cursor mode on, and while any mode is on MATLAB turns zoom and
            % pan off in every plot of the window.
            if isempty(ax) || ~isgraphics(ax); return; end
            try
                axtoolbar(ax, {'export', 'datacursor', 'pan', 'zoomin', 'zoomout', 'restoreview'});
                enableDefaultInteractivity(ax);
            catch
                % (older releases: no axes toolbar)
            end
        end

        function relayout(fig)
            % Lay out every GUI Layout Toolbox panel of the window again, once it is shown.
            %
            % THIS IS THE "TITLES COVER THE CONTENT" BUG. uix.Panel sizes its contents to the
            % panel's inner height (the panel minus its title), measured when it first draws.
            % On R2025a+ that happens before the title's height is known, so the contents are
            % made too tall, and they stick out at the top, under the title. A panel only
            % re-measures when its own size changes, so panels whose size never changes
            % (most of them) kept the hidden header rows for good. Nudging the padding makes
            % each panel lay out again, now with the right title height.
            if isempty(fig) || ~isgraphics(fig); return; end
            drawnow;
            panels = findall(fig, '-isa', 'uix.Panel');
            for k = 1:numel(panels)
                p = panels(k).Padding;
                panels(k).Padding = p + 1;
                panels(k).Padding = p;
            end
            drawnow;
        end

        function h = sectionHeight(contentHeight)
            % Total height of a section whose rows (with their spacing) add up to contentHeight.
            h = contentHeight + MainGuiStyle.TITLE_H + 2*MainGuiStyle.PAD + 2;
        end

        function h = rowsHeight(heights)
            % Height of a VBox with these row heights and GAP spacing.
            h = sum(heights) + MainGuiStyle.GAP * max(numel(heights) - 1, 0);
        end

        function box = hbox(parent, spacing)
            if nargin < 2; spacing = MainGuiStyle.GAP; end
            box = uix.HBox('Parent', parent, 'Spacing', spacing, ...
                'BackgroundColor', MainGuiStyle.col().panel);
        end

        function box = vbox(parent, spacing)
            if nargin < 2; spacing = MainGuiStyle.GAP; end
            box = uix.VBox('Parent', parent, 'Spacing', spacing, ...
                'BackgroundColor', MainGuiStyle.col().panel);
        end

        function g = grid(parent, spacing)
            if nargin < 2; spacing = 4; end
            g = uix.Grid('Parent', parent, 'Spacing', spacing, ...
                'BackgroundColor', MainGuiStyle.col().panel);
        end

        function e = empty(parent)
            e = uix.Empty('Parent', parent);
            try; e.BackgroundColor = MainGuiStyle.col().panel; catch; end
        end

        %% ---- text
        function h = label(parent, str, align, varargin)
            % Plain label: secondary ink on the panel background, centred vertically in its
            % row (a text uicontrol draws at the top of its box, so on its own it sits higher
            % than the field next to it).
            if nargin < 3 || isempty(align); align = 'left'; end
            C = MainGuiStyle.col();
            holder = MainGuiStyle.centred(parent, 16);
            h = uicontrol('Parent', holder, 'Style', 'text', 'String', str, ...
                'FontSize', MainGuiStyle.FONT, 'ForegroundColor', C.ink2, ...
                'BackgroundColor', C.panel, 'HorizontalAlignment', align, varargin{:});
            MainGuiStyle.finishCentred(holder, h);
        end

        function holder = centred(parent, textHeight)
            % A box that keeps one text control in the middle of the row's height
            holder = uix.VBox('Parent', parent, 'BackgroundColor', MainGuiStyle.col().panel);
            uix.Empty('Parent', holder);
            setappdata(holder, 'textHeight', textHeight);
        end

        function finishCentred(holder, h)
            uix.Empty('Parent', holder);
            th = getappdata(holder, 'textHeight');
            fs = get(h, 'FontSize');
            th = max(th, ceil(fs * 1.75));
            holder.Heights = [-1 th -1];
            holder.MinimumHeights = [0 0 0];
            try; holder.BackgroundColor = get(h, 'BackgroundColor'); catch; end
        end

        function h = header(parent, str, align, varargin)
            % Column header: small, bold, secondary ink.
            if nargin < 3 || isempty(align); align = 'center'; end
            C = MainGuiStyle.col();
            h = uicontrol('Parent', parent, 'Style', 'text', 'String', str, ...
                'FontSize', MainGuiStyle.FONT - 1, 'FontWeight', 'bold', ...
                'ForegroundColor', C.ink3, 'BackgroundColor', C.panel, ...
                'HorizontalAlignment', align, varargin{:});
        end

        function h = hint(parent, str, align, varargin)
            % A one-line explanation or summary, in tertiary ink.
            if nargin < 3 || isempty(align); align = 'left'; end
            C = MainGuiStyle.col();
            holder = MainGuiStyle.centred(parent, 15);
            h = uicontrol('Parent', holder, 'Style', 'text', 'String', str, ...
                'FontSize', MainGuiStyle.FONT - 1, 'ForegroundColor', C.ink3, ...
                'BackgroundColor', C.panel, 'HorizontalAlignment', align, varargin{:});
            MainGuiStyle.finishCentred(holder, h);
        end

        function h = axisTag(parent, letter)
            % The X / Y / Z at the start of a row.
            C = MainGuiStyle.col();
            holder = MainGuiStyle.centred(parent, 19);
            h = uicontrol('Parent', holder, 'Style', 'text', 'String', upper(letter), ...
                'FontSize', MainGuiStyle.FONT + 2, 'FontWeight', 'bold', ...
                'ForegroundColor', C.ink, 'BackgroundColor', C.panel, ...
                'HorizontalAlignment', 'center');
            MainGuiStyle.finishCentred(holder, h);
        end

        %% ---- inputs
        function h = edit(parent, str, varargin)
            % A number the user types. Right-aligned, like numbers in a table.
            if nargin < 2; str = ''; end
            C = MainGuiStyle.col();
            h = uicontrol('Parent', parent, 'Style', 'edit', 'String', str, ...
                'FontSize', MainGuiStyle.FONT_EDIT, 'ForegroundColor', C.ink, ...
                'BackgroundColor', C.field, 'HorizontalAlignment', 'right', varargin{:});
        end

        function h = readout(parent, str, varargin)
            % A value shown, not typed (the stage position, a computed step). Looks like a
            % field, but grey and not editable.
            if nargin < 2; str = ''; end
            C = MainGuiStyle.col();
            h = uicontrol('Parent', parent, 'Style', 'edit', 'String', str, ...
                'Enable', 'inactive', 'FontSize', MainGuiStyle.FONT_EDIT, ...
                'ForegroundColor', C.ink, 'BackgroundColor', C.readonly, ...
                'HorizontalAlignment', 'right', varargin{:});
        end

        function h = checkbox(parent, str, value, varargin)
            if nargin < 3; value = false; end
            C = MainGuiStyle.col();
            h = uicontrol('Parent', parent, 'Style', 'checkbox', 'String', str, ...
                'Value', value, 'FontSize', MainGuiStyle.FONT, ...
                'ForegroundColor', C.ink, 'BackgroundColor', C.panel, varargin{:});
        end

        function h = popup(parent, items, varargin)
            C = MainGuiStyle.col();
            h = uicontrol('Parent', parent, 'Style', 'popupmenu', 'String', items, ...
                'FontSize', MainGuiStyle.FONT, 'ForegroundColor', C.ink, ...
                'BackgroundColor', C.field, varargin{:});
        end

        function h = slider(parent, varargin)
            h = uicontrol('Parent', parent, 'Style', 'slider', ...
                'BackgroundColor', MainGuiStyle.col().panel, varargin{:});
        end

        %% ---- buttons
        function h = button(parent, str, kind, varargin)
            % kind: 'normal' (default), 'primary' (blue, the section's main action),
            % 'stop' (brick red), 'quiet' (small, for secondary actions).
            if nargin < 3 || isempty(kind); kind = 'normal'; end
            h = uicontrol('Parent', parent, 'Style', 'pushbutton', 'String', str, varargin{:});
            MainGuiStyle.paintButton(h, kind);
        end

        function h = toggle(parent, str, varargin)
            h = uicontrol('Parent', parent, 'Style', 'togglebutton', 'String', str, ...
                'FontSize', MainGuiStyle.FONT, varargin{:});
        end

        function paintButton(h, kind)
            C = MainGuiStyle.col();
            switch kind
                case 'primary'
                    set(h, 'BackgroundColor', C.accent, 'ForegroundColor', [1 1 1], ...
                        'FontWeight', 'bold', 'FontSize', MainGuiStyle.FONT + 2);
                case 'stop'
                    set(h, 'BackgroundColor', C.stop, 'ForegroundColor', [1 1 1], ...
                        'FontWeight', 'bold', 'FontSize', MainGuiStyle.FONT + 2);
                case 'pending'
                    set(h, 'BackgroundColor', C.accent, 'ForegroundColor', [1 1 1], ...
                        'FontWeight', 'bold', 'FontSize', MainGuiStyle.FONT);
                case 'quiet'
                    set(h, 'BackgroundColor', [0.94 0.94 0.94], 'ForegroundColor', C.ink2, ...
                        'FontWeight', 'normal', 'FontSize', MainGuiStyle.FONT - 1);
                otherwise
                    set(h, 'BackgroundColor', [0.94 0.94 0.94], 'ForegroundColor', C.ink, ...
                        'FontWeight', 'normal', 'FontSize', MainGuiStyle.FONT);
            end
        end

        %% ---- states
        function setFieldEnabled(h, tf)
            % Grey out an edit field that does not apply right now (it keeps its value).
            C = MainGuiStyle.col();
            for k = 1:numel(h)
                if tf
                    set(h(k), 'Enable', 'on', 'BackgroundColor', C.field, 'ForegroundColor', C.ink);
                else
                    set(h(k), 'Enable', 'off', 'BackgroundColor', C.off, 'ForegroundColor', C.ink3);
                end
            end
        end

        function setPending(h, tf)
            % A target typed but not yet applied: pale amber, bold. Back to a plain field
            % when cleared.
            C = MainGuiStyle.col();
            for k = 1:numel(h)
                if tf
                    set(h(k), 'BackgroundColor', C.pending, 'ForegroundColor', C.pendingInk, ...
                        'FontWeight', 'bold');
                else
                    set(h(k), 'BackgroundColor', C.field, 'ForegroundColor', C.ink, ...
                        'FontWeight', 'normal');
                end
            end
        end

        %% ---- numbers
        function s = num(x, digits)
            % A coordinate or setting for a field: up to `digits` decimals, no trailing
            % zeros, no exponent for ordinary numbers.
            if nargin < 2; digits = 3; end
            if isempty(x) || ~isfinite(x)
                s = '';
                return
            end
            s = sprintf(sprintf('%%.%df', digits), x);
            if contains(s, '.')
                s = regexprep(s, '0+$', '');
                s = regexprep(s, '\.$', '');
            end
            if strcmp(s, '-0'); s = '0'; end
        end

        function [s, scale, unit] = signal(v, baseUnit, digits)
            % A detector value with a readable unit: 0.0523 V -> '52.3 mV'.
            % baseUnit: 'V' (photodiode) or 'kcps' (counter; never rescaled).
            if nargin < 3; digits = 4; end
            scale = 1; unit = baseUnit;
            if strcmp(baseUnit, 'V') && isfinite(v)
                [scale, unit] = MainGuiStyle.voltScale(abs(v));
            end
            if ~isfinite(v)
                s = ['—  ' unit];
                return
            end
            s = [sprintf('%.*g', digits, v * scale) '  ' unit];
        end

        function [scale, unit] = voltScale(maxAbs)
            % Pick V / mV / µV so that the numbers sit between 1 and 1000.
            if isempty(maxAbs) || maxAbs >= 1 || maxAbs == 0 || ~isfinite(maxAbs)
                scale = 1; unit = 'V';
            elseif maxAbs >= 1e-3
                scale = 1e3; unit = 'mV';
            else
                scale = 1e6; unit = [char(181) 'V'];
            end
        end

        function s = duration(seconds)
            % '42 s', '3:05 min', '1:12 h'
            if ~isfinite(seconds) || seconds < 0
                s = '—';
            elseif seconds < 60
                s = sprintf('%.0f s', seconds);
            elseif seconds < 3600
                s = sprintf('%d:%02d min', floor(seconds/60), floor(mod(seconds, 60)));
            else
                s = sprintf('%d:%02d h', floor(seconds/3600), floor(mod(seconds, 3600)/60));
            end
        end

        function u = micron()
            u = [char(181) 'm'];
        end
    end
end
