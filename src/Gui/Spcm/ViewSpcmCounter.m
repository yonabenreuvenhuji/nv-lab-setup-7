classdef ViewSpcmCounter < ViewVBox & EventListener
    %VIEWSPCMCOUNTER The live detector reading: big number(s), the trace, the controls.
    %
    %   SIGNAL                    REFERENCE                  (two only in "signal & reference")
    %   52.31 mV                  80.04 mV
    %   mean 52.18 · sd 0.21      mean 80.02 · sd 0.09       signal / reference = 0.6535
    %   +------------------------------------------------------------------------------+
    %   |  trace, in the fit figures' style: blue signal, amber reference               |
    %   +------------------------------------------------------------------------------+
    %   [ Start ]  [Reset]  [Save]  [Pop out]       Integration (ms) [100]
    %   Input [ Signal (input 1) ▾]                 [x] Show last [ 50 ] points
    %
    % The readout is the LAST reading, with the mean and standard deviation of the points
    % shown in the trace. With "signal & reference" both inputs are read at once, so there
    % are two readouts (in the two series colours) and their ratio. Volts are shown as V,
    % mV or µV, whichever reads best; the trace's axis uses the same unit.
    %
    % Lifetime and g2 views (time-tagger setups) keep their own plots.

    properties
        vAxes           % axes view, to use for the plotting

        btnStartStop
        btnReset
        btnSave
        btnPopout
        edtIntegrationTime

        wrap            % positive integer, how many records in plot
        cbxUsingWrap
        edtWrap

        readouts        % struct array: name, value, stats (uicontrols) per series
        tvRatio         % ratio line (signal & reference)
        vInput          % ViewPhotodiodeControl, if any

        leftLabel
        isStandalone = false;
        popoutFcn = []
    end

    properties (Access = private)
        hSeries = gobjects(0)   % errorbar per series in the trace
        lastLimits = []         % [xlim ylim] this view set last; different now = the user zoomed / panned
    end

    properties (Constant)
        BOTTOM_LABEL = 'time (s)'; % Text for horiz. axis
        DEFAULT_WRAP_VALUE = 50;    % Value of wrap set in initiation
        DEFAULT_USING_WRAP = true;  % boolean, does this window uses wrap
        READOUT_H = 62
    end

    methods
        function obj = ViewSpcmCounter(parent, controller, varargin)
            obj@ViewVBox(parent, controller, 0, MainGuiStyle.GAP);
            obj@EventListener(SpcmCounter.NAME);
            C = MainGuiStyle.col();
            S = MainGuiStyle;
            obj.wrap = obj.DEFAULT_WRAP_VALUE;
            if ~isempty(varargin) && strcmp(varargin{1}, 'isStandalone')
                obj.isStandalone = varargin{2};
            end
            if obj.isStandalone
                obj.component.Padding = 10;
                fig = ancestor(obj.component, 'figure');
                set(fig, 'Color', C.panel);
            end
            obj.component.BackgroundColor = C.panel;
            spcm = getObjByName(Spcm.NAME);
            obj.leftLabel = Spcm.signalUnits();     % 'kcps', or 'V' with the photodiode

            %%% Readouts %%%
            readoutRow = MainGuiStyle.hbox(obj.component, 10);
            obj.readouts = struct('box', {}, 'name', {}, 'value', {}, 'stats', {});
            for k = 1:3     % (3: the ratio, in "signal & reference")
                box = MainGuiStyle.vbox(readoutRow, 0);
                r.box = box;
                r.name = uicontrol('Parent', box, 'Style', 'text', 'String', '', ...
                    'FontSize', S.FONT - 1, 'FontWeight', 'bold', 'BackgroundColor', C.panel, ...
                    'HorizontalAlignment', 'left');
                r.value = uicontrol('Parent', box, 'Style', 'text', 'String', '—', ...
                    'FontSize', S.FONT_BIG + 4*obj.isStandalone, 'FontWeight', 'bold', ...
                    'ForegroundColor', C.ink, 'BackgroundColor', C.panel, 'HorizontalAlignment', 'left');
                r.stats = MainGuiStyle.hint(box, '');
                box.Heights = [14 -1 15];
                obj.readouts(k) = r;
            end
            obj.tvRatio = obj.readouts(3).value;
            set(obj.readouts(3).value, 'FontSize', S.FONT_BIG - 4 + 4*obj.isStandalone);
            readoutRow.Widths = [-1 0 0];

            %%% Trace %%%
            plotBox = uipanel('Parent', obj.component, 'BorderType', 'none', 'BackgroundColor', C.surface);
            obj.vAxes = axes('Parent', plotBox, 'ActivePositionProperty', 'outerposition');
            obj.styleAxes();
            MainGuiStyle.makeInteractive(obj.vAxes);
            axes('Visible', 'off', 'Position', [0 0 eps eps]);     % a spare current axes, so stray plot() calls don't land here

            %%% Controls %%%
            row1 = MainGuiStyle.hbox(obj.component, S.GAP);
            obj.btnStartStop = ButtonStartStop(row1);
            obj.btnStartStop.colors = {C.accent, C.stop};
            obj.btnStartStop.isRunning = false;
            set(obj.btnStartStop.component, 'FontSize', S.FONT + 2, 'TooltipString', 'Read the detector continuously');
            obj.btnStartStop.startCallback = @obj.btnStartCallback;
            obj.btnStartStop.stopCallback = @obj.btnStopCallback;
            obj.btnReset = MainGuiStyle.button(row1, 'Reset', 'normal', 'Callback', @obj.btnResetCallback, ...
                'TooltipString', 'Clear the trace');
            obj.btnSave = MainGuiStyle.button(row1, 'Save', 'normal', 'Callback', @obj.btnSaveCallback, ...
                'TooltipString', 'Save the trace');
            if obj.isStandalone
                MainGuiStyle.empty(row1);
            else
                obj.btnPopout = MainGuiStyle.button(row1, 'Pop out', 'quiet', 'Callback', @obj.btnPopoutCallback, ...
                    'TooltipString', 'Open the reading in its own, larger window');
            end
            MainGuiStyle.label(row1, 'Integration (ms)', 'right');
            obj.edtIntegrationTime = MainGuiStyle.edit(row1, num2str(SpcmCounter.INTEGRATION_TIME_DEFAULT_MILLISEC), ...
                'Callback', @obj.edtIntegrationTimeCallback, 'TooltipString', 'Time per reading, in ms');
            row1.Widths = [72 56 52 60 -1 54];

            row2 = MainGuiStyle.hbox(obj.component, S.GAP);
            hasInput = spcm.hasPhotodiode();
            if hasInput
                obj.vInput = ViewPhotodiodeControl(struct('component', row2), controller);
            else
                MainGuiStyle.empty(row2);
            end
            obj.cbxUsingWrap = MainGuiStyle.checkbox(row2, 'Show last', obj.DEFAULT_USING_WRAP, ...
                'Callback', @obj.cbxUsingWrapCallback, 'TooltipString', 'Show only the latest points (off: all of them)');
            obj.edtWrap = MainGuiStyle.edit(row2, num2str(obj.DEFAULT_WRAP_VALUE), 'Callback', @obj.edtWrapCallback);
            MainGuiStyle.label(row2, 'points', 'left');
            row2.Widths = [-1 72 44 40];

            % Time-tagger setups: counter selection, lifetime, g2
            shouldHaveControl = (spcm.hasAltCount() || spcm.hasG2() || spcm.hasLifetime()) && ~obj.isStandalone;
            heights = [obj.READOUT_H, -1, S.ROW_H, S.ROW_H];
            if shouldHaveControl
                spcmControlView = ViewSpcmControl(obj.component, controller);
                heights(end + 1) = spcmControlView.height;
            end
            obj.setHeights(heights);
            obj.height = 520;
            obj.width = 430;
            if obj.isStandalone
                obj.height = 560;
                obj.width = 900;
            end

            SpcmCounter.init();
            counter = getObjByName(SpcmCounter.NAME);
            counter.addGraphicAxes(obj.vAxes); % So that experiment could plot on it
            obj.update;     % There might already be records in the counter
        end

        function tf = isUsingWrap(obj)
            tf = obj.cbxUsingWrap.Value;
        end

        function refresh(obj)
            % Just uicontrols, not axes
            spcmCount = getObjByName(SpcmCounter.NAME);
            if isempty(spcmCount)
                obj.btnStartStop.isRunning = false;
            else
                obj.edtIntegrationTime.String = num2str(spcmCount.integrationTimeMillisec);
                obj.btnStartStop.isRunning = spcmCount.isRunning;
            end
            if ~isempty(obj.vInput)
                obj.vInput.update();
            end
        end

        function update(obj)
            % Axes AND uicontrols
            counter = getObjByName(SpcmCounter.NAME);
            if isempty(counter); return; end
            counter.addGraphicAxes(obj.vAxes); % So that experiment could plot on it

            switch counter.currentType
                case 1 % Counter / photodiode reading
                    if obj.isUsingWrap
                        [time, values, sterr, realTime] = counter.getRecords(obj.wrap);
                    else
                        [time, values, sterr, realTime] = counter.getRecords;
                    end
                    % getRecords pads a short history up to the window length, ending
                    % with a 0 at the last time (to fix the time axis). Those are not
                    % readings: not drawn, not in the readout.
                    padding = isnan(realTime);
                    values(:, padding) = NaN;
                    sterr(:, padding) = NaN;
                    obj.drawTrace(time, values, sterr);
                case 2 % Lifetime
                    obj.hSeries = gobjects(0);
                    spcm = getObjByName(Spcm.NAME);
                    binWidth = spcm.binWidth/1000; % in ns
                    bins = binWidth:binWidth:spcm.nBins*binWidth;
                    AxesHelper.fill(obj.vAxes, spcm.lastTimeHist, 1, bins, [], 'time (ns)', 'counts (kcps)');
                case 3 % G2
                    obj.hSeries = gobjects(0);
                    spcm = getObjByName(Spcm.NAME);
                    binWidth = spcm.binWidth/1000; % in ns
                    bins = -spcm.nBins*binWidth:binWidth:spcm.nBins*binWidth;
                    AxesHelper.fill(obj.vAxes, spcm.lastTimeG2, 1, bins, [], 'delay (ns)', 'g2');
            end
            obj.refresh;
        end
    end

    %% Drawing
    methods (Access = private)
        function styleAxes(obj)
            C = ExpPlotStyle.col();
            ax = obj.vAxes;
            set(ax, 'FontSize', 9 + 4*obj.isStandalone, 'XColor', C.ink2, 'YColor', C.ink2, 'TickDir', 'out', ...
                'Color', C.surface, 'Box', 'off', 'GridLineStyle', '-', 'GridAlpha', 0.12, 'Layer', 'bottom');
            grid(ax, 'on');
            xlabel(ax, obj.BOTTOM_LABEL, 'Color', C.ink2);
        end

        function [names, unitBase] = seriesNames(obj, nSeries)
            % What each row of the counter's records is
            spcm = getObjByName(Spcm.NAME);
            unitBase = obj.leftLabel;
            if spcm.hasPhotodiode() && isprop(spcm, 'currentChannelGUI')
                switch spcm.currentChannelGUI
                    case 1; names = {'Signal'};
                    case 2; names = {'Reference'};
                    case 3; names = {'Signal', 'Reference'};
                    case 4; names = {['Signal ' char(8722) ' reference']};
                    case 5; names = {['Signal / ' char(10216) 'reference' char(10217)]};
                    case 6; names = {'Laser block'};
                    otherwise; names = {'Reading'};
                end
            elseif nSeries == 2
                names = {'SPCM 1', 'SPCM 2'};
            else
                names = {'Counts'};
            end
            names = names(1:min(numel(names), max(nSeries, 1)));
        end

        function drawTrace(obj, time, values, sterr)
            C = ExpPlotStyle.col();
            ax = obj.vAxes;
            % Rows that hold data (the second only when two values are read)
            nSeries = 1;
            if size(values, 1) >= 2 && any(isfinite(values(2, :)) & values(2, :) ~= 0)
                nSeries = 2;
            end
            [names, unitBase] = obj.seriesNames(nSeries);
            nSeries = numel(names);

            % One unit for the axis and the readouts
            shown = values(1:nSeries, :);
            scale = 1; unit = unitBase;
            if strcmp(unitBase, 'V')
                [scale, unit] = MainGuiStyle.voltScale(max(abs(shown(isfinite(shown))), [], 'all'));
            end

            % (Re)create the series if needed
            if numel(obj.hSeries) ~= nSeries || ~all(isgraphics(obj.hSeries))
                cla(ax);
                set(ax, 'NextPlot', 'add');
                obj.hSeries = gobjects(1, nSeries);
                slots = {C.s1, C.s2};
                for k = 1:nSeries
                    h = errorbar(ax, NaN, NaN, NaN);
                    set(h, 'Color', slots{k}, 'LineWidth', 1.1, 'Marker', 'o', 'MarkerSize', 4 + 2*obj.isStandalone, ...
                        'MarkerFaceColor', slots{k}, 'MarkerEdgeColor', C.surface, ...
                        'DisplayName', names{k});
                    if isprop(h, 'CapSize'); h.CapSize = 0; end
                    try     % data tips: time, value, error
                        h.DataTipTemplate.DataTipRows(1).Label = 'time (s)';
                        h.DataTipTemplate.DataTipRows(2).Label = names{k};
                    catch
                    end
                    obj.hSeries(k) = h;
                end
                obj.styleAxes();
                MainGuiStyle.makeInteractive(ax);
                obj.lastLimits = [];
                if nSeries > 1
                    lg = legend(ax, obj.hSeries, names, 'Location', 'northwest', 'Orientation', 'horizontal');
                    set(lg, 'Box', 'off', 'TextColor', C.ink2, 'FontSize', 8 + 3*obj.isStandalone);
                else
                    legend(ax, 'off');
                end
            end

            % Zoomed or panned by the user (limits are not the ones set here last time):
            % leave the view alone until it is restored (double-click / toolbar) or Reset
            userView = ~isempty(obj.lastLimits) && ...
                (any(abs(ax.XLim - obj.lastLimits(1:2)) > 1e-12) || any(abs(ax.YLim - obj.lastLimits(3:4)) > 1e-12));
            for k = 1:nSeries
                e = sterr(k, :) * scale;
                e(~isfinite(e)) = 0;
                set(obj.hSeries(k), 'XData', time, 'YData', values(k, :) * scale, ...
                    'YNegativeDelta', e, 'YPositiveDelta', e, 'DisplayName', names{k});
            end
            if nSeries == 1
                ylabel(ax, sprintf('%s (%s)', lower(names{1}), unit), 'Color', C.ink2);
            else
                ylabel(ax, sprintf('photodiode (%s)', unit), 'Color', C.ink2);
            end
            if ~userView
                t = time(isfinite(time));
                if numel(t) >= 2 && max(t) > min(t)
                    xlim(ax, [min(t), max(t)]);
                end
                ylim(ax, 'auto');
                if nSeries > 1
                    yl = ylim(ax);
                    ylim(ax, [yl(1), yl(2) + 0.2*(yl(2) - yl(1))]);     % room for the legend
                end
                obj.lastLimits = [ax.XLim, ax.YLim];
            end

            obj.showReadouts(names, values(1:nSeries, :), scale, unit);
            drawnow;    % (a full drawnow: it is also where the reading loop sees Stop)
        end

        function showReadouts(obj, names, values, scale, unit)
            C = ExpPlotStyle.col();
            slots = {C.s1, C.s2};
            n = numel(names);
            for k = 1:2
                r = obj.readouts(k);
                if k > n
                    set([r.name r.value r.stats], 'String', '');
                    continue
                end
                v = values(k, :);
                v = v(isfinite(v));
                set(r.name, 'String', upper(names{k}), 'ForegroundColor', slots{k});
                if isempty(v)
                    set(r.value, 'String', ['—  ' unit]);
                    set(r.stats, 'String', 'no reading yet');
                else
                    set(r.value, 'String', sprintf('%.4g  %s', v(end) * scale, unit));
                    set(r.stats, 'String', sprintf('mean %.4g · sd %.2g', mean(v) * scale, std(v) * scale), ...
                        'TooltipString', sprintf('Mean and standard deviation of the %d readings in the trace', numel(v)));
                end
            end
            row = obj.readouts(1).box.Parent;
            r = obj.readouts(3);
            if n == 2
                sig = values(1, :); ref = values(2, :);
                ok = isfinite(sig) & isfinite(ref) & ref ~= 0;
                set(r.name, 'String', 'SIGNAL / REF', 'ForegroundColor', C.ink2);
                if any(ok)
                    last = find(ok, 1, 'last');
                    set(r.value, 'String', sprintf('%.4f', sig(last)/ref(last)));
                    set(r.stats, 'String', sprintf('mean %.4f', mean(sig(ok)./ref(ok))));
                else
                    set([r.value r.stats], 'String', '');
                end
                row.Widths = [-1 -1 -0.8];
            else
                set([r.name r.value r.stats], 'String', '');
                row.Widths = [-1 0 0];
            end
        end
    end

    %% Callbacks
    methods
        function cbxUsingWrapCallback(obj, ~, ~)
            MainGuiStyle.setFieldEnabled(obj.edtWrap, obj.isUsingWrap);
            obj.update;
        end
        function edtWrapCallback(obj, ~, ~)
            if ~ValidationHelper.isValuePositiveInteger(obj.edtWrap.String)
                EventStation.anonymousWarning('The number of points needs to be a positive integer! Reverting.')
                obj.edtWrap.String = num2str(obj.wrap);
            end
            obj.wrap = str2double(obj.edtWrap.String);
            obj.update;
        end
        function btnStartCallback(obj, ~, ~)
            spcmCount = obj.getCounter;
            spcmCount.run;
        end
        function btnStopCallback(obj, ~, ~)
            spcmCount = obj.getCounter;
            spcmCount.pause;
        end
        function btnResetCallback(obj, ~ ,~)
            spcmCount = obj.getCounter;
            spcmCount.resetHistory;
        end
        function btnSaveCallback(obj, ~ ,~)
            spcmCount = obj.getCounter;
            spcmCount.save;
        end
        function btnPopoutCallback(obj, ~, ~)
            if ~isempty(obj.popoutFcn)
                obj.popoutFcn();
            end
        end
        function edtIntegrationTimeCallback(obj, ~, ~)
            spcmCount = obj.getCounter;
            spcmCount.integrationTimeMillisec = str2double(obj.edtIntegrationTime.String);
            % The counter will take care of the rest
        end
    end

    methods (Static)
        function spcmCounter = getCounter
            spcmCounter = getObjByName(SpcmCounter.NAME);
            if isempty(spcmCounter)
                % If we don't already have one, we want to create it
                spcmCounter = SpcmCounter;
            end
        end
    end

    %% overridden from EventListener
    methods
        function onEvent(obj, event)
            if isfield(event.extraInfo, SpcmCounter.EVENT_DATA_UPDATED)   % event = update
                obj.update;
            else
                obj.refresh;
            end
            if isfield(event.extraInfo, SpcmCounter.EVENT_SPCM_COUNTER_RESET)    % event = reset
                for h = obj.hSeries
                    if isgraphics(h)
                        set(h, 'XData', NaN, 'YData', NaN, 'YNegativeDelta', NaN, 'YPositiveDelta', NaN);
                    end
                end
                obj.lastLimits = [];     % Reset also gives the automatic view back
                obj.update;
            end
        end
    end
end
