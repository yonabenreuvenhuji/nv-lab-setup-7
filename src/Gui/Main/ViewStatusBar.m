classdef ViewStatusBar < GuiComponent & EventListener
    %VIEWSTATUSBAR One line at the bottom of the main window: the last error or warning, and
    %what sent it. Errors are red, warnings amber. The message stays (it used to vanish after
    %a few seconds, often before it was read) until the next one, or until "Clear"; the time
    %it arrived is shown, so an old message is not mistaken for a new one. The full text is
    %in the tooltip, and in the command window as always.

    properties (Access = protected)
        tvIcon
        tvMsg
        tvFrom
        btnClear
    end

    methods
        function obj = ViewStatusBar(parent, controller)
            obj@GuiComponent(parent, controller);
            obj@EventListener();
            obj.setListeningForAll(true);
            C = MainGuiStyle.col();

            obj.component = uix.HBox('Parent', parent.component, 'Spacing', 8, 'Padding', 4, ...
                'BackgroundColor', C.readonly);
            obj.tvIcon = uicontrol('Parent', obj.component, 'Style', 'text', 'String', '', ...
                'FontSize', MainGuiStyle.FONT + 2, 'FontWeight', 'bold', ...
                'BackgroundColor', C.readonly, 'HorizontalAlignment', 'center');
            obj.tvMsg = uicontrol('Parent', obj.component, 'Style', 'text', 'String', 'Ready', ...
                'FontSize', MainGuiStyle.FONT + 1, 'ForegroundColor', C.ink3, ...
                'BackgroundColor', C.readonly, 'HorizontalAlignment', 'left');
            obj.tvFrom = uicontrol('Parent', obj.component, 'Style', 'text', 'String', '', ...
                'FontSize', MainGuiStyle.FONT, 'ForegroundColor', C.ink3, ...
                'BackgroundColor', C.readonly, 'HorizontalAlignment', 'right');
            obj.btnClear = MainGuiStyle.button(obj.component, 'Clear', 'quiet', 'Visible', 'off', ...
                'Callback', @(h,e) obj.clear());
            obj.component.Widths = [20 -1 260 50];

            obj.width = 800;
            obj.height = 30;
        end

        function clear(obj)
            C = MainGuiStyle.col();
            set(obj.tvMsg, 'String', 'Ready', 'ForegroundColor', C.ink3, 'TooltipString', '');
            set(obj.tvIcon, 'String', '');
            set(obj.tvFrom, 'String', '');
            obj.btnClear.Visible = 'off';
        end
    end

    methods
        function out = onEvent(obj, event)
            out = false;
            if ~(event.isError || event.isWarning)
                return
            end
            C = MainGuiStyle.col();
            msg = event.extraInfo.(Event.ERROR_MSG);
            msg = strrep(msg, '\n', newline);    % (some senders pass an unformatted '\n')
            oneLine = regexprep(strtrim(msg), '\s*\n\s*', '  ');
            if event.isError
                color = C.err; icon = '!';
            else
                color = C.warn; icon = '!';
            end
            set(obj.tvIcon, 'String', icon, 'ForegroundColor', color);
            set(obj.tvMsg, 'String', oneLine, 'ForegroundColor', color, 'TooltipString', msg);
            set(obj.tvFrom, 'String', sprintf('%s  ·  %s', event.creator.name, char(datetime('now', 'Format', 'HH:mm:ss'))));
            obj.btnClear.Visible = 'on';
            out = true;
        end
    end
end
