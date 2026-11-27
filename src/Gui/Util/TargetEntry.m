classdef TargetEntry
    %TARGETENTRY The "now | target" pattern for setting several values at once.
    %
    % Each value has a read-only field showing what it is NOW, and next to it an always
    % visible TARGET field. Typing a target (and Enter) makes it pending: the field turns
    % amber and the apply button lights up. Apply sets every pending value at once and clears
    % the targets; Clear drops them. An empty target means "leave this one as it is".
    %
    % It replaces "click the value, type into the blue box that appears, press Change To
    % Blue": the box only appeared after a click, and its pale blue was easy to miss.
    % Used by the Helmholtz current / field panels; the main window's stage Position section
    % works the same way.

    methods (Static)
        function styleNow(h)
            % The present value: grey, not editable
            C = MainGuiStyle.col();
            set(h, 'Enable', 'inactive', 'BackgroundColor', C.readonly, 'ForegroundColor', C.ink, ...
                'TooltipString', 'The present value');
        end

        function styleTarget(h, tip)
            if nargin < 2
                tip = 'Type a target and press Enter (it turns amber); then press the Set button. Empty = no change.';
            end
            C = MainGuiStyle.col();
            set(h, 'Enable', 'on', 'BackgroundColor', C.field, 'ForegroundColor', C.ink, ...
                'FontWeight', 'normal', 'TooltipString', tip);
        end

        function showPending(h, tf)
            MainGuiStyle.setPending(h, tf);
        end

        function paintApply(btn, active)
            % The apply button: blue and bold when there is something to apply
            if active
                MainGuiStyle.paintButton(btn, 'pending');
                btn.Enable = 'on';
            else
                MainGuiStyle.paintButton(btn, 'normal');
                btn.Enable = 'off';
            end
        end

        function box = header(parent, widths, nowColumn, spacing)
            % A header row over a grid: 'now' and 'target' over the column that holds them
            if nargin < 4; spacing = 5; end
            box = uix.HBox('Parent', parent, 'Spacing', spacing);
            for k = 1:numel(widths)
                if k == nowColumn
                    pair = uix.HBox('Parent', box, 'Spacing', 0);
                    uicontrol('Parent', pair, 'Style', 'text', 'String', 'now', 'FontSize', 8, ...
                        'FontWeight', 'bold', 'ForegroundColor', MainGuiStyle.col().ink3);
                    uicontrol('Parent', pair, 'Style', 'text', 'String', 'target', 'FontSize', 8, ...
                        'FontWeight', 'bold', 'ForegroundColor', MainGuiStyle.col().pendingInk);
                else
                    uix.Empty('Parent', box);
                end
            end
            box.Widths = widths;
        end
    end
end
