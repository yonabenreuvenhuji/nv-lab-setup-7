classdef HelmholtzAxesLayout
    %HELMHOLTZAXESLAYOUT The layout shared by the three value sections of the Helmholtz window
    %(spherical B, cartesian B, coil currents):
    %
    %   <title>
    %              Min     −    Now     +     Max    Target   |  Step (G)  [ 1 ]
    %   |B| (G)   2.55    [−]  45.208  [+]   87.87  [      ]  |  [ Set to target ]
    %   θ (deg)   0.27    [−]  33.705  [+]   87.31  [      ]  |  [ Clear targets ]
    %   φ (deg)   ...
    %
    % Min / Max: how far that value can go with the other two held (shown, not typed).
    % Now: the present value. Target: type one (Enter) -- it turns amber -- then "Set to
    % target" sets every pending value at once (TargetEntry).
    %
    % The title is a text row (MainGuiStyle.section), not a panel title: on R2025a+ the
    % GUI Layout Toolbox put a titled panel's first rows UNDER its title (the half-hidden
    % "Step:" rows).

    properties (Constant)
        WIDTHS = [64 46 28 74 28 46 74]     % label | min | − | now | + | max | target
    end

    methods (Static)
        function h = build(parent, title, labels, tips, stepUnit, stepTip)
            S = MainGuiStyle;
            C = MainGuiStyle.col();
            [h.panel, inner, h.title] = MainGuiStyle.section(parent, title);

            top = MainGuiStyle.hbox(inner, 14);
            grid = MainGuiStyle.grid(top, 4);
            % column 1: what each row is, with its unit
            MainGuiStyle.empty(grid);
            for i = 1:3
                MainGuiStyle.label(grid, labels{i}, 'left', 'FontWeight', 'bold', ...
                    'ForegroundColor', C.ink, 'TooltipString', tips{i});
            end
            % column 2: lowest reachable value
            MainGuiStyle.header(grid, 'Min', 'right');
            h.lowerLimLbl = gobjects(1, 3);
            for i = 1:3
                h.lowerLimLbl(i) = MainGuiStyle.hint(grid, '', 'right', ...
                    'TooltipString', 'Lowest value it can have with the other two held');
            end
            % column 3: one step down
            MainGuiStyle.empty(grid);
            h.btnMoveLeft = gobjects(1, 3);
            for i = 1:3
                h.btnMoveLeft(i) = MainGuiStyle.button(grid, char(8722), 'normal', 'FontSize', 12, ...
                    'TooltipString', 'One step down');
            end
            % column 4: now
            MainGuiStyle.header(grid, 'Now');
            h.tvCurB = gobjects(1, 3);
            for i = 1:3
                h.tvCurB(i) = MainGuiStyle.readout(grid, '');
            end
            % column 5: one step up
            MainGuiStyle.empty(grid);
            h.btnMoveRight = gobjects(1, 3);
            for i = 1:3
                h.btnMoveRight(i) = MainGuiStyle.button(grid, '+', 'normal', 'FontSize', 12, ...
                    'TooltipString', 'One step up');
            end
            % column 6: highest reachable value
            MainGuiStyle.header(grid, 'Max', 'left');
            h.upperLimLbl = gobjects(1, 3);
            for i = 1:3
                h.upperLimLbl(i) = MainGuiStyle.hint(grid, '', 'left', ...
                    'TooltipString', 'Highest value it can have with the other two held');
            end
            % column 7: target
            uicontrol('Parent', grid, 'Style', 'text', 'String', 'Target', ...
                'FontSize', S.FONT - 1, 'FontWeight', 'bold', 'ForegroundColor', C.pendingInk, ...
                'BackgroundColor', C.panel, 'HorizontalAlignment', 'center');
            h.edtCurB = gobjects(1, 3);
            for i = 1:3
                h.edtCurB(i) = MainGuiStyle.edit(grid, '');
            end
            set(grid, 'Widths', HelmholtzAxesLayout.WIDTHS, 'Heights', [S.HEAD_H, S.ROW_H * ones(1, 3)]);
            gridWidth = sum(HelmholtzAxesLayout.WIDTHS) + 4*numel(HelmholtzAxesLayout.WIDTHS);

            % right: step, apply, clear
            right = MainGuiStyle.vbox(top, S.GAP);
            MainGuiStyle.empty(right);      % (level with the header row)
            stepRow = MainGuiStyle.hbox(right, 4);
            MainGuiStyle.label(stepRow, sprintf('Step (%s)', stepUnit), 'right');
            h.edtSteps = MainGuiStyle.edit(stepRow, '', 'TooltipString', stepTip);
            stepRow.Widths = [-1 56];
            h.btnMoveToBlue = MainGuiStyle.button(right, 'Set to target', 'normal', 'Enable', 'off', ...
                'TooltipString', 'Set every value that has a target (amber) to it, all at once');
            h.btnClearTargets = MainGuiStyle.button(right, 'Clear targets', 'quiet', 'Enable', 'off', ...
                'TooltipString', 'Forget the targets');
            right.Heights = [S.HEAD_H - S.GAP + 4, S.ROW_H, S.ROW_H, S.ROW_H];
            top.Widths = [gridWidth, -1];

            inner.Heights = S.HEAD_H + 3*S.ROW_H + 3*4;
            h.height = MainGuiStyle.sectionHeight(inner.Heights);
            h.width = gridWidth + 14 + 150 + 2*S.PAD;
            h.hboxCurB = gobjects(1, 3);    % (now and target are separate columns now)
        end

        function s = fmtLim(v)
            if isnan(v); s = '--'; else; s = sprintf('%.2f', v); end
        end
    end
end
