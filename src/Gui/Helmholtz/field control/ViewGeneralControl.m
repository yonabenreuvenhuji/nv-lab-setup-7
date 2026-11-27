classdef ViewGeneralControl < GuiComponent & EventListener
    %VIEWGENERALCONTROL The "Supplies" section of the Helmholtz window:
    %
    %   Supplies
    %   Control [ ON ][ AUTO ][ OFF ]    ( OUTPUT OFF )    [ ] Flip direction
    %
    %   ON    the supplies' output is on
    %   AUTO  on only while an experiment runs
    %   OFF   off
    % The pill shows whether the output is on right now. Flip direction is greyed out when
    % the setup has no flip switch (flip_direction.available in setupInfo.json).

    properties (Constant)
        CONTROL_STATE = {'ON', 'AUTO', 'OFF'};
        TIPS = {'Supplies on', 'On only while an experiment runs', 'Supplies off'};
    end

    properties
        btnState        % 1x3 toggle buttons, in the order of CONTROL_STATE
        dspOutput       % output indicator
        cbxFlipDir      % checkbox
        helm
    end

    methods
        function obj = ViewGeneralControl(parent, controller)
            expNames = Experiment.getExperimentNames();
            namesToListenTo = [expNames(:)', {Helmholtz.NAME}];

            obj@GuiComponent(parent, controller);
            obj@EventListener(namesToListenTo);
            obj.helm = getObjByName('Helmholtz');
            S = MainGuiStyle;

            [panel, inner] = MainGuiStyle.section(parent.component, 'Supplies');
            obj.component = panel;
            row = MainGuiStyle.hbox(inner, 6);
            MainGuiStyle.label(row, 'Control', 'left');
            obj.btnState = gobjects(1, 3);
            for k = 1:3
                obj.btnState(k) = MainGuiStyle.toggle(row, obj.CONTROL_STATE{k}, ...
                    'TooltipString', obj.TIPS{k}, 'FontWeight', 'bold', ...
                    'Callback', @(h,e) obj.setControlState(k));
            end
            MainGuiStyle.empty(row);
            obj.dspOutput = uicontrol('Parent', row, 'Style', 'edit', 'Enable', 'inactive', 'String', 'OUTPUT OFF', ...
                'FontSize', S.FONT + 1, 'FontWeight', 'bold', 'ForegroundColor', [1 1 1], ...
                'HorizontalAlignment', 'center', 'TooltipString', 'Whether the supplies'' output is on now');
            MainGuiStyle.empty(row);
            obj.cbxFlipDir = MainGuiStyle.checkbox(row, 'Flip direction', false, ...
                'Callback', @(h,e) obj.cbxFlipDirCallback());
            if obj.helm.canFlip
                obj.cbxFlipDir.TooltipString = 'Reverse the field with the flip switch';
            else
                set(obj.cbxFlipDir, 'Enable', 'off', 'TooltipString', ...
                    'This setup has no flip switch (flip_direction.available is false in setupInfo.json)');
            end
            row.Widths = [50 50 58 50 -1 110 -1 110];
            inner.Heights = S.ROW_H;

            obj.height = MainGuiStyle.sectionHeight(S.ROW_H);
            obj.width = 520;
            obj.refresh;
        end

        function refresh(obj)
            C = MainGuiStyle.col();
            k = find(strcmp(obj.helm.helmControl, obj.helm.CONTROL_STATE), 1);
            for j = 1:3
                if j == k
                    set(obj.btnState(j), 'Value', 1, 'BackgroundColor', C.accent, 'ForegroundColor', [1 1 1]);
                else
                    set(obj.btnState(j), 'Value', 0, 'BackgroundColor', [0.94 0.94 0.94], 'ForegroundColor', C.ink);
                end
            end
            if obj.helm.helmOn
                set(obj.dspOutput, 'String', 'OUTPUT ON', 'BackgroundColor', C.ok);
            else
                set(obj.dspOutput, 'String', 'OUTPUT OFF', 'BackgroundColor', C.stop);
            end
            obj.cbxFlipDir.Value = obj.helm.flip;
        end

        function setControlState(obj, k)
            obj.helm.control(obj.helm.CONTROL_STATE{k});
            obj.refresh();     % (also when the state did not change: puts the toggles back)
        end

        function cbxFlipDirCallback(obj)
            obj.helm.flipDirection(obj.cbxFlipDir.Value);
        end
    end

    %% overridden from EventListener
    methods
        function onEvent(obj, event)
            if event.isError ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_OUTPUT_STATE_CHANGED) ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_FLIP_STATE_CHANGED) ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_CONTROL_STATE_CHANGED)
                obj.refresh();
            end
        end
    end
end
