classdef ViewLaser < GuiComponent
    %VIEWLASER The section for one laser:
    %   [x] Light on (fast AOM switch)    [x] Source on      (what the laser supports)
    %   Source  [x]  [ 50 % ]  [=====|=====]
    %   AOM     [x]  [ 80 % ]  [=======|===]

    methods
        function obj = ViewLaser(parent, controller, laserGate)
            obj@GuiComponent(parent, controller);
            sourceAvailBool = laserGate.isSourceAvail();
            aomAvailBool = laserGate.isAomAvail();
            S = MainGuiStyle;
            C = MainGuiStyle.col();

            [panel, inner] = MainGuiStyle.section(parent.component, sprintf('Laser %s %s', char(183), laserGate.name));
            obj.component = panel;
            innerView = struct('component', inner);

            switchRow = ViewBox(innerView, controller, @uix.HBox, 'Spacing', 10, 'BackgroundColor', C.panel);
            cbxFastOn = ViewBooleanSwitch(switchRow, controller, laserGate.aomSwitch, 'Light on (fast AOM switch)');
            set(cbxFastOn.cbxEnabled, 'BackgroundColor', C.panel, 'FontSize', S.FONT, 'ForegroundColor', C.ink, ...
                'TooltipString', 'The pulse generator''s gate on the AOM. Off: no light reaches the sample.');
            widths = -1;
            if sourceAvailBool && laserGate.source.canSetEnabled && ~laserGate.source.canSetValue
                % only an on/off switch for the source
                cbxSourceOn = ViewBooleanSwitch(switchRow, controller, laserGate.source, 'Source on');
                set(cbxSourceOn.cbxEnabled, 'BackgroundColor', C.panel, 'FontSize', S.FONT, 'ForegroundColor', C.ink);
                widths = [-1 -1];
            end
            switchRow.component.Widths = widths;
            heights = S.ROW_H;

            if sourceAvailBool
                sourceView = ViewLaserPart(innerView, controller, laserGate.source, 'Source');
                heights(end + 1) = sourceView.height;
            end
            if aomAvailBool
                aom = laserGate.aom;
                if isa(aom, 'AomDoubleNiDaqControlled')
                    aomView = ViewDoubleAom(innerView, controller, aom);
                else
                    aomView = ViewLaserPart(innerView, controller, aom, 'AOM power');
                end
                heights(end + 1) = aomView.height;
            end
            inner.Heights = heights;
            obj.height = MainGuiStyle.sectionHeight(MainGuiStyle.rowsHeight(heights));
            obj.width = 420;
        end
    end
end
