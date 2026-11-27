classdef ViewStagePanelGeneral < GuiComponent
    %VIEWSTAGEPANELGENERAL The "Stage" section: reconnect the controller, and the closed-loop
    %switch (for stages that have both loop modes).

    properties
        btnResetStage   % button
        cbxClosedLoop   % checkbox

        stageName       % string
    end

    methods
        function obj = ViewStagePanelGeneral(parent, controller, stage)
            obj@GuiComponent(parent, controller);
            obj.stageName = stage.name;

            [panel, inner] = MainGuiStyle.section(parent.component, 'Stage');
            obj.component = panel;
            row = MainGuiStyle.hbox(inner, 8);

            hasClosedLoop = stage.hasClosedLoop;
            hasOpenLoop = stage.hasOpenLoop;
            obj.btnResetStage = MainGuiStyle.button(row, 'Reconnect', 'normal', ...
                'Callback', @obj.btnResetStageCallback, ...
                'TooltipString', 'Reconnect to the stage controller (what "Reset Stage" did)');
            obj.cbxClosedLoop = MainGuiStyle.checkbox(row, 'Closed loop', hasClosedLoop, ...
                'Enable', BooleanHelper.boolToOnOff(hasClosedLoop && hasOpenLoop), ...
                'Callback', @obj.cbxClosedLoopCallback);
            if hasClosedLoop && hasOpenLoop
                obj.cbxClosedLoop.TooltipString = 'Closed loop: position feedback. Open loop: the stage can only be moved by the joystick/scan';
            else
                obj.cbxClosedLoop.TooltipString = 'This stage has only one loop mode';
            end
            row.Widths = [90 -1];
            inner.Heights = MainGuiStyle.ROW_H;

            obj.height = MainGuiStyle.sectionHeight(MainGuiStyle.ROW_H);
            obj.width = 220;
        end

        %%%% Callbcaks %%%%
        function btnResetStageCallback(obj,~,~)
            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
            stage.Reconnect;
            stage.sendEventPositionChanged;     % It might have moved
        end

        function cbxClosedLoopCallback(obj,~,~)
            stage = getObjByName(obj.stageName);
            if isempty(stage); throwBaseObjException(obj.stageName); end
            mode = BooleanHelper.ifTrueElse(obj.cbxClosedLoop.Value, 'Closed', 'Open');
            stage.ChangeLoopMode(mode);
        end
    end
end
