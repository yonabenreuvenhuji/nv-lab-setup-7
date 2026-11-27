classdef ViewStagePanelTrack < GuiComponent & EventListener
    %VIEWSTAGEPANELTRACK The "Tracking" section: find the brightest point near the stage
    %position (opens the tracker window), once or continuously.

    properties
        btnTrack        % button. Starts tracker GUI
        cbxContinuous   % checkbox

        mStageName      % string. Name of stage being tracked
    end

    methods
        function obj = ViewStagePanelTrack(parent, controller, stage, laser) %#ok<INUSD>
            obj@GuiComponent(parent, controller);
            obj@EventListener(TrackablePosition.NAME)
            obj.mStageName = stage.name;

            [panel, inner] = MainGuiStyle.section(parent.component, 'Tracking');
            obj.component = panel;
            row = MainGuiStyle.hbox(inner, 8);
            obj.btnTrack = MainGuiStyle.button(row, 'Track', 'normal', ...
                'FontWeight', 'bold', 'Callback', @obj.btnTrackCallback, ...
                'TooltipString', 'Move to the brightest point around the current position');
            obj.cbxContinuous = MainGuiStyle.checkbox(row, 'Keep tracking', false, ...
                'Callback', @obj.cbxContinuousCallback, ...
                'TooltipString', 'Track again and again (continuous tracking)');
            row.Widths = [80 -1];
            inner.Heights = MainGuiStyle.ROW_H;

            obj.height = MainGuiStyle.sectionHeight(MainGuiStyle.ROW_H);
            obj.width = 220;
        end

        function refresh(obj)
            trackablePos = getObjByName(TrackablePosition.NAME);
            if isempty(trackablePos); throwBaseObjException(TrackablePosition.NAME); end
            obj.cbxContinuous.Value = trackablePos.isRunningContinuously;
        end

        %%%% Callbcaks %%%%
        function btnTrackCallback(obj,~,~)
            % We need to first make sure we have a working experiment, and
            % that it is running the right stage and laser
            trackablePos = getObjByName(TrackablePosition.NAME);
            if isempty(trackablePos)
                % We create the tracker anew, with this stage
                trackablePos = TrackablePosition.getInstance(obj.mStageName);
            else
                % Just in case the stage is not the right one,
                trackablePos.mStageName = obj.mStageName;
            end
            trackablePos.resetTrack;
            trackablePos.startTrack;
        end

        function cbxContinuousCallback(obj,~,~)
            trackablePos = getObjByName(TrackablePosition.NAME);
            if isempty(trackablePos); throwBaseObjException(TrackablePosition.NAME); end
            trackablePos.isRunningContinuously = obj.cbxContinuous.Value;
        end
    end

    %% overridden from EventListener
    methods
        function onEvent(obj, event)
            if (isfield(event.extraInfo, Trackable.EVENT_CONTINUOUS_TRACKING_CHANGED) ...
                    && strcmp(event.creator.NAME, TrackablePosition.NAME)) ...
                    || event.isError
                obj.refresh();
            end
        end
    end
end
