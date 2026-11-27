classdef ViewStage < ViewVBox
    %VIEWSTAGE Everything about one stage, as a column of sections:
    %
    %   Position     where the stage is, moving it, saved positions
    %   Scan area    what to scan: per axis scan / from / to / points / centre
    %   Scan         the Scan button, progress, pixel time, line order and pattern
    %   MW contrast  (scannable stages)
    %   Tracking | Stage
    %   Soft limits
    %   Tilt         (stages that have it)
    %
    % Each section view sets its own height; this view stacks them.

    methods
        function obj = ViewStage(parent, controller, stage, camera)
            if ~exist('camera', 'var'); camera = []; end
            obj@ViewVBox(parent, controller, 0, MainGuiStyle.GAP);
            C = MainGuiStyle.col();
            obj.component.BackgroundColor = C.window;

            % Title row: the stage's name
            title = uicontrol('Parent', obj.component, 'Style', 'text', ...
                'String', stage.name, 'FontSize', MainGuiStyle.FONT + 2, 'FontWeight', 'bold', ...
                'ForegroundColor', C.ink, 'BackgroundColor', C.window, 'HorizontalAlignment', 'left');
            heights = 20;
            width = 0;

            sections = {};
            sections{end + 1} = ViewStagePanelMovementControl(obj, controller, stage);
            if stage.isScannable
                sections{end + 1} = ViewStagePanelScanParams(obj, controller, stage);
                sections{end + 1} = ViewStagePanelScan(obj, controller, stage.name);
            end

            % Tracking (or the focal distance of camera setups) next to the stage's own settings
            row = ViewHBox(obj, controller, 0, MainGuiStyle.GAP);
            row.component.BackgroundColor = C.window;
            if stage.isScannable
                if isempty(camera)
                    vTrack = ViewStagePanelTrack(row, controller, stage);
                else
                    vTrack = ViewFocalDistance(row, controller);
                end
                vGeneral = ViewStagePanelGeneral(row, controller, stage);
                row.setWidths([-1 -1]);
                row.height = max(vTrack.height, vGeneral.height);
            else
                vGeneral = ViewStagePanelGeneral(row, controller, stage);
                row.setWidths(-1);
                row.height = vGeneral.height;
            end
            sections{end + 1} = row;

            sections{end + 1} = ViewStagePanelLimits(obj, controller, stage);
            if stage.isScannable
                % Rarely used (scans with the MW on), so last
                sections{end + 1} = ViewStagePanelMWContrast(obj, controller, stage);
            end
            if stage.tiltAvailable
                sections{end + 1} = ViewStagePanelTiltCorrection(obj, controller, stage);
            end

            for k = 1:numel(sections)
                heights(end + 1) = sections{k}.height; %#ok<AGROW>
                width = max(width, sections{k}.width);
            end
            obj.setHeights(heights);
            obj.height = MainGuiStyle.rowsHeight(heights);
            obj.width = width;
            set(title, 'TooltipString', sprintf('%s (axes: %s)', stage.name, upper(stage.availableAxes)));
        end
    end
end
