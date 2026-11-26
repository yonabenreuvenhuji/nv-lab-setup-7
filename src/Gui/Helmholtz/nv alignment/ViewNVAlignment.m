classdef ViewNVAlignment < GuiComponent & EventListener
    %VIEWNVALIGNMENT Put B along an NV axis: one row per signed NV direction.
    %
    % With a coil calibration loaded, B is in the NV crystal frame and every NV axis is at
    % theta 54.7356 (or 125.2644), phi +-45 / +-135 -- so aligning is a setting, not a scan.
    % Each row says how far along that direction the supplies reach (they are unipolar, so
    % usually only one or two directions are reachable at all), and "Go" sets |B| along it.

    properties
        helm
        edtRho          % |B| to align at, G
        lblDir          % 1x8
        lblAngles       % 1x8
        lblReach        % 1x8
        btnGo           % 1x8
        btnNearest
        lblStatus
    end

    methods
        function obj = ViewNVAlignment(parent, controller)
            expNames = Experiment.getExperimentNames();
            namesToListenTo = {expNames{:}, Helmholtz.NAME}; %#ok<CCAT>
            obj@GuiComponent(parent, controller);
            obj@EventListener(namesToListenTo);
            obj.helm = getObjByName('Helmholtz');

            panelMain = uix.Panel('Parent', parent.component, 'Title', ...
                'Align B with an NV axis', 'Padding', 5);
            vbox = uix.VBox('Parent', panelMain, 'Spacing', 4);

            hb = uix.HBox('Parent', vbox, 'Spacing', 5);
            uicontrol(obj.PROP_LABEL{:}, 'Parent', hb, 'String', '|B| (G)');
            obj.edtRho = uicontrol(obj.PROP_EDIT{:}, 'Parent', hb, ...
                'String', sprintf('%.2f', obj.helm.Brho), ...
                'TooltipString', 'field magnitude to use for Go');
            obj.btnNearest = uicontrol(obj.PROP_BUTTON{:}, 'Parent', hb, ...
                'String', 'B to nearest NV', ...
                'TooltipString', 'turn the present field onto the closest NV axis, keeping |B|');
            hb.Widths = [70 70 -1];

            grid = uix.Grid('Parent', vbox, 'Spacing', 3);
            dirs = Helmholtz.nvDirections();
            obj.lblDir = gobjects(1, 8); obj.lblAngles = gobjects(1, 8);
            obj.lblReach = gobjects(1, 8); obj.btnGo = gobjects(1, 8);
            for k = 1:8
                obj.lblDir(k) = uicontrol(obj.PROP_TEXT_NORMAL{:}, 'Parent', grid, ...
                    'String', sprintf('[%+d %+d %+d]', round(dirs(k, :)*sqrt(3))), ...
                    'FontName', 'FixedWidth');
            end
            for k = 1:8
                [~, t, p] = Helmholtz.cart2sphPhys(dirs(k, :));
                obj.lblAngles(k) = uicontrol(obj.PROP_TEXT_NORMAL{:}, 'Parent', grid, ...
                    'String', sprintf('theta %.2f  phi %.0f', t, p), 'FontName', 'FixedWidth');
            end
            for k = 1:8
                obj.lblReach(k) = uicontrol(obj.PROP_TEXT_NORMAL{:}, 'Parent', grid, ...
                    'String', '', 'FontName', 'FixedWidth', ...
                    'TooltipString', '|B| range the supplies can make along this direction');
            end
            for k = 1:8
                obj.btnGo(k) = uicontrol(obj.PROP_BUTTON_SMALL{:}, 'Parent', grid, ...
                    'String', 'Go', 'Callback', @(h,e) obj.goCallback(k));
            end
            set(grid, 'Widths', [80 135 115 40], 'Heights', 20 * ones(1, 8));

            obj.lblStatus = uicontrol(obj.PROP_TEXT_NORMAL{:}, 'Parent', vbox, ...
                'String', '', 'HorizontalAlignment', 'left');
            vbox.Heights = [26, 8*20 + 7*3, 32];

            obj.btnNearest.Callback = @(h,e) obj.helm.BtoNV();
            obj.edtRho.Callback = @(h,e) obj.refresh();

            obj.height = 26 + 8*23 + 32 + 4*4 + 30;
            obj.width = 400;
            obj.refresh();
        end

        function refresh(obj)
            rho = str2double(obj.edtRho.String);
            if ~isfinite(rho) || rho < 0
                rho = obj.helm.Brho;
                obj.edtRho.String = sprintf('%.2f', rho);
            end
            T = obj.helm.nvTable();
            for k = 1:8
                if isnan(T(k).rhoMin)
                    obj.lblReach(k).String = 'not reachable';
                    ok = false;
                else
                    obj.lblReach(k).String = sprintf('%.1f - %.1f G', T(k).rhoMin, T(k).rhoMax);
                    ok = rho >= T(k).rhoMin - 1e-9 && rho <= T(k).rhoMax + 1e-9;
                end
                obj.btnGo(k).Enable = BooleanHelper.boolToOnOff(ok);
                if T(k).angleToB < 0.05 && obj.helm.Brho > 0
                    obj.lblDir(k).BackgroundColor = [0.75 1 0.75];
                else
                    obj.lblDir(k).BackgroundColor = [1 1 1];
                end
            end
            [kN, c] = obj.helm.nearestNV();
            status = sprintf('|B| = %.2f G, %.2f deg from [%+d %+d %+d]', obj.helm.Brho, ...
                acosd(min(c, 1)), round(T(kN).dir));
            if obj.helm.isCalibrated
                obj.lblStatus.String = {status, 'angles in the NV crystal frame'};
                obj.lblStatus.ForegroundColor = [0 0 0];
            else
                obj.lblStatus.String = {status, ...
                    'UNCALIBRATED: nominal coil frame, alignment only approximate'};
                obj.lblStatus.ForegroundColor = [0.8 0 0];
            end
        end

        function goCallback(obj, k)
            rho = str2double(obj.edtRho.String);
            obj.helm.alignToNV(k, rho);
        end
    end

    methods
        function onEvent(obj, event)
            if event.isError ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_MAGNETIC_FIELD_CHANGED) ...
                    || isfield(event.extraInfo, Helmholtz.EVENT_CALIBRATION_CHANGED)
                obj.refresh();
            end
        end
    end
end
