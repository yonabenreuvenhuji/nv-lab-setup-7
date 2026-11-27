classdef ViewNVAlignment < GuiComponent & EventListener
    %VIEWNVALIGNMENT Choose a direction and where the dips should land; Go sets the field.
    %
    % DIRECTION  one of the 8 signed NV axes (one orientation aligned: four lines) or the 6
    %            cube axes (all four NV orientations equivalent: two lines, 4x deep). Each is
    %            listed with the |B| range the supplies reach exactly along it.
    % |B| AND THE DIPS ARE ONE INPUT. Type |B|, or where the outer pair should land (left or
    %            right), or the inner pair, and the rest follow -- through the exact
    %            Hamiltonian, so the off-axis lines' second-order shift is included. Along a
    %            cube axis there is only one pair, so the inner boxes are off.
    % PREVIEW    nothing moves until Go. Meanwhile the expected-ESR plot shows the target
    %            spectrum dashed over the present one, and the bar shows where the target |B|
    %            sits in the direction's reachable range.
    % OUT OF EXACT REACH  the supplies are unipolar, so near a cube axis they can need a few
    %            mA of negative current. Then Go sets the closest reachable field instead, if
    %            it is within 1 deg, and the status line says how far off it is.

    properties (Constant)
        MAX_DEG = 1;            % the furthest "closest reachable" field Go will accept
        DIPS = {'outerLo', 'outerHi', 'innerLo', 'innerHi'};
    end

    properties
        helm
        popDir          % direction menu
        edtRho          % |B|
        edtDip          % 1x4: outer left, outer right, inner left, inner right
        btnGo
        btnClear
        lblStatus
        axBar           % the reach bar
        targets         % Helmholtz.targetTable(), as last listed
        rho             % the target |B|, G
        userTarget = false  % true once a target was typed or picked; until then the
                            % panel simply follows the present field
    end

    methods
        function obj = ViewNVAlignment(parent, controller)
            expNames = Experiment.getExperimentNames();
            namesToListenTo = {expNames{:}, Helmholtz.NAME}; %#ok<CCAT>
            obj@GuiComponent(parent, controller);
            obj@EventListener(namesToListenTo);
            obj.helm = getObjByName('Helmholtz');
            C = Helmholtz.plotColours();

            panelMain = uix.Panel('Parent', parent.component, 'Title', ...
                'Align B: direction, |B| and where the dips land', 'Padding', 5);
            vbox = uix.VBox('Parent', panelMain, 'Spacing', 5);

            hb = uix.HBox('Parent', vbox, 'Spacing', 5);
            uicontrol(obj.PROP_LABEL{:}, 'Parent', hb, 'String', 'Direction');
            obj.popDir = uicontrol(obj.PROP_POPUP{:}, 'Parent', hb, 'String', {' '}, ...
                'FontName', 'FixedWidth', 'TooltipString', ...
                'NV axes: one orientation aligned.  +-x,y,z: all four orientations equal.');
            hb.Widths = [85 -1];

            barBox = uicontainer('Parent', vbox, 'BackgroundColor', C.surface);
            obj.axBar = axes('Parent', barBox, 'Units', 'normalized', ...
                'Position', [0.04 0.62 0.92 0.3]);
            try
                obj.axBar.Toolbar.Visible = 'off';
                disableDefaultInteractivity(obj.axBar);
            catch
            end

            grid = uix.Grid('Parent', vbox, 'Spacing', 4);
            uicontrol(obj.PROP_LABEL{:}, 'Parent', grid, 'String', '|B| (G)');
            uicontrol(obj.PROP_LABEL{:}, 'Parent', grid, 'String', 'outer dips (MHz)', ...
                'TooltipString', 'the widest pair: the aligned orientation along an NV axis');
            uicontrol(obj.PROP_LABEL{:}, 'Parent', grid, 'String', 'inner dips (MHz)', ...
                'TooltipString', 'the pair nearest D: the three other orientations');
            obj.edtRho = uicontrol(obj.PROP_EDIT{:}, 'Parent', grid);
            obj.edtDip = gobjects(1, 4);
            obj.edtDip(1) = uicontrol(obj.PROP_EDIT{:}, 'Parent', grid, 'TooltipString', 'left (low)');
            obj.edtDip(3) = uicontrol(obj.PROP_EDIT{:}, 'Parent', grid, 'TooltipString', 'left (low)');
            uix.Empty('Parent', grid);
            obj.edtDip(2) = uicontrol(obj.PROP_EDIT{:}, 'Parent', grid, 'TooltipString', 'right (high)');
            obj.edtDip(4) = uicontrol(obj.PROP_EDIT{:}, 'Parent', grid, 'TooltipString', 'right (high)');
            set(grid, 'Widths', [130 -1 -1], 'Heights', [26 26 26]);

            hbGo = uix.HBox('Parent', vbox, 'Spacing', 8);
            obj.btnGo = uicontrol(obj.PROP_BUTTON_BIG_GREEN{:}, 'Parent', hbGo, 'String', 'Go', ...
                'TooltipString', 'set the field to the target above');
            obj.btnClear = uicontrol(obj.PROP_BUTTON{:}, 'Parent', hbGo, 'String', 'Reset to present', ...
                'TooltipString', 'forget the target: back to the present field and the nearest direction');
            hbGo.Widths = [-1 150];

            obj.lblStatus = uicontrol(obj.PROP_TEXT_NORMAL{:}, 'Parent', vbox, 'String', '', ...
                'HorizontalAlignment', 'left');
            vbox.Heights = [28, 58, 86, 38, 46];

            obj.popDir.Callback = @(h,e) obj.directionPicked();
            obj.edtRho.Callback = @(h,e) obj.rhoTyped();
            for k = 1:4
                obj.edtDip(k).Callback = @(h,e) obj.dipTyped(k);
            end
            obj.btnGo.Callback = @(h,e) obj.goCallback();
            obj.btnClear.Callback = @(h,e) obj.resetToPresent();

            obj.height = 28 + 58 + 86 + 38 + 46 + 5*5 + 30;
            obj.width = 420;
            obj.resetToPresent();
        end

        function resetToPresent(obj)
            % Target = the present |B| along the nearest listed direction.
            obj.listTargets();
            [~, k] = min([obj.targets.angleToB]);
            obj.popDir.Value = k;
            obj.rho = obj.helm.Brho;
            obj.userTarget = false;
            obj.update();
        end

        function directionPicked(obj)
            obj.userTarget = true;
            obj.update();
        end

        function listTargets(obj)
            obj.targets = obj.helm.targetTable();
            labels = cell(1, numel(obj.targets));
            for k = 1:numel(obj.targets)
                t = obj.targets(k);
                if isnan(t.rhoMin)
                    r = 'closest only';
                else
                    r = sprintf('%.1f-%.1f G', t.rhoMin, t.rhoMax);
                end
                if strcmp(t.kind, 'NV')
                    name = strrep(t.label, 'NV ', '');
                else
                    name = [strtok(t.label) ' all equal'];
                end
                labels{k} = sprintf('%-13s %s', name, r);
            end
            v = min(max(obj.popDir.Value, 1), numel(labels));
            obj.popDir.String = labels;
            obj.popDir.Value = v;
        end

        function t = target(obj)
            t = obj.targets(obj.popDir.Value);
        end

        function rhoTyped(obj)
            v = str2double(obj.edtRho.String);
            if ~isfinite(v) || v < 0 || v > 120
                EventStation.anonymousWarning('|B| must be a number from 0 to 120 G.');
            else
                obj.rho = v;
                obj.userTarget = true;
            end
            obj.update();
        end

        function dipTyped(obj, k)
            f = str2double(obj.edtDip(k).String);
            t = obj.target();
            r = NaN;
            if isfinite(f); r = Helmholtz.rhoForDip(t.u, obj.helm.D, obj.DIPS{k}, f); end
            if isnan(r)
                EventStation.anonymousWarning('%s', sprintf(['No |B| up to 120 G along %s puts ' ...
                    'that dip at %s MHz.'], t.label, obj.edtDip(k).String));
            else
                obj.rho = r;
                obj.userTarget = true;
            end
            obj.update();
        end

        function update(obj)
            % Fill every box from (direction, |B|), work out reachability, set the preview.
            t = obj.target();
            D = obj.helm.D;
            d = Helmholtz.dipsAlong(t.u, obj.rho, D);
            vals = [d.outer, d.inner];
            obj.edtRho.String = sprintf('%.3f', obj.rho);
            for k = 1:4; obj.edtDip(k).String = sprintf('%.2f', vals(k)); end
            isEqual = strcmp(t.kind, 'equal');
            set(obj.edtDip(3:4), 'Enable', BooleanHelper.boolToOnOff(~isEqual));

            Bt = obj.rho * t.u;
            exact = obj.rho >= t.rhoMin - 1e-9 && obj.rho <= t.rhoMax + 1e-9;
            ang = 0; Bset = Bt;
            if ~exact
                [Bset, ~, ang] = obj.helm.closestReachable(Bt);
            end
            ok = exact || ang <= obj.MAX_DEG;
            obj.btnGo.Enable = BooleanHelper.boolToOnOff(ok);
            if ok
                obj.edtRho.ForegroundColor = [0 0 0];
            else
                obj.edtRho.ForegroundColor = [0.8 0 0];
            end

            lines = {};
            if exact
                lines{end+1} = sprintf('%.2f G along %s is reachable exactly.', obj.rho, t.label);
            elseif ok
                lines{end+1} = sprintf(['Not exactly reachable: Go sets the closest field, ' ...
                    '%.2f G, %.2f deg off.'], norm(Bset), ang);
            else
                lines{end+1} = sprintf('Out of reach: the closest field is %.1f deg off.', ang);
            end
            if isEqual
                lines{end+1} = sprintf('All 4 orientations equal: one pair, %.2f / %.2f MHz.', d.outer);
            end
            lines{end+1} = sprintf('Present field: |B| = %.2f G, %.2f deg from this direction.', ...
                obj.helm.Brho, t.angleToB);
            if ~obj.helm.isCalibrated
                lines{end+1} = 'UNCALIBRATED: directions are in the nominal coil frame.';
                obj.lblStatus.ForegroundColor = [0.8 0 0];
            else
                obj.lblStatus.ForegroundColor = [0 0 0];
            end
            obj.lblStatus.String = lines;
            obj.drawBar(t, ok);

            % the preview: only when it differs from the field already set
            if norm(Bset - obj.helm.B) > 0.01
                obj.helm.setPreview(Bset);
            elseif ~isempty(obj.helm.previewB)
                obj.helm.setPreview([]);
            end
        end

        function drawBar(obj, t, ok)
            % The direction's exact reach (green band), the target |B| (purple) and the
            % present |B| if the field already points this way (orange).
            C = Helmholtz.plotColours();
            ax = obj.axBar;
            cla(ax); hold(ax, 'on');
            his = [obj.targets.rhoMax];
            xmax = max([his(isfinite(his)), obj.rho, obj.helm.Brho, 10]) * 1.05;
            plot(ax, [0 xmax], [0 0], '-', 'Color', C.band, 'LineWidth', 6);
            if ~isnan(t.rhoMin)
                plot(ax, [t.rhoMin t.rhoMax], [0 0], '-', 'Color', C.ok, 'LineWidth', 6);
            end
            if t.angleToB < 0.5
                plot(ax, obj.helm.Brho, 0, 'o', 'MarkerSize', 8, 'LineWidth', 2, ...
                     'Color', C.s1, 'MarkerFaceColor', C.s1);
            end
            mk = C.s2; if ~ok; mk = [0.8 0 0]; end
            plot(ax, obj.rho, 0, 'v', 'MarkerSize', 10, 'Color', mk, 'MarkerFaceColor', mk);
            hold(ax, 'off');
            set(ax, 'XLim', [0 xmax], 'YLim', [-1 1], 'YColor', 'none', 'Color', 'none', ...
                'XColor', C.ink2, 'FontSize', 8, 'TickDir', 'out', 'Box', 'off');
            xlabel(ax, '|B| (G):  green = reachable exactly,  \nabla target,  \bullet present', ...
                   'FontSize', 8, 'Color', C.ink2);
        end

        function goCallback(obj)
            t = obj.target();
            if obj.rho >= t.rhoMin - 1e-9 && obj.rho <= t.rhoMax + 1e-9
                obj.helm.setB(obj.rho * t.u);
            else
                obj.helm.setDirection(t.u, obj.rho, obj.MAX_DEG);
            end
            obj.helm.setPreview([]);
            obj.userTarget = false;
            obj.listTargets();
            obj.update();
        end
    end

    methods
        function onEvent(obj, event)
            if event.isError; return; end
            if isfield(event.extraInfo, Helmholtz.EVENT_CALIBRATION_CHANGED) || ...
                    (isfield(event.extraInfo, Helmholtz.EVENT_MAGNETIC_FIELD_CHANGED) && ~obj.userTarget)
                % nothing typed: follow the present field (a load, the arrows, a script)
                obj.resetToPresent();
            elseif isfield(event.extraInfo, Helmholtz.EVENT_MAGNETIC_FIELD_CHANGED)
                % keep what was typed; refresh the angles, the reach and the status
                obj.listTargets();
                obj.update();
            end
        end
    end
end
