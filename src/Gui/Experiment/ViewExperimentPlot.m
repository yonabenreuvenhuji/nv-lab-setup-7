classdef ViewExperimentPlot < ViewVBox & EventListener
    %VIEWEXPERIMENTPLOT The experiment window: the live plot, its fit, and the run's controls.
    %
    %   +--------------------------------------------+--------------+
    %   |                                            | Fit          |
    %   |   live plot, with the fit's curve and      |  [Fit][Clear]|
    %   |   the fit's numbers as its title           |  options     |
    %   |                                            |  Auto-fit    |
    %   |                                            |  results     |
    %   +--------------------------------------------+--------------+
    %   | averages: n of N, time per average, left   | this average |
    %   | [Start] [Restart] (View) [Open In Figure] [Save] status   |
    %   | [                  Halt Experiment                       ] |
    %
    % THE FIT is the experiment's robust fitter (esrFitRobust, rabiFitRobust, ...), run on what
    % the plot shows: a.fit from a script does the same thing as the Fit button, and Open In
    % Figure is the fitter's own figure. The options are the experiment's (ExpFit, fitSpec);
    % an experiment without one shows no panel.
    %
    % SAVE writes the run now -- over the file it was saved to, else where the autosave would
    % put it -- fits included. The status beside it follows the run: not saved, saved, changed
    % since (a new fit, more averages), or what the autosave will do while running.

    properties (Access = private)
        expName

        vAxes
        progressbarAverages
        progressbarParameters
        btnStartStop
        btnRestart
        btnOpenInFig
        btnSave
        txtSaveStatus
        bgDisplay
        radioNormalDisplay
        radioAlterDisplay
        radioRawDisplay
        btnEmergencyStop

        % fit panel
        hasFitPanel = false
        btnFit
        btnClearFit
        txtFitWhat
        fitControls         % struct array: the spec's control, plus its handle (h)
        popupAutoFit
        txtFitResult
        btnFitCopy
        btnFitWorkspace

        % timing, for the time left
        resumeTic = []
        resumeIter = 0
    end

    properties (Constant, Access = private)
        AUTOFIT_VALUES = {'off', 'end', 'average'}
        AUTOFIT_ITEMS = {'Off', 'When the run ends', 'After every average'}
        FIT_PANEL_WIDTH = 270
        COLOR_SAVED = [0.80 0.93 0.80]
        COLOR_UNSAVED = [1.00 0.86 0.60]
        COLOR_RUNNING = [0.86 0.91 0.98]
        COLOR_EMPTY = [0.92 0.92 0.92]
    end

    methods
        function obj = ViewExperimentPlot(expName, parent, controller)
            padding = 8;
            spacing = 6;
            obj@ViewVBox(parent, controller, padding, spacing);
            obj@EventListener({expName, SaveLoadCatExp.NAME});
            obj.expName = expName;
            exp = getExpByName(expName);

            fig = obj.component;    % for brevity

            %%% Plot, and the fit panel beside it
            boxMain = uix.HBox('Parent', fig, 'Spacing', 8);
            obj.vAxes = axes('Parent', uicontainer('Parent', boxMain), ...
                'NextPlot', 'replacechildren', ...
                'OuterPosition', [-0.02 0 1.06 1], ...
                'UserData', 'experiment');     % To be plotted on by the Experiment
            obj.hasFitPanel = ~isempty(exp.fitSpec());
            if obj.hasFitPanel
                obj.buildFitPanel(boxMain, exp);
                boxMain.Widths = [-1, obj.FIT_PANEL_WIDTH];
            end

            %%% Progress: averages, and the current average's points
            boxProgress = uix.HBox('Parent', fig, 'Spacing', 12);
            obj.progressbarAverages = progressbar(boxProgress, 0, 'Averages');
            obj.progressbarParameters = progressbar(boxProgress, 0, 'This average');
            boxProgress.Widths = [-3, -2];
            axes();  % So as not to accidently overwrite on these axes

            %%% Run, view, figure, save
            hboxControls = uix.HBox('Parent', fig, 'Spacing', 8);
            obj.btnStartStop = ButtonStartStop(hboxControls, 'Start', 'Pause');
                obj.btnStartStop.startCallback = @obj.btnStartCallback;
                obj.btnStartStop.stopCallback  = @obj.btnStopCallback;
            obj.btnRestart = uicontrol(obj.PROP_BUTTON_BIG_BLUE{:}, ...
                'Parent', hboxControls, ...
                'String', 'Restart', ...
                'Callback', @obj.btnRestartCallback);
            wView = obj.buildViewSelector(hboxControls, exp);
            if isMATLABReleaseOlderThan('R2025a')
                openLabel = '<html><center>Open In<br />Figure</center></html>';  % Java-based uicontrols render HTML
            else
                openLabel = sprintf('Open In\nFigure');     % R2025a+ uicontrols show HTML as raw text
            end
            obj.btnOpenInFig = uicontrol(obj.PROP_BUTTON{:}, ...
                'Parent', hboxControls, ...
                'String', openLabel, ...
                'Tooltip', ['The fit''s own figure when the plot has a fit (the robust fitter''s look), ' ...
                            'else a copy of the plot'], ...
                'Callback', @obj.btnOpenInFigureCallback);
            obj.buildSave(hboxControls);
            hboxControls.Widths = [-1, 90, wView, 80, 230];
            hboxControls.MinimumWidths = [150, 90, wView, 80, 230];
            wControls = sum(hboxControls.MinimumWidths) + 4*hboxControls.Spacing + 2*padding;

            obj.btnEmergencyStop = uicontrol(obj.PROP_BUTTON_BIG_RED{:}, ...
                'Parent', fig, ...
                'String', 'Halt Experiment', ...
                'Callback', @obj.btnEmergencyStopCallback);

            fig.Heights = [-1, 42, 50, 36];

            obj.height = 640;
            obj.width = max([wControls, 760 + obj.hasFitPanel * (obj.FIT_PANEL_WIDTH + 8)]);

            exp.addGraphicAxes(obj.vAxes); % So that experiment could plot on it
            exp.plotResults;
            obj.refresh;
            obj.update;
        end
    end

    %% Building
    methods (Access = private)
        function wTotal = buildViewSelector(obj, parent, exp)
            % Normal | the experiment's alternate view | the raw detection windows
            obj.bgDisplay = uibuttongroup( ...
                'Parent', parent, ...
                'Title', 'View', ...
                'SelectionChangedFcn', @obj.callbackRadioSelection);
            strings = {exp.displayType1, Spcm.localiseUnits(exp.displayType2)};   % older experiments say 'kcps'
            tags = {'normal', 'alternate'};
            if exp.hasRawView()
                strings{end+1} = exp.displayTypeRaw;
                tags{end+1} = 'raw';
            end
            x = 8;
            radios = gobjects(1, numel(strings));
            for i = 1:numel(strings)
                w = 28 + 7*numel(strings{i});
                radios(i) = uicontrol(obj.PROP_RADIO{:}, 'Parent', obj.bgDisplay, ...
                    'String', strings{i}, ...
                    'Position', [x, 6, w, 20], ...    % [fromLeft, fromBottom, width, height]
                    'UserData', tags{i});
                x = x + w + 6;
            end
            obj.radioNormalDisplay = radios(1);
            obj.radioAlterDisplay = radios(2);
            if numel(radios) > 2
                obj.radioRawDisplay = radios(3);
            end
            wTotal = x + 4;
        end

        function buildSave(obj, parent)
            box = uix.VBox('Parent', parent, 'Spacing', 2);
            obj.btnSave = uicontrol(obj.PROP_BUTTON{:}, ...
                'Parent', box, ...
                'String', 'Save', ...
                'Tooltip', ['Save this run now, fits included: over the file it was saved to, else ' ...
                            'where the autosave puts it. Right-click: Save as, show the file.'], ...
                'Callback', @(~, ~) obj.saveCallback(''));
            obj.txtSaveStatus = uicontrol(obj.PROP_TEXT_NO_BG{:}, ...
                'Parent', box, ...
                'HorizontalAlignment', 'center', ...
                'String', '');
            box.Heights = [-1, 18];

            figHandle = ancestor(parent, 'figure');
            cm = uicontextmenu(figHandle);
            uimenu(cm, 'Label', 'Save as...', 'Callback', @(~, ~) obj.saveAsCallback());
            uimenu(cm, 'Label', 'Show the file in Explorer', 'Callback', @(~, ~) obj.showFileCallback());
            uimenu(cm, 'Label', 'Copy the file''s path', 'Callback', @(~, ~) obj.copyPathCallback());
            obj.btnSave.UIContextMenu = cm;
            obj.txtSaveStatus.UIContextMenu = cm;
        end

        function buildFitPanel(obj, parent, exp)
            spec = exp.fitSpec();
            % No panel Title: on R2025a+ figures it is drawn over the panel's content, and the
            % buttons ended up under it. The heading is the box's first row instead.
            panel = uipanel('Parent', parent, 'BorderType', 'line', ...
                'HighlightColor', [0.75 0.75 0.75]);
            box = uix.VBox('Parent', panel, 'Padding', 6, 'Spacing', 4);
            uicontrol(obj.PROP_TEXT_NO_BG{:}, 'Parent', box, 'FontSize', 9, 'FontWeight', 'bold', ...
                'String', sprintf('Fit  (%s)', spec.fitter));
            heights = 18;

            row = uix.HBox('Parent', box, 'Spacing', 6);
            obj.btnFit = uicontrol(obj.PROP_BUTTON{:}, 'Parent', row, 'String', 'Fit', ...
                'Tooltip', ['Fit what the plot shows with the robust fitter and draw it. During a ' ...
                            'run, after the current average. The same as a.fit in a script.'], ...
                'Callback', @(~, ~) obj.fitCallback());
            obj.btnClearFit = uicontrol(obj.PROP_BUTTON{:}, 'Parent', row, 'String', 'Clear', ...
                'Tooltip', 'Remove this view''s fit (a.clearFit)', ...
                'Callback', @(~, ~) obj.clearFitCallback());
            row.Widths = [-2, -1];
            heights(end+1) = 30;

            obj.txtFitWhat = uicontrol(obj.PROP_TEXT_NO_BG{:}, 'Parent', box, 'String', '', ...
                'ForegroundColor', [0.35 0.35 0.35]);
            heights(end+1) = 28;

            % The fitter's options: the ones that change the fit, then the figure's
            ctl = spec.controls;
            ctl = [ctl(~strcmp({ctl.scope}, 'figure')), ctl(strcmp({ctl.scope}, 'figure'))];
            ctlOut = ctl;
            [ctlOut.h] = deal([]);
            shownFigureHeader = false;
            for i = 1:numel(ctl)
                c = ctl(i);
                if strcmp(c.scope, 'figure') && ~shownFigureHeader
                    uicontrol(obj.PROP_TEXT_NO_BG{:}, 'Parent', box, ...
                        'String', 'In the fit figure (Open In Figure):', 'ForegroundColor', [0.35 0.35 0.35]);
                    heights(end+1) = 16; %#ok<AGROW>
                    shownFigureHeader = true;
                end
                cb = @(src, ~) obj.fitOptionCallback(src);
                switch c.style
                    case 'check'
                        h = uicontrol(obj.PROP_CHECKBOX{:}, 'FontSize', 9, 'Parent', box, ...
                            'String', c.label, 'Callback', cb);
                    otherwise
                        r = uix.HBox('Parent', box, 'Spacing', 4);
                        uicontrol(obj.PROP_TEXT_NO_BG{:}, 'Parent', r, 'String', c.label, 'FontSize', 9);
                        if strcmp(c.style, 'popup')
                            h = uicontrol(obj.PROP_POPUP{:}, 'FontSize', 9, 'Parent', r, ...
                                'String', c.items, 'Callback', cb);
                        else
                            h = uicontrol(obj.PROP_EDIT_SMALL{:}, 'FontSize', 9, 'Parent', r, ...
                                'String', '', 'Callback', cb);
                        end
                        r.Widths = [105, -1];
                end
                set(h, 'Tooltip', c.tip, 'UserData', i);
                ctlOut(i).h = h;
                heights(end+1) = 22; %#ok<AGROW>
            end
            obj.fitControls = ctlOut;

            r = uix.HBox('Parent', box, 'Spacing', 4);
            uicontrol(obj.PROP_TEXT_NO_BG{:}, 'Parent', r, 'String', 'Auto-fit', 'FontSize', 9, ...
                'FontWeight', 'bold');
            obj.popupAutoFit = uicontrol(obj.PROP_POPUP{:}, 'FontSize', 9, 'Parent', r, ...
                'String', obj.AUTOFIT_ITEMS, ...
                'Tooltip', ['Fit without pressing Fit. "When the run ends" fits before the ' ...
                            'autosave, so the saved file and its .png have the fit (a.autoFit).'], ...
                'Callback', @(src, ~) obj.autoFitCallback(src));
            r.Widths = [105, -1];
            heights(end+1) = 24;

            obj.txtFitResult = uicontrol('Style', 'edit', 'Parent', box, 'Max', 2, 'Min', 0, ...
                'Enable', 'inactive', 'HorizontalAlignment', 'left', ...
                'FontName', 'FixedWidth', 'FontSize', 8, 'BackgroundColor', [0.985 0.985 0.98], ...
                'String', '');
            heights(end+1) = -1;

            row = uix.HBox('Parent', box, 'Spacing', 6);
            obj.btnFitCopy = uicontrol(obj.PROP_BUTTON_SMALL{:}, 'Parent', row, 'String', 'Copy', ...
                'Tooltip', 'Copy the results above to the clipboard', ...
                'Callback', @(~, ~) obj.copyFitCallback());
            obj.btnFitWorkspace = uicontrol(obj.PROP_BUTTON_SMALL{:}, 'Parent', row, ...
                'String', 'To workspace (o)', ...
                'Tooltip', ['Put the fit''s out struct in the base workspace as o -- the same o as ' ...
                            'o = esrFitRobust(...) in a script. Overwrites any o there.'], ...
                'Callback', @(~, ~) obj.fitToWorkspaceCallback());
            row.Widths = [-1, -2];
            heights(end+1) = 24;

            box.Heights = heights;
        end
    end

    %% Callbacks
    methods (Access = private)
        function btnOpenInFigureCallback(obj, ~, ~)
            exp = getObjByName(obj.expName);
            if isempty(exp)
                EventStation.anonymousWarning('%s Experiment does not exist!', obj.expName)
                return
            end
            view = exp.currentView();
            rec = exp.fitOf(view);
            if ~isempty(rec) && numel(rec.traces) == 1 && exp.canFit() && ~exp.isRunning
                % The fitter's own figure: the same call as the fit on the plot, with plotIt
                % on. With several traces (raw windows) a copy of the plot, which has every
                % trace and its fit on one axes, compares them better than one figure each.
                try
                    exp.fit(struct('plotIt', true, 'verbose', false));
                    return
                catch err
                    EventStation.anonymousWarning('The fit figure failed (%s). Showing a copy of the plot.', err.message);
                end
            end
            isVisible = true;
            fig = AxesHelper.copyToNewFigure(exp.gAxes, isVisible);
            ExpPlotStyle.styleFigure(fig);
            expNameFixed = strrep(obj.expName, '_', ' ');
            set(fig, 'Name', expNameFixed, 'NumberTitle', 'off');
            if isempty(rec)
                % With a fit, the title is the fit's numbers; keep them.
                ax = findobj(fig, 'Type', 'axes', 'UserData', 'experiment');   % the plot, not a top/right axis
                if isempty(ax); ax = gca; end
                ExpPlotStyle.styleTitle(ax(1), expNameFixed);
            end
        end

        function btnStopCallback(obj, ~, ~)
            exp = getExpByName(obj.expName);
            if isempty(exp)
                EventStation.anonymousWarning('%s Experiment does not exist!', obj.expName)
                return
            end
            exp.pause;
            obj.refresh;
            drawnow
        end

        function btnStartCallback(obj, ~, ~)
            exp = getExpByName(obj.expName);
            exp.run;
            obj.refresh;
        end

        function btnRestartCallback(obj, ~, ~)
            strQuestion = sprintf('Do you want to restart the Experiment?');
            strTitle = 'Restart Experiment';
            if QuestionUserYesNo(strTitle, strQuestion)
                exp = getExpByName(obj.expName);
                exp.restart();
                obj.refresh;
            end
        end

        function callbackRadioSelection(obj, ~, event)
            exp = getExpByName(obj.expName);
            view = event.NewValue.UserData;
            exp.isPlotAlternate = strcmp(view, 'alternate');
            exp.isPlotRaw = strcmp(view, 'raw');
            exp.plotResultsWhenSwitchingViews;
            obj.refreshFit;
            drawnow
        end

        function btnEmergencyStopCallback(obj, ~, ~)
            exp = getObjByName(obj.expName);
            if isempty(exp)
                EventStation.anonymousWarning('%s Experiment does not exist!', obj.expName)
                return
            end
            if exp.isRunning
                exp.emergencyStop;
            elseif exp.currParamIter > 0
                exp.clearUncompleteRun;
            end
        end

        %%% Fit %%%
        function fitCallback(obj)
            exp = getObjByName(obj.expName);
            set(obj.btnFit, 'String', 'Fitting...', 'Enable', 'off');
            drawnow
            try
                exp.requestFit();
            catch err
                % Shown in the panel (exp.fitError); the command window gets the stack.
                fprintf(2, '%s fit failed: %s\n', obj.expName, err.message);
            end
            obj.refreshFit;
        end

        function clearFitCallback(obj)
            exp = getObjByName(obj.expName);
            exp.clearFit();
        end

        function fitOptionCallback(obj, src)
            % A control changed: store it in exp.fitOptions, and refit if this view has a fit.
            exp = getObjByName(obj.expName);
            c = obj.fitControls(src.UserData);
            opts = exp.fitOptions;
            switch c.style
                case 'popup'
                    v = c.values{src.Value};
                case 'check'
                    v = logical(src.Value);
                    if isequal(v, c.default); v = []; end     % the fitter's default: pass nothing
                case 'edit'
                    v = sscanf(strrep(src.String, ',', ' '), '%f')';
                    if isempty(strtrim(src.String)); v = []; end
            end
            if isempty(v)
                if isfield(opts, c.key); opts = rmfield(opts, c.key); end
            else
                opts.(c.key) = v;
            end
            exp.fitOptions = opts;
            obj.syncFitControls(exp);
            rec = exp.fitOf(exp.currentView());
            if ~isempty(rec) && ~strcmp(c.scope, 'figure') && ~exp.isRunning
                obj.fitCallback();
            end
        end

        function autoFitCallback(obj, src)
            exp = getObjByName(obj.expName);
            exp.autoFit = obj.AUTOFIT_VALUES{src.Value};
        end

        function copyFitCallback(obj)
            txt = obj.txtFitResult.String;
            if iscell(txt); txt = strjoin(txt, newline); end
            clipboard('copy', txt);
        end

        function fitToWorkspaceCallback(obj)
            exp = getObjByName(obj.expName);
            rec = exp.fitOf(exp.currentView());
            if isempty(rec); return; end
            if numel(rec.traces) == 1
                o = rec.traces(1).out;
            else
                o = {rec.traces.out};
            end
            assignin('base', 'o', o);
            fprintf('o = the %s fit of %s (%s view). The saved copy: no function handles, no experiment.\n', ...
                rec.fitter, obj.expName, rec.view);
        end

        %%% Save %%%
        function saveCallback(obj, path)
            exp = getObjByName(obj.expName);
            try
                exp.saveToFile(path);
            catch err
                EventStation.anonymousWarning('Saving %s failed: %s', obj.expName, err.message);
            end
            obj.refreshSave;
        end

        function saveAsCallback(obj)
            exp = getObjByName(obj.expName);
            if exp.currIter == 0; return; end
            start = exp.savedFile;
            if isempty(start); start = SaveLoad.PATH_DEFAULT_SAVE; end
            [file, folder] = uiputfile('*.mat', sprintf('Save %s as', obj.expName), start);
            if isequal(file, 0); return; end
            obj.saveCallback(fullfile(folder, file));
        end

        function showFileCallback(obj)
            exp = getObjByName(obj.expName);
            if ~isempty(exp.savedFile) && isfile(exp.savedFile)
                system(sprintf('explorer /select,"%s" &', exp.savedFile));
            else
                folder = SaveLoad.PATH_DEFAULT_AUTO_SAVE;
                if isfolder(folder); winopen(folder); end
            end
        end

        function copyPathCallback(obj)
            exp = getObjByName(obj.expName);
            clipboard('copy', exp.savedFile);
        end
    end

    %% Updating display
    methods
        function update(obj)
            exp = getObjByName(obj.expName);
            if isempty(exp) || ~isvalid(obj); return; end

            % Averages, with the time per average and the time left while running
            nDone = exp.currIter;  % The number of the current iteration is also the number of iterations done
            frac = nDone / max(exp.averages, 1);
            str = sprintf('%d of %d averages (%.0f%%)', nDone, exp.averages, frac*100);
            if exp.isRunning && ~isempty(obj.resumeTic) && nDone > obj.resumeIter
                per = toc(obj.resumeTic) / (nDone - obj.resumeIter);
                left = per * (exp.averages - nDone);
                str = sprintf('%s   |   %s per average, %s left', str, ...
                    ViewExperimentPlot.duration(per), ViewExperimentPlot.duration(left));
            end
            progressbar(obj.progressbarAverages, frac, str);
            updateParameterProgressbar(obj);
        end

        function updateParameterProgressbar(obj)
            exp = getObjByName(obj.expName);
            if isempty(exp) || ~isvalid(obj); return; end
            nDone = exp.currParamIter;
            nAll = exp.getTotalNumberOfParams;
            frac = nDone / max(nAll, 1);
            str = sprintf('this average: %d of %d %s', nDone, nAll, exp.parameterName);
            progressbar(obj.progressbarParameters, frac, str);
        end

        function refresh(obj)
            if ~isvalid(obj); return; end
            obj.refreshControls;
            obj.refreshSave;
            obj.refreshFit;
        end

        function refreshControls(obj)
            % Start/Stop status
            if ~isvalid(obj); return; end
            exp = getExpByName(obj.expName);
            inTheMiddleOfRun = exp.currParamIter > 0 || exp.currIter > 0;

            % Take care of Start/Stop button
            if ~exp.restartFlag && inTheMiddleOfRun
                obj.btnStartStop.startString = 'Resume';
            else
                obj.btnStartStop.startString = 'Start';
            end
            if exp.stopFlag
                obj.btnStartStop.stopString = 'Pausing...';
            else
                obj.btnStartStop.stopString = 'Pause';
            end

            % Take care of Emergency Stop button
            bEnable = exp.isRunning || (exp.currParamIter > 0 && exp.currParamIter < exp.getTotalNumberOfParams);
            obj.btnEmergencyStop.Enable = BooleanHelper.boolToOnOff(bEnable);
            if exp.isRunning || exp.currParamIter == 0
                obj.btnEmergencyStop.String = 'Halt Experiment';
            else
                obj.btnEmergencyStop.String = 'Clear Last Average';
            end

            % Take care of Restart button
            obj.btnStartStop.isRunning = exp.isRunning;
            obj.btnRestart.Enable = BooleanHelper.boolToOnOff(~exp.isRunning && inTheMiddleOfRun);

            % Display type
            obj.radioAlterDisplay.Enable = BooleanHelper.boolToOnOff(exp.isPlotAlternateAvailable);
            switch exp.currentView()
                case 'raw'
                    selected = obj.radioRawDisplay;
                case 'alternate'
                    selected = obj.radioAlterDisplay;
                otherwise
                    selected = obj.radioNormalDisplay;
            end
            if ~isempty(selected) && ~isequal(obj.bgDisplay.SelectedObject, selected)
                obj.bgDisplay.SelectedObject = selected;
            end
        end

        function refreshSave(obj)
            if ~isvalid(obj) || isempty(obj.txtSaveStatus) || ~isvalid(obj.txtSaveStatus); return; end
            exp = getObjByName(obj.expName);
            if isempty(exp); return; end
            [state, detail] = exp.saveStatus();
            switch state
                case 'saved'
                    [~, f, ext] = fileparts(detail);
                    str = ['Saved: ' f ext];
                    bg = obj.COLOR_SAVED;
                    tip = detail;
                case 'modified'
                    str = 'Changed since saved';
                    bg = obj.COLOR_UNSAVED;
                    tip = detail;
                case 'unsaved'
                    str = 'Not saved';
                    bg = obj.COLOR_UNSAVED;
                    tip = 'This run is not in any file yet. Save writes it to _AutoSave.';
                case 'running'
                    if exp.shouldAutosave
                        str = 'Autosaves at the end';
                        bg = obj.COLOR_RUNNING;
                    else
                        str = 'Autosave OFF: Save after';
                        bg = obj.COLOR_UNSAVED;
                    end
                    tip = detail;
                otherwise
                    str = 'No data';
                    bg = obj.COLOR_EMPTY;
                    tip = detail;
            end
            set(obj.txtSaveStatus, 'String', str, 'BackgroundColor', bg, 'Tooltip', tip);
            canSave = ~exp.isRunning && exp.currIter > 0;
            set(obj.btnSave, 'Enable', BooleanHelper.boolToOnOff(canSave));
            if any(strcmp(state, {'saved'}))
                set(obj.btnSave, 'String', 'Save again');
            else
                set(obj.btnSave, 'String', 'Save');
            end
        end

        function refreshFit(obj)
            if ~obj.hasFitPanel || ~isvalid(obj) || ~isvalid(obj.btnFit); return; end
            exp = getObjByName(obj.expName);
            if isempty(exp); return; end
            spec = exp.fitSpec();
            view = exp.currentView();
            rec = exp.fitOf(view);
            available = ExpFit.isAvailable(spec);
            canFit = available && exp.canFit() && exp.currIter > 0;

            % What a press of Fit will do
            if ~available
                what = sprintf('%s is not on the path: run setup7AddPath.', spec.fitter);
            elseif ~exp.canFit()
                what = 'This experiment draws its own plot: fit it from a script.';
            elseif exp.fitRequested
                what = 'Fit asked for: it runs when this average ends.';
            else
                switch view
                    case 'normal'
                        what = 'Fits the normalised data, as the fitter reads it from the experiment.';
                    case 'alternate'
                        what = sprintf('Fits each trace of "%s" on its own.', Spcm.localiseUnits(exp.displayType2));
                    otherwise
                        what = 'Fits each raw window on its own (counts scaled to ~1).';
                end
            end
            set(obj.txtFitWhat, 'String', what);

            if isempty(rec); fitLabel = 'Fit'; else; fitLabel = 'Refit'; end
            set(obj.btnFit, 'String', fitLabel, 'Enable', BooleanHelper.boolToOnOff(canFit && ~exp.fitRequested));
            set(obj.btnClearFit, 'Enable', BooleanHelper.boolToOnOff(~isempty(rec)));
            set([obj.btnFitCopy, obj.btnFitWorkspace], 'Enable', BooleanHelper.boolToOnOff(~isempty(rec)));

            obj.syncFitControls(exp);
            obj.popupAutoFit.Value = max(1, find(strcmp(obj.AUTOFIT_VALUES, exp.autoFit), 1));

            % The results
            if ~isempty(exp.fitError)
                txt = {'Fit failed:', exp.fitError};
            elseif isempty(rec)
                if exp.currIter == 0
                    txt = {'No data yet.'};
                else
                    txt = {'No fit on this view.'};
                end
            else
                txt = ExpFit.summary(rec, exp.currIter);
            end
            set(obj.txtFitResult, 'String', txt);
        end

        function syncFitControls(obj, exp)
            % The controls show exp.fitOptions, or the fitter's default for an unset option.
            view = exp.currentView();
            opts = exp.fitOptions;
            for i = 1:numel(obj.fitControls)
                c = obj.fitControls(i);
                has = isfield(opts, c.key);
                switch c.style
                    case 'popup'
                        k = 1;
                        if has
                            k = find(cellfun(@(v) isequal(v, opts.(c.key)), c.values), 1);
                            if isempty(k); k = 1; end
                        end
                        c.h.Value = k;
                    case 'check'
                        v = c.default;
                        if has; v = logical(opts.(c.key)); end
                        c.h.Value = double(v);
                    case 'edit'
                        s = '';
                        if has; s = strtrim(sprintf('%g ', opts.(c.key))); end
                        c.h.String = s;
                end
                inScope = ~strcmp(c.scope, 'normal') || strcmp(view, 'normal');
                c.h.Enable = BooleanHelper.boolToOnOff(inScope);
            end
        end
    end

    methods (Static, Access = private)
        function s = duration(sec)
            if ~isfinite(sec); s = '?'; return; end
            sec = round(sec);
            if sec < 60
                s = sprintf('%d s', sec);
            elseif sec < 3600
                s = sprintf('%d:%02d min', floor(sec/60), mod(sec, 60));
            else
                s = sprintf('%d:%02d h', floor(sec/3600), floor(mod(sec, 3600)/60));
            end
        end
    end

    %% overridden from EventListener
    methods
        % When events happen, this function jumps.
        % event is the event sent from the EventSender
        function onEvent(obj, event)
            if ~isvalid(obj); return; end
            info = event.extraInfo;
            if ~isstruct(info); return; end
            if strcmp(event.creator.name, SaveLoadCatExp.NAME)
                obj.refreshSave;    % saved, loaded, deleted
                return
            end
            if isfield(info, Experiment.EVENT_EXP_RESUMED)
                exp = getObjByName(obj.expName);
                obj.resumeTic = tic;
                obj.resumeIter = exp.currIter;
            end
            if isfield(info, Experiment.EVENT_DATA_UPDATED) ...
                    || isfield(info, Experiment.EVENT_EXP_RESUMED) ...
                    || isfield(info, Experiment.EVENT_EXP_PAUSED)
                obj.refresh;
                obj.update;
            elseif isfield(info, Experiment.EVENT_PARAM_ITERATION_DONE)
                obj.refreshControls;
                obj.updateParameterProgressbar;
            end
            if isfield(info, Experiment.EVENT_FIT_UPDATED)
                obj.refreshFit;
            end
            if isfield(info, Experiment.EVENT_SAVE_STATE)
                obj.refreshSave;
            end
        end
    end
end
