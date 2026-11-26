function p = setup7Paths(which)
%SETUP7PATHS The one place that knows where Setup 7 keeps things.
%
%   p = setup7Paths                  % the NEWEST dated results folder
%   p = setup7Paths('2026-08-14')    % a named one
%   p = setup7Paths('list')          % print every results folder and return them
%
% Lives in Setup 7\functions\, and that is not incidental: the Setup 7 root is derived
% from this file's own location, so there is no drive letter and no folder name written
% down anywhere. It works whether the tree is reached as G:\My Drive\NV Lab\Setup 7 or
% through the .shortcut-targets-by-id path, and it keeps working if the whole tree moves.
%
% ===================================================================
% WHY THIS EXISTS
% ===================================================================
%
% Before this, twelve separate files each carried their own copy of
%
%     'G:\My Drive\NV Lab\Setup 7\calibrations\calibration results\20250505 - with switches\expProperties.json'
%
% as a default. Moving the results folder therefore meant finding and editing twelve
% files, and missing one would not raise an error -- it would quietly load the OLD
% calibration into one experiment while everything else used the new one. Two
% calibrations live at once, and the only symptom is data that does not reproduce.
%
% A hardcoded path repeated N times is not N times the work to change. It is one
% opportunity per copy to change it wrongly.
%
% ===================================================================
% OUTPUT
% ===================================================================
%
%   p.root          ...\Setup 7
%   p.functions     ...\Setup 7\functions
%   p.experiments   ...\Setup 7\experiments
%   p.calibrations  ...\Setup 7\calibrations
%   p.results       ...\Setup 7\calibrations\results          (the parent)
%   p.folder        ...\Setup 7\calibrations\results\<name>   (the one selected)
%   p.name          <name>, e.g. '2026-08-14'
%   p.json          ...\<name>\expProperties.json
%   p.helmholtz     the same file as p.json: the Helmholtz state (currents, field, coil
%                   calibration) is the "Helmholtz" block of expProperties.json now.
%                   h.load(p.helmholtz) / h.save(p.helmholtz) work as before.
%   p.helmholtzLegacy ...\<name>\helmholtzState.mat -- older folders; Helmholtz.load falls
%                   back to it when the json has no "Helmholtz" block
%   p.all           every results folder name found, newest last
%
% The dated folder names sort correctly as strings because they are YYYY-MM-DD. If anyone
% ever names one '14-08-2026' the newest-folder logic silently picks the wrong one, so the
% chosen folder is always printed on the first call of a session.
%
% ===================================================================
% USAGE
% ===================================================================
%
%   In a script:      jsonLocation = setup7Paths().json;
%   In a calibration:  the default is '' and it resolves here, so nothing has to pass it.
%   Pinning a run to an older calibration, deliberately:
%                      calibB1Rabi(struct('jsonLocation', setup7Paths('2026-08-14').json))
%
% Base MATLAB only. Reads the filesystem, writes nothing.

    if ~exist('which', 'var'); which = ''; end

    % ---- the root, from this file's location -----------------------------
    % setup7Paths.m sits in Setup 7\functions\, so two levels up is Setup 7. Deriving it
    % rather than writing it down is what makes this survive the tree moving, and it is
    % also the only way one file can be correct for both the My Drive path and the
    % shortcut-targets path that the same folder answers to.
    here = fileparts(mfilename('fullpath'));
    p.root        = fileparts(here);
    p.functions   = fullfile(p.root, 'functions');
    p.experiments = fullfile(p.root, 'experiments');
    p.calibrations = fullfile(p.root, 'calibrations');
    p.results     = fullfile(p.calibrations, 'results');

    if ~isfolder(p.calibrations)
        error(['setup7Paths: expected to find a calibrations folder at\n    %s\n' ...
               'This file derives the Setup 7 root from its OWN location, so it must live ' ...
               'in Setup 7\\functions\\. It is currently in\n    %s\n' ...
               'Move it there, or the derivation is wrong for every caller.'], ...
              p.calibrations, here);
    end
    if ~isfolder(p.results)
        error(['setup7Paths: no results folder at\n    %s\n' ...
               'Calibration outputs (expProperties.json) are expected in dated subfolders ' ...
               'of it.'], p.results);
    end

    % ---- which results folder -------------------------------------------
    d = dir(p.results);
    names = {d([d.isdir]).name};
    names = names(~ismember(names, {'.', '..'}));
    % Only folders that actually contain the json count as results folders. A stray
    % subfolder should not be able to become "the newest calibration" by being named late
    % in the alphabet.
    hasJson = cellfun(@(n) isfile(fullfile(p.results, n, 'expProperties.json')), names);
    p.all = sort(names(hasJson));

    if isempty(p.all)
        error(['setup7Paths: no results folder under\n    %s\n' ...
               'contains an expProperties.json. Folders present: %s'], ...
              p.results, strjoin(sort(names), ', '));
    end

    if strcmpi(which, 'list')
        fprintf('\nresults folders under %s:\n', p.results);
        for i = 1:numel(p.all)
            tag = '';
            if i == numel(p.all); tag = '   <- newest, the default'; end
            fprintf('  %s%s\n', p.all{i}, tag);
        end
        p.name = p.all{end};
    elseif isempty(which)
        p.name = p.all{end};
    else
        j = find(strcmp(p.all, which), 1);
        if isempty(j)
            error(['setup7Paths: no results folder named "%s". Available: %s\n' ...
                   'Pass one of those, or omit the argument for the newest.'], ...
                  which, strjoin(p.all, ', '));
        end
        p.name = p.all{j};
    end

    p.folder    = fullfile(p.results, p.name);
    p.json      = fullfile(p.folder, 'expProperties.json');
    p.helmholtz = p.json;
    p.helmholtzLegacy = fullfile(p.folder, 'helmholtzState.mat');

    % The Helmholtz state is not required to exist -- level 0 to 2 run with the field off
    % and never touch it -- so this warns rather than erroring. A calibration that needs it
    % will fail on its own terms, with a message about the field rather than about paths.
    hasBlock = contains(fileread(p.json), '"Helmholtz"');
    if ~hasBlock && ~isfile(p.helmholtzLegacy)
        warning('setup7Paths:noHelmholtzState', ...
                ['no Helmholtz state in %s (no "Helmholtz" block in expProperties.json, no ' ...
                 'helmholtzState.mat). Levels 0-2 do not need it; levels 3 and 4 do. Save ' ...
                 'one with h.save(setup7Paths().json).'], p.folder);
    end

    % ---- say which one, once per session ---------------------------------
    % Announced rather than silent, because "which calibration was that run against" is a
    % question that comes up weeks later and the answer has to be in the log.
    persistent announced
    if isempty(announced) || ~strcmp(announced, p.folder)
        fprintf('Setup 7 calibration: %s\n', p.folder);
        if numel(p.all) > 1 && isempty(which)
            fprintf('  (newest of %d; setup7Paths(''list'') shows them all)\n', numel(p.all));
        end
        announced = p.folder;
    end
end
