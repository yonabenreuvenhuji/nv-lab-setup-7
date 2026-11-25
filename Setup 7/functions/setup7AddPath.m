function added = setup7AddPath(opts)
%SETUP7ADDPATH Put every Setup 7 folder on the MATLAB path. Run once per session.
%
%   setup7AddPath                                  % add them, report what happened
%   setup7AddPath(struct('save', true))            % and make it permanent (savepath)
%   setup7AddPath(struct('check', true))           % report only, change nothing
%
% The tree is now nine folders deep in places, and a calibration in `calibrations\level 3\`
% calls `esrMultiLineFit` in `functions\`, which means folders have to be added rather than
% just cd'd into. Same derivation as setup7Paths: this file lives in Setup 7\functions\, so
% two levels up is the root and no drive letter is written down anywhere.
%
% ===================================================================
% WHY THIS ALSO CHECKS FOR SHADOWING
% ===================================================================
%
% Adding folders to the path is easy; the failure mode is that a name already exists
% somewhere else on it. MATLAB then silently runs whichever copy comes first, and the
% symptom is not an error -- it is an edit that appears to have no effect, or a fixed bug
% that comes back. Two live examples in this tree:
%
%   ratioAcrossAverages.m   was byte-identical in prod\Util\ and
%                           prod\Experiment\Experiments\
%   extractESRdata.m        differs from extractESRData.m only in one letter's case, and
%                           Windows does not distinguish them
%
% So after adding, this looks up every function Setup 7 owns and reports any that resolve
% to more than one file. A clean run prints nothing about shadowing.
%
% ===================================================================
% OPTIONS
% ===================================================================
%
%   save    call savepath so it survives restarting MATLAB (default false). Off by default
%           because savepath writes to a shared pathdef.m, which on a lab machine other
%           people also use.
%   check   report only, add nothing (default false)
%   quiet   suppress the folder list, keep the warnings (default false)

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('save', false, 'check', false, 'quiet', false), ...
                        'setup7AddPath');

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    if ~isfolder(fullfile(root, 'calibrations'))
        error(['setup7AddPath: no calibrations folder under\n    %s\n' ...
               'This file derives the Setup 7 root from its own location, so it has to ' ...
               'live in Setup 7\\functions\\. It is currently in\n    %s'], root, here);
    end

    % The order is deliberate: functions FIRST, so that if a name in it ever collides with
    % something in a level folder, the library copy wins. Nothing should collide, and the
    % shadow check below says so, but the ordering makes the intended winner explicit.
    wanted = {fullfile(root, 'functions'), ...
              fullfile(root, 'experiments'), ...
              fullfile(root, 'experiments', '_diagnostics'), ...
              fullfile(root, 'calibrations'), ...
              fullfile(root, 'calibrations', 'level 0'), ...
              fullfile(root, 'calibrations', 'level 1'), ...
              fullfile(root, 'calibrations', 'level 2'), ...
              fullfile(root, 'calibrations', 'level 3'), ...
              fullfile(root, 'calibrations', 'level 4')};

    % 'old files' and 'calibrations\results' are deliberately NOT on the path. The archive
    % holds superseded copies of things that are still on the path under their old names --
    % putting it on the path is how a deleted calibration comes back to life.
    added = {};
    p = strsplit(path, pathsep);
    for i = 1:numel(wanted)
        w = wanted{i};
        if ~isfolder(w)
            fprintf('  MISSING  %s\n', w);
            continue
        end
        onPath = any(strcmpi(p, w));
        if onPath
            if ~opts.quiet; fprintf('  already  %s\n', shortName(w, root)); end
        elseif opts.check
            fprintf('  WOULD ADD  %s\n', shortName(w, root));
        else
            addpath(w);
            added{end+1} = w;                                               %#ok<AGROW>
            if ~opts.quiet; fprintf('  added    %s\n', shortName(w, root)); end
        end
    end

    if opts.save && ~opts.check && ~isempty(added)
        try
            savepath;
            fprintf('  path saved (savepath). It will survive restarting MATLAB.\n');
        catch err
            fprintf(['  savepath FAILED (%s). The folders are on the path for THIS session\n' ...
                     '  only. That usually means pathdef.m is read-only or shared.\n'], ...
                    err.message);
        end
    end

    checkShadows(root);
end

% ======================================================================

function checkShadows(root)
% Every function name Setup 7 owns, looked up to see whether more than one file answers.
    names = {};
    folders = {fullfile(root, 'functions'), fullfile(root, 'experiments'), ...
               fullfile(root, 'experiments', '_diagnostics')};
    for lev = 0:4
        folders{end+1} = fullfile(root, 'calibrations', sprintf('level %d', lev)); %#ok<AGROW>
    end
    for i = 1:numel(folders)
        if ~isfolder(folders{i}); continue; end
        d = dir(fullfile(folders{i}, '*.m'));
        for k = 1:numel(d)
            [~, n] = fileparts(d(k).name);
            names{end+1} = n;                                               %#ok<AGROW>
        end
    end

    % which(...,'-all') is MATLAB only -- Octave implements which but not the flag, and
    % returns just the first hit, which would make this report "no shadowing" without having
    % looked. Detected and said out loud rather than passing silently.
    canCheck = true;
    try
        w = which('setup7Paths', '-all');
        canCheck = iscell(w);
    catch
        canCheck = false;
    end
    if ~canCheck
        fprintf(['\nShadow check SKIPPED: this MATLAB/Octave does not support ' ...
                 'which(...,''-all''),\nso duplicate function names cannot be detected here.\n']);
        return
    end

    shadowed = {};
    for i = 1:numel(names)
        try
            w = which(names{i}, '-all');
        catch
            continue
        end
        if ~iscell(w); w = {w}; end
        % Drop MATLAB's own "method of class" and "built-in" entries -- those are not file
        % shadows, and reporting them would bury the ones that are.
        w = w(cellfun(@(s) ischar(s) && endsWith(lower(s), '.m'), w));
        if numel(w) > 1
            shadowed{end+1} = struct('name', names{i}, 'where', {w});       %#ok<AGROW>
        end
    end

    if isempty(shadowed)
        fprintf('\nNo shadowing: every Setup 7 function name resolves to exactly one file.\n');
        return
    end

    fprintf(['\n*** %d SETUP 7 NAME(S) RESOLVE TO MORE THAN ONE FILE. MATLAB runs whichever\n' ...
             'comes first on the path, so an edit to the other copy has no effect and a\n' ...
             'fixed bug can come back. Delete or rename all but one of each.\n'], ...
            numel(shadowed));
    for i = 1:numel(shadowed)
        fprintf('\n  %s\n', shadowed{i}.name);
        for k = 1:numel(shadowed{i}.where)
            marker = '   ';
            if k == 1; marker = ' ->'; end     % the one that actually runs
            fprintf('  %s %s\n', marker, shadowed{i}.where{k});
        end
    end
end

function s = shortName(full, root)
    s = strrep(full, root, '<Setup 7>');
end

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d); given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('%s: unknown option "%s". Valid: %s.', who, bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || isempty(o.(known{i}))
            o.(known{i}) = d.(known{i});
        end
    end
end
