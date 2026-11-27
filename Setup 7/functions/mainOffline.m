function mainOffline(opts)
%MAINOFFLINE The main setup window on dummy hardware -- no stage, detector or C:\lab needed.
%
%   mainOffline                                  % photodiode (setup 7) on a dummy xyz stage
%   mainOffline(struct('detector', 'counter'))   % a photon counter instead of the photodiode
%
% For trying out the GUI on any PC (it needs only the GUI Layout Toolbox). Nothing is read
% from or written to C:\lab: the configuration is the repo's own Setup 7\jsons\setupInfo.json
% with every device swapped for its dummy, and saves go to a temp folder. The dummy
% photodiode reads a few bright spots at the dummy stage's position, so a scan draws an image.
%
% Run it in a fresh MATLAB, not next to the live setup: the setup's devices are created once
% per session. The override is removed when the window closes.

    if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'detector'); opts.detector = 'photodiode'; end

    here = fileparts(mfilename('fullpath'));
    s7 = fileparts(here);
    repo = fileparts(s7);
    if isempty(which('Setup'))
        addpath(genpath(fullfile(repo, 'src')));
    end
    if isempty(which('fitStyle'))
        addpath(here);
    end
    % Next to the live setup: refuse (its devices exist already, and would be used)
    live = false;
    try
        live = ~isempty(getObjByName('spcm')) && isempty(getappdata(groot, 'mainOfflineActive'));
    catch
    end
    if live
        error('mainOffline:live', ['The setup''s devices already exist in this MATLAB session. ' ...
            'Run mainOffline in a fresh MATLAB.']);
    end

    S = jsondecode(fileread(fullfile(s7, 'jsons', 'setupInfo.json')));
    S.showEvents = false;
    S.niDaq.dummy = true;
    S.pulseGenerator = struct('type', 'dummy', 'channelNames', {S.pulseGenerator.channelNames}, ...
        'channelValues', S.pulseGenerator.channelValues);
    switch lower(opts.detector)
        case 'photodiode'
            S.spcm = struct('classname', 'photodiodeDummy', 'Switch', S.spcm.Switch);
        otherwise
            S.spcm = struct('classname', 'dummy', 'Switch', S.spcm.Switch);
    end
    S.stages = struct('type', 'Dummy', 'name', 'Stage (dummy)', 'axes', 'xyz', 'tilt_available', false);
    for k = 1:numel(S.lasers)
        S.lasers(k).aomControl.classname = 'dummy';
    end
    fg = S.frequencyGenerators(1);
    fg.type = 'AWGdummy';
    fg.name = 'AWGdummy';
    S.frequencyGenerators = fg;
    S.savesFolder = fullfile(tempdir, 'nv-lab-offline-saves', filesep);
    if ~exist(S.savesFolder, 'dir'); mkdir(S.savesFolder); end

    setappdata(groot, 'JsonInfoReaderOverride', S);
    setappdata(groot, 'mainOfflineActive', true);
    fprintf('mainOffline: dummy hardware, saves go to %s\n', S.savesFolder);

    Setup.init;

    % A scan area around the dummy photodiode's bright spots (the dummy stage's range is
    % +-6.5 mm, which would make a first scan useless)
    stage = ClassStage.getStages;
    stage = stage{1};
    p = stage.scanParams;
    p.from = [-40 -40 -10];
    p.to = [40 40 10];
    p.numPoints = [41 41 5];
    p.fixedPos = [0 0 0];
    p.isFixed = [false false true];
    p.pixelTime = 0.002;
    p.autoSave = false;
    stage.sendEventScanParamsChanged;

    GuiControllerImage().start;
    % The override stays for this MATLAB session (objects created from it, like the event
    % station, read it again later). Restart MATLAB to go back to the real setup.
end
