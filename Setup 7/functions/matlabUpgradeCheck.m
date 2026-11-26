function report = matlabUpgradeCheck(opts)
%MATLABUPGRADECHECK Read-only check that setup 7's code is ready for the running MATLAB release.
%
%   matlabUpgradeCheck                               % full check, prints a report
%   r = matlabUpgradeCheck(struct('deps', false))    % skip the ~1 min code analysis
%
% Written for the R2023b -> R2026b upgrade. Run it in the NEW release on the lab PC
% before connecting to any hardware: it touches no instrument and writes no files.
% (If src\ is not on the path yet it adds it for this session only.)
%
% ===================================================================
% WHAT IT CHECKS
% ===================================================================
%
%   1. Code. Runs this release's Code Compatibility Analyzer over every file setup 7
%      runs: the dependency closure of main.m, ExperimentAWG and the AWG experiments,
%      the digitizer calibrations, ExpT1 and everything in Setup 7\. Errors and
%      Warnings are listed; findings in code setup 7 never runs are listed separately
%      as "expected" (see KNOWN below).
%   2. Toolboxes. Checks the eight toolboxes setup 7's own code calls (see NEEDED for
%      what each is used for). Products referenced only by other setups' code are
%      listed for information.
%   3. Classes. Loads every experiment class the GUI lists and the setup-7 hardware
%      classes, and reports any that fail to load.
%   4. Environment. GUI Layout Toolbox (the GUI is built on uix.*; use 2.4 or newer),
%      a VISA backend for visadev (Tektronix AWG), NI-DAQmx (driver DLL, the C header
%      and a C compiler for LoadNIDAQmx's loadlibrary, and every member of the
%      undocumented daq.ni.NIDAQmx class the NiDaq wrappers call), the libximc driver
%      in <repo>\drivers\, the CompuScope MEX, and whether Java is available. R2026b installs no Java; the code then uses its
%      Java-free paths (PulseStreamerHttpNoJava for the Pulse Streamer, a keyHash
%      digest in DataHash), so missing Java is reported for information only.
%
% A clean result means "ready to try on hardware", not "verified on hardware".

    if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'deps'); opts.deps = true; end

    here = fileparts(mfilename('fullpath'));
    s7   = fileparts(here);
    repo = fileparts(s7);
    src  = fullfile(repo, 'src');
    if ~isfolder(src) || ~isfolder(fullfile(s7, 'calibrations'))
        error('matlabUpgradeCheck must live in <repo>\\Setup 7\\functions\\. It is in %s', here);
    end
    if isempty(which('Setup'))
        addpath(genpath(src));
        fprintf('(added %s to the path for this session)\n', src);
    end

    report = struct('release', version('-release'), 'fail', {{}}, 'warn', {{}}, 'expected', {{}}, 'info', {{}});
    fprintf('\n=== Setup 7 upgrade check on MATLAB R%s ===\n', report.release);

    % Findings in files setup 7 does not run (other setups, or code paths it never takes).
    KNOWN = { ...
        'FrequencyGeneratorTektronixAWG.m', 'legacy visa() branch, used only if useVisadev = false'; ...
        'FrequencyGeneratorRigolAWG.m',     'Rigol AWG (not on setup 7)'; ...
        'FrequencyGeneratorSGT100A.m',      'R&S SGT100A (not on setup 7)'; ...
        'SerialControlled.m',               'legacy serial(): Windfreak, laser sources, Arduino (not on setup 7)'; ...
        'ViewStagePanelTiltCorrection.m',   'tilt panel, shown only for tiltable stages (not the Standa)'; ...
        'DataHash.m',                       'Java digest, used only when a JRE is installed; base64 mode never used'; ...
        'urlread2.m',                       'Java HTTP, used only when a JRE is installed'; ...
        'RegionofInteresst.m',              'camera GUI (not on setup 7)'; ...
        'ramseyWalkDiag.m',                 'legend name-value call that already means what it says'};

    % Toolboxes setup 7's own code calls, and what for (from a per-file dependency scan).
    NEEDED = { ...
        'Instrument Control Toolbox',            'visadev for the Tektronix AWG; calibResonator, scopeProbe'; ...
        'Data Acquisition Toolbox',              'NiDaq talks to NI-DAQmx through daq.ni.NIDAQmx (needs the NI-DAQmx support package too)'; ...
        'Statistics and Machine Learning Toolbox', 'range() in FrequencyGeneratorTektronixAWG; fitStd'; ...
        'Curve Fitting Toolbox',                 'fit() in ESRfit, rabiFit, calibLaserSaturation'; ...
        'Signal Processing Toolbox',             'findpeaks in rabiFit'; ...
        'Mapping Toolbox',                       'extractfield in PathHelper, which builds the GUI experiment list at startup'; ...
        'Phased Array System Toolbox',           'rotx/roty/rotz in Helmholtz'; ...
        'Image Processing Toolbox',              'getrect in ImageScanResult (selecting a region of a scan)'};

    %% 1. Code Compatibility Analyzer --------------------------------------------
    files = {};
    if opts.deps
        fprintf('\n[1] Code analysis (dependency scan takes about a minute)...\n');
        s7files = dir(fullfile(s7, '**', '*.m'));
        s7files = s7files(~contains({s7files.folder}, [filesep 'results']) & ...
                          ~contains({s7files.folder}, [filesep 'old files']));
        addpath(strjoin(unique({s7files.folder}), pathsep));
        entry = [{fullfile(src, 'main.m'), fullfile(src, 'Experiment', 'ExperimentAWG.m'), ...
                  fullfile(src, 'Experiment', 'Experiments', 'ExpT1.m')}, ...
                 fullfile({s7files.folder}, {s7files.name})];
        pats = {fullfile(src, 'Experiment', 'Experiments', 'AWG*.m'), ...
                fullfile(src, 'Experiment', 'Experiments', 'AWGEXPDRDOIDS', '*.m'), ...
                fullfile(src, 'Experiment', 'Experiments', 'Calibrations', 'AWG*.m'), ...
                fullfile(src, 'Experiment', 'Experiments', 'Calibrations', 'digitizer calibration', '*.m')};
        for k = 1:numel(pats)
            d = dir(pats{k});
            entry = [entry, fullfile({d.folder}, {d.name})]; %#ok<AGROW>
        end
        [files, products] = matlab.codetools.requiredFilesAndProducts(unique(entry));
        files = unique([files, entry]);
        files = files(endsWith(files, '.m'));
        r = analyzeCodeCompatibility(files).Recommendations;
        fprintf('    %d files analyzed, %d findings (%d info-level, not listed)\n', ...
                numel(files), height(r), sum(string(r.Severity) == "Info"));
        for k = 1:height(r)
            sev = string(r.Severity(k));
            if sev == "Info"; continue; end
            [~, fn, fe] = fileparts(char(r.File(k)));
            line = sprintf('%s %s:%d  %s', sev, [fn fe], r.LineNumber(k), char(r.Description(k)));
            j = find(strcmp(KNOWN(:, 1), [fn fe]), 1);
            if ~isempty(j)
                report.expected{end+1} = sprintf('%s  [%s]', line, KNOWN{j, 2});
            elseif sev == "Error"
                report.fail{end+1} = line;
            else
                report.warn{end+1} = line;
            end
        end

        others = setdiff({products.Name}, [{'MATLAB'}, NEEDED(:, 1)']);
        if ~isempty(others)
            report.info{end+1} = ['Also referenced (other setups'' code, or uncertain name matches): ' ...
                                  strjoin(others, ', ')];
        end
    end

    %% 2. Toolboxes ---------------------------------------------------------------
    fprintf('\n[2] Toolboxes\n');
    v = ver;
    installed = {v.Name};
    for k = 1:size(NEEDED, 1)
        if any(strcmp(installed, NEEDED{k, 1}))
            fprintf('    %-40s installed\n', NEEDED{k, 1});
        else
            report.fail{end+1} = sprintf('%s not installed: %s', NEEDED{k, 1}, NEEDED{k, 2});
        end
    end

    %% 3. Class loading -----------------------------------------------------------
    fprintf('\n[3] Class loading\n');
    folders = {fullfile(src, 'Experiment', 'Experiments'), ...
               fullfile(src, 'Experiment', 'Experiments', 'Calibrations'), ...
               fullfile(src, 'Tracker', 'Trackables')};
    nExp = 0;
    for f = folders
        d = dir(fullfile(f{1}, '*.m'));
        for k = 1:numel(d)
            [~, cls] = fileparts(d(k).name);
            nExp = nExp + 1;
            try
                getConstPropertyfromString(cls, 'NAME');
            catch err
                msg = sprintf('experiment class %s does not load: %s', cls, err.message);
                if any(strcmp(cls, {'ExpESR_CS_scan', 'ExpFluorescenceMagnetScan'}))
                    report.expected{end+1} = [msg '  [pre-existing, not setup 7]'];
                else
                    report.fail{end+1} = msg;
                end
            end
        end
    end
    hw = {'Setup', 'Experiment', 'ExperimentAWG', 'PulseStreamerNewClass', 'PulseStreamer', 'NiDaq', ...
          'PhotoDiodeDigitizerNiDaqControlled', 'Digitizer', 'ClassStanda8SMC4_stages', ...
          'FrequencyGeneratorTektronixAWG', 'FrequencyGeneratorSRS', 'LaserGate', 'AomNiDaqControlled', ...
          'SwitchPgControlled', 'MW', 'Helmholtz', 'powerSupply_E3632A', 'TrackablePosition', ...
          'StageScanner', 'GuiControllerImage', 'ViewPhotodiodeControl'};
    for k = 1:numel(hw)
        try
            meta.class.fromName(hw{k});
        catch err
            report.fail{end+1} = sprintf('class %s does not load: %s', hw{k}, strtok(err.message, newline));
        end
    end
    fprintf('    %d experiment classes and %d setup-7 classes checked\n', nExp, numel(hw));

    %% 4. Environment -------------------------------------------------------------
    fprintf('\n[4] Environment\n');
    if isempty(which('uix.HBox'))
        report.fail{end+1} = ['GUI Layout Toolbox not installed (the GUI is built on uix.*). ' ...
                              'Install version 2.4 or newer from the Add-On Explorer.'];
    else
        glt = '';
        try
            a = matlab.addons.installedAddons;
            glt = char(a.Version(contains(a.Name, 'GUI Layout', 'IgnoreCase', true)));
        catch
        end
        fprintf('    GUI Layout Toolbox %s: %s\n', glt, fileparts(which('uix.HBox')));
        if ~isempty(glt) && ~isMATLABReleaseOlderThan('R2025a') && verLessThan_(glt, '2.4')
            report.warn{end+1} = sprintf(['GUI Layout Toolbox %s is older than 2.4, the first ' ...
                                          'version for R2025a+ graphics. Update it.'], glt);
        end
    end
    try
        list = visadevlist; %#ok<NASGU>
        fprintf('    VISA backend found (visadevlist works)\n');
    catch err
        report.fail{end+1} = ['No VISA backend for visadev (Tektronix AWG): ' strtok(err.message, newline) ...
                              ' Install NI-VISA or Keysight IO Libraries (not both with TekVISA).'];
    end
    % NI-DAQmx is reached two ways: LoadNIDAQmx.m loadlibrary's nicaiu.dll with the
    % nidaqmx.h header (needs a C compiler), and the DAQmx* wrappers call the Data
    % Acquisition Toolbox's undocumented daq.ni.NIDAQmx class.
    cc = mex.getCompilerConfigurations('C', 'Selected');
    if isempty(cc)
        report.fail{end+1} = ['No C compiler selected (run mex -setup C). LoadNIDAQmx.m calls loadlibrary ' ...
                              'with the nidaqmx.h header, which needs one (MinGW-w64 or Visual C++).'];
    else
        fprintf('    C compiler: %s\n', cc(1).Name);
    end
    if ~isfile('C:\Windows\System32\nicaiu.dll')
        report.fail{end+1} = 'NI-DAQmx driver (C:\Windows\System32\nicaiu.dll) not found.';
    end
    loader = fileread(fullfile(src, 'Physical', 'NiDaq', 'nidaqmx', 'LoadNIDAQmx.m'));
    hdr = regexp(loader, '(?m)^\s*headerPath\s*=\s*''([^'']+)''', 'tokens', 'once');
    if ~isempty(hdr) && ~isfile(hdr{1})
        report.fail{end+1} = sprintf('NI-DAQmx C header not found at %s (path set in LoadNIDAQmx.m).', hdr{1});
    end
    used = {};
    nd = dir(fullfile(src, 'Physical', 'NiDaq', '**', '*.m'));
    for k = 1:numel(nd)
        t = regexp(fileread(fullfile(nd(k).folder, nd(k).name)), 'daq\.ni\.NIDAQmx\.(\w+)', 'tokens');
        used = [used, cellfun(@(c) c{1}, t, 'UniformOutput', false)]; %#ok<AGROW>
    end
    used = unique(used);
    mc = meta.class.fromName('daq.ni.NIDAQmx');
    if isempty(mc)
        report.fail{end+1} = sprintf(['daq.ni.NIDAQmx not found (used %d ways by the NiDaq wrappers). Install ' ...
                                      'the Data Acquisition Toolbox and its NI-DAQmx support package.'], numel(used));
    else
        have = [{mc.MethodList.Name}, {mc.PropertyList.Name}];
        missing = setdiff(used, have);
        if isempty(missing)
            fprintf('    daq.ni.NIDAQmx: all %d members the NiDaq wrappers use exist\n', numel(used));
        else
            report.fail{end+1} = sprintf(['daq.ni.NIDAQmx (an undocumented class) lacks %d member(s) the NiDaq ' ...
                                          'wrappers use: %s'], numel(missing), strjoin(missing, ', '));
        end
    end
    ximc = fullfile(repo, 'drivers', 'ximc-2.10.5', 'ximc', 'win64', 'libximc.dll');
    if ~isfile(ximc)
        report.fail{end+1} = sprintf('Standa libximc driver not found at %s (copy ximc-2.10.5 into drivers\\).', ximc);
    end
    if isempty(which('CsMl_Initialize'))
        report.warn{end+1} = 'CompuScope MATLAB SDK (CsMl_*) not on the path; the digitizer needs it.';
    end
    if usejava('jvm')
        report.info{end+1} = 'Java available: the Pulse Streamer uses urlread2 and DataHash uses MD5, as before.';
    else
        report.info{end+1} = ['No Java (normal from R2026b): the Pulse Streamer uses PulseStreamerHttpNoJava ' ...
                              'and DataHash a keyHash digest. To use the original Java paths instead, ' ...
                              'install "MATLAB Support for OpenJDK" and run jenv -clear.'];
    end

    %% Summary --------------------------------------------------------------------
    printList('FAIL - fix before running', report.fail);
    printList('WARN - check', report.warn);
    printList('expected - not on setup 7''s path', report.expected);
    printList('info', report.info);
    fprintf('\n=== %d fail, %d warn, %d expected ===\n', numel(report.fail), numel(report.warn), numel(report.expected));
    if isempty(report.fail)
        fprintf('Ready to try on hardware (start the GUI, then run the calibration levels).\n');
    end
    report.files = files;
end

function tf = verLessThan_(versionText, minimum)
    % '2.3.6' < '2.4' ? Compare numerically, field by field.
    a = str2double(strsplit(versionText, '.'));
    b = str2double(strsplit(minimum, '.'));
    n = max(numel(a), numel(b));
    a(end+1:n) = 0;
    b(end+1:n) = 0;
    k = find(a ~= b, 1);
    tf = ~isempty(k) && a(k) < b(k);
end

function printList(title, items)
    if isempty(items); return; end
    fprintf('\n%s (%d):\n', title, numel(items));
    for k = 1:numel(items)
        fprintf('  - %s\n', items{k});
    end
end
