function [h, rig] = helmholtzOfflineRig(jsonFile, opts)
%HELMHOLTZOFFLINERIG A real Helmholtz object on dummy supplies, for offline graders.
%
%   [h, rig] = helmholtzOfflineRig;                 % nominal model, no window
%   [h, rig] = helmholtzOfflineRig(jsonFile);       % then h.load(jsonFile)
%   [h, rig] = helmholtzOfflineRig(jsonFile, struct('alongsideLive', true))
%   ... use h ...
%   clear rig                                       % closes h, removes the stub
%
% In a fresh MATLAB (no Helmholtz object yet): h is named 'Helmholtz', so the GUI and
% getObjByName see it, and nothing is read from or written to C:\lab -- for as long as `rig`
% exists, a stub JsonInfoReader (in a temp folder, first on the path) serves the repo's own
% Setup 7\jsons\setupInfo.json with the three supplies swapped for "type": "dummy".
%
% Next to the LIVE setup (a Helmholtz object already exists): REFUSES, unless
% alongsideLive is true. Then h is a second object named 'Helmholtz_offlineRig' on its own
% dummy supplies: getObjByName('Helmholtz') and the GUI keep the live one, its events reach
% no listener, and no stub is installed (the live session's JsonInfoReader stays). The live
% object, its supplies and its currents are never touched. Only for code that takes the
% Helmholtz handle as an argument (calibB1Rabi's 'helmholtz' option) -- anything that looks
% the object up by name would find the live one.

    if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'alongsideLive'); opts.alongsideLive = false; end
    RIG_NAME = 'Helmholtz_offlineRig';

    live = [];
    try; live = getObjByName('Helmholtz'); catch; end
    if ~isempty(live) && ~opts.alongsideLive
        error(['helmholtzOfflineRig: a Helmholtz object already exists in this session. ' ...
               'Run this grader in a fresh MATLAB, not alongside the setup.']);
    end
    here = fileparts(mfilename('fullpath'));
    s7 = fileparts(fileparts(here));
    repo = fileparts(s7);
    if isempty(live)
        addpath(genpath(fullfile(repo, 'src')));
        addpath(fullfile(s7, 'functions'));
    end
    % a rig left behind by a grader that crashed before its cleanup ran
    old = [];
    try; old = getObjByName(RIG_NAME); catch; end
    if ~isempty(old); try; old.close; catch; end; end

    S = struct();
    S.flip_direction = struct('available', false);
    ax = 'XYZ'; ratio = [5.382 8.034 13.649]; offset = [-0.373 0.179 2.121];
    addr = 'DUMMY';
    if ~isempty(live); addr = 'OFFLINERIG'; end    % own supply names, even if the live one is dummy
    for k = 1:3
        S.(ax(k)) = struct('type', 'dummy', 'name', ax(k), 'address', [addr ax(k)], ...
            'range', 'P15V', 'max_current', 5.2, 'magnetic_field_ratio', ratio(k), ...
            'magnetic_field_offset', offset(k));
    end

    stubDir = '';
    if isempty(live)
        setupJson = jsondecode(fileread(fullfile(s7, 'jsons', 'setupInfo.json')));
        setupJson.Helmholtz = S;
        setupJson.showEvents = false;
        setappdata(groot, 'JsonInfoReaderStub', setupJson);
        stubDir = [tempname '_jsonstub'];
        mkdir(stubDir);
        writeStub(fullfile(stubDir, 'JsonInfoReader.m'));
        addpath(stubDir, '-begin');
        clear('JsonInfoReader');
        h = Helmholtz(struct('setup', S, 'gui', false));
    else
        h = Helmholtz(struct('setup', S, 'gui', false, 'name', RIG_NAME));
    end
    rig = onCleanup(@() teardown(h, stubDir));
    if exist('jsonFile', 'var') && ~isempty(jsonFile)
        h.load(jsonFile);
    end
end

function teardown(h, stubDir)
    try; if isvalid(h); h.close; end; catch; end
    if isempty(stubDir); return; end         % next to a live session: no stub was installed
    try; rmpath(stubDir); catch; end
    try; rmdir(stubDir, 's'); catch; end
    try; rmappdata(groot, 'JsonInfoReaderStub'); catch; end
    clear('JsonInfoReader');
end

function writeStub(file)
    L = {
'classdef JsonInfoReader'
'    % TEST STUB written by testHelmholtzOffline -- serves getappdata(groot, ''JsonInfoReaderStub'').'
'    properties (Constant)'
'        JSON_LOCATION_DIR = '''''
'        JSON_LOCATION_FILENAME = ''stub'''
'        JSON_PARAMS_LOCATION_FILENAME = ''stub'''
'    end'
'    methods (Static)'
'        function s = getJson()'
'            s = getappdata(groot, ''JsonInfoReaderStub'');'
'            if ~isfield(s, ''debugMode''); s.debugMode = false; end'
'        end'
'        function p = getParams(jsonLocation)'
'            p = jsondecode(fileread(jsonLocation));'
'        end'
'        function setParams(varargin)'
'            error(''JsonInfoReader stub: setParams is not available in the offline test'');'
'        end'
'        function n = setupNumber()'
'            s = JsonInfoReader.getJson();'
'            n = ''7'';'
'            if isfield(s, ''setupNumber''); n = s.setupNumber; if ~ischar(n); n = num2str(n); end; end'
'        end'
'        function o = getDefaultObject(listName, varargin)'
'            s = JsonInfoReader.getJson();'
'            o = s.(listName);'
'            if iscell(o); o = o{1}; elseif numel(o) > 1; o = o(1); end'
'        end'
'    end'
'end'};
    fid = fopen(file, 'w');
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end
