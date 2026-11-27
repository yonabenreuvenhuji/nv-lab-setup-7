function [h, rig] = helmholtzOfflineRig(jsonFile)
%HELMHOLTZOFFLINERIG A real Helmholtz object on dummy supplies, for offline graders.
%
%   [h, rig] = helmholtzOfflineRig;                 % nominal model, no window
%   [h, rig] = helmholtzOfflineRig(jsonFile);       % then h.load(jsonFile)
%   ... use h ...
%   clear rig                                       % closes h, removes the stub
%
% Nothing is read from or written to C:\lab: for as long as `rig` exists, a stub
% JsonInfoReader (in a temp folder, first on the path) serves the repo's own
% Setup 7\jsons\setupInfo.json with the three supplies swapped for "type": "dummy".
% REFUSES if a Helmholtz object already exists -- it must never touch a live session.

    live = [];
    try; live = getObjByName('Helmholtz'); catch; end
    if ~isempty(live)
        error(['helmholtzOfflineRig: a Helmholtz object already exists in this session. ' ...
               'Run offline graders in a fresh MATLAB, not alongside the setup.']);
    end
    here = fileparts(mfilename('fullpath'));
    s7 = fileparts(fileparts(here));
    repo = fileparts(s7);
    addpath(genpath(fullfile(repo, 'src')));
    addpath(fullfile(s7, 'functions'));

    S = struct();
    S.flip_direction = struct('available', false);
    ax = 'XYZ'; ratio = [5.382 8.034 13.649]; offset = [-0.373 0.179 2.121];
    for k = 1:3
        S.(ax(k)) = struct('type', 'dummy', 'name', ax(k), 'address', ['DUMMY' ax(k)], ...
            'range', 'P15V', 'max_current', 5.2, 'magnetic_field_ratio', ratio(k), ...
            'magnetic_field_offset', offset(k));
    end
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
    rig = onCleanup(@() teardown(h, stubDir));
    if exist('jsonFile', 'var') && ~isempty(jsonFile)
        h.load(jsonFile);
    end
end

function teardown(h, stubDir)
    try; if isvalid(h); h.close; end; catch; end
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
