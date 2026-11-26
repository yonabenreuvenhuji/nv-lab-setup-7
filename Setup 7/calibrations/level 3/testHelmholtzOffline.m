function out = testHelmholtzOffline(opts)
%TESTHELMHOLTZOFFLINE Drive a real Helmholtz object (and its window) on dummy supplies.
%
%   testHelmholtzOffline                       % with the GUI
%   testHelmholtzOffline(struct('gui', false))
%
% No hardware, and nothing read from or written to C:\lab: for the duration of the test a
% stub JsonInfoReader (written to a temp folder put first on the path, and removed at the
% end) serves the repo's own Setup 7\jsons\setupInfo.json with the three supplies swapped
% for "type": "dummy" ones (powerSupply_dummy). Checks, on the real class:
%
%   - the nominal model before any calibration (B = diag(ratio)*I + offset, coil frame)
%   - a legacy helmholtzState.mat comes back as the CURRENTS that made its field
%   - a coil calibration changes what B is, not the currents
%   - alignToNV puts B exactly on theta 54.7356 / phi 45, and BtoNV snaps to the nearest
%   - out-of-reach settings are refused, and limits() is right at its own edges
%   - save to a json (a temp copy) and load back: same currents, same calibration, and
%     the other keys of the json untouched
%   - with the GUI: every refresh that the above triggers runs without an error
%
% Run it in a FRESH MATLAB (it refuses if a Helmholtz object exists). Leaves no Helmholtz
% object behind. Needs the GUI Layout Toolbox (uix) for the window; gui = false does not.
% On a machine without C:\lab, MATLAB prints EventListener destructor warnings when it
% exits after a gui = true run -- the stub is gone by then; they are harmless.

    if ~exist('opts', 'var'); opts = struct(); end
    if ~isfield(opts, 'gui'); opts.gui = true; end

    % NEVER on top of a live session: it would have to close the real Helmholtz object,
    % i.e. the supplies' connections. Refuse instead.
    live = [];
    try; live = getObjByName('Helmholtz'); catch; end
    if ~isempty(live)
        error(['testHelmholtzOffline: a Helmholtz object already exists in this session. ' ...
               'Run the test in a fresh MATLAB, not alongside the setup.']);
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

    % ---- the stub setup json -------------------------------------------------------
    setupJson = jsondecode(fileread(fullfile(s7, 'jsons', 'setupInfo.json')));
    setupJson.Helmholtz = S;
    setupJson.showEvents = false;
    setappdata(groot, 'JsonInfoReaderStub', setupJson);
    stubDir = [tempname '_jsonstub'];
    mkdir(stubDir);
    writeStub(fullfile(stubDir, 'JsonInfoReader.m'));
    addpath(stubDir, '-begin');
    clear('JsonInfoReader');
    unstub = onCleanup(@() removeStub(stubDir));

    nFail = 0;
    fprintf('\n==================================================================\n');
    fprintf('testHelmholtzOffline -- the Helmholtz class on dummy supplies\n');
    fprintf('==================================================================\n');
    lastwarn('');
    h = Helmholtz(struct('setup', S, 'gui', opts.gui));
    cleanup = onCleanup(@() closeIt(h));

    % ---- nominal model ----------------------------------------------------------
    ok = ~h.isCalibrated && norm(h.B - (ratio .* h.I + offset)) < 1e-12 && ...
         strcmp(h.frameName, 'nominal coil frame');
    nFail = nFail + report('nominal model: B = ratio.*I + offset, coil frame', ok);

    % ---- legacy .mat -------------------------------------------------------------
    % the newest results folder that still has one (2026-11-23's was moved into its json)
    d = dir(fullfile(s7, 'calibrations', 'results', '*', 'helmholtzState.mat'));
    legacy = '';
    if ~isempty(d)
        [~, iNew] = max(cellfun(@(f) datenum(f(end-9:end), 'yyyy-mm-dd'), {d.folder})); %#ok<DATNM>
        legacy = fullfile(d(iNew).folder, d(iNew).name);
    end
    if ~isempty(legacy) && isfile(legacy)
        L = load(legacy);
        h.load(legacy);
        ok = norm(h.B - L.expStruct.B) < 1e-9;
        nFail = nFail + report('legacy helmholtzState.mat: same field via the currents', ok);
    else
        fprintf('  (no legacy helmholtzState.mat in any results folder to test with)\n');
    end

    % ---- a calibration -----------------------------------------------------------
    R = rotAxis([0.3 -0.5 1], 6);
    Mc = R * diag([5.2 8.4 12.0]) + [0 0.1 0.2; -0.1 0 0.15; 0.05 -0.2 0];
    cal = struct('date', '2026-11-27 10:00', 'frame', 'NV', 'matrix', Mc, ...
                 'offset', [0.18 -0.27 0.41], 'D', 2870.4, 'rmsMHz', 0.2, 'nSpectra', 20, ...
                 'currentMin', [0.3 0.3 0.3], 'currentMax', [5 5 5]);
    I0 = h.I;
    h.setCoilCalibration(cal);
    ok = isequal(h.I, I0) && norm(h.B - (Mc * I0' + cal.offset')') < 1e-12 && h.isCalibrated;
    nFail = nFail + report('calibration changes B, not the currents', ok);

    % ---- NV alignment ------------------------------------------------------------
    T = h.nvTable();
    reachable = find(~isnan([T.rhoMin]));
    ok = ~isempty(reachable);
    nFail = nFail + report(sprintf('some NV direction is reachable (%d of 8)', numel(reachable)), ok);
    k1 = find([T.k] == 1);
    if ~isnan(T(k1).rhoMin)
        rho = min(max(20, T(k1).rhoMin + 1), T(k1).rhoMax - 1);
        h.alignToNV([1 1 1], rho);
        ok = abs(h.Btheta - 54.7356103) < 1e-6 && abs(h.Bphi - 45) < 1e-9 && ...
             abs(h.Brho - rho) < 1e-9 && norm(cross(h.B, [1 1 1])) < 1e-9;
        nFail = nFail + report(sprintf('alignToNV([1 1 1], %.1f): theta 54.7356, phi 45 exactly', rho), ok);
        h.setBsphere(rho, 52, 47);
        h.BtoNV();
        ok = abs(h.Btheta - 54.7356103) < 1e-6 && abs(h.Bphi - 45) < 1e-9;
        nFail = nFail + report('BtoNV from 2-3 deg off snaps onto [1 1 1]', ok);
    else
        fprintf('  [1 1 1] is not reachable with this test calibration; skipped\n');
    end

    % ---- limits and refusal --------------------------------------------------------
    allOk = true;
    for w = {'Bx', 'By', 'Bz', 'Brho', 'Btheta', 'Bphi'}
        [lo, hi] = h.limits(w{1});
        if isnan(lo); allOk = false; continue; end
        inside = h.checkLimits(w{1}, [lo, hi], 0);
        span = hi - lo; eps0 = max(1e-3, 1e-4 * abs(span));
        beyond = true;
        if ~(strcmp(w{1}, 'Bphi') && span >= 359.99) && ...
                ~(strcmp(w{1}, 'Btheta') && (lo <= 1e-9 || hi >= 180 - 1e-9)) && ...
                ~(strcmp(w{1}, 'Brho') && lo <= 1e-9)
            beyond = ~h.checkLimits(w{1}, hi + eps0, 0);
        end
        allOk = allOk && inside && beyond;
    end
    nFail = nFail + report('every limit is reachable at its edge and not past it', allOk);
    threw = false;
    try
        h.setB(h.B + [0 0 500]);
    catch
        threw = true;
    end
    nFail = nFail + report('a field out of reach is refused', threw);
    h.setI([2 3 1.5]);
    ok = abs(h.Btheta - acosd(h.B(3)/norm(h.B))) < 1e-9 && abs(h.Bphi - atan2d(h.B(2), h.B(1))) < 1e-9;
    nFail = nFail + report('Btheta = polar (acos Bz/|B|), Bphi = azimuth (atan2 By, Bx)', ok);

    % ---- json round trip -----------------------------------------------------------
    src = fullfile(s7, 'calibrations', 'results', '2026-11-23', 'expProperties.json');
    tmp = [tempname '.json'];
    copyfile(src, tmp);
    J0 = jsondecode(fileread(src));
    Isaved = h.I;
    h.save(tmp);
    h.clearCoilCalibration();
    h.setI([1 1 1]);
    h.load(tmp);
    J1 = jsondecode(fileread(tmp));
    others = setdiff(fieldnames(J0), {'Helmholtz'});
    same = all(cellfun(@(f) isequaln(J0.(f), J1.(f)), others));
    ok = norm(h.I - Isaved) < 1e-9 && h.isCalibrated && norm(h.coilMatrix - Mc) < 1e-9 && same;
    nFail = nFail + report('save to json and load back: currents, calibration, other keys untouched', ok);
    delete(tmp);
    d = dir([tmp '.bak.*']); for i = 1:numel(d); delete(fullfile(d(i).folder, d(i).name)); end

    % ---- no listener blew up -------------------------------------------------------
    [msg, id] = lastwarn;
    ok = ~(contains(msg, 'onEventTask') || contains(id, 'onEvent'));
    if ~ok; fprintf('    last warning: %s\n', msg); end
    nFail = nFail + report(sprintf('no GUI refresh raised an error (gui = %d)', opts.gui), ok);

    fprintf('\n------------------------------------------------------------------\n');
    if nFail == 0; fprintf('ALL CHECKS PASSED\n'); else; fprintf('%d CHECK(S) FAILED\n', nFail); end
    out = struct('nFail', nFail);
end

function n = report(name, ok)
    if ok; s = 'ok'; else; s = 'FAIL'; end
    fprintf('  %-72s %s\n', name, s);
    n = double(~ok);
end

function R = rotAxis(k, deg)
    k = k(:) / norm(k);
    K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    R = eye(3) + sind(deg) * K + (1 - cosd(deg)) * K * K;
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

function removeStub(stubDir)
    try; rmpath(stubDir); catch; end
    try; rmdir(stubDir, 's'); catch; end
    try; rmappdata(groot, 'JsonInfoReaderStub'); catch; end
    clear('JsonInfoReader');
end

function closeIt(h)
    try
        f = findall(0, 'Type', 'figure', 'Name', 'Helmholtz_GUI');
        if ~isempty(f); delete(f); end
    catch
    end
    try; if isvalid(h); h.close; end; catch; end
end
