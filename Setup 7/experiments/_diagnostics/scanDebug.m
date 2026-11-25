function out = scanDebug(opts)
%SCANDEBUG  What does the Setup 7 stage scan really measure, and when?
%
% Uses the scan settings in the GUI (Stage - 8SMC4 scan parameters, pixel
% time, photodiode option) and only public methods of the stage, the
% photodiode, the digitizer and the pulse generator. Changes no code and
% no settings; the stage, the green laser, the 'detector' channel and the
% scan acquisition mode are restored at the end (also after an error).
%
% Tests (each one runs even if an earlier one failed):
%   A  Digitizer trigger timing, stage still at the scan centre. One
%      segment per trigger, exactly as a scan line, while the green laser
%      is switched on for odd and off for even triggers. If the segments
%      are taken at the triggers, the values alternate. Done with the
%      scan's channels and with signal & reference (the reference
%      photodiode always sees the laser). Digitizer status is logged
%      before, during and after every trigger. Each trigger is a short
%      pulse that ends while the segment is still being recorded, as in
%      the scan code; holdTrigger = true keeps the line high until the
%      segment has ended (the first version, which gives a second
%      segment per pulse on Setup 7).
%   A2 The same laser on/off pattern with forced single reads
%      (readFromTime): shows whether the scan channel sees the laser at
%      the centre at all.
%   B  Stage moves: stage.checkMoveWait on both scan axes.
%   C  The middle line along the fast axis, four ways:
%      C1 script: move, read back position, settle, one triggered segment
%         per pixel (raw data kept);
%      C2 script: move, read back, settle, forced read per pixel (the
%         method the scan used before the fix);
%      C3 the real scan code (ScanX/ScanZ + photodiode), trigger mode;
%      C4 the real scan code, software mode.
%   D  The middle line along the slow axis, C1 and C2.
%
% Result: printed verdicts, a figure, and scanDebug_<time>.mat next to
% this file (send that file, or its name, for analysis).
%
%   out = scanDebug;                               % GUI settings
%   out = scanDebug(struct('nTrig', 12));
%   out = scanDebug(struct('replay', out));        % re-analyse, no hardware
%
% Options: nTrig (10), laserChannel ('greenLaser'), triggerChannel
% ('detector'), laserSettle (0.03 s), settle (stage.scanSettleTime),
% margin (0.003 s), maxSamples (1e5), tests ('A A2 B C D'), plot (true),
% save (true), tPixel ([] = GUI pixel time), holdTrigger (false).

if nargin < 1; opts = struct(); end
if isfield(opts, 'replay')
    out = analyse(opts.replay);
    printReport(out);
    if getOpt(opts, 'plot', true); plotReport(out); end
    return
end

o.nTrig = getOpt(opts, 'nTrig', 10);
o.laserChannel = getOpt(opts, 'laserChannel', 'greenLaser');
o.triggerChannel = getOpt(opts, 'triggerChannel', 'detector');
o.laserSettle = getOpt(opts, 'laserSettle', 0.03);
o.margin = getOpt(opts, 'margin', 0.003);
o.maxSamples = getOpt(opts, 'maxSamples', 1e5);
o.tests = getOpt(opts, 'tests', 'A A2 B C D');
o.doPlot = getOpt(opts, 'plot', true);
o.doSave = getOpt(opts, 'save', true);
o.tPixel = getOpt(opts, 'tPixel', []);
o.holdTrigger = getOpt(opts, 'holdTrigger', false);

%% objects
ctx.stage = getObjByName('Stage - 8SMC4');
ctx.spcm = getObjByName(Spcm.NAME);
ctx.digi = getObjByName(Digitizer.NAME);
ctx.pg = getObjByName(PulseGenerator.NAME);
if isempty(ctx.stage) || isempty(ctx.spcm) || isempty(ctx.digi) || isempty(ctx.pg)
    error('scanDebug needs the Standa stage, the photodiode, the digitizer and the pulse generator (start the Setup 7 GUI first).');
end
ctx.o = o;

%% geometry, from the GUI scan parameters
sp = ctx.stage.scanParams;
g = geometry(sp);
if ~isempty(o.tPixel); g.tPixel = o.tPixel; end
ctx.g = g;
o.settle = getOpt(opts, 'settle', ctx.stage.scanSettleTime);
ctx.o = o;

out = struct();
out.when = datestr(now, 'yyyy-mm-dd HH:MM:SS');
out.opts = o;
out.geometry = g;
out.info = collectInfo(ctx);
fprintf('\nscanDebug %s\n  fast axis %s: %d points %.2f..%.2f um, slow axis %s: %d points %.2f..%.2f um, fixed %s = %.2f um, pixel time %.4g s\n', ...
    out.when, g.letterA, numel(g.vecA), g.vecA(1), g.vecA(end), g.letterB, numel(g.vecB), g.vecB(1), g.vecB(end), ...
    g.letterC, g.posC, g.tPixel);
fprintf('  photodiode option %d, scan channels %s, scan mode ''%s''\n', out.info.channelOption, mat2str(out.info.scanChannels), out.info.scanMode);

%% state to restore
state.pos = ctx.stage.Pos(ctx.stage.availableAxes);
state.axes = ctx.stage.availableAxes;
state.laserOn = logical(ctx.pg.isOn(o.laserChannel));
state.mode = ctx.spcm.scanAcquisitionMode;
out.startState = state;
cleanupObj = onCleanup(@() restoreState(ctx, state)); %#ok<NASGU>
if ~state.laserOn
    fprintf('  note: %s was off; it is switched on for the tests and off again at the end.\n', o.laserChannel);
end
ctx.spcm.setSPCMEnable(true);
setLaser(ctx, true);
pause(0.3);                                         % photodiode supply / laser settle

%% tests
has = @(name) any(strcmp(strsplit(o.tests), name));
if has('A')
    out.A1 = runSafe(@() testTriggerTiming(ctx, out.info.scanChannels), 'A (scan channels)');
    if ~isequal(out.info.scanChannels, [1 2])
        out.A2ch = runSafe(@() testTriggerTiming(ctx, [1 2]), 'A (signal & reference)');
    end
end
if has('A2')
    out.A2 = runSafe(@() testSoftwareToggle(ctx), 'A2 (forced reads, laser on/off)');
end
if has('B')
    out.B = runSafe(@() testMoves(ctx), 'B (stage moves)');
end
if has('C')
    out.C1 = runSafe(@() testLineScript(ctx, 'A', 'trigger'), 'C1 (fast line, script, triggered)');
    out.C2 = runSafe(@() testLineScript(ctx, 'A', 'software'), 'C2 (fast line, script, forced reads)');
    out.C3 = runSafe(@() testLineClasses(ctx, 'trigger'), 'C3 (fast line, scan code, trigger mode)');
    out.C4 = runSafe(@() testLineClasses(ctx, 'software'), 'C4 (fast line, scan code, software mode)');
end
if has('D')
    out.D1 = runSafe(@() testLineScript(ctx, 'B', 'trigger'), 'D1 (slow line, script, triggered)');
    out.D2 = runSafe(@() testLineScript(ctx, 'B', 'software'), 'D2 (slow line, script, forced reads)');
end
[out.lastWarning, out.lastWarningId] = lastwarn;

%% analyse, save, show
out = analyse(out);
if o.doSave
    folder = fileparts(mfilename('fullpath'));
    out.file = fullfile(folder, sprintf('scanDebug_%s.mat', datestr(now, 'yyyymmdd_HHMMSS')));
    try
        save(out.file, 'out', '-v7');
    catch err
        fprintf('  could not save: %s\n', err.message);
    end
end
printReport(out);
if o.doPlot
    try
        plotReport(out);
    catch err
        fprintf('  (no figure: %s)\n', err.message);
    end
end
end

%% ===================================================================
%% tests
%% ===================================================================
function r = testTriggerTiming(ctx, chans)
% n triggered segments at the scan centre, laser on for odd, off for even.
o = ctx.o; g = ctx.g; digi = ctx.digi; pg = ctx.pg;
n = o.nTrig;
moveTo(ctx, [g.letterA g.letterB], [g.centerA g.centerB]);
r.chans = chans;
r.laser = mod(1:n, 2) == 1;
tUs = g.tPixel * digi.UNIT_CONV;
rate = pickRate(tUs, digi, o.maxSamples);
pg.Off(o.triggerChannel);
setLaser(ctx, true);
digi.prepareAcquisition(chans, n, tUs, 1, rate);
r.config = digiConfig(digi);
r.board = boardConfig(digi, chans);
segS = digi.segmentPostTrigSize / digi.sampleRate / digi.UNIT_CONV;
r.segSeconds = segS;
t0 = tic;
digi.startExperiment;
r.tArmed = toc(t0);
r.preTrigger = pollStatus(digi, 0.1, t0);        % [t status] changes, no trigger sent yet
r.statusBefore = nan(1, n); r.statusAfter = nan(1, n); r.statusMid = nan(1, n);
r.tOnCall = nan(1, n); r.tOnReturn = nan(1, n); r.tOffReturn = nan(1, n);
r.poll = cell(1, n);
for k = 1:n
    setLaser(ctx, r.laser(k));
    pause(o.laserSettle);
    r.statusBefore(k) = gageStatus(digi);
    r.tOnCall(k) = toc(t0);
    pg.On(o.triggerChannel);
    r.tOnReturn(k) = toc(t0);
    if ~o.holdTrigger
        pg.Off(o.triggerChannel);               % short pulse: low again while the segment is recorded
        r.tOffReturn(k) = toc(t0);
    end
    remaining = segS + o.margin - (toc(t0) - r.tOnReturn(k));
    if k <= 3
        r.poll{k} = pollStatus(digi, max(remaining, 0), t0);
    else
        r.statusMid(k) = gageStatus(digi);
        waitSeconds(remaining);
    end
    r.statusAfter(k) = gageStatus(digi);
    if o.holdTrigger
        pg.Off(o.triggerChannel);
        r.tOffReturn(k) = toc(t0);
    end
end
r.holdTrigger = o.holdTrigger;
setLaser(ctx, true);
r.statusEnd = waitNotBusy(digi, 1);
r = readSegments(ctx, r, chans);
end

function r = testSoftwareToggle(ctx)
% Forced reads (readFromTime) with the laser on/off, at the scan centre.
o = ctx.o; g = ctx.g; spcm = ctx.spcm;
n = o.nTrig;
moveTo(ctx, [g.letterA g.letterB], [g.centerA g.centerB]);
r.laser = mod(1:n, 2) == 1;
spcm.prepareReadByTime(g.tPixel);
r.v = []; r.e = [];
for k = 1:n
    setLaser(ctx, r.laser(k));
    pause(o.laserSettle);
    [v, e] = spcm.readFromTime();
    r.v(:, k) = v(:); r.e(:, k) = e(:);
end
setLaser(ctx, true);
spcm.clearTimeRead();
end

function r = testMoves(ctx)
g = ctx.g;
moveTo(ctx, [g.letterA g.letterB], [g.centerA g.centerB]);
r.A = ctx.stage.checkMoveWait(g.letterA, abs(g.pitchA), min(10, numel(g.vecA) - 1));
moveTo(ctx, [g.letterA g.letterB], [g.centerA g.centerB]);
r.B = ctx.stage.checkMoveWait(g.letterB, abs(g.pitchB), min(5, numel(g.vecB) - 1));
end

function r = testLineScript(ctx, which, how)
% Middle line, moved and measured by this script, with full logging.
o = ctx.o; g = ctx.g; stage = ctx.stage; digi = ctx.digi; pg = ctx.pg; spcm = ctx.spcm;
if strcmp(which, 'A')
    L = g.letterA; vec = g.vecA; otherL = g.letterB; otherV = g.centerB;
else
    L = g.letterB; vec = g.vecB; otherL = g.letterA; otherV = g.centerA;
end
n = numel(vec);
r.axis = L; r.vec = vec; r.how = how; r.other = {otherL, otherV};
lead = stage.scanLeadIn;
moveTo(ctx, [otherL L], [otherV, vec(1) - sign(vec(end) - vec(1)) * lead]);
chans = ctx.spcm.activeChannelsGUI;
r.chans = chans;
r.tMoveStart = nan(1, n); r.tMoveEnd = nan(1, n); r.posRead = nan(1, n);
r.statusBefore = nan(1, n); r.statusAfter = nan(1, n); r.tOnReturn = nan(1, n); r.tRead = nan(1, n);
t0 = tic;
if strcmp(how, 'trigger')
    tUs = g.tPixel * digi.UNIT_CONV;
    rate = pickRate(tUs, digi, o.maxSamples);
    pg.Off(o.triggerChannel);
    digi.prepareAcquisition(chans, n, tUs, 1, rate);
    r.config = digiConfig(digi);
    segS = digi.segmentPostTrigSize / digi.sampleRate / digi.UNIT_CONV;
    r.segSeconds = segS;
    digi.startExperiment;
else
    spcm.prepareReadByTime(g.tPixel);
    r.v = nan(1, n);
end
for i = 1:n
    r.tMoveStart(i) = toc(t0);
    stage.Move(L, vec(i));
    r.tMoveEnd(i) = toc(t0);
    r.posRead(i) = stage.Pos(L);
    pause(o.settle);
    if strcmp(how, 'trigger')
        r.statusBefore(i) = gageStatus(digi);
        pg.On(o.triggerChannel);
        r.tOnReturn(i) = toc(t0);
        if ~o.holdTrigger; pg.Off(o.triggerChannel); end
        waitSeconds(segS + o.margin - (toc(t0) - r.tOnReturn(i)));
        r.statusAfter(i) = gageStatus(digi);
        if o.holdTrigger; pg.Off(o.triggerChannel); end
    else
        [v, ~] = spcm.readFromTime();
        r.v(i) = v(1);
        r.tRead(i) = toc(t0);
    end
end
if strcmp(how, 'trigger')
    r.statusEnd = waitNotBusy(digi, 1);
    r = readSegments(ctx, r, chans);
    if isfield(r, 'value'); r.v = r.value(1, :); end
else
    spcm.clearTimeRead();
end
r.tTotal = toc(t0);
end

function r = testLineClasses(ctx, mode)
% Middle fast line through the real scan code (same calls as StageScanner.scan1D).
g = ctx.g; stage = ctx.stage; spcm = ctx.spcm;
L = g.letterA; vec = g.vecA; n = numel(vec);
xyz = num2cell(g.fixedXYZ);
xyz{g.idxB} = g.centerB;
xyz{g.idxA} = vec;
spcm.scanAcquisitionMode = mode;
r.mode = mode; r.vec = vec; r.axis = L;
t0 = tic;
feval(['PrepareScan' upper(L)], stage, xyz{1}, xyz{2}, xyz{3}, 0, 0, g.tPixel);
spcm.prepareCountByStage(stage.name, n, 2 * n * g.tPixel, false, g.tPixel);
spcm.startScanCount();
r.tPrepared = toc(t0);
feval(['Scan' upper(L)], stage, xyz{1}, xyz{2}, xyz{3}, 0, 0, g.tPixel);
r.tScanned = toc(t0);
r.statusBeforeRead = gageStatus(ctx.digi);
try
    [v, e] = spcm.readFromScan();
    r.v = v; r.e = e;
catch err
    r.readError = err.message;
end
r.tRead = toc(t0);
spcm.clearScanRead();
stage.AbortScan();
r.posEnd = stage.Pos(L);
end

%% ===================================================================
%% analysis
%% ===================================================================
function out = analyse(out)
v = struct();
% A: do the segments follow the triggers?
for f = {'A1', 'A2ch'}
    if isfield(out, f{1}) && ~isfield(out.(f{1}), 'error') && isfield(out.(f{1}), 'segMean')
        r = out.(f{1});
        rr = struct();
        rr.chans = r.chans;
        for c = 1:size(r.segMean, 1)
            m = r.segMean(c, :);
            rr.z(c) = onOffZ(m, r.laser);
            rr.spreadRel(c) = std(m) / max(abs(mean(m)), eps);
        end
        rr.preTriggerStatuses = unique(r.preTrigger(:, 2))';
        rr.statusBefore = r.statusBefore; rr.statusAfter = r.statusAfter; rr.statusEnd = r.statusEnd;
        rr.earlyComplete = find(r.statusAfter(1:end-1) == 0, 1);
        v.(f{1}) = rr;
    end
end
if isfield(out, 'A2') && ~isfield(out.A2, 'error')
    for c = 1:size(out.A2.v, 1)
        v.A2.z(c) = onOffZ(out.A2.v(c, :), out.A2.laser);
    end
end
% lines
lines = {'C1', 'C2', 'C3', 'C4', 'D1', 'D2'};
for k = 1:numel(lines)
    f = lines{k};
    if isfield(out, f) && ~isfield(out.(f), 'error') && isfield(out.(f), 'v')
        y = out.(f).v(:)';
        lv.structure = (max(y) - min(y)) / max(noiseOf(y), eps);
        lv.range = [min(y) max(y)];
        if isfield(out.(f), 'posRead')
            lv.maxPosErr = max(abs(out.(f).posRead - out.(f).vec));
        end
        if isfield(out.(f), 'statusAfter')
            lv.earlyComplete = find(out.(f).statusAfter(1:end-1) == 0, 1);
        end
        v.(f) = lv;
        clear lv
    end
end
pairs = {'C1', 'C2'; 'C3', 'C2'; 'C4', 'C2'; 'C3', 'C4'; 'D1', 'D2'};
for k = 1:size(pairs, 1)
    a = pairs{k, 1}; b = pairs{k, 2};
    if isfield(v, a) && isfield(v, b)
        ya = out.(a).v(:); yb = out.(b).v(:);
        if numel(ya) == numel(yb) && std(ya) > 0 && std(yb) > 0
            cc = corrcoef(ya, yb);
            v.corr.([a '_' b]) = cc(1, 2);
        end
    end
end
out.verdict = v;
end

function printReport(out)
v = out.verdict;
fprintf('\n================ scanDebug report ================\n');
names = {'A1', 'A2ch'};
for k = 1:2
    f = names{k};
    if isfield(out, f)
        if isfield(out.(f), 'error')
            fprintf('A  %s: ERROR %s\n', mat2str(out.(f).chans), out.(f).error);
        elseif isfield(v, f)
            r = v.(f);
            fprintf('A  channels %s: laser on/off contrast z = %s (|z| > 5: segments follow the triggers)\n', mat2str(r.chans), mat2str(round(r.z * 10) / 10));
            fprintf('   status before any trigger: %s (expect 1 = waiting); status after each segment: %s (must not be 0 before the last)\n', ...
                mat2str(r.preTriggerStatuses), mat2str(r.statusAfter));
            if ~isempty(r.earlyComplete)
                fprintf('   >>> digitizer reported COMPLETE after trigger %d of %d: it does not wait for the triggers\n', r.earlyComplete, numel(r.statusAfter));
            end
        end
    end
end
if isfield(out, 'A2')
    if isfield(out.A2, 'error')
        fprintf('A2 ERROR %s\n', out.A2.error);
    else
        fprintf('A2 forced reads, laser on/off contrast z = %s (does the scan channel see the laser here?)\n', mat2str(round(v.A2.z * 10) / 10));
    end
end
if isfield(out, 'B')
    if isfield(out.B, 'error')
        fprintf('B  ERROR %s\n', out.B.error);
    else
        for ax = {'A', 'B'}
            rep = out.B.(ax{1});
            fprintf('B  axis %s: new wait %.1f ms median, old wait %.1f ms, stage short at old "done": %d of %d moves\n', ...
                ax{1}, 1e3 * median([rep.tNewDone]), 1e3 * median([rep.tOldDone]), sum([rep.shortAtOldDoneUm] > 0.5), numel(rep));
        end
    end
end
lines = {'C1', 'C2', 'C3', 'C4', 'D1', 'D2'};
for k = 1:numel(lines)
    f = lines{k};
    if ~isfield(out, f); continue; end
    if isfield(out.(f), 'error')
        fprintf('%s ERROR %s\n', f, out.(f).error);
    elseif isfield(out.(f), 'readError')
        fprintf('%s read ERROR %s\n', f, out.(f).readError);
    elseif isfield(v, f)
        s = sprintf('%s structure %.1f (range %.4g..%.4g)', f, v.(f).structure, v.(f).range(1), v.(f).range(2));
        if isfield(v.(f), 'maxPosErr'); s = [s sprintf(', max position error %.3f um', v.(f).maxPosErr)]; end
        if isfield(v.(f), 'earlyComplete') && ~isempty(v.(f).earlyComplete)
            s = [s sprintf(', DIGITIZER COMPLETE AFTER PIXEL %d', v.(f).earlyComplete)];
        end
        fprintf('%s\n', s);
    end
end
if isfield(v, 'corr')
    fn = fieldnames(v.corr);
    for k = 1:numel(fn)
        fprintf('   correlation %s: %.3f\n', strrep(fn{k}, '_', ' vs '), v.corr.(fn{k}));
    end
end
fprintf('--------------------------------------------------\n');
yn = {'NO', 'yes'};
if isfield(v, 'A2ch') || isfield(v, 'A1')
    if isfield(v, 'A2ch'); r = v.A2ch; zz = r.z(end); else; r = v.A1; zz = max(abs(r.z)); end
    fprintf('Digitizer takes one segment per trigger (laser on/off visible):  %s\n', yn{1 + (abs(zz) > 5 && isempty(r.earlyComplete))});
end
if isfield(v, 'corr')
    show = {'C1_C2', 'Triggered line by this script matches forced reads:'; ...
            'C3_C2', 'Scan code, trigger mode, matches forced reads:      '; ...
            'C4_C2', 'Scan code, software mode, matches forced reads:     '; ...
            'D1_D2', 'Slow axis: triggered line matches forced reads:     '};
    for k = 1:size(show, 1)
        if isfield(v.corr, show{k, 1})
            fprintf('%s %s (r = %.3f)\n', show{k, 2}, yn{1 + (v.corr.(show{k, 1}) > 0.9)}, v.corr.(show{k, 1}));
        end
    end
end
if isfield(v, 'C2')
    fprintf('Forced-read line shows structure (range/noise %.0f): %s\n', v.C2.structure, yn{1 + (v.C2.structure > 10)});
end
if isfield(out, 'file'); fprintf('saved: %s\n', out.file); end
fprintf('==================================================\n');
end

function plotReport(out)
figure('Name', 'scanDebug', 'Color', 'w');
nr = 2; nc = 3; p = 0;
for f = {'A1', 'A2ch'}
    p = p + 1;
    subplot(nr, nc, p);
    if isfield(out, f{1}) && isfield(out.(f{1}), 'segMean')
        r = out.(f{1});
        m = r.segMean ./ max(abs(mean(r.segMean, 2)), eps);
        plot(1:size(m, 2), m', '-o'); hold on;
        yl = ylim; on = find(r.laser);
        plot(on, yl(2) * ones(size(on)), 'g^');
        title(sprintf('A %s: segment means (green = laser on)', mat2str(r.chans)));
        xlabel('trigger'); ylabel('relative');
    end
end
p = p + 1; subplot(nr, nc, p);
if isfield(out, 'A1') && isfield(out.A1, 'trace')
    tr = squeeze(out.A1.trace(1, :, :));
    if isvector(tr); tr = tr(:)'; end
    plot(tr'); title('A scan channel: decimated trace of each segment'); xlabel('sample block');
end
p = p + 1; subplot(nr, nc, p); hold on;
leg = {};
for f = {'C1', 'C2', 'C3', 'C4'}
    if isfield(out, f{1}) && isfield(out.(f{1}), 'v')
        plot(out.(f{1}).vec, out.(f{1}).v, '-o'); leg{end+1} = f{1}; %#ok<AGROW>
    end
end
if ~isempty(leg); legend(leg); end
title(sprintf('middle line along %s', out.geometry.letterA)); xlabel('um');
p = p + 1; subplot(nr, nc, p); hold on;
leg = {};
for f = {'D1', 'D2'}
    if isfield(out, f{1}) && isfield(out.(f{1}), 'v')
        plot(out.(f{1}).vec, out.(f{1}).v, '-o'); leg{end+1} = f{1}; %#ok<AGROW>
    end
end
if ~isempty(leg); legend(leg); end
title(sprintf('middle line along %s', out.geometry.letterB)); xlabel('um');
p = p + 1; subplot(nr, nc, p); hold on;
for f = {'C1', 'C2', 'D1', 'D2'}
    if isfield(out, f{1}) && isfield(out.(f{1}), 'posRead')
        plot(out.(f{1}).posRead - out.(f{1}).vec, '-o');
    end
end
title('read-back minus commanded position'); ylabel('um'); xlabel('pixel');
end

%% ===================================================================
%% helpers
%% ===================================================================
function g = geometry(sp)
from = sp.from; to = sp.to; np = sp.numPoints; fixedPos = sp.fixedPos; isFixed = logical(sp.isFixed);
scanIdx = find(~isFixed);
if numel(scanIdx) ~= 2
    error('scanDebug expects a 2D scan in the GUI (exactly two axes not fixed).');
end
% same choice of fast axis as StageScanner.optimize2dScanDirections (no chunking for the Standa)
if np(scanIdx(2)) <= np(scanIdx(1))
    g.idxA = scanIdx(1); g.idxB = scanIdx(2);
else
    g.idxA = scanIdx(2); g.idxB = scanIdx(1);
end
g.idxC = setdiff(1:3, scanIdx);
axesLetters = 'xyz';
g.letterA = axesLetters(g.idxA); g.letterB = axesLetters(g.idxB); g.letterC = axesLetters(g.idxC);
vec = @(i) linspace(from(i), to(i), np(i));
g.vecA = vec(g.idxA); g.vecB = vec(g.idxB);
g.pitchA = g.vecA(2) - g.vecA(1); g.pitchB = g.vecB(2) - g.vecB(1);
g.centerA = g.vecA(ceil(end / 2)); g.centerB = g.vecB(ceil(end / 2));
g.posC = fixedPos(g.idxC);
g.fixedXYZ = fixedPos;
g.tPixel = sp.pixelTime;
end

function info = collectInfo(ctx)
info.channelOption = ctx.spcm.currentChannelGUI;
info.scanChannels = ctx.spcm.activeChannelsGUI;
info.scanMode = ctx.spcm.scanAcquisitionMode;
info.triggerMargin = ctx.spcm.scanTriggerMargin;
info.stageBidirectional = ctx.stage.scanBidirectional;
info.stageLeadIn = ctx.stage.scanLeadIn;
info.stageSettle = ctx.stage.scanSettleTime;
info.digiTriggerTimeout = ctx.digi.acqInfo.TriggerTimeout;
info.pgDetectorOn = ctx.pg.isOn(ctx.o.triggerChannel);
end

function moveTo(ctx, letters, pos)
ctx.stage.Move(letters, pos);
end

function setLaser(ctx, on)
if on
    ctx.pg.On(ctx.o.laserChannel);
else
    ctx.pg.Off(ctx.o.laserChannel);
end
end

function s = gageStatus(digi)
s = CsMl_QueryStatus(digi.handle);
end

function log = pollStatus(digi, duration, t0)
% Status changes over 'duration' s: rows [time since t0, status].
tStart = tic;
last = -99; log = zeros(0, 2);
while toc(tStart) < duration
    s = gageStatus(digi);
    if s ~= last
        log(end+1, :) = [toc(t0), s]; %#ok<AGROW>
        last = s;
    end
end
log(end+1, :) = [toc(t0), gageStatus(digi)];
end

function s = waitNotBusy(digi, maxSeconds)
t = tic;
s = gageStatus(digi);
while s == 2 && toc(t) < maxSeconds
    s = gageStatus(digi);
end
end

function r = readSegments(ctx, r, chans)
digi = ctx.digi;
if r.statusEnd ~= 0
    digi.stopAcquisition;
    r.readError = sprintf('digitizer status %d after the last trigger (not complete), nothing read', r.statusEnd);
    return
end
digi.waitForDeviceReady;
raw = digi.readExperimentData(chans, true, true);      % [chan x segment x sample]
nS = size(raw, 3);
r.nSamples = nS;
r.segMean = reshape(mean(raw, 3, 'omitnan'), size(raw, 1), size(raw, 2));
r.segStd = reshape(std(raw, 0, 3), size(raw, 1), size(raw, 2));
r.nanFrac = reshape(mean(isnan(raw), 3), size(raw, 1), size(raw, 2));
nb = min(150, nS);
edges = round(linspace(0, nS, nb + 1));
tr = zeros(size(raw, 1), size(raw, 2), nb);
for b = 1:nb
    tr(:, :, b) = mean(raw(:, :, edges(b)+1:edges(b+1)), 3, 'omitnan');
end
r.trace = tr;
try
    [r.value, r.valueErr] = ctx.spcm.reduceSegments(raw);   % same reduction as the scan
catch
end
end

function c = digiConfig(digi)
c.acqInfo = digi.acqInfo; c.trigInfo = digi.trigInfo; c.sampleRateMHz = digi.sampleRate;
c.segmentPostTrigSize = digi.segmentPostTrigSize; c.segmentPreTrigSize = digi.segmentPreTrigSize;
end

function b = boardConfig(digi, chans)
try
    b = digi.readFromBoard(chans(1));
catch err
    b = struct('error', err.message);
end
end

function rate = pickRate(tUs, digi, maxSamples)
rates = digi.AVALIABLE_SAMPLE_RATES(digi.AVALIABLE_SAMPLE_RATES <= 1);
ok = rates(tUs .* rates <= maxSamples);
if isempty(ok); ok = rates(end); end
whole = abs(tUs .* ok - round(tUs .* ok)) < 1e-6 & round(tUs .* ok) >= 1;
if any(whole); rate = ok(find(whole, 1)); else; rate = ok(1); end
end

function waitSeconds(s)
t = tic;
if s > 0.005; pause(s - 0.003); end
while toc(t) < s
end
end

function z = onOffZ(m, laser)
a = m(laser); b = m(~laser);
if numel(a) < 2 || numel(b) < 2; z = NaN; return; end
se = sqrt(var(a) / numel(a) + var(b) / numel(b));
z = (mean(a) - mean(b)) / max(se, eps * max(abs([a b])) + realmin);
end

function s = noiseOf(y)
d = diff(y(:));
s = 1.4826 * median(abs(d - median(d))) / sqrt(2);
end

function r = runSafe(f, name)
fprintf('  running %s ...\n', name);
t = tic;
try
    r = f();
    r.seconds = toc(t);
catch err
    r = struct('error', err.message);
    if ~isempty(err.stack)
        r.where = sprintf('%s line %d', err.stack(1).name, err.stack(1).line);
    end
    fprintf('    ERROR in %s: %s\n', name, err.message);
end
end

function restoreState(ctx, state)
try; ctx.pg.Off(ctx.o.triggerChannel); catch; end
try; ctx.digi.stopAcquisition; catch; end
try; ctx.spcm.scanAcquisitionMode = state.mode; catch; end
try; ctx.stage.AbortScan(); catch; end
try; ctx.stage.Move(state.axes, state.pos); catch err; fprintf('  could not move back to the start position: %s\n', err.message); end
try; setLaser(ctx, state.laserOn); catch; end
try; ctx.spcm.setSPCMEnable(false); catch; end
end

function v = getOpt(s, f, d)
if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
