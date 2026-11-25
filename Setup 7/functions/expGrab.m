function g = expGrab(e, tag)
%EXPGRAB Copy everything worth keeping out of an Experiment into a plain struct.
%
%   g = expGrab(r)
%   g = expGrab(r, 'walk')
%
% ===================================================================
% WHY THIS IS NECESSARY, AND IT IS NOT OBVIOUS
% ===================================================================
%
% Experiment objects are SINGLETONS KEYED BY CLASS NAME. Experiment's constructor does
%
%       oldObjOrNan = removeObjIfExists(name);
%       if isa(oldObjOrNan, 'Experiment')
%           oldObjOrNan.delete();
%           fprintf('Creating a new %s Experiment. All of its parameters were reset...')
%       end
%
% so constructing a second AWGExpRamsey DELETES the first one. Not shadows -- deletes. The
% variable holding it is still there, but it is a dead handle, and the next thing you do with
% it is
%
%       Invalid or deleted object.
%
% That is why comparing two runs of the same experiment class cannot be done by keeping two
% objects around:
%
%       r  = AWGExpRamsey;  r.run;          % walk
%       rc = AWGExpRamsey;                  % <-- r is dead from here
%       rc.tau = r.tau;                     % Invalid or deleted object
%
% The data has to come OUT of the object before the next one is built. That is all this does.
% The returned struct is plain data and survives anything that happens to the object.
%
%       W = expGrab(r, 'walk');       % after the first run
%       rc = AWGExpRamsey;  ...       % kills r, W is untouched
%       C = expGrab(rc, 'classic');
%
% The same trap applies to a quick-look run: `rq = AWGExpRamsey` in a scratch cell silently
% destroys the r you were about to plot.
%
% WHAT IT COPIES
%   the x axis, the two ratio traces with their errors, the raw per-window signal and sterr,
%   and the settings needed to know whether two runs are comparable at all. Everything is
%   copied, not referenced, and missing fields come back empty rather than raising -- a grab
%   that throws on a half-finished run is useless exactly when you want it.

    if ~exist('tag', 'var'); tag = ''; end
    g = struct('tag', tag, 'name', '', 'when', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
               'x', [], 'value', [], 'sterr', [], 'value2', [], 'sterr2', [], ...
               'signal', [], 'rawSterr', [], 'currIter', NaN, 'settings', struct());  %#ok<TNOW1,DATST>

    g.name = tryGet(e, 'NAME', '');

    % ---- the x axis. tau is the property; mCurrentXAxisParam.value is what the plot used.
    g.x = tryGet(e, 'tau', []);
    ax = tryGet(e, 'mCurrentXAxisParam', []);
    if ~isempty(ax)
        v = tryGet(ax, 'value', []);
        if ~isempty(v); g.x = v; end
    end
    g.x = g.x(:)';

    % ---- the reduced traces
    [g.value,  g.sterr]  = paramData(e, 'signalParam');
    [g.value2, g.sterr2] = paramData(e, 'signalParam2');

    % ---- the raw per-window arrays, trimmed to the averages actually done
    g.currIter = tryGet(e, 'currIter', NaN);
    sig = tryGet(e, 'signal', []);
    ste = tryGet(e, 'sterr', []);
    n = g.currIter;
    if ~isempty(sig) && isfinite(n) && n >= 1 && ndims(sig) == 3 && size(sig, 3) >= n
        g.signal   = sig(:, :, 1:n);
        if ~isempty(ste) && size(ste, 3) >= n; g.rawSterr = ste(:, :, 1:n); end
    else
        g.signal = sig; g.rawSterr = ste;
    end

    % ---- the settings that decide whether two runs are comparable
    keep = {'frequency', 'amplitude', 'amplitudeDiff', 'channelPhase', 'timeDelay', ...
            'nChannels', 'halfPiTime', 'piTime', 'detectionDuration', ...
            'referenceDetectionDuration', 'laserInitializationDuration', 'lastDelay', ...
            'repeats', 'averages', 'constantTime', 'balancedSequence', 'doubleMeasurement', ...
            'smallDelay', 'isTracking', 'useSequenceWalk', 'passesPerRun', ...
            'passesPerRunUsed', 'walkWarmupTriggers', 'resyncEveryPass', ...
            'maxDigitizerSegments', 'cycles', 'pulsePhase', ...
            'mwGateMode', 'mwHoldUs', 'mwGateOpenUs', 'mwHoldUsed'};
    % mwGateMode and friends are AWGExpT1's. They are in this list because the DOSE an arm was
    % taken at is not recoverable from the traces -- two t1GateTest arms differ in nothing else,
    % so a grab that drops the gate mode is a grab of an unlabelled measurement. isprop guards
    % every name here, so classes without them are unaffected.
    for i = 1:numel(keep)
        if isprop(e, keep{i})
            g.settings.(keep{i}) = tryGet(e, keep{i}, []);
        end
    end
end

% ======================================================================

function [v, s] = paramData(e, which)
    v = []; s = [];
    p = tryGet(e, which, []);
    if isempty(p); return; end
    v = tryGet(p, 'value', []);
    s = tryGet(p, 'sterr', []);
    v = v(:)'; s = s(:)';
end

function v = tryGet(o, name, dflt)
% Deliberately forgiving. A grab is most useful on a run that stopped early or was
% interrupted, which is exactly when half the fields are missing -- so nothing here raises.
    v = dflt;
    try
        if isobject(o) || isstruct(o)
            v = o.(name);
        end
    catch
    end
    if isempty(v) && ~isempty(dflt); v = dflt; end
end
