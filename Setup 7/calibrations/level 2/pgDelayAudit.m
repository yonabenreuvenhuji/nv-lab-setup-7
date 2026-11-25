function out = pgDelayAudit(gap, opts)
%PGDELAYAUDIT Why fixSequenceDelays produces a negative duration, without running anything.
%
%   pgDelayAudit                 % audit the whole calibLastDelay gap scan
%   pgDelayAudit(0.005)          % just that gap, verbosely
%   out = pgDelayAudit(0.005, struct('mwBurst', 1, 'detectionDuration', 2.76))
%
% NO HARDWARE. Reads the channel delays out of setupInfo.json and replays the arithmetic in
% PulseGenerator.fixSequenceDelays on paper. Safe to run any time.
%
% ===================================================================
% WHAT THIS IS FOR
% ===================================================================
%
% calibLastDelay dies at small gaps with
%
%   Error using Pulse/set.duration
%   Duration must not be negative! Requested: -2.174000e-01
%   ... in PulseGenerator/fixSequenceDelays (line 167)
%
% and then retries ten times and emergency-stops, which is a lot of output saying nothing.
% The failure is NOT in the sequence as written -- every event in it is positive. It is in
% the DELAY COMPENSATION applied on top, and it depends on the gap, which is why the scan
% dies at 0.005 us while ordinary Rabi runs happily at lastDelay = 0.55.
%
% ===================================================================
% THE ARITHMETIC, WHICH IS THE WHOLE ANSWER
% ===================================================================
%
% fixSequenceDelays takes each channel's pulse edges and moves them so the LIGHT (or the
% microwave, or the trigger) arrives when the sequence said it should:
%
%       start := start - onDelay
%       end   := end   - offDelay
%
% and then folds anything that went negative around to the end of the sequence. Two things
% fall out of that, and both can bite:
%
%   1. A pulse's duration becomes  D + onDelay - offDelay.  If a channel's offDelay exceeds
%      its onDelay by more than the pulse is long, the corrected duration is NEGATIVE. Short
%      pulses on a channel with an offDelay bigger than its onDelay are the danger, and the
%      scanned gap is exactly such a pulse -- it starts at 5 ns.
%
%   2. A pulse whose start goes negative is split at zero and the front half wrapped to the
%      end of the sequence. That is fine on its own, but the split is done with
%          insert = @(a,x,n) cat(2, x(1:n), a, x(n+1:end))
%      and if TWO pulses on one channel straddle zero, indSplit is a 2-vector and x(1:n)
%      silently uses only its first element. The second straddling pulse is then never
%      split, its start is wrapped while its end is not, and the duration comes out as
%      end - start - sequenceDuration: a large negative number. On this setup the AWG
%      trigger has onDelay 1.3676 us, far larger than any other channel, so it is the one
%      that gets wrapped.
%
% This function prints, per channel, the raw pulses, the corrected pulses and every
% resulting duration, and marks the negative one. That turns "-0.2174" into a channel name
% and a number you can act on.
%
% OPTIONS (struct fields)
%   gaps                the gap list to audit (default 0.005:0.025:0.78, the level-2 scan)
%   mwBurst             MW pulse length, us (default 1)
%   detectionDuration   us (default 2.76)
%   referenceDetection  us (default = detectionDuration)
%   laserInit           us (default 70.5536)
%   trigDuration        us (default 0.1, FrequencyGeneratorTektronixAWG.TRIG_DURATION)
%   nChannels           1 or 2 (default 2)
%   verbose             per-pulse table (default true for one gap, false for a scan)
%
% OUTPUT
%   out.delays      the channel delays it read
%   out.minSafeGap  smallest gap in the list with no negative duration, NaN if none
%   out.perGap      struct array, one per gap: .gap .ok .worst .channel

    if ~exist('gap', 'var'); gap = []; end
    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('gaps', 0.005:0.025:0.78, 'mwBurst', 1, ...
                        'detectionDuration', 2.76, 'referenceDetection', [], ...
                        'laserInit', 70.5536, 'trigDuration', 0.1, 'nChannels', 2, ...
                        'verbose', [], 'expect', []), 'pgDelayAudit');
    if ~isempty(gap); opts.gaps = gap; end
    if isempty(opts.referenceDetection); opts.referenceDetection = opts.detectionDuration; end
    if isempty(opts.verbose); opts.verbose = numel(opts.gaps) == 1; end

    D = readDelays();
    out = struct('delays', D, 'minSafeGap', NaN, 'perGap', []);

    fprintf('\n=== pgDelayAudit ===\n');
    fprintf('channel delays from setupInfo.json (us):\n');
    fprintf('   %-12s %-10s %-10s %s\n', 'channel', 'onDelay', 'offDelay', 'off - on');
    names = fieldnames(D);
    worstOffMinusOn = -Inf; worstName = '';
    for i = 1:numel(names)
        d = D.(names{i});
        fprintf('   %-12s %-10.4f %-10.4f %+.4f\n', names{i}, d(1), d(2), d(2)-d(1));
        if d(2)-d(1) > worstOffMinusOn; worstOffMinusOn = d(2)-d(1); worstName = names{i}; end
    end
    fprintf(['\n   The rightmost column is what matters for mechanism 1: a pulse shorter than\n' ...
             '   (off - on) on that channel becomes negative. Worst here is %s at %+.4f us.\n'], ...
            worstName, worstOffMinusOn);
    maxOn = 0; maxOnName = '';
    for i = 1:numel(names)
        if D.(names{i})(1) > maxOn; maxOn = D.(names{i})(1); maxOnName = names{i}; end
    end
    fprintf(['   And the largest onDelay is %s at %.4f us, which is what decides mechanism 2:\n' ...
             '   any pulse starting less than %.4f us into the sequence gets wrapped to the end.\n'], ...
            maxOnName, maxOn, maxOn);

    % ---- per gap ----
    per = struct('gap', {}, 'ok', {}, 'worst', {}, 'channel', {});
    for g = opts.gaps(:)'
        % BOTH merge models. Whether perChannelEdgeTimes fuses back-to-back events on
        % one channel into a single pulse decides whether the scanned gap survives as a
        % pulse of its own -- and a pulse only 5 ns long is exactly what mechanism 1
        % needs. Sequence.m is not something this file can read, so it tries both and
        % reports whichever produces a negative.
        [ok, worst, ch, detail] = auditOne(g, D, opts, true);
        if ok
            [ok2, worst2, ch2, detail2] = auditOne(g, D, opts, false);
            if ~ok2
                ok = ok2; worst = worst2; ch = ch2; detail = detail2;
                if opts.verbose
                    fprintf(['\n   (the MERGED model is clean; what follows is the UNMERGED one,\n' ...
                             '    i.e. each addEvent kept as its own pulse)\n']);
                end
            end
        end
        per(end+1) = struct('gap', g, 'ok', ok, 'worst', worst, 'channel', ch); %#ok<AGROW>
        if opts.verbose
            fprintf('\n--- gap %.4f us ---\n', g);
            fprintf('   %-12s %-9s %-20s %-20s %s\n', 'channel', 'pulse', 'raw [start end]', ...
                    'corrected [start end]', 'duration');
            for k = 1:numel(detail)
                dd = detail(k);
                mark = ''; if dd.dur < 0; mark = '   <== NEGATIVE'; end
                fprintf('   %-12s %-9d [%7.4f %8.4f]  [%8.4f %9.4f]  %+9.4f%s\n', ...
                        dd.channel, dd.index, dd.rawStart, dd.rawEnd, ...
                        dd.start, dd.end, dd.dur, mark);
            end
        end
    end
    out.perGap = per;

    okGaps = [per([per.ok]).gap];
    if isempty(okGaps)
        fprintf('\n   EVERY gap in the list produces a negative duration.\n');
    else
        out.minSafeGap = min(okGaps);
        bad = [per(~[per.ok]).gap];
        if isempty(bad)
            fprintf('\n   Every gap in the list is safe. The crash is not explained by this model --\n');
            fprintf('   see the note at the end.\n');
        else
            fprintf('\n   %d of %d gaps fail. Smallest SAFE gap: %.4f us.\n', ...
                    numel(bad), numel(per), out.minSafeGap);
            fprintf('   Failing gaps: %s\n', numList(bad));
            fprintf(['\n   So the scan can be run as-is from %.4f us upward:\n' ...
                     '       calibLastDelay(struct(''gaps'', %.4f:0.025:%.3f, ...))\n'], ...
                    out.minSafeGap, out.minSafeGap, max(opts.gaps));
        end
    end

    % ---- back-solve against the number MATLAB actually reported ----
    if ~isempty(opts.expect) && isfinite(opts.expect)
        fprintf('\n--- back-solving for the reported %.4f us ---\n', opts.expect);
        fprintf(['   For each channel and each pulse in the sequence, the offDelay - onDelay\n' ...
                 '   that would produce exactly that duration. Compare against the measured\n' ...
                 '   column above: whichever row matches is the delay that is wrong.\n\n']);
        fprintf('   %-12s %-10s %-14s %-14s %s\n', 'channel', 'pulse len', 'needs off-on', ...
                'actual off-on', 'match?');
        for g = opts.gaps(:)'
            for mg = [true false]
                [~, ~, ~, det] = auditOne(g, D, opts, mg);
                for k = 1:numel(det)
                    len = det(k).rawEnd - det(k).rawStart;
                    need = len - opts.expect;
                    act  = D.(det(k).channel)(2) - D.(det(k).channel)(1);
                    if abs(need - act) < 5e-4
                        fprintf('   %-12s %-10.4f %-14.4f %-14.4f  <== YES (gap %.4f)\n', ...
                                det(k).channel, len, need, act, g);
                    end
                end
            end
            if g ~= opts.gaps(1); break; end   % one gap's worth of rows is enough
        end
        fprintf(['\n   No YES row means the sequence modelled here is not the one that ran.\n' ...
                 '   Send this whole printout and the failing command.\n']);
    end

    fprintf(['\n   IF THE NEGATIVE DURATION THIS PREDICTS DOES NOT MATCH THE ONE MATLAB\n' ...
             '   REPORTED, the delays above are not the ones in force -- setupInfo.json has\n' ...
             '   been edited since main last ran, or the sequence differs from the one modelled\n' ...
             '   here. Restart main and re-run this, and compare the two numbers directly.\n']);
    fprintf('\n');
end

% ======================================================================

function [ok, worst, worstCh, detail] = auditOne(gap, D, o, merge)
% Build AWGExpLastDelayDuration's sequence for this gap, apply fixSequenceDelays, and
% return every corrected pulse. Event list transcribed from that class's prepare().
    initDur = o.laserInit - o.detectionDuration;
    if o.nChannels == 2; mw = {'MW','MW2'}; else; mw = {'MW'}; end

    ev = {};   % {duration, {channels}}
    ev = addEv(ev, o.trigDuration, [{'trigger'}, mw]);
    ev = addEv(ev, o.mwBurst,      mw);
    ev = addEv(ev, gap,            mw);                       % <-- the scanned gap
    ev = addEv(ev, o.detectionDuration,  {'greenLaser','detector'});
    ev = addEv(ev, initDur,        {'greenLaser'});
    ev = addEv(ev, o.trigDuration, {});
    ev = addEv(ev, max(o.gaps),    {});                       % tau(end)
    ev = addEv(ev, gap,            {});
    ev = addEv(ev, o.referenceDetection, {'greenLaser','detector'});
    ev = addEv(ev, initDur,        {'greenLaser'});

    % absolute edges and total duration
    t = 0; edges = zeros(numel(ev), 2);
    for i = 1:numel(ev); edges(i,:) = [t, t + ev{i}{1}]; t = t + ev{i}{1}; end
    Dseq = t;

    chNames = fieldnames(D);
    detail = struct('channel', {}, 'index', {}, 'rawStart', {}, 'rawEnd', {}, ...
                    'start', {}, 'end', {}, 'dur', {});
    worst = Inf; worstCh = '';
    for c = 1:numel(chNames)
        name = chNames{c};
        % contiguous events on this channel merge into one pulse, as perChannelEdgeTimes does
        raw = zeros(0,2);
        for i = 1:numel(ev)
            if any(strcmp(name, ev{i}{2}))
                if merge && ~isempty(raw) && abs(raw(end,2) - edges(i,1)) < 1e-12
                    raw(end,2) = edges(i,2);
                else
                    raw(end+1,:) = edges(i,:); %#ok<AGROW>
                end
            end
        end
        if isempty(raw); continue; end
        rawKeep = raw;

        on = D.(name)(1); off = D.(name)(2);
        st = raw(:,1)' - on;
        en = raw(:,2)' - off;

        % the split, reproduced INCLUDING the single-element behaviour of x(1:n)
        idx = find(st < 0 & en > 0);
        if ~isempty(idx)
            i1 = idx(1);                      % x(1:n) uses only n(1)
            st = [st(1:i1), 0, st(i1+1:end)];
            en = [en(1:i1-1), 0, en(i1:end)];
            rawKeep = [rawKeep(1:i1,:); rawKeep(i1,:); rawKeep(i1+1:end,:)];
        end
        st = st + (st < 0)  * Dseq;
        en = en + (en <= 0) * Dseq;

        for j = 1:numel(st)
            dur = en(j) - st(j);
            detail(end+1) = struct('channel', name, 'index', j, ...
                                   'rawStart', rawKeep(min(j,size(rawKeep,1)),1), ...
                                   'rawEnd',   rawKeep(min(j,size(rawKeep,1)),2), ...
                                   'start', st(j), 'end', en(j), 'dur', dur); %#ok<AGROW>
            if dur < worst; worst = dur; worstCh = name; end
        end
    end
    ok = worst >= 0;
end

function ev = addEv(ev, dur, channels)
    ev{end+1} = {dur, channels};
end

function D = readDelays()
% Every pulse-generator channel's on/off delay, from the same places calibChannelDelays
% looks: lasers keep theirs under 'switch', microwave switches under 'microwave', the AWG
% trigger under the frequency generator as onDelay1/offDelay1.
    D = struct();
    try
        j = JsonInfoReader.getJson();
    catch err
        error('pgDelayAudit: could not read the setup json (%s). Run main first.', err.message);
    end
    try
        for L = asList(j.lasers)
            LL = L{1};
            if isstruct(LL) && isfield(LL, 'switch')
                D.(vname(LL.switch.switchChannelName)) = ...
                    [LL.switch.onDelay, LL.switch.offDelay];
            end
        end
    catch
    end
    try
        f = fieldnames(j.microwave);
        for i = 1:numel(f)
            m = j.microwave.(f{i});
            D.(vname(m.switchChannelName)) = [m.onDelay, m.offDelay];
        end
    catch
    end
    try
        for G = asList(j.frequencyGenerators)
            g = G{1};
            if isfield(g, 'trigger1ChannelName')
                offd = 0;
                if isfield(g, 'offDelay1'); offd = g.offDelay1; end
                D.(vname(g.trigger1ChannelName)) = [g.onDelay1, offd];
            end
        end
    catch
    end
    % Channels with no json entry get 0/0, which is what channelName2Delays returns.
    for extra = {'detector'}
        if ~isfield(D, extra{1}); D.(extra{1}) = [0 0]; end
    end
    if isempty(fieldnames(D))
        error('pgDelayAudit: found no channel delays in the setup json.');
    end
end

function n = vname(s)
    n = regexprep(char(s), '\W', '_');
    if isempty(regexp(n, '^[A-Za-z]', 'once')); n = ['x' n]; end
end

function c = asList(x)
    if iscell(x); c = x; elseif isstruct(x); c = num2cell(x); else; c = {x}; end
end

function s = numList(v)
% Comma-separated numbers. Not compose(): that is R2016b+ and absent in Octave, and this
% file is checked offline.
    p = cell(1, numel(v));
    for i = 1:numel(v); p{i} = sprintf('%.4f', v(i)); end
    s = strjoin(p, ', ');
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
        if ~isfield(o, known{i}); o.(known{i}) = d.(known{i}); end
    end
end
