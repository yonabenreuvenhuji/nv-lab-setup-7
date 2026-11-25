function out = detectionNullTest(opts)
%DETECTIONNULLTEST Does the detection itself make the two ESR windows read differently?
%
%   out = detectionNullTest                       % all conditions, ESR timing
%   out = detectionNullTest(struct('window', 125, 'averages', 10))
%
% NO MICROWAVE, NO ESR. Runs ExpDetectionNullTest: the green laser is on the whole
% time and the detector opens and closes two windows per repeat with the CW-ESR timing
% (window, gap, window, gap). The light never changes, so both windows must read the
% same and every window must be flat. Anything else is the detection chain.
%
% CONDITIONS (all run by default)
%   A  detection only, ESR gap (10 us)
%   B  detection only, long gap (50 us)    -> does window 2 depend on how soon it
%                                             follows window 1 (digitizer re-arm)?
%   C  as A, plus the MW switches switched ON during window 1, as in the ESR, with the
%      AWG STOPPED so there is no microwave behind them
%                                          -> electrical pickup from the switch lines?
%
% WHAT IT REPORTS, per condition
%   offset = (window 1)/(window 2) - 1, from one ratio per repeat (laser drift cancels)
%   start  = (first 1 us of a window)/(rest of it) - 1, for each window
% and compares them with the ESR's own scales: the reference dip (~0.6% of the level)
% and the contrast (~5%).
%
% WHAT IT CAN AND CANNOT SHOW
%   A detection offset is the same at every MW frequency, so it can shift the whole
%   normalised ESR baseline but can NOT make a dip shaped like the ESR line. If all three
%   conditions are clean, the readout is ruled out at the level that matters.
%
% OPTIONS (struct fields, all optional)
%   window          us, detection window (default 125, the ESR window)
%   gaps            us, [ESR gap, long gap] (default [10 50])
%   switchChannels  channels toggled in condition C (default {'MW','MW2'})
%   conditions      which to run, any of 'ABC' (default 'ABC')
%   repeats         per pg.run (default 200)
%   averages        default 20
%   plot            default true

    if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
    opts = withDefaults(opts, struct('window', 125, 'gaps', [10 50], ...
        'switchChannels', {{'MW', 'MW2'}}, 'conditions', 'ABC', 'repeats', 200, ...
        'averages', 20, 'plot', true), 'detectionNullTest');

    conds = struct('id', {}, 'label', {}, 'gap', {}, 'toggle', {});
    if any(opts.conditions == 'A')
        conds(end+1) = struct('id', 'A', 'label', sprintf('detection only, gap %g us', opts.gaps(1)), ...
                              'gap', opts.gaps(1), 'toggle', {{}});
    end
    if any(opts.conditions == 'B')
        conds(end+1) = struct('id', 'B', 'label', sprintf('detection only, gap %g us', opts.gaps(2)), ...
                              'gap', opts.gaps(2), 'toggle', {{}});
    end
    if any(opts.conditions == 'C')
        conds(end+1) = struct('id', 'C', 'label', sprintf('MW switches on in window 1 (no MW), gap %g us', ...
                              opts.gaps(1)), 'gap', opts.gaps(1), 'toggle', {opts.switchChannels});
    end

    fprintf('\n=== detectionNullTest ===\n');
    fprintf('laser on throughout, %g us windows, %d repeats x %d averages per condition\n', ...
            opts.window, opts.repeats, opts.averages);
    stopAwg();          % no microwave can be produced, whatever the switches do

    results = cell(1, numel(conds));
    for k = 1:numel(conds)
        fprintf('\n--- condition %s: %s ---\n', conds(k).id, conds(k).label);
        e = ExpDetectionNullTest;
        e.detectionDuration = opts.window;
        e.referenceDetectionDuration = opts.window;
        e.gap = conds(k).gap;
        e.toggleChannels = conds(k).toggle;
        e.repeats = opts.repeats;
        e.averages = opts.averages;
        e.run;
        r = e.analyse(true);
        r.id = conds(k).id;
        r.label = conds(k).label;
        results{k} = r;
    end
    results = [results{:}];

    out = struct('results', results, 'opts', opts);
    out.verdict = verdict(results);
    report(results, out.verdict);
    if opts.plot; plotIt(results); end
end

% =========================================================================

function stopAwg()
    awg = [];
    try; awg = getObjByName('TektronixAWG'); catch; end
    if isempty(awg)
        try
            fgs = FrequencyGenerator.getFG();
            for i = 1:numel(fgs)
                if contains(fgs{i}.name, 'TektronixAWG'); awg = fgs{i}; break; end
            end
        catch
        end
    end
    if isempty(awg)
        fprintf('NOTE no AWG object found -- make sure nothing else is producing microwave.\n');
        return
    end
    awg.sendCommand('awgcontrol:stop:immediate');
    fprintf(['AWG stopped: no microwave is produced in any condition. (The next ESR re-arms it ' ...
             'by itself.)\n']);
end

% =========================================================================

function v = verdict(R)
    REF_DIP = 6e-3;     % the ESR reference dip, as a fraction of the level
    SMALL   = 0.1 * REF_DIP;
    lines = {};
    ok = true;
    transient = false;
    for k = 1:numel(R)
        r = R(k);
        sig = abs(r.offset) / max(r.ratioSterr, eps);
        if abs(r.offset) < SMALL || sig < 3
            lines{end+1} = sprintf('%s: windows agree (offset %+.1e +/- %.0e, under 10%% of the ESR reference dip).', ...
                                   r.id, r.offset, r.ratioSterr); %#ok<AGROW>
        else
            ok = false;
            lines{end+1} = sprintf(['%s: windows DIFFER by %+.2e +/- %.0e (%.0f%% of the ESR reference ' ...
                                    'dip). Same at every MW frequency, so it shifts the ESR baseline ' ...
                                    'but cannot make the dip.'], r.id, r.offset, r.ratioSterr, ...
                                   100 * abs(r.offset) / REF_DIP); %#ok<AGROW>
        end
        for w = 1:2
            s = r.(sprintf('start%d', w));
            if abs(s) > 1e-3
                transient = true;
                lines{end+1} = sprintf(['%s: window %d starts %+.1e off its own level in the first ' ...
                                        '1 us -- a gate/arming transient.'], r.id, w, s); %#ok<AGROW>
            end
        end
    end
    ids = [R.id];
    if all(ismember('AB', ids))
        a = R(ids == 'A'); b = R(ids == 'B');
        d = a.offset - b.offset; e = hypot(a.ratioSterr, b.ratioSterr);
        if abs(d) > 3 * e && abs(d) > SMALL
            ok = false;
            lines{end+1} = sprintf(['A vs B: window 2 depends on the gap (%+.1e +/- %.0e): the second ' ...
                                    'segment is affected by how soon it follows the first.'], d, e);
        end
    end
    if all(ismember('AC', ids))
        a = R(ids == 'A'); c = R(ids == 'C');
        d = c.offset - a.offset; e = hypot(a.ratioSterr, c.ratioSterr);
        if abs(d) > 3 * e && abs(d) > SMALL
            ok = false;
            lines{end+1} = sprintf(['C vs A: switching the MW switches changes window 1 by %+.1e +/- %.0e ' ...
                                    'with no microwave: electrical pickup into the detection.'], d, e);
        else
            lines{end+1} = sprintf('C vs A: no pickup from the MW switch lines (%+.1e +/- %.0e).', d, e);
        end
    end
    if ok && transient
        lines{end+1} = ['OVERALL: the two windows read the same, so the detection cannot be the ' ...
                        'source of the "without MW" dip. The start-of-window transient above is in ' ...
                        'both windows alike and cancels in the ratio.'];
    elseif ok
        lines{end+1} = ['OVERALL: the detection is clean. It cannot be the source of the ' ...
                        '"without MW" dip, nor of any baseline shift at the ESR''s level.'];
    else
        lines{end+1} = 'OVERALL: the detection is NOT perfectly symmetric -- see above for size and cause.';
    end
    v = strjoin(lines, newline);
end

% =========================================================================

function report(R, v)
    fprintf('\n=== detectionNullTest summary ===\n');
    fprintf('%-3s %-46s %12s %10s %11s %11s\n', '', 'condition', 'w1/w2 - 1', '+/-', 'start w1', 'start w2');
    for k = 1:numel(R)
        r = R(k);
        fprintf('%-3s %-46s %+12.2e %10.1e %+11.1e %+11.1e\n', r.id, r.label, r.offset, ...
                r.ratioSterr, r.start1, r.start2);
    end
    fprintf('(for scale: ESR reference dip ~6e-3, ESR contrast ~5e-2)\n\n%s\n\n', v);
end

% =========================================================================

function plotIt(R)
    figure('Name', 'detectionNullTest', 'Color', 'w');
    n = numel(R);
    for k = 1:n
        r = R(k);
        subplot(n + 1, 1, k);
        ref = mean(r.trace2);
        sm = max(1, round(numel(r.t) / 500));          % light smoothing for display only
        plot(r.t, movmean(r.trace1 / ref - 1, sm), 'b-'); hold on;
        plot(r.t, movmean(r.trace2 / ref - 1, sm), 'r-');
        yline(0, ':k');
        ylabel('relative'); grid on;
        title(sprintf('%s: %s   (w1/w2 - 1 = %+.1e +/- %.0e)', r.id, r.label, r.offset, r.ratioSterr));
        if k == 1; legend({'window 1', 'window 2'}, 'Location', 'best'); end
    end
    xlabel('time in window (us)');
    subplot(n + 1, 1, n + 1);
    errorbar(1:n, [R.offset], [R.ratioSterr], 'o', 'LineWidth', 1.5); hold on;
    yline(6e-4, '--k'); yline(-6e-4, '--k'); yline(0, ':k');
    set(gca, 'XTick', 1:n, 'XTickLabel', {R.id});
    xlim([0.5, n + 0.5]);
    ylabel('w1/w2 - 1'); grid on;
    title('window offsets (dashed: 10% of the ESR reference dip)');
end

% =========================================================================

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
