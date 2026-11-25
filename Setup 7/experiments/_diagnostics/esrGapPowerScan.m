function out = esrGapPowerScan(opts)
%ESRGAPPOWERSCAN Short, strong laser pulses between the ESR windows: which gap length is best?
%
%   out = esrGapPowerScan()                          % windows at -0.25, gaps at full power
%   out = esrGapPowerScan(struct('gaps', [5 10 20])) % your own gap lengths
%   out = esrGapPowerScan(struct('data', out))       % re-analyse a saved run, no hardware
%
% DIAGNOSTIC. Writes nothing. Uses AWGExpESR_gapAom, so it needs that class's cable:
% PulseStreamer output gapAomPgChannel (default 1) -> DAQ gapAomDaqTerminal (default PFI0).
%
% WHAT IT DOES
%   The ESR reference dip comes from NVs the MW knocks out and the laser resets only slowly
%   (~2.5 ms at the usual power). A short burst of full laser power between the windows may
%   reset them -- charge-state conversion speeds up roughly with power squared. This scan
%   keeps the detection windows at the usual power and puts full power in the gaps, for
%   several gap lengths, and measures for each:
%     - the signal and reference dips (as esrRefDipScan does: per average, projected on the
%       baseline line shape),
%     - the ESR contrast  = signal dip - reference dip,
%     - the SNR per unit time = contrast x sqrt(window / (window + gap)), relative to the
%       usual ESR. That is the number that decides: a longer gap only pays if the contrast it
%       recovers outweighs the time it adds.
%
%   Conditions (zero field, 2870 +/- 25 MHz, walk mode):
%     1. usual ESR: windows and gaps at the usual power, gap 10 us        (the reference)
%     2. full power in the gaps, for each gap in opts.gaps
%     3. control: the longest gap at the USUAL power -- separates "longer gap" from "more power"
%
%   The gap power is applied from the start of each gap until gapAomReturnLead (default 2 us)
%   before the next window, so the full-power time is about gap - 2 us (minus ~1 us DAQ/AOM
%   response). The shortest allowed gap is ~3 us.
%
% OPTIONS
%   windowAom   usual power, laser panel units (default -0.25)
%   gapAom      full power, panel units (default []: the top of the AOM range in the json,
%               i.e. 0 on Setup 7. If -0.6 is your BRIGHTEST setting, pass gapAom = -0.6.)
%   gaps        gap lengths at full power, us (default [4 6 10 20 50 100])
%   baseGap     gap of the usual ESR, us (default 10)
%   control     true (default): add the longest gap at usual power
%   window      detection window, us (default 100)
%   freqs, amplitude, repeats, averages, init, gapAomPgChannel, gapAomDaqTerminal,
%   gapAomReturnLead, helmholtzOff, jsonLocation   -- passed to esrRefDipScan
%   plot, data
%
% OUTPUT
%   out.scan   table rows: gap, gapAom, sig, ref, ratio, contrast, snrRel (+ errors), timeX
%   out.best   the recommended row (or empty if nothing beats the usual ESR)
%   out.raw    the esrRefDipScan result (re-analysable on its own)

    if ~exist('opts', 'var'); opts = struct(); end
    d = struct('data', [], 'windowAom', -0.25, 'gapAom', [], 'gaps', [4 6 10 20 50 100], ...
               'baseGap', 10, 'control', true, 'window', 100, 'plot', true);
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(opts, f{i}); opts.(f{i}) = d.(f{i}); end
    end

    if ~isempty(opts.data)
        raw = opts.data.raw;
        raw = esrRefDipScan(struct('data', raw, 'plot', false));
    else
        if isempty(opts.gapAom)
            part = AWGExpESR_gapAom.gapAomPart();
            opts.gapAom = part.maxValue;
            fprintf(['Full power taken as %g (top of the AOM range). If your brightest setting is %g, ' ...
                     'stop and rerun with gapAom = %g.\n'], opts.gapAom, part.minValue, part.minValue);
        end
        W = opts.window;
        C = {struct('label', sprintf('usual ESR: gap %g us, power %g', opts.baseGap, opts.windowAom), ...
                    'role', 'base', 'gap', opts.baseGap, 'window', W)};
        for g = opts.gaps(:)'
            C{end+1} = struct('label', sprintf('gap %g us at full power (%g)', g, opts.gapAom), ...
                              'role', 'gapAom', 'gap', g, 'window', W, 'gapAom', opts.gapAom); %#ok<AGROW>
        end
        if opts.control && max(opts.gaps) ~= opts.baseGap
            C{end+1} = struct('label', sprintf('control: gap %g us at usual power', max(opts.gaps)), ...
                              'role', 'gapControl', 'gap', max(opts.gaps), 'window', W);
        end
        pass = {'freqs', 'amplitude', 'repeats', 'averages', 'init', 'gapAomPgChannel', ...
                'gapAomDaqTerminal', 'gapAomReturnLead', 'helmholtzOff', 'jsonLocation'};
        ro = struct('conditions', {C}, 'windowAom', opts.windowAom, 'plot', false);
        for i = 1:numel(pass)
            if isfield(opts, pass{i}); ro.(pass{i}) = opts.(pass{i}); end
        end
        raw = esrRefDipScan(ro);
    end

    out = struct('raw', raw, 'opts', rmfieldSafe(opts, 'data'));
    [out.scan, out.best, out.verdict] = evaluate(raw);
    report(out);
    if isempty(opts.data)
        try
            fn = fullfile(setup7Paths().folder, 'esrGapPowerScan.mat');
            save(fn, 'out');
            fprintf('saved to %s   (re-analyse: esrGapPowerScan(struct(''data'', out)))\n\n', fn);
        catch
        end
    end
    if opts.plot; plotScan(out); end
end

% ======================================================================

function [T, best, v] = evaluate(raw)
    Cd = raw.cond;
    ib = find(strcmp({Cd.role}, 'base'), 1);
    if isempty(ib); error('esrGapPowerScan: no usual-ESR condition in the data.'); end
    b = Cd(ib);
    snr = @(c) (c.aSig - c.aRef) * sqrt(c.window / (c.window + c.gap));
    snrErr = @(c) sqrt(c.aSigErr^2 + c.aRefErr^2) * sqrt(c.window / (c.window + c.gap));
    s0 = snr(b); e0 = snrErr(b);
    T = struct([]);
    for i = 1:numel(Cd)
        c = Cd(i);
        ga = []; if isfield(c, 'gapAom'); ga = c.gapAom; end
        r = struct('label', c.label, 'role', c.role, 'gap', c.gap, 'gapAom', ga, ...
                   'sig', c.aSig, 'sigErr', c.aSigErr, 'ref', c.aRef, 'refErr', c.aRefErr, ...
                   'ratio', c.ratio, 'ratioErr', c.ratioErr, ...
                   'contrast', c.aSig - c.aRef, 'contrastErr', sqrt(c.aSigErr^2 + c.aRefErr^2));
        r.snrRel = snr(c) / s0;
        r.snrRelErr = abs(r.snrRel) * sqrt((snrErr(c)/max(abs(snr(c)), eps))^2 + (e0/max(abs(s0), eps))^2);
        r.timeX = 1 / max(r.snrRel, eps)^2;    % averaging time for the same SNR, vs usual ESR
        if isempty(T); T = r; else; T(end+1) = r; end %#ok<AGROW>
    end

    v = {};
    ih = find(strcmp({T.role}, 'gapAom'));
    best = [];
    if isempty(ih)
        v{end+1} = 'No full-power conditions in the data.';
        return
    end
    % best = highest SNR per time; among those within 1 sigma of it, the SHORTEST gap
    [~, k] = max([T(ih).snrRel]);
    top = T(ih(k));
    near = ih([T(ih).snrRel] >= top.snrRel - top.snrRelErr);
    [~, kk] = min([T(near).gap]);
    cand = T(near(kk));

    v{end+1} = sprintf(['Usual ESR: contrast %.2f%% (signal dip %.2f%%, reference dip %.2f%%, ratio %.3f).'], ...
                       100*T(ib).contrast, 100*T(ib).sig, 100*T(ib).ref, T(ib).ratio);
    gain = cand.snrRel - 1; gErr = cand.snrRelErr;
    if gain > max(2*gErr, 0.03)
        best = cand;
        v{end+1} = sprintf(['BEST: gap %g us at full power. ESR contrast %.2f%% (ratio %.3f), SNR per time ' ...
                            '%.2f +/- %.2f x the usual ESR, i.e. %.2fx the averaging time for the same SNR.'], ...
                           cand.gap, 100*cand.contrast, cand.ratio, cand.snrRel, cand.snrRelErr, cand.timeX);
        v{end+1} = sprintf(['Use it: a = AWGExpESR_gapAom; a.delayBetweenDetections = %g; a.gapAomValue = %g; ' ...
                            '(laser panel at the usual power for the windows).'], cand.gap, cand.gapAom);
    elseif cand.snrRel > 1 - max(2*gErr, 0.03)
        v{end+1} = sprintf(['No clear gain: the best full-power gap (%g us) gives %.2f +/- %.2f x the usual SNR per ' ...
                            'time. Not worth the extra setup; keep the usual ESR.'], cand.gap, cand.snrRel, cand.snrRelErr);
    else
        v{end+1} = sprintf(['Full power in the gaps does not pay: the best (%g us) gives only %.2f +/- %.2f x the ' ...
                            'usual SNR per time. Keep the usual ESR.'], cand.gap, cand.snrRel, cand.snrRelErr);
    end
    % does the strong pulse remove the reference dip, and is it the power or the length?
    okr = ih(isfinite([T(ih).ratio]));
    if ~isempty(okr)
        [~, km] = min([T(okr).ratio]);
        m = T(okr(km));
        v{end+1} = sprintf('Smallest reference dip: ratio %.3f +/- %.3f at gap %g us full power (usual ESR %.3f).', ...
                           m.ratio, m.ratioErr, m.gap, T(ib).ratio);
    end
    ic = find(strcmp({T.role}, 'gapControl'), 1);
    if ~isempty(ic)
        ihc = ih([T(ih).gap] == T(ic).gap);
        if ~isempty(ihc) && isfinite(T(ic).ratio) && isfinite(T(ihc).ratio)
            h = T(ihc); c = T(ic);
            dr = c.ratio - h.ratio; de = hypot(c.ratioErr, h.ratioErr);
            if dr > 2*de
                v{end+1} = sprintf(['At %g us, full power gives ratio %.3f against %.3f at the usual power: it is the ' ...
                                    'POWER that clears the reference dip, not the length.'], h.gap, h.ratio, c.ratio);
            else
                v{end+1} = sprintf(['At %g us, full power (%.3f) and the usual power (%.3f) give the same ratio: the ' ...
                                    'extra power does not help -- only gap length does.'], h.gap, h.ratio, c.ratio);
            end
        end
    end
end

function report(out)
    T = out.scan;
    fprintf('\n=== esrGapPowerScan ===\n\n');
    fprintf('%-38s %8s %8s %7s %9s %13s %8s\n', 'condition', 'signal', 'ref', 'ratio', 'contrast', ...
            'SNR/time (x)', 'time (x)');
    for i = 1:numel(T)
        r = T(i);
        if isfinite(r.ratio); rs = sprintf('%.3f', r.ratio); else; rs = '--'; end
        fprintf('%-38s %7.3f%% %7.3f%% %7s %8.3f%% %6.2f +/- %4.2f %7.2f\n', r.label, 100*r.sig, 100*r.ref, ...
                rs, 100*r.contrast, r.snrRel, r.snrRelErr, r.timeX);
    end
    fprintf('\nWHAT IT SAYS\n');
    for i = 1:numel(out.verdict)
        fprintf('   %s\n', out.verdict{i});
    end
    fprintf('\n');
end

function plotScan(out)
    T = out.scan;
    ib = find(strcmp({T.role}, 'base'), 1);
    ih = find(strcmp({T.role}, 'gapAom'));
    ic = find(strcmp({T.role}, 'gapControl'), 1);
    [~, o] = sort([T(ih).gap]); ih = ih(o);
    g = [T(ih).gap];
    blue = [0.15 0.42 0.69]; grey = [0.45 0.45 0.45]; red = [0.80 0.25 0.15];
    figure('Color', 'w', 'Position', [100 100 1000 420]);
    for p = 1:2
        ax = subplot(1, 2, p); hold(ax, 'on');
        set(ax, 'XScale', 'log', 'Box', 'off', 'TickDir', 'out', 'XGrid', 'on', 'YGrid', 'on', ...
            'GridAlpha', 0.12, 'FontSize', 9);
        if p == 1
            y = [T(ih).ratio]; e = [T(ih).ratioErr]; y0 = T(ib).ratio;
            yl = 'reference dip / signal dip'; tl = 'reference dip left';
        else
            y = [T(ih).snrRel]; e = [T(ih).snrRelErr]; y0 = 1;
            yl = 'SNR per time  (usual ESR = 1)'; tl = 'what it buys';
        end
        plot(ax, [min(g)/1.3, max(g)*1.3], [y0 y0], '--', 'Color', grey);
        h1 = errorbar(ax, g, y, e, 'o-'); styleH(h1, blue, 1.2);
        hs = h1; ls = {'full power in the gaps'};
        if ~isempty(ic)
            c = T(ic);
            if p == 1; yc = c.ratio; ec = c.ratioErr; else; yc = c.snrRel; ec = c.snrRelErr; end
            hs(end+1) = errorbar(ax, c.gap, yc, ec, 's'); styleH(hs(end), red, 1);
            ls{end+1} = 'usual power (control)';
        end
        if p == 2 && ~isempty(out.best)
            plot(ax, out.best.gap, out.best.snrRel, 'p', 'MarkerSize', 14, 'Color', [0.1 0.1 0.1], ...
                 'MarkerFaceColor', [1 0.85 0.2]);
        end
        xlim(ax, [min(g)/1.3, max(g)*1.3]);
        xlabel(ax, 'gap (delayBetweenDetections), us'); ylabel(ax, yl);
        title(ax, tl, 'FontWeight', 'bold');
        if p == 1; legend(ax, hs, ls, 'Location', 'best', 'Box', 'off'); end
    end
end

function styleH(h, col, lw)
    try
        set(h, 'Color', col, 'MarkerFaceColor', col, 'LineWidth', lw);
    catch
        try; set(h, 'Color', col); catch; end
    end
end

function s = rmfieldSafe(s, f)
    if isfield(s, f); s = rmfield(s, f); end
end
