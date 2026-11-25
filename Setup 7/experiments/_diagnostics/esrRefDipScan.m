function out = esrRefDipScan(opts)
%ESRREFDIPSCAN Why does the "without MW" (reference) ESR trace dip at resonance? Scan and decide.
%
%   out = esrRefDipScan()                              % the default set, zero field
%   out = esrRefDipScan(struct('conditions', {{...}})) % your own set
%   out = esrRefDipScan(struct('data', out))           % re-analyse a saved run, no hardware
%
% DIAGNOSTIC. Writes nothing. Runs AWGExpESR at zero field over a narrow sweep around 2870 MHz
% under several conditions, measures the dip in BOTH traces each time, and prints which
% explanation the pattern supports and what to do.
%
% ===================================================================
% WHAT IS KNOWN (2026-11-04)
% ===================================================================
%   - The reference trace is the real averaged reference level (processData does not mix it
%     with the signal), and it dips at the NV line by ~10% of the signal dip.
%   - Detector memory: ruled out (detectorMemoryTest).
%   - Spin memory with the level-2 time constants: ruled out. Gap 100 us + window 200 us should
%     have shrunk the dip ~20x and did not change it.
%   - The fluorescence-after-dark effect: cannot make a resonance-shaped dip (frequency
%     independent), and both ESR windows share the same history anyway.
%
% ===================================================================
% WHAT EACH CONDITION TESTS
% ===================================================================
%   baseline          gap 10 us, window 100 us             the reference value
%   switches closed   MW switches never opened             no MW at the NVs: any dip left is
%                                                          leakage (both traces) or processing
%                                                          (reference only)
%   gap 300, 2000 us  laser on through the gap             memory up to ~ms: shrinks if the laser
%                                                          erases it; two gaps give its time constant
%   long window       window 200 us                        dilution of a memory at window start
%   MW -10 dB         lower MW power                       ratio GROWS at low power if a small,
%                                                          easily saturated sub-population makes
%                                                          it (NVs the laser barely re-polarises,
%                                                          or charge-state memory)
%   classic           not walk mode                        walk ordering artefacts / memory
%                                                          across blocks
%
% The AWG burst is set to window + 5 us (laserInitializationDuration = 5) in every condition, so
% it always ends before the reference window: MW overhang is excluded by construction.
%
% WHAT IT COSTS. The reference dip lowers the ESR ratio contrast by about (ratio) but does NOT
% move the line. Line positions and fields (level 3) are unaffected. Fix it only if contrast
% matters.
%
% OPTIONS
%   freqs        MHz (default 2870 + (-25:1:25))
%   amplitude    dBm (default: the json value)
%   repeats, averages   walk passes and averages (default 100, 4)
%   init         laserInitializationDuration for every condition, us (default 5)
%   conditions   cell of structs: label, and any of gap, window, ampShift (dB), switches (logical),
%                walk (logical), role ('base','noMW','gap','window','power','classic'),
%                gapAom (AOM value during the gaps; runs AWGExpESR_gapAom -- needs its cable)
%   windowAom    AOM value (laser panel units) for the detection windows during the whole scan;
%                [] = leave the panel as it is. Restored at the end.
%   gapAomPgChannel, gapAomDaqTerminal, gapAomReturnLead   passed to AWGExpESR_gapAom ([] = its default)
%   helmholtzOff default true: h.control('off') first (zero field)
%   jsonLocation, plot, data
%
% OUTPUT
%   out.cond(i)  label, role, settings, freqs, sig, sigErr, ref, refErr, aSig, aRef, ratio, ratioErr
%   out.verdict  the conclusion lines that were printed

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('data', [], 'jsonLocation', '', ...
        'freqs', 2870 + (-25:1:25), 'amplitude', [], 'repeats', 100, 'averages', 4, ...
        'init', 5, 'conditions', [], 'helmholtzOff', true, 'plot', true, ...
        'windowAom', [], 'gapAomPgChannel', [], 'gapAomDaqTerminal', '', 'gapAomReturnLead', []));

    if ~isempty(opts.data)
        out = opts.data;
        fprintf('\n=== esrRefDipScan (REPLAY of saved data, no hardware) ===\n');
    else
        if isempty(opts.jsonLocation)
            try; opts.jsonLocation = setup7Paths().json; catch; opts.jsonLocation = ''; end
        end
        conds = opts.conditions;
        if isempty(conds)
            conds = { ...
                struct('label', 'baseline: gap 10, window 100',        'role', 'base'), ...
                struct('label', 'switches never opened (no MW at NVs)', 'role', 'noMW', 'switches', false), ...
                struct('label', 'gap 300 us, laser on',                 'role', 'gap', 'gap', 300), ...
                struct('label', 'long gap: 2000 us, laser on',          'role', 'gap', 'gap', 2000), ...
                struct('label', 'long window: 200 us',                  'role', 'window', 'window', 200), ...
                struct('label', 'MW 10 dB lower',                       'role', 'power', 'ampShift', -10), ...
                struct('label', 'classic (not walk)',                   'role', 'classic', 'walk', false)};
        end
        if opts.helmholtzOff
            try
                h = getObjByName('Helmholtz'); h.control('off');
                fprintf('Helmholtz off (zero field).\n');
            catch
                fprintf('Could not switch the Helmholtz off -- make sure the field is off.\n');
            end
        end
        fprintf('\n=== esrRefDipScan: %d conditions, %d frequencies, %d passes x %d averages ===\n', ...
                numel(conds), numel(opts.freqs), opts.repeats, opts.averages);
        restoreAom = [];
        if ~isempty(opts.windowAom)
            part = AWGExpESR_gapAom.gapAomPart();
            old = part.valueInternal; if isempty(old); old = part.value; end
            part.value = opts.windowAom;
            restoreAom = onCleanup(@() setAomSafe(part, old));    % back on exit, error or not
            fprintf('Laser panel (AOM) set to %g for the windows; it goes back to %g at the end.\n', ...
                    opts.windowAom, old);
        end
        out = struct('cond', struct([]), 'opts', rmfieldSafe(opts, 'data'));
        for i = 1:numel(conds)
            c = withDefaults(conds{i}, struct('gap', 10, 'window', 100, 'ampShift', 0, ...
                'switches', true, 'walk', true, 'role', '', 'gapAom', []));
            fprintf('\n--- [%d/%d] %s ---\n', i, numel(conds), c.label);
            if isempty(c.gapAom)
                a = AWGExpESR;
            else
                a = AWGExpESR_gapAom;
            end
            if ~isempty(opts.jsonLocation); a.loadJsonParams('all', opts.jsonLocation); end
            if ~isempty(c.gapAom)
                a.gapAomValue = c.gapAom;
                if ~isempty(opts.gapAomPgChannel);   a.gapAomPgChannel = opts.gapAomPgChannel;     end
                if ~isempty(opts.gapAomDaqTerminal); a.gapAomDaqTerminal = opts.gapAomDaqTerminal; end
                if ~isempty(opts.gapAomReturnLead);  a.gapAomReturnLead = opts.gapAomReturnLead;   end
            end
            amp0 = a.amplitude;
            if ~isempty(opts.amplitude); amp0 = opts.amplitude; end
            if numel(amp0) > 1; amp0 = amp0(1); end
            a.frequency = opts.freqs;
            a.repeats = opts.repeats;  a.averages = opts.averages;
            a.amplitude = amp0 + c.ampShift;
            a.delayBetweenDetections = c.gap;
            a.detectionDuration = c.window;
            a.referenceDetectionDuration = c.window;
            a.laserInitializationDuration = opts.init;     % burst ends before the reference window
            a.isTracking = 0;
            a.useSequenceWalk = c.walk;
            a.useSequenceJump = false;
            a.useFgen = false;
            a.driveMWSwitches = c.switches;
            try; a.newOrderEveryAverage = c.walk; catch; end
            a.run;
            K  = a.currIter;
            S1 = reshape(a.signal(1, :, 1:K), [], K);
            S2 = reshape(a.signal(2, :, 1:K), [], K);
            cc = struct('label', c.label, 'role', c.role, 'gap', c.gap, 'window', c.window, ...
                        'amplitude', a.amplitude(1), 'switches', c.switches, 'walk', c.walk, ...
                        'freqs', opts.freqs(:), 'sigAvg', S1, 'refAvg', S2);
            cc.gapAom = c.gapAom;                     % [] = usual power in the gaps
            if isempty(out.cond); out.cond = cc; else; out.cond(end+1) = cc; end
        end
        clear restoreAom                              % panel back to its value now
    end

    out = analyse(out);
    report(out);
    if isempty(opts.data)
        try
            fn = fullfile(setup7Paths().folder, 'esrRefDipScan.mat');
            save(fn, 'out');
            fprintf('saved to %s   (re-analyse: esrRefDipScan(struct(''data'', out)))\n\n', fn);
        catch
        end
    end
    if opts.plot; plotAll(out); end
end

% ======================================================================
%  ANALYSIS: both traces projected onto the baseline condition's dip shape
% ======================================================================

function out = analyse(out)
% The dip is formed PER AVERAGE (each average against its own wing baseline) and only then
% averaged: the laser level drifts between averages by far more than the reference dip, and
% that drift is common to every frequency, so it cancels inside one average but not across them.
    n = numel(out.cond);
    f = out.cond(1).freqs(:);
    fc = 0.5*(min(f) + max(f));
    wing = abs(f - fc) >= 0.6*max(abs(f - fc));
    iT = find(strcmp({out.cond.role}, 'base'), 1); if isempty(iT); iT = 1; end
    [dT, ~] = dipPerAvg(out.cond(iT).sigAvg, f, wing);
    T = smooth3(dT); T(wing) = 0;
    [mx, k] = max(T); T = T / max(mx, eps);
    out.template = T; out.lineCentre = f(k);
    for i = 1:n
        c = out.cond(i);
        [d1, s1] = dipPerAvg(c.sigAvg, f, wing);
        [d2, s2] = dipPerAvg(c.refAvg, f, wing);
        [a1, ea1] = project(d1, s1, T);
        [a2, ea2] = project(d2, s2, T);
        r = a2 / a1;
        er = abs(r) * sqrt((ea1/max(abs(a1), eps))^2 + (ea2/max(abs(a2), eps))^2);
        out.cond(i).aSig = a1;
        out.cond(i).aSig = a1; out.cond(i).aSigErr = ea1;
        out.cond(i).aRef = a2; out.cond(i).aRefErr = ea2;
        out.cond(i).ratio = r; out.cond(i).ratioErr = er;
        out.cond(i).dSig = d1; out.cond(i).dRef = d2;
    end
    % no real signal dip (significant AND >10% of the baseline's): the ratio means nothing
    aT = out.cond(iT).aSig;
    for i = 1:n
        if out.cond(i).aSig < max(3*out.cond(i).aSigErr, 0.10*aT)
            out.cond(i).ratio = NaN; out.cond(i).ratioErr = NaN;
        end
    end
    out.verdict = verdict(out);
end

function [d, s] = dipPerAvg(Y, f, wing)
% Fractional dip of every average against its own straight-line wing baseline; mean and
% standard error across averages.
    K = size(Y, 2); D = zeros(size(Y));
    for k = 1:K
        p = polyfit(f(wing) - mean(f), Y(wing, k), 1);
        D(:, k) = 1 - Y(:, k) ./ polyval(p, f - mean(f));
    end
    d = mean(D, 2);
    if K >= 2
        s = std(D, 0, 2) / sqrt(K);
    else
        s = std(d(wing)) * ones(size(d));
    end
    s(~isfinite(s) | s <= 0) = median(s(isfinite(s) & s > 0));
end

function [a, ea] = project(d, s, T)
% Amplitude of the template in d (weighted least squares), with its error scaled by the
% residual scatter so it is never more optimistic than the data.
    w = 1 ./ s.^2;
    a = sum(w .* T .* d) / sum(w .* T.^2);
    ea = 1 / sqrt(sum(w .* T.^2));
    r = d - a*T; chi = sum(w .* r.^2) / max(numel(d) - 1, 1);
    ea = ea * sqrt(max(chi, 1));
end

function v = verdict(out)
    v = {};
    C = out.cond; roles = {C.role};
    ib = find(strcmp(roles, 'base'), 1);
    if isempty(ib)
        v{end+1} = 'No baseline condition, so only the table applies.'; return
    end
    b = C(ib);
    v{end+1} = sprintf(['Baseline: signal dip %.2f%%, reference dip %.2f%%, ratio %.3f +/- %.3f. It lowers the ' ...
                        'ESR contrast by ~%.0f%% and does not move the line.'], 100*b.aSig, 100*b.aRef, ...
                        b.ratio, b.ratioErr, 100*b.ratio);
    rel = @(c) c.ratio / b.ratio;
    relErr = @(c) abs(c.ratio / b.ratio) * sqrt((c.ratioErr/max(abs(c.ratio), eps))^2 + (b.ratioErr/max(abs(b.ratio), eps))^2);
    vague = @(c) relErr(c) > 0.25;
    inconc = @(c, what) sprintf(['%s: %.0f%% +/- %.0f%% of baseline -- too noisy to tell. Raise averages (or the MW ' ...
                                 'power for a deeper dip) and rerun this condition.'], what, 100*rel(c), 100*relErr(c));
    i = find(strcmp(roles, 'noMW'), 1);
    if ~isempty(i)
        c = C(i);
        % a dip counts only if it is significant AND a real fraction of the baseline dip
        dS = c.aSig > max(3*c.aSigErr, 0.10*b.aSig);
        dR = c.aRef > max(3*c.aRefErr, 0.25*b.aRef);
        if dS && dR
            v{end+1} = sprintf(['Switches never opened, and BOTH traces still dip (%.2f%% / %.2f%%): MW reaches the NVs ' ...
                                'through the closed switches. Leakage is the cause -- check switch isolation.'], ...
                               100*c.aSig, 100*c.aRef);
        elseif dR
            v{end+1} = sprintf(['Switches never opened: no signal dip, but the REFERENCE still dips (%.2f%%). Not NV ' ...
                                'physics -- the reference value is being contaminated in the processing/walk reduction.'], ...
                               100*c.aRef);
        else
            v{end+1} = 'Switches never opened: no dip in either trace. The reference dip needs MW at the NVs during the signal window.';
        end
    end
    ig = find(strcmp(roles, 'gap') & arrayfun(@(c) isfinite(c.ratio), C));
    if ~isempty(ig)
        [~, o] = sort([C(ig).gap]); ig = ig(o);
        gone = arrayfun(@(c) ~vague(c) && rel(c) < 0.5, C(ig));
        if any(gone)
            % time constant from every gap point that still has a measurable dip, plus baseline
            gp = [b.gap, C(ig).gap]; rp = [b.ratio, C(ig).ratio]; ep = [b.ratioErr, C(ig).ratioErr];
            use = rp > 2*ep;
            if sum(use) >= 2
                w = (rp(use) ./ ep(use)).^2;
                pf = wpolyfit(gp(use), log(rp(use)), w);
                tau = -1/pf(1);
                v{end+1} = sprintf(['A laser-on gap shrinks it (%s): a memory the laser erases, tau ~ %.0f us. ' ...
                                    'Removing it needs delayBetweenDetections ~ %.0f us (3 tau).'], gapList(C(ig), rel), tau, 3*tau);
                v{end+1} = costLine(3*tau, b);
            else
                g1 = C(ig(find(gone, 1))).gap;
                v{end+1} = sprintf(['It is gone by a %g us laser-on gap (%s): a memory the laser erases, tau below ~%.0f us. ' ...
                                    'Use delayBetweenDetections ~ %g us, or scan shorter gaps to pin tau.'], ...
                                   g1, gapList(C(ig), rel), g1/3, g1);
                v{end+1} = costLine(g1, b);
            end
        elseif all(arrayfun(@(c) ~vague(c) && abs(rel(c) - 1) < max(2*relErr(c), 0.2), C(ig)))
            v{end+1} = sprintf(['A laser-on gap does not change it (%s): whatever carries the dip survives ' ...
                                'milliseconds of laser (a very long gap, ~10 ms, has removed it before, so tau is ' ...
                                'several ms). Too slow to fix with a gap.'], gapList(C(ig), rel));
            v{end+1} = costLine(10000, b);
        elseif any(arrayfun(@(c) vague(c), C(ig)))
            v{end+1} = inconc(C(ig(end)), sprintf('Gap %g us', C(ig(end)).gap));
        else
            v{end+1} = sprintf('A laser-on gap changes it partly (%s): part memory, part something slower.', gapList(C(ig), rel));
        end
    end
    i = find(strcmp(roles, 'window'), 1);
    if ~isempty(i) && isfinite(C(i).ratio)
        c = C(i);
        if vague(c)
            v{end+1} = inconc(c, sprintf('Window %g us', c.window));
        elseif rel(c) < 0.7 && (1 - rel(c)) > 2*relErr(c)
            v{end+1} = sprintf('Window %g us dilutes it to %.0f%%: it sits at the start of the reference window.', c.window, 100*rel(c));
        else
            v{end+1} = sprintf('Window %g us does not dilute it (%.0f%%): it is spread through the reference window, not at its start.', ...
                               c.window, 100*rel(c));
        end
    end
    i = find(strcmp(roles, 'power'), 1);
    if ~isempty(i) && isfinite(C(i).ratio)
        c = C(i);
        if vague(c) && abs(rel(c) - 1) < 3*relErr(c)
            v{end+1} = inconc(c, 'MW 10 dB lower');
        elseif rel(c) > 1.3 && (rel(c) - 1) > 2*relErr(c)
            v{end+1} = sprintf(['At lower MW power the ratio GROWS (%.0f%% of baseline): a small sub-population that saturates ' ...
                                'easily -- NVs the laser barely re-polarises (out of focus) or a charge-state memory. A ' ...
                                'property of the sample/optics: a smaller pinhole, better focus or more laser power reduce it.'], ...
                               100*rel(c));
        elseif rel(c) < 0.7 && (1 - rel(c)) > 2*relErr(c)
            v{end+1} = sprintf('At lower MW power the ratio SHRINKS (%.0f%% of baseline): it scales faster than the signal with power.', 100*rel(c));
        else
            v{end+1} = sprintf('MW power does not change the ratio (%.0f%% of baseline).', 100*rel(c));
        end
    end
    i = find(strcmp(roles, 'classic'), 1);
    if ~isempty(i) && isfinite(C(i).ratio)
        c = C(i);
        if vague(c)
            v{end+1} = inconc(c, 'Classic mode');
        elseif abs(rel(c) - 1) > max(2*relErr(c), 0.2)
            v{end+1} = sprintf(['Classic mode gives %.0f%% of the walk value: the walk ordering matters, so the memory spans ' ...
                                'blocks (frequencies) -- something slower than one repeat.'], 100*rel(c));
        else
            v{end+1} = 'Classic and walk agree: not an artefact of the walk ordering.';
        end
    end
end

% ======================================================================

function report(out)
    fprintf('\n%-40s %12s %12s %16s\n', 'condition', 'signal dip', 'ref dip', 'ref / signal');
    for i = 1:numel(out.cond)
        c = out.cond(i);
        if isfinite(c.ratio); rs = sprintf('%.3f +/- %.3f', c.ratio, c.ratioErr); else; rs = '--'; end
        fprintf('%-40s %10.3f%%  %10.3f%%  %16s\n', c.label, 100*c.aSig, 100*c.aRef, rs);
    end
    fprintf('\nWHAT IT SAYS\n');
    for i = 1:numel(out.verdict); fprintf('   %s\n', out.verdict{i}); end
    fprintf('\n');
end

function plotAll(out)
    cols = [0.15 0.42 0.69; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56; 0.10 0.10 0.10; ...
            0.30 0.75 0.93; 0.64 0.08 0.18];
    figure('Color', 'w', 'Position', [80 80 1250 460]);
    for p = 1:2
        ax = subplot(1, 3, p); hold(ax, 'on'); set(ax, 'Box', 'off', 'TickDir', 'out', ...
            'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.12, 'FontSize', 9);
        hs = []; ls = {};
        for i = 1:numel(out.cond)
            c = out.cond(i); col = cols(1 + mod(i-1, size(cols, 1)), :);
            if p == 1; y = c.dSig; else; y = c.dRef; end
            hs(end+1) = plot(ax, c.freqs, 100*y, '.-', 'Color', col); %#ok<AGROW>
            ls{end+1} = c.label; %#ok<AGROW>
        end
        xlabel(ax, 'frequency (MHz)'); ylabel(ax, 'dip (% of baseline)');
        if p == 1
            title(ax, 'signal (with MW) dip', 'FontWeight', 'bold');
        else
            title(ax, 'reference (without MW) dip', 'FontWeight', 'bold');
            legend(ax, hs, ls, 'Location', 'best', 'Box', 'off');
        end
    end
    ax = subplot(1, 3, 3); hold(ax, 'on'); set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 9, ...
        'XGrid', 'on', 'GridAlpha', 0.12);
    n = numel(out.cond);
    for i = 1:n
        c = out.cond(i); col = cols(1 + mod(i-1, size(cols, 1)), :);
        if isfinite(c.ratio)
            barh(ax, n + 1 - i, c.ratio, 0.6, 'FaceColor', col, 'EdgeColor', 'none');
            plot(ax, c.ratio + [-1 1]*c.ratioErr, (n + 1 - i)*[1 1], 'k-', 'LineWidth', 1.2);
        end
    end
    set(ax, 'YTick', 1:n, 'YTickLabel', fliplr(1:n));
    xlabel(ax, 'reference dip / signal dip'); ylabel(ax, 'condition # (legend order)');
    title(ax, 'ratio per condition', 'FontWeight', 'bold');
end

% ======================================================================

function s = costLine(gap, b)
% What a gap costs: each repeat is 2*(window + gap) long but only the window measures, so the
% averaging time for the same SNR scales with (window + gap).
    f = (b.window + gap) / (b.window + b.gap);
    if f > 3
        loss = 'the small';
        if isfinite(b.ratio); loss = sprintf('the ~%.0f%%', 100*b.ratio); end
        s = sprintf(['A %.0f us gap makes every ESR ~%.0fx slower for the same SNR -- not worth it. The dip ' ...
                     'does not move the line: accept %s contrast loss, or reduce it optically (tighter ' ...
                     'pinhole, better focus, more laser power: the slow NVs are the ones that see little laser).'], ...
                    gap, f, loss);
    else
        s = sprintf('A %.0f us gap costs only ~%.1fx in averaging time: worth using.', gap, f);
    end
end

function s = gapList(C, rel)
    p = arrayfun(@(c) sprintf('%g us: %.0f%%', c.gap, 100*rel(c)), C, 'UniformOutput', false);
    s = strjoin(p, ', ');
end

function p = wpolyfit(x, y, w)
    x = x(:); y = y(:); w = w(:);
    A = [x, ones(size(x))] .* sqrt(w);
    p = (A \ (y .* sqrt(w)))';
end

function y = smooth3(x)
    x = x(:); y = x;
    for k = 2:numel(x)-1; y(k) = mean(x(k-1:k+1)); end
end

function setAomSafe(part, v)
    try
        part.value = v;
        fprintf('Laser panel (AOM) back to %g.\n', v);
    catch err
        fprintf('Could not set the laser panel back to %g (%s) -- set it by hand.\n', v, err.message);
    end
end

function s = rmfieldSafe(s, f)
    if isfield(s, f); s = rmfield(s, f); end
end

function o = withDefaults(o, d)
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(o, f{i}) || (isempty(o.(f{i})) && ~any(strcmp(f{i}, {'data', 'conditions', 'amplitude', 'gapAom', 'windowAom'})))
            o.(f{i}) = d.(f{i});
        end
    end
end
