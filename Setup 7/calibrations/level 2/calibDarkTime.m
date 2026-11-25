function d = calibDarkTime(opts)
%CALIBDARKTIME Contrast vs dark time before the readout, under several hardware conditions.
%
%   d = calibDarkTime()                                  % the default diagnostic set
%   d = calibDarkTime(struct('conditions', {{...}}))     % your own set
%   d = calibDarkTime(struct('data', d))                 % refit / replot, no hardware
%
% LEVEL 2 DIAGNOSTIC (level2.m cell 2.1b). Not a calibration: it writes nothing and returns no
% json value. It exists because of what the long lastDelay scans showed on 2026-11-04:
%
%   on resonance   contrast 1.12% -> 1.77% as the gap goes 0 -> 150 us
%   off resonance  contrast 0.14% -> 0.72%, same shape          (2600 MHz, no spin contrast)
%
% So the slow rise needs no spins: the real spin contrast (on - off) is flat at ~1.05% after the
% MW step. What rises is a HARDWARE ARTEFACT that exists only in the signal half, and it grows
% with how long the MW switches were open before the window -- they open at the trigger and,
% in the normal sequence, close exactly when the window opens. A switch control line coupling
% into a 22 mV photodiode signal through an AC-coupled path gives exactly that shape.
%
% Each condition below changes ONE thing about the signal half (AWGExpLastDelayDuration's
% isolation switches) and reruns the same gap scan. Which conditions remove the rise says what
% causes it:
%
%   off resonance                      the artefact alone (baseline for the others)
%   + mwOff       AWG plays zeros      still rises  -> switching/triggering, not MW power
%   + gapSwitch=0 switch shuts after   flat         -> the switch edge at the window is it
%                 the burst
%   + refSwitch   reference half       flat         -> it cancels when both halves switch
%                 switches identically                 alike: THE FIX for every pulsed class
%   ON resonance + refSwitch           the true spin contrast with the fix: should be flat
%
% It matters beyond lastDelay: in every balanced pulsed sequence the switches open only in the
% signal half, so the same artefact sits on the contrast of Rabi, Ramsey, Echo and T1, sized by
% how long the switch is open there -- i.e. it changes with tau.
%
% OPTIONS
%   gaps          dark times to scan, us (default: 0.2 ... 200, log-spaced, 18 points)
%   conditions    cell of structs: field 'label' plus any calibLastDelay option to override
%                 (frequency, mwOff, gapSwitch, refSwitch, noTrigger, ...).
%   preset        used when conditions is empty:
%                   'balanced' (default) off resonance UNBALANCED, off resonance balanced, ON
%                              resonance balanced. Round 3, see below.
%                   'trigger'  AWG not triggered; and no trigger + reference switches too.
%                   'hardware' the five conditions above (round 1).
%   balancedSequence  default 1 for every condition unless the condition overrides it.
%
% ROUND 1 RESULT (2026-11-04): all four off-resonance conditions rose IDENTICALLY, 0.17 -> 0.85%
% (t63 ~29 us): MW off, switch shut early and reference switching too changed nothing. So it is
% neither the MW power nor the MW switch. On resonance (with the AWG fix) the contrast is now
% 2.07% at short gaps, about double the 1.12% before -- the alternating-sequence bug was real.
% ROUND 2: AWG never triggered -- still 0.565%. No trigger + reference switching -- still 0.560%.
%
% THE ANSWER: none of it was hardware. Experiment's balancedSequence defaults to 0, the lastDelay
% scan never set it, and the UNBALANCED sequence puts the reference window straight after
% continuous laser while the signal window follows TRIG + burst + gap of dark. Fluorescence after
% a dark period is lower than under steady light (0.14% after ~1 us, ~0.7% after 100 us), and
% that read as contrast. refSwitch did nothing because it only acts in the balanced branch -- so
% round 2's "halves identical" were not identical. calibLastDelay is now balanced by default;
% round 3 confirms the rise disappears. The real experiments are not affected: Rabi, Ramsey and
% Echo run balanced, and T1's contrast is |0> - |1> after the same dark.
%   frequency     on-resonance frequency, MHz (default 2870)
%   offFrequency  MHz (default 2600)
%   repeats, averages (1000, 10), detectionDuration (3 us), laserInitializationDuration (100 us)
%   fitFrom       gaps below this are shown but not fitted, so the MW step stays out (0.3 us)
%   data          a previous result: refit and replot only
%   jsonLocation, plot
%
% OUTPUT
%   d.cond(i)     label, gaps, contrast, contrastErr, fit (c0, rise, riseErr, t63, t90, flat)
%   d.verdict     cell of the conclusion lines that were printed

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('data', [], 'jsonLocation', '', ...
        'gaps', [0.2 0.5 1 2 3 5 7 10 15 20 30 40 60 80 100 130 160 200], ...
        'conditions', [], 'frequency', 2870, 'offFrequency', 2600, ...
        'repeats', 1000, 'averages', 10, 'detectionDuration', 3, ...
        'laserInitializationDuration', 100, 'fitFrom', 0.3, 'preset', 'balanced', ...
        'balancedSequence', 1, 'plot', true));

    if ~isempty(opts.data)
        d = opts.data;
        fprintf('\n=== calibDarkTime (REPLAY of saved data, no hardware) ===\n');
    else
        conds = opts.conditions;
        off = opts.offFrequency; on = opts.frequency;
        if isempty(conds) && strcmp(opts.preset, 'balanced')
            % Round 3: the cause was the UNBALANCED sequence (reference after laser light, signal
            % after dark). Old scan next to the balanced one, and the spin contrast balanced.
            conds = { ...
                struct('label', 'off resonance, UNBALANCED (the old scan)', 'frequency', off, ...
                       'balancedSequence', 0), ...
                struct('label', 'off resonance, balanced', 'frequency', off, 'balancedSequence', 1), ...
                struct('label', sprintf('ON resonance (%g MHz), balanced', on), 'frequency', on, ...
                       'balancedSequence', 1)};
        elseif isempty(conds) && strcmp(opts.preset, 'trigger')
            % Round 2 (after 2026-11-04 round 1 ruled out MW power and the switch): the only
            % signal-half-only thing left is the AWG trigger -- or the readout itself.
            conds = { ...
                struct('label', 'AWG not triggered (no MW at all)', 'frequency', off, 'noTrigger', true), ...
                struct('label', 'halves identical: no trigger, reference switches too', ...
                       'frequency', off, 'noTrigger', true, 'refSwitch', true)};
        elseif isempty(conds)
            conds = { ...
                struct('label', sprintf('off resonance (%g MHz)', off), 'frequency', off), ...
                struct('label', 'off res, MW off (AWG zeros)', 'frequency', off, 'mwOff', true), ...
                struct('label', 'off res, switch shuts after burst', 'frequency', off, 'gapSwitch', false), ...
                struct('label', 'off res, reference switches too', 'frequency', off, 'refSwitch', true), ...
                struct('label', sprintf('ON resonance (%g MHz), reference switches too', on), ...
                       'frequency', on, 'refSwitch', true)};
        end
        g = opts.gaps(:)';
        perRep = 2*(1.1 + g + opts.laserInitializationDuration);          % us per repeat
        estMin = (sum(perRep)*opts.repeats*1e-6 + 0.5*numel(g)) * opts.averages / 60;
        fprintf('\n=== calibDarkTime ===\n');
        fprintf('%d conditions x %d gaps (%.2g to %.3g us), %d x %d. Roughly %.0f min per condition.\n', ...
                numel(conds), numel(g), g(1), g(end), opts.repeats, opts.averages, estMin);

        d = struct('cond', struct('label', {}, 'overrides', {}, 'gaps', {}, 'contrast', {}, ...
                                  'contrastErr', {}, 'contrastPerAvg', {}), 'opts', opts);
        for i = 1:numel(conds)
            cnd = conds{i};
            fprintf('\n--- [%d/%d] %s ---\n', i, numel(conds), cnd.label);
            o = struct('gaps', g, 'repeats', opts.repeats, 'averages', opts.averages, ...
                       'detectionDuration', opts.detectionDuration, ...
                       'laserInitializationDuration', opts.laserInitializationDuration, ...
                       'jsonLocation', opts.jsonLocation, 'frequency', opts.frequency, ...
                       'balancedSequence', opts.balancedSequence, 'analyse', false, 'plot', false);
            f = fieldnames(cnd);
            for j = 1:numel(f)
                if ~strcmp(f{j}, 'label'); o.(f{j}) = cnd.(f{j}); end
            end
            q = calibLastDelay(o);
            pa = [];
            if isfield(q, 'contrastPerAvg'); pa = q.contrastPerAvg; end
            d.cond(i) = struct('label', cnd.label, 'overrides', rmfieldSafe(cnd, 'label'), ...
                               'gaps', q.gaps(:), 'contrast', q.contrast(:), ...
                               'contrastErr', q.contrastErr(:), 'contrastPerAvg', pa);
        end
    end

    % ---- fit every curve ----------------------------------------------------
    for i = 1:numel(d.cond)
        cd = d.cond(i);
        k = cd.gaps >= opts.fitFrom & isfinite(cd.contrast) & isfinite(cd.contrastErr);
        fit = satFit(cd.gaps(k), cd.contrast(k), cd.contrastErr(k));
        fits(i) = fit; %#ok<AGROW>
    end
    for i = 1:numel(d.cond); d.cond(i).fit = fits(i); end

    % ---- table --------------------------------------------------------------
    fprintf('\n%-46s %10s %10s %16s %8s %8s\n', 'condition', 'C(first)', 'C(last)', ...
            'rise', 't63 us', 't90 us');
    for i = 1:numel(d.cond)
        F = d.cond(i).fit;
        if F.flat; tag = '  flat'; else; tag = ''; end
        fprintf('%-46s %9.3f%% %9.3f%% %+8.3f +/- %.3f%% %8s %8s%s\n', d.cond(i).label, ...
                100*F.cFirst, 100*F.cLast, 100*F.rise, 100*F.riseErr, num(F.t63), num(F.t90), tag);
    end

    % ---- verdict ------------------------------------------------------------
    d.verdict = verdict(d);
    fprintf('\nWHAT IT SAYS\n');
    for i = 1:numel(d.verdict); fprintf('   %s\n', d.verdict{i}); end
    fprintf('\n');

    if opts.plot; plotAll(d, opts); end
end

% ======================================================================

function v = verdict(d)
% One line per condition, plus the conclusion the combination supports.
    v = {};
    n = numel(d.cond);
    A = repmat(struct('off', false, 'bal', false, 'iso', false, 'noMW', false, 'onRes', false), 1, n);
    defBal = 0;     % runs saved before 2026-11-04 evening were unbalanced and do not record it
    if isfield(d, 'opts') && isfield(d.opts, 'balancedSequence'); defBal = d.opts.balancedSequence; end
    for i = 1:n
        ov = d.cond(i).overrides;
        A(i).bal = defBal;
        if isfield(ov, 'balancedSequence'); A(i).bal = logical(ov.balancedSequence); end
        % an isolation switch counts only if it is set AWAY from its normal value
        A(i).iso = any(cellfun(@(f) isfield(ov, f) && logical(ov.(f)) ~= strcmp(f, 'gapSwitch'), ...
                   {'mwOff', 'refSwitch', 'noTrigger', 'gapSwitch'}));
        A(i).noMW = (isfield(ov, 'noTrigger') && ov.noTrigger) || (isfield(ov, 'mwOff') && ov.mwOff);
        A(i).onRes = isfield(ov, 'frequency') && abs(ov.frequency - d.opts.frequency) < 50 && ~A(i).noMW;
        A(i).off = ~A(i).onRes;
    end
    if any(~[A.bal])
        v{end+1} = ['UNBALANCED conditions (reference window straight after laser light): any rise there ' ...
                    'is the lower fluorescence after a dark period, not the microwave.'];
    end
    iUB = find([A.off] & ~[A.bal] & ~[A.iso], 1);
    iB  = find([A.off] &  [A.bal] & ~[A.iso], 1);
    iOn = find([A.onRes] & [A.bal] & ~[A.iso], 1);
    if ~isempty(iUB)
        F = d.cond(iUB).fit;
        v{end+1} = sprintf('Off resonance, unbalanced: rises by %.3f%% (t63 %s us) -- the fluorescence-after-dark effect.', ...
                           100*F.rise, num(F.t63));
    end
    if ~isempty(iB)
        F = d.cond(iB).fit;
        if F.flat
            v{end+1} = ['Off resonance, balanced: FLAT. The slow rise was the unbalanced sequence, not hardware. ' ...
                        'Rabi/Ramsey/Echo (balanced) and T1 (|0>-|1>) are not affected.'];
        else
            v{end+1} = sprintf('Off resonance, balanced: still rises by %.3f%% -- a real signal-half artefact remains.', 100*F.rise);
        end
    end
    if ~isempty(iOn)
        F = d.cond(iOn).fit;
        if F.flat
            v{end+1} = sprintf(['On resonance, balanced: spin contrast %.3f%%, flat in dark time. lastDelay is set by ' ...
                                'the MW step alone.'], 100*F.c0);
        else
            v{end+1} = sprintf(['On resonance, balanced: the contrast changes by %.3f%% across the scan -- the spin ' ...
                                'contrast itself depends on the dark time.'], 100*F.rise);
        end
    end
    % isolation conditions, kept for completeness
    for i = 1:n
        if ~A(i).iso; continue; end
        ov = d.cond(i).overrides; F = d.cond(i).fit;
        tag = ''; if ~A(i).bal; tag = ' (unbalanced: refSwitch has no effect there)'; end
        what = strjoin(fieldnames(rmfieldSafe(rmfieldSafe(ov, 'frequency'), 'balancedSequence'))', '+');
        if F.flat
            v{end+1} = sprintf('%s: flat.%s', what, tag); %#ok<AGROW>
        else
            v{end+1} = sprintf('%s: rises by %.3f%%.%s', what, 100*F.rise, tag); %#ok<AGROW>
        end
    end
    if isempty(v); v = {'No condition here that the rules know; read the table.'}; end
end

function tf = isOffPlain(ov)
% The plain off-resonance condition: a frequency override and nothing else.
    f = fieldnames(ov);
    tf = numel(f) == 1 && strcmp(f{1}, 'frequency');
end

% ======================================================================

function F = satFit(g, C, dC)
% Flat, or c + sum of 1-2 exponentials in g (a saturating rise or fall), chosen by BIC.
% Time constants by fminsearch, everything else linear (variable projection).
    g = g(:); C = C(:); w = 1 ./ dC(:).^2; N = numel(g);
    gmin = min(g); gmax = max(g);
    [~, c0] = wls(ones(N, 1), C, w);
    best = struct('n', 0, 'tau', [], 'bic', c0 + log(N), 'chi', c0);
    o = optimset('Display', 'off', 'MaxFunEvals', 2000, 'MaxIter', 2000, 'TolX', 1e-6, 'TolFun', 1e-6);
    lo = log(max(0.3*min(diff(sort(g))), 0.02)); hi = log(5*gmax);
    for n = 1:2
        if N < 2*n + 3; break; end
        if n == 1
            starts = {log(gmax/30), log(gmax/10), log(gmax/3)};
        else
            starts = {log([gmax/50 gmax/5]), log([gmax/20 gmax/2]), log([gmax/100 gmax/10])};
        end
        f = @(th) vpChi(th, g, C, w, lo, hi);
        for s = 1:numel(starts)
            th = starts{s}; th = min(max(th, lo + 0.01), hi - 0.01);
            th = fminsearch(f, th, o);
            c2 = f(th);
            bic = c2 + (2*n + 1)*log(N);
            if bic < best.bic
                best = struct('n', n, 'tau', sort(exp(th)), 'bic', bic, 'chi', c2);
            end
        end
    end
    B = [ones(N, 1), exp(-g ./ best.tau(:)')];
    [x, c2] = wls(B, C, w);
    dof = max(N - (2*best.n + 1), 1);
    s2  = max(c2/dof, 1);
    covx = pinv((B .* sqrt(w))' * (B .* sqrt(w))) * s2;
    model = @(gg) [ones(numel(gg), 1), exp(-gg(:) ./ best.tau(:)')] * x;
    dvec = [0, exp(-gmax ./ best.tau(:)') - exp(-gmin ./ best.tau(:)')];   % d(rise)/dx
    rise = dvec * x;  riseErr = sqrt(max(dvec * covx * dvec', 0));
    if best.n == 0; rise = 0; riseErr = sqrt(2)*median(dC); end
    % times to 63% / 90% of the rise across the scan
    gg = linspace(gmin, gmax, 4000)'; mv = model(gg);
    t63 = crossT(gg, mv, mv(1) + 0.63*(mv(end) - mv(1)));
    t90 = crossT(gg, mv, mv(1) + 0.90*(mv(end) - mv(1)));
    flat = best.n == 0 || abs(rise) < max(3*riseErr, 3*median(dC));
    F = struct('n', best.n, 'tau', best.tau, 'coef', x, 'chi2nu', c2/dof, 'model', model, ...
               'cFirst', model(gmin), 'cLast', model(gmax), 'c0', mean(C), ...
               'rise', rise, 'riseErr', riseErr, 't63', t63, 't90', t90, 'flat', flat);
    if flat; F.t63 = NaN; F.t90 = NaN; end
end

function c2 = vpChi(th, g, C, w, lo, hi)
    if any(th < lo | th > hi); c2 = Inf; return; end
    if numel(th) > 1 && abs(th(1) - th(2)) < 0.05; c2 = Inf; return; end   % two equal taus
    [~, c2] = wls([ones(numel(g), 1), exp(-g ./ exp(th(:)'))], C, w);
end

function [x, c2] = wls(A, y, w)
    sw = sqrt(w(:)); Aw = A .* sw; yw = y(:) .* sw;
    x = Aw \ yw; r = yw - Aw*x; c2 = r'*r;
end

function t = crossT(gg, mv, level)
    if mv(end) >= mv(1); k = find(mv >= level, 1); else; k = find(mv <= level, 1); end
    if isempty(k); t = NaN; else; t = gg(k); end
end

function s = num(v)
    if isfinite(v); s = sprintf('%.1f', v); else; s = '--'; end
end

% ======================================================================

function plotAll(d, opts)
    cols = [0.15 0.42 0.69; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56; 0.10 0.10 0.10; ...
            0.30 0.75 0.93; 0.64 0.08 0.18];
    figure('Color', 'w', 'Position', [80 80 1250 480]);
    ax = subplot(1, 3, [1 2]); hold(ax, 'on');
    set(ax, 'XScale', 'log', 'Box', 'off', 'TickDir', 'out', 'XGrid', 'on', 'YGrid', 'on', ...
            'GridAlpha', 0.12, 'FontSize', 9);
    hs = []; ls = {};
    for i = 1:numel(d.cond)
        c = d.cond(i); col = cols(1 + mod(i-1, size(cols, 1)), :);
        h = errorbar(ax, c.gaps, 100*c.contrast, 100*c.contrastErr, '.');
        set(h, 'Color', col, 'MarkerSize', 12);
        try; set(h, 'CapSize', 0); catch; end
        gg = exp(linspace(log(max(min(c.gaps(c.gaps >= opts.fitFrom)), 1e-3)), log(max(c.gaps)), 300))';
        plot(ax, gg, 100*c.fit.model(gg), '-', 'Color', col, 'LineWidth', 1.2);
        hs(end+1) = h; ls{end+1} = c.label; %#ok<AGROW>
    end
    xlabel(ax, 'dark time before the window (gap, us)'); ylabel(ax, 'contrast (%)');
    legend(ax, hs, ls, 'Location', 'best', 'Box', 'off');
    title(ax, 'contrast vs dark time: which condition removes the rise?', 'FontWeight', 'bold');

    ax2 = subplot(1, 3, 3); hold(ax2, 'on');
    set(ax2, 'Box', 'off', 'TickDir', 'out', 'XGrid', 'on', 'GridAlpha', 0.12, 'FontSize', 9);
    n = numel(d.cond);
    for i = 1:n
        F = d.cond(i).fit; col = cols(1 + mod(i-1, size(cols, 1)), :);
        barh(ax2, n + 1 - i, 100*F.rise, 0.6, 'FaceColor', col, 'EdgeColor', 'none');
        plot(ax2, 100*(F.rise + [-1 1]*F.riseErr), (n + 1 - i)*[1 1], 'k-', 'LineWidth', 1.2);
    end
    set(ax2, 'YTick', 1:n, 'YTickLabel', fliplr(1:n));
    plot(ax2, [0 0], [0.4 n + 0.6], 'k-');
    xlabel(ax2, 'rise across the scan (% contrast)'); ylabel(ax2, 'condition # (legend order)');
    title(ax2, 'rise per condition', 'FontWeight', 'bold');
end

% ======================================================================

function s = rmfieldSafe(s, f)
    if isfield(s, f); s = rmfield(s, f); end
end

function o = withDefaults(o, d)
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(o, f{i}) || (isempty(o.(f{i})) && ~any(strcmp(f{i}, {'data', 'conditions'})))
            o.(f{i}) = d.(f{i});
        end
    end
end
