function out = calibAmplitudeDiff(opts)
%CALIBAMPLITUDEDIFF Measure both channels' gain curves, so amplitudeDiff is known at any drive.
%
%   out = calibAmplitudeDiff
%   out = calibAmplitudeDiff(struct('amplitudes', -30:3:-9))
%
% LEVEL 4, after calibMwPair. Coils on, at the field calibB1Rabi chose.
%
% calibMwPair measures amplitudeDiff at ONE amplitude, and that is all it can honestly do:
% each channel has its own amplifier, the two compress at different points, so a delta
% measured at -30 dBm need not hold at -16. This scans the amplitude and produces the whole
% relationship, which amplitudeDiffAt then evaluates at any drive.
%
% Writes results\<date>\amplitudeDiffCurve.mat. Nothing else is changed.
%
% ===================================================================
% amplitudeDiff IS A HORIZONTAL OFFSET, NOT A VERTICAL GAP
% ===================================================================
%
% This is the whole reason the scan measures two CURVES rather than one delta per point.
%
% What we want is the drive offset that makes the two channels deliver the same field:
%
%       Omega1(a)  =  Omega2(a + delta)
%
% Read that carefully -- Omega2 is evaluated at a DIFFERENT amplitude. The obvious
% alternative, taking 20*log10(Omega1(a)/Omega2(a)) at each point and calling it delta, is
% the vertical gap between the curves, and the two agree only when both curves have unit
% slope. Unit slope is exactly what stops being true once an amplifier compresses.
%
% Writing V for the vertical gap and s2 for the local slope of channel 2 in dB per dB,
%
%       delta  ~=  V / s2
%
% so at s2 = 0.7 a 0.9 dB gap needs 1.29 dB of offset. That 0.4 dB error is eight times the
% 0.05 dB the Rabi measurement is actually capable of. The vertical gap is not a slightly
% worse estimate of delta; below the compression knee it is a different quantity.
%
% So: measure Omega1(a) and Omega2(a) over a range, hand both curves to amplitudeDiffAt,
% and let it solve. Nothing here averages a per-point ratio.
%
% ===================================================================
% WHY EVERY POINT COSTS A FULL UPLOAD
% ===================================================================
%
% On this AWG the output amplitude is NOT a front-panel setting. setAmplitude pins the
% instrument at MAX_VPP and folds the requested power into MaxNormPowerAllowed, which is
% the factor multiplying the cosine when the waveforms are built:
%
%       Power_Volt = obj.dBmToVpp(Power);
%       obj.MaxNormPowerAllowed = Power_Volt/obj.MAX_VPP;      % scales the SAMPLES
%
% So amplitude lives in the samples, every point re-uploads, and there is no way to scan it
% cheaply. It also means amplitude MUST stay in the load-cache fingerprint: dropping it
% looks like an easy speed-up and would silently reuse waveforms carrying the old drive
% level while the code believed the new one.
%
% What can be saved is the tau list. The scan walks from HIGH amplitude to LOW and sizes
% each sweep from the previous point's Rabi frequency, so a fast Rabi gets a short sweep
% and a small upload, and a slow one gets a long enough sweep to show a period.
%
% ===================================================================
% THE BOTTOM OF THE RANGE WILL FAIL, AND THAT IS PHYSICS
% ===================================================================
%
% Omega scales as 10^(a/20), so from -16 dBm (about 8 MHz per channel here) down to -40 is
% a factor of 16 -- roughly 0.5 MHz, a 2 us period. That is comparable to the Rabi decay,
% so less than one oscillation survives and no fit can recover a frequency from it.
%
% Expect the lowest one or two points to be rejected. They are reported as "too slow to
% fit at this drive" rather than silently dropped, because a rejected point at the bottom
% is expected while a rejected point in the middle means something is wrong.
%
% ===================================================================
% THE COMBINED RUN, AND WHY ITS BALANCE DOES NOT MATTER
% ===================================================================
%
% A third run per point drives both channels and gives
%
%       R(a) = Omega_both / (Omega1 + Omega2)
%
% which is 1 for a linear network at the optimum phase and falls once something shared
% compresses. That is what fixes the top of the usable range, and it cannot be inferred
% from the single-channel curves -- each amplifier alone may be perfectly happy.
%
% The combined run needs SOME delta, and the properly solved one is not available until the
% whole scan is finished. It uses the vertical-gap value, and that is not an approximation
% worth worrying about: at the optimum phase Omega_both = Omega1 + Omega2 exactly, whatever
% the imbalance, so R is completely insensitive to how well the two are matched. The
% balance affects how much drive you get, not whether the two add.
%
% ===================================================================
% OPTIONS
% ===================================================================
%
%   amplitudes    dBm for channel 1, scanned HIGH to LOW whatever order you give
%                 (default linspace(-40, -10, 10))
%   channelPhase  degrees on channel 2 for the combined run (default [] = from the json,
%                 which calibMwPair wrote)
%   doCombined    the third run per point, for R(a) (default true)
%   nRepsPerPoint how many times to measure each channel at each amplitude (default 1;
%                 raise to 2 if the drift warning fires)
%   offAmplitude  dBm for the channel held off (default -70, i.e. 54 dB down = 0.2 % in
%                 amplitude)
%   tauMaxCap     longest sweep the scan will build, us (default 4). The low-amplitude
%                 points want a long sweep and there is no point going past the decay.
%   nTau          points per sweep (default 60). Fewer than the standard 99 because the
%                 upload scales with it and only the frequency is wanted here.
%   save          write amplitudeDiffCurve.mat into the results folder (default true)
%
%   repeats, averages, rabiDetectionMHz, nFreq, minContrast, maxFErrFrac, isTracking,
%   plotIt, jsonLocation      as calibMwPair
%   newOrderEveryAverage  walk-mode Rabi runs: a new random tau order before every average
%                  (true/false). Default [] = AWGExpRabi's own setting (false). Costs one
%                  sequence upload + arm per average.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct( ...
        'amplitudes', [], 'channelPhase', [], 'doCombined', true, ...
        'nRepsPerPoint', 1, 'offAmplitude', -70, 'tauMaxCap', 4, 'nTau', 60, ...
        'save', true, 'repeats', 300, 'averages', 2, 'rabiDetectionMHz', 0.3, ...
        'nFreq', 1, 'minContrast', 0.005, 'maxFErrFrac', 0.05, 'isTracking', 0, ...
        'useSequenceWalk', true, 'passesPerRun', [], 'walkWarmupTriggers', 1, ...
        'resyncEveryPass', true, 'newOrderEveryAverage', [], ...
        'plotIt', true, 'jsonLocation', ''), ...
        'calibAmplitudeDiff');

    if isempty(opts.amplitudes); opts.amplitudes = linspace(-40, -10, 10); end
    if isempty(opts.jsonLocation); opts.jsonLocation = setup7Paths().json; end

    % HIGH to LOW. Each point sizes its sweep from the previous Rabi frequency, and the
    % prediction only works downhill: a faster Rabi than expected still fits on a long
    % sweep, a slower one on a short sweep shows no oscillation at all.
    amps = sort(round(opts.amplitudes(:)'*100)/100, 'descend');

    out = struct('ok', false, 'opts', opts, 'amplitudes', amps);
    nRun = numel(amps)*(2*opts.nRepsPerPoint + double(opts.doCombined));

    fprintf('\n=== calibAmplitudeDiff ===\n');
    fprintf('%d amplitudes, %.1f to %.1f dBm. %d Rabi runs.\n', ...
            numel(amps), min(amps), max(amps), nRun);
    fprintf(['Every one re-uploads: on this AWG the amplitude is folded into the waveform\n' ...
             'samples, not set on the instrument, so there is no cache to hit. Budget about\n' ...
             '%.0f minutes.\n'], 1.4*nRun);
    fprintf('The field is not touched.\n');

    phi = opts.channelPhase;
    if isempty(phi)
        try
            j = jsondecode(fileread(opts.jsonLocation));
            if isfield(j, 'channelPhase') && numel(j.channelPhase) >= 2
                phi = j.channelPhase(2);
            end
        catch
        end
    end
    if isempty(phi); phi = 0; end
    out.channelPhase = phi;
    if opts.doCombined
        fprintf('Combined runs at channelPhase = %.2f deg (from the json).\n', phi);
    end

    % ---------------- the scan ----------------
    P = [];
    lastOmega = [];                 % drives the next point's sweep length
    printHeader();
    for i = 1:numel(amps)
        a = amps(i);
        p = struct('amplitude', a, 'omega1', NaN, 'omega1Err', NaN, ...
                   'omega2', NaN, 'omega2Err', NaN, 'omegaBoth', NaN, 'omegaBothErr', NaN, ...
                   'deltaVertical', NaN, 'R', NaN, 'Rerr', NaN, ...
                   'tauMax', NaN, 'why1', '', 'why2', '', 'whyBoth', '', 'good', false);

        tauList = sweepFor(lastOmega, opts);
        p.tauMax = max(tauList);

        [m1, m2] = deal([]);
        for r = 1:opts.nRepsPerPoint
            order = [1 2]; if mod(r,2) == 0; order = [2 1]; end
            for k = order
                amp = opts.offAmplitude*ones(1,2); amp(k) = a;
                m = runRabi(amp, 0, tauList, opts);
                if k == 1; m1 = [m1 m]; else; m2 = [m2 m]; end     %#ok<AGROW>
            end
        end
        [p.omega1, p.omega1Err, p.why1] = combine(m1);
        [p.omega2, p.omega2Err, p.why2] = combine(m2);

        if isfinite(p.omega1) && isfinite(p.omega2)
            p.deltaVertical = 20*log10(p.omega1/p.omega2);
            lastOmega = max(p.omega1, p.omega2);
            if opts.doCombined
                % The vertical delta is fine here -- see the header. At the optimum phase
                % Omega_both = Omega1 + Omega2 regardless of imbalance, so R does not care.
                mb = runRabi([a, a + p.deltaVertical], phi, tauList, opts);
                p.omegaBoth = mb.f; p.omegaBothErr = mb.fErr; p.whyBoth = mb.note;
                if mb.ok
                    S = p.omega1 + p.omega2;
                    p.R = mb.f/S;
                    p.Rerr = p.R*hypot(mb.fErr/mb.f, hypot(p.omega1Err, p.omega2Err)/S);
                end
            end
            p.good = true;
        end
        printRow(p);
        P = [P p];                                                 %#ok<AGROW>
    end
    out.points = P;

    good = logical([P.good]);
    fprintf('\n%d of %d amplitudes usable.\n', sum(good), numel(P));
    listRejects(P);
    if sum(good) < 4
        error(['calibAmplitudeDiff: only %d usable amplitudes -- not enough for a curve. ' ...
               'If the failures are all at the bottom of the range that is expected; ' ...
               're-run with amplitudes starting higher.'], sum(good));
    end

    % ---------------- the curves ----------------
    G = P(good);
    curve = struct();
    curve.amplitudes = [G.amplitude];
    curve.omega1     = [G.omega1];      curve.omega1Err = [G.omega1Err];
    curve.omega2     = [G.omega2];      curve.omega2Err = [G.omega2Err];
    curve.omegaBoth  = [G.omegaBoth];
    curve.R          = [G.R];           curve.Rerr      = [G.Rerr];
    curve.channelPhase = phi;
    curve.measuredAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');        %#ok<TNOW1,DATST>
    curve.jsonLocation = opts.jsonLocation;
    curve.opts = opts;
    % Sorted ascending in amplitude, because every consumer interpolates on it and an
    % interpolant built on a descending abscissa is a class of bug nobody enjoys.
    [~, ord] = sort(curve.amplitudes);
    for f = {'amplitudes','omega1','omega1Err','omega2','omega2Err','omegaBoth','R','Rerr'}
        curve.(f{1}) = curve.(f{1})(ord);
    end
    out.curve = curve;

    reportCurves(curve);
    out.linearityCeiling = reportLinearity(curve);
    curve.linearityCeiling = out.linearityCeiling;
    out.curve = curve;

    % ---------------- save ----------------
    if opts.save
        try
            P7 = setup7Paths();
            f = fullfile(P7.folder, 'amplitudeDiffCurve.mat');
            save(f, 'curve', '-mat');
            out.savedTo = f;
            fprintf('\nSaved to %s\n', f);
            fprintf('amplitudeDiffAt reads it from there, so nothing has to pass it around.\n');
        catch err
            fprintf('\nCould not save the curve (%s). It is in out.curve.\n', err.message);
        end
    end

    if opts.plotIt; drawIt(curve); end
    out.ok = true;

    fprintf('\n--- using it ---\n');
    fprintf('  d = amplitudeDiffAt(-16)          %% [0 delta] at -16 dBm\n');
    fprintf(['  expPropsWrite(struct(''amplitudeDiff'', amplitudeDiffAt(-16)))   %% writes it into\n' ...
             '                                    %% expProperties.json (backup + verify)\n']);
    fprintf('  [d, w, piNs] = amplitudeDiffAt(-16)   %% and the Rabi frequency to expect\n');
end

% ======================================================================

function tauList = sweepFor(lastOmega, opts)
% Sweep length from the previous point's Rabi frequency: about 3.5 periods, capped.
% Without this the low-amplitude points get a sweep far shorter than one period and fail
% for a reason that has nothing to do with the drive.
    if isempty(lastOmega) || ~isfinite(lastOmega) || lastOmega <= 0
        tauMax = 1.0;                                   % the standard sweep, for point one
    else
        % The next point is lower in amplitude, so slower. Allow generously for that rather
        % than predicting the drop: a sweep that is too long costs upload time, one that is
        % too short costs the point.
        tauMax = min(3.5/max(lastOmega/2, 0.05), opts.tauMaxCap);
    end
    tauMax = max(tauMax, 0.3);
    dtau = max(round(tauMax/opts.nTau, 3), 0.005);
    tauList = round(dtau : dtau : tauMax, 3);
end

function m = runRabi(amps, phaseDeg, tauList, opts)
    m = struct('ok', false, 'note', '', 'f', NaN, 'fErr', NaN, 'contrast', NaN, ...
               'r2', NaN, 'piNs', NaN, 'seconds', NaN);
    b = AWGExpRabi;
    b.loadJsonParams('all', opts.jsonLocation);
    b.nChannels = 2;                        % before amplitude, always
    b.amplitudeDiff = [0 0];                % explicit vector below, so no double-adding
    b.amplitude = amps;
    b.timeDelay = 0;
    b.channelPhase = [0 phaseDeg];
    b.repeats = opts.repeats;   b.averages = opts.averages;
    b.isTracking = opts.isTracking;
    b.balancedSequence = 1;
    b.detectionDuration = opts.rabiDetectionMHz;
    b.referenceDetectionDuration = b.detectionDuration;
    b.shouldAutosave = false;
    b.tau = tauList;
    b.useSequenceWalk = opts.useSequenceWalk;
    b.passesPerRun = opts.passesPerRun;
    b.walkWarmupTriggers = opts.walkWarmupTriggers;
    b.resyncEveryPass = opts.resyncEveryPass;
    if ~isempty(opts.newOrderEveryAverage)
        b.newOrderEveryAverage = logical(opts.newOrderEveryAverage);
    end

    txt = ''; tId = tic;
    try
        txt = evalc('b.run');
    catch err
        m.seconds = toc(tId); m.note = ['run failed: ' err.message];
        printIndented(txt); return
    end
    m.seconds = toc(tId);
    try
        o = rabiFitRobust(struct('exp', b, 'plotIt', false, 'verbose', false, ...
                                 'nFreq', opts.nFreq));
    catch err
        m.note = ['fit threw: ' err.message]; return
    end
    if ~isstruct(o) || ~isfield(o,'ok') || ~o.ok
        m.note = 'fit failed'; return
    end
    m.f = o.f; m.fErr = o.fErr; m.contrast = o.contrast; m.r2 = o.r2;
    m.piNs = 1e3*o.pulses(2);

    % A frequency below about one period per sweep is not measured, it is extrapolated from
    % a fraction of a cycle. Rejected explicitly, with the reason, because at the bottom of
    % an amplitude scan this is the EXPECTED outcome and should not look like a fault.
    span = max(tauList) - min(tauList);
    if m.f*span < 1
        m.note = sprintf(['too slow to fit at this drive: %.3f MHz is %.2f of a period over ' ...
                          'a %.1f us sweep'], m.f, m.f*span, span);
        return
    end
    if m.contrast < opts.minContrast
        m.note = sprintf('contrast %.2f %% below %.2f %%', 100*m.contrast, 100*opts.minContrast);
        return
    end
    if m.fErr/max(m.f,eps) > opts.maxFErrFrac
        m.note = sprintf('frequency only %.0f %% determined', 100*(1-m.fErr/max(m.f,eps)));
        return
    end
    m.ok = true; m.note = 'ok';
end

function [v, e, why] = combine(ms)
    v = NaN; e = NaN; why = '';
    if isempty(ms); why = 'no runs'; return; end
    ok = logical([ms.ok]);
    if ~any(ok); why = ms(1).note; return; end
    f = [ms(ok).f]; s = [ms(ok).fErr];
    s(~isfinite(s) | s <= 0) = median(s(isfinite(s) & s > 0));
    if all(~isfinite(s)); v = mean(f); e = std(f)/sqrt(numel(f)); return; end
    w = 1./s.^2;  v = sum(w.*f)/sum(w);  e = 1/sqrt(sum(w));
    if numel(f) > 1; e = max(e, std(f)/sqrt(numel(f))); end
end

% ======================================================================

function printHeader()
    fprintf('\n%8s %10s %10s %9s %11s %8s  %s\n', 'ampl/dBm', 'O1/MHz', 'O2/MHz', ...
            'gap/dB', 'Oboth/MHz', 'R', 'note');
end

function printRow(p)
    note = 'ok';
    if ~p.good
        note = 'REJECTED';
        if ~isempty(p.why1); note = ['ch1: ' p.why1]; end
        if ~isempty(p.why2); note = [note '  ch2: ' p.why2]; end
    elseif ~isempty(p.whyBoth) && ~strcmp(p.whyBoth, 'ok')
        note = ['both: ' p.whyBoth];
    end
    fprintf('%8.2f %10.4f %10.4f %9.3f %11.4f %8.4f  %s\n', ...
            p.amplitude, p.omega1, p.omega2, p.deltaVertical, p.omegaBoth, p.R, note);
end

function listRejects(P)
    bad = P(~logical([P.good]));
    if isempty(bad); return; end
    fprintf('\n--- rejected amplitudes ---\n');
    lowest = min([P.amplitude]);
    for i = 1:numel(bad)
        tag = '';
        if bad(i).amplitude <= lowest + 6
            tag = '   (bottom of the range -- expected, see the header)';
        end
        fprintf('  %+7.2f dBm: %s %s%s\n', bad(i).amplitude, bad(i).why1, bad(i).why2, tag);
    end
end

function reportCurves(c)
    fprintf('\n--- the two gain curves ---\n');
    fprintf('%10s %10s %10s %9s %9s\n', 'ampl/dBm', 'O1/MHz', 'O2/MHz', 'slope1', 'slope2');
    s1 = localSlope(c.amplitudes, 20*log10(c.omega1));
    s2 = localSlope(c.amplitudes, 20*log10(c.omega2));
    for i = 1:numel(c.amplitudes)
        fprintf('%10.2f %10.4f %10.4f %9.3f %9.3f\n', ...
                c.amplitudes(i), c.omega1(i), c.omega2(i), s1(i), s2(i));
    end
    fprintf(['\nSlope is dB of Rabi per dB of drive. 1.00 is a linear amplifier; below that\n' ...
             'is compression, and it is what makes amplitudeDiff amplitude-dependent -- if\n' ...
             'both slopes were 1.00 everywhere, one number would do.\n']);
    lo = min(min(s1), min(s2));
    if lo > 0.95
        fprintf(['\nBoth channels are linear to within 5 %% across this whole range, so one\n' ...
                 'amplitudeDiff really would have been enough. Worth knowing, and it means the\n' ...
                 'curve will be nearly flat.\n']);
    end
end

function ceilingA = reportLinearity(c)
% Where the pair stops adding. Reported as the highest amplitude at which R is still
% consistent with 1, which is the top of the usable operating range.
    ceilingA = NaN;
    if all(~isfinite(c.R)); fprintf('\nNo combined runs, so no linearity ceiling.\n'); return; end
    fprintf('\n--- R = Omega_both / (Omega1 + Omega2) ---\n');
    for i = 1:numel(c.amplitudes)
        if ~isfinite(c.R(i)); continue; end
        flag = '';
        if c.R(i) < 1 - max(3*c.Rerr(i), 0.02); flag = '  <- below 1, the pair is compressing'; end
        fprintf('%10.2f dBm   R = %.4f +/- %.4f%s\n', c.amplitudes(i), c.R(i), c.Rerr(i), flag);
    end
    okR = isfinite(c.R) & (c.R > 1 - max(3*c.Rerr, 0.02));
    if any(okR)
        ceilingA = max(c.amplitudes(okR));
        fprintf(['\nLinear up to %.2f dBm. Above that the two channels together deliver less\n' ...
                 'than the sum of their parts, which the single-channel curves cannot show --\n' ...
                 'each amplifier alone is still fine there.\n'], ceilingA);
    else
        fprintf('\nR is below 1 everywhere measured. Something shared is compressing at every drive.\n');
    end
end

function s = localSlope(a, G)
% dG/da by central differences, one-sided at the ends. Not a fit: a fitted slope would
% impose a shape, and the point of the table is to show the shape.
    s = nan(size(a));
    for i = 1:numel(a)
        lo = max(i-1,1); hi = min(i+1,numel(a));
        if hi > lo; s(i) = (G(hi)-G(lo))/(a(hi)-a(lo)); end
    end
end

function drawIt(c)
    C = struct('surface',[252 252 251]/255,'ink',[11 11 11]/255,'ink2',[82 81 78]/255, ...
               'data',[42 120 214]/255,'fit',[235 104 52]/255,'ink3',[140 138 133]/255);
    fig = figure('Name','calibAmplitudeDiff','NumberTitle','off','Color',C.surface);

    ax1 = subplot(3,1,1); hold(ax1,'on');
    plot(ax1, c.amplitudes, 20*log10(c.omega1), 'o-', 'LineWidth',1.6,'MarkerSize',5, ...
         'Color',C.data,'MarkerFaceColor',C.surface);
    plot(ax1, c.amplitudes, 20*log10(c.omega2), 's-', 'LineWidth',1.6,'MarkerSize',5, ...
         'Color',C.fit,'MarkerFaceColor',C.surface);
    % A unit-slope reference through the lowest point. Compression is a departure from a
    % straight line, and a straight line is very hard to judge by eye without one drawn.
    ref = 20*log10(c.omega1(1)) + (c.amplitudes - c.amplitudes(1));
    plot(ax1, c.amplitudes, ref, ':', 'LineWidth',1.2,'Color',C.ink3);
    legend(ax1, {'channel 1','channel 2','slope 1 reference'}, 'Location','northwest', ...
           'Box','off','TextColor',C.ink2);
    style(ax1,C); ylabel(ax1,'20 log_{10}\Omega  (dB)','Color',C.ink2);
    title(ax1,'gain curves','Color',C.ink,'FontWeight','normal');

    ax2 = subplot(3,1,2); hold(ax2,'on');
    d = nan(size(c.amplitudes));
    for i = 1:numel(c.amplitudes)
        try; d(i) = solveDelta(c, c.amplitudes(i)); catch; end
    end
    plot(ax2, c.amplitudes, 20*log10(c.omega1./c.omega2), 's--', 'LineWidth',1.2, ...
         'MarkerSize',5,'Color',C.ink3,'MarkerFaceColor',C.surface);
    plot(ax2, c.amplitudes, d, 'o-', 'LineWidth',1.8,'MarkerSize',5, ...
         'Color',C.data,'MarkerFaceColor',C.surface);
    legend(ax2, {'vertical gap (wrong)','solved offset (used)'}, 'Location','best', ...
           'Box','off','TextColor',C.ink2);
    style(ax2,C); ylabel(ax2,'amplitudeDiff  (dB)','Color',C.ink2);

    ax3 = subplot(3,1,3); hold(ax3,'on');
    if any(isfinite(c.R))
        hE = errorbar(ax3, c.amplitudes, c.R, c.Rerr);
        set(hE,'LineStyle','none','Color',C.data,'LineWidth',0.75);
        try; set(hE,'CapSize',4); catch; end
        plot(ax3, c.amplitudes, c.R, 'o','MarkerSize',5,'LineWidth',1, ...
             'Color',C.data,'MarkerFaceColor',C.surface);
        plot(ax3, [min(c.amplitudes) max(c.amplitudes)], [1 1], '-', ...
             'LineWidth',1.2,'Color',C.fit);
    end
    style(ax3,C); ylabel(ax3,'R','Color',C.ink2);
    xlabel(ax3,'channel 1 amplitude  (dBm)','Color',C.ink2);
    try; set(fig,'Position',[200 80 760 820]); catch; end
end

function style(ax,C)
    grid(ax,'on');
    set(ax,'Color',C.surface,'XColor',C.ink2,'YColor',C.ink2,'Box','off', ...
           'GridLineStyle','-','Layer','top','TickDir','out','FontSize',10);
    try; set(ax,'GridAlpha',0.18); catch; end
end

function d = solveDelta(c, a)
% Local copy of the solve, only so the plot can draw it. amplitudeDiffAt is the one callers
% should use -- it validates the range and reports why a solve failed.
    G1 = 20*log10(c.omega1); G2 = 20*log10(c.omega2);
    target = interp1(c.amplitudes, G1, a, 'pchip');
    if ~all(diff(G2) > 0); d = NaN; return; end
    if target < min(G2) || target > max(G2); d = NaN; return; end
    d = interp1(G2, c.amplitudes, target, 'pchip') - a;
end

function printIndented(txt)
    if isempty(txt); return; end
    L = strsplit(txt, newline);
    for i = 1:numel(L)
        if ~isempty(strtrim(L{i})); fprintf('      | %s\n', L{i}); end
    end
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
        if ~isfield(o, known{i}) || isempty(o.(known{i}))
            o.(known{i}) = d.(known{i});
        end
    end
end
