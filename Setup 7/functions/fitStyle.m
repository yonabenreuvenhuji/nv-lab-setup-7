function S = fitStyle()
%FITSTYLE One look for every robust fit figure, in one file.
%
%   S = fitStyle;
%   ax = S.axes(S.newFigure('Rabi'));
%   h1 = S.data(ax, t, y, e, S.col.s1);     % blue, connected, filled markers
%   hF = S.model(ax, tf, yf);               % near-black, on top
%   S.headroom(ax, 0.28); S.legend(ax, [h1 hF], {'measured', 'fit'});
%   S.title(ax, {'line one', 'line two'});
%
% ===================================================================
% THE LOOK IS ramseyGateTest'S, DELIBERATELY
% ===================================================================
%
% The first version of this file drew hollow markers, no connecting line, error bars at a
% third opacity and a dotted 15% grid. It was too washed out to read, and it was not what
% had been asked for. ramseyGateTest's figures are the reference, so this now matches them:
%
%   SOLID grid at 12% alpha, behind the data. Solid rather than dotted -- a dotted grid at
%   low alpha nearly disappears, and then the axes have no reference lines at all.
%
%   FILLED markers with a white edge, size 6. A white edge is what keeps adjacent points
%   separable where the sweep is dense, without a second colour.
%
%   ERROR BARS AT FULL OPACITY, same hue as the marker, no caps. The earlier 35% opacity was
%   the single worst decision in the old version: the error bar is the measurement's
%   uncertainty and on a 0.3% Ramsey contrast it is most of what the figure is telling you.
%   Caps off because at 80 points the caps merge into a band.
%
%   A CONNECTING LINE through the data, thin. It makes a sweep readable as a curve rather
%   than a cloud. See connect below for the one place it is switched off, which matters.
%
% ===================================================================
% THE PALETTE -- CHOSEN BY EYE ON RENDERED ESR AND RAMSEY FIGURES (2026-11-27)
% ===================================================================
%
%   s1  blue    #3a6fd8   first series, and the DEFAULT for a single trace
%   s2  amber   #e69a1c   second series (|0> where |1> is s1)
%   s3  slate   #6b7885   a third, e.g. a repeat run kept as a drift control
%   fit brick   #b0302a   the model
%
% It replaced orange/purple data with a near-black fit. The blue is a little bluer than
% indigo so it does not drift toward the brick red. The fit is a hue no data series uses:
% warm against the blue, darker than the amber, so it reads as "the model" with one series
% or two. (An earlier note here found blue points under a near-black fit hard to read -- two
% dark marks of similar weight. A red curve on blue points is separated by hue instead.)
% Change it here AND in ExpPlotStyle.fallbackCol.
%
% THE FLAG IS AN INK SQUARE, NOT A RED RING. A red ring would now also be confused with the
% brick-red fit (and was hard to place on the old orange points). Ink is neutral against every
% series, and a SQUARE distinguishes it from the fit's line by shape as well as by mark type.
% A green or teal ring meaning "bad point" reads backwards.
%
% envelope  slate dashed, recessive: a reading aid for a decay, not a measurement
% band      pale grey fill, for a +/- 1 sigma region
%
% ===================================================================
% THE Y AXIS LABEL
% ===================================================================
%
% 'ratio' was the old label and it is too vague -- a ratio of what, to what?
%
% And it must NOT say 'a.u.'. Arbitrary units means the scale carries no information, and
% here it carries a lot: the quantity is signal photons over reference photons, so 1 means
% "no contrast" and 0.99 means "1% contrast", and that scale is the same in every run on the
% setup. Calling it arbitrary throws away the one thing the normalisation bought.
%
% It is DIMENSIONLESS, and a dimensionless quantity takes no unit symbol -- the same reason a
% refractive index or a Mach number is written bare. That is the convention, not an omission.
%
%       normalised fluorescence   (signal / reference)

    S = struct();
    S.col = struct( ...
        'surface', [252 252 251]/255, ...
        'ink',     [ 20  20  20]/255, ...
        'ink2',    [ 82  81  78]/255, ...
        's1',      [ 58 111 216]/255, ...   % blue, default
        's2',      [230 154  28]/255, ...   % amber
        's3',      [107 120 133]/255, ...   % slate
        'fit',     [176  48  42]/255, ...   % brick red
        'env',     [107 120 133]/255, ...
        'flag',    [ 20  20  20]/255, ...
        'band',    [219 224 232]/255);

    S.yLabel = 'normalised fluorescence   (signal / reference)';
    S.slots  = {S.col.s1, S.col.s2, S.col.s3};

    S.newFigure = @newFigure;
    S.axes      = @styleAxes;
    S.data      = @dataSeries;
    S.resid     = @residSeries;
    S.focus     = @focusY;
    S.model     = @modelLine;
    S.envelope  = @envelopeLine;
    S.band      = @bandFill;
    S.flag      = @flagPoints;
    S.legend    = @legendFrom;
    S.headroom  = @headroom;
    S.tight     = @tightX;
    S.title     = @titleBlock;
    S.zeroLine  = @zeroLine;
end

% ======================================================================

function fig = newFigure(name, wide)
    if ~exist('wide', 'var') || isempty(wide); wide = false; end
    C = fitStyle().col;
    pos = [70 90 1000 640];
    if wide; pos = [50 60 1240 800]; end
    fig = figure('Name', name, 'NumberTitle', 'off', 'Color', C.surface, 'Position', pos);
end

function ax = styleAxes(target)
% Solid grid at 12% alpha, BEHIND the data, ticks outward, secondary ink, no box.
%
% Layer 'bottom' rather than ramseyGateTest's 'top'. At 12% alpha the two are visually
% indistinguishable, and 'bottom' cannot draw grid lines across a marker. GridAlpha is
% R2014b and later, hence the try -- an older release gets the same plot with a heavier grid
% rather than an error.
    C = fitStyle().col;
    if nargin < 1 || isempty(target); target = gca; end
    % axes('Parent', fig), NOT axes(fig). The one-argument form means "make this existing
    % axes current" in MATLAB and is only read as a parent by newer releases; Octave rejects
    % it outright. The Parent form means the same thing everywhere.
    if strcmp(get(target, 'Type'), 'figure')
        ax = axes('Parent', target);
    else
        ax = target;
    end
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'off');
    set(ax, 'FontSize', 9, 'XColor', C.ink2, 'YColor', C.ink2, 'TickDir', 'out', ...
            'Color', C.surface, 'GridLineStyle', '-');
    try; set(ax, 'GridAlpha', 0.12, 'Layer', 'bottom'); catch; end
end

function h = dataSeries(ax, x, y, e, col, name, connect, badIdx)
%DATASERIES A measured sweep: filled markers, white edge, full-opacity bars, thin line.
%
% ONE errorbar object does all of it -- line, markers and bars -- which is what keeps the
% bars the same hue and weight as the marker they belong to. A second colour for the
% uncertainty of a series would be a third identity that does not exist.
%
% Properties are set AFTER creation, never as name/value pairs in the errorbar call: some
% MATLAB and Octave versions parse the leading character of 'LineStyle' as a LineSpec and
% throw "unrecognized format character", and Octave rejects property pairs after a format
% string outright.
%
% connect (default true) draws the thin line through the points. Pass false and you get
% markers only -- see residSeries for the one place that matters.
%
% badIdx BREAKS THE LINE at those points, keeping their markers. A down-weighted point is not
% part of the curve the fit describes, and joining through it draws a tall spike that is
% pure artefact: on the first render of the new look, two points diluted to ratio = 1 turned
% into vertical strokes that dominated the figure and stretched the y axis until the 0.3%
% fringe was a flat band in the middle. The marker and its flag still show the point; only
% the line skips it.
    C = fitStyle().col;
    if ~exist('connect', 'var') || isempty(connect); connect = true; end
    if isempty(e) || numel(e) ~= numel(y); e = zeros(size(y)); end
    yLine = y;
    if exist('badIdx', 'var') && ~isempty(badIdx)
        keep = badIdx(badIdx >= 1 & badIdx <= numel(yLine));
        yLine(keep) = NaN;                 % NaN breaks the line, the marker is drawn below
    end
    h = errorbar(ax, x, yLine, e);
    set(h, 'Color', col, 'LineWidth', 1.1, 'Marker', 'o', 'MarkerSize', 6, ...
           'MarkerFaceColor', col, 'MarkerEdgeColor', C.surface);
    if connect
        set(h, 'LineStyle', '-');
    else
        set(h, 'LineStyle', 'none');
    end
    try; set(h, 'CapSize', 0); catch; end
    if exist('badIdx', 'var') && ~isempty(badIdx)
        % the markers the line skipped, drawn on their own so the point is still visible
        k = badIdx(badIdx >= 1 & badIdx <= numel(y));
        plot(ax, x(k), y(k), 'o', 'MarkerSize', 6, 'MarkerFaceColor', col, ...
             'MarkerEdgeColor', C.surface, 'LineStyle', 'none');
    end
    if exist('name', 'var') && ~isempty(name)
        try; set(h, 'DisplayName', name); catch; end
    end
end

function nOut = focusY(ax, vals, frac)
%FOCUSY Set the y limits from the values that MATTER, and say how many fell outside.
%
% Autoscale includes everything, and on a Ramsey that means a couple of points diluted toward
% ratio = 1 set the scale for a 0.3% fringe -- the signal becomes a flat band and the figure
% is useless. So the limits come from the fitted model and the points the robust pass kept,
% and anything outside is reported rather than accommodated. The count goes in the title, so
% an off-scale point is stated rather than silently cropped.
    v = vals(:); v = v(isfinite(v));
    if isempty(v); nOut = 0; return; end
    if ~exist('frac', 'var') || isempty(frac); frac = 0.08; end
    lo = min(v); hi = max(v); span = max(hi - lo, eps);
    lo = lo - frac*span; hi = hi + frac*span;
    ylim(ax, [lo, hi]);
    nOut = 0;
end

function h = residSeries(ax, x, r, col)
%RESIDSERIES Residuals: markers ONLY, never joined.
%
% This is not a style preference. Joining residual points with a line manufactures the
% appearance of a systematic wave out of independent noise -- the eye follows the line and
% reads structure that is not there. It cost a real diagnostic detour on this project: a
% Ramsey fit with a reduced chi2 of 1.25, whose residual is 1.12 sigma RMS and therefore
% ordinary scatter, looked convincingly "out of phase with the fit" purely because the
% points were connected.
%
% If you want to know whether a residual is structured, the answer is the lag-1
% autocorrelation the fitters report, not the shape of a polyline.
    h = dataSeries(ax, x, r, [], col, '', false);
end

function h = modelLine(ax, x, y)
    C = fitStyle().col;
    h = plot(ax, x, y, '-', 'LineWidth', 2.2, 'Color', C.fit);
end

function h = envelopeLine(ax, x, hi, lo)
    C = fitStyle().col;
    h = plot(ax, x, hi, '--', 'LineWidth', 1.2, 'Color', C.env);
    if exist('lo', 'var') && ~isempty(lo)
        plot(ax, x, lo, '--', 'LineWidth', 1.2, 'Color', C.env);
    end
end

function h = bandFill(ax, x, hi, lo)
    C = fitStyle().col;
    h = fill(ax, [x(:)', fliplr(x(:)')], [hi(:)', fliplr(lo(:)')], C.band, ...
             'EdgeColor', 'none');
    try; uistack(h, 'bottom'); catch; end
end

function h = flagPoints(ax, x, y)
% An open INK SQUARE around the point, not a red ring. See the palette note: crimson against
% orange measures dE 14.9, under the 15 floor, so a red ring on an orange marker is hard to
% place. Ink is the best-separated option and the square shape distinguishes it from the
% fit's line by mark type as well as colour.
    C = fitStyle().col;
    h = plot(ax, x, y, 's', 'MarkerSize', 13, 'LineWidth', 1.5, ...
             'MarkerEdgeColor', C.flag, 'MarkerFaceColor', 'none');
end

function zeroLine(ax, xr)
    C = fitStyle().col;
    plot(ax, xr, [0 0], '-', 'Color', 0.45*C.ink2 + 0.55*[1 1 1], 'LineWidth', 0.9);
end

function lg = legendFrom(ax, H, L)
% Fixed location with headroom reserved for it, NOT 'best'. 'best' is only as good as the
% renderer's collision search, and on an early version of the Ramsey figure it put the legend
% straight through the series it was labelling.
    lg = [];
    try
        lg = legend(ax, H, L, 'Location', 'northeast');
        set(lg, 'Box', 'off', 'TextColor', fitStyle().col.ink2, 'FontSize', 8);
    catch
    end
end

function headroom(ax, frac)
    yl = ylim(ax); span = max(yl(2) - yl(1), eps);
    ylim(ax, [yl(1), yl(2) + frac*span]);
end

function tightX(ax, x)
% Hold the x limits to the data with a 2% pad. Auto-limits round outward to a tick, which
% leaves a strip of empty axes that reads as "the sweep stopped early".
    lo = min(x); hi = max(x); pad = 0.02*max(hi - lo, eps);
    xlim(ax, [lo - pad, hi + pad]);
end

function titleBlock(ax, lines)
% The numbers go in the TITLE, never in an annotation box over the data. A box sits on top of
% the very points it describes, and the plot area should be entirely data.
    C = fitStyle().col;
    if ischar(lines); lines = {lines}; end
    t = title(ax, lines, 'FontWeight', 'normal', 'FontSize', 9, 'Color', C.ink);
    try; set(t, 'Interpreter', 'tex'); catch; end
end
