function res = fitPinholeFocus(S, x, y, opts)
%FITPINHOLEFOCUS  Find the stage position that puts a hole at the laser focus.
%
%   res = fitPinholeFocus(S, x, y)
%   res = fitPinholeFocus(S, x, y, opts)
%
%   A hole on piezo stages is scanned across the beam in one transverse axis
%   (x, or z) and along the beam axis (y); a fixed photodiode behind it reads
%   the transmitted power. The whole 2D scan is fitted at once with
%
%     S(y,x) = B + A * T( d ; a, w(y) ),  d = sqrt((x - xc(y))^2 + z(y)^2)
%
%     T(d; a, w)  fraction of a Gaussian beam (1/e^2 intensity radius w)
%                 that passes a round hole of radius a whose centre is a
%                 distance d from the beam axis. Exact in every regime:
%                 hole >> beam gives the erf-up / erf-down plateau whose edges
%                 are sharpest at the focus; hole << beam gives a Gaussian
%                 whose peak is highest and narrowest at the focus.
%     w(y)  = w0 * sqrt(1 + ((y - y0)/zR)^2)    Gaussian-beam caustic
%     xc(y) = xc + tilt*(y - y0)                beam not exactly along the
%                                               stage y axis
%     z(y)  = zOff + zTilt*(y - y0)             offset in the unscanned
%                                               axis; fitted only when the
%                                               data need it (opts.otherAxis)
%
%   The answer is res.best = [xc y0]: hole centre on the beam, at the focus.
%   Readings far off the model (> 10 robust sigma and > 30% of A, at most
%   1% of the points, e.g. readings taken before the stage arrived) are
%   ignored and listed in res.outliers.
%
%   INPUTS
%     S   photodiode readings. S(i,j) is at y(i), x(j): the layout of
%         imagesc(x, y, S) and of [X,Y] = meshgrid(x, y). A matrix in the
%         other layout ([numel(x) numel(y)], or ndgrid grids) is transposed
%         automatically. NaN marks missing points.
%     x   transverse stage positions (vector or meshgrid/ndgrid matrix)
%     y   beam-axis stage positions  (vector or meshgrid/ndgrid matrix)
%
%   OPTIONS (all optional)
%     opts.sigma         noise of S, scalar or matrix like S. [] (default):
%                        uniform, errors scaled by the residual scatter.
%     opts.rowAmplitude  false (default): one A for the whole scan.
%                        true: a separate A for every y row (laser drift).
%                        With false, the per-row amplitudes are checked
%                        anyway; if they drift by > 3% (and > 4 sigma) the
%                        fit is redone with true, a warning says so, and the
%                        single-amplitude result is kept in
%                        res.singleAmplitudeFit.
%     opts.fix           parameters to hold fixed, e.g. opts.fix.a = 12.5
%                        (known hole radius), opts.fix.tilt = 0.
%                        Fields: xc, tilt, a, w0, y0, zR, zOff, zTilt.
%     opts.init          starting values overriding the automatic ones
%                        (same fields).
%     opts.otherAxis     'auto' (default): also fit the offset and tilt of
%                        the unscanned axis, and keep that fit if it is
%                        significantly better (F > 20). true: always.
%                        false: never (z(y) = 0).
%     opts.otherAxisName name of the unscanned axis for messages (default
%                        'z' for an x scan, 'x' for a z scan).
%     opts.plot          true (default): data, model, residual, per-row
%                        widths and centres, profile at the focus.
%     opts.verbose       true (default): print the result.
%     opts.axisNames     {'x','y'} (default). Use {'z','y'} for a zy scan.
%     opts.units         'um' (default); only used in labels.
%
%   OUTPUT (struct)
%     res.best              [xc y0]
%     res.xc,  res.xcErr    hole centre (transverse) at the focus
%     res.y0,  res.y0Err    focus position along the beam
%     res.w0,  res.zR, res.a, res.tilt  (+ ...Err), res.A, res.B
%     res.otherAxisFitted   true if z(y) was fitted; then
%     res.zOff, res.zTilt   offset of the unscanned axis at y0 and its slope
%                           (the overall sign is not determined by this scan)
%     res.zCross            y where the beam crosses the hole centre in the
%                           unscanned axis (+ ...Err)
%     res.wAt(y), res.xcAt(y), res.zAt(y), res.halfChordAt(y)  function handles
%     res.model, res.residual, res.R2, res.chi2red, res.rmsResidual
%     res.rows              per-row check fits (a and B held at the global
%                           values, amplitude free per row):
%                           y, xc, xcErr, w, wErr, A, ok
%     res.outliers          [y x S] of the ignored readings
%     res.notes             findings (unscanned-axis offset, outliers)
%     res.warnings          cell array of plain-language warnings
%     res.singleAmplitudeFit  only if a laser drift was detected (see
%                           opts.rowAmplitude)
%     res.p, res.pErr, res.cov, res.names, res.free
%
%   NOTES
%   - An offset z0 in the unscanned axis makes this scan see a chord of
%     half-width sqrt(a^2 - z0^2) and wider edges. A constant offset does not
%     bias xc or y0, but if the beam is tilted in z relative to the stage y
%     travel, z0 changes along y and a simple fit mistakes the place where
%     the beam crosses the z centre for the focus. That is why z(y) is
%     fitted when the data need it. Do xy, then zy at the new y0.
%   - With a hole larger than the beam, the xy scan finds (mostly) the focus
%     of the beam's x width and the zy scan that of its z width. If the two
%     y0 differ by more than their errors the beam is astigmatic; the best
%     compromise is between them.
%   - Bidirectional raster scans with a lag show up as a shift between even
%     and odd rows; this is detected and reported.
%   - The photodiode must collect everything that passes the hole.
%   - No toolboxes needed (own Levenberg-Marquardt). Runs in Octave too.
%
%   Example
%     res = fitPinholeFocus(S, xPos, yPos, struct('axisNames', {{'z','y'}}));
%     moveTo(res.best)

if nargin < 4 || isempty(opts), opts = struct(); end
opts = setDefaults(opts);
warnState = quietSingular(); %#ok<NASGU>  restored when the function exits

% ------------------------------------------------------------------ input
[S, x, y, transposed] = arrangeInput(S, x, y);
[Ny, Nx] = size(S);
if Ny < 3 || Nx < 4
    error('fitPinholeFocus:size', 'Need at least 3 y rows and 4 x points (got %d x %d).', Ny, Nx);
end
M = isfinite(S);
if isempty(opts.sigma)
    W = ones(Ny, Nx);
else
    sg = double(opts.sigma);
    if isscalar(sg)
        W = ones(Ny, Nx) / sg;
    else
        if transposed, sg = sg.'; end
        if ~isequal(size(sg), [Ny Nx])
            error('fitPinholeFocus:sigma', 'opts.sigma must be a scalar or the same size as S.');
        end
        W = 1 ./ sg;
    end
    M = M & isfinite(W) & W > 0;
end
W(~M) = 0;
if nnz(M) < 20, error('fitPinholeFocus:data', 'Fewer than 20 valid points.'); end

Lx = x(end) - x(1);  Ly = y(end) - y(1);
dx = min(diff(x));   dy = min(diff(y));

P = buildP(S, M, W, x, y, opts.rowAmplitude);

% parameters: 1 xc, 2 tilt, 3 log a, 4 log w0, 5 y0, 6 log zR,
%             7 zOff (offset in the unscanned axis at y0), 8 zTilt
names = {'xc', 'tilt', 'a', 'w0', 'y0', 'zR', 'zOff', 'zTilt'};
isLog = logical([0 0 1 1 0 1 0 0]);
typ = [Lx, Lx/Ly, 1, 1, Ly, 1, Lx, Lx/Ly];
lb  = [x(1) - Lx, -5*Lx/Ly, log(dx/20), log(dx/50), y(1) - 3*Ly, log(dy/50), -20*Lx, -5*Lx/Ly];
ub  = [x(end) + Lx, 5*Lx/Ly, log(20*Lx), log(20*Lx), y(end) + 3*Ly, log(1000*Ly), 20*Lx, 5*Lx/Ly];

% --------------------------------------------- starting values (row fits)
v = sort(S(M));
B0 = v(max(1, round(0.02*numel(v))));
rowPk = zeros(Ny, 1);
for k = 1:Ny
    sk = S(k, M(k,:));
    if numel(sk) >= 3
        rowPk(k) = max(conv(sk - B0, [1 1 1]/3, 'same'));
    end
end
if max(rowPk) <= 0
    error('fitPinholeFocus:signal', 'No signal above the background was found.');
end
sel = find(rowPk > 0.2*max(rowPk));
if numel(sel) > 25, sel = sel(round(linspace(1, numel(sel), 25))); end

ri = struct('y', y(sel), 'xc', NaN(numel(sel),1), 'a', NaN(numel(sel),1), ...
            'w', NaN(numel(sel),1), 'pk', rowPk(sel));
for i = 1:numel(sel)
    k = sel(i);  m = M(k,:);
    xr = P.x(m);  sr = S(k,m);  wr = W(k,m);
    above = (sr - B0) > 0.5*rowPk(k);
    if ~any(above), continue; end
    xa = xr(above);
    fw = max(xa) - min(xa) + dx;
    starts = [fw/2 fw/8; fw/2.4 fw/3; fw/5 fw/1.1];
    bestSse = Inf;
    for s = 1:size(starts, 1)
        q0 = [(max(xa) + min(xa))/2, log(starts(s,1)), log(starts(s,2))];
        f = @(q) rowResid(q, xr, sr, wr, B0, NaN, 0, P);
        [q, sse] = lmFit(f, q0, true(1,3), lb([1 3 4]), ub([1 3 4]), typ([1 3 4]), 80);
        if sse < bestSse
            bestSse = sse;
            ri.xc(i) = q(1);  ri.a(i) = exp(q(2));  ri.w(i) = exp(q(3));
        end
    end
end
okr = isfinite(ri.w);
if ~any(okr), error('fitPinholeFocus:init', 'Could not fit any row; check the data.'); end

q0 = zeros(1, 6);
% beam axis: weighted line through the row centres
yr = ri.y(okr);  xcr = ri.xc(okr);  wt = ri.pk(okr);
if numel(yr) >= 2 && max(yr) > min(yr)
    c = wlinfit(yr, xcr, wt, 1);
    q0(2) = c(2);
else
    c = [xcr(1); 0];
end
% caustic: w^2 is quadratic in y
wOk = ri.w(okr);
wr2 = wOk.^2;
[~, kmin] = min(wOk);
fallback = [yr(kmin), log(wOk(kmin)), log(Ly/2)];
quadOK = false;
if numel(yr) >= 3
    ym = mean(y);
    cq = wlinfit((yr - ym)/Ly, wr2, wt, 2);
    if cq(3) > 0
        y0q = ym - Ly*cq(2)/(2*cq(3));
        w02 = cq(1) - cq(2)^2/(4*cq(3));
        if w02 > 0 && y0q > y(1) - Ly && y0q < y(end) + Ly
            quadOK = true;
            w0q = sqrt(w02);
            q0([5 4 6]) = [y0q, log(w0q), log(w0q*Ly/sqrt(cq(3)))];
        end
    end
end
if ~quadOK, q0([5 4 6]) = fallback; end
q0(1) = c(1) + c(2)*q0(5);
near = ri.w(okr) <= 1.5*min(ri.w(okr));
ar = ri.a(okr);
q0(3) = log(median(ar(near)));
q0(7:8) = 0;

% offset in the unscanned axis: it shrinks the chord seen in this scan,
% a_row^2 = a^2 - (zTilt*(y - yCross))^2, a downward parabola in y
zStart = [];
if numel(yr) >= 4
    cz = wlinfit((yr - mean(y))/Ly, ar.^2, wt, 2);
    if cz(3) < 0
        aZ2 = cz(1) - cz(2)^2/(4*cz(3));
        if aZ2 > 0
            yCr = mean(y) - Ly*cz(2)/(2*cz(3));
            tz = sqrt(-cz(3))/Ly;
            zStart = [sqrt(aZ2), tz, yCr];   % [a, |zTilt|, y of the crossing]
        end
    end
end

% user overrides and fixed parameters
free = true(1, 8);
for j = 1:8
    if isfield(opts.init, names{j}), q0(j) = toQ(opts.init.(names{j}), isLog(j)); end
    if isfield(opts.fix, names{j})
        q0(j) = toQ(opts.fix.(names{j}), isLog(j));
        free(j) = false;
    end
end
q0 = min(max(q0, lb), ub);
zFree = free(7:8);
if isequal(opts.otherAxis, false) || (ischar(opts.otherAxis) && strcmpi(opts.otherAxis, 'off'))
    zFree(:) = false;
end

% ------------------------------------------------- global fit, z held fixed
freeB = free;  freeB(7:8) = false;
qAlt = [];          % second start (narrowest row) if the first looks doubtful
if quadOK && ~isfield(opts.init, 'y0') && free(5)
    qAlt = q0;  qAlt([5 4 6]) = fallback;
    qAlt(1) = c(1) + c(2)*qAlt(5);
    qAlt(~free) = q0(~free);
    qAlt = min(max(qAlt, lb), ub);
end
[qBest, sseBest, nIt, startUsed] = basicFit(q0, P, freeB, lb, ub, typ, qAlt);

% gross outliers (e.g. readings taken before the stage arrived)
outl = zeros(0, 3);
[~, lin, F] = projResid(qBest, P);
Rw = P.W .* (S - F);  Rw(~P.M) = NaN;
rv = Rw(P.M);
sc = 1.4826*median(abs(rv - median(rv)));
bad = P.M & abs(Rw) > 10*sc & abs(S - F) > 0.3*median(abs(lin(2:end)));
if any(bad(:)) && nnz(bad) <= max(3, ceil(0.01*nnz(P.M)))
    [iy, ix] = find(bad);
    outl = [P.y(iy), P.x(ix).', S(bad)];
    M = P.M & ~bad;
    P = buildP(S, M, W, x, y, opts.rowAmplitude);
    [qBest, sseBest, nIt, startUsed] = basicFit(qBest, P, freeB, lb, ub, typ, []);
end

% ----------------------------- add the offset/tilt of the unscanned axis
useZ = false;  Fz = NaN;
if any(zFree)
    zs = {};
    qa = qBest;  aB = exp(qBest(3));
    if ~isempty(zStart) && zStart(1) >= 0.8*aB
        qa(3) = log(zStart(1));
        if free(8), qa(8) = zStart(2); end
        if free(7), qa(7) = zStart(2)*(qBest(5) - zStart(3)); end
        if abs(qa(7)) < 0.05*zStart(1), qa(7) = 0.05*zStart(1); end
        zs{end+1} = qa;
    end
    qa = qBest;
    if free(7), qa(7) = 0.3*aB; end
    qa(3) = log(sqrt(aB^2 + qa(7)^2));
    zs{end+1} = qa;
    sseZ = Inf;
    for i = 1:numel(zs)
        [qq, ss, ni] = lmFit(@(q) projResid(q, P), min(max(zs{i}, lb), ub), free, lb, ub, typ, 300);
        if ss < sseZ, qZ = qq;  sseZ = ss;  nItZ = ni; end
    end
    nParZ = nnz(free) + 2;
    Fz = ((sseBest - sseZ)/nnz(zFree)) / (sseZ/max(nnz(P.M) - nParZ, 1));
    if isequal(opts.otherAxis, true) || (Fz > 20 && sseZ < sseBest)
        useZ = true;
        qBest = qZ;  sseBest = sseZ;  nIt = nItZ;
    end
end
if ~useZ, free = freeB; end
q = qBest;
[~, lin, F, T] = projResid(q, P);
M = P.M;

% ------------------------------------------------------- uncertainties
nf = find(free);
nObs = nnz(M);
if ~P.rowAmp
    B = lin(1);  A = lin(2);
    J = zeros(nObs, numel(nf) + 2);
    for j = 1:numel(nf)
        h = 1e-5*typ(nf(j));
        qp = q;  qp(nf(j)) = qp(nf(j)) + h;
        qm = q;  qm(nf(j)) = qm(nf(j)) - h;
        Tp = modelT(qp, P);  Tm = modelT(qm, P);
        J(:,j) = -A * P.Wm .* (Tp(M) - Tm(M)) / (2*h);
    end
    J(:, end-1) = -P.Wm;
    J(:, end)   = -P.Wm .* T(M);
    nPar = numel(nf) + 2;
else
    B = lin(1);  A = lin(2:end);
    J = zeros(nObs, numel(nf));
    for j = 1:numel(nf)
        h = 1e-5*typ(nf(j));
        qp = q;  qp(nf(j)) = qp(nf(j)) + h;
        qm = q;  qm(nf(j)) = qm(nf(j)) - h;
        J(:,j) = (projResid(qp, P) - projResid(qm, P)) / (2*h);
    end
    nPar = numel(nf) + 1 + nnz(sum(abs(T) .* M, 2) > 0);
end
dof = max(nObs - nPar, 1);
chi2red = sseBest / dof;
if isempty(opts.sigma), s2 = chi2red; else, s2 = 1; end
C = covFromJ(J) * s2;

% to physical units
pq  = q;  pq(isLog) = exp(q(isLog));
G = ones(1, numel(nf));  G(isLog(nf)) = pq(nf(isLog(nf)));
Cq = C(1:numel(nf), 1:numel(nf)) .* (G.' * G);
pErr = zeros(1, 8);
pErr(nf) = sqrt(max(diag(Cq), 0)).';

res = struct();
res.otherAxisFitted = useZ;
res.zOff = pq(7);  res.zOffErr = pErr(7);
res.zTilt = pq(8);  res.zTiltErr = pErr(8);
if useZ && abs(pq(8)) > 0
    res.zCross = pq(5) - pq(7)/pq(8);
    gz = [-1/pq(8), pq(7)/pq(8)^2];          % d zCross / d(zOff, zTilt), plus d/dy0 = 1
    iz = [find(nf == 5), find(nf == 7), find(nf == 8)];
    if numel(iz) == 3
        gg = [1 gz];
        res.zCrossErr = sqrt(max(gg*Cq(iz, iz)*gg.', 0));
    else
        res.zCrossErr = NaN;
    end
else
    res.zCross = NaN;  res.zCrossErr = NaN;
end
res.otherAxisF = Fz;
res.outliers = outl;
res.xc = pq(1);  res.xcErr = pErr(1);
res.y0 = pq(5);  res.y0Err = pErr(5);
res.best = [res.xc res.y0];
res.w0 = pq(4);  res.w0Err = pErr(4);
res.zR = pq(6);  res.zRErr = pErr(6);
res.a  = pq(3);  res.aErr  = pErr(3);
res.tilt = pq(2);  res.tiltErr = pErr(2);
res.A = A;  res.B = B;
if ~P.rowAmp
    linErr = sqrt(max(diag(C(end-1:end, end-1:end)), 0));
    res.BErr = linErr(1);  res.AErr = linErr(2);
end
res.wAt  = @(yy) res.w0 * sqrt(1 + ((yy - res.y0)/res.zR).^2);
res.xcAt = @(yy) res.xc + res.tilt*(yy - res.y0);
res.zAt  = @(yy) res.zOff + res.zTilt*(yy - res.y0);
res.halfChordAt = @(yy) sqrt(max(res.a^2 - res.zAt(yy).^2, 0));
res.model = F;
res.residual = S - F;
res.residual(~M) = NaN;
res.R2 = 1 - sseBest / sum(P.Wm.^2 .* (S(M) - mean(S(M))).^2);
res.chi2red = chi2red;
res.rmsResidual = sqrt(mean(res.residual(M).^2));
res.names = names;  res.free = free;
res.p = pq;  res.pErr = pErr;
covFull = zeros(8);  covFull(nf, nf) = Cq;
res.cov = covFull;
res.x = P.x;  res.y = P.y;  res.data = S;
res.fitInfo = struct('iterations', nIt, 'start', startUsed, 'init', toPhys(q0, isLog), ...
                     'transposedInput', transposed);

% --------------------------------------------- per-row check fits
rows = struct('y', P.y, 'xc', NaN(Ny,1), 'xcErr', NaN(Ny,1), 'w', NaN(Ny,1), ...
              'wErr', NaN(Ny,1), 'A', NaN(Ny,1), 'AErr', NaN(Ny,1), 'ok', false(Ny,1));
Amed = median(abs(A));
for k = 1:Ny
    if rowPk(k) < 0.1*max(rowPk), continue; end
    m = M(k,:);
    if nnz(m) < 4, continue; end
    xr = P.x(m);  sr = S(k,m);  wr = W(k,m);
    qr0 = [res.xcAt(P.y(k)), log(res.wAt(P.y(k)))];
    zk = res.zAt(P.y(k));
    f = @(qq) rowResid(qq, xr, sr, wr, B, res.a, zk, P);
    [qr, sser] = lmFit(f, qr0, true(1,2), lb([1 4]), ub([1 4]), typ([1 4]), 60);
    Jr = zeros(numel(sr), 2);
    typR = typ([1 4]);
    for j = 1:2
        h = 1e-5*typR(j);
        qp = qr;  qp(j) = qp(j) + h;  qm = qr;  qm(j) = qm(j) - h;
        Jr(:,j) = (f(qp) - f(qm)) / (2*h);
    end
    dofr = max(numel(sr) - 3, 1);
    if isempty(opts.sigma), s2r = sser/dofr; else, s2r = 1; end
    Cr = covFromJ(Jr) * s2r;
    [~, Ar, ttr] = rowResid(qr, xr, sr, wr, B, res.a, zk, P);
    rows.xc(k) = qr(1);  rows.xcErr(k) = sqrt(max(Cr(1,1), 0));
    rows.w(k)  = exp(qr(2));  rows.wErr(k) = rows.w(k)*sqrt(max(Cr(2,2), 0));
    rows.A(k)  = Ar;  rows.AErr(k) = sqrt(s2r/ttr);
    rows.ok(k) = abs(Ar) > 0.05*Amed && isfinite(rows.wErr(k));
end
res.rows = rows;

% -------------------------------------------------------------- warnings
wn = {};
xn = opts.axisNames{1};  yn = opts.axisNames{2};  un = opts.units;
if res.y0 < y(1) || res.y0 > y(end)
    wn{end+1} = sprintf(['The focus (%s0 = %.4g) is outside the scanned %s range [%.4g, %.4g]: ' ...
        'it is an extrapolation. Rescan around it.'], yn, res.y0, yn, y(1), y(end));
end
if res.xc < x(1) || res.xc > x(end)
    wn{end+1} = sprintf('The hole centre (%s = %.4g) is outside the scanned %s range.', xn, res.xc, xn);
end
if free(5) && res.y0Err > Ly/4
    wn{end+1} = sprintf(['%s0 is poorly determined (+-%.3g %s). The beam width barely changes over ' ...
        'the scan: scan a longer %s range.'], yn, res.y0Err, un, yn);
end
yFar = max(abs([y(1) y(end)] - min(max(res.y0, y(1)), y(end))));
if res.wAt(res.y0 + yFar) < 1.1*res.w0
    wn{end+1} = sprintf(['The beam width changes by less than 10%% over the scanned %s range ' ...
        '(zR = %.3g %s): the focus position rests on a small effect.'], yn, res.zR, un);
end
yc = min(max(res.y0, y(1)), y(end));
Tend = modelT(q, setfield(P, 'y', yc)); %#ok<SFLD>
if max(Tend([1 end])) > 0.1
    wn{end+1} = sprintf(['At the focus the signal does not fall back to background at the ends of ' ...
        'the %s scan: background and amplitude are partly degenerate. Scan a wider %s range.'], xn, xn);
end
if res.a < res.w0
    wn{end+1} = sprintf(['The hole (radius %.3g) is smaller than the waist (w0 = %.3g %s): the signal ' ...
        'is Gaussian, not a plateau, so a and A are correlated. xc and %s0 are still fine.'], ...
        res.a, res.w0, un, yn);
end
if any(A < 0)
    wn{end+1} = 'The fitted amplitude is negative: is the photodiode signal inverted?';
end
hitB = free & (abs(q - lb) < 1e-9*max(typ, 1) | abs(q - ub) < 1e-9*max(typ, 1));
if any(hitB)
    wn{end+1} = sprintf('The fit hit a bound on: %s. Treat the result with care.', strjoin(names(hitB), ', '));
end
if ~isempty(opts.sigma) && chi2red > 3
    wn{end+1} = sprintf('chi2/dof = %.2g: the model does not describe the data to within opts.sigma.', chi2red);
end
% findings worth knowing that are not problems
notes = {};
zn = opts.otherAxisName;
if useZ
    tiltSig = abs(res.zTilt) > 3*res.zTiltErr;
    crossIn = tiltSig && isfinite(res.zCross) && res.zCross >= y(1) && res.zCross <= y(end);
    if ~tiltSig
        crossTxt = sprintf('and the same all along the scan (no significant %s tilt)', zn);
    elseif crossIn
        crossTxt = sprintf(['changing by %s %s per %s %s, so the beam crosses the hole''s %s centre ' ...
            'at %s = %s'], pm(abs(res.zTilt)*1000, res.zTiltErr*1000), un, '1000', un, zn, yn, ...
            pm(res.zCross, res.zCrossErr));
    else
        crossTxt = sprintf(['changing by %s %s per %s %s; it would cross the hole''s %s centre at ' ...
            '%s = %.4g, outside the scan'], pm(abs(res.zTilt)*1000, res.zTiltErr*1000), un, '1000', ...
            un, zn, yn, res.zCross);
    end
    notes{end+1} = sprintf(['The beam is also off the hole''s centre in %s: by %s %s at the focus, %s. ' ...
        'The fit includes this (an offset narrows the plateau and widens the %s edges, which a ' ...
        'simpler fit mistakes for defocus). The sign of the %s offset cannot be told from this ' ...
        'scan; the %s%s scan will show it.'], zn, pm(abs(res.zOff), res.zOffErr), un, crossTxt, ...
        xn, zn, zn, yn);
    if crossIn && abs(res.zCross - res.y0) < 2*hypot(res.zCrossErr, res.y0Err)
        wn{end+1} = sprintf(['The beam crosses the %s centre close to the focus, so the two effects ' ...
            'are hard to separate in this scan; check %s0 with the %s%s scan.'], zn, yn, zn, yn);
    end
end
if ~isempty(res.outliers)
    tx = sprintf('(%s,%s) = ', xn, yn);
    for i = 1:size(res.outliers, 1)
        tx = [tx sprintf('(%.4g, %.4g) ', res.outliers(i,2), res.outliers(i,1))]; %#ok<AGROW>
    end
    notes{end+1} = sprintf('Ignored %d reading(s) far off the model: %s', size(res.outliers, 1), tx);
end
res.notes = notes;

% alternate rows shifted (bidirectional raster with lag)
okr = rows.ok;
if nnz(okr) >= 6
    dxc = rows.xc - res.xcAt(rows.y);
    ev = okr & mod((1:Ny).', 2) == 0;  od = okr & mod((1:Ny).', 2) == 1;
    if nnz(ev) >= 3 && nnz(od) >= 3
        dd = mean(dxc(ev)) - mean(dxc(od));
        se = sqrt(var(dxc(ev))/nnz(ev) + var(dxc(od))/nnz(od));
        if abs(dd) > 4*se && abs(dd) > 0.2*dx
            wn{end+1} = sprintf(['Even and odd rows are shifted by %.3g %s in %s. A bidirectional scan ' ...
                'with lag? xc is the average of the two.'], dd, un, xn);
        end
    end
end
% laser drift: per-row amplitudes (A free per row) trending along y
if ~P.rowAmp && nnz(okr) >= 5
    ya = rows.y(okr);  Aa = rows.A(okr);  ea = rows.AErr(okr);
    X = [ones(size(ya)), (ya - mean(ya))/Ly];
    Wa = 1 ./ ea.^2;
    Hn = X.' * (X .* Wa);
    ca = Hn \ (X.' * (Aa .* Wa));
    chiA = sum(Wa .* (Aa - X*ca).^2) / max(numel(ya) - 2, 1);
    Ca = inv(Hn) * max(chiA, 1);
    dA = ca(2);  dAErr = sqrt(Ca(2,2));
    if abs(dA) > 0.03*abs(A) && abs(dA) > 4*dAErr
        oAlt = opts;
        oAlt.rowAmplitude = true;  oAlt.plot = false;  oAlt.verbose = false;
        if ~isempty(opts.sigma), oAlt.sigma = 1 ./ P.W; end   % already in this layout
        alt = fitPinholeFocus(S, P.x, P.y, oAlt);
        single = struct('xc', res.xc, 'xcErr', res.xcErr, 'y0', res.y0, 'y0Err', res.y0Err, ...
            'w0', res.w0, 'zR', res.zR, 'a', res.a, 'A', A, 'R2', res.R2);
        note = sprintf(['The amplitude changes by %.1f%% across the %s scan (laser drift?). A single ' ...
            'amplitude would bias the result (it gave %s = %s, %s = %s), so the result above uses ' ...
            'one amplitude per row (opts.rowAmplitude = true). The single-amplitude fit is in ' ...
            'res.singleAmplitudeFit.'], 100*dA/abs(A), yn, xn, pm(single.xc, single.xcErr), ...
            yn, pm(single.y0, single.y0Err));
        res = alt;
        res.singleAmplitudeFit = single;
        res.fitInfo.transposedInput = transposed;
        res.warnings = [{note}, res.warnings];
        P.rowAmp = true;
        wn = res.warnings;
    end
end
res.warnings = wn;

% ------------------------------------------------------------- report
if opts.verbose
    fprintf('\nfitPinholeFocus: hole at the focus at  %s = %s,  %s = %s %s\n', ...
        xn, pm(res.xc, res.xcErr), yn, pm(res.y0, res.y0Err), un);
    fprintf('  waist w0 = %s %s (1/e^2 radius), zR = %s %s\n', pm(res.w0, res.w0Err), un, pm(res.zR, res.zRErr), un);
    if res.otherAxisFitted
        aTxt = sprintf('(%s offset included)', opts.otherAxisName);
    else
        aTxt = '(as seen in this scan)';
    end
    fprintf('  hole radius a = %s %s %s, tilt d%s/d%s = %s\n', ...
        pm(res.a, res.aErr), un, aTxt, xn, yn, pm(res.tilt, res.tiltErr));
    if P.rowAmp
        fprintf('  background B = %.4g, amplitude per row (median %.4g)\n', res.B, median(res.A));
    else
        fprintf('  background B = %.4g, full-beam amplitude A = %.4g\n', res.B, res.A);
    end
    fprintf('  R^2 = %.5f, rms residual = %.3g (%.2g%% of A)', res.R2, res.rmsResidual, ...
        100*res.rmsResidual/abs(median(res.A)));
    if ~isempty(opts.sigma), fprintf(', chi2/dof = %.3g', res.chi2red); end
    fprintf('\n');
    for i = 1:numel(res.notes), fprintf('  NOTE: %s\n', res.notes{i}); end
    for i = 1:numel(wn), fprintf('  WARNING: %s\n', wn{i}); end
end

if opts.plot
    plotResult(res, P, opts);
end
end

% ======================================================================
%   model
% ======================================================================
function T = apertureT(D, a, s, xg, wg)
% Fraction of a 2D Gaussian (std s in each axis, s = w/2) inside a disk of
% radius a whose centre is at distance D from the beam axis. D is a matrix,
% s a scalar or a column matching the rows of D. The Gaussian is integrated
% analytically along the offset direction (erf over each chord) and with
% Gauss-Legendre in the perpendicular direction, v = a*sin(th), which keeps
% the integrand smooth in all regimes (accurate to ~1e-11 with 20 nodes).
thMax = asin(min(1, 9*s/a));
r2s = 1 ./ (sqrt(2)*s);
if isscalar(s) && isrow(D)
    % one row: all nodes at once (avoids loop overhead in the row fits)
    th = 0.5*thMax*(xg + 1);
    c  = a*cos(th);
    g  = (0.5*thMax*wg) .* exp(-0.5*(a*sin(th)/s).^2) / (sqrt(2*pi)*s) .* c;
    T  = g.' * (erf((c - D)*r2s) + erf((c + D)*r2s));
    return
end
T = zeros(size(D));
for k = 1:numel(xg)
    th = 0.5*thMax*(xg(k) + 1);
    c  = a*cos(th);
    g  = (0.5*thMax*wg(k)) .* exp(-0.5*(a*sin(th)./s).^2) ./ (sqrt(2*pi)*s) .* c;
    T  = T + g .* (erf((c - D).*r2s) + erf((c + D).*r2s));
end
end

function T = modelT(q, P)
w = exp(q(4)) * sqrt(1 + ((P.y - q(5))/exp(q(6))).^2);
D = P.x - (q(1) + q(2)*(P.y - q(5)));
zz = q(7) + q(8)*(P.y - q(5));          % offset in the unscanned axis
if any(zz ~= 0)
    D = sqrt(D.^2 + zz.^2);
else
    D = abs(D);
end
T = apertureT(D, exp(q(3)), w/2, P.xg, P.wg);
end

function P = buildP(S, M, W, x, y, rowAmp)
W(~M) = 0;
P.x  = x(:).';  P.y = y(:);
P.S  = S;  P.M = M;  P.W = W;
P.Sz = S;  P.Sz(~M) = 0;
P.Wm = W(M);  P.Ws = W(M) .* S(M);
P.rowAmp = logical(rowAmp);
[P.xg, P.wg] = gaussLegendre(20);
end

function [q, sse, nIt, startUsed] = basicFit(q0, P, free, lb, ub, typ, qAlt)
% Global fit; a second start (qAlt) if the first result looks doubtful.
f = @(qq) projResid(qq, P);
[q, sse, nIt] = lmFit(f, q0, free, lb, ub, typ, 300);
startUsed = 1;
sv = P.S(P.M);
r2 = 1 - sse / sum(P.Wm.^2 .* (sv - mean(sv)).^2);
if ~isempty(qAlt) && (r2 < 0.99 || q(5) < P.y(1) || q(5) > P.y(end))
    [q2, sse2, n2] = lmFit(f, qAlt, free, lb, ub, typ, 300);
    if sse2 < sse
        q = q2;  sse = sse2;  nIt = n2;  startUsed = 2;
    end
end
end

function [r, lin, F, T] = projResid(q, P)
% Residual with the linear parameters (B, A or B, A_k) solved exactly.
T = modelT(q, P);
if ~P.rowAmp
    X = [P.Wm, P.Wm .* T(P.M)];
    lin = solveSmall(X, P.Ws);
    F = lin(1) + lin(2)*T;
else
    qm = P.W .* T;  sm = P.W .* P.Sz;  em = P.W;
    qq = sum(qm.^2, 2);  qs = sum(qm.*sm, 2);  qe = sum(qm.*em, 2);
    ok = qq > 1e-12*max(qq);
    cs = zeros(size(qq));  ce = cs;
    cs(ok) = qs(ok)./qq(ok);  ce(ok) = qe(ok)./qq(ok);
    u = sm - qm.*cs;  v = em - qm.*ce;
    Bh = sum(u(:).*v(:)) / sum(v(:).^2);
    Ak = cs - Bh*ce;
    F = Bh + Ak.*T;
    lin = [Bh; Ak];
end
r = P.Wm .* (P.S(P.M) - F(P.M));
end

function [r, A, tt] = rowResid(q, xr, sr, wr, B, aFixed, zk, P)
% One row: q = [xc, log a, log w] (aFixed NaN) or [xc, log w] (a fixed);
% zk = offset of the row in the unscanned axis.
if isnan(aFixed)
    a = exp(q(2));  w = exp(q(3));
else
    a = aFixed;  w = exp(q(2));
end
D = abs(xr - q(1));
if zk ~= 0, D = sqrt(D.^2 + zk^2); end
T = apertureT(D, a, w/2, P.xg, P.wg);
t = wr .* T;  s = wr .* (sr - B);
tt = t*t.';
if tt > 0, A = (t*s.')/tt; else, A = 0; end
r = (s - A*t).';
end

% ======================================================================
%   numerics
% ======================================================================
function [q, sse, it] = lmFit(fun, q, free, lb, ub, typ, maxIter)
% Levenberg-Marquardt with Marquardt scaling, forward-difference Jacobian
% and box bounds (clamping).
nf = find(free);
r = fun(q);  sse = r.'*r;
it = 0;
if isempty(nf), return; end
lam = 1e-3;  nSmall = 0;
for it = 1:maxIter
    J = zeros(numel(r), numel(nf));
    for j = 1:numel(nf)
        h = 1e-6*typ(nf(j));
        qh = q;  qh(nf(j)) = qh(nf(j)) + h;
        J(:,j) = (fun(qh) - r) / h;
    end
    g = J.'*r;  H = J.'*J;
    dg = diag(H);  dg = max(dg, 1e-12*max(dg) + realmin);
    accepted = false;
    while lam <= 1e12
        dq = -(H + lam*diag(dg)) \ g;
        qn = q;  qn(nf) = min(max(q(nf) + dq.', lb(nf)), ub(nf));
        rn = fun(qn);  ssen = rn.'*rn;
        if isfinite(ssen) && ssen < sse
            accepted = true;  break;
        end
        lam = lam*10;
    end
    if ~accepted, break; end
    gain = (sse - ssen) / sse;
    step = max(abs(qn(nf) - q(nf)) ./ typ(nf));
    q = qn;  r = rn;  sse = ssen;
    lam = max(lam/10, 1e-10);
    if gain < 1e-10 || step < 1e-10
        nSmall = nSmall + 1;
        if nSmall >= 2, break; end
    else
        nSmall = 0;
    end
end
end

function C = covFromJ(J)
% inv(J'J) with column scaling, so that pinv does not drop directions just
% because the parameters have very different units.
cn = sqrt(sum(J.^2, 1));
cn(cn == 0) = 1;
Js = J ./ cn;
C = pinv(Js.' * Js) ./ (cn.' * cn);
end

function c = solveSmall(X, s)
H = X.'*X;
c = (H + 1e-14*trace(H)*eye(size(H))) \ (X.'*s);
end

function c = wlinfit(t, v, w, deg)
% weighted polynomial fit, coefficients from constant term up
t = t(:);  v = v(:);  w = sqrt(max(w(:), 0));
X = ones(numel(t), deg + 1);
for d = 1:deg, X(:, d+1) = t.^d; end
c = (X .* w) \ (v .* w);
end

function [xg, wg] = gaussLegendre(n)
k = 1:n-1;
b = k ./ sqrt(4*k.^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[xg, i] = sort(diag(D));
wg = 2*V(1, i).'.^2;
end

% ======================================================================
%   input, options, output helpers
% ======================================================================
function [S, x, y, transposed] = arrangeInput(S, x, y)
S = double(S);  x = double(x);  y = double(y);
transposed = false;
if ~isvector(x) || ~isvector(y)
    % grids: find which dimension each one varies along
    if ~isequal(size(x), size(S)) || ~isequal(size(y), size(S))
        error('fitPinholeFocus:grid', 'Grid inputs x and y must have the same size as S.');
    end
    if all(all(x == x(1,:))) && all(all(y == y(:,1)))      % meshgrid: rows = y
        x = x(1,:);  y = y(:,1);
    elseif all(all(x == x(:,1))) && all(all(y == y(1,:)))  % ndgrid: rows = x
        x = x(:,1);  y = y(1,:);  S = S.';  transposed = true;
    else
        error('fitPinholeFocus:grid', 'x and y grids are not regular (meshgrid or ndgrid).');
    end
else
    nx = numel(x);  ny = numel(y);
    if isequal(size(S), [ny nx])
        % expected layout
    elseif isequal(size(S), [nx ny])
        S = S.';  transposed = true;
    else
        error('fitPinholeFocus:size', 'size(S) = [%d %d] matches neither [numel(y) numel(x)] = [%d %d] nor its transpose.', ...
            size(S,1), size(S,2), ny, nx);
    end
end
x = x(:);  y = y(:);
[x, ix] = sort(x);  [y, iy] = sort(y);
S = S(iy, ix);
if any(diff(x) <= 0) || any(diff(y) <= 0)
    error('fitPinholeFocus:grid', 'x and y positions must be distinct.');
end
end

function opts = setDefaults(opts)
d = struct('sigma', [], 'rowAmplitude', false, 'fix', struct(), 'init', struct(), ...
           'plot', true, 'verbose', true, 'axisNames', {{'x', 'y'}}, 'units', 'um', ...
           'otherAxis', 'auto', 'otherAxisName', '');
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = d.(fn{i}); end
end
if ischar(opts.axisNames), opts.axisNames = {opts.axisNames(1), opts.axisNames(end)}; end
if isempty(opts.otherAxisName)
    if strcmpi(opts.axisNames{1}, 'x'), opts.otherAxisName = 'z';
    elseif strcmpi(opts.axisNames{1}, 'z'), opts.otherAxisName = 'x';
    else, opts.otherAxisName = 'other axis';
    end
end
if isnumeric(opts.otherAxis) || islogical(opts.otherAxis)
    opts.otherAxis = logical(opts.otherAxis);
end
end

function v = toQ(val, isLogPar)
if isLogPar, v = log(val); else, v = val; end
end

function p = toPhys(q, isLog)
p = q;  p(isLog) = exp(q(isLog));
end

function s = pm(v, e)
if e == 0
    s = sprintf('%.5g (fixed)', v);
elseif ~isfinite(e)
    s = sprintf('%.5g', v);
else
    nd = min(max(0, 1 - floor(log10(e))), 10);   % 2 significant digits on the error
    s = sprintf(sprintf('%%.%df +- %%.%df', nd, nd), v, e);
end
end

function c = quietSingular()
ids = {'MATLAB:nearlySingularMatrix', 'MATLAB:singularMatrix', 'MATLAB:rankDeficientMatrix', ...
       'Octave:nearly-singular-matrix', 'Octave:singular-matrix'};
st = {};
for i = 1:numel(ids)
    try   % MATLAB rejects the Octave IDs and vice versa
        s = warning('query', ids{i});
        warning('off', ids{i});
        st{end+1} = s; %#ok<AGROW>
    catch
    end
end
c = onCleanup(@() cellfun(@(s) warning(s.state, s.identifier), st));
end

% ======================================================================
%   plot
% ======================================================================
function plotResult(res, P, opts)
xn = opts.axisNames{1};  yn = opts.axisNames{2};  un = opts.units;
x = P.x;  y = P.y;
yy = linspace(y(1), y(end), 200).';
figure('Name', 'fitPinholeFocus', 'Color', 'w');

lim = [min(res.data(:)) max(res.data(:))];
if ~(lim(2) > lim(1)), lim = lim + [-1 1]; end

subplot(2, 3, 1);
imagesc(x, y, res.data, lim);  axis xy;  colorbar;  hold on;
hc = res.halfChordAt(yy);
plot(res.xcAt(yy), yy, 'w-', res.xcAt(yy) - hc, yy, 'w:', res.xcAt(yy) + hc, yy, 'w:');
if ~isempty(res.outliers)
    plot(res.outliers(:,2), res.outliers(:,1), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(res.xc, res.y0, 'w+', 'MarkerSize', 16, 'LineWidth', 2);
plot(res.xc, res.y0, 'ko', 'MarkerSize', 10, 'LineWidth', 1);
xlabel(sprintf('%s (%s)', xn, un));  ylabel(sprintf('%s (%s)', yn, un));
title(sprintf('Data: best %s = %.4g, %s = %.4g', xn, res.xc, yn, res.y0));

subplot(2, 3, 2);
imagesc(x, y, res.model, lim);  axis xy;  colorbar;
xlabel(sprintf('%s (%s)', xn, un));  ylabel(sprintf('%s (%s)', yn, un));
title('Global fit');

subplot(2, 3, 3);
rr = res.residual;  m = max(abs(rr(:)));  if ~(m > 0), m = 1; end
imagesc(x, y, rr, [-m m]);  axis xy;  colorbar;
xlabel(sprintf('%s (%s)', xn, un));  ylabel(sprintf('%s (%s)', yn, un));
title(sprintf('Residual (rms %.2g%% of A)', 100*res.rmsResidual/abs(median(res.A))));

ok = res.rows.ok;
subplot(2, 3, 4);
errorbar(res.rows.y(ok), res.rows.w(ok), res.rows.wErr(ok), 'o');  hold on;
plot(yy, res.wAt(yy), 'r-', 'LineWidth', 1.5);
yl = ylim;  plot([res.y0 res.y0], yl, 'k--');  ylim(yl);
xlabel(sprintf('%s (%s)', yn, un));  ylabel(sprintf('beam radius w (%s)', un));
title(sprintf('w_0 = %.3g, z_R = %.3g %s', res.w0, res.zR, un));
legend('rows (A free per row)', 'global fit', 'Location', 'best');

subplot(2, 3, 5);
errorbar(res.rows.y(ok), res.rows.xc(ok), res.rows.xcErr(ok), 'o');  hold on;
plot(yy, res.xcAt(yy), 'r-', 'LineWidth', 1.5);
yl = ylim;  plot([res.y0 res.y0], yl, 'k--');  ylim(yl);
xlabel(sprintf('%s (%s)', yn, un));  ylabel(sprintf('hole centre %s (%s)', xn, un));
title(sprintf('tilt d%s/d%s = %.3g', xn, yn, res.tilt));

subplot(2, 3, 6);
[~, k] = min(abs(y - res.y0));
xf = linspace(x(1), x(end), 400);
Pk = P;  Pk.x = xf;  Pk.y = y(k);
q = [res.xc res.tilt log(res.a) log(res.w0) res.y0 log(res.zR) res.zOff res.zTilt];
Tk = modelT(q, Pk);
if numel(res.A) > 1, Ak = res.A(k); else, Ak = res.A; end
plot(x, res.data(k,:), 'o', xf, res.B + Ak*Tk, 'r-', 'LineWidth', 1.5);
xlabel(sprintf('%s (%s)', xn, un));  ylabel('signal');
title(sprintf('Row nearest the focus (%s = %.4g)', yn, y(k)));
end
