function d = coilCurrentDesign(opts)
%COILCURRENTDESIGN The coil currents at which calibCoilMatrix measures.
%
%   d = coilCurrentDesign(struct('M', M, 'b0', b0, 'Imin', [0 0 0], 'Imax', [5.2 5.2 5.2]))
%
% Space-filling (a Halton sequence) over the current box, kept only where the field the
% CURRENT model predicts is inside Brange. Why each choice:
%
%   SPREAD, NOT A STAR. The fit is linear in the currents, so leverage is the spread of the
%   currents, and a design that varies all three coils together also separates the
%   off-diagonal terms from the diagonal ones. Small steps around one point would leave
%   M poorly determined in every direction except the one stepped.
%
%   A FIELD SHELL. Below ~12 G the eight lines crowd inside two linewidths of D and the
%   per-spectrum field solve has little to work with. Above Bmax the outer lines leave the
%   sweep -- and esrFieldSolve's own warning applies: without the outer pair, |B| is weakly
%   constrained. The window is judged with the model in use (nominal or the last
%   calibration), so it is padded by the margin below.
%
%   QUASI-RANDOM ORDER. Consecutive points are far apart in current, so a slow drift (D
%   with temperature, the coils warming) does not line up with any one coil's current.
%
% OPTIONS
%   M, b0        the model the design is judged with (3x3 G/A, 1x3 G)
%   Imin, Imax   1x3 current limits, A
%   frac         [lo hi] fraction of each coil's range to use (default [0.05 1])
%   nPoints      (default 16)
%   Brange       [min max] |B|, Gauss (default [12 40])
%   skip         Halton points skipped at the start (default 17)
%
% OUTPUT
%   d.I          nPoints x 3, A
%   d.Bpred      nPoints x 3, the field the model predicts there
%   d.BmagPred   nPoints x 1
%   d.nCandidates, d.opts

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('M', eye(3), 'b0', [0 0 0], 'Imin', [0 0 0], ...
        'Imax', [1 1 1], 'frac', [0.05 1], 'nPoints', 16, 'Brange', [12 40], 'skip', 17));

    Imin = opts.Imin(:)'; Imax = opts.Imax(:)';
    lo = Imin + opts.frac(1) * (Imax - Imin);
    hi = Imin + opts.frac(2) * (Imax - Imin);
    nCand = 4000;
    U = [halton(nCand + opts.skip, 2), halton(nCand + opts.skip, 3), ...
         halton(nCand + opts.skip, 5)];
    U = U(opts.skip+1:end, :);
    Ic = lo + U .* (hi - lo);
    Bc = (opts.M * Ic' + opts.b0(:))';
    Bm = sqrt(sum(Bc.^2, 2));
    keep = find(Bm >= opts.Brange(1) & Bm <= opts.Brange(2));
    keep = keep(1:min(opts.nPoints, numel(keep)));

    d = struct();
    d.I = Ic(keep, :);
    d.Bpred = Bc(keep, :);
    d.BmagPred = Bm(keep);
    d.nCandidates = nCand;
    d.opts = opts;
    if numel(keep) < opts.nPoints
        warning('coilCurrentDesign:short', ['only %d of %d points have a predicted |B| ' ...
            'in [%g %g] G inside the current box. Widen Brange or frac.'], numel(keep), ...
            opts.nPoints, opts.Brange(1), opts.Brange(2));
    end
end

function x = halton(n, base)
% The first n points of the van der Corput sequence in this base.
    x = zeros(n, 1);
    for i = 1:n
        f = 1; r = 0; k = i;
        while k > 0
            f = f / base;
            r = r + f * mod(k, base);
            k = floor(k / base);
        end
        x(i) = r;
    end
end

function o = withDefaults(o, d)
    if isempty(o); o = struct(); end
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(o, fn{i}) || isempty(o.(fn{i})); o.(fn{i}) = d.(fn{i}); end
    end
end
