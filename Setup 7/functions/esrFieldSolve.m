function out = esrFieldSolve(centres, centreErr, opts)
%ESRFIELDSOLVE The field vector from ESR line positions, by fitting the exact Hamiltonian.
%
%   out = esrFieldSolve(centres)                       % MHz
%   out = esrFieldSolve(centres, centreErr, struct('D0', 2870))
%
% Four parameters -- Bx, By, Bz in the crystal frame and D -- predict ALL EIGHT line
% positions exactly. Fit them to the measured centres and everything else follows: |B|,
% the tilt from the nearest axis, the four projections, and a prediction for where every
% line sits at any other field setting.
%
% ===================================================================
% WHY THIS REPLACES PAIRING-AND-ARITHMETIC
% ===================================================================
%
% The previous analysis paired lines by mirror symmetry about a single D and then took
% splittings. That is exact ONLY at first order, and it fails in a specific, measurable
% way that showed up the moment the analysis was run on a synthetic spectrum with a known
% field:
%
%   at 45 G and 7 degrees of tilt, the ALIGNED orientation has almost no transverse
%   field, so its pair straddles D = 2870.6. The three off-axis orientations have 10 to 18
%   G along the axis and 41 to 44 G across it, so their pairs are pushed UP by 7.1 to 7.9
%   MHz by the second-order transverse shift. Their means are 2877.8, 2877.3, 2878.4.
%
% There is no single D that pairs all four orientations. Estimating one from the data
% gives 2878, because three of the four pairs agree on it -- and the aligned pair, whose
% mean is the physical D and whose splitting is the most reliable number in the whole
% spectrum, is then declared to have no mirror partner AND DISCARDED. That is exactly
% what happened: the 2745.45/2995.94 pair, correct to 0.1 MHz, was thrown out as "not an
% NV line" while the shifted pairs were kept.
%
% Fitting the Hamiltonian has no such failure mode: the shift is part of the model, so
% each orientation is allowed its own pair mean and the physics supplies the relationship
% between them instead of the analysis having to assume they are equal.
%
% ===================================================================
% HOW THE FIT IS MADE ROBUST
% ===================================================================
%
% ASSIGNMENT, not ordering. Lines get missed -- blended, too shallow, outside the sweep --
% so the k-th measured line is not the k-th predicted one. Each MEASURED line is matched
% to its nearest PREDICTED line, so a missing line weakens the fit rather than corrupting it.
%
% A PREDICTED LINE INSIDE THE SWEEP WITH NOTHING THERE COSTS SOMETHING (missWeight). Letting
% unmatched predictions cost nothing was the first version, and it has a specific failure:
% if the sweep does not contain the outer pair, HALF the field fits the remaining lines
% about as well as the true field does, because at half |B| both pairs move inside the
% window. Tested by clipping a 45 G synthetic to 2790-2950: the fit reported 23.06 G and a
% 41 degree tilt, with no complaint. The true field predicts 6 lines in that window and
% all 6 are seen; the 23 G solution predicts 8 and only 6 are seen, so charging for the two
% absent predictions separates them. The charge is BOUNDED and modest, so a genuinely weak
% or blended line that the fitter missed costs a little and does not overturn a good fit.
%
% A BOUNDED COST. The residual for each line is min(d^2/sigma^2, huber^2), so a feature
% that is not an NV line at all -- the resonator dip, which is what the discarded mirror
% test was really trying to catch -- contributes a fixed penalty instead of dragging the
% field to fit it. Such lines come back in out.unmatched, named and quantified, which is
% strictly more useful than silently dropping them.
%
% MULTI-START. The cost has local minima: a tilt in the wrong azimuth can match three of
% four orientations. So the fit starts from a grid of tilts and azimuths, with |B| seeded
% from the largest observed splitting, and keeps the best. 4 parameters and a 3x3
% eigenproblem make this cheap -- a few hundred Nelder-Mead runs, well under a second.
%
% NELDER-MEAD, from base MATLAB. fminsearch, not lsqnonlin: no Optimization Toolbox on
% this setup, and the assignment step makes the cost piecewise anyway, so a
% derivative-free method is the honest choice.
%
% WHAT IT CANNOT DO, MEASURED
%
% Below about 2.5 degrees of tilt the three off-axis orientations are no longer resolved --
% their splittings differ by less than a linewidth -- so line POSITIONS carry almost no
% information about the tilt and this fit cannot recover it. Graded against synthetic
% spectra at 45 G with known fields: 10 deg -> 9.98, 7 -> 6.96, 4 -> 3.63, 2 -> 1.95,
% 1 -> 0.10, 0 -> 2.75. Note the last two: near alignment the reported tilt is not merely
% imprecise, it can be wrong in EITHER direction. |B| stays good throughout (45.03 to 45.06
% G in every one of those cases), and so does D.
%
% So use this for |B|, for D, for the line assignment, and for the tilt when the tilt is
% more than a couple of degrees. Do NOT use it as the endgame alignment metric: the off-axis
% group's WIDTH still responds where its positions do not, which is what
% calibFieldAlign's 'innerspan' metric measures.
%
% OPTIONS
%   D0, dTol    zero-field splitting start and how far it may move (default 2870, 12 MHz)
%   huberMHz    a line further than this from any prediction is treated as an outlier
%               rather than pulling the field (default 6)
%   sweepRange  [lo hi] MHz. Predicted lines outside it are known to be unobservable and
%               are not held against the fit (default: from the measured lines, padded)
%   nAz         azimuth start points (default 12)
%   tiltStarts  degrees (default [0 2 5 8 12 18 25])
%   missWeight  cost charged for a predicted line inside the sweep with nothing there,
%               in units of one huber residual (default 0.5). This is what stops half the
%               true |B| fitting a sweep that clips the outer pair -- see above.
%   verbose     print the report (default true)
%
% OUTPUT
%   ok, why
%   B, Bvec, tiltDeg, D          |B| in Gauss, the vector in the crystal frame, the angle
%                                to the NEAREST <111> axis, and the fitted D
%   projections                  signed B.a_i, Gauss, sorted by |.| descending
%   splittings, means            predicted, per orientation, MHz
%   predicted                    all 8 predicted line positions, MHz
%   matched                      struct array: measured centre, predicted, axis, residual
%   unmatched                    measured lines further than huberMHz from any prediction
%   rms, chi2                    residual of the matched lines
%   outerVisible                 false when the aligned pair falls outside the sweep, in
%                                which case |B| is weakly constrained -- see below
%   predictLines                 function handle: f = out.predictLines(Bscale) gives the
%                                line positions with |B| scaled, direction held
%   fieldForLineAt               function handle: [scale, which] = fieldForLineAt(target)
%                                solves for the |B| scaling that puts a line on target

    if ~exist('centreErr', 'var') || isempty(centreErr)
        centreErr = ones(size(centres));
    end
    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('D0', 2870, 'dTol', 12, 'huberMHz', 6, ...
        'sweepRange', [], 'nAz', 12, 'tiltStarts', [0 2 5 8 12 18 25], 'verbose', true, ...
        'missWeight', 0.5));

    GAMMA = 2.8025;
    c = sort(centres(:))';
    ce = centreErr(:)';
    ce(~isfinite(ce) | ce <= 0) = 1;
    if numel(ce) ~= numel(c); ce = ones(size(c)); end

    out = struct('ok', false, 'why', '', 'B', NaN, 'tiltDeg', NaN, 'D', NaN);
    if numel(c) < 3
        out.why = sprintf('%d lines is not enough to determine a field vector', numel(c));
        return
    end

    if isempty(opts.sweepRange)
        pad = 0.1*(max(c) - min(c)) + 10;
        rng = [min(c) - pad, max(c) + pad];
    else
        rng = opts.sweepRange;
    end

    Ax = [ 1  1  1; 1 -1 -1; -1  1 -1; -1 -1  1] / sqrt(3);

    % |B| start: the largest splitting belongs to the best-aligned orientation, so
    % max(splitting)/(2*gamma) is a lower bound on |B| and usually within a percent.
    spanB = (max(c) - min(c)) / (2*GAMMA);

    % COARSE SCAN FIRST, then optimise only the promising starts. Running fminsearch from
    % every one of ~90 start points takes 11 s, and calibFieldAlign calls this once per
    % scan point -- 28 times in a two-round alignment. Evaluating the cost is a handful of
    % 3x3 eigenproblems, so scanning is nearly free and only the few starts that already
    % look right are worth the several hundred function evaluations of a full descent.
    fun = @(p) costOf(p, c, ce, Ax, GAMMA, opts, rng);
    starts = []; costs = [];
    a1 = Ax(1,:)';
    u = [0;0;1] - (a1'*[0;0;1])*a1;
    if norm(u) < 1e-9; u = [0;1;0] - (a1'*[0;1;0])*a1; end
    u = u/norm(u); w = cross(a1, u);
    for tiltDeg = opts.tiltStarts
        for iaz = 1:opts.nAz
            az = 2*pi*(iaz-1)/opts.nAz;
            th = tiltDeg*pi/180;
            nHat = cos(th)*a1 + sin(th)*(cos(az)*u + sin(az)*w);
            % |B| STARTS MUST GO WELL ABOVE THE OBSERVED SPAN. max(splitting)/(2*gamma) is a
            % lower bound on |B| ONLY IF the aligned pair is in the sweep. If it is not --
            % the sweep was too narrow, or the outer lines are too shallow to fit -- the
            % widest visible splitting belongs to an off-axis orientation and is between 1/3
            % and 1/2 of the aligned one, so the true |B| is 2x to 3x this seed and was never
            % in the start set. That is exactly how a clipped 45 G sweep came back as 20 G
            % with a 31 degree tilt: the right answer wins the cost function easily, it was
            % simply never tried.
            for bs = [0.9 1.0 1.1 1.5 2.0 3.0]
                p0 = [bs*spanB*nHat(:)', opts.D0];
                starts(end+1, :) = p0;                                  %#ok<AGROW>
                costs(end+1) = fun(p0);                                 %#ok<AGROW>
            end
        end
    end
    [~, ord] = sort(costs);
    nTry = min(12, numel(costs));
    best = []; bestCost = inf;
    o = optimset('Display', 'off', 'MaxFunEvals', 4000, 'MaxIter', 4000, ...
                 'TolX', 1e-4, 'TolFun', 1e-6);
    for k = ord(1:nTry)
        try
            [p, fv] = fminsearch(fun, starts(k, :), o);
        catch
            continue
        end
        if fv < bestCost; bestCost = fv; best = p; end
    end
    if isempty(best)
        out.why = 'no start point converged'; return
    end

    % A second polish from the winner, in case the best start was still on a slope.
    o = optimset('Display', 'off', 'MaxFunEvals', 8000, 'MaxIter', 8000);
    best = fminsearch(@(p) costOf(p, c, ce, Ax, GAMMA, opts, rng), best, o);

    Bvec = best(1:3)'; D = best(4);
    [pred, spl, mns, proj] = predictAll(Bvec, D, Ax, GAMMA);

    % ---- assignment and residuals -------------------------------------------
    M = struct('centre', {}, 'predicted', {}, 'axis', {}, 'resid', {});
    unmatched = []; matchedIdx = [];
    for k = 1:numel(c)
        [d, j] = min(abs(pred - c(k)));
        if d <= opts.huberMHz
            M(end+1) = struct('centre', c(k), 'predicted', pred(j), ...
                              'axis', ceil(j/2), 'resid', c(k) - pred(j));   %#ok<AGROW>
            matchedIdx(end+1) = k;                                            %#ok<AGROW>
        else
            unmatched(end+1) = c(k);                                          %#ok<AGROW>
        end
    end
    out.matched = M; out.unmatched = unmatched;
    inRange = pred(pred >= rng(1) & pred <= rng(2));
    miss = [];
    for j = 1:numel(inRange)
        if min(abs(c - inRange(j))) > opts.huberMHz
            miss(end+1) = inRange(j);                                         %#ok<AGROW>
        end
    end
    out.missing = miss;
    if numel(M) < 3
        out.why = sprintf(['only %d of %d lines could be matched to any predicted ' ...
            'transition -- the spectrum is not an NV pattern, or D is far off'], numel(M), numel(c));
        return
    end
    res = [M.resid];
    out.rms  = sqrt(mean(res.^2));
    % Index-based, NOT ismember on the centres: two lines can share a value after a refit
    % and ismember would then return a mask of the wrong length, which is an error thrown
    % from a line that has nothing to do with the cause.
    out.chi2 = mean((res ./ ce(matchedIdx)).^2);

    [~, iBig] = max(abs(proj));
    out.B = norm(Bvec);
    out.Bvec = Bvec;
    out.D = D;
    out.tiltDeg = acosd(min(abs(proj(iBig))/max(norm(Bvec), eps), 1));
    [~, ord] = sort(abs(proj), 'descend');
    out.projections = proj(ord);
    out.splittings = spl(ord);
    out.means = mns(ord);
    out.predicted = pred;
    out.alignedAxis = iBig;
    out.ok = true;

    % IS THE ALIGNED PAIR ACTUALLY IN THE SWEEP? If it is not, |B| is barely constrained,
    % and this has to be said out loud rather than left for the caller to discover.
    %
    % The four orientations are not equal witnesses. The aligned one has the largest
    % projection, so its splitting is the most direct measure of |B| there is; the other
    % three have projections between B/3 and B/2 and, crucially, their line positions trade
    % off against the field DIRECTION -- a larger |B| tilted differently reproduces them
    % almost exactly. Measured on a 45 G synthetic clipped to 2790-2950, which keeps the six
    % off-axis lines and loses the aligned pair: fitting the exact line positions recovers
    % 45.00 G and 7.00 deg, but with the ordinary sub-MHz scatter of a real fit on those same
    % six lines, a 55.3 G solution fits to 0.03 MHz rms. Ten Gauss for three hundredths of a
    % MHz. The fit is not wrong, the data does not contain the answer.
    %
    % So this is not a defect to be corrected by a better optimiser or a tighter penalty. The
    % only cure is a sweep wide enough to hold the outer pair, and the only honest thing the
    % code can do is refuse to present |B| as measured when it is not.
    alignedLines = pred(2*iBig-1 : 2*iBig);
    out.outerVisible = all(alignedLines >= rng(1) & alignedLines <= rng(2));
    out.alignedLines = alignedLines;
    if ~out.outerVisible
        out.why = sprintf(['the aligned pair (%.1f and %.1f MHz) falls outside the sweep, ' ...
            'so |B| is weakly constrained'], alignedLines(1), alignedLines(2));
    end

    dirHat = Bvec/max(norm(Bvec), eps);
    out.predictLines = @(scale) predictAll(scale*norm(Bvec)*dirHat, D, Ax, GAMMA);
    out.fieldForLineAt = @(target) solveScale(target, norm(Bvec), dirHat, D, Ax, GAMMA);
    out.alignedAxisIdx = iBig;

    if opts.verbose
        report(out, c, opts, GAMMA, rng(1), rng(2));
    end
end

% ======================================================================

function J = costOf(p, c, ce, Ax, GAMMA, opts, rng)
    Bvec = p(1:3)'; D = p(4);
    if abs(D - opts.D0) > opts.dTol
        J = 1e6*(1 + (abs(D - opts.D0) - opts.dTol)); return
    end
    pred = predictAll(Bvec, D, Ax, GAMMA);
    pred = pred(pred >= rng(1) & pred <= rng(2));
    if isempty(pred); J = 1e6; return; end
    h = opts.huberMHz;
    J = 0;
    for k = 1:numel(c)
        d = min(abs(pred - c(k)));
        J = J + min((d/ce(k))^2, (h/ce(k))^2);
    end
    % predictions inside the sweep with no measured line near them
    ceTyp = median(ce);
    for j = 1:numel(pred)
        if min(abs(c - pred(j))) > h
            J = J + opts.missWeight*(h/ceTyp)^2;
        end
    end
    J = J / numel(c);
end

function [pred, spl, mns, proj] = predictAll(Bvec, D, Ax, GAMMA)
% All eight transition frequencies, from the exact 3x3 Hamiltonian per orientation.
    Bvec = Bvec(:);
    Sz = diag([1 0 -1]);
    Sx = [0 1 0; 1 0 1; 0 1 0]/sqrt(2);
    proj = (Ax*Bvec)';
    B = norm(Bvec);
    pred = zeros(1, 8); spl = zeros(1,4); mns = zeros(1,4);
    for k = 1:4
        bz = proj(k);
        bperp = sqrt(max(B^2 - bz^2, 0));
        H = D*(Sz*Sz) + GAMMA*(bperp*Sx + bz*Sz);
        ev = sort(real(eig((H + H')/2)));
        pred(2*k-1) = ev(2) - ev(1);
        pred(2*k)   = ev(3) - ev(1);
        spl(k) = pred(2*k) - pred(2*k-1);
        mns(k) = (pred(2*k) + pred(2*k-1))/2;
    end
end

function [scale, which, fAll, sols] = solveScale(target, B0, dirHat, D, Ax, GAMMA)
% The |B| scaling that puts SOME line on target, by bisection on each line in turn.
%
% Solved rather than extrapolated: the splitting is only linear in |B| while the transverse
% shift is quadratic, so scaling a splitting linearly is wrong by several MHz at the fields
% in use here -- and this is the number the coil gets set from.
%
% ALL solutions are returned, not just the closest: which line you want on the target is a
% physics choice, not a numerical one. For Rabi you want the ALIGNED orientation, whose
% line may well need a different |B| from whichever line happens to reach the target first.
    scale = NaN; which = NaN; fAll = [];
    sols = struct('scale', {}, 'line', {}, 'axis', {}, 'lines', {});
    lo = 0.05; hi = 4;
    for j = 1:8
        f = @(s) pickLine(predictAll(s*B0*dirHat, D, Ax, GAMMA), j) - target;
        a = lo; b = hi; fa = f(a); fb = f(b);
        if ~isfinite(fa) || ~isfinite(fb) || fa*fb > 0; continue; end
        for it = 1:80
            m = 0.5*(a+b); fm = f(m);
            if fa*fm <= 0; b = m; fb = fm; else; a = m; fa = fm; end
        end
        sj = 0.5*(a+b);
        [pj, ~, ~, projj] = predictAll(sj*B0*dirHat, D, Ax, GAMMA);
        sols(end+1) = struct('scale', sj, 'line', j, 'axis', ceil(j/2), ...
                             'lines', sort(pj));                            %#ok<AGROW>
        if isnan(scale) || abs(sj - 1) < abs(scale - 1)
            scale = sj; which = j; fAll = sort(pj);
        end
        if false; disp(projj); end
    end
end

function v = pickLine(pred, j)
    v = pred(j);
end

function report(out, c, opts, GAMMA, rngLo, rngHi)
    fprintf('\n--- field from the exact Hamiltonian (4 parameters, %d lines) ---\n', numel(c));
    fprintf('  |B|      = %.3f G\n', out.B);
    fprintf('  D        = %.3f MHz\n', out.D);
    fprintf('  tilt from the nearest <111> axis = %.3f deg\n', out.tiltDeg);
    fprintf('  projections/G = %s\n', sprintf('%+9.3f', out.projections));
    fprintf('  residual rms = %.3f MHz over %d matched lines\n', out.rms, numel(out.matched));
    if isfield(out, 'outerVisible') && ~out.outerVisible
        fprintf(['\n  *** |B| IS NOT MEASURED HERE. The aligned pair is predicted at %.1f and\n' ...
                 '      %.1f MHz, both outside the swept range %.0f-%.0f. Without it only the\n' ...
                 '      off-axis lines constrain the field, and their positions trade off\n' ...
                 '      against its DIRECTION: on a clipped 45 G test spectrum a 55 G solution\n' ...
                 '      fitted to 0.03 MHz rms. Widen the sweep to include %.0f-%.0f MHz before\n' ...
                 '      believing the |B| above. The tilt and D are far less affected.\n'], ...
                out.alignedLines(1), out.alignedLines(2), rngLo, rngHi, ...
                min(out.alignedLines)-15, max(out.alignedLines)+15);
    end
    fprintf('\n%4s %12s %12s %10s   %s\n', 'axis', 'measured', 'predicted', 'resid', 'note');
    for k = 1:numel(out.matched)
        m = out.matched(k);
        fprintf('%4d %12.3f %12.3f %+10.3f\n', m.axis, m.centre, m.predicted, m.resid);
    end
    if ~isempty(out.missing)
        fprintf(['\n  %d predicted transition(s) inside the sweep were NOT seen:\n'], ...
                numel(out.missing));
        fprintf('    %.2f MHz\n', out.missing);
        fprintf(['    Usually too shallow or merged into a neighbour. If several are missing,\n' ...
                 '    the field solution is suspect -- check that the sweep contains the OUTER\n' ...
                 '    pair, without which half the true |B| fits the rest almost as well.\n']);
    end
    if ~isempty(out.unmatched)
        fprintf(['\n  %d line(s) match NO predicted transition within %.1f MHz:\n'], ...
                numel(out.unmatched), opts.huberMHz);
        fprintf('    %.2f MHz\n', out.unmatched);
        fprintf(['    These are not NV transitions for this field. The resonator dip is the\n' ...
                 '    usual candidate (it moves: seen near 2820 and near 2763), though\n' ...
                 '    esrFitRobust normally fits it separately -- see out.resonance.\n']);
    end
    seen = unique([out.matched.axis]);
    if numel(seen) < 4
        fprintf(['\n  only %d of the 4 orientations were seen. The fit still closes on the\n' ...
                 '  ones present -- the <111> geometry fixes the rest -- but |B| and the tilt\n' ...
                 '  are less constrained than the residual alone suggests.\n'], numel(seen));
    end
    fprintf('\n  predicted splittings/MHz = %s\n', sprintf('%9.3f', out.splittings));
    fprintf('  predicted pair means/MHz = %s\n', sprintf('%9.3f', out.means));
    fprintf(['  (the off-axis means sit %.2f MHz above D: the second-order transverse\n' ...
             '   shift, which is why a single mirror D cannot pair every orientation)\n'], ...
            max(out.means) - out.D);
end

function o = withDefaults(o, d)
    if isempty(o); o = struct(); end
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(o, fn{i}) || isempty(o.(fn{i}))
            o.(fn{i}) = d.(fn{i});
        end
    end
end
