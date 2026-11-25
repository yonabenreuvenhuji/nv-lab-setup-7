function out = testFieldPipeline(opts)
%TESTFIELDPIPELINE Grade the whole field-survey chain against a KNOWN field. No hardware.
%
%   testFieldPipeline
%   testFieldPipeline(struct('tilts', [0 2 4 7 10], 'B', 45))
%
% Builds synthetic ESR sweeps from the exact NV Hamiltonian with esrSynthSpectrum, pushes
% each one through calibFieldSurvey in replay mode, and compares what comes out with what
% went in. Nothing here touches the AWG, the magnet or the counter.
%
% WHY THIS EXISTS. Every earlier change to this analysis was tested by running the magnet,
% which is slow, ties up the setup, and -- the real problem -- provides no truth to check
% against: when the survey reports 44.8 G and a 6.5 degree tilt there is nothing to compare
% it with, so a fitter bug and a real measurement look identical. Two defects that had been
% live for several runs were found within minutes of having a known answer to grade against:
%
%   1. Seeded refits started every line's depth at the noise floor, because the amplitude
%      start point came from a 7-point running maximum -- a 10 MHz window, narrower than the
%      lines themselves, so the "local maximum" was a point on the dip. Refits therefore
%      began from near-zero depths and could not beat the unseeded fit, which is why a
%      blended line never split however it was seeded.
%
%   2. The mirror-partner classifier discarded the ALIGNED pair, the most reliable pair in
%      the spectrum. There is no single D that pairs all four orientations once the
%      second-order transverse shift is included, and a D estimated from the data lands on
%      the value shared by the three off-axis pairs.
%
% Neither was visible in a live run. Both are obvious here.
%
% TESTS
%   A  field recovery across a range of tilts: |B|, tilt and D against truth
%   B  a simulated theta scan: does the alignment metric's parabola vertex land on the
%      known optimum, with a known residual left in the other angle
%   C  robustness: a spurious non-NV line, and a sweep that clips the outer pair. The
%      second is not a test that |B| survives -- it cannot, the aligned pair is the only
%      strong constraint on |B| -- but that the code reports the result as unconstrained
%      instead of presenting a confident wrong number.
%
% Returns a struct of results and prints a PASS/FAIL summary. Tolerances are set from what
% the method can actually do, not from what would be nice: 0.2 G on |B|, 0.5 deg on tilt
% ABOVE 2.5 degrees (below that, line positions genuinely stop carrying the tilt -- see the
% table in test A), and 0.3 deg on the scan vertex.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('tilts', [0 1 2 4 7 10], 'B', 45, 'az', 20, ...
        'phiRes', 1.5, 'scanOffsets', -3:1:3, 'tolB', 0.2, 'tolTilt', 0.5, ...
        'tolVertex', 0.3, 'verbose', false));

    out = struct();
    nFail = 0;

    % ---- A: field recovery ------------------------------------------------
    fprintf('\n==================================================================\n');
    fprintf('TEST A: field recovery from a known synthetic spectrum\n');
    fprintf('==================================================================\n');
    fprintf('%9s %9s %9s %9s %9s %8s %8s %8s\n', ...
            'tilt in', '|B| out', 'tilt out', 'D out', 'span', 'nLines', 'rms', 'verdict');
    A = struct('tiltIn', {}, 'B', {}, 'tilt', {}, 'D', {}, 'span', {}, 'nLines', {}, 'rms', {});
    for tIn = opts.tilts
        s = esrSynthSpectrum(struct('B', opts.B, 'tiltDeg', tIn, 'tiltAzDeg', opts.az, 'verbose', false));
        r = runQuiet(s, opts.verbose);
        [sp, w0] = spanOfFit(r);
        bErr = abs(r.B - opts.B);
        tErr = abs(r.tiltDeg - tIn);
        % Above 2.5 deg the tilt must be right. Below it, only |B| is graded: the three
        % off-axis orientations are then unresolved and their POSITIONS carry no tilt
        % information at all, which is exactly why the scan metric is a width and not this.
        if tIn > 2.5
            ok = bErr < opts.tolB && tErr < opts.tolTilt;
        else
            ok = bErr < opts.tolB;
        end
        nFail = nFail + ~ok;
        fprintf('%9.2f %9.3f %9.3f %9.2f %9.2f %8d %8.2f %8s\n', ...
                tIn, r.B, r.tiltDeg, r.D, sp, r.nLines, r.rms, verdict(ok));
        A(end+1) = struct('tiltIn', tIn, 'B', r.B, 'tilt', r.tiltDeg, 'D', r.D, ...
                          'span', sp, 'nLines', r.nLines, 'rms', r.rms);      %#ok<AGROW>
        if isfinite(w0); end
    end
    out.A = A;
    fprintf(['\n  Read the tilt column against the tilt-in column. It tracks down to about\n' ...
             '  2.5 deg and then collapses toward zero -- that is physics, not a bug: below\n' ...
             '  that the off-axis splittings differ by less than a linewidth. The span column\n' ...
             '  keeps moving all the way to zero tilt, which is why the alignment scan is run\n' ...
             '  on the span and not on the tilt.\n']);

    % ---- B: simulated alignment scan --------------------------------------
    fprintf('\n==================================================================\n');
    fprintf('TEST B: alignment scan, %g deg deliberately left in the other angle\n', opts.phiRes);
    fprintf('==================================================================\n');
    fprintf('%9s %10s %9s %11s\n', 'd theta', 'true tilt', 'span', 'metric');
    mv = []; xs = [];
    for d = opts.scanOffsets
        tt = sqrt(d^2 + opts.phiRes^2);
        az = atan2d(opts.phiRes, d); if d == 0; az = 90; end
        s = esrSynthSpectrum(struct('B', opts.B, 'tiltDeg', tt, 'tiltAzDeg', az, 'verbose', false));
        r = runQuiet(s, opts.verbose);
        [sp, w0] = spanOfFit(r);
        if isfinite(sp) && isfinite(w0); m = max(sp - w0, 0)^2; else; m = NaN; end
        fprintf('%9.1f %10.3f %9.2f %11.2f\n', d, tt, sp, m);
        mv(end+1) = m; xs(end+1) = d;                                         %#ok<AGROW>
    end
    g = isfinite(mv);
    vert = NaN;
    if sum(g) >= 4
        x = xs(g)'; y = mv(g)';
        b = [ones(numel(x),1), x, x.^2] \ y;
        if b(3) > 0; vert = -b(2)/(2*b(3)); end
    end
    okB = isfinite(vert) && abs(vert) < opts.tolVertex;
    nFail = nFail + ~okB;
    fprintf('\n  parabola vertex at d theta = %+.3f deg, truth 0.000  -> %s\n', vert, verdict(okB));
    out.B_vertex = vert;

    % ---- C: robustness ----------------------------------------------------
    fprintf('\n==================================================================\n');
    fprintf('TEST C: a spurious non-NV line, and a clipped sweep\n');
    fprintf('==================================================================\n');
    s = esrSynthSpectrum(struct('B', opts.B, 'tiltDeg', 7, 'tiltAzDeg', opts.az, 'verbose', false));
    % A resonator-like dip: broad, shallow, and at no NV frequency.
    gRes = 12;
    s.signal = s.signal - 0.006 * gRes^2 ./ ((s.freqs - 2800).^2 + gRes^2);
    r = runQuiet(s, opts.verbose);
    okC1 = abs(r.B - opts.B) < 3*opts.tolB;
    fprintf('  with a spurious dip at 2800: |B| %.3f (truth %.1f), tilt %.2f -> %s\n', ...
            r.B, opts.B, r.tiltDeg, verdict(okC1));

    % Clip the sweep so the outer pair falls outside it. |B| should then be reported from
    % the off-axis lines alone, which is much weaker -- the point of the test is that it
    % does not silently report a confident wrong answer.
    s2 = esrSynthSpectrum(struct('B', opts.B, 'tiltDeg', 7, 'tiltAzDeg', opts.az, ...
                                 'freqs', 2790:1.5:2950, 'verbose', false));
    r2 = runQuiet(s2, opts.verbose);
    % The pass condition is NOT that |B| comes back right. With the aligned pair outside the
    % sweep the data genuinely does not determine |B| -- a 55 G solution fits those six lines
    % to 0.03 MHz -- so demanding 45 G would be demanding that the code guess well. What is
    % required is that it SAYS SO: outerVisible must come back false.
    flagged = isstruct(r2.solve) && isfield(r2.solve, 'outerVisible') && ~r2.solve.outerVisible;
    okC2 = ~r2.ok || flagged;
    fprintf('  outer pair clipped out:      |B| %.3f (truth %.1f), tilt %.2f, rms %.2f\n', ...
            r2.B, opts.B, r2.tiltDeg, r2.rms);
    fprintf('    flagged as unconstrained: %d -> %s\n', flagged, verdict(okC2));
    nFail = nFail + ~okC1 + ~okC2;
    out.C = struct('spurious', r, 'clipped', r2);

    fprintf('\n==================================================================\n');
    if nFail == 0
        fprintf('ALL CHECKS PASSED\n');
    else
        fprintf('%d CHECK(S) FAILED -- do not trust a live run until these are clean\n', nFail);
    end
    fprintf('==================================================================\n\n');
    out.nFail = nFail;
end

% ======================================================================

function r = runQuiet(s, verbose)
% One survey in replay mode, reduced to the numbers the tests grade.
    o = struct('replay', s, 'plotIt', false, 'verboseRefine', false);
    if verbose
        sv = calibFieldSurvey(o);
    else
        txt = evalc('sv = calibFieldSurvey(o);');
        if isempty(txt); end
    end
    r = struct('ok', false, 'B', NaN, 'tiltDeg', NaN, 'D', NaN, 'rms', NaN, ...
               'nLines', 0, 'fit', [], 'solve', []);
    if ~isstruct(sv); return; end
    if isfield(sv, 'fit') && isstruct(sv.fit) && isfield(sv.fit, 'nLines')
        r.nLines = sv.fit.nLines; r.fit = sv.fit;
    end
    if isfield(sv, 'solve') && isstruct(sv.solve) && isfield(sv.solve, 'ok') && sv.solve.ok
        r.ok = true; r.B = sv.solve.B; r.tiltDeg = sv.solve.tiltDeg;
        r.D = sv.solve.D; r.rms = sv.solve.rms; r.solve = sv.solve;
    end
end

function [span, w0] = spanOfFit(r)
% The same off-axis span calibFieldAlign minimises, computed the same way: the ALIGNED pair
% is identified by the field fit's axis assignment, not by taking the lowest and highest
% line. Position-based selection breaks whenever a spurious line lands outside the real
% outer pair -- the true outer line is then counted as part of the off-axis group and the
% span jumps to the full aligned splitting. One such point in a scan moves the parabola
% vertex by degrees.
    span = NaN; w0 = NaN;
    if ~isstruct(r.fit) || ~isfield(r.fit, 'lines') || isempty(r.fit.lines); return; end
    cc = [r.fit.lines.centre]; ww = [r.fit.lines.fwhm];
    if numel(cc) < 3; return; end
    D = r.D; if ~isfinite(D); D = (min(cc)+max(cc))/2; end
    isO = [];
    if isstruct(r.solve) && isfield(r.solve, 'matched') && ~isempty(r.solve.matched) ...
            && isfield(r.solve, 'alignedAxisIdx')
        mc = [r.solve.matched.centre]; ma = [r.solve.matched.axis];
        isO = false(size(cc)); known = false(size(cc));
        for k = 1:numel(cc)
            [dd, j] = min(abs(mc - cc(k)));
            if dd < 1e-6
                known(k) = true;
                isO(k) = (ma(j) == r.solve.alignedAxisIdx);
            end
        end
        if any(isO) && any(known & ~isO)
            cc = cc(known); ww = ww(known); isO = isO(known);
        else
            isO = [];
        end
    end
    if isempty(isO); isO = (cc == min(cc)) | (cc == max(cc)); end
    w0 = min(ww(isO));
    sides = [NaN NaN];
    for k = 1:2
        if k == 1; m = ~isO & cc < D; else; m = ~isO & cc > D; end
        if any(m); sides(k) = max(cc(m)+ww(m)/2) - min(cc(m)-ww(m)/2); end
    end
    sides = sides(isfinite(sides));
    if ~isempty(sides); span = min(sides); end
end

function s = verdict(ok)
    if ok; s = 'ok'; else; s = 'FAIL'; end
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
