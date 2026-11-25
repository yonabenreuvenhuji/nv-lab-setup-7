function out = esrSynthSpectrum(opts)
%ESRSYNTHSPECTRUM An NV ESR sweep with KNOWN field, from the exact spin-1 Hamiltonian.
%
%   s = esrSynthSpectrum(struct('B', 45, 'tiltDeg', 7))
%   sv = calibFieldSurvey(struct('replay', s, 'plotIt', false));
%
% This exists because every change to the fitter has so far been tested by running the
% magnet. That is slow, it ties up the setup, and -- the real problem -- a live run gives
% no ground truth: when the analysis reports |B| = 44.8 G and a 6.5 degree tilt there is
% nothing to check it against, so a fitter bug and a genuine measurement look identical.
% Here the field is an input, so every number the survey prints can be graded.
%
% THE LINE POSITIONS ARE EXACT, not the usual D +/- gamma*B_par. For each orientation the
% 3x3 Hamiltonian
%
%       H = D*Sz^2 + gamma*(Bx*Sx + By*Sy + Bz*Sz)        (z along that NV axis)
%
% is diagonalised and the two transition frequencies are taken as E2-E1 and E3-E1. That
% automatically contains the second-order transverse-field shift -- the effect that moves
% an off-axis pair's MEAN up by ~7 MHz at 45 G while leaving its SPLITTING almost
% untouched, and therefore the effect the survey's own consistency check is testing. A
% synthetic built from the first-order formula would make that check fail on physically
% perfect data, which is exactly what happened the first time this was tried.
%
% OPTIONS
%   B          field magnitude, Gauss (default 45)
%   tiltDeg    angle between B and the NEAREST <111> axis (default 7)
%   tiltAzDeg  which way it is tilted, degrees about that axis (default 20). Changes how
%              the three off-axis orientations spread; 0 and 60 are the symmetric cases
%              that keep two of them degenerate, so the default is deliberately neither.
%   D          zero-field splitting, MHz (default 2870.5)
%   freqs      sweep, MHz (default 2725:1.5:3015 -- the grid actually in use)
%   fwhm       intrinsic linewidth, MHz (default 8)
%   depth      contrast of ONE line, in fluorescence units (default 0.010)
%   noise      point-to-point noise, same units (default 0.0009, as measured)
%   sterrScale sterr is reported this many times the true noise (default 8 -- AWGExpESR's
%              sterr is drift-dominated and measured far too large, and the fitter's
%              Birge scaling has to cope with that, so the synthetic reproduces it)
%   resDepth   depth of a broad resonator dip (default 0.004)
%   resCentre  where that dip sits, MHz (default 2820)
%   resWidth   its Gaussian sigma, MHz (default 90)
%   seed       RNG seed (default 12345)
%   verbose    print the truth table (default true; testFieldPipeline turns it off)
%
% OUTPUT is a struct that calibFieldSurvey accepts as 'replay', plus the truth:
%   freqs, signal, sterr        the sweep
%   trueLines                   every transition frequency, MHz
%   trueSplittings, trueMeans   per orientation
%   trueProjections             B.a_i in Gauss, signed
%   B, tiltDeg, D               what went in

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('B', 45, 'tiltDeg', 7, 'tiltAzDeg', 20, ...
        'D', 2870.5, 'freqs', 2725:1.5:3015, 'fwhm', 8, 'depth', 0.010, ...
        'noise', 0.0009, 'sterrScale', 8, 'resDepth', 0.004, 'resCentre', 2820, ...
        'resWidth', 90, 'seed', 12345, 'verbose', true));

    GAMMA = 2.8025;                 % MHz/G

    % The four <111> axes, normalised. Their sum is zero, which is the identity the
    % survey's projection solver relies on, so it is worth having the synthetic built
    % from the same vectors rather than from a paraphrase of them.
    Ax = [ 1  1  1; 1 -1 -1; -1  1 -1; -1 -1  1] / sqrt(3);

    % Field direction: start along axis 1, tilt by tiltDeg about a perpendicular
    % direction chosen by tiltAzDeg.
    a1 = Ax(1,:)';
    u  = [0;0;1] - (a1'*[0;0;1])*a1;  if norm(u) < 1e-9; u = [0;1;0] - (a1'*[0;1;0])*a1; end
    u  = u/norm(u);
    w  = cross(a1, u);
    th = opts.tiltDeg*pi/180; az = opts.tiltAzDeg*pi/180;
    nHat = cos(th)*a1 + sin(th)*(cos(az)*u + sin(az)*w);
    nHat = nHat/norm(nHat);
    Bvec = opts.B * nHat;

    proj = (Ax * Bvec)';            % signed projections, Gauss

    % ---- exact transition frequencies, one orientation at a time -------------
    Sz = diag([1 0 -1]);
    Sx = [0 1 0; 1 0 1; 0 1 0]/sqrt(2);
    Sy = [0 -1i 0; 1i 0 -1i; 0 1i 0]/sqrt(2);
    lines = []; splittings = zeros(1,4); means = zeros(1,4);
    for k = 1:4
        bz = proj(k);
        bperp = sqrt(max(opts.B^2 - bz^2, 0));
        H = opts.D*(Sz*Sz) + GAMMA*(bperp*Sx + bz*Sz);
        ev = sort(real(eig((H + H')/2)));
        fLo = ev(2) - ev(1);
        fHi = ev(3) - ev(1);
        lines = [lines fLo fHi];                              %#ok<AGROW>
        splittings(k) = fHi - fLo;
        means(k)      = (fHi + fLo)/2;
    end
    lines = sort(lines);

    % ---- build the sweep -----------------------------------------------------
    f = opts.freqs(:);
    y = 1 - opts.resDepth*exp(-((f - opts.resCentre)/opts.resWidth).^2) ...
          - 2e-5*(f - mean(f))/10;
    g = opts.fwhm/2;
    for k = 1:numel(lines)
        y = y - opts.depth * g^2 ./ ((f - lines(k)).^2 + g^2);
    end
    try
        rs = RandStream('mt19937ar', 'Seed', opts.seed);
        y = y + opts.noise*randn(rs, size(f));
    catch
        randn('seed', opts.seed);                                          %#ok<RAND>
        y = y + opts.noise*randn(size(f));
    end

    out = struct();
    out.freqs  = f;
    out.signal = y;
    out.sterr  = repmat(opts.noise*opts.sterrScale, size(f));
    out.trueLines      = lines;
    out.trueSplittings = splittings;
    out.trueMeans      = means;
    out.trueProjections = proj;
    out.B = opts.B; out.tiltDeg = opts.tiltDeg; out.D = opts.D;
    out.opts = opts;

    if ~opts.verbose; return; end
    fprintf('\n=== esrSynthSpectrum: |B| = %.2f G, tilt %.2f deg (az %.0f), D = %.1f ===\n', ...
            opts.B, opts.tiltDeg, opts.tiltAzDeg, opts.D);
    fprintf('projections/G : %s\n', sprintf('%+8.3f', proj));
    fprintf('%4s %10s %10s %11s %10s\n', 'axis', 'f-/MHz', 'f+/MHz', 'splitting', 'mean');
    for k = 1:4
        fprintf('%4d %10.3f %10.3f %11.3f %10.3f\n', k, means(k)-splittings(k)/2, ...
                means(k)+splittings(k)/2, splittings(k), means(k));
    end
    fprintf(['NOTE the off-axis means sit ABOVE D by up to %.2f MHz: that is the ' ...
             'second-order\n     transverse shift, and it is in the data on purpose.\n'], ...
            max(means) - opts.D);
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
