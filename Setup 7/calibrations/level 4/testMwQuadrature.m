function out = testMwQuadrature(opts)
%TESTMWQUADRATURE Grade the both-on balance test against a pair with a known imbalance.
%
%   out = testMwQuadrature                    % all cases, ~10 s
%   out = testMwQuadrature(struct('nPull', 40))
%
% No instruments. Plays AWGExpMwQuadrature's OWN configuration list and pulse definitions
% (configList, segmentsFor) through a two-level NV ensemble -- the 14N triplet at
% -2.16/0/+2.16 MHz, a Gaussian inhomogeneous width, and the two PCBs' fields falling off in
% opposite directions through the depth -- reads it out as 1 - contrast * P(ms=-1) plus
% noise, and hands the per-average ratios to mwQuadratureFit, exactly as calibMwQuadrature
% does. The truth is a channel-2 shortfall of epsDb at the played level; the test must
% recover delta + d* = delta + epsDb.
%
% CASES
%   A  the basic case at -32 dBm drive (1.72 MHz per channel)
%   B  a 5 degree channelPhase error               (must not move the answer)
%   C  +/-10% opposite field gradients from the PCBs (must not move it)
%   D  the other sign of imbalance, pulse times 10% long
%   E  strong drive (8 MHz per channel)
%   F  the opposite phase sign convention          (the zero crossing must not move)
%   G  pulls: the error bar covers the truth at the rate it claims
%
% The hyperfine leaves a small known bias at weak drive (+0.03-0.05 dB here; it falls to
% <0.01 dB at 8 MHz), so the tolerance on the cases is max(3 sigma, 0.08 dB).

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('nPull', 30, 'seed', 3, 'verbose', true), 'testMwQuadrature');
    % The generator is named: a bare rng(seed) errors once anything in the session has
    % called rand('state', k), and rng(seed, 'twister') draws the same numbers.
    rngOld = rng;  rng(opts.seed, 'twister');
    cleanup = onCleanup(@() rng(rngOld)); %#ok<NASGU>

    base = struct('omega1', 1.72, 'epsDb', 0.5, 'phaseErrDeg', 0, 'grad', 0, ...
                  'pulseFactor', 1, 'convention', -1, 'offsets', [-1 0 1], ...
                  'contrast', 0.03, 'noise', 2e-4, 'nAvg', 20, 'hyperfine', true);
    cases = {
        'A basic, -32 dBm',              struct()
        'B channelPhase error 5 deg',    struct('phaseErrDeg', 5)
        'C opposite 10% gradients',      struct('grad', 0.10)
        'D other sign, pulses 10% long', struct('epsDb', -0.7, 'pulseFactor', 1.1)
        'E strong drive 8 MHz',          struct('omega1', 8, 'epsDb', 0.3)
        'F opposite phase convention',   struct('convention', +1)
        };
    out = struct('pass', true, 'cases', []);
    fprintf('\n============================================================\n');
    fprintf('testMwQuadrature\n');
    fprintf('============================================================\n');
    for i = 1:size(cases, 1)
        s = merge(base, cases{i, 2});
        [f, truth] = simulate(s);
        tol = max(3 * f.dStarErr, 0.08);
        ok = f.ok && abs(f.dStar - truth) <= tol;
        fprintf('  %-32s truth %+.3f  got %+.3f +/- %.3f dB  slope %+.4f/dB of contrast  %s\n', ...
                cases{i, 1}, truth, f.dStar, f.dStarErr, f.slopeNorm, tf(ok));
        if ~ok && ~f.ok; fprintf('      (%s)\n', f.why); end
        out.cases(i).name = cases{i, 1};  out.cases(i).ok = ok;  out.cases(i).fit = f;
        out.pass = out.pass && ok;
    end

    % G: pulls -- is the quoted error honest?
    pulls = nan(1, opts.nPull);
    for k = 1:opts.nPull
        s = base;  s.epsDb = 0.6 * (2 * rand - 1);
        [f, truth] = simulate(s);
        if f.ok; pulls(k) = (f.dStar - truth - 0.04) / f.dStarErr; end   % 0.04: the hyperfine bias
    end
    pulls = pulls(isfinite(pulls));
    sp = std(pulls);
    okG = numel(pulls) >= 0.9 * opts.nPull && sp > 0.6 && sp < 1.5;
    fprintf('  %-32s %d runs: pull mean %+.2f, std %.2f (want ~0, ~1)  %s\n', ...
            'G pulls', numel(pulls), mean(pulls), sp, tf(okG));
    out.pulls = pulls;  out.pass = out.pass && okG;

    fprintf('\n%s\n', ternary(out.pass, 'ALL TESTS PASS', '*** SOME TESTS FAILED ***'));
end

% ======================================================================

function [f, truth] = simulate(s)
    cfg = AWGExpMwQuadrature.configList(s.offsets);
    % Fields at the NVs, per ensemble member: channel 1 as set, channel 2 short by epsDb.
    [D, h1, h2, wgt] = ensemble(s);
    O1 = s.omega1;
    O2 = s.omega1 * 10^(-s.epsDb / 20);
    piPair = 1 / (2 * (O1 + O2));                       % us, the in-phase pair's pi
    T2 = s.pulseFactor * piPair / 2;
    T1 = s.pulseFactor * sqrt(2) * piPair / 2;
    phErr = deg2rad(s.phaseErrDeg);
    P = zeros(numel(cfg), 1);
    for i = 1:numel(cfg)
        segs1 = AWGExpMwQuadrature.segmentsFor(cfg(i), 1, T1, T2);
        segs2 = AWGExpMwQuadrature.segmentsFor(cfg(i), 2, T1, T2);
        % The hardware sees channel c's waveform cos(w t - phase - extra) as a drive of
        % amplitude k*Omega_c and phase convention*(phase + extra); channel 2 carries the
        % residual channelPhase error.
        edges = unique([segs1(:, 1); segs1(:, 2); segs2(:, 1); segs2(:, 2)])';
        a = ones(size(D));  b = zeros(size(D));          % state |0> = (a, b)
        for e = 1:numel(edges) - 1
            t0 = edges(e);  t1 = edges(e + 1);  dt = t1 - t0;
            Om = zeros(size(D));
            Om = Om + O1 * h1 .* segAmp(segs1, t0, s.convention, 0);
            Om = Om + O2 * h2 .* segAmp(segs2, t0, s.convention, phErr);
            [a, b] = evolve(a, b, Om, D, dt);
        end
        P(i) = sum(wgt .* abs(b).^2);
    end
    nA = s.nAvg;
    R = 1 - s.contrast * P + s.noise * randn(numel(cfg), nA);
    f = mwQuadratureFit(cfg, R);
    truth = s.epsDb;
end

function A = segAmp(segs, t0, conv, phErr)
    % Complex drive factor of the segment that is on at t0 (k * exp(conv*i*(extra + err))).
    A = 0;
    for s = 1:size(segs, 1)
        if t0 >= segs(s, 1) - 1e-12 && t0 < segs(s, 2) - 1e-12
            A = segs(s, 3) * exp(1i * conv * (segs(s, 4) + phErr));
        end
    end
end

function [D, h1, h2, wgt] = ensemble(s)
    % Deterministic quadrature: hyperfine lines x Gauss-Hermite detuning x depth points.
    if s.hyperfine; hf = [-2.16 0 2.16]; else; hf = 0; end
    xg = [-2.0201828705 -0.9585724646 0 0.9585724646 2.0201828705];     % Gauss-Hermite, 5 pts
    wg = [0.0199532421 0.3936193232 0.9453087205 0.3936193232 0.0199532421];
    sig = 1.0 / 2.355;                                                  % 1 MHz FWHM
    zz = linspace(-1, 1, 5);
    [H, X, Z] = ndgrid(hf, xg, zz);
    [~, WX, ~] = ndgrid(hf, wg, zz);
    D = H(:) + sqrt(2) * sig * X(:);
    wgt = WX(:);  wgt = wgt / sum(wgt);
    h1 = 1 + s.grad * Z(:);
    h2 = 1 - s.grad * Z(:);
end

function [a, b] = evolve(a, b, Om, D, t)
    % Two-level rotation for Rabi frequency Om (complex, MHz) and detuning D (MHz), t in us.
    ax = real(Om);  ay = imag(Om);  az = D;
    W = sqrt(ax.^2 + ay.^2 + az.^2) + 1e-15;
    th = pi * W * t;
    c = cos(th);  sn = sin(th);
    nx = ax ./ W;  ny = ay ./ W;  nz = az ./ W;
    a2 = (c - 1i * sn .* nz) .* a + (-1i * sn .* (nx - 1i * ny)) .* b;
    b2 = (-1i * sn .* (nx + 1i * ny)) .* a + (c + 1i * sn .* nz) .* b;
    a = a2;  b = b2;
end

function s = merge(a, b)
    s = a;
    f = fieldnames(b);
    for i = 1:numel(f); s.(f{i}) = b.(f{i}); end
end

function s = tf(v)
    if v; s = 'ok'; else; s = '*** FAIL ***'; end
end

function o = ternary(c, a, b)
    if c; o = a; else; o = b; end
end

function o = withDefaults(o, d, who)
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d); given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('%s: unknown option "%s". Valid: %s.', who, bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || isempty(o.(known{i})); o.(known{i}) = d.(known{i}); end
    end
end
