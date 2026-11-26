function out = testCoilMatrixFit(opts)
%TESTCOILMATRIXFIT Grade the coil-matrix calibration against KNOWN coils. No hardware.
%
%   testCoilMatrixFit                         % line-level tests, ~20 s
%   testCoilMatrixFit(struct('fullChain', true))   % also synthetic spectra -> esrFitRobust
%
% Builds a coil system with everything the real one might have wrong -- the diamond turned
% relative to the coils, coil axes a few degrees off perpendicular, sensitivities 10-20% off
% the nominal ratios, and a background field -- then generates what calibCoilMatrix would
% measure at the design currents, and checks the fit recovers it.
%
% The number that matters is NOT how well M's entries come back, it is the one the
% calibration exists for: SET THE FIELD ALONG AN NV AXIS USING THE FITTED MODEL, THEN ASK
% THE TRUE COILS WHERE THE FIELD ACTUALLY POINTS. That angle is what alignment costs you,
% and it is gauge-free -- it does not care which of the 48 equivalent crystal frames the
% fit chose.
%
% TESTS
%   A  diamond ~5 deg off the coils: M, b0, D against truth, and the gauge -- with the
%      nominal (diagonal) matrix as reference the NV frame must come out as the coil-ish
%      one, i.e. g = identity
%   B  diamond turned 40 deg about a skew axis, nominal ratios 25% wrong: the start comes
%      from the data, not from the nominal model
%   C  gauge continuity: with a previous calibration as gaugeRef, the same frame comes back
%   D  only 8 spectra: the |B|^2 fit is skipped and the nominal start is used
%   E  (fullChain) synthetic SPECTRA through calibCoilMatrix's replay path: line finding by
%      esrFitRobust, the resonator dip, blends, everything
%
% Each measured spectrum is degraded the way real ones are: 0.15 MHz centre noise, lines
% closer than 4 MHz merged into one, lines outside the sweep lost, a line dropped at
% random in a third of the spectra, and a spurious line (the resonator) in a quarter.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('fullChain', false, 'nPoints', 16, ...
        'noiseMHz', 0.15, 'tolAlignDeg', 0.15, 'tolBmagG', 0.08, 'verbose', false));

    here = fileparts(mfilename('fullpath'));
    fn = fullfile(fileparts(fileparts(here)), 'functions');
    if isfolder(fn) && ~contains(path, fn); addpath(fn); end
    % the Helmholtz class (static methods only) for calibCoilMatrix's reach table
    src = fullfile(fileparts(fileparts(fileparts(here))), 'src');
    for p = {fullfile(src, 'Physical', 'MagneticStage'), fullfile(src, '_BaseClasses'), ...
             fullfile(src, '_BaseClasses', 'BaseObject'), fullfile(src, 'ClassPretty')}
        if isfolder(p{1}) && ~contains(path, p{1}); addpath(p{1}); end
    end

    Mnom = diag([5.382 8.034 13.649]);
    b0nom = [-0.373 0.179 2.121];
    Imax = [5.2 5.2 5.2];
    nFail = 0;
    out = struct();

    fprintf('\n==================================================================\n');
    fprintf('testCoilMatrixFit -- coil matrix from NV lines, against known coils\n');
    fprintf('==================================================================\n');

    % ---- A ----------------------------------------------------------------------
    T = makeTruth(struct('crystalAxis', [0.3 -0.5 1], 'crystalDeg', 5, 'sensErr', 0.12, ...
                         'skewDeg', 2, 'seed', 1));
    [r, f] = runCase('A  diamond ~5 deg off the coils', T, Mnom, b0nom, Imax, Mnom, opts);
    nFail = nFail + f;
    [gStar, ~] = bestGauge(r.M, T.M);
    ok = isequal(gStar, eye(3));
    nFail = nFail + ~ok;
    fprintf('    gauge: NV frame named like the coil frame (g = identity)   %s\n', okStr(ok));
    out.A = r;

    % ---- B ----------------------------------------------------------------------
    T = makeTruth(struct('crystalAxis', [1 2 -0.7], 'crystalDeg', 40, 'sensErr', 0.15, ...
                         'skewDeg', 3, 'seed', 2));
    [r, f] = runCase('B  diamond turned 40 deg, nominal ratios 25% off', T, Mnom*1.25, ...
                     b0nom, Imax, Mnom*1.25, opts);
    nFail = nFail + f; out.B = r;

    % ---- C ----------------------------------------------------------------------
    T = makeTruth(struct('crystalAxis', [-1 0.4 0.2], 'crystalDeg', 25, 'sensErr', 0.1, ...
                         'skewDeg', 2, 'seed', 3));
    g0 = [0 1 0; 0 0 -1; -1 0 0];                    % some arbitrary earlier naming
    [r, f] = runCase('C  gauge continuity with a previous calibration', T, Mnom, b0nom, ...
                     Imax, g0 * T.M, opts);
    nFail = nFail + f;
    ok = norm(r.M - g0 * T.M) / norm(T.M) < 0.01;
    nFail = nFail + ~ok;
    fprintf('    frame matches the previous calibration''s naming            %s\n', okStr(ok));
    out.C = r;

    % ---- D ----------------------------------------------------------------------
    T = makeTruth(struct('crystalAxis', [0.2 1 0.3], 'crystalDeg', 12, 'sensErr', 0.1, ...
                         'skewDeg', 2, 'seed', 4));
    o2 = opts; o2.nPoints = 8;
    [r, f] = runCase('D  only 8 spectra (nominal start)', T, Mnom, b0nom, Imax, Mnom, o2);
    nFail = nFail + f; out.D = r;

    % ---- E ----------------------------------------------------------------------
    if opts.fullChain
        T = makeTruth(struct('crystalAxis', [0.5 0.5 1], 'crystalDeg', 8, 'sensErr', 0.1, ...
                             'skewDeg', 2, 'seed', 5));
        f = runFullChain(T, Mnom, b0nom, Imax, opts);
        nFail = nFail + f;
    else
        fprintf('\nE  full chain (synthetic spectra -> esrFitRobust) skipped; run with\n');
        fprintf('   testCoilMatrixFit(struct(''fullChain'', true)), about 1-2 minutes.\n');
    end

    fprintf('\n------------------------------------------------------------------\n');
    if nFail == 0
        fprintf('ALL CHECKS PASSED\n');
    else
        fprintf('%d CHECK(S) FAILED -- do not trust a live calibration until these are clean\n', nFail);
    end
    out.nFail = nFail;
end

% ======================================================================

function T = makeTruth(o)
% A coil system: coil k makes field sens(k)*c_k per Ampere, c_k a few degrees off the
% coil-frame axis; the crystal frame is the coil frame turned by crystalDeg about
% crystalAxis; plus a background field.
    rs = RandStream('mt19937ar', 'Seed', o.seed);
    nomSens = [5.382 8.034 13.649];
    sens = nomSens .* (1 + o.sensErr * (2*rand(rs, 1, 3) - 1));
    C = eye(3);
    for k = 1:3
        v = C(:, k) + tand(o.skewDeg) * randn(rs, 3, 1) / sqrt(3);
        C(:, k) = sens(k) * v / norm(v);
    end
    ax = o.crystalAxis(:) / norm(o.crystalAxis);
    R = axisAngle(ax, o.crystalDeg);
    T.M = R' * C;                            % coil field expressed in the crystal frame
    T.b0 = (R' * [0.18; -0.27; 0.41])';      % ~0.5 G, the earth-sized background
    T.D = 2870.6;
end

function [r, nFail] = runCase(name, T, Mnom, b0nom, Imax, gaugeRef, opts)
    fprintf('\n%s\n', name);
    d = coilCurrentDesign(struct('M', Mnom, 'b0', b0nom, 'Imin', [0 0 0], 'Imax', Imax, ...
                                 'nPoints', opts.nPoints));
    data = simulate(T, d.I, opts.noiseMHz, 100 + round(1000*T.M(1)));
    r = coilMatrixFit(data, struct('nominalM', Mnom, 'nominalB0', b0nom, ...
                                   'gaugeRef', gaugeRef, 'verbose', opts.verbose));
    nFail = 0;
    if ~r.ok
        fprintf('    fit FAILED: %s\n', r.why);
        nFail = 1; return
    end
    [gStar, Mg] = bestGauge(r.M, T.M);
    b0g = (gStar * r.b0(:))';
    errM = max(abs(Mg(:) - T.M(:)));
    errB0 = max(abs(b0g - T.b0));
    fprintf('    %d spectra, %d lines, %d outliers, rms %.3f MHz, chi2 %.2f\n', ...
            r.nSpectra, r.nLines, r.nOutliers, r.rmsMHz, r.chi2Red);
    ok = errM < 0.05;
    fprintf('    max |dM|   %.4f G/A  (quoted max err %.4f)                %s\n', errM, ...
            max(r.Merr(:)), okStr(ok));
    nFail = nFail + ~ok;
    ok = errB0 < 0.15;
    fprintf('    max |db0|  %.4f G                                       %s\n', errB0, okStr(ok));
    nFail = nFail + ~ok;
    ok = abs(r.D - T.D) < 0.1;
    fprintf('    D          %.3f vs %.3f MHz                          %s\n', r.D, T.D, okStr(ok));
    nFail = nFail + ~ok;

    % THE NUMBER THAT MATTERS: aim the field along an NV axis with the fitted model.
    [worstAng, worstMag, nAx] = alignmentCheck(r, T, Imax);
    ok = nAx > 0 && worstAng < opts.tolAlignDeg && worstMag < opts.tolBmagG;
    fprintf(['    aimed at %d reachable NV axes (20 G): true field off-axis by <= %.3f deg,\n' ...
             '      |B| off by <= %.3f G                                       %s\n'], ...
            nAx, worstAng, worstMag, okStr(ok));
    nFail = nFail + ~ok;
    r.alignDeg = worstAng;
end

function [worstAng, worstMag, nAx] = alignmentCheck(r, T, Imax)
    AX = [1 1 1; 1 -1 -1; -1 1 -1; -1 -1 1] / sqrt(3);
    worstAng = 0; worstMag = 0; nAx = 0;
    for k = 1:4
        for s = [1 -1]
            Bt = 20 * s * AX(k, :)';
            Iq = r.M \ (Bt - r.b0(:));
            if any(Iq < 0) || any(Iq' > Imax); continue; end
            Btrue = T.M * Iq + T.b0(:);
            cosMax = max(abs(AX * Btrue)) / norm(Btrue);
            worstAng = max(worstAng, acosd(min(cosMax, 1)));
            worstMag = max(worstMag, abs(norm(Btrue) - 20));
            nAx = nAx + 1;
        end
    end
end

function data = simulate(T, Icur, noise, seed)
    rs = RandStream('mt19937ar', 'Seed', seed);
    sweep = [2740 3000];
    data = struct('I', {}, 'centres', {}, 'centreErr', {}, 'sweep', {});
    for j = 1:size(Icur, 1)
        B = T.M * Icur(j, :)' + T.b0(:);
        f = sort(nvLines(B, T.D));
        f = f(f > sweep(1) + 3 & f < sweep(2) - 3);
        % blends: lines closer than 4 MHz come out as one line at their mean
        merged = [];
        k = 1;
        while k <= numel(f)
            grp = f(k);
            while k < numel(f) && f(k+1) - f(k) < 4
                k = k + 1; grp(end+1) = f(k);                             %#ok<AGROW>
            end
            merged(end+1) = mean(grp);                                    %#ok<AGROW>
            k = k + 1;
        end
        f = merged + noise * randn(rs, size(merged));
        if rand(rs) < 0.33 && numel(f) > 4
            f(randi(rs, numel(f))) = [];
        end
        if rand(rs) < 0.25
            f(end+1) = 2763 + 2*randn(rs);                                %#ok<AGROW>
        end
        data(j) = struct('I', Icur(j, :), 'centres', sort(f), ...
                         'centreErr', noise * ones(size(f)), 'sweep', sweep);
    end
end

function nFail = runFullChain(T, Mnom, b0nom, Imax, opts)
    fprintf('\nE  full chain: synthetic spectra -> esrFitRobust -> coilMatrixFit\n');
    d = coilCurrentDesign(struct('M', Mnom, 'b0', b0nom, 'Imin', [0 0 0], 'Imax', Imax, ...
                                 'nPoints', opts.nPoints));
    freqs = 2740:1:3000;
    rep = struct('I', {}, 'freqs', {}, 'signal', {}, 'sterr', {});
    for j = 1:size(d.I, 1)
        B = T.M * d.I(j, :)' + T.b0(:);
        s = esrSynthSpectrum(struct('Bvec', B, 'D', T.D, 'freqs', freqs, 'seed', 50 + j, ...
                                    'resCentre', 2763, 'resWidth', 30, 'resDepth', 0.006, ...
                                    'verbose', false));
        rep(j) = struct('I', d.I(j, :), 'freqs', s.freqs, 'signal', s.signal, 'sterr', s.sterr);
    end
    t0 = tic;
    r = calibCoilMatrix(struct('replay', struct('spectra', rep), 'nominalM', Mnom, ...
                               'nominalB0', b0nom, 'plotIt', false, 'verbose', false));
    fprintf('    %.0f s for %d spectra\n', toc(t0), numel(rep));
    nFail = 0;
    if ~r.ok
        fprintf('    FAILED: %s\n', r.why); nFail = 1; return
    end
    fit = r.fit;
    fprintf('    %d spectra, %d lines, %d outliers, rms %.3f MHz\n', fit.nSpectra, ...
            fit.nLines, fit.nOutliers, fit.rmsMHz);
    [worstAng, worstMag, nAx] = alignmentCheck(fit, T, Imax);
    ok = nAx > 0 && worstAng < 2*opts.tolAlignDeg && worstMag < 2*opts.tolBmagG;
    fprintf(['    aimed at %d reachable NV axes: true field off-axis by <= %.3f deg, ' ...
             '|B| off by <= %.3f G   %s\n'], nAx, worstAng, worstMag, okStr(ok));
    nFail = nFail + ~ok;
end

% ======================================================================

function f = nvLines(B, D)
    GAMMA = 2.8025;
    AX = [1 1 1; 1 -1 -1; -1 1 -1; -1 -1 1] / sqrt(3);
    f = zeros(1, 8);
    for k = 1:4
        bz = AX(k, :) * B;
        bp = sqrt(max(B'*B - bz^2, 0));
        H = [D + GAMMA*bz, GAMMA*bp/sqrt(2), 0; GAMMA*bp/sqrt(2), 0, GAMMA*bp/sqrt(2); ...
             0, GAMMA*bp/sqrt(2), D - GAMMA*bz];
        ev = sort(eig(H));
        f(2*k-1) = ev(2) - ev(1); f(2*k) = ev(3) - ev(1);
    end
end

function [g, Mg] = bestGauge(M, Mtrue)
% The signed permutation mapping the fitted frame onto the true one.
    P = perms(1:3); best = inf; g = eye(3);
    for i = 1:6
        Pm = eye(3); Pm = Pm(P(i, :), :);
        for s = 0:7
            gi = diag(1 - 2*[bitget(s,1) bitget(s,2) bitget(s,3)]) * Pm;
            e = norm(gi * M - Mtrue, 'fro');
            if e < best; best = e; g = gi; end
        end
    end
    Mg = g * M;
end

function R = axisAngle(k, deg)
    K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    R = eye(3) + sind(deg) * K + (1 - cosd(deg)) * K * K;
end

function s = okStr(ok)
    if ok; s = 'ok'; else; s = 'FAIL'; end
end

function o = withDefaults(o, d)
    if isempty(o); o = struct(); end
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(o, fn{i}) || isempty(o.(fn{i})); o.(fn{i}) = d.(fn{i}); end
    end
end
