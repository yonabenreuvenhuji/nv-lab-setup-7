function out = testHelmholtzLimits(opts)
%TESTHELMHOLTZLIMITS Grade the Helmholtz field limits and angle convention. No hardware.
%
%   testHelmholtzLimits
%   testHelmholtzLimits(struct('nCases', 200))
%
% The GUI's "how far can Bx (or Brho, Btheta, Bphi) go with the other two held" bounds used
% to be computed as if the reachable fields were a box in B -- which they are not: they are
% the parallelepiped M*[current box] + b0, and with a full coil matrix not even aligned with
% the axes. Helmholtz.fieldLimits now solves it exactly. This checks that against brute
% force: for random coil matrices, offsets and working points, every coordinate is swept on
% a fine grid, each grid value is tested for reachability directly (are the currents inside
% the box?), and the connected run containing the present value is compared with the
% closed-form interval.
%
% Also checks the angle convention -- theta POLAR from +z, phi AZIMUTH from +x -- on the
% one case that matters most: theta = 54.7356, phi = 45 must be the [1 1 1] NV axis.
%
% Uses only Helmholtz's static methods, so no supplies, no GUI, no setupInfo.json.

    if ~exist('opts', 'var'); opts = struct(); end
    if ~isfield(opts, 'nCases'); opts.nCases = 60; end
    if ~isfield(opts, 'tolDeg'); opts.tolDeg = 0.02; end

    here = fileparts(mfilename('fullpath'));
    repo = fileparts(fileparts(fileparts(here)));
    src = fullfile(repo, 'src');
    for p = {fullfile(src, 'Physical', 'MagneticStage'), fullfile(src, '_BaseClasses'), ...
             fullfile(src, '_BaseClasses', 'BaseObject'), fullfile(src, 'ClassPretty')}
        if isfolder(p{1}) && ~contains(path, p{1}); addpath(p{1}); end
    end

    nFail = 0;
    fprintf('\n==================================================================\n');
    fprintf('testHelmholtzLimits -- exact field limits vs brute force\n');
    fprintf('==================================================================\n');

    % ---- the convention ---------------------------------------------------------
    u = Helmholtz.sph2cartPhys(1, Helmholtz.NV_THETA, 45);
    ok = norm(u - [1 1 1]/sqrt(3)) < 1e-12;
    fprintf('theta 54.7356 (polar), phi 45 (azimuth) is [1 1 1]/sqrt(3)        %s\n', okStr(ok));
    nFail = nFail + ~ok;
    [r, t, p] = Helmholtz.cart2sphPhys([0 0 5]);
    ok = abs(r - 5) < 1e-12 && abs(t) < 1e-12;
    fprintf('B along +z has theta = 0                                          %s\n', okStr(ok));
    nFail = nFail + ~ok;
    [~, t, p] = Helmholtz.cart2sphPhys([0 3 0]);
    ok = abs(t - 90) < 1e-12 && abs(p - 90) < 1e-12;
    fprintf('B along +y has theta = 90, phi = 90                               %s\n', okStr(ok));
    nFail = nFail + ~ok;
    dirs = Helmholtz.nvDirections();
    [~, t, p] = Helmholtz.cart2sphPhys(dirs(1:4, :));
    ok = all(abs(t - 54.7356) < 1e-3) && all(abs(abs(p) - 45) < 1e-9 | abs(abs(p) - 135) < 1e-9);
    fprintf('the four upward NV axes all have theta 54.7356, phi +-45 / +-135  %s\n', okStr(ok));
    nFail = nFail + ~ok;
    f = Helmholtz.nvLines([0 0 0], 2870);
    ok = all(abs(f - 2870) < 1e-6);
    fprintf('zero field: all 8 lines at D                                       %s\n', okStr(ok));
    nFail = nFail + ~ok;
    f = Helmholtz.nvLines(30 * [1 1 1]/sqrt(3), 2870);
    ok = abs((f(end) - f(1)) - 2*Helmholtz.GAMMA*30) < 0.5;
    fprintf('30 G along [111]: outer splitting 2*gamma*30 (to 2nd order)        %s\n', okStr(ok));
    nFail = nFail + ~ok;

    % ---- limits vs brute force ----------------------------------------------------
    rs = RandStream('mt19937ar', 'Seed', 11);
    whats = {'Bx', 'By', 'Bz', 'Brho', 'Btheta', 'Bphi'};
    worst = struct();
    for w = whats; worst.(w{1}) = 0; end
    nEmpty = 0;
    for n = 1:opts.nCases
        % a coil matrix: nominal-ish sensitivities, a random crystal rotation, some skew
        sens = [5.4 8.0 13.6] .* (1 + 0.2*(2*rand(rs, 1, 3) - 1));
        C = diag(sens) + 0.4 * randn(rs, 3);
        q = randn(rs, 4, 1); q = q / norm(q);
        R = quat(q);
        if n <= 5; R = eye(3); end
        M = R * C;
        b0 = 0.5 * randn(rs, 1, 3);
        Imin = [0 0 0]; Imax = [5.2 5.2 5.2];
        Inow = Imin + rand(rs, 1, 3) .* (Imax - Imin);
        if n == 6; Inow(2) = Imax(2); end            % on a limit
        if n == 7; Inow(1) = Imin(1); end
        Bnow = (M * Inow' + b0')';
        Minv = M \ eye(3);
        reach = @(B) all((Minv * (B' - b0'))' >= Imin - 1e-6 & (Minv * (B' - b0'))' <= Imax + 1e-6, 2);
        [rho, th, ph] = Helmholtz.cart2sphPhys(Bnow);
        for w = whats
            [lo, hi] = Helmholtz.fieldLimits(M, b0, Imin, Imax, Bnow, w{1});
            % brute force: a fine grid, the connected run of reachable values holding now
            switch w{1}
                case {'Bx', 'By', 'Bz'}
                    c = w{1}(2) - 'w';
                    g = Bnow(c) + linspace(-150, 150, 60001);
                    Bs = repmat(Bnow, numel(g), 1); Bs(:, c) = g;
                    now = Bnow(c);
                case 'Brho'
                    g = linspace(0, 200, 40001);
                    Bs = Helmholtz.sph2cartPhys(g, th*ones(size(g)), ph*ones(size(g)));
                    now = rho;
                case 'Btheta'
                    g = linspace(0, 180, 36001);
                    Bs = Helmholtz.sph2cartPhys(rho*ones(size(g)), g, ph*ones(size(g)));
                    now = th;
                case 'Bphi'
                    g = ph + linspace(-360, 360, 144001);
                    Bs = Helmholtz.sph2cartPhys(rho*ones(size(g)), th*ones(size(g)), g);
                    now = ph;
            end
            okv = reach(Bs);
            [~, i0] = min(abs(g - now));
            if ~okv(i0)
                % the present point is on the edge and the nearest grid value fell outside
                if i0 > 1 && okv(i0 - 1); i0 = i0 - 1; elseif i0 < numel(g) && okv(i0 + 1); i0 = i0 + 1; end
            end
            if ~okv(i0)
                if ~isnan(lo); nEmpty = nEmpty + 1; end
                continue
            end
            iL = i0; while iL > 1 && okv(iL - 1); iL = iL - 1; end
            iR = i0; while iR < numel(g) && okv(iR + 1); iR = iR + 1; end
            bLo = g(iL); bHi = g(iR);
            step = g(2) - g(1);
            if strcmp(w{1}, 'Bphi') && bHi - bLo >= 359
                % the whole circle is reachable
                err = 0;
                if ~(lo <= -179.99 && hi >= 179.99) && ~(hi - lo >= 359.9); err = inf; end
            else
                err = max(abs(lo - bLo), abs(hi - bHi)) - step;
            end
            worst.(w{1}) = max(worst.(w{1}), err);
        end
    end
    tol = struct('Bx', 0.01, 'By', 0.01, 'Bz', 0.01, 'Brho', 0.01, 'Btheta', opts.tolDeg, ...
                 'Bphi', opts.tolDeg);
    fprintf('\n%d random coil systems and working points, each coordinate swept:\n', opts.nCases);
    for w = whats
        ok = worst.(w{1}) <= tol.(w{1});
        nFail = nFail + ~ok;
        fprintf('  %-7s worst disagreement with brute force %9.4f (tol %.3f)     %s\n', w{1}, ...
                max(worst.(w{1}), 0), tol.(w{1}), okStr(ok));
    end
    ok = nEmpty == 0;
    fprintf('  intervals reported where brute force found none: %d                %s\n', ...
            nEmpty, okStr(ok));
    nFail = nFail + ~ok;

    fprintf('\n------------------------------------------------------------------\n');
    if nFail == 0
        fprintf('ALL CHECKS PASSED\n');
    else
        fprintf('%d CHECK(S) FAILED\n', nFail);
    end
    out = struct('nFail', nFail, 'worst', worst);
end

function R = quat(q)
    w = q(1); x = q(2); y = q(3); z = q(4);
    R = [1-2*(y^2+z^2), 2*(x*y-z*w),   2*(x*z+y*w);
         2*(x*y+z*w),   1-2*(x^2+z^2), 2*(y*z-x*w);
         2*(x*z-y*w),   2*(y*z+x*w),   1-2*(x^2+y^2)];
end

function s = okStr(ok)
    if ok; s = 'ok'; else; s = 'FAIL'; end
end
