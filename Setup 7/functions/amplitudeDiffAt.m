function [amplitudeDiff, omega, piNs, info] = amplitudeDiffAt(a, opts)
%AMPLITUDEDIFFAT The channel-2 offset that balances the pair at a given drive.
%
%   d = amplitudeDiffAt(-16)                 % [0 delta], ready for the json or b.amplitudeDiff
%   [d, w, piNs] = amplitudeDiffAt(-16)      % and the Rabi frequency and pi time to expect
%   d = amplitudeDiffAt(-16, struct('curve', c))   % from a curve in hand, not the saved one
%
% Reads amplitudeDiffCurve.mat from the current results folder (setup7Paths) and solves for
% the offset at whatever amplitude you ask for. calibAmplitudeDiff measures the curve.
%
% ===================================================================
% WHAT IT SOLVES, AND WHY NOT JUST INTERPOLATE THE MEASURED DELTAS
% ===================================================================
%
% amplitudeDiff is defined by
%
%       Omega1(a)  =  Omega2(a + delta)
%
% -- channel 2 at a DIFFERENT amplitude, producing the field channel 1 produces at a. So
% this takes channel 1's curve, reads off the field it delivers at the requested drive, and
% asks channel 2's curve what drive it needs to match it. Two interpolations and an
% inversion, not an interpolation of deltas.
%
% Interpolating a table of per-point deltas would be interpolating the VERTICAL gap between
% the curves, which equals the horizontal offset only where both curves have unit slope.
% Under compression, at channel-2 slope s2, the two differ by a factor 1/s2: a 0.9 dB gap
% at s2 = 0.7 means 1.29 dB of offset. calibAmplitudeDiff's plot draws both so the size of
% the difference on your hardware is visible.
%
% ===================================================================
% SHAPE-PRESERVING INTERPOLATION, AND NO EXTRAPOLATION AT ALL
% ===================================================================
%
% pchip, not spline. A cubic spline through ten measured points will overshoot between
% them, and an overshoot in a gain curve is a non-monotone patch -- which makes the
% inversion ambiguous, and does it silently. pchip cannot overshoot, so a curve that was
% measured monotone stays monotone.
%
% Outside the measured amplitudes this REFUSES rather than extrapolating. That is not
% caution for its own sake: the entire reason the curve exists is compression, and
% compression is precisely the region where a polynomial fitted below the knee predicts
% the wrong answer with a confident-looking number. If you need -5 dBm, measure to -5 dBm.
%
% One more refusal worth knowing about. If channel 1 at your requested drive delivers more
% field than channel 2 reached ANYWHERE in the scan, there is no offset that balances them
% -- channel 2 saturates first. That is a real hardware finding rather than an interpolation
% failure, and it is reported as such.
%
% OPTIONS
%   curve       a curve struct to use instead of the saved one (default [] = load it)
%   file        an explicit .mat to load (default [] = results folder via setup7Paths)
%   verbose     print what it did (default true)
%   quietRange  do not warn when the drive is above the linearity ceiling (default false)
%
% OUTPUT
%   amplitudeDiff  [0 delta], the vector the json and b.amplitudeDiff both want
%   omega          the Rabi frequency per channel expected at that drive, MHz
%   piNs           the pi time expected with BOTH channels, ns. Uses the measured R from
%                  the curve where it is available, so it is what you will actually get
%                  rather than the ideal both-channels-add value.
%   info           .delta .slope1 .slope2 .withinLinear .curveDate .range

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('curve', [], 'file', [], 'verbose', true, ...
                                     'quietRange', false), 'amplitudeDiffAt');

    if ~isnumeric(a) || ~isscalar(a) || ~isfinite(a)
        error('amplitudeDiffAt: give one amplitude in dBm, e.g. amplitudeDiffAt(-16).');
    end

    % ---- the curve ------------------------------------------------------
    c = opts.curve;
    if isempty(c)
        f = opts.file;
        if isempty(f)
            P = setup7Paths();
            f = fullfile(P.folder, 'amplitudeDiffCurve.mat');
        end
        if ~isfile(f)
            error(['amplitudeDiffAt: no amplitudeDiffCurve.mat at\n    %s\n' ...
                   'Run calibAmplitudeDiff to measure it. Until then the single value in ' ...
                   'expProperties.json is all there is, and it is only right at the ' ...
                   'amplitude calibMwPair used.'], f);
        end
        S = load(f);
        if ~isfield(S, 'curve')
            error('amplitudeDiffAt: %s does not contain a "curve" variable.', f);
        end
        c = S.curve;
    end
    for need = {'amplitudes','omega1','omega2'}
        if ~isfield(c, need{1})
            error('amplitudeDiffAt: the curve has no "%s" field.', need{1});
        end
    end

    A  = c.amplitudes(:)';
    G1 = 20*log10(c.omega1(:)');
    G2 = 20*log10(c.omega2(:)');
    keep = isfinite(A) & isfinite(G1) & isfinite(G2);
    A = A(keep); G1 = G1(keep); G2 = G2(keep);
    [A, ord] = sort(A); G1 = G1(ord); G2 = G2(ord);
    if numel(A) < 3
        error('amplitudeDiffAt: the curve has only %d usable points.', numel(A));
    end

    info = struct('delta', NaN, 'slope1', NaN, 'slope2', NaN, 'withinLinear', true, ...
                  'curveDate', '', 'range', [min(A) max(A)]);
    if isfield(c, 'measuredAt'); info.curveDate = c.measuredAt; end

    % ---- in range? ------------------------------------------------------
    if a < min(A) || a > max(A)
        error(['amplitudeDiffAt: %.2f dBm is outside the measured range %.2f to %.2f dBm, ' ...
               'and this will not extrapolate a compression curve -- below the knee a fit ' ...
               'predicts the wrong answer confidently. Re-run calibAmplitudeDiff over the ' ...
               'range you need.'], a, min(A), max(A));
    end

    % ---- channel 2 must be usable as an inverse -------------------------
    if ~all(diff(G2) > 0)
        bad = find(diff(G2) <= 0);
        error(['amplitudeDiffAt: channel 2''s gain curve is not monotone -- it does not ' ...
               'increase between %.2f and %.2f dBm. An inverse is ambiguous there. Either ' ...
               'those points are noise (re-measure them) or the amplifier is doing ' ...
               'something worth understanding before this number is trusted.'], ...
              A(bad(1)), A(bad(1)+1));
    end

    % ---- the solve ------------------------------------------------------
    target = interp1(A, G1, a, 'pchip');
    if target > max(G2)
        error(['amplitudeDiffAt: at %.2f dBm channel 1 delivers %.2f dB, and channel 2 did ' ...
               'not reach that anywhere up to %.2f dBm (its best was %.2f dB). No offset ' ...
               'balances them here -- channel 2 saturates first. That is a hardware limit, ' ...
               'not a fitting problem.'], a, target, max(A), max(G2));
    end
    if target < min(G2)
        error(['amplitudeDiffAt: at %.2f dBm channel 1 delivers %.2f dB, below anything ' ...
               'channel 2 was measured at (its lowest was %.2f dB at %.2f dBm). Extend the ' ...
               'scan downward.'], a, target, min(G2), min(A));
    end
    a2 = interp1(G2, A, target, 'pchip');
    info.delta = a2 - a;

    % local slopes, for the report and because they are the reason this file exists
    h = max(0.25, 0.02*(max(A)-min(A)));
    info.slope1 = numSlope(@(x) interp1(A, G1, x, 'pchip'), a,  h, [min(A) max(A)]);
    info.slope2 = numSlope(@(x) interp1(A, G2, x, 'pchip'), a2, h, [min(A) max(A)]);

    amplitudeDiff = [0, info.delta];
    omega = 10^(target/20);
    % Balanced and IN PHASE the pair delivers 2*omega, so pi is 500/(2*omega) in ns. But
    % measured R is 0.90 to 0.97 on this setup rather than 1, so quoting the ideal would
    % over-promise the pi time by 3 to 10 %. Where the curve carries R, it is used.
    info.R = NaN;
    Reff = 1;
    if isfield(c, 'R') && any(isfinite(c.R))
        rr = c.R(:)'; aa = c.amplitudes(:)';
        k = isfinite(rr) & isfinite(aa);
        if sum(k) >= 2
            info.R = interp1(aa(k), rr(k), a, 'pchip', 'extrap');
            Reff = max(min(info.R, 1), 0.1);
        end
    end
    info.omegaBoth = 2*omega*Reff;
    piNs = 500/max(info.omegaBoth, eps);

    % ---- is this drive inside the linear region? ------------------------
    if isfield(c, 'linearityCeiling') && isfinite(c.linearityCeiling)
        info.withinLinear = a <= c.linearityCeiling + 1e-9;
        if ~info.withinLinear && ~opts.quietRange
            warning('amplitudeDiffAt:aboveLinear', ...
                    ['%.2f dBm is above the %.2f dBm where the pair last added linearly ' ...
                     '(R < 1 beyond it). The balance below is still right; you simply get ' ...
                     'less total drive than Omega1 + Omega2 would suggest.'], ...
                    a, c.linearityCeiling);
        end
    end

    if opts.verbose
        fprintf('\namplitudeDiff at %.2f dBm:  [0 %+.3f]\n', a, info.delta);
        fprintf('  channel 2 wants %.2f dBm to match channel 1''s %.3f MHz\n', a2, omega);
        fprintf('  local slopes: ch1 %.3f, ch2 %.3f dB/dB\n', info.slope1, info.slope2);
        vgap = interp1(A, G1 - G2, a, 'pchip');
        fprintf('  the vertical gap here is %+.3f dB -- the offset is %+.3f, a %.3f dB\n', ...
                vgap, info.delta, abs(info.delta - vgap));
        fprintf('  difference, which is what taking the gap as the answer would have cost\n');
        if isfinite(info.R)
            fprintf('  expect %.3f MHz per channel; both give %.3f MHz (R = %.3f), pi %.1f ns\n', ...
                    omega, info.omegaBoth, info.R, piNs);
            fprintf('  the ideal R = 1 would be %.1f ns -- the shortfall is the phase, not the drive\n', ...
                    500/(2*omega));
        else
            fprintf('  expect %.3f MHz per channel, %.1f ns pi with both (assumes R = 1)\n', ...
                    omega, piNs);
        end
        if ~isempty(info.curveDate)
            fprintf('  curve measured %s\n', info.curveDate);
        end
    end
end

% ======================================================================

function s = numSlope(fn, x, h, lims)
% Central difference, pulled inside the measured range at the ends so the interpolant is
% never asked to extrapolate just to report a slope.
    lo = max(x-h, lims(1)); hi = min(x+h, lims(2));
    if hi <= lo; s = NaN; return; end
    s = (fn(hi) - fn(lo))/(hi - lo);
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
