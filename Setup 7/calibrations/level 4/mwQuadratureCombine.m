function c = mwQuadratureCombine(qa, qb)
%MWQUADRATURECOMBINE Is there a channel skew? If not, the both-on balance.  (level 4, 4.2b)
%
%   qc = mwQuadratureCombine(q, q3)     q: calibMwQuadrature at quadAngle 0.5 (pi/2)
%                                       q3: the same at quadAngle 1.5 (3pi/2), either order
%
% No hardware. A delay between the two channels' pulse envelopes at the NVs reads as an
% imbalance in calibMwQuadrature: one channel plays alone at each pulse edge. With the
% quadrature pulse at 3pi/2 instead of pi/2, a real imbalance reads the same and the skew
% bias changes sign. Simulated (testMwQuadrature's ensemble, 6.7 MHz per channel):
%
%   no skew, pulse times within +/-10%:   d1 and d3 agree within 0.07 dB
%   2 ns skew:                            d1 - d3 = 1.9 to 2.7 dB
%
% So the DIFFERENCE is a robust skew test. The MEAN is not a robust skew correction: how
% far the 3pi/2 bias flips depends on the pulse length, and at 2 ns the mean was still off
% by -0.3 to +0.47 dB for pulses 10% short or long. Hence:
%
%   |d1 - d3| > 3 sigma + skewTolDb   'skew'   a timing offset between the channels is
%                                              present; no balance is given (deltaBoth NaN).
%                                              Its size is reported, roughly, in ns.
%   otherwise                          the balance is the mean, delta + (d1 + d3)/2, with
%                                      calibMwQuadrature's own rule: 'agree' when within
%                                      3 sigma + tolDb of the delta played, else 'disagree';
%                                      'inconclusive' when the error exceeds maxErrDb.
%
% Also 'inconclusive' when either run has no zero crossing. The ns figure assumes ~0.082 dB
% of skew bias per ns per MHz of per-channel drive in EACH run, so |d1 - d3| / (2 * 0.082 *
% Omega); it is +/-20% at best, and an overestimate when channel 2 compresses. Which channel
% is late is not reported: the sign depends on the AWG's phase convention.

    if ~(isstruct(qa) && isstruct(qb) && isfield(qa, 'fit') && isfield(qb, 'fit'))
        error('mwQuadratureCombine: give two calibMwQuadrature outputs.');
    end
    angA = angleOf(qa);  angB = angleOf(qb);
    if abs(angA - angB) < 1e-9
        error(['mwQuadratureCombine: both runs are at quadAngle %g. One must be 0.5 (pi/2) and ' ...
               'the other 1.5 (3pi/2).'], angA);
    end
    if angA < angB; q1 = qa; q3 = qb; else; q1 = qb; q3 = qa; end

    % Same pair, same point: otherwise the two crossings measure different things.
    same = @(x, y, tol) abs(x - y) <= tol;
    checks = { ...
        'amplitude', same(q1.amplitude, q3.amplitude, 0.01); ...
        'delta',     same(q1.delta, q3.delta, 0.001); ...
        'phiMax',    same(q1.phiMax, q3.phiMax, 0.01); ...
        'frequency', same(q1.frequency, q3.frequency, 1e-6); ...
        'offsets',   isequal(q1.fit.offsets, q3.fit.offsets)};
    bad = checks(~[checks{:, 2}], 1);
    if ~isempty(bad)
        error(['mwQuadratureCombine: the two runs differ in %s -- run both at the same ' ...
               'amplitude, delta, phase, frequency and offsets.'], strjoin(bad', ', '));
    end

    tolDb = 0.1;  maxErrDb = 0.3;  skewTolDb = 0.15;
    if isfield(q1, 'opts') && isfield(q1.opts, 'tolDb'); tolDb = q1.opts.tolDb; end
    if isfield(q1, 'opts') && isfield(q1.opts, 'maxErrDb'); maxErrDb = q1.opts.maxErrDb; end

    c = struct('verdict', 'inconclusive', 'why', '', 'delta', q1.delta, ...
               'deltaErr', getDeltaErr(q1), 'amplitude', q1.amplitude, 'phiMax', q1.phiMax, ...
               'dHalf', q1.dStar, 'dHalfErr', q1.dStarErr, ...
               'dThreeHalves', q3.dStar, 'dThreeHalvesErr', q3.dStarErr, ...
               'diffDb', NaN, 'diffErr', NaN, 'skewTolDb', skewTolDb, 'skewNsRough', NaN, ...
               'imbalanceDb', NaN, 'imbalanceErr', NaN, 'deltaBoth', NaN, 'deltaBothErr', NaN, ...
               'halfWhen', q1.when, 'threeHalvesWhen', q3.when);

    if ~(q1.fit.ok && q3.fit.ok && all(isfinite([q1.dStar q3.dStar q1.dStarErr q3.dStarErr])))
        c.why = 'one of the two runs has no zero crossing (see its own report)';
        report(c);
        return
    end
    c.diffDb  = q1.dStar - q3.dStar;
    c.diffErr = hypot(q1.dStarErr, q3.dStarErr);
    % Per-channel Rabi frequency from the PAIR's pi time: pi = 1/(2*(O1+O2)) = 1/(4*O).
    if isfield(q1, 'piTime') && isfinite(q1.piTime) && q1.piTime > 0
        omega = 1 / (4 * q1.piTime);
        c.skewNsRough = abs(c.diffDb) / (2 * 0.082 * omega);
    end

    if abs(c.diffDb) > 3 * c.diffErr + skewTolDb
        c.verdict = 'skew';
        c.why = 'the pi/2 and 3pi/2 crossings differ: a timing offset between the channels';
    else
        e = c.diffErr / 2;
        c.imbalanceDb = (q1.dStar + q3.dStar) / 2;   c.imbalanceErr = e;
        c.deltaBoth   = q1.delta + c.imbalanceDb;    c.deltaBothErr = e;
        if e > maxErrDb
            c.why = sprintf('combined error %.2f dB > maxErrDb %.2f', e, maxErrDb);
        elseif abs(c.imbalanceDb) <= 3 * hypot(e, c.deltaErr) + tolDb
            c.verdict = 'agree';
        else
            c.verdict = 'disagree';
        end
    end
    report(c);
end

% ======================================================================

function report(c)
    fprintf('\n============ BOTH CHANNELS ON: pi/2 and 3pi/2 together ============\n');
    fprintf('  zero crossings   pi/2 %+.3f +/- %.3f dB    3pi/2 %+.3f +/- %.3f dB\n', ...
            c.dHalf, c.dHalfErr, c.dThreeHalves, c.dThreeHalvesErr);
    if isfinite(c.diffDb)
        fprintf('  difference       %+.3f +/- %.3f dB  (skew test: %.2f dB allowed)', c.diffDb, ...
                c.diffErr, 3 * c.diffErr + c.skewTolDb);
        if isfinite(c.skewNsRough) && strcmp(c.verdict, 'skew')
            fprintf('  ~ %.1f ns between the channels (rough)', c.skewNsRough);
        end
        fprintf('\n');
    end
    if isfinite(c.imbalanceDb)
        fprintf('  balance with both on:  delta = %+.3f +/- %.3f dB   (single-channel: %+.3f +/- %.3f)\n', ...
                c.deltaBoth, c.deltaBothErr, c.delta, c.deltaErr);
    end
    switch c.verdict
        case 'agree'
            fprintf(['  AGREE: no channel skew, and the single-channel balance holds with both ' ...
                     'on. 4.3 may write it.\n']);
        case 'disagree'
            fprintf(['  DISAGREE by %.1f sigma, with no channel skew: the pair really is balanced ' ...
                     'differently\n  when both play. level4.m prints the next round.\n'], ...
                    abs(c.imbalanceDb) / max(hypot(c.imbalanceErr, c.deltaErr), eps));
        case 'skew'
            fprintf(['  SKEW: the two channels'' pulses are offset in time. This test then reads ' ...
                     'the skew, not the\n  balance -- no both-on balance is given. ' ...
                     'level4.m 4.3 says what that leaves.\n']);
        otherwise
            fprintf('  INCONCLUSIVE: %s\n', c.why);
    end
    fprintf('===================================================================\n');
end

function a = angleOf(q)
    a = 0.5;
    if isfield(q, 'quadAngle') && ~isempty(q.quadAngle); a = q.quadAngle; end
end

function e = getDeltaErr(q)
    e = 0;
    if isfield(q, 'deltaErr') && ~isempty(q.deltaErr) && isfinite(q.deltaErr); e = q.deltaErr; end
end
