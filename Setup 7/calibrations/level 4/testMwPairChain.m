function out = testMwPairChain(opts)
%TESTMWPAIRCHAIN Grade the calibMwPair arithmetic against a known two-channel pair.
%
%   out = testMwPairChain
%   out = testMwPairChain(struct('nTrial', 300))
%
% No instruments. Builds synthetic Rabi frequencies from a pair whose true amplitudes and
% true relative phase are known, pushes them through the same arithmetic calibMwPair uses,
% and checks that what comes back is what went in.
%
% This exists mostly for the SIGNS. The physics of the calibration is easy; the places it
% goes wrong are whether amplitudeDiff comes out positive when channel 2 is the weak one,
% whether dB is 20*log10 or 10*log10 of the thing being compared, and whether the fitted
% phase is the maximum or the null. None of those show up as a failure on the bench -- they
% show up as a setup that has been running with the channels 180 degrees apart for a month.
%
% TESTS
%   A  the dB round trip: apply the fitted amplitudeDiff and check the channels match
%   B  the phase fit recovers phi0, and phiMax is the MAXIMUM (not the null)
%   C  the closure ratio is 1 for a linear pair and detects a shared compression
%   D  the antiphase residual maps back to the same imbalance stage 1 measured
%   E  error bars on the recovered delta cover the truth at the rate they claim

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('omega1', 8.4, 'omega2', 7.6, 'phi0', 37, ...
                                     'relNoise', 0.004, 'nTrial', 200, 'seed', 7), ...
                        'testMwPairChain');

    O1 = opts.omega1; O2 = opts.omega2; phi0 = opts.phi0;
    trueDb = 20*log10(O1/O2);
    out = struct('pass', true, 'trueDb', trueDb, 'phi0', phi0);

    fprintf('\n============================================================\n');
    fprintf('testMwPairChain\n');
    fprintf('  truth: Omega1 %.3f, Omega2 %.3f MHz  ->  amplitudeDiff [0 %+.4f] dB\n', ...
            O1, O2, trueDb);
    fprintf('         relative phase of maximum drive %.1f deg\n', phi0);
    fprintf('         Omega1+Omega2 %.3f,  |Omega1-Omega2| %.3f MHz\n', O1+O2, abs(O1-O2));
    fprintf('============================================================\n');
    randn('state', opts.seed); rand('state', opts.seed);                    %#ok<RAND>

    % ---------------- TEST A: the dB round trip ----------------
    fprintf('\n--- TEST A: does the fitted dB offset actually equalise them? ---\n');
    m1 = O1*(1 + opts.relNoise*randn); m2 = O2*(1 + opts.relNoise*randn);
    delta = 20*log10(m1/m2);
    % Applying delta dB to channel 2 multiplies its AMPLITUDE, and therefore its Rabi
    % frequency, by 10^(delta/20). If the sign or the factor of 20 is wrong this diverges
    % instead of converging, which is the whole point of checking it.
    O2after = O2 * 10^(delta/20);
    a = struct();
    a.sign    = grade('delta > 0 when channel 2 is weaker', delta > 0, O2 < O1);
    a.value   = gradeNear('delta', delta, trueDb, 0.08, 'dB');
    a.equalise = gradeNear('Omega2 after applying it', O2after, O1, 0.06, 'MHz');
    fprintf('    channel 2 went from %.4f to %.4f MHz against channel 1 at %.4f\n', ...
            O2, O2after, O1);
    out.A = a; out.pass = out.pass && allOk(a);

    % ---------------- TEST B: the phase, and which extremum ----------------
    fprintf('\n--- TEST B: the phase fit ---\n');
    pd = 0:45:315;
    w  = pairOmega(O1, O2, pd, phi0);
    e  = opts.relNoise*w;
    pf = phaseScanFit(pd, w + e.*randn(size(w)), e, struct('verbose', false));
    b = struct();
    b.ok      = grade('fit ok', pf.ok, true);
    b.phi     = gradeNear('phiMax', angDiff(pf.phiMax, phi0), 0, 1.0, 'deg');
    % phiMax must be the MAXIMUM. Confusing it with the null is a silent, total failure --
    % the setup would run with the two channels cancelling -- so it is checked explicitly
    % against the model rather than trusted from the formula.
    wAtMax  = pairOmega(O1, O2, pf.phiMax,  phi0);
    wAtNull = pairOmega(O1, O2, pf.phiNull, phi0);
    b.isMax   = grade('phiMax really is the maximum', wAtMax > wAtNull, true);
    b.nullOpp = gradeNear('phiNull is phiMax + 180', angDiff(pf.phiNull, pf.phiMax+180), 0, 0.01, 'deg');
    b.sum     = gradeNear('Omega1+Omega2 from the curve', pf.sumOmega, O1+O2, 0.06, 'MHz');
    b.imb     = gradeNear('imbalance from the curve', pf.imbalanceDb, trueDb, 0.15, 'dB');
    fprintf('    Omega at phiMax %.3f, at phiNull %.3f MHz\n', wAtMax, wAtNull);
    out.B = b; out.pass = out.pass && allOk(b);

    % ---------------- TEST C: closure, linear and compressed ----------------
    fprintf('\n--- TEST C: the closure ratio R ---\n');
    S = O1 + O2;
    Rlin = pairOmega(O1, O2, phi0, phi0)/S;
    % A compression both channels share multiplies the summed field by a common factor. It
    % must show up in R and must NOT change the balance -- that is the claim calibMwPair
    % makes in its header, so it is worth testing rather than asserting.
    comp = 0.90;
    Rcomp = comp*pairOmega(O1, O2, phi0, phi0)/S;
    dbUnderComp = 20*log10((comp*O1)/(comp*O2));
    c = struct();
    c.linear   = gradeNear('R, linear pair',        Rlin,  1.00, 0.005, '');
    c.comp     = gradeNear('R, 10 % compression',   Rcomp, 0.90, 0.005, '');
    c.balSafe  = gradeNear('imbalance under compression', dbUnderComp, trueDb, 1e-9, 'dB');
    fprintf(['    a shared compression changes R but leaves the dB offset identical, which is\n' ...
             '    why it costs drive without invalidating amplitudeDiff\n']);
    % A phase error also lowers R, so the two causes must be distinguishable. They are: the
    % phase error is found and removed by TEST B, so a shortfall AT THE FITTED MAXIMUM
    % cannot be phase.
    dphi = 25;
    Rphase = pairOmega(O1, O2, phi0 + dphi, phi0)/S;
    fprintf('    for comparison, a %d deg phase error alone gives R = %.4f\n', dphi, Rphase);
    c.phaseCost = grade('a 25 deg phase error is a small R shortfall', Rphase > 0.94, true);
    out.C = c; out.pass = out.pass && allOk(c);

    % ---------------- TEST D: antiphase residual -> imbalance ----------------
    fprintf('\n--- TEST D: the antiphase residual maps back to the imbalance ---\n');
    resid = pairOmega(O1, O2, phi0 + 180, phi0);
    impliedDb = 20*log10((S + resid)/(S - resid));
    d = struct();
    d.resid   = gradeNear('residual at antiphase', resid, abs(O1-O2), 1e-9, 'MHz');
    d.implied = gradeNear('imbalance implied by it', impliedDb, trueDb, 1e-6, 'dB');
    % And the floor a phase error puts under that null, which is the reason stage 4 cannot
    % be trusted below a certain imbalance.
    for pe = [1 3 10]
        floorOmega = pairOmega(O1, O1, 180 + pe, 0);      % equal channels, phase off by pe
        fprintf('    with equal channels and the phase %2d deg off the null, residual %.4f MHz\n', ...
                pe, floorOmega);
    end
    fprintf(['    so a null measurement claims an imbalance that is really phase error unless\n' ...
             '    the phase is known well -- which is why stage 4 prints that floor alongside.\n']);
    out.D = d; out.pass = out.pass && allOk(d);

    % ---------------- TEST E: are the delta error bars honest? ----------------
    fprintf('\n--- TEST E: coverage of the delta error bar, %d trials ---\n', opts.nTrial);
    nIn = 0; devs = [];
    for k = 1:opts.nTrial
        s1 = O1*(1 + opts.relNoise*randn); s2 = O2*(1 + opts.relNoise*randn);
        dk = 20*log10(s1/s2);
        % Same propagation as calibMwPair: d(20 log10(x))/dx = 20/(ln10 * x), so the
        % fractional errors on the two frequencies add in quadrature and scale by 20/ln10.
        ek = (20/log(10))*hypot(opts.relNoise, opts.relNoise);
        devs = [devs dk - trueDb];                                          %#ok<AGROW>
        if abs(dk - trueDb) <= ek; nIn = nIn + 1; end
    end
    ek = (20/log(10))*hypot(opts.relNoise, opts.relNoise);
    e = struct();
    fprintf('  scatter of delta   %.4f dB\n', std(devs));
    fprintf('  propagated error   %.4f dB\n', ek);
    fprintf('  inside 1 sigma     %.1f %%  (want about 68)\n', 100*nIn/opts.nTrial);
    e.scatter = gradeNear('scatter vs propagated', std(devs)/ek, 1.0, 0.2, '');
    e.cover   = gradeNear('1-sigma coverage', 100*nIn/opts.nTrial, 68, 10, '%');
    e.bias    = gradeNear('bias', mean(devs), 0, 0.02, 'dB');
    fprintf(['  At %.1f %% Rabi errors the balance lands at %.3f dB, which is the number worth\n' ...
             '  remembering: one-at-a-time with Rabi is a hundredth-of-a-dB measurement, and no\n' ...
             '  null on a decaying oscillation gets close to that.\n'], 100*opts.relNoise, ek);
    out.E = e; out.pass = out.pass && allOk(e);

    fprintf('\n============================================================\n');
    if out.pass
        fprintf('ALL TESTS PASS.\n');
    else
        fprintf('SOMETHING FAILED -- read the FAIL lines above.\n');
    end
    fprintf('============================================================\n');
end

% ======================================================================

function w = pairOmega(O1, O2, phaseDeg, phi0)
% The exact two-drive interference. This is the model, and the test's only physics.
    w = sqrt(O1^2 + O2^2 + 2*O1*O2*cosd(phaseDeg - phi0));
end

function d = angDiff(a, b)
% Signed difference wrapped into (-180, 180], so 359 vs 1 is 2 degrees and not 358.
    d = mod(a - b + 180, 360) - 180;
end

function ok = grade(name, got, want)
    ok = isequal(logical(got), logical(want));
    fprintf('  %-38s got %d, want %d   %s\n', name, logical(got), logical(want), tf(ok));
end

function ok = gradeNear(name, got, want, tol, unit)
    ok = isfinite(got) && abs(got - want) <= tol;
    fprintf('  %-38s %.4f vs %.4f %s (tol %g)   %s\n', name, got, want, unit, tol, tf(ok));
end

function s = tf(v)
    if v; s = 'ok'; else; s = 'FAIL'; end
end

function ok = allOk(s)
    ok = true; f = fieldnames(s);
    for i = 1:numel(f)
        v = s.(f{i});
        if islogical(v) || (isnumeric(v) && isscalar(v)); ok = ok && logical(v); end
    end
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
