function out = testB1RabiChain(opts)
%TESTB1RABICHAIN Grade the calibB1Rabi maths against a resonator whose answer is known.
%
%   out = testB1RabiChain
%   out = testB1RabiChain(struct('nTrial', 200))
%
% No instruments. This builds synthetic scans from a known resonator and a known coil, runs
% them through the two pieces that turn measurements into an answer -- resonatorLorentzianFit
% and coilLineRelation -- and checks that the answer comes back. It exists because the two
% previous calibrations in this series each shipped with a bug that only a known-truth test
% could have found: the multi-line fitter could not split a blend because its seeded
% amplitudes came from the wrong envelope, and calibFieldSurvey threw away correctly aligned
% pairs because the mirror-partner test assumed a single D that does not exist. Both were
% invisible on real data, where there is nothing to compare against.
%
% WHAT EACH TEST IS FOR
%
%   A  a clean wide scan. Does the chain recover f0, Q, the on-resonance Rabi frequency and
%      the optimum coil setting at all? If A fails nothing else matters.
%   B  scan width swept from too-narrow to generous. Establishes WHERE the honest boundary
%      is: at which span the Lorentzian stops being determined, that the code says so, and
%      that the parabola fallback covers the gap. A calibration that silently degrades is
%      worse than one that refuses.
%   C  a scan sitting entirely on one wing. The peak is not in the data. The only correct
%      behaviour is to refuse, and the fit must not return a number dressed up with an
%      error bar.
%   D  a rogue point tracking a different NV orientation, the exact failure that produced
%      the 0.13*Brho + 9.7 coil fit on a live run. coilLineRelation has to notice from the
%      slope alone, because that is all it is given.
%   E  many noise realisations of the same scan. Are the quoted error bars honest -- does
%      the truth fall inside one sigma about 68% of the time? An error bar that is wrong by
%      3x is what turns a settled calibration into an argument.
%
% Every test prints what it got, what it wanted, and ok or FAIL.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('f0', 2823.0, 'Q', 94, 'wMax', 20.0, 'floorMHz', 1.5, ...
                                     'coilK', 0.974, 'coilM', -0.813, 'D0', 2870.5, ...
                                     'relNoise', 0.02, 'lineNoise', 0.15, ...
                                     'nTrial', 120, 'seed', 20260814), 'testB1RabiChain');

    GAMMA = 2.8025;
    W  = opts.f0/opts.Q;
    cT = opts.floorMHz^2;
    AT = opts.wMax^2 - cT;
    truth = struct('f0', opts.f0, 'fwhm', W, 'Q', opts.Q, 'wMax', opts.wMax, ...
                   'A', AT, 'c', cT, 'coilK', opts.coilK, 'coilM', opts.coilM);
    out = struct('truth', truth, 'pass', true);

    % The field that puts the line exactly on resonance, worked out from the truth so the
    % test knows the answer it is grading against.
    %   f = D0 - gamma*(k*rho + m)  ->  rho = (D0 - m*gamma - f)/(gamma*k)
    rhoTrue = (opts.D0 - GAMMA*opts.coilM - opts.f0)/(GAMMA*opts.coilK);
    truth.rhoOptimum = rhoTrue;
    out.truth = truth;

    fprintf('\n============================================================\n');
    fprintf('testB1RabiChain\n');
    fprintf('  truth: f0 %.2f MHz, FWHM %.2f, Q %.0f, Rabi at peak %.2f MHz\n', ...
            opts.f0, W, opts.Q, opts.wMax);
    fprintf('         coil |B| = %.3f*Brho %+.3f  ->  optimum Brho %.3f G\n', ...
            opts.coilK, opts.coilM, rhoTrue);
    fprintf('============================================================\n');

    rand('state', opts.seed); randn('state', opts.seed);                    %#ok<RAND>

    % ---------------- TEST A ----------------
    fprintf('\n--- TEST A: clean wide scan, 13 fields over 10 G ---\n');
    S = synth(linspace(rhoTrue - 5, rhoTrue + 5, 13), opts, truth, GAMMA, true);
    fit = resonatorLorentzianFit(S.fLine, S.fRabi, S.fRabiErr, struct('verbose', false));
    rel = coilLineRelation(S.rho, S.fLine, S.fLineErr, struct('verbose', false, 'gamma', GAMMA));
    [rhoStar, rhoErr] = rel.invert(fit.f0, fit.f0Err);
    a = struct();
    a.usable  = grade('fit usable',        fit.usable, true);
    a.f0      = gradeNear('f0',            fit.f0,   opts.f0,   1.0, 'MHz');
    a.Q       = gradeNear('Q',             fit.Q,    opts.Q,    0.25*opts.Q, '');
    a.wMax    = gradeNear('Rabi at peak',  fit.wMax, opts.wMax, 0.5, 'MHz');
    a.slope   = grade('coil slope sane',   rel.slopeSane, true);
    a.rho     = gradeNear('optimum Brho',  rhoStar,  rhoTrue,   0.25, 'G');
    fprintf('    (quoted errors: f0 %.2f, Brho %.3f)\n', fit.f0Err, rhoErr);
    out.A = a; out.pass = out.pass && allOk(a);

    % ---------------- TEST B ----------------
    fprintf('\n--- TEST B: how narrow can the scan get before it stops working ---\n');
    fprintf('    A REFUSAL IS A PASS HERE. What must not happen is a confident wrong number.\n');
    fprintf('%8s %10s %10s %8s %8s %10s %10s  %s\n', 'span/MHz', 'lorentz f0', 'FWHM', ...
            'usable', 'parab f0', 'best avail', 'err/MHz', 'verdict');
    b = struct('rows', [], 'pass', true);
    for span = [8 12 18 25 35 50 70]
        rhoSpan = span/(GAMMA*opts.coilK);
        S = synth(linspace(rhoTrue - rhoSpan/2, rhoTrue + rhoSpan/2, 11), opts, truth, GAMMA, true);
        fit = resonatorLorentzianFit(S.fLine, S.fRabi, S.fRabiErr, struct('verbose', false));
        % What the calibration would actually hand over, by the same rule calibB1Rabi uses.
        if fit.usable
            fBest = fit.f0; eBest = fit.f0Err; src = 'lorentz';
        elseif fit.parabola.ok && fit.parabola.inRange
            fBest = fit.parabola.f0; eBest = fit.parabola.f0Err; src = 'parabola';
        else
            fBest = NaN; eBest = NaN; src = 'refused';
        end
        % Honest means: either it refused, or the number it gave is right.
        good = strcmp(src, 'refused') || abs(fBest - opts.f0) < max(3*eBest, 1.0);
        if ~good; b.pass = false; end
        fprintf('%8.0f %10.2f %10.2f %8d %8.2f %10.2f %10.2f  %s %s\n', ...
                span, fit.f0, fit.fwhm, fit.usable, fit.parabola.f0, fBest, eBest, ...
                src, BooleanHelper2(good, 'ok', 'FAIL'));
        b.rows = [b.rows; span, fit.f0, fit.fwhm, double(fit.usable), fBest, eBest, double(good)];
    end
    fprintf('  %s\n', BooleanHelper2(b.pass, ...
        'ok -- every span either refused or gave a right answer', ...
        'FAIL -- some span gave a confident wrong centre'));
    out.B = b; out.pass = out.pass && b.pass;

    % ---------------- TEST C ----------------
    fprintf('\n--- TEST C: scan entirely on the low wing, peak not in the data ---\n');
    rhoLo = (opts.D0 - GAMMA*opts.coilM - (opts.f0 - 45))/(GAMMA*opts.coilK);
    S = synth(linspace(rhoLo - 3, rhoLo + 3, 9), opts, truth, GAMMA, true);
    fit = resonatorLorentzianFit(S.fLine, S.fRabi, S.fRabiErr, struct('verbose', false));
    c = struct();
    c.notUsable = grade('fit refuses (usable == false)', fit.usable, false);
    handedOver = fit.usable || (fit.parabola.ok && fit.parabola.inRange);
    c.noHandover = grade('nothing handed over', handedOver, false);
    fprintf('    scan covered %.1f to %.1f MHz, true peak at %.1f -- %.0f MHz outside.\n', ...
            min(S.fLine), max(S.fLine), opts.f0, opts.f0 - max(S.fLine));
    out.C = c; out.pass = out.pass && allOk(c);

    % ---------------- TEST D ----------------
    fprintf('\n--- TEST D: one point tracking a different NV orientation ---\n');
    S = synth(linspace(rhoTrue - 5, rhoTrue + 5, 13), opts, truth, GAMMA, true);
    Sbad = S;
    % A different orientation's lower line sits well above the aligned one, and its
    % position moves with a smaller projection -- so one substituted point both jumps and
    % breaks the slope. This is what the live failure looked like.
    Sbad.fLine(4) = Sbad.fLine(4) + 28;
    relClean = coilLineRelation(S.rho,    S.fLine,    S.fLineErr,    struct('verbose', false, 'gamma', GAMMA));
    relBad   = coilLineRelation(Sbad.rho, Sbad.fLine, Sbad.fLineErr, struct('verbose', false, 'gamma', GAMMA));
    d = struct();
    d.cleanSane = grade('clean scan: slope sane',    relClean.slopeSane, true);
    d.badLoud   = grade('rogue point: chi2 > 10',    relBad.chi2Reduced > 10, true);
    fprintf('    clean: slope %.3f, chi2 %.1f.  With the rogue point: slope %.3f, chi2 %.1f.\n', ...
            relClean.a, relClean.chi2Reduced, relBad.a, relBad.chi2Reduced);
    fprintf(['    NOTE the rogue slope is still "sane" (%d) -- one bad point in thirteen bends\n' ...
             '    the line without breaking it, so the SLOPE test alone does not catch this and\n' ...
             '    reduced chi2 is what does. Both are printed by coilLineRelation for that\n' ...
             '    reason; the upstream orientation gate is the real defence.\n'], relBad.slopeSane);
    out.D = d; out.pass = out.pass && allOk(d);

    % ---------------- TEST E ----------------
    fprintf('\n--- TEST E: are the error bars honest? %d noise realisations ---\n', opts.nTrial);
    nIn1 = 0; nIn2 = 0; nOk = 0; devs = []; pulls = [];
    for k = 1:opts.nTrial
        S = synth(linspace(rhoTrue - 5, rhoTrue + 5, 13), opts, truth, GAMMA, true);
        fit = resonatorLorentzianFit(S.fLine, S.fRabi, S.fRabiErr, struct('verbose', false));
        if ~fit.usable || ~isfinite(fit.f0Err) || fit.f0Err <= 0; continue; end
        nOk = nOk + 1;
        dev = fit.f0 - opts.f0;
        devs = [devs dev];                                                  %#ok<AGROW>
        pulls = [pulls dev/fit.f0Err];                                      %#ok<AGROW>
        if abs(dev) <= fit.f0Err;   nIn1 = nIn1 + 1; end
        if abs(dev) <= 2*fit.f0Err; nIn2 = nIn2 + 1; end
    end
    e = struct();
    if nOk < 20
        fprintf('  only %d of %d trials fitted -- not enough to judge the error bars.\n', nOk, opts.nTrial);
        e.enough = grade('enough successful trials', false, true);
    else
        f1 = nIn1/nOk; f2 = nIn2/nOk;
        pullSd = std(pulls);
        fprintf('  scatter of f0        %.3f MHz\n', std(devs));
        fprintf('  mean quoted error    %.3f MHz\n', mean(abs(devs))/max(mean(abs(pulls)), eps));
        fprintf('  bias                 %+.3f MHz  (a systematic offset, not noise)\n', mean(devs));
        fprintf('  inside 1 sigma       %.1f %%   (want about 68)\n', 100*f1);
        fprintf('  inside 2 sigma       %.1f %%   (want about 95)\n', 100*f2);
        fprintf('  sd of the pull       %.2f     (want about 1: >1 errors too small)\n', pullSd);
        e.enough  = grade('enough successful trials', true, true);
        e.oneSig  = gradeNear('1-sigma coverage', 100*f1, 68, 15, '%');
        e.twoSig  = gradeNear('2-sigma coverage', 100*f2, 95, 8,  '%');
        e.pull    = gradeNear('pull width',       pullSd, 1.0, 0.5, '');
        e.bias    = gradeNear('bias',             mean(devs), 0, 0.3, 'MHz');
    end
    out.E = e; out.pass = out.pass && allOk(e);

    fprintf('\n============================================================\n');
    if out.pass
        fprintf('ALL TESTS PASS. The maths behind calibB1Rabi recovers a known resonator,\n');
        fprintf('refuses when the peak is not in the data, and quotes error bars that cover\n');
        fprintf('the truth as often as they claim to.\n');
    else
        fprintf('SOMETHING FAILED. Read the FAIL lines above -- each one names the quantity,\n');
        fprintf('what it should have been and what it was.\n');
    end
    fprintf('============================================================\n');
end

% ======================================================================

function S = synth(rhos, opts, truth, GAMMA, withNoise)
% One synthetic scan: fields in, measured line positions and Rabi frequencies out.
%
% Built from the SAME physics the fit assumes, which is the point -- this tests the
% estimator, not the model. A test that also changed the model would not be able to say
% which of the two failed.
    rhos = rhos(:)';
    Btrue = truth.coilK*rhos + truth.coilM;
    fLine = opts.D0 - GAMMA*Btrue;

    w2 = truth.A ./ (1 + 4*(fLine - truth.f0).^2 / truth.fwhm^2) + truth.c;
    fRabi = sqrt(w2);

    fRabiErr = opts.relNoise * fRabi;       % a Rabi fit's error scales with the frequency
    fLineErr = opts.lineNoise * ones(size(fLine));
    if withNoise
        fRabi = fRabi + fRabiErr .* randn(size(fRabi));
        fLine = fLine + fLineErr .* randn(size(fLine));
    end
    S = struct('rho', rhos, 'fLine', fLine, 'fLineErr', fLineErr, ...
               'fRabi', fRabi, 'fRabiErr', fRabiErr, 'Btrue', Btrue);
end

function ok = grade(name, got, want)
    ok = isequal(logical(got), logical(want));
    fprintf('  %-32s got %d, want %d   %s\n', name, logical(got), logical(want), ...
            BooleanHelper2(ok, 'ok', 'FAIL'));
end

function ok = gradeNear(name, got, want, tol, unit)
    ok = isfinite(got) && abs(got - want) <= tol;
    fprintf('  %-32s %.3f vs %.3f %s (tol %.3f)   %s\n', name, got, want, unit, tol, ...
            BooleanHelper2(ok, 'ok', 'FAIL'));
end

function tf = allOk(s)
    tf = true;
    f = fieldnames(s);
    for i = 1:numel(f)
        v = s.(f{i});
        if islogical(v) || (isnumeric(v) && isscalar(v)); tf = tf && logical(v); end
    end
end

function s = BooleanHelper2(tf, a, b)
% Local, so this test does not depend on the lab's BooleanHelper being on the path.
    if tf; s = a; else; s = b; end
end

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d);
    given = fieldnames(o);
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
