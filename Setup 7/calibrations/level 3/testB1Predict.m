function out = testB1Predict(opts)
%TESTB1PREDICT Grade calibB1Rabi's predict mode on a simulated rig. No hardware.
%
%   testB1Predict
%   testB1Predict(struct('plotIt', true, 'compareMeasure', true))
%
% A real Helmholtz object on dummy supplies (helmholtzOfflineRig) carries the 2026-11-23
% coil calibration -- a separate one when the setup is running, so this is safe in the lab
% session; the ESR survey and the Rabi are simulated (calibB1Rabi's surveyFn and
% rabiFn hooks). The simulation has the two things the 2026-11-27 lab run showed:
%
%   - the coil calibration is slightly off: the TRUE |B| is 0.99*Brho - 0.10 G
%   - the resonator (f0 2763.75 MHz, 31.6 MHz FWHM) pulls a SINGLE-LINE fit of the lower
%     aligned line by up to 3 MHz when the line is on it; |B| and D from the whole spectrum
%     are not pulled (they come mostly from lines far from the resonator)
%
% and the Rabi frequency follows the TRUE line: Omega^2 = Lorentzian(f_line) + floor, plus
% the detuning of the drive from the true line, sqrt(Omega^2 + Delta^2).
%
% Checks: 3 surveys in all (2 anchors + verify) instead of ~23; the anchors recover the true
% |B| law; the resonator centre comes back; at the chosen field the TRUE line sits on f0;
% and the frequency handed to the json is the true line. compareMeasure also runs the old
% per-point ESR mode on the same rig for comparison (not graded).

    if ~exist('opts', 'var'); opts = struct(); end
    if ~isfield(opts, 'plotIt'); opts.plotIt = false; end
    if ~isfield(opts, 'compareMeasure'); opts.compareMeasure = false; end

    here = fileparts(mfilename('fullpath'));
    s7 = fileparts(fileparts(here));
    src = fullfile(s7, 'calibrations', 'results', '2026-11-23', 'expProperties.json');
    json = [tempname '.json'];
    copyfile(src, json);
    % Safe next to the live setup: alongside a live session h is a SECOND object on dummy
    % supplies (named 'Helmholtz_offlineRig'), handed to calibB1Rabi through its 'helmholtz'
    % option. The live object, its supplies and its currents are never touched.
    [h, rig] = helmholtzOfflineRig(json, struct('alongsideLive', true));   % delete(rig) closes h

    TRUE = struct('a', 0.99, 'b', -0.10, 'D', 2870.20, 'f0', 2763.75, 'W', 31.6, ...
                  'A', 219, 'floor', 18.7);
    rs = RandStream('mt19937ar', 'Seed', 3);
    nSurvey = 0; nRabi = 0;

    fprintf('\n==================================================================\n');
    fprintf('testB1Predict -- calibB1Rabi predict mode on a simulated rig\n');
    fprintf('==================================================================\n');

    base = struct('rhoCoarse', 34:1:46, 'settle', 0, 'saveRaw', false, ...
        'plotIt', opts.plotIt, 'jsonLocation', json, 'surveyFn', @fakeSurvey, ...
        'rabiFn', @fakeRabi, 'controlAuto', false, 'helmholtz', h);
    o = base; o.esrMode = 'predict';
    res = calibB1Rabi(o);

    nFail = 0;
    ok = strcmp(res.esrMode, 'predict') && nSurvey == 3;
    nFail = nFail + report(sprintf('predict mode, %d ESR surveys in all (2 anchors + verify)', nSurvey), ok);
    A = res.anchor;
    err = max(abs((A.a*[34 46] + A.b) - (TRUE.a*[34 46] + TRUE.b)));
    ok = A.n == 2 && err < 0.15;
    nFail = nFail + report(sprintf('anchors recover the true |B| law (worst %.3f G)', err), ok);
    ok = abs(res.resonantFrequency - TRUE.f0) < 0.5;
    nFail = nFail + report(sprintf('resonator centre %.2f vs %.2f MHz', res.resonantFrequency, TRUE.f0), ok);
    fTrue = trueLine();
    ok = abs(fTrue - TRUE.f0) < 0.5;
    nFail = nFail + report(sprintf('at the chosen field the TRUE line is on f0 (%+.2f MHz)', fTrue - TRUE.f0), ok);
    ok = isfield(res, 'driveFrequency') && abs(res.driveFrequency - fTrue) < 0.3;
    nFail = nFail + report(sprintf('drive frequency for the json is the true line (%+.2f MHz)', ...
                           res.driveFrequency - fTrue), ok);
    out.predict = res;

    if opts.compareMeasure
        nSurvey = 0;
        h.setBrho(38.19);
        o = base; o.esrMode = 'measure';
        resM = calibB1Rabi(o);
        fprintf(['\n  (measure mode on the same rig: %d surveys, f0 %.2f MHz, true line at its\n' ...
                 '   chosen field %+.2f MHz from f0 -- not graded)\n'], nSurvey, ...
                resM.resonantFrequency, trueLine() - TRUE.f0);
        out.measure = resM;
    end

    fprintf('\n------------------------------------------------------------------\n');
    if nFail == 0; fprintf('ALL CHECKS PASSED\n'); else; fprintf('%d CHECK(S) FAILED\n', nFail); end
    out.nFail = nFail;
    % The fake-hardware handles are nested functions: anything holding one keeps this whole
    % workspace -- and so the rig and its Helmholtz object -- alive. Strip them, then close.
    for f = fieldnames(out)'
        if isstruct(out.(f{1})) && isfield(out.(f{1}), 'opts')
            out.(f{1}).opts.surveyFn = []; out.(f{1}).opts.rabiFn = [];
        end
    end
    delete(rig);
    delete(json);

    % ---------------------------------------------------------------------------
    function f = trueLine()
        % the aligned lower line at the coils' present setting, with the TRUE |B|
        u = h.B / h.Brho;
        f = min(Helmholtz.nvLines((TRUE.a*h.Brho + TRUE.b) * u, TRUE.D));
    end

    function sv = fakeSurvey(freqs, ~)
        nSurvey = nSurvey + 1;
        u = h.B / h.Brho;
        Bt = (TRUE.a*h.Brho + TRUE.b) * u;
        L4 = Helmholtz.nvLinesByAxis(Bt, TRUE.D);         % [lower upper] per orientation
        [~, kAl] = max(abs(Helmholtz.NV_AXES * Bt(:)));
        lo = min(freqs) + 10; hi = max(freqs) - 10;
        lines = struct('centre', {}, 'centreErr', {}, 'fwhm', {}, 'depth', {});
        M = struct('centre', {}, 'predicted', {}, 'axis', {}, 'resid', {});
        for k = 1:4
            for s = 1:2
                f = L4(k, s);
                if f < lo || f > hi; continue; end
                c = f + 0.1*randn(rs);
                if k == kAl && s == 1                       % the resonator pulls this one
                    c = c + 3*sin((f - TRUE.f0)/6) * exp(-((f - TRUE.f0)/15)^2);
                end
                lines(end+1) = struct('centre', c, 'centreErr', 0.2, 'fwhm', 20, 'depth', 0.01); %#ok<AGROW>
                M(end+1) = struct('centre', c, 'predicted', f, 'axis', k, 'resid', c - f);      %#ok<AGROW>
            end
        end
        al = L4(kAl, :);
        visible = all(al >= lo & al <= hi);
        sv = struct();
        sv.fit = struct('lines', lines, 'ok', true);
        sv.solve = struct('ok', visible, 'B', norm(Bt) + 0.04*randn(rs), ...
            'D', TRUE.D + 0.05*randn(rs), 'tiltDeg', 0.2, 'outerVisible', visible, ...
            'alignedAxisIdx', kAl, 'alignedLines', al, 'matched', M);
        sv.B = sv.solve.B;
    end

    function r = fakeRabi(freq, ~)
        nRabi = nRabi + 1;
        fl = trueLine();
        w2 = TRUE.A ./ (1 + ((fl - TRUE.f0)/(TRUE.W/2)).^2) + TRUE.floor;
        w = sqrt(w2 + (freq - fl)^2) * (1 + 0.005*randn(rs));
        r = struct('ok', true, 'why', '', 'f', w, 'fErr', 0.01*w, 'contrast', 0.02, ...
                   'r2', 0.95, 'r2Ceiling', 0.97, 'r2Deficit', 0.02, 'chi2', 1, ...
                   'residAutocorr', 0, 'piNs', 500/w, 'piErrNs', 1, 'seconds', 0, ...
                   'shapeWarn', false);
    end
end

function n = report(name, ok)
    if ok; s = 'ok'; else; s = 'FAIL'; end
    fprintf('  %-72s %s\n', name, s);
    n = double(~ok);
end
