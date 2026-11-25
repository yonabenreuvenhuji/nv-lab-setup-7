function testWalkCompareMetrics()
%TESTWALKCOMPAREMETRICS Offline check that the agreement test behaves.
%
% The point of this test is the case that burned us: two CORRECT but NOISY
% measurements of the same spectrum must PASS even though their raw correlation is
% well below 0.9. A flat 0.9 threshold failed walk mode repeatedly while its graph
% visibly matched the reference.

    fprintf('\n=== testWalkCompareMetrics ===\n');
    nPass = 0; nFail = 0;
    function check(name, cond, detail)
        if cond, fprintf('  PASS  %s\n', name); nPass=nPass+1;
        else,    fprintf('  FAIL  %s   %s\n', name, detail); nFail=nFail+1; end
    end

    randn('seed', 3);
    nF = 24;
    k = (1:nF).';
    truth = 1 - 0.022*exp(-((k-12)/2.2).^2);      % 2.2% dip, a few points wide

    fprintf('\n  %-34s %-8s %-8s %-8s %-9s %s\n', ...
            'scenario', 'r', 'rMax', 'rRatio', 'sigmaD', 'verdict');

    function m = scenario(name, noise, transform)
        e1 = noise*ones(nF,1); e2 = noise*ones(nF,1);
        ref  = truth + noise*randn(nF,1);
        walk = truth + noise*randn(nF,1);
        walk = transform(walk);
        m = walkCompareMetrics(ref, e1, walk, e2);
        fprintf('  %-34s %-+8.3f %-8.3f %-8.2f %-9.2f %s\n', name, ...
                m.correlation, m.rMax, m.rRatio, m.sigmaDiff, ...
                BooleanHelper_ifTrueElse(m.pass, 'PASS', 'FAIL'));
    end

    % 1. clean, low noise: should pass with high r
    m = scenario('correct, low noise (repeats~100)', 0.0008, @(x) x);
    check('1 correct + low noise passes', m.pass, 'rejected a correct pair');

    % 2. THE CASE THAT MATTERED: correct but noisy. r is low, must still pass.
    m = scenario('correct, high noise (repeats~20)', 0.0035, @(x) x);
    check(sprintf('2 correct + high noise passes despite r = %.2f', m.correlation), ...
          m.pass, 'this is the false failure that stalled walk mode');
    check('2b rRatio recognises the noise ceiling', m.rRatio > 0.7, ...
          sprintf('rRatio %.2f', m.rRatio));

    % 3. shifted by one point: must fail, and must SAY the shift is 1
    m = scenario('shifted by +1 point', 0.0008, @(x) circshift(x, 1));
    check('3 a one-point shift fails', ~m.pass, 'accepted a shifted spectrum');
    % Sign convention: bestShift is the shift APPLIED TO walk that realigns it, so
    % a walk lagging ref by +1 point is corrected by -1 and reports -1.
    check(sprintf('3b reports a one-point offset (bestShift = %+d)', m.bestShift), ...
          abs(m.bestShift) == 1, sprintf('got %+d', m.bestShift));

    % 3c THE OTHER FALSE FAILURE: a noisy but correct pair must NOT be labelled
    % shifted just because one of the 11 candidate shifts beats zero by 0.0003.
    % That is what turned the 96-point walk run into a FAIL.
    nFsave = nF; nF = 96; k = (1:nF).';
    truth = 1 - 0.019*exp(-((k-48)/8).^2);
    m = scenario('correct, noisy, 96 points', 0.0022, @(x) x);
    check(sprintf('3c no spurious shift on a noisy correct pair (raw argmax %+d, reported %+d)', ...
                  m.rawBestShift, m.bestShift), ...
          m.bestShift == 0, sprintf('reported a shift of %+d on correct data', m.bestShift));
    check('3d and it passes', m.pass, sprintf('rRatio %.2f, sigmaDiff %.2f', m.rRatio, m.sigmaDiff));

    % 3e a REAL one-point shift at 96 points must still be caught -- at the SNR
    % the real sweep actually has. The measured run was contrast 1.84% with a
    % median error near 4e-4, i.e. a ~48 sigma dip, so that is what is used here.
    %
    % The limitation this makes explicit: whether a one-point shift is detectable
    % by correlation depends on SNR and on how many points span the dip. At
    % 48 sigma it moves r far more than the noise on r and is caught. At the much
    % noisier 3c settings, or against a dip spanning many points, a single-point
    % shift moves r by less than the noise on r -- no correlation test can
    % separate it from agreement there, and it also barely distorts the spectrum,
    % which is why that blind spot is acceptable and why firstPointSigma is
    % reported as its own number rather than folded into r.
    truth = 1 - 0.019*exp(-((k-48)/2.2).^2);
    m = scenario('shifted by +1, 96 points, real SNR', 0.0004, @(x) circshift(x, 1));
    check(sprintf('3e a real shift still detected at 96 points (%+d)', m.bestShift), ...
          abs(m.bestShift) == 1, sprintf('got %+d, margin %.3f', m.bestShift, m.shiftMargin));

    nF = nFsave; k = (1:nF).';
    truth = 1 - 0.022*exp(-((k-12)/2.2).^2);

    % 4. scrambled: must fail (this is what jump mode produced)
    m = scenario('scrambled frequencies', 0.0008, @(x) x(randperm(nF)));
    check('4 scrambled fails', ~m.pass, 'accepted scrambled data');

    % 5. stuck on one frequency: flat, must fail
    m = scenario('stuck on one frequency', 0.0008, @(x) mean(x)*ones(nF,1));
    check('5 flat/stuck fails', ~m.pass, 'accepted a flat spectrum');

    fprintf('\n  %d passed, %d failed\n\n', nPass, nFail);
    if nFail > 0, error('testWalkCompareMetrics: %d check(s) failed', nFail); end
end

function s = BooleanHelper_ifTrueElse(c, a, b)
    if c, s = a; else, s = b; end
end
