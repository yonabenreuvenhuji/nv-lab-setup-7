function testWalkReduce()
%TESTWALKREDUCE Offline check that walkReduce matches Experiment.processData.
%
% No hardware. Run with:  testWalkReduce
%
% The walk acquisition mode reduces its data with walkReduce instead of
% processData, because the frequency axis is interleaved into the stream. If the
% two ever disagree, the walk path would produce a subtly different spectrum from
% the classic path and the difference would be very hard to spot by eye. So this
% builds the same measurement both ways and asserts they agree to machine
% precision.
%
% Covers: the photodiode branch (reference sterr identically zero), the plain
% counting branch (kcps normalisation), the gated-integrator branch, and the
% discard-first-repeat threshold at exactly 10 and 11 passes.

    fprintf('\n=== testWalkReduce ===\n');
    nPass = 0; nFail = 0;
    function check(name, cond, detail)
        if cond, fprintf('  PASS  %s\n', name); nPass=nPass+1;
        else,    fprintf('  FAIL  %s   %s\n', name, detail); nFail=nFail+1; end
    end

    rand('seed', 11); randn('seed', 11);

    configs = { struct('hasPhotodiode',true,  'detectionWithGI',false, 'GIgain',1, ...
                       'detectionDuration',100, 'referenceDetectionDuration',100), ...
                struct('hasPhotodiode',false, 'detectionWithGI',false, 'GIgain',1, ...
                       'detectionDuration',100, 'referenceDetectionDuration',100), ...
                struct('hasPhotodiode',false, 'detectionWithGI',false, 'GIgain',1, ...
                       'detectionDuration',300, 'referenceDetectionDuration',500), ...
                struct('hasPhotodiode',true,  'detectionWithGI',true,  'GIgain',7.5, ...
                       'detectionDuration',100, 'referenceDetectionDuration',100) };
    names = {'photodiode', 'counting', 'counting, unequal windows', 'photodiode + GI'};

    for c = 1:length(configs)
        opts = configs{c};
        for nPasses = [5, 10, 11, 60]
            nFreq = 7;
            raw = zeros(2, nFreq, nPasses);
            for k = 1:nFreq
                level = 100 + 3*k;
                dip   = 1 - 0.08*exp(-((k-4)/1.3)^2);
                raw(2,k,:) = level  * (1 + 0.02*randn(1, nPasses));   % reference
                raw(1,k,:) = level*dip*(1 + 0.02*randn(1, nPasses));  % signal
            end

            [sigWalk, errWalk] = walkReduce(raw, 2, nFreq, opts);

            % Same data, reduced one frequency at a time through a replica of
            % Experiment.processData -- exactly how the classic path would see it.
            sigRef = zeros(2, nFreq); errRef = zeros(2, nFreq);
            for k = 1:nFreq
                stream = reshape(raw(:, k, :), 2*nPasses, 1);   % [s r s r ...]
                [sg, er] = processDataReplica(stream, nPasses, opts);
                sigRef(:, k) = sg(:); errRef(:, k) = er(:);
            end

            dSig = max(abs(sigWalk(:) - sigRef(:)));
            dErr = max(abs(errWalk(:) - errRef(:)));
            check(sprintf('%-26s nPasses=%2d : signal matches (max dev %.2e)', ...
                          names{c}, nPasses, dSig), dSig < 1e-9, 'signal mismatch');
            check(sprintf('%-26s nPasses=%2d : sterr  matches (max dev %.2e)', ...
                          names{c}, nPasses, dErr), dErr < 1e-9, 'sterr mismatch');
        end
    end

    % The photodiode branch replaces the reference column with a constant, so its
    % error must be negligible -- that is what makes the downstream ratio error
    % exact. Not tested against exactly zero: averaging N identical values is not
    % bit-exact for N that is not a power of two, so both this and the real
    % processData leave residue around 1e-14. Assert it is negligible relative to
    % the level instead.
    opts = configs{1};
    raw  = 100 + randn(2, 5, 40);
    [sig, err] = walkReduce(raw, 2, 5, opts);
    relRefErr = max(err(2,:) ./ abs(sig(2,:)));
    check(sprintf('photodiode reference sterr is negligible (%.1e relative)', relRefErr), ...
          relRefErr < 1e-12, sprintf('relative %g', relRefErr));

    % REGRESSION: a single pass arrives as [2 x nFreq], not [2 x nFreq x 1],
    % because MATLAB drops trailing singleton dimensions. This is what broke the
    % first version of walk mode on its very first pass.
    opts = configs{1};
    onePass = 100 + randn(2, 9);            % note: 2-D, no third dimension
    okOnePass = true;
    try
        [s1, e1] = walkReduce(onePass, 2, 9, opts);
        okOnePass = isequal(size(s1), [2 9]) && isequal(size(e1), [2 9]);
    catch
        okOnePass = false;
    end
    check('single pass arriving as 2-D [2 x nFreq] is handled', okOnePass, ...
          'walkReduce rejected a dropped trailing singleton');

    fprintf('\n  %d passed, %d failed\n\n', nPass, nFail);
    if nFail > 0, error('testWalkReduce: %d check(s) failed', nFail); end
end


function [signal, sterrv] = processDataReplica(rawData, repeats, opts)
% Faithful replica of the parts of Experiment.processData that matter here.
    n = repeats;
    m = length(rawData)/n;
    if mod(10, m)
        m = 2;                          % detectionPeriodsPerRepeat
        n = floor(length(rawData)/m);
        rawData = rawData(1:m*n);
    end
    s = (reshape(rawData, m, n))';

    timeNormalization = [opts.detectionDuration opts.referenceDetectionDuration]*1e-6;

    start = 1; if size(s,1) > 10, start = 2; end

    if opts.hasPhotodiode
        refSignal = s(:, 1:2:end) ./ s(:, 2:2:end);
        ref = mean(s(start:end, 2:2:end));
        s(:, 1:2:end) = ref .* refSignal;
        s(:, 2:2:end) = ref +  refSignal * 0;
    end

    signal = mean(s(start:end, :));
    sterrv = ste(s(start:end, :));

    if ~opts.hasPhotodiode
        signal = signal./timeNormalization/1e3;
        sterrv = sterrv./timeNormalization/1e3;
    else
        if opts.detectionWithGI
            tn = timeNormalization / 1e-6;
            signal = signal./tn./opts.GIgain;
            sterrv = sterrv./tn./opts.GIgain;
        end
    end
end
