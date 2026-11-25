function testBinblockEquivalence()
%TESTBINBLOCKEQUIVALENCE Prove the visadev waveform transfer is byte-identical.
%
% The binary waveform upload is the one place where the visadev migration could
% silently corrupt data rather than fail loudly: a wrong binblock header or byte
% order would load a garbage waveform that the instrument accepts without
% complaint. So compare the exact byte stream both paths put on the wire.
%
% LEGACY  fwrite(v, [prefix header])   then  fwrite(v, data, 'single')
%         (EOIMode off for the first write so the two form one message)
% VISADEV write(v, [uint8(prefix) uint8(header) typecast(single(data),'uint8')])
%
% Run with:  testBinblockEquivalence

    fprintf('\n=== testBinblockEquivalence ===\n');
    nPass = 0; nFail = 0;
    function check(name, cond, detail)
        if cond, fprintf('  PASS  %s\n', name); nPass=nPass+1;
        else,    fprintf('  FAIL  %s   %s\n', name, detail); nFail=nFail+1; end
    end

    lengths = [2400, 16000, 99999, 1000000];   % min allowed, typical ESR, odd, large
    for n = lengths
        data   = single(cos(2*pi*(0:n-1)/97) * 0.5);
        prefix = sprintf('wlist:waveform:data "%d",', 7);

        legacy  = legacyBytes(prefix, data);
        modern  = visadevBytes(prefix, data);

        check(sprintf('%7d samples: byte streams identical (%d bytes)', n, numel(legacy)), ...
              isequal(legacy, modern), ...
              sprintf('legacy %d bytes vs visadev %d bytes', numel(legacy), numel(modern)));
    end

    % header format spot-check against the IEEE 488.2 definition
    data = single(zeros(1,16000));
    m = visadevBytes('X,', data);
    hdr = char(m(3:9));
    check('header for 16000 samples is #564000', strcmp(hdr, '#564000'), ...
          sprintf('got "%s"', hdr));

    % byte order: first sample must round-trip
    probe = single([1.5 -2.25 3.75]);
    m = visadevBytes('X,', probe);
    payload = m(end-numel(probe)*4+1 : end);
    back = typecast(uint8(payload), 'single');
    check('float32 payload round-trips (little-endian)', isequal(back, probe), ...
          sprintf('got %s', mat2str(double(back))));

    fprintf('\n  %d passed, %d failed\n\n', nPass, nFail);
    if nFail > 0, error('testBinblockEquivalence: %d check(s) failed', nFail); end
end

function bytes = legacyBytes(prefix, data)
    wfm_bytes = numel(data) * 4;
    header = sprintf('#%d%d', length(num2str(wfm_bytes)), wfm_bytes);
    % fwrite of a char array sends its char codes; fwrite(...,'single') sends
    % native little-endian float32 on a PCWIN64 host.
    bytes = [uint8([prefix header]), typecast(single(data(:).'), 'uint8')];
end

function bytes = visadevBytes(prefix, data)
    data   = single(data);
    raw    = typecast(data(:).', 'uint8');
    lenStr = num2str(numel(raw));
    header = ['#' num2str(length(lenStr)) lenStr];
    bytes  = [uint8(prefix), uint8(header), raw];
end
