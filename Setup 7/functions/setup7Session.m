function [jsonLocation, h, awg, P] = setup7Session(resultsFolder)
%SETUP7SESSION Start a Setup 7 session: main, calibration paths, Helmholtz (state and coil
% calibration from expProperties.json), AWG settings.
%
%   [jsonLocation, h, awg, P] = setup7Session;                % newest calibration folder
%   [jsonLocation, h, awg, P] = setup7Session('2026-08-14');  % pin an older one
%
% Re-run it whenever you re-run main: main rebuilds the device objects, which drops the
% AWG settings below. setup7Paths('list') shows the available calibration folders.
%
% awg.keepOn = true leaves the MW outputs on at disconnect, so an experiment that re-runs
% without reloading the AWG still emits MW. A stopgap: set it false once every class
% re-enables its outputs (workChannels) on every run.

    if nargin < 1; resultsFolder = ''; end

    main;

    if isempty(resultsFolder); P = setup7Paths; else; P = setup7Paths(resultsFolder); end
    jsonLocation = P.json;

    h = Helmholtz;
    h.load(P.json);     % "Helmholtz" block of expProperties.json: currents + coil calibration
    h.control('auto');

    awg = getObjByName('TektronixAWG');
    if isempty(awg)
        fgs = FrequencyGenerator.getFG();
        for iFg = 1:numel(fgs)
            if contains(fgs{iFg}.name, 'TektronixAWG'); awg = fgs{iFg}; break; end
        end
    end
    if isempty(awg)
        error('setup7Session: no TektronixAWG found after main.');
    end

    awg.useVisadev           = true;    % ~2 ms queries instead of ~1.2 s
    awg.commDelayTime        = 0;       % no 5 ms busy-wait after every SCPI write
    awg.skipRedundantDeletes = true;
    awg.armWithoutChecking   = false;   % the check costs ~2 ms with visadev
    awg.keepOn               = true;
    awg.resetCommandCount;

    fprintf('AWG session: visadev %d, commDelayTime %g s, skipRedundantDeletes %d, keepOn %d\n', ...
            awg.useVisadev, awg.commDelayTime, awg.skipRedundantDeletes, awg.keepOn);
end
