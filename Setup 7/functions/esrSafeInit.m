function initUs = esrSafeInit(detectionUs, gapUs, jsonInitUs)
%ESRSAFEINIT The largest laserInitializationDuration a CW-ESR window can survive.
%
%   a.laserInitializationDuration = esrSafeInit(a.detectionDuration, ...
%                                              a.delayBetweenDetections, ...
%                                              a.laserInitializationDuration);
%
% One line, called wherever a CW ESR window is set. See esrWindowPick for the full argument;
% the short version:
%
% In CW ESR the AWG loops a short microwave burst for  detectionDuration +
% laserInitializationDuration,  while one repeat of the sequence lasts  2 *
% (detectionDuration + delayBetweenDetections).  The AWG cannot accept the next trigger
% while it is still playing, so if the burst is longer than a repeat, repeats get NO
% MICROWAVE and the spectrum comes back as a comb -- alternate points at exactly ratio 1
% with the dip only on the others. Nothing errors and nothing warns.
%
% laserInitializationDuration is NOT an initialisation in CW ESR. That sequence has no init
% event; the green laser is on throughout. The field only sets how much extra microwave the
% AWG plays past the window, as padding. AWGExpESR's constructor uses 10 us for that job.
% expProperties.json holds a much larger value because the PULSED experiments need a real
% initialisation, and loadJsonParams copies it in regardless -- so a level-2 edit to the json
% silently breaks every CW ESR that was working. Which is what happened: a 10 us window was
% fine until laserInitializationDuration went from 26.5 to 70.6.
%
% This caps the padding at 75% of one repeat, leaving at least 1.33x headroom at any window,
% and never returns more than the value it was given -- so on a healthy json it is a no-op.

    if ~exist('gapUs', 'var') || isempty(gapUs); gapUs = 10; end
    if ~exist('jsonInitUs', 'var') || isempty(jsonInitUs); jsonInitUs = 10; end
    if ~isfinite(detectionUs) || detectionUs <= 0
        error('esrSafeInit: detectionDuration must be a positive number of microseconds.');
    end

    repeatUs = 2 * (detectionUs + gapUs);
    initMax  = 0.75 * repeatUs - detectionUs;      % burst <= 75% of one repeat
    initUs   = max(1, min(jsonInitUs, initMax));

    if initUs < jsonInitUs - 1e-9
        warning('esrSafeInit:capped', ...
                ['laserInitializationDuration capped %g -> %.3f us for a %g us CW window ' ...
                 '(one repeat is %.1f us). The json value would have made the AWG burst ' ...
                 'longer than a repeat, dropping triggers and producing a comb spectrum.'], ...
                jsonInitUs, initUs, detectionUs, repeatUs);
    end
end
