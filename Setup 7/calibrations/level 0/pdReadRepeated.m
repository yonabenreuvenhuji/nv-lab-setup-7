function out = pdReadRepeated(spcm, integrationTime, nRepeats, channelConfig)
%PDREADREPEATED Read the photodiode N times and report the level with an honest error.
%
%   out = pdReadRepeated(spcm, integrationTime, nRepeats)
%   out = pdReadRepeated(spcm, integrationTime, nRepeats, channelConfig)
%
% Shared by the Level-0 calibrations. Returns, per digitizer input:
%   mean        mean over all repeats
%   sterrWithin the error the digitizer itself reports, averaged in quadrature.
%               This is photon/electronic noise within one integration.
%   sterrAcross the scatter OF THE REPEAT MEANS. This includes drift.
%   sterr       the larger of the two -- see below
%   drifted     true when sterrAcross is more than 2x sterrWithin, AND a within
%               estimate actually exists. On Setup 7 it does not: readFromTime
%               computes its error over the samples in one read, prepareReadByTime
%               asks the digitizer for a single segment, so ste() of one number is
%               exactly 0 and every read reported "drifting". A comparison against
%               zero is not a test.
%
% WHY TWO ERRORS. Averaging N reads and quoting the within-read error assumes the
% level is constant, which for a laser and a photodiode it is not. If the level
% drifts during the measurement, the within-read error stays small and lies to
% you while the scatter across repeats grows. Reporting both, and taking the
% larger, means a drifting measurement declares itself instead of looking precise.
% Nothing in the existing calibrations records any uncertainty at all, so this is
% the piece everything else here is built on.
%
% channelConfig follows spcm.currentChannelGUI:
%   1 fluorescence (input_1)   2 laser monitor (input_2)   3 both
%   4 differential             5 referenced               6 laser block
% Default 3, which gives fluorescence and laser monitor from the same read, so
% the two are measured simultaneously rather than sequentially.

    if ~exist('channelConfig', 'var') || isempty(channelConfig); channelConfig = 3; end
    if ~exist('nRepeats', 'var') || isempty(nRepeats); nRepeats = 10; end

    savedConfig = spcm.currentChannelGUI;
    spcm.currentChannelGUI = channelConfig;

    try
        spcm.prepareReadByTime(integrationTime);
        % Take value and error from the SAME read. An earlier version called
        % readFromTime twice for k = 1 -- once for the value, once for the error --
        % which wasted a read and paired a value with a different read's error.
        [v1, e1] = spcm.readFromTime();
        nCh = numel(v1);

        vals  = nan(nRepeats, nCh);
        errs  = nan(nRepeats, nCh);
        vals(1, :) = v1(:)';
        errs(1, :) = e1(:)';

        for k = 2:nRepeats
            [v, e] = spcm.readFromTime();
            vals(k, :) = v(:)';
            errs(k, :) = e(:)';
        end
        spcm.clearTimeRead();
    catch err
        try; spcm.clearTimeRead(); catch; end
        spcm.currentChannelGUI = savedConfig;
        rethrow(err);
    end

    spcm.currentChannelGUI = savedConfig;

    out.raw          = vals;
    out.nRepeats     = nRepeats;
    out.mean         = mean(vals, 1, 'omitnan');
    % Within-read error, combined over repeats: errors add in quadrature and the
    % mean of N of them is sqrt(sum(e^2))/N.
    out.sterrWithin  = sqrt(sum(errs.^2, 1, 'omitnan')) ./ sum(~isnan(errs), 1);
    if nRepeats > 1
        out.sterrAcross = std(vals, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(vals), 1));
    else
        out.sterrAcross = nan(1, nCh);
    end
    out.sterr    = max(out.sterrWithin, out.sterrAcross);
    out.hasWithin = out.sterrWithin > 0;
    % Relative scatter, which is what a human should judge drift by. 12 sigma on a
    % 50 uV error bar is 0.6 mV; whether that matters depends on the size of the
    % signal, and sigma alone hides that.
    out.relScatter = out.sterrAcross ./ max(abs(out.mean), eps);
    % DRIFT is judged on relative scatter, NOT against the within-read error.
    % Two attempts to use the within-read error both misfired: it comes out at or
    % near zero on this setup (readFromTime averages a read down to one number, so
    % its ste has nothing to work with), and "bigger than twice zero" is true for
    % any real measurement. So every single line printed "drifting", including runs
    % that were stable to 4 microvolts. A threshold on the actual relative scatter
    % says something; a comparison against zero never did.
    out.drifted  = out.relScatter > 0.01;      % 1% of the level
end
