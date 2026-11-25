function results = calibLaserSaturation(aomValues, integrationTime, nRepeats, dark, referenceAom)
%CALIBLASERSATURATION NV fluorescence versus green power, and the saturation point.
%
%   results = calibLaserSaturation()
%   results = calibLaserSaturation(-0.6:0.04:0, 0.1, 8)
%   results = calibLaserSaturation([], [], [], darkResults)   % reuse calibDarkLevel
%   results = calibLaserSaturation([], [], [], d, -0.2)       % pick the reference point
%
% LEVEL 0 CALIBRATION. Depends only on calibDarkLevel (for the zero) and on the
% spot already being on an NV. Nothing about microwave, timing or field enters.
%
% WHY THE X AXIS IS NOT THE AOM VOLTAGE. The obvious sweep is over the AOM control
% voltage, but an AOM's transmission is not linear in its drive, so fitting a
% saturation curve against control volts fits the AOM's response and the NV's
% together and the result is not P_sat in any meaningful sense.
%
% This setup does not need that compromise: digitizer input_2 is the laser
% sampling channel, a direct monitor of the light actually delivered. So the AOM
% voltage is only the knob that is turned, and the fit is fluorescence against
% MEASURED power. Both come from the same read (currentChannelGUI = 3), so they
% are simultaneous rather than sequential and drift affects them together.
%
% THE MODEL
%
%   F(P) = Fsat * P / (P + Psat)  +  c * P
%
% The first term is the two-level saturation of the NV. The second is everything
% that scales linearly with power and does not saturate: scattered green reaching
% the detector, substrate fluorescence, diamond Raman. Leaving it out is the
% classic way to get a wrong Psat -- an unsaturating linear background masquerades
% as a larger Fsat and pushes Psat high. It is fitted here, and reported, so you
% can see how much of your signal is not NV.
%
% ORDER AND DRIFT. The power points are visited in RANDOM order, and the first
% point is re-measured at the end. Random order stops a monotonic drift (heating,
% laser warm-up, focus creep) from being absorbed into the fitted curve as
% curvature; the repeat gives a direct measure of how much the setup moved during
% the scan. A saturation curve taken in ascending order on a drifting setup looks
% perfectly convincing and is wrong.
%
% AND AN INTERLEAVED REFERENCE, which is the part that actually helps. A first
% version relied on averaging, and on this setup that was the wrong tool: the read
% noise is ~50 uV on a 20 mV signal while the SAME point moved 5% across one scan.
% Drift was 20x the noise, so more repeats bought nothing and cost time, which
% made the drift worse.
%
% So every power point is now bracketed by a measurement at one fixed reference
% AOM setting, and the ratio is what gets fitted. Slow changes in collection
% efficiency -- focus creep, position drift, laser warm-up -- are common to both
% and divide out. It is the same trick the ESR sequence already uses with its
% signal and reference windows, for the same reason. The reference series is also
% reported on its own, so you can see how much the setup moved and when.
%
% WHAT TO DO WITH Psat. It is a measurement, not a recommendation. Higher power
% gives more photons but also more heating, more ionisation to NV0, and generally
% worse ESR contrast. Labs commonly work at or below Psat. This function reports
% Psat, the AOM voltage that reaches it, and the voltages for a few fractions of
% it, and leaves the choice to you.
%
% It does NOT write to the json.

    spcm  = getObjByName(Spcm.NAME);
    laser = getObjByName(LaserGate.GREEN_LASER_NAME);
    if isempty(spcm);  throwBaseObjException(Spcm.NAME); end
    if isempty(laser); throwBaseObjException(LaserGate.GREEN_LASER_NAME); end
    if ~laser.isAomAvail || ~laser.aom.canSetValue
        error('This laser has no settable AOM, so power cannot be swept from software.');
    end

    aom = laser.aom;
    if ~exist('aomValues', 'var') || isempty(aomValues)
        % Full allowed span. NOTE the AOM value is a VOLTAGE here, not a percent:
        % AomNiDaqControlled passes NiDaq.UNITS, and on Setup 7 the json gives
        % minVal -0.6, maxVal 0. Which end is full power is a property of the
        % hardware, not of the numbers, so the sweep covers the whole range and
        % the laser monitor tells us which way round it is.
        aomValues = linspace(aom.minValue, aom.maxValue, 13);
    end
    if ~exist('integrationTime', 'var') || isempty(integrationTime); integrationTime = 0.1; end
    if ~exist('nRepeats', 'var') || isempty(nRepeats); nRepeats = 4; end
    if ~exist('dark', 'var'); dark = []; end
    if ~exist('referenceAom', 'var') || isempty(referenceAom)
        % Mid-range by default. It only has to be a setting that gives a decent
        % signal and stays fixed; its absolute value never enters the result.
        referenceAom = mean([aom.minValue, aom.maxValue]);
    end

    aomValues = aomValues(:).';
    nPts = numel(aomValues);

    fprintf('\n=== calibLaserSaturation: %d power points, %g s x %d per point ===\n', ...
            nPts, integrationTime, nRepeats);
    % NOTE the units string says '%' but the values are VOLTS. LaserPartAbstract's
    % constructor is (name, minVal, maxVal, func, units) and AomNiDaqControlled
    % calls it with 4 arguments, so NiDaq.UNITS lands in the 'func' slot and 'units'
    % falls back to its default of '%'. The numbers are niDaq volts.
    fprintf('AOM range [%g, %g] V  (the class labels this "%s" -- see the code comment)\n', ...
            aom.minValue, aom.maxValue, aom.units);

    % ---- zero point ------------------------------------------------------
    if isempty(dark)
        fprintf('No dark measurement passed in -- taking one now (laser gate closed).\n');
        savedOn = laser.isOn;
        laser.isOn = false;
        pause(0.5);
        darkRead = pdReadRepeated(spcm, integrationTime, nRepeats, 3);
        laser.isOn = savedOn;
    elseif isstruct(dark) && isfield(dark, 'darkWithDetector')
        darkRead = dark.darkWithDetector;        % straight from calibDarkLevel
        fprintf('Using the dark level from calibDarkLevel.\n');
    else
        darkRead = dark;
    end
    fprintf('dark: fluorescence %+.6f V, monitor %+.6f V\n', darkRead.mean(1), darkRead.mean(2));

    results = struct('aomValues', aomValues, 'integrationTime', integrationTime, ...
                     'nRepeats', nRepeats, 'dark', darkRead);

    savedAom = aom.value;
    savedOn  = laser.isOn;

    fl   = nan(1, nPts);  flE  = nan(1, nPts);
    pw   = nan(1, nPts);  pwE  = nan(1, nPts);
    ref  = nan(1, nPts);  refE = nan(1, nPts);      % interleaved reference

    try
        laser.isOn = true;

        % Random visit order -- see the header.
        order = randperm(nPts);
        fprintf('\nreference AOM setting: %.4f V (measured before every point)\n', referenceAom);
        fprintf('\n%-6s %-10s %-13s %-13s %-13s %s\n', ...
                'point', 'AOM (V)', 'monitor (V)', 'fl (V)', 'ref (V)', 'fl/ref');
        for k = 1:nPts
            i = order(k);

            % Reference first, then the point. Same read structure for both.
            aom.value = referenceAom;
            pause(0.3);
            rr = pdReadRepeated(spcm, integrationTime, nRepeats, 3);
            ref(i)  = rr.mean(1) - darkRead.mean(1);
            refE(i) = sqrt(rr.sterr(1)^2 + darkRead.sterr(1)^2);

            aom.value = aomValues(i);
            pause(0.3);                          % analog settle
            r = pdReadRepeated(spcm, integrationTime, nRepeats, 3);
            fl(i)  = r.mean(1) - darkRead.mean(1);
            flE(i) = sqrt(r.sterr(1)^2 + darkRead.sterr(1)^2);
            pw(i)  = r.mean(2) - darkRead.mean(2);
            pwE(i) = sqrt(r.sterr(2)^2 + darkRead.sterr(2)^2);

            fprintf('%-6d %-10.4f %+11.6f %+11.6f %+11.6f  %.4f\n', ...
                    i, aomValues(i), pw(i), fl(i), ref(i), fl(i)/max(ref(i), eps));
        end

        % Normalised signal and its error. Ratio errors add in quadrature
        % relatively, which is the honest propagation for fl/ref.
        results.reference      = ref;
        results.referenceSterr = refE;
        norm  = fl ./ ref;
        normE = abs(norm) .* sqrt((flE./max(abs(fl),eps)).^2 + (refE./max(abs(ref),eps)).^2);
        results.normalised      = norm;
        results.normalisedSterr = normE;

        % How much did the reference itself move? This is the drift, measured
        % directly rather than inferred.
        refSpread = (max(ref) - min(ref)) / max(abs(mean(ref)), eps);
        results.referenceSpread = refSpread;
        fprintf('\nreference moved %.1f%% peak-to-peak across the scan\n', 100*refSpread);
        if refSpread > 0.03
            fprintf(['   That drift is why averaging harder does not help: it is a real\n' ...
                     '   change in the setup, not read noise. Dividing by the reference is\n' ...
                     '   what removes it. If it is much larger than a few percent, consider\n' ...
                     '   turning tracking on or letting the laser settle first.\n']);
        end

        % Drift check: re-measure the first point visited.
        iFirst = order(1);
        aom.value = aomValues(iFirst);
        pause(0.3);
        rEnd = pdReadRepeated(spcm, integrationTime, nRepeats, 3);
        flEnd = rEnd.mean(1) - darkRead.mean(1);
        dSig  = (flEnd - fl(iFirst)) / max(sqrt(rEnd.sterr(1)^2 + flE(iFirst)^2), eps);
        % As a FRACTION OF THE SCAN'S FULL SWING, not of the point itself.
        %
        % Dividing by the point is only meaningful when the point carries signal, and with a
        % random visit order it often does not. A run whose first visited point was AOM
        % -0.5 V -- where the fluorescence is -71 uV, i.e. zero -- reported "+35.3% drift"
        % from a change of 25 uV, and level 0 then repeated that 35.3% as though the setup
        % had fallen apart. It had not: the same run was +1.1 sigma, which is no drift.
        % The swing is a denominator that is always the same size as the measurement.
        swing = max(fl) - min(fl);
        dPct  = 100 * (flEnd - fl(iFirst)) / max(swing, eps);
        results.driftSigma       = dSig;
        results.driftPercent     = dPct;
        results.driftAbs         = flEnd - fl(iFirst);
        results.driftPointAom    = aomValues(iFirst);
        % Could the re-measured point have shown drift at all?
        results.driftPointUsable = abs(fl(iFirst)) > 0.05 * max(abs(fl));

    catch err
        aom.value = savedAom; laser.isOn = savedOn;
        fprintf('\n  calibLaserSaturation aborted: %s\n', err.message);
        rethrow(err);
    end

    aom.value = savedAom;
    laser.isOn = savedOn;

    results.fluorescence = fl;  results.fluorescenceSterr = flE;
    results.power        = pw;  results.powerSterr        = pwE;

    % Report BOTH. Sigma alone is misleading when the error bars are tens of
    % microvolts: 12 sigma sounded catastrophic on a run where the actual change
    % was 5%, which is a real effect but not a ruined measurement.
    fprintf('\ndrift check: AOM %.3f V re-measured at the end, changed by %+.0f uV', ...
            results.driftPointAom, 1e6*results.driftAbs);
    fprintf('  = %+.1f%% of the scan''s full swing (%+.1f sigma)\n', dPct, dSig);
    if ~results.driftPointUsable
        fprintf(['   The random visit order put a near-dark point first, so this particular\n' ...
                 '   check is weak -- there is little signal there to drift. The interleaved\n' ...
                 '   reference above is the stronger drift measure and covers the whole scan.\n']);
    end
    % Needs BOTH a real size and real significance. A run reported "-2.0% (-0.1
    % sigma)" and still printed the warning, because the percentage test fired
    % alone -- but at 0.1 sigma there is no drift to warn about.
    if abs(dPct) > 2 && abs(dSig) > 3
        fprintf(['   A few percent of drift over the scan. Random visit order spreads it\n' ...
                 '   across the curve instead of letting it masquerade as curvature, so the\n' ...
                 '   shape is still meaningful -- but do not quote Psat to better than this.\n']);
    end

    % ---- is the AOM monotone? -------------------------------------------
    [~, srt] = sort(aomValues);
    pwSorted = pw(srt);
    results.monotone = all(diff(pwSorted) > -3*max(pwE));
    if ~results.monotone
        fprintf(['\n   NOTE the laser monitor is not monotone in the AOM voltage. That is a\n' ...
                 '   property of the AOM, not an error, and it is exactly why the fit uses\n' ...
                 '   measured power rather than AOM voltage as its x axis.\n']);
    end
    if pwSorted(end) < pwSorted(1)
        fprintf('   Full power is at the LOW end of the AOM range (%g V).\n', aom.minValue);
    else
        fprintf('   Full power is at the HIGH end of the AOM range (%g V).\n', aom.maxValue);
    end

    % ---- fit -------------------------------------------------------------
    % ---- re-zero the monitor IN SITU -------------------------------------
    % The laser-gate-closed dark is the wrong zero for the monitor channel, and the
    % data says so directly: at the lowest AOM setting the fluorescence is zero --
    % no light reaches the NV -- yet the dark-subtracted monitor reads about
    % -1.2 mV. A power axis cannot be negative. The detector baseline also wandered
    % from +4.8 mV to -0.3 mV between two sessions, which is larger than the
    % monitor's entire ~3.5 mV swing, so an absolute dark taken at a different
    % moment simply is not a usable zero for this channel.
    %
    % So the zero is taken from the sweep itself: the monitor reading at the
    % settings where the NV sees no light. Same channel, same session, same
    % baseline, and it makes the x axis physical.
    dark0 = norm < 0.02 * max(norm);        % points where the NV is dark
    if sum(dark0) >= 2
        monitorZero = mean(pw(dark0));
        fprintf('\n--- monitor zero ---\n');
        fprintf('   from the laser-off dark : 0 (by construction)\n');
        fprintf('   from %d zero-fluorescence points : %+.6f V\n', sum(dark0), monitorZero);
        fprintf('   difference %+.6f V, i.e. %.0f%% of the monitor''s full swing\n', ...
                -monitorZero, 100*abs(monitorZero)/max(max(pw)-min(pw), eps));
        pw = pw - monitorZero;
        results.monitorZero = monitorZero;
        results.powerRezeroed = pw;
    else
        fprintf(['\n   Fewer than 2 points had zero fluorescence, so the monitor cannot be\n' ...
                 '   re-zeroed from the sweep. Extend the sweep to lower power.\n']);
        results.monitorZero = 0;
    end

    % ---- drop the AOM roll-over region -----------------------------------
    % Above its efficiency peak the AOM delivers less light, so fluorescence falls
    % while the monitor barely moves -- on this setup fluorescence drops 17% over a
    % span where the monitor changes 2%. Those two cannot both be right, so that
    % region says nothing about NV saturation and is excluded from the fit rather
    % than being allowed to bend it.
    [~, iPeak] = max(norm);
    [aomSortVal, aomSortIdx] = sort(aomValues);
    peakPos = find(aomSortIdx == iPeak);
    keep = false(1, nPts);
    keep(aomSortIdx(1:peakPos)) = true;      % up to and including the peak
    nDropped = nPts - sum(keep);
    if nDropped > 0
        fprintf('\n   Fluorescence peaks at AOM %.3f V. Dropping the %d point(s) above it:\n', ...
                aomValues(iPeak), nDropped);
        fprintf('   past the AOM''s efficiency peak, more drive means less light.\n');
    end
    results.keptForFit = keep;

    % ---- is the monitor fit to be an x axis at all? ------------------------
    % Check this BEFORE fitting, and refuse to fit if it fails. A fit against a
    % power axis that does not track power produces confident-looking coefficients
    % with meaningless values -- on this setup it returned Psat = 7.3 with an error
    % bar of 3e7, next to an R^2 of 0.82. Better to say the axis is unusable.
    %
    % The test: over the kept region, does the monitor RANK the points the same way
    % the fluorescence does? Fluorescence is monotone in delivered power (below
    % saturation it rises, at saturation it flattens -- it never falls), so if the
    % monitor disagrees about the ORDER, it is not measuring power.
    kf = find(keep & isfinite(norm) & isfinite(pw));
    if numel(kf) >= 5
        rMon = rankOf(pw(kf));
        rFl  = rankOf(norm(kf));
        rankAgreement = corrCoef(rMon, rFl);
        monScatter = median(abs(diff(pw(kf)))) / max(max(pw(kf)) - min(pw(kf)), eps);
        results.monitorRankAgreement = rankAgreement;
        results.monitorRelScatter    = monScatter;
        fprintf('\n--- is the monitor usable as a power axis? ---\n');
        fprintf('   rank agreement with fluorescence : %+.3f  (need > 0.9)\n', rankAgreement);
        fprintf('   point-to-point scatter           : %.0f%% of its full swing\n', 100*monScatter);
        results.monitorUsable = rankAgreement > 0.9 && monScatter < 0.15;
    else
        results.monitorUsable = false;
        results.monitorRankAgreement = NaN;
        fprintf('\n   Too few kept points to judge the monitor channel.\n');
    end

    if ~results.monitorUsable
        fprintf(['\n   NOT FITTING. The monitor does not order the points the way the\n' ...
                 '   fluorescence does, so it is not tracking the light delivered to the\n' ...
                 '   diamond, and any Psat from it would be meaningless however good the\n' ...
                 '   R^2 looked. Fluorescence never DECREASES with more power, so a\n' ...
                 '   disagreement in ordering is the monitor''s fault, not the NV''s.\n' ...
                 '\n' ...
                 '   Two causes, opposite fixes. If the sampler photodiode is noise-limited,\n' ...
                 '   more light on it helps. If it is compressing -- likely, since it sees\n' ...
                 '   the excitation beam, far brighter than the fluorescence, yet reports a\n' ...
                 '   SMALLER voltage than the main diode -- then more light makes it worse\n' ...
                 '   and it needs attenuation. Test by halving the light into the sampler\n' ...
                 '   and seeing whether its reading halves.\n' ...
                 '\n' ...
                 '   Free first step either way: max_voltage is +/-2 V for a signal of a few\n' ...
                 '   mV, about 1%% of full scale. A smaller input range costs nothing and\n' ...
                 '   gains resolution everywhere.\n' ...
                 '\n' ...
                 '   The brightest-setting result below needs none of this and is reported\n' ...
                 '   regardless.\n']);
        results = reportBrightest(aomValues, norm, normE, refSpreadPct(results), results);
        results = finishPlots(results, aomValues, pw, pwE, norm, normE, aom);
        results = laserVerdict(results, aom);
        return
    end

    % Fit the REFERENCE-NORMALISED signal. Fsat then comes out in units of the
    % reference level rather than volts, which is fine: Psat is what we are after
    % and it lives on the x axis.
    ok = isfinite(norm) & isfinite(pw) & pw > 0 & keep;
    if sum(ok) < 4
        fprintf('\n   Only %d usable points; not fitting. Check the monitor channel.\n', sum(ok));
        return
    end
    P = pw(ok).';  F = norm(ok).';  W = 1 ./ max(normE(ok).', eps).^2;

    ft = fittype('Fsat*P/(P+Psat) + c*P', 'independent', 'P', 'coefficients', {'Fsat','Psat','c'});
    % fitoptions(ft) already returns a curvefit.nlsqoptions for a nonlinear
    % fittype, and its Method is read-only -- setting it errors. Do not add it back.
    opts = fitoptions(ft);
    opts.Weights    = W;
    % c >= 0: a non-saturating background cannot be negative. Letting it go
    % negative gave the fit another direction to wander in and was part of why
    % repeated runs landed on different answers.
    opts.Lower      = [0, 0, 0];
    opts.StartPoint = [max(F)*1.5, median(P), 0];
    opts.Display    = 'off';
    [fitres, gof] = fit(P, F, ft, opts);
    ci = confint(fitres);
    co = coeffvalues(fitres);

    results.fit      = fitres;
    results.gof      = gof;
    results.Fsat     = co(1);   results.FsatSterr = (ci(2,1)-ci(1,1))/2;
    results.Psat     = co(2);   results.PsatSterr = (ci(2,2)-ci(1,2))/2;
    results.linear   = co(3);   results.linearSterr = (ci(2,3)-ci(1,3))/2;

    fprintf('\n--- fit  F = Fsat*P/(P+Psat) + c*P ---\n');
    fprintf('   Fsat  %.6f  +/- %.6f  V   (saturated NV fluorescence)\n', results.Fsat, results.FsatSterr);
    fprintf('   Psat  %.6f  +/- %.6f  V   (monitor units, half-saturation power)\n', results.Psat, results.PsatSterr);
    fprintf('   c     %.6f  +/- %.6f      (non-saturating background per unit power)\n', results.linear, results.linearSterr);
    fprintf('   R^2   %.4f\n', gof.rsquare);

    % Is the fit even well posed? Psat is only constrained if the data straddles
    % it. All-high-power data makes the saturating term flat, Fsat and c then trade
    % off freely, and the fit lands somewhere different on every run -- which is
    % exactly the instability seen here, and no amount of averaging cures it.
    results.PsatInRange = results.Psat > min(P) && results.Psat < max(P);
    if ~results.PsatInRange
        if results.Psat >= max(P)
            fprintf(['\n   Psat = %.6f is AT OR ABOVE the highest power available (%.6f).\n' ...
                     '   That is a physical result, not a fitting failure: the NV never\n' ...
                     '   saturates anywhere in this AOM range, so you are operating on the\n' ...
                     '   near-linear part of the curve throughout. Psat cannot be measured\n' ...
                     '   without more power at the sample -- more points or more averaging\n' ...
                     '   will not produce it. What you CAN take from this is the setting that\n' ...
                     '   delivers the most light, which is the fluorescence peak above.\n'], ...
                    results.Psat, max(P));
        else
            fprintf(['\n   Psat = %.6f lies BELOW the lowest measured power (%.6f), so every\n' ...
                     '   point sits on the flat part of the curve and Fsat trades off freely\n' ...
                     '   against c. Extend the sweep to lower power until the curve bends.\n'], ...
                    results.Psat, min(P));
        end
    end
    if results.PsatSterr > 0.5*results.Psat
        fprintf(['   The uncertainty on Psat is more than half its value, so it is not a\n' ...
                 '   measurement yet whatever the R^2 says.\n']);
    end

    % ---- linear regime, and why c is NOT the background here ---------------
    % Below saturation the NV's OWN fluorescence is linear in power, approximately
    % Fsat*P/Psat. So when Psat is far above the accessible power the two terms of
    % the model are degenerate -- a straight line can be built entirely from
    % either -- and c cannot be read as "background".
    %
    % An earlier version did read it that way and told the user that 86% of their
    % signal was background, and to go and check their collection optics and
    % centring. That was wrong: the fit had assigned everything to c because the
    % DATA is linear, which says nothing at all about how much light is not NV.
    linearOnly = ~results.PsatInRange || results.PsatSterr > 0.5*results.Psat;
    results.linearRegime = linearOnly;
    if linearOnly
        slope = results.linear + results.Fsat/max(results.Psat, eps);
        % Do NOT claim which coefficient is well determined. Across runs the line
        % was absorbed sometimes by c (c = 398 +/- 96, Fsat unconstrained) and
        % sometimes by Fsat/Psat (Fsat = 10.8 +/- 270, c = 0.0005 +/- 5610). Both
        % describe the same straight line. The SLOPE is what is determined; the
        % split between the terms is not, and an earlier version asserted the
        % opposite on a run where c had a relative error of 10^7.
        fprintf(['\n   LINEAR REGIME. The individual coefficients carry errors comparable to\n' ...
                 '   or larger than themselves -- which of them absorbs the line varies from\n' ...
                 '   run to run -- but their combination is well determined: the data is a\n' ...
                 '   straight line, F = %.0f * P. The curve never bends over the power\n' ...
                 '   available, so Psat is not measurable here, and that is a statement\n' ...
                 '   about the power range rather than about the noise.\n'], slope);
        fprintf(['\n   Note c must NOT be read as a background fraction in this regime. Below\n' ...
                 '   saturation the NV''s own response is linear too, so the saturating and\n' ...
                 '   linear terms are degenerate and the split between them is arbitrary.\n' ...
                 '   To measure the real background, park on a spot with no NV at the same\n' ...
                 '   power and compare -- that separates them properly, and this scan\n' ...
                 '   cannot.\n']);
    else
        % Only meaningful once the curve genuinely bends.
        atPsat    = results.Fsat/2;
        bgAtPsat  = results.linear * results.Psat;
        if atPsat + bgAtPsat > 0
            results.backgroundFractionAtPsat = bgAtPsat / (atPsat + bgAtPsat);
            fprintf('   At Psat, %.0f%% of the signal is the non-saturating background.\n', ...
                    100*results.backgroundFractionAtPsat);
            if results.backgroundFractionAtPsat > 0.3
                fprintf(['   That is a lot, and here it IS interpretable because the curve\n' ...
                         '   bends. Worth checking scattered green and the centring, since\n' ...
                         '   this background dilutes every contrast measured later.\n']);
            end
        end
    end

    % ---- AOM voltage for a few working points ---------------------------
    % Interpolate the measured monitor-vs-AOM curve rather than assuming a form.
    if results.linearRegime
        fprintf(['\n   Skipping the "AOM voltage for a fraction of Psat" table: Psat is not\n' ...
                 '   measurable in the linear regime, so those settings would be arbitrary.\n']);
        results = reportBrightest(aomValues, norm, normE, refSpreadPct(results), results);
        results = finishPlots(results, aomValues, pw, pwE, norm, normE, aom);
        results = laserVerdict(results, aom);
        return
    end
    fprintf('\n--- AOM voltage for a given fraction of Psat ---\n');
    fractions = [0.25 0.5 1.0 2.0];
    results.workingPoints = struct('fraction', {}, 'power', {}, 'aomValue', {});
    aomSorted = aomValues(srt);
    for f = fractions
        target = f * results.Psat;
        if target < min(pwSorted) || target > max(pwSorted)
            fprintf('   %4.2f x Psat = %.6f V  -- outside the measured range\n', f, target);
            continue
        end
        % Monotone interpolation on whichever direction the monitor runs.
        [pu, iu] = unique(pwSorted);
        v = interp1(pu, aomSorted(iu), target, 'linear');
        fprintf('   %4.2f x Psat = %.6f V  ->  AOM %.4f V\n', f, target, v);
        results.workingPoints(end+1) = struct('fraction', f, 'power', target, 'aomValue', v); %#ok<AGROW>
    end
    results = reportBrightest(aomValues, norm, normE, refSpreadPct(results), results);

    fprintf(['\n   Psat is a measurement, not a recommendation. More power means more\n' ...
             '   photons but also heating, NV0 ionisation and usually worse ESR contrast.\n' ...
             '   Pick the working point deliberately; nothing here writes it anywhere.\n']);

    % ---- plot ------------------------------------------------------------
    figure;
    subplot(2,1,1); hold on
    errorbar(pw, norm, normE, '.', 'MarkerSize', 12);
    Pfine = linspace(0, max(pw)*1.05, 200);
    plot(Pfine, fitres(Pfine), '-');
    plot(Pfine, results.linear*Pfine, '--');
    xlabel('measured power (monitor V, dark-subtracted)');
    ylabel('fluorescence / reference');
    legend('data', 'fit', 'linear background', 'Location', 'southeast');
    title(sprintf('Psat = %.4f +/- %.4f monitor V,  R^2 = %.4f', ...
                  results.Psat, results.PsatSterr, gof.rsquare));

    subplot(2,1,2);
    errorbar(aomValues, pw, pwE, '.-');
    xlabel(sprintf('AOM value (%s)', aom.units));
    ylabel('measured power (monitor V)');
    title('AOM response -- the knob, not the physics');

    results = laserVerdict(results, aom);
    fprintf('\n');
end

function results = reportBrightest(aomValues, norm, normE, refSpread, results)
    % The one number that is trustworthy whatever the fit does: where the most light
    % actually arrives. Straight from the normalised data, no model, and independent of the
    % monitor channel entirely.
    %
    % BUT NOT THE ARGMAX. An earlier version returned max(norm) and its AOM value, and the
    % comment recorded that "across three runs this came back at -0.15 or -0.16 V every
    % time". A later run put it at -0.10, and the two are not in conflict -- they are two
    % samples from a FLAT TOP. On that run -0.15 and -0.10 differed by 0.7% while the
    % interleaved reference itself moved 1.6% across the scan, so the peak position was not
    % resolved and picking the argmax was picking noise.
    %
    % So report the plateau: every setting that is not distinguishable from the best, and
    % recommend its midpoint. The midpoint is also the robust choice operationally -- it
    % puts ordinary AOM drift along the flat top instead of down a shoulder.
    [peakVal, iPeak] = max(norm);
    fprintf('\n--- brightest setting (model-free) ---\n');
    fprintf('   AOM %.3f V gives the highest normalised fluorescence (%.4f)\n', ...
            aomValues(iPeak), peakVal);

    % Indistinguishable = within the combined error bars, OR within the amount the
    % reference itself wandered. The second usually dominates and is the honest limit:
    % differences smaller than the drift were not measured.
    tolStat  = 2*sqrt(normE.^2 + normE(iPeak)^2);
    tolDrift = (max(refSpread, 0)/100) * peakVal;
    near     = (peakVal - norm) <= max(tolStat, tolDrift);

    % Contiguous run of near-peak settings AROUND the peak, in AOM order, so a stray
    % near-peak point at the far end of the range cannot widen the plateau.
    [aSort, ord] = sort(aomValues);
    nearSort = near(ord);
    jPeak = find(ord == iPeak, 1);
    lo = jPeak; while lo > 1 && nearSort(lo-1); lo = lo - 1; end
    hi = jPeak; while hi < numel(nearSort) && nearSort(hi+1); hi = hi + 1; end
    plateau = aSort(lo:hi);
    recommended = mean([min(plateau), max(plateau)]);

    results.brightestAom  = aomValues(iPeak);
    results.brightestNorm = peakVal;
    results.plateauAom    = [min(plateau) max(plateau)];
    results.recommendedAom = recommended;

    if numel(plateau) > 1
        fprintf(['   The top is FLAT from %.3f to %.3f V -- those settings are not\n' ...
                 '   distinguishable from the best given %.1f%% reference movement.\n'], ...
                min(plateau), max(plateau), refSpread);
        fprintf('   USE AOM %.3f V, the middle of that plateau.\n', recommended);
        fprintf(['   Not the argmax: on a flat top the argmax moves run to run, and the\n' ...
                 '   midpoint keeps AOM drift on the plateau instead of down a shoulder.\n']);
    else
        fprintf('   USE AOM %.3f V. No neighbouring setting matched it, so the peak is resolved.\n', ...
                recommended);
    end
    fprintf('   Measured, not fitted, and independent of the monitor channel.\n');
end

function results = finishPlots(results, aomValues, pw, pwE, norm, normE, aom)
% The no-usable-Psat figure, reached from TWO different branches, which is why the titles
% used to be wrong. The old top title read "no fit: the monitor is not tracking delivered
% power" unconditionally -- so a run that reached here through the LINEAR-REGIME branch,
% with rank agreement 0.991 and a perfectly good monitor, was told its monitor was broken.
% Say which of the two actually happened.
    linear = isfield(results, 'linearRegime') && results.linearRegime;
    figure;
    subplot(2,1,1); hold on
    errorbar(pw, norm, normE, '.', 'MarkerSize', 12);
    xlabel('monitor power (re-zeroed, V)');
    ylabel('fluorescence / reference');
    if linear
        % A straight line through the origin is the whole model here, so draw it: the
        % point of the panel is that the data does NOT bend.
        s = sum(pw(:).*norm(:)) / max(sum(pw(:).^2), eps);
        Pf = linspace(0, max(pw)*1.05, 50);
        plot(Pf, s*Pf, '-');
        legend('data', sprintf('straight line, F = %.0f P', s), 'Location', 'southeast');
        title('linear regime: no saturation over the available power');
    else
        title('no fit: the monitor is not tracking delivered power');
    end

    subplot(2,1,2);
    errorbar(aomValues, norm, normE, '.-');
    if isfield(results, 'plateauAom') && isfield(results, 'recommendedAom')
        yl = ylim;
        pa = results.plateauAom;
        patch([pa(1) pa(2) pa(2) pa(1)], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.9 1], ...
              'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
        set(gca, 'Children', flipud(get(gca, 'Children')));   % band behind the data
        try
            xline(results.recommendedAom, '--', sprintf('use %.3f V', results.recommendedAom));
        catch
            plot(results.recommendedAom*[1 1], yl, '--', 'HandleVisibility', 'off');
        end
        ylim(yl);
    end
    xlabel(sprintf('AOM value (%s)', aom.units));
    ylabel('fluorescence / reference');
    title('this is the curve that IS reproducible -- shaded band is the flat top');
end

function results = laserVerdict(results, aom)
% Say plainly what to set, where it lives, and what the run means. Reached from all three
% exit paths, so the instruction does not depend on whether the fit happened to work.
    fprintf('\n============================================================\n');
    fprintf('WHAT TO DO WITH THIS\n');
    fprintf('============================================================\n');

    % ---- 1. the setting ----
    if isfield(results, 'recommendedAom') && isfinite(results.recommendedAom)
        v = results.recommendedAom;
        fprintf('\n1. SET THE AOM TO %.3f V. In code, in this session:\n\n', v);
        fprintf('       laser = getObjByName(LaserGate.GREEN_LASER_NAME);\n');
        fprintf('       laser.aom.value = %.3f;\n\n', v);
        fprintf(['   or move the Green Laser slider in the GUI to the same number. This is\n' ...
                 '   the ONE actionable output of level 0.2; everything else is diagnosis.\n']);
        if isfield(results, 'plateauAom') && diff(results.plateauAom) > 0
            fprintf('   Anywhere in %.3f to %.3f V is equivalent; the middle is the safe pick.\n', ...
                    results.plateauAom(1), results.plateauAom(2));
        end
    else
        fprintf('\n1. No usable brightest setting was found -- see the messages above.\n');
    end

    % ---- 2. where it does and does not live ----
    fprintf(['\n2. IT DOES NOT GO IN expProperties.json. There is no laserPower field there\n' ...
             '   and nothing reads one -- do not add it. The AOM value lives on the LaserGate\n' ...
             '   object, and the lab''s own save-state carries it (LaserGate stores it under\n' ...
             '   its aom_value property), so setting it once and saving the setup state is\n' ...
             '   what makes it stick. Set it at the top of a session otherwise.\n']);

    % ---- 3. do not run the AOM at full drive ----
    if isfield(results, 'brightestAom') && isfinite(results.brightestAom) && ...
       isfield(results, 'aomValues') && ~isempty(results.aomValues)
        [~, iMaxDrive] = max(results.aomValues);      % 0 V is maximum drive on Setup 7
        if isfield(results, 'normalised') && numel(results.normalised) >= iMaxDrive && ...
           isfinite(results.normalised(iMaxDrive)) && isfield(results, 'brightestNorm')
            lossPct = 100*(results.brightestNorm - results.normalised(iMaxDrive)) / ...
                      max(results.brightestNorm, eps);
            if lossPct > 2
                fprintf(['\n3. DO NOT RUN THE AOM AT %.3f V. Maximum drive is not maximum light:\n' ...
                         '   it gives %.0f%% LESS fluorescence than the plateau, because the AOM is\n' ...
                         '   past its own diffraction-efficiency peak. Anyone reaching for "full\n' ...
                         '   power" by pushing the knob to the end of its range loses %.0f%% of the\n' ...
                         '   signal and has no way to see it.\n'], ...
                        results.aomValues(iMaxDrive), lossPct, lossPct);
            end
        end
    end

    % ---- 4. what Psat did or did not say ----
    fprintf('\n4. Psat. ');
    if isfield(results, 'linearRegime') && results.linearRegime
        frac = NaN;
        if isfield(results, 'Psat') && isfinite(results.Psat) && results.Psat > 0 && ...
           isfield(results, 'power') && any(isfinite(results.power))
            frac = 100*max(results.power)/results.Psat;
        end
        fprintf(['NOT MEASURED, and that is a statement about the laser rather than\n' ...
                 '   about the fit. The response never bends over the power available, so the\n' ...
                 '   NV does not come close to saturating.']);
        if isfinite(frac)
            fprintf(' At full plateau power you are at roughly\n   %.0f%% of Psat.\n', frac);
        else
            fprintf('\n');
        end
        fprintf(['   Record NOTHING from the fit -- not Psat, not Fsat, not c. In the linear\n' ...
                 '   regime the saturating and linear terms are degenerate and which of them\n' ...
                 '   absorbs the line is arbitrary, which is why their error bars exceed their\n' ...
                 '   values while R^2 sits at 0.998.\n']);
        fprintf(['\n   The useful consequence: you are photon-starved, so more green at the\n' ...
                 '   sample would improve every SNR above this level, and polarisation is\n' ...
                 '   slower than it would otherwise be. Cross-check that against\n' ...
                 '   laserInitializationDuration in expProperties.json -- if it has had to be\n' ...
                 '   raised well above the 10 us the experiment classes default to, that is\n' ...
                 '   the same finding measured a second way, and the fix is more power (or\n' ...
                 '   better coupling), not a longer init.\n']);
    elseif isfield(results, 'Psat') && isfinite(results.Psat)
        fprintf('measured at %.6f monitor V. Use it as a scale, not a target.\n', results.Psat);
    else
        fprintf('not available on this run.\n');
    end

    % ---- 5. was the run itself sound ----
    fprintf('\n5. WAS THIS RUN SOUND?\n');
    if isfield(results, 'referenceSpread')
        fprintf('   reference moved %.1f%% across the scan  --  %s\n', 100*results.referenceSpread, ...
                pick(results.referenceSpread < 0.03, 'fine', 'high; let the laser settle and repeat'));
    end
    if isfield(results, 'driftSigma')
        fprintf('   end-of-scan repeat %+.1f sigma  --  %s\n', results.driftSigma, ...
                pick(abs(results.driftSigma) < 3, 'no significant drift', 'the setup moved; repeat'));
    end
    if isfield(results, 'monitorUsable')
        fprintf('   monitor usable as a power axis: %s\n', pick(results.monitorUsable, 'yes', 'no'));
    end
    fprintf(['\n   If all three are fine, level 0 is done: set the AOM and go to level 1.\n' ...
             '   Nothing from level 0 needs to be written to a file.\n']);
    fprintf('\n============================================================\n');
end

function s = pick(cond, a, b)
    if cond; s = a; else; s = b; end
end

function p = refSpreadPct(results)
% The reference's peak-to-peak movement as a PERCENT. Stored as a fraction, printed as a
% percent elsewhere, so convert in exactly one place rather than at each call site.
    p = 0;
    if isfield(results, 'referenceSpread') && isfinite(results.referenceSpread)
        p = 100 * results.referenceSpread;
    end
end

function r = rankOf(x)
    [~, idx] = sort(x);
    r = zeros(size(x));
    r(idx) = 1:numel(x);
end

function c = corrCoef(a, b)
    a = a(:) - mean(a(:));
    b = b(:) - mean(b(:));
    d = sqrt(sum(a.^2) * sum(b.^2));
    if d == 0; c = 0; else; c = sum(a.*b)/d; end
end
