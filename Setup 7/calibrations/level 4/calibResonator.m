function results = calibResonator(opts)
%CALIBRESONATOR Reflected microwave versus frequency, both channels driven together.
%
%   results = calibResonator()
%   results = calibResonator(struct('freqs', 2.6e9:5e6:3.0e9))
%
% LEVEL 3. Measures the microwave chain's frequency response from the circulators'
% reverse ports -- no NVs, no field, no laser.
%
% ===================================================================
% WHY BOTH CHANNELS AT ONCE
% ===================================================================
%
% Each circulator sends port 1 to port 2 and port 2 to port 3, so the amplifier's
% output goes to the antenna and whatever the antenna sends back comes out of port 3,
% where the 20 dB attenuator sits. Maximum transfer into the PCB therefore means
% MINIMUM at port 3: the resonance is a dip in what this measures.
%
% The first version drove one channel and left the other dark. That was wrong. The two
% PCBs sit either side of the same diamond and are mutually coupled, so the wave
% returning to circulator 1 is channel 1's own reflection PLUS whatever crosses over
% from channel 2. With channel 2 off, its amplifier presents some arbitrary off-state
% impedance and the boundary condition is not the one the experiments run in.
%
% So both channels are driven, at the balance and phase already calibrated, and both
% reverse ports are read simultaneously. What that yields is the ACTIVE reflection --
% how much power fails to get in when the pair is operating normally -- which is the
% quantity that matters, though it is not a one-port S11 and should not be called one.
%
% PHASE CAVEAT. channelPhase was calibrated at 2870 MHz. The electrical phase between
% the two paths drifts with frequency, so a fixed value is only correct near there.
% Across a several-hundred-MHz sweep the relative phase is uncontrolled: good enough
% to locate features, not to quote numbers far from 2870.
%
% ===================================================================
% WHAT IS AND IS NOT REPORTED
% ===================================================================
%
% Reported: the reflected amplitude per channel, each channel's reflection relative to
% its OWN worst point in dB, the summed reflected power, and where the minima are.
% All of that is model-free.
%
% NOT reported: |Gamma|, Q, or a flattening table. Earlier versions produced all
% three and every one was an artefact. Separating the antenna's reflection from the
% chain's gain needs a known reflection standard -- a short on the antenna port -- and
% the PCBs cannot be disconnected. Without it, a peak in the reading may be poor
% matching or may be the amplifier being strong there, and nothing in the data
% distinguishes them. Fitting a smooth "baseline" through the wings assumed a shape
% the data does not have: on Setup 7 there are two sharp peaks at 2.825 and 3.225 GHz,
% the baseline fit ran straight through them, and the result was |Gamma| > 1 clipped to
% 1 and a table of +156 dB corrections. A fitted Q of 12 with r^2 = 0.32 was the fit
% saying it did not believe itself.
%
% If a connector between amplifier and PCB ever becomes reachable, one sweep with a
% short there gives the reference and all of the above becomes available honestly.
%
% ===================================================================
% TWO SETTINGS THAT ARE NOT OPTIONAL
% ===================================================================
%
% SETTLE TIME. Three seconds after each frequency change, and it is not padding: at
% 0.3 s roughly half the points came back as flat traces. The chain needs seconds to
% settle after an FGEN frequency change.
%
% SYNCHRONISATION. Every :DIGitize is followed by *OPC?. :DIGitize does not reliably
% block on this instrument with an untriggered acquisition, and reading the waveform
% straight after is a race whose losses look like a flat trace with its spectral peak
% at DC. The per-point "hits" column counts how many repeats actually contained the
% tone; below 5/5 means the synchronisation is failing again.
%
% ===================================================================
% OPTIONS (struct fields, all optional)
%   channels       AWG channels to drive (default [1 2])
%   scopeChannels  scope inputs carrying their reverse ports (default [1 4])
%   scopeResource  VISA resource (default: the EXR/MXR/DSO found by visadevlist,
%                  USB preferred over LAN. Inputs are forced to 50 ohm and the
%                  setting is read back -- see setInputFifty.)
%   jsonLocation   expProperties.json, for amplitudeDiff and channelPhase
%   amplitudeDiff  dB per AWG channel, OVERRIDES the json when non-empty
%                  (default []). Element k is channel k, same as the json.
%   channelPhase   degrees per AWG channel, OVERRIDES the json when non-empty
%                  (default []). Use these two to sweep the same frequency grid
%                  at a deliberately wrong balance or phase: if the PCBs are
%                  coupled, the reflection minimum MOVES, and if it does not move
%                  then the crossover term is negligible and each reverse port is
%                  reading its own one-port reflection.
%   freqs          Hz (default: 25 MHz steps over 2.5-3.3 GHz, with 5 MHz steps
%                  through 2.79-2.91 and 3.14-3.26 where the sharp features are)
%   vpp            channel 1 amplitude, volts peak-to-peak (default 0.01)
%   nAvg           acquisitions per frequency, median-combined (default 5)
%   settle         seconds after a frequency change (default 3)
%   recordTime     seconds per acquisition (default 200e-9)
%   range          scope full-scale volts to start from (default 0.5)
%   maxVolts       abort threshold at the scope, volts peak (default 2)
%   peakTol        fractional tolerance for "the peak is at the drive" (default 0.03)
%   randomize      visit the frequencies in random order (default true). Costs
%                  nothing and decorrelates slow drift from frequency: in
%                  sequential order a warming amplifier looks like a real slope.
%   nSweeps        repeat the whole sweep this many times (default 1). More than 1
%                  gives an error bar per point that INCLUDES drift, which is what
%                  the Lorentzian fit needs to mean anything.
%   fitDip         locate the dip centre with a local parabola (default true)
%   fitSpan        Hz either side of the minimum included in that parabola
%                  (default 12e6 -- narrow on purpose: a parabola is only valid
%                  where nothing but the turning point is inside the window)
%
% Writes nothing. Restores the bench on exit, including on error.

    if ~exist('opts', 'var'); opts = struct(); end
    defFreqs = unique([2.5e9:25e6:3.3e9, 2.79e9:5e6:2.91e9, 3.14e9:5e6:3.26e9]);
    opts = withDefaults(opts, struct( ...
        'channels', [1 2], 'scopeChannels', [1 4], 'scopeResource', '', ...
        'jsonLocation', '', ...
        'freqs', defFreqs, 'vpp', 0.01, 'nAvg', 5, 'settle', 3, ...
        'recordTime', 200e-9, 'range', 0.5, 'maxVolts', 2, 'peakTol', 0.03, ...
        'randomize', true, 'nSweeps', 1, 'fitDip', true, 'fitSpan', 12e6, ...
        'amplitudeDiff', [], 'channelPhase', []));

    % ---- where the calibration lives ------------------------------------
    % Resolved by setup7Paths rather than written down here -- see setup7Paths for why a
    % path repeated in twelve files is worse than N times the work to change.
    %
    % TOLERANT ON PURPOSE. This function can do useful work with no json at all: a replay
    % or dry run never reads it, and the offline graders run with no Setup 7 tree on the
    % path. So a failure to resolve leaves jsonLocation EMPTY, which every use site below
    % already guards with ~isempty, and the calibration falls back to class defaults. An
    % eager error here broke testFieldPipeline, which fits saved spectra and touches no
    % hardware and no json.
    if isempty(opts.jsonLocation)
        try
            opts.jsonLocation = setup7Paths().json;
        catch
            opts.jsonLocation = '';
        end
    end

    ch  = opts.channels(:)';
    sch = opts.scopeChannels(:)';
    if numel(ch) ~= numel(sch)
        error('channels and scopeChannels must have the same length.');
    end
    f  = opts.freqs(:)';
    nF = numel(f);
    nC = numel(ch);

    % ---- calibrated balance and phase, from the json, or overridden --------
    cal = readCal(opts.jsonLocation);
    [ampDiff, cal] = overrideCal(cal, 'amplitudeDiff', opts.amplitudeDiff, 'dB',  max(ch));
    [phase,   cal] = overrideCal(cal, 'channelPhase',  opts.channelPhase,  'deg', max(ch));
    % amplitudeDiff is in dB of POWER, added to the commanded dBm by ExperimentAWG's
    % amplitude setter. In volts that is 10^(dB/20).
    vppCh = opts.vpp * 10.^(ampDiff(ch)/20);

    fprintf('\n=== calibResonator ===\n');
    fprintf('%d frequencies, %.3f to %.3f GHz\n', nF, f(1)/1e9, f(end)/1e9);
    for i = 1:nC
        fprintf('  AWG ch%d -> circulator -> scope ch%d : %.4g Vpp, phase %g deg\n', ...
                ch(i), sch(i), vppCh(i), phase(ch(i)));
    end
    fprintf('%s\n', cal.source);
    fprintf(['Both channels driven together: the PCBs are coupled, so channel 1''s reverse\n' ...
             '   port carries its own reflection PLUS the crossover from channel 2. Driving\n' ...
             '   one alone measures a boundary condition the experiments never use.\n']);
    fprintf(['Port 3 is a MINIMUM where transfer into the PCB is maximum, so the resonance\n' ...
             '   shows up as a dip.\n']);
    fprintf(['No |Gamma|, Q or flattening table: separating the antenna from the amplifier\n' ...
             '   needs a reflection standard and the PCBs cannot be disconnected. Everything\n' ...
             '   below is relative.\n']);

    % ---- instruments -------------------------------------------------------
    awg = findAwg();
    sc  = openScope(opts.scopeResource);
    pg  = getObjByName(PulseGenerator.NAME);
    cleanup = onCleanup(@() shutDown(awg, pg));                          %#ok<NASGU>

    setupScope(sc, sch, opts);

    % Both switches open. Whether one call can hold two channels is not something the
    % abstract PulseGenerator interface promises, so open each and then verify.
    mwNames = cell(1, nC);
    for i = 1:nC
        mwNames{i} = BooleanHelper.ifTrueElse(ch(i) == 1, 'MW', 'MW2');
        pg.on(mwNames{i});
    end
    for i = 1:nC
        try
            if ~pg.isOn(mwNames{i})
                fprintf(['WARNING switch "%s" reads OFF after pg.on. If only one switch can be\n' ...
                         '   held at a time, that channel contributes nothing and this is a\n' ...
                         '   single-channel measurement after all.\n'], mwNames{i});
            end
        catch
        end
    end
    fprintf('switches held open: %s\n', strjoin(mwNames, ', '));

    % FGEN mode, run and outputs established ONCE through awgFgenTone, which is where
    % the validated start ordering lives.
    awgFgenTone(struct('channels', ch, 'freq', f(1)*ones(1,nC), 'vpp', vppCh, ...
                       'type', 'SINE', 'phase', phase(ch)));

    % ---- sweep -------------------------------------------------------------
    Vall = nan(nC, nF, opts.nSweeps);      % median amplitude per point per sweep
    Wall = nan(nC, nF, opts.nSweeps);      % scatter WITHIN a point, across nAvg
    hits  = zeros(nC, nF);
    valid = true(nC, nF);
    Vrms  = nan(nC, nF);

    hdr = sprintf('%-6s %-10s', 'sweep', 'GHz');
    for i = 1:nC
        hdr = [hdr sprintf('%-14s%-8s', sprintf('ch%d mV pk', sch(i)), 'hits')]; %#ok<AGROW>
    end
    fprintf('\n%s %s\n', hdr, 'note');

    for sw = 1:opts.nSweeps
        % Randomised visit order. Slow drift -- an amplifier warming up, the AWG
        % settling -- is indistinguishable from a real frequency dependence when the
        % sweep runs monotonically, and it costs nothing to remove that confusion.
        order = 1:nF;
        if opts.randomize; order = randperm(nF); end

        for kk = order
            st = setTone(awg, ch, f(kk), vppCh, phase(ch));
            pause(opts.settle);
            d = acquireAll(sc, sch, f(kk), opts);

            ln = sprintf('%-6d %-10.4f', sw, f(kk)/1e9);
            note = '';
            if ~st.ok
                note = sprintf('[GENERATOR state %s, outputs %s]', st.runState, st.outState);
            end
            for i = 1:nC
                Vall(i,kk,sw) = d(i).A;
                Wall(i,kk,sw) = d(i).Astd;
                Vrms(i,kk)    = d(i).vrms;
                hits(i,kk)    = d(i).hits;
                okPt = st.ok && d(i).hits >= ceil(0.6*opts.nAvg) && ...
                       abs(d(i).peakFreq - f(kk)) <= opts.peakTol*f(kk) && ...
                       d(i).samplesPerCycle >= 2.5;
                valid(i,kk) = valid(i,kk) && okPt;
                ln = [ln sprintf('%-14.4g%-8s', 1e3*d(i).A, ...
                      sprintf('%d/%d', d(i).hits, opts.nAvg))];               %#ok<AGROW>
                if d(i).fill > 0.85
                    note = [note sprintf(' [ch%d near clipping]', sch(i))];   %#ok<AGROW>
                end
                if d(i).A > opts.maxVolts
                    error(['Scope channel %d reads %.2f V, over the %.2f V abort threshold. ' ...
                           'Reduce opts.vpp or add attenuation.'], sch(i), d(i).A, opts.maxVolts);
                end
            end
            if sw == 1 && kk == order(1)
                fprintf('   [%.3g GS/s, %.1f samples/cycle, %d points]\n', ...
                        d(1).sampleRate/1e9, d(1).samplesPerCycle, d(1).n);
            end
            fprintf('%s %s\n', ln, note);
        end
    end

    % ---- combine, with an error bar that does not lie ----------------------
    % Two independent estimates, and the larger wins. The WITHIN estimate is the
    % scatter of the nAvg acquisitions at one setting, which is electronic noise. The
    % ACROSS estimate is the scatter of the sweep medians, which also contains drift.
    % Quoting only the within estimate assumes the chain is stationary over the tens
    % of minutes a multi-sweep run takes, and it is not -- the same reasoning as
    % pdReadRepeated at Level 0.
    V = median(Vall, 3, 'omitnan');
    Wwithin = median(Wall, 3, 'omitnan') / sqrt(max(opts.nAvg, 1));
    if opts.nSweeps > 1
        Wacross = std(Vall, 0, 3, 'omitnan') / sqrt(opts.nSweeps);
    else
        Wacross = nan(size(V));
    end
    Verr = max(Wwithin, Wacross);
    Verr(~isfinite(Verr)) = Wwithin(~isfinite(Verr));

    results = struct('freqs', f(:), 'channels', ch, 'scopeChannels', sch, ...
                     'amplitude', V, 'amplitudeErr', Verr, 'amplitudeAll', Vall, ...
                     'vrms', Vrms, 'hits', hits, 'valid', valid, ...
                     'errWithin', Wwithin, 'errAcross', Wacross, ...
                     'vppPerChannel', vppCh, 'phase', phase(ch), 'opts', opts, 'cal', cal);

    % ---- model-free reading ------------------------------------------------
    fprintf('\n--- coverage ---\n');
    for i = 1:nC
        fprintf('  scope ch%d: %d of %d points contain the tone\n', sch(i), sum(valid(i,:)), nF);
    end
    if ~any(valid(:))
        error('No point on any channel contained the tone. Run scopeProbe.');
    end

    % Reflected POWER in detector units. Summing across channels assumes the two
    % detection paths have the same attenuation and response, which they do not exactly
    % -- both are nominally 20 dB but nothing here measures that. So the sum is
    % indicative of total loss; the per-channel curves are the trustworthy part.
    P = V.^2;
    Ptot = sum(P, 1, 'omitnan');
    results.reflectedPower = P;
    results.reflectedPowerTotal = Ptot(:);
    results.minFreq = nan(1, nC);

    fprintf('\n--- where reflection is lowest (best transfer into the PCB) ---\n');
    for i = 1:nC
        vi = valid(i,:);
        if ~any(vi); continue; end
        iv = find(vi);
        [~, im] = min(V(i,vi)); im = iv(im);
        [~, iM] = max(V(i,vi)); iM = iv(iM);
        results.minFreq(i) = f(im);
        fprintf('  ch%d: min at %.4f GHz (%.4g mV), worst at %.4f GHz (%.4g mV), %.1f dB apart\n', ...
                sch(i), f(im)/1e9, 1e3*V(i,im), f(iM)/1e9, 1e3*V(i,iM), ...
                20*log10(V(i,iM)/max(V(i,im), eps)));
        % A minimum defined by one point is a sampling accident until it is resolved,
        % so show its neighbours rather than quoting the frequency alone.
        if im > 1 && im < nF
            step = min(f(im)-f(im-1), f(im+1)-f(im));
            fprintf('        (local step %.0f MHz; neighbours %.4g and %.4g mV)\n', ...
                    step/1e6, 1e3*V(i,im-1), 1e3*V(i,im+1));
        end
    end
    vt = all(valid, 1);
    if any(vt)
        ivt = find(vt);
        [~, it] = min(Ptot(vt)); it = ivt(it);
        fprintf('  combined reflected power is lowest at %.4f GHz\n', f(it)/1e9);
        results.minFreqCombined = f(it);
    end

    % ---- centre of the dip, by parabola ------------------------------------
    % A PARABOLA over a narrow window, not a Lorentzian.
    %
    % The Lorentzian was wrong and said so: reduced chi2 of 16662 on channel 1, an
    % error bar on the centre of 240 MHz for a dip 20 MHz wide, and a fitted FWHM of
    % 43.7 +/- 1033.9 MHz. The reason is that over +/-80 MHz this data is not one
    % Lorentzian on a smooth background -- channel 1 has a second broad dip near 2.71
    % and a shoulder around 2.80, and forcing a single-mode model across all of it
    % produces a fit dominated by the parts that do not belong to the mode.
    %
    % A parabola makes no claim about the shape. Every smooth minimum is locally
    % quadratic, so over a window narrow enough that only the turning point is inside
    % it, the vertex IS the centre, and the fit has three parameters instead of five.
    % What a parabola cannot give is a width -- it has none -- so the width below is
    % measured directly from the data by interpolation rather than fitted.
    results.fit = struct([]);
    if opts.fitDip
        fprintf('\n--- dip centre by parabola (window +/-%.0f MHz) ---\n', opts.fitSpan/1e6);
        for i = 1:nC
            vi = valid(i,:);
            if sum(vi) < 5; continue; end
            iv = find(vi);
            [~, im] = min(V(i,vi)); im = iv(im);
            win = vi & abs(f - f(im)) <= opts.fitSpan;
            if sum(win) < 4
                fprintf('  ch%d: only %d points within the window, not fitting\n', ...
                        sch(i), sum(win));
                continue
            end
            x = (f(win) - f(im))/1e6;              % MHz, for conditioning
            y = V(i,win).^2;                       % power
            ey = 2 * V(i,win) .* max(Verr(i,win), eps);
            w  = 1 ./ max(ey, eps).^2;

            % Weighted quadratic, closed form. Nothing to converge, no start point.
            X = [ones(numel(x),1), x(:), x(:).^2];
            W = diag(w(:));
            b = (X'*W*X) \ (X'*W*y(:));
            r = y(:) - X*b;
            dof = max(numel(x) - 3, 1);
            chi2 = (r'*W*r)/dof;
            C = inv(X'*W*X) * max(chi2, 0.09);     % two-sided Birge scaling, floored

            if b(3) <= 0
                fprintf('  ch%d: window curves the wrong way, no minimum inside it\n', sch(i));
                continue
            end
            xv  = -b(2)/(2*b(3));                  % vertex, MHz from f(im)
            g2  = -1/(2*b(3)); g3 = b(2)/(2*b(3)^2);
            xve = sqrt(max(g2^2*C(2,2) + g3^2*C(3,3) + 2*g2*g3*C(2,3), 0));
            f0  = f(im) + xv*1e6;
            f0e = xve*1e6;

            % Width WITHOUT a model: interpolate where the power rises to halfway
            % between the minimum and the nearest local maximum on each side.
            [wLo, wHi, fwhmMeas] = halfDepthWidth(f, V(i,:).^2, vi, im);

            fprintf('  ch%d: centre %.4f +/- %.4f GHz  (%d points, reduced chi2 %.2f)\n', ...
                    sch(i), f0/1e9, f0e/1e9, sum(win), chi2);
            if isfinite(fwhmMeas)
                fprintf('        measured half-depth width %.1f MHz (%.4f to %.4f GHz), Q ~ %.0f\n', ...
                        fwhmMeas/1e6, wLo/1e9, wHi/1e9, f0/fwhmMeas);
            else
                fprintf('        half-depth width not measurable inside the swept range\n');
            end
            if chi2 > 5
                fprintf(['        the window is still too wide for a parabola -- reduce\n' ...
                         '        opts.fitSpan until this falls near 1\n']);
            end

            results.fit(i).channel = sch(i);
            results.fit(i).f0 = f0;      results.fit(i).f0Err = f0e;
            results.fit(i).chi2 = chi2;
            results.fit(i).coeff = b(:)';
            results.fit(i).xOrigin = f(im);
            results.fit(i).window = f(win);
            results.fit(i).widthMeasured = fwhmMeas;
            results.fit(i).widthEdges = [wLo wHi];
        end
    end

    fprintf(['\n   Read these as RELATIVE. A peak may be poor matching or may be the amplifier\n' ...
             '   being strong there; without a reflection standard nothing in this data\n' ...
             '   separates the two. The frequency of the MINIMUM is the useful output, because\n' ...
             '   both explanations put the best transfer at the same place.\n']);
    fprintf('   Nothing was written.\n\n');

    % ---- plots -------------------------------------------------------------
    figure
    subplot(2,1,1); hold on
    for i = 1:nC
        errorbar(f/1e9, 1e3*V(i,:), 1e3*Verr(i,:), 'o-');
    end
    for i = 1:numel(results.fit)
        if ~isfield(results.fit, 'coeff') || isempty(results.fit(i).coeff); continue; end
        wf = results.fit(i).window;
        xr = linspace(min(wf), max(wf), 200);
        xm = (xr - results.fit(i).xOrigin)/1e6;
        b  = results.fit(i).coeff;
        yp = b(1) + b(2)*xm + b(3)*xm.^2;
        plot(xr/1e9, 1e3*sqrt(max(yp, 0)), 'k-', 'LineWidth', 1.2);
        line([results.fit(i).f0 results.fit(i).f0]/1e9, ylim, 'Color', [0.4 0.4 0.4], ...
             'LineStyle', ':');
    end
    ylabel('reflected (mV pk)');
    legend(arrayfun(@(c) sprintf('scope ch%d', c), sch, 'UniformOutput', false), ...
           'Location', 'best');
    title('reflected amplitude -- low means power is getting into the PCB');

    subplot(2,1,2); hold on
    for i = 1:nC
        ref = max(V(i, valid(i,:)));
        plot(f/1e9, 20*log10(V(i,:)/max(ref, eps)), 'o-');
    end
    xlabel('frequency (GHz)'); ylabel('dB below that channel''s worst');
    title('relative only -- no reflection standard, so this is not |\Gamma|');
end

% ======================================================================

function cal = readCal(jsonLocation)
% amplitudeDiff and channelPhase out of expProperties.json. Both default to zero with
% a loud note rather than a silent assumption, because running with the wrong balance
% or phase measures a configuration the experiments do not use.
    cal = struct('amplitudeDiff', [0 0], 'channelPhase', [0 0], ...
                 'source', 'defaults [0 0] and [0 0] -- json NOT read');
    if isempty(jsonLocation); return; end
    try
        j = jsondecode(fileread(jsonLocation));
    catch err
        fprintf('WARNING could not read %s\n  (%s)\n', jsonLocation, err.message);
        return
    end
    got = {};
    if isfield(j, 'amplitudeDiff')
        cal.amplitudeDiff = j.amplitudeDiff(:)';
        got{end+1} = sprintf('amplitudeDiff %s dB', mat2str(cal.amplitudeDiff));
    else
        got{end+1} = 'amplitudeDiff absent -- using [0 0]';
    end
    if isfield(j, 'channelPhase')
        cal.channelPhase = j.channelPhase(:)';
        got{end+1} = sprintf('channelPhase %s deg', mat2str(cal.channelPhase));
    else
        got{end+1} = 'channelPhase ABSENT -- using [0 0]';
    end
    if isfield(j, 'timeDelay') && j.timeDelay ~= 0
        got{end+1} = sprintf(['WARNING timeDelay = %g us is still non-zero. It is the other, ' ...
                              'coarsely quantised phase mechanism and it ADDS to channelPhase'], ...
                             j.timeDelay);
    end
    cal.source = ['from the json: ' strjoin(got, '; ')];
end

function [v, cal] = overrideCal(cal, field, override, unit, maxCh)
% opts.amplitudeDiff / opts.channelPhase beat the json when non-empty. Same
% convention as the json: element k is AWG channel k, so [0 -1.4] means channel 1
% at reference and channel 2 1.4 dB down. Missing trailing entries pad with zero.
%
% The point of the override is to answer "does the reflection minimum move when the
% relative phase moves?" -- if the two PCBs are genuinely coupled it must, and if it
% does not then each port is seeing its own one-port reflection and the crossover
% term is negligible. Sweeping at channelPhase [0 26.4] and again at [0 116.4] is
% that test.
    v = cal.(field);
    if ~isempty(override)
        v = override(:)';
        cal.(field) = v;
        cal.source  = [cal.source sprintf('\n  OVERRIDDEN %s = %s %s (opts, not the json)', ...
                                          field, mat2str(v), unit)];
    end
    if numel(v) < maxCh; v(end+1:maxCh) = 0; end
    cal.(field) = v;
end

function st = setTone(awg, chans, freqHz, vppCh, phaseCh)
% New frequency on every channel, then the full bring-up, unconditionally.
%
% Not conditional on a state query, because this instrument reports a STALE run state
% for about a second in BOTH directions: 0 for a second after a successful start, and
% still 2 for a while after a parameter change has stopped it. A version that
% restarted only when RSTATE? read 0 therefore never restarted, and the whole sweep
% measured an unpowered chain.
    cmds = {};
    for i = 1:numel(chans)
        c = chans(i);
        cmds{end+1} = sprintf('FGEN:CHAN%d:TYPE SINE', c);                 %#ok<AGROW>
        cmds{end+1} = sprintf('FGEN:CHAN%d:FREQ %.9g', c, freqHz);         %#ok<AGROW>
        cmds{end+1} = sprintf('FGEN:CHAN%d:AMPLITUDE %.6g', c, vppCh(i));  %#ok<AGROW>
        cmds{end+1} = sprintf('FGEN:CHAN%d:PHASE %g', c, phaseCh(i));      %#ok<AGROW>
    end
    awg.sendCommandBatch(cmds);

    st = struct('runState', '0', 'outState', '', 'ok', false);
    awg.sendCommand('awgcontrol:run:immediate');
    t0 = tic;
    while toc(t0) < 10
        st.runState = strtrim(awg.queryRaw('AWGCONTROL:RSTATE?'));
        if ~strcmp(st.runState, '0'); break; end
        pause(0.05);
    end
    awg.workChannels = chans;                 % outputs on, AFTER the run
    outs = cell(1, numel(chans)); allOn = true;
    for i = 1:numel(chans)
        outs{i} = strtrim(awg.queryRaw(sprintf('OUTPUT%d:STATE?', chans(i))));
        allOn = allOn && any(strcmp(outs{i}, {'1','ON'}));
    end
    st.outState = strjoin(outs, ',');
    st.ok = ~strcmp(st.runState, '0') && allOn;
end

function d = acquireAll(sc, sch, freqHz, opts)
% Range each scope channel, then nAvg simultaneous acquisitions of all of them.
    q = [];
    for i = 1:numel(sch)
        for attempt = 1:5
            rng = str2double(query1(sc, sprintf(':CHANnel%d:RANGe?', sch(i))));
            q = acquireOnce(sc, sch, freqHz);
            % Every range change goes through setRange, which reports what the
            % instrument actually accepted. The EXR clamps the 50 ohm input well below
            % the 8 V the old scope allowed, so comparing the REQUEST against the old
            % value would loop forever against a ceiling; comparing the ACCEPTED value
            % exits as soon as the scope stops moving.
            if q(i).fill > 0.85
                newRng = setRange(sc, sch(i), min(rng*2.5, 8));
                if newRng <= rng*1.01; break; end
            elseif q(i).fill < 0.06 && rng > 0.05
                newRng = setRange(sc, sch(i), max(rng/2.5, 0.008));
                if newRng >= rng*0.99; break; end
            else
                break
            end
        end
    end

    A = nan(numel(sch), opts.nAvg); vr = A; pk = A; fl = A;
    for j = 1:opts.nAvg
        q = acquireOnce(sc, sch, freqHz);
        for i = 1:numel(sch)
            A(i,j) = q(i).A; vr(i,j) = q(i).vrms;
            pk(i,j) = q(i).peakFreq; fl(i,j) = q(i).fill;
        end
    end
    d = struct('A', {}, 'Astd', {}, 'vrms', {}, 'peakFreq', {}, 'fill', {}, ...
               'hits', {}, 'sampleRate', {}, 'samplesPerCycle', {}, 'n', {});
    for i = 1:numel(sch)
        d(i).A = median(A(i,:), 'omitnan');
        % Robust scatter of the repeats, converted to a standard-deviation equivalent.
        % MAD rather than std so one glitched acquisition cannot inflate the error bar
        % that the Lorentzian fit will be weighted by.
        d(i).Astd = 1.4826 * median(abs(A(i,:) - d(i).A), 'omitnan');
        d(i).vrms = median(vr(i,:), 'omitnan');
        d(i).peakFreq = median(pk(i,:), 'omitnan');
        d(i).fill = median(fl(i,:), 'omitnan');
        d(i).hits = sum(abs(pk(i,:) - freqHz) <= opts.peakTol*freqHz);
        d(i).sampleRate = q(i).sampleRate;
        d(i).samplesPerCycle = q(i).samplesPerCycle;
        d(i).n = q(i).n;
    end
end

function q = acquireOnce(sc, sch, freqHz)
% Digitise ALL channels in one command so they are simultaneous, wait for completion,
% then read each. Simultaneity matters: the two reverse ports are being compared with
% each other, and acquiring them at different moments would fold any drift between
% them into the comparison.
    lst = strjoin(arrayfun(@(c) sprintf('CHANnel%d', c), sch, 'UniformOutput', false), ',');
    writeline(sc, [':DIGitize ' lst]);
    try
        query1(sc, '*OPC?');            % blocks until the acquisition is complete
    catch
        pause(0.2);
    end
    q = struct('vrms', {}, 'A', {}, 'fill', {}, 'sampleRate', {}, ...
               'samplesPerCycle', {}, 'n', {}, 'peakFreq', {});
    for i = 1:numel(sch)
        c = sch(i);
        % Fill fraction is computed from VOLTS against the channel's full-scale range,
        % not from raw ADC counts. The old code divided by 32768, which assumed the
        % DSO90404A's 8-bit ADC mapped to the full signed-16-bit span. The EXR is a
        % 10-bit part and can run in high-resolution modes deeper than that, so the
        % count-based fraction would be wrong by whatever factor the mapping happens
        % to use -- and it feeds the auto-ranging, so a wrong fill means a wrong
        % vertical scale on every point.
        rngV = NaN;
        try; rngV = str2double(query1(sc, sprintf(':CHANnel%d:RANGe?', c))); catch; end
        vrms = NaN;
        if vrmsSupported()
            vv = str2double(queryOpt(sc, sprintf(':MEASure:VRMS? DISPlay,AC,CHANnel%d', c)));
            if isfinite(vv) && abs(vv) < 1e30; vrms = vv; end
        end
        writeline(sc, sprintf(':WAVeform:SOURce CHANnel%d', c));
        pre = str2double(strsplit(strtrim(query1(sc, ':WAVeform:PREamble?')), ','));
        xinc = pre(5); xorg = pre(6); xref = pre(7);
        yinc = pre(8); yorg = pre(9); yref = pre(10);
        writeline(sc, ':WAVeform:DATA?');
        raw = double(readbinblock(sc, 'int16'));
        try; readline(sc); catch; end
        v = (raw - yref)*yinc + yorg;
        n = numel(v);
        t = ((0:n-1) - xref)*xinc + xorg;
        w = 0.5*(1 - cos(2*pi*(0:n-1)'/(n-1)));
        q(i).A = 2*abs(sum(v(:).*w.*exp(-2i*pi*freqHz*t(:))))/sum(w);
        if ~isfinite(vrms); vrms = sqrt(mean(v.^2)); end
        q(i).vrms = vrms;
        if isfinite(rngV) && rngV > 0
            q(i).fill = 2*max(abs(v))/rngV;
        else
            q(i).fill = max(abs(raw))/32768;    % last resort if RANGe? did not answer
        end
        q(i).sampleRate = 1/xinc;
        q(i).samplesPerCycle = q(i).sampleRate/freqHz;
        q(i).n = n;
        Vf = abs(fft(v(:).*w)); Vf = Vf(1:floor(n/2));
        [~, ip] = max(Vf(2:end)); ip = ip + 1;
        q(i).peakFreq = (ip-1)*(q(i).sampleRate/n);
    end
end

function setupScope(sc, sch, opts)
    cmds = { ':SYSTem:HEADer OFF', ':ACQuire:MODE RTIMe', ':ACQuire:AVERage OFF', ...
             ':ACQuire:COMPlete 100', ':TRIGger:SWEep AUTO', ...
             sprintf(':TIMebase:SCALe %g', opts.recordTime/10), ...
             ':TIMebase:POSition 0', ':WAVeform:FORMat WORD', ...
             ':WAVeform:BYTeorder LSBFirst' };
    for i = 1:numel(cmds); writeline(sc, cmds{i}); end
    for c = sch
        writeline(sc, sprintf(':CHANnel%d:DISPlay ON', c));
        setInputFifty(sc, c);                  % before the range: the two interact
        writeline(sc, sprintf(':CHANnel%d:OFFSet 0', c));
        setRange(sc, c, opts.range);
    end
    % Checked, so a mnemonic the EXR does not have is named here rather than surfacing as
    % an anonymous -113 in the error queue three hundred acquisitions later.
    refused = {};
    opt = { ':ACQuire:SRATe:ANALog MAX'
            ':ACQuire:POINts:ANALog AUTO'
            ':ACQuire:INTerpolate OFF'
            ':WAVeform:STReaming OFF' };
    % :WAVeform:UNSigned is NOT in that list. The EXR rejects it (-113, no such setting)
    % and it is not needed: scopeProbe's ASCII cross-check read 70.66 mV against the WORD
    % path's 70.66 mV, identical to four figures, which settles the int16 interpretation
    % and the preamble scaling by measurement instead of by command.
    for i = 1:numel(opt)
        if ~sendChecked(sc, opt{i}); refused{end+1} = opt{i}; end   %#ok<AGROW>
    end
    if ~isempty(refused)
        fprintf('  (this instrument refused: %s)\n', strjoin(refused, ', '));
    end

    % Probe the RMS measurement ONCE. See vrmsSupported for why this cannot be left to
    % fail per acquisition.
    vrmsSupported(true);
    probe = queryOpt(sc, sprintf(':MEASure:VRMS? DISPlay,AC,CHANnel%d', sch(1)));
    if strcmp(probe, '(not supported)') || strcmp(probe, '(empty response)') || ...
       ~isfinite(str2double(probe))
        vrmsSupported(false);
        fprintf(['  :MEASure:VRMS? not usable on this instrument -- RMS will be computed\n' ...
                 '   from the trace instead. No effect on the measurement.\n']);
    end

    % Sample rate is not assumed from the datasheet, it is asked for. The EXR404A is
    % 16 GSa/s, which at 2.87 GHz is 5.6 samples per cycle -- comfortable. But some
    % scopes halve the rate as more channels are enabled, and this is the run where
    % channel count goes up, so the number is measured and complained about rather
    % than trusted.
    sr = NaN;
    try; sr = str2double(query1(sc, ':ACQuire:SRATe:ANALog?')); catch; end
    if isfinite(sr) && sr > 0
        fMax = max(opts.freqs);
        spc  = sr / fMax;
        fprintf('  sample rate %.3g GSa/s -> %.2f samples/cycle at %.3f GHz\n', ...
                sr/1e9, spc, fMax/1e9);
        if spc < 2.5
            warning(['calibResonator:sampleRate', newline, ...
                     'Only %.2f samples per cycle at the top of the sweep. Below about ' ...
                     '2.5 the tone is too close to Nyquist for the DFT amplitude to be ' ...
                     'trustworthy. Enable fewer scope channels or shorten the frequency ' ...
                     'range.'], spc);
        end
    end
    % Drain quietly. The queue usually holds stale "please turn on a waveform of the
    % type required" errors from measurement queries issued before a channel was
    % displayed -- twenty of them appeared in one run and none of them mattered.
    for k = 1:40
        r = strtrim(query1(sc, ':SYSTem:ERRor? STRing'));
        if isempty(r) || startsWith(r, '0') || startsWith(r, '+0'); break; end
    end
end

function trySend(sc, cmd)
    try; writeline(sc, cmd); catch; fprintf('  (scope did not accept: %s)\n', cmd); end
end

function s = query1(sc, cmd)
    writeline(sc, cmd);
    s = char(readline(sc));
end

function ok = sendChecked(sc, cmd)
% Write a command and immediately ask whether the instrument objected. Returns false if
% it did, naming the command. Attribution matters: an error queue drained at the end of
% a 101-point sweep says something was refused, not which line to fix.
    ok = true;
    try
        writeline(sc, cmd);
    catch
        fprintf('  (could not send: %s)\n', cmd);
        ok = false; return
    end
    e = queryOpt(sc, ':SYSTem:ERRor? STRing');
    if ~isempty(e) && ~startsWith(e, '0') && ~startsWith(e, '+0') && ...
       ~strcmp(e, '(not supported)') && ~strcmp(e, '(empty response)')
        fprintf('  REFUSED %s\n     -> %s\n', cmd, e);
        ok = false;
    end
end

function tf = vrmsSupported(setTo)
% Whether :MEASure:VRMS? works on this instrument, decided ONCE and remembered.
%
% This exists to prevent a specific disaster. The measurement query lives inside
% acquireOnce, which runs nAvg times per frequency per sweep. If the instrument does
% not recognise the mnemonic it returns nothing and the read blocks for the full VISA
% timeout -- 20 s, times 5 repeats, times 2 channels, times 101 frequencies, is over
% five hours of pure waiting for a number that is only ever printed. So it is probed
% once during setup and switched off for the rest of the run if it fails; the RMS is
% then computed from the trace, which is what the fallback always did anyway.
    persistent flag
    if nargin > 0; flag = setTo; end
    if isempty(flag); flag = true; end
    tf = flag;
end

function s = queryOpt(sc, cmd)
% A query the instrument is allowed to refuse. An unsupported mnemonic on the Infiniium
% family produces no response at all -- the command goes to the error queue and the read
% blocks -- so the timeout drops to 2 s for the attempt and the connection is flushed
% afterwards. Without the flush the next legitimate read would collect this query's
% late answer and every reading after it would be shifted by one.
    old = sc.Timeout;
    s = '(not supported)';
    try
        sc.Timeout = 2;
        writeline(sc, cmd);
        s = strtrim(char(readline(sc)));
        if isempty(s); s = '(empty response)'; end
    catch
        try; flush(sc); catch; end
        try
            for k = 1:20
                writeline(sc, ':SYSTem:ERRor? STRing');
                r = strtrim(char(readline(sc)));
                if isempty(r) || startsWith(r, '0') || startsWith(r, '+0'); break; end
            end
        catch
        end
    end
    try; sc.Timeout = old; catch; end
end

function sc = openScope(resource)
% Keysight EXR404A (and the MXR siblings) run the Infiniium software, so the SCPI
% below is the same dialect the DSO90404A used. Model matching covers both families
% and prefers USB over LAN when the same instrument appears twice in visadevlist.
    if isempty(resource)
        L = visadevlist;
        M = upper(string(L.Model));
        hit = startsWith(M, "EXR") | startsWith(M, "MXR") | contains(M, "DSO");
        if ~any(hit)
            error(['No Keysight scope found in visadevlist (looked for EXR*, MXR*, ' ...
                   '*DSO*). Pass opts.scopeResource explicitly.']);
        end
        idx = find(hit);
        isUsb = startsWith(upper(string(L.ResourceName(idx))), "USB");
        if any(isUsb); idx = idx(find(isUsb, 1)); else; idx = idx(1); end
        resource = L.ResourceName(idx);
        fprintf('scope resource: %s (%s)\n', resource, L.Model(idx));
    end
    sc = visadev(resource);
    sc.Timeout = 20;
    fprintf('scope: %s\n', strtrim(query1(sc, '*IDN?')));
end

function ok = setInputFifty(sc, c)
% 50 ohm termination, VERIFIED. The EXR inputs are switchable 50 ohm / 1 Mohm, which
% the DSO90404A was not -- it was 50 ohm and nothing else, so this never had to be
% set. It does now, and it is not a nicety: the 1 Mohm path is bandwidth limited to a
% few hundred MHz, so a 2.87 GHz tone into a 1 Mohm input reads as almost nothing.
% That failure looks exactly like the microvolt readings we already chased once, and
% it would be silent. Hence the readback rather than a fire-and-forget write.
%
% Two command spellings are tried because the Infiniium and InfiniiVision families
% disagree: :CHANnel<n>:INPut DC50 versus :CHANnel<n>:IMPedance FIFTy. Whichever the
% instrument actually honours is the one kept.
    ok = false;
    attempts = { sprintf(':CHANnel%d:INPut DC50', c),      sprintf(':CHANnel%d:INPut?', c); ...
                 sprintf(':CHANnel%d:IMPedance FIFTy', c), sprintf(':CHANnel%d:IMPedance?', c) };
    for k = 1:size(attempts, 1)
        % queryOpt, not query1: an instrument that rejects the first spelling returns
        % nothing at all, so a plain readback is a full 20 s VISA stall per attempt.
        % The EXR takes :INPut DC50 on the first try and never reaches the second.
        try; writeline(sc, attempts{k,1}); catch; end
        r = upper(queryOpt(sc, attempts{k,2}));
        if contains(r, "DC50") || contains(r, "DCF") || contains(r, "FIFT") || contains(r, "50")
            fprintf('  scope ch%d input: %s (50 ohm confirmed)\n', c, r);
            ok = true; return
        end
    end
    warning(['calibResonator:inputImpedance', newline, ...
             'Could not confirm 50 ohm on scope channel %d. If this input is on 1 Mohm ' ...
             'the 2.87 GHz tone is being rolled off by the input path and every ' ...
             'amplitude in this sweep is meaningless. Set it by hand on the front panel ' ...
             'before trusting anything.'], c);
end

function actual = setRange(sc, c, want)
% Write the full-scale range and read back what the instrument actually took. The EXR
% clamps the 50 ohm input to a smaller maximum range than the old scope allowed, so
% asking for 8 V silently lands somewhere else; the ranging loop needs to know the
% real value or it will keep trying to climb past a ceiling that is already reached.
    try
        writeline(sc, sprintf(':CHANnel%d:RANGe %g', c, want));
        actual = str2double(query1(sc, sprintf(':CHANnel%d:RANGe?', c)));
    catch
        actual = NaN;
    end
    if ~isfinite(actual); actual = want; end
end

function shutDown(awg, pg)
% Microwave off, switches closed, AWG back in the mode the experiments expect. Runs
% even if the sweep errors, which matters because an abort would otherwise leave both
% amplifiers driven.
    try; pg.off; catch; end
    try; awg.workChannels = []; catch; end
    try; awg.setValue('mode', 'AWG'); awg.invalidateLoaded; catch; end
    fprintf('\n(bench restored: switches off, AWG outputs off, AWG mode)\n');
end

function awg = findAwg()
    try
        fgs = FrequencyGenerator.getFG();
        for i = 1:numel(fgs)
            if contains(fgs{i}.name, 'TektronixAWG') || contains(fgs{i}.name, 'AWGdummy')
                awg = fgs{i}; return
            end
        end
    catch
    end
    error('No AWG found via FrequencyGenerator.getFG().');
end

function [fLo, fHi, width] = halfDepthWidth(f, P, valid, im)
% Width at half depth, measured by interpolation rather than fitted.
%
% "Half depth" is relative to the nearest local MAXIMUM on each side, not to some
% global baseline -- there is no trustworthy global baseline here, and the shoulders
% are what the dip actually sits between. Walk outward from the minimum until the
% power stops rising, take that as the local shoulder, and interpolate where the
% curve crosses halfway between the minimum and the shoulder.
    fLo = NaN; fHi = NaN; width = NaN;
    n = numel(f);
    P(~valid) = NaN;
    if ~isfinite(P(im)); return; end

    % walk left to the first turning point
    iL = im;
    while iL > 1 && (~isfinite(P(iL-1)) || P(iL-1) >= P(iL)); iL = iL - 1; if iL == 1; break; end; end
    iR = im;
    while iR < n && (~isfinite(P(iR+1)) || P(iR+1) >= P(iR)); iR = iR + 1; if iR == n; break; end; end
    if iL == im || iR == im; return; end

    halfL = (P(im) + P(iL))/2;
    halfR = (P(im) + P(iR))/2;
    fLo = crossBetween(f, P, im, iL, halfL);
    fHi = crossBetween(f, P, im, iR, halfR);
    if isfinite(fLo) && isfinite(fHi); width = abs(fHi - fLo); end
end

function fx = crossBetween(f, P, iFrom, iTo, level)
% First crossing of 'level' walking from iFrom towards iTo, linearly interpolated.
    fx = NaN;
    stepDir = sign(iTo - iFrom);
    k = iFrom;
    while k ~= iTo
        k2 = k + stepDir;
        if isfinite(P(k)) && isfinite(P(k2)) && ...
                ((P(k) - level)*(P(k2) - level) <= 0) && P(k2) ~= P(k)
            fr = (level - P(k))/(P(k2) - P(k));
            fx = f(k) + fr*(f(k2) - f(k));
            return
        end
        k = k2;
    end
end

function o = withDefaults(o, d)
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(o, fn{i}) || isempty(o.(fn{i})); o.(fn{i}) = d.(fn{i}); end
    end
end
