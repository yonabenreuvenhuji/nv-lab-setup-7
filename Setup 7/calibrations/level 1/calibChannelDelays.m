function results = calibChannelDelays(channelSet, opts)
%CALIBCHANNELDELAYS Absolute on/off delay of each pulse-generator channel, from ONE fit
%to the whole fluorescence trace.
%
%   results = calibChannelDelays({'greenLaser'})
%   results = calibChannelDelays({'greenLaser','MW','trigger'}, opts)
%   results = calibChannelDelays([], struct('data', r))      % refit an earlier result
%
% LEVEL 1 CALIBRATION. Depends on Level 0 (a working detector and a chosen laser power).
% Everything above it depends on this: if the channel delays are wrong, every detection
% window at Level 2 is measured against a shifted time origin.
%
% ===================================================================
% fixDelays MUST BE 0
% ===================================================================
%
% fixDelays makes the pulse generator apply the stored onDelay/offDelay from
% setupInfo.json before uploading. With fixDelays = 0 this measures the ABSOLUTE delay,
% which is the value to store. With 1 it would measure the RESIDUAL left after the stored
% compensation -- and a residual stored as a delay collapses the correction towards zero.
% This function forces fixDelays = 0.
%
% ===================================================================
% THE MEASUREMENT
% ===================================================================
%
% ExpAWGDigitizerCalibrateDelays plays a staircase: the listed channels switch on one at
% a time, channelDelay apart, then off again in reverse order, while the digitizer
% records the whole window. With N channels and offset = detectionStartStopOffset:
%
%   channel k switches ON  at  offset + (k-1)*channelDelay
%   channel k switches OFF at  offset + (2N-k)*channelDelay
%
% The trace is taken through ExpCalibrateDelaysTrace (Setup 7\experiments), which is the
% same experiment with plain averaging of the raw trace. See that file for why: the
% inherited processData's glitch filter is what made every run with channelDelay > 5 us
% come back empty (and crash at the old line 223).
%
% ===================================================================
% THE FIT -- one function for the whole trace
% ===================================================================
%
%   F(t) = dark + light * L(t) * (1 - sum_j C_j * D_j(t))
%
% Every switching moment is a Heaviside step at its own time t0, followed by an
% exponential response, and the delay is where that Heaviside sits:
%
%   S(t; t0, tau) = H(t - t0) * (1 - exp(-(t - t0)/tau))
%
% THE LASER   L(t) = S(t; tOn, tauRise) - S(t; tOff, tauFall)
%             Both steps are smoothed by a Gaussian of width sigma (one sigma for the
%             trace): the AOM and the detector make the edges S-shaped, and a sharp
%             corner forced onto an S-shaped edge drags t0 around. sigma -> 0 is the bare
%             Heaviside. In the full model the light also starts ~1-5% high after turn-on
%             and settles over ~0.1-0.5 us (metastable singlet filling): L is multiplied by
%             1 + K exp(-(t - tOn)/tauSettle).
%
% THE MW      D_j(t) = drop(t; tStart) - recovery(t; tEnd), one per window where
%             microwave reaches the NVs (see WHICH EDGES ARE VISIBLE). C_j is the drop as a
%             fraction of the light -- a few percent, so the level it drops to stays close
%             to the maximum.
%   drop      full model: Rabi nutation,
%                 H(x) * (1 - exp(-x/tau) (cos Wx + sin Wx/(W tau))),   x = t - tStart
%             which starts with ZERO slope (the microwave moves the population
%             quadratically at first), overshoots, rings at W and settles. W -> 0 is the
%             critically damped step. A plain exponential start puts tStart ~20 ns late
%             on a Rabi-shaped edge.
%   recovery  full model: (1-w) S(t; tEnd, tauFast) + w S(t; tEnd, tauSlow). The NV is
%             re-polarised fast (~0.5 us) but part of the light comes back over several
%             us. Fitting this with ONE exponential (or reading a 90% crossing off it, as
%             the old version did) puts the switch hundreds of ns to 1.5 us late.
%
% opts.model = 'simple' drops the extras: smoothed Heaviside-times-exponential everywhere
% and nothing else (7 parameters for the laser, 5 per MW window). 'full' (default) has
% 9 + 8. Tested on synthetic traces with Rabi-shaped drops and two-part recoveries: 'full'
% returns the laser edges to < 1 ns and the MW start to ~2 ns; 'simple' is 20 ns late on
% the MW start and 300 ns early on the MW end.
%
% Weighted by the scatter across averages. Levenberg-Marquardt written here (base MATLAB,
% no toolbox), started from the data: levels from medians, laser edges from half-level
% crossings, MW edges and the Rabi frequency from chi-square scans. Errors are from the
% covariance, inflated by sqrt(chi2/dof) when that is above 1 -- it will be, because the
% Rabi ringing beats (four NV orientations, three hyperfine lines) and one frequency
% cannot follow it after the first dip. That does not move the start time.
%
% ===================================================================
% WHERE THE DELAY IS READ  (opts.delayAt)
% ===================================================================
%
% 0.9 (default)      delay read where the fitted step is 90% complete -- 90% of the level
%                    on a turn-on, 10% left on a turn-off. This is the Setup 7 convention:
%                    the values in setupInfo.json and Level 2's light-on point use it. It
%                    is now read off the FITTED curve instead of interpolated between noisy
%                    samples.
% 'onset'            delay = t0 - commanded, the position of the Heaviside: the moment the
%                    hardware started to switch, whatever the response after it.
%
% Both are always in results.edges (t0 and t90) and both are printed. For the laser they
% differ by ~2.3*tauRise ~ 65 ns by definition alone.
%
% ONE EXCEPTION: an edge whose response is slow (tau > 0.2 us -- the MW recovery, and the
% end of the AWG waveform) is read at its onset under either convention. Its 90% point is
% the NV being re-polarised, ~1.3 us after the switch; the old version reported exactly
% that as an MW off delay. Such edges are marked in the table.
%
% ===================================================================
% WHICH EDGES ARE VISIBLE
% ===================================================================
%
% Microwave reaches the NVs only while BOTH its switch (MW or MW2) is open AND the AWG is
% playing. The AWG starts on trigger ON and plays a fixed-length waveform that LoadAWG
% sizes to end at the commanded trigger OFF; the trigger's own falling edge does nothing.
% So the ORDER of the list decides which edges exist:
%
%   {greenLaser, MW, trigger}  switch already open, the window is the AWG itself.
%                              Seen: trigger ON, and the end of the AWG waveform,
%                              reported as trigger OFF = trigger ON + waveform length --
%                              a second look at the trigger-on delay, never compared with
%                              a stored off delay.
%   {greenLaser, trigger, MW}  AWG already playing, the window is the switch.
%                              Seen: MW ON and MW OFF.
%   MW without trigger, or trigger without an MW switch: no microwave, nothing to see.
%
% Edges that cannot be seen are listed with the reason instead of being fitted to noise.
% MW work needs zero field (all four NV orientations at 2870 MHz). This function does NOT
% touch the Helmholtz -- run h.control('off') first.
%
% ===================================================================
% OPTIONS (all optional, struct fields)
%   channelDelay  us per staircase step (default 1; use >= 3 with the trigger)
%   offset        detectionStartStopOffset in us (default 1)
%   repeats, averages   default 100, 100
%   frequency     MHz for the AWG (default 2870, the zero-field line)
%   amplitude     dBm (default -10)
%   jsonLocation  expProperties.json for loadJsonParams (default: newest results folder)
%   delayAt       0.9 (default) or another completion fraction, or 'onset' -- see above
%   model         'full' (default) or 'simple' -- see THE FIT
%   plot          default true
%   data          an earlier results struct (anything with time, trace, trace_err):
%                 refit it without touching hardware. Its channelSet, channelDelay and
%                 offset are used when it has them.
%
% Writes nothing. Prints the measured delays next to the stored ones and what, if
% anything, to change.

    if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
    opts = withDefaults(opts, struct('channelDelay', 1, 'offset', 1, ...
                                     'repeats', 100, 'averages', 100, ...
                                     'frequency', 2870, 'amplitude', -10, ...
                                     'jsonLocation', '', 'delayAt', 0.9, ...
                                     'model', 'full', 'plot', true, 'data', []));
    if ~any(strcmpi(opts.model, {'full', 'simple'}))
        error('calibChannelDelays: model must be ''full'' or ''simple''.');
    end
    opts.model = lower(opts.model);

    replay = ~isempty(opts.data);
    if replay
        D = opts.data;
        if ~isstruct(D) || ~isfield(D, 'time') || ~isfield(D, 'trace')
            error('calibChannelDelays: opts.data must be a struct with fields time and trace.');
        end
        if (~exist('channelSet', 'var') || isempty(channelSet)) && isfield(D, 'channelSet')
            channelSet = D.channelSet;
        end
        if isfield(D, 'channelDelay'); opts.channelDelay = D.channelDelay; end
        if isfield(D, 'offset');       opts.offset       = D.offset;       end
    end
    if ~exist('channelSet', 'var') || isempty(channelSet)
        channelSet = {'greenLaser'};
    end
    if ~iscell(channelSet); channelSet = {channelSet}; end
    if ~(isnumeric(opts.channelDelay) && isscalar(opts.channelDelay) && opts.channelDelay > 0)
        error('calibChannelDelays: channelDelay must be one positive number of us.');
    end
    delayAt = checkDelayAt(opts.delayAt);

    channelSet = channelSet(:)';
    names = cellfun(@charOf, channelSet, 'UniformOutput', false);
    cd = opts.channelDelay;
    of = opts.offset;

    fprintf('\n=== calibChannelDelays: {%s} ===\n', strjoin(names, ', '));
    fprintf('staircase: %g us per step, %g us offset each side\n', cd, of);
    if ischar(delayAt)
        fprintf('delay = fitted switch time (the position of each step''s Heaviside).\n');
    else
        fprintf('delay = time the fitted step is %.0f%% complete.\n', 100*delayAt);
    end

    % The staircase step has to be longer than the longest delay, or an edge lands in the
    % next step and cannot be attributed. The AWG trigger's delay is ~1.4 us on Setup 7.
    if any(strcmpi(names, 'trigger')) && cd < 3
        fprintf(['\nWARNING channelDelay = %g us but the AWG trigger delay is ~1.4 us. Use\n' ...
                 '        opts.channelDelay = 3 or more for any set containing the trigger.\n'], cd);
    end
    if any(strcmpi(names, 'trigger')) || any(strncmpi(names, 'MW', 2))
        fprintf(['NOTE this set involves microwave. Run it at ZERO FIELD (all four NV\n' ...
                 '     orientations at 2870 MHz). This function does not touch the Helmholtz.\n' ...
                 '     Using %g MHz.\n'], opts.frequency);
    end

    % Which edges can be seen, and which window each MW dip belongs to. Done before the
    % run, so a set with nothing to see fails in a second rather than after the averages.
    M = buildModel(names, cd, of);

    % ---- the trace ----------------------------------------------------------
    if replay
        t = D.time(:);
        y = D.trace(:);
        e = [];
        if isfield(D, 'trace_err'); e = D.trace_err(:); end
        nAvg = NaN;
        if isfield(D, 'nAverages'); nAvg = D.nAverages; end
        fprintf('REPLAY: refitting %d samples. No hardware touched.\n', numel(t));
    else
        % Resolved by setup7Paths. Tolerant on purpose: an empty jsonLocation just
        % leaves the class defaults.
        if isempty(opts.jsonLocation)
            try
                opts.jsonLocation = setup7Paths().json;
            catch
                opts.jsonLocation = '';
            end
        end

        c = ExpCalibrateDelaysTrace;
        if ~isempty(opts.jsonLocation)
            c.loadJsonParams('all', opts.jsonLocation);
        end
        c.listOfChannels           = channelSet;
        c.channelDelay             = cd;
        c.detectionStartStopOffset = of;
        c.repeats                  = opts.repeats;
        c.averages                 = opts.averages;
        c.frequency                = opts.frequency;
        c.amplitude                = opts.amplitude;
        % The json turns the glitch filter on. On a trace it discards whole repeats for
        % one outlying sample out of thousands -- see ExpCalibrateDelaysTrace.
        c.checkGlitchInData        = false;
        % NOT negotiable. See the header.
        c.fixDelays = 0;
        fprintf('fixDelays = 0, so this measures ABSOLUTE delays, not residuals.\n');

        c.run;
        [t, y, e, nAvg] = traceFromExperiment(c);
    end

    results = struct('channelSet', {channelSet}, 'time', t, 'trace', y, 'trace_err', e, ...
                     'channelDelay', cd, 'offset', of, 'fixDelays', 0, ...
                     'nAverages', nAvg, 'delayAt', {opts.delayAt});

    % ---- the fit ------------------------------------------------------------
    [fitR, edges] = fitTrace(t, y, e, M, cd, delayAt, opts.model);
    results.edges = edges;
    results.fit   = fitR;

    printEdges(edges, fitR, delayAt);

    stored = storedDelays();
    compareStored(edges, stored, delayAt);

    if opts.plot
        plotFit(fitR, edges, names, delayAt, cd);
    end
    results.actions = delayVerdict(edges, stored, delayAt);
end

% ======================================================================
% WHAT CAN BE SEEN
% ======================================================================

function M = buildModel(names, cd, of)
    N = numel(names);
    onT  = of + (0:N-1) * cd;
    offT = of + (2*N - (1:N)) * cd;

    % Edge order as before: every ON in list order, then every OFF in reverse.
    edges = repmat(blankEdge('', '', NaN), 1, 0);
    for k = 1:N
        edges(end+1) = blankEdge(names{k}, 'on',  onT(k));   %#ok<AGROW>
    end
    for k = N:-1:1
        edges(end+1) = blankEdge(names{k}, 'off', offT(k));  %#ok<AGROW>
    end
    onIdx  = @(k) k;
    offIdx = @(k) 2*N - k + 1;

    iL = find(strcmpi(names, 'greenLaser'), 1);
    if isempty(iL)
        error(['calibChannelDelays: greenLaser must be in the channel set. Without light ' ...
               'there is no fluorescence and no edge to see.']);
    end
    M.laser = struct('tOnCmd', onT(iL), 'tOffCmd', offT(iL), ...
                     'onEdge', onIdx(iL), 'offEdge', offIdx(iL));
    edges(onIdx(iL)).role  = 'laserOn';
    edges(offIdx(iL)).role = 'laserOff';

    iT = find(strcmpi(names, 'trigger'), 1);
    iS = find(strncmpi(names, 'MW', 2));
    M.dips = struct('switchName', {}, 'tStartCmd', {}, 'tEndCmd', {}, 'startEdge', {}, ...
                    'endEdge', {}, 'endIsAwg', {}, 'wf', {});

    if isempty(iT) && ~isempty(iS)
        for s = iS(:)'
            edges(onIdx(s)).note  = 'no trigger in this set, so the AWG never plays';
            edges(offIdx(s)).note = edges(onIdx(s)).note;
        end
    elseif ~isempty(iT) && isempty(iS)
        edges(onIdx(iT)).note  = 'no MW switch in this set, so no microwave reaches the NVs';
        edges(offIdx(iT)).note = edges(onIdx(iT)).note;
    elseif ~isempty(iT)
        % The AWG waveform length, exactly as ExpAWGDigitizerCalibrateDelays.LoadAWG sets
        % it: from trigger ON to the commanded trigger OFF.
        if iT == N
            wf = cd;
        else
            wf = (2*(N - iT) + 1) * cd;
        end
        awgOn  = onT(iT);
        awgOff = awgOn + wf;
        for s = iS(:)'
            a = max(onT(s), awgOn);
            b = min(offT(s), awgOff);
            if b <= a
                edges(onIdx(s)).note  = 'switch never open while the AWG plays';
                edges(offIdx(s)).note = edges(onIdx(s)).note;
                continue
            end
            if a < onT(iL) || b > offT(iL)
                edges(onIdx(s)).note  = 'MW window is outside the laser window -- dark, so invisible';
                edges(offIdx(s)).note = edges(onIdx(s)).note;
                continue
            end
            if onT(s) > awgOn;  eStart = onIdx(s);  else; eStart = onIdx(iT);  end
            if offT(s) < awgOff
                eEnd = offIdx(s);  endIsAwg = false;
            else
                eEnd = offIdx(iT); endIsAwg = true;
            end
            j = numel(M.dips) + 1;
            M.dips(j) = struct('switchName', names{s}, 'tStartCmd', a, 'tEndCmd', b, ...
                               'startEdge', eStart, 'endEdge', eEnd, 'endIsAwg', endIsAwg, ...
                               'wf', wf);
            edges(eStart).role = 'mwStart';
            edges(eEnd).role   = 'mwEnd';
            edges(eStart).dip  = j;
            edges(eEnd).dip    = j;
            if endIsAwg
                edges(eEnd).compare = false;
                edges(eEnd).note = sprintf(['end of the %g us AWG waveform = trigger ON + ' ...
                                            'waveform; checks the trigger-on delay'], wf);
            end
        end
    end

    % Two switches open over the SAME window (both opened before the trigger, say) give
    % one dip whose two contrasts cannot be told apart. Fit them as one window.
    j = 2;
    while j <= numel(M.dips)
        same = find([M.dips(1:j-1).startEdge] == M.dips(j).startEdge & ...
                    [M.dips(1:j-1).endEdge]   == M.dips(j).endEdge, 1);
        if ~isempty(same)
            M.dips(same).switchName = [M.dips(same).switchName '+' M.dips(j).switchName];
            M.dips(j) = [];
            for i = 1:numel(edges)
                if edges(i).dip == j; edges(i).dip = same; elseif edges(i).dip > j; edges(i).dip = edges(i).dip - 1; end
            end
            % the switches' own edges now carry no window of their own
        else
            j = j + 1;
        end
    end

    % Say why everything else is invisible.
    for i = 1:numel(edges)
        if ~isempty(edges(i).role) || ~isempty(edges(i).note); continue; end
        ch = edges(i).channel;
        if strncmpi(ch, 'MW', 2)
            edges(i).note = 'switch moves while the AWG is silent';
        elseif strcmpi(ch, 'trigger')
            edges(i).note = 'no MW switch is open at this moment';
        elseif strcmpi(ch, 'detector')
            edges(i).note = 'the acquisition gate itself';
        else
            edges(i).note = 'no optical effect expected';
        end
    end
    M.edges = edges;
end

% ======================================================================
% THE TRACE FROM THE EXPERIMENT
% ======================================================================

function [t, y, e, K] = traceFromExperiment(c)
% Straight from c.signal (one trace per average), NOT from summaryData: summaryData is
% whatever the plot drew, and when the trace is all NaN the plot draws a placeholder --
% which is how the old version ended up indexing a short vector with a long mask.
    K = c.currIter;
    if K < 1
        error('calibChannelDelays:noAverage', ['calibChannelDelays: no average completed -- ' ...
              'every trial of the first average failed.\nThe reason is in the "Experiment ' ...
              'failed at trial ..." warnings printed above.']);
    end
    t = c.time(:);
    S = reshape(c.signal(1, :, 1:K), [], K);
    E = reshape(c.sterr(1, :, 1:K), [], K);
    n = min(numel(t), size(S, 1));
    t = t(1:n);  S = S(1:n, :);  E = E(1:n, :);

    ok = isfinite(S);
    k  = sum(ok, 2);
    S0 = S;  S0(~ok) = 0;
    y  = sum(S0, 2) ./ k;
    if K >= 3
        % Scatter of the averages about their mean: includes slow drifts, which the
        % within-average error cannot see.
        d = S - y;  d(~ok) = 0;
        e = sqrt(sum(d.^2, 2) ./ max(k - 1, 1)) ./ sqrt(max(k, 1));
    else
        E0 = E;  E0(~isfinite(E0)) = 0;
        e = sqrt(sum(E0.^2, 2) ./ max(k, 1)) ./ sqrt(max(k, 1));
    end

    bad = ~isfinite(y);
    if all(bad)
        error(['calibChannelDelays: %d averages completed but every sample is NaN.\n' ...
               'The digitizer returned no usable data -- check that the photodiode input ' ...
               'is connected and that the detector channel triggers the digitizer.'], K);
    end
    fprintf('\n%d averages, %d samples of %.0f ns', K, n, 1e3*median(diff(t)));
    if any(bad)
        fprintf(' (%d NaN samples ignored)', nnz(bad));
    end
    fprintf('\n');
end

% ======================================================================
% THE FIT
% ======================================================================

function [R, edges] = fitTrace(t, y, e, M, cd, delayAt, modelKind)
    edges = M.edges;
    t = t(:);  y = y(:);
    if isempty(e); e = NaN(size(y)); end
    e = e(:);
    keep = isfinite(t) & isfinite(y);
    t = t(keep);  y = y(keep);  e = e(keep);
    if numel(t) < 20
        error('calibChannelDelays: only %d usable samples -- nothing to fit.', numel(t));
    end
    dt = median(diff(t));

    % Per-sample error. A floor keeps one lucky sample from taking all the weight.
    eOK = isfinite(e) & e > 0;
    if nnz(eOK) < 0.5*numel(e)
        e(:) = robustNoise(y);
    else
        floorE = 0.25 * median(e(eOK));
        e(~eOK | e < floorE) = floorE;
    end

    % ---- starting values, read off the data -------------------------------
    L = M.laser;
    pre  = t < L.tOnCmd - 0.02;
    post = t > L.tOffCmd + 0.6;
    if nnz(pre) >= 5
        dark0 = median(y(pre));
    elseif nnz(post) >= 5
        dark0 = median(y(post));
    else
        ys = sort(y);  dark0 = ys(max(1, round(0.02*numel(ys))));
    end
    tPlatEnd = L.tOffCmd;
    if ~isempty(M.dips); tPlatEnd = min([M.dips.tStartCmd]); end
    plat = t > L.tOnCmd + 0.3 & t < tPlatEnd;
    if nnz(plat) < 5; plat = t > L.tOnCmd + 0.3 & t < L.tOffCmd; end
    top0 = median(y(plat));
    lev0 = top0 - dark0;
    if ~(lev0 > 5*median(e))
        error(['calibChannelDelays: no laser step in the trace (light minus dark = %.3g, ' ...
               'noise %.3g). Is the laser reaching the sample and the photodiode?'], ...
              lev0, median(e));
    end

    ys = movmean(y, 3);
    half = dark0 + lev0/2;
    tOn0  = crossTime(t, ys, half, L.tOnCmd  - 0.1*cd, +1);
    tOff0 = crossTime(t, ys, half, L.tOffCmd - 0.1*cd, -1);
    if isnan(tOn0);  tOn0  = L.tOnCmd  + 0.2; end
    if isnan(tOff0); tOff0 = L.tOffCmd + 0.2; end

    % ---- parameter vector -------------------------------------------------
    % Times in us; every tau and sigma as a logarithm, so it stays positive. P collects
    % (name, kind, start value); I says where each one sits in the vector.
    full = strcmp(modelKind, 'full');
    P = struct('name', {}, 'kind', {}, 'v', {}, 'lo', {}, 'hi', {});
    [P, I.b]   = addPar(P, 'dark',        'lin',  dark0);
    [P, I.A]   = addPar(P, 'light',       'lin',  lev0);
    [P, I.on]  = addPar(P, 'tOn',         'time', tOn0 - 0.01);
    [P, I.tr]  = addPar(P, 'ln tauRise',  'log',  log(0.03), log(1e-3), log(2));
    [P, I.off] = addPar(P, 'tOff',        'time', tOff0 - 0.01);
    [P, I.tf]  = addPar(P, 'ln tauFall',  'log',  log(0.03), log(1e-3), log(2));
    [P, I.sg]  = addPar(P, 'ln sigma',    'log',  log(0.01), log(1e-4), log(0.2));
    I.K = 0;  I.tK = 0;
    if full
        % After the laser comes on the light starts ~1% high and settles over a few
        % hundred ns (the metastable singlet filling, charge-state equilibration). Its
        % time constant is kept above 0.1 us so it cannot borrow the edge's own shape.
        seg = t > tOn0 + 0.15 & t < tOn0 + 0.4;
        K0 = 0.005;
        if nnz(seg) > 3; K0 = min(max((mean(y(seg)) - top0) / lev0, 0.001), 0.2); end
        [P, I.K]  = addPar(P, 'K settle',      'frac', K0, 0, 0.5);
        [P, I.tK] = addPar(P, 'ln tauSettle',  'log',  log(0.5), log(0.1), log(5));
    end

    nd = numel(M.dips);
    I.seen = false(1, nd);
    [I.C, I.ta, I.tda, I.W, I.tb, I.tdb, I.w, I.tdb2] = deal(zeros(1, nd));
    timePar = zeros(1, numel(edges));     % an edge shared by two MW windows gets ONE time
    scan = zeros(0, 4);                   % [parameter, lo, hi, isEnd] to scan first
    dipInfo = repmat(struct('z', NaN, 'C0', NaN), 1, nd);
    for j = 1:nd
        a = M.dips(j).tStartCmd;  b = M.dips(j).tEndCmd;
        inDip = t > a + 0.6*cd & t < b;
        if nnz(inDip) < 5; inDip = t > a + 0.5*(b - a) & t < b; end
        z = (median(y(plat)) - median(y(inDip))) / ...
            sqrt(var(y(plat))/nnz(plat) + var(y(inDip))/nnz(inDip));
        C0 = (top0 - median(y(inDip))) / lev0;
        dipInfo(j) = struct('z', z, 'C0', C0);
        if ~(abs(z) >= 5)
            msg = sprintf(['no MW dip seen (%.1f sigma) -- AWG on resonance? field off? ' ...
                           'switch cabled?'], abs(z));
            edges(M.dips(j).startEdge).note = msg;
            edges(M.dips(j).endEdge).note   = msg;
            continue
        end
        I.seen(j) = true;
        [P, I.C(j)] = addPar(P, sprintf('C%d', j), 'frac', C0, -1, 1);
        es = M.dips(j).startEdge;
        if timePar(es) == 0
            [P, timePar(es)] = addPar(P, sprintf('tStart%d', j), 'time', a + 0.1);
            scan(end+1, :) = [timePar(es), a - 0.1*cd, a + 0.9*cd, 0];   %#ok<AGROW>
        end
        I.ta(j) = timePar(es);
        [P, I.tda(j)] = addPar(P, sprintf('ln tauDrop%d', j), 'log', log(0.03), log(1e-3), log(5));
        if full
            % The drop is Rabi nutation: it starts with zero slope (the microwave moves the
            % population quadratically at first), overshoots, rings at the Rabi frequency
            % and settles. W in rad/us; W = 0 is the critically damped limit.
            [P, I.W(j)] = addPar(P, sprintf('rabi%d', j), 'rabi', 0, 0, 2*pi*40);
        end
        ee = M.dips(j).endEdge;
        if timePar(ee) == 0
            [P, timePar(ee)] = addPar(P, sprintf('tEnd%d', j), 'time', b + 0.1);
            scan(end+1, :) = [timePar(ee), b - 0.1*cd, b + 0.9*cd, 1];   %#ok<AGROW>
        end
        I.tb(j) = timePar(ee);
        [P, I.tdb(j)] = addPar(P, sprintf('ln tauRecovery%d', j), 'log', log(0.5), log(1e-3), log(20));
        if full
            % The recovery is not one exponential: a fast re-polarisation and a slow tail
            % that has not finished microseconds later.
            [P, I.w(j)]    = addPar(P, sprintf('w slow%d', j),           'frac', 0.4, 0, 1);
            [P, I.tdb2(j)] = addPar(P, sprintf('ln tauRecSlow%d', j),    'log',  log(3), log(0.2), log(50));
        end
    end
    p0 = [P.v]';
    pnames = {P.name};
    kinds  = {P.kind};

    % ---- MW edges: scan first ----------------------------------------------
    % A dip edge is a few percent of the light; started blind, the optimiser can settle
    % on the wrong side of it. A chi-square scan along the edge's own step, everything
    % else held, finds it; the fit then refines everything at once.
    p = p0;
    for s = 1:size(scan, 1)
        k = scan(s, 1);  lo = scan(s, 2);  hi = scan(s, 3);
        % A start is judged on its own step only (the end of the window may still be
        % misplaced); an end also needs the slow recovery after it.
        if scan(s, 4)
            win = t >= lo & t <= min(hi + 1, L.tOffCmd);
        else
            win = t >= lo & t <= hi;
        end
        grid = lo : dt : hi;
        chi = zeros(size(grid));
        for g = 1:numel(grid)
            q = p;  q(k) = grid(g);
            r = (y(win) - modelEval(q, t(win), I)) ./ e(win);
            chi(g) = r' * r;
        end
        [~, m] = min(chi);
        p(k) = grid(m);
    end

    % Rabi frequency and damping of each MW start: a coarse grid, then the start time
    % again with them in place. The chi-square is flat in W at W = 0 (the response is
    % even in W), so the optimiser alone could never leave it.
    for j = find(I.seen & I.W > 0)
        ta = p(I.ta(j));
        win = t >= ta - 0.1 & t <= ta + min(1.5, 0.9*cd);
        best = Inf;
        for Wg = 2*pi*(0:0.5:20)
            for tg = [0.02 0.04 0.08 0.15 0.3 0.6]
                q = p;  q(I.W(j)) = Wg;  q(I.tda(j)) = log(tg);
                r = (y(win) - modelEval(q, t(win), I)) ./ e(win);
                if r' * r < best
                    best = r' * r;  bW = Wg;  bT = tg;
                end
            end
        end
        p(I.W(j)) = bW;  p(I.tda(j)) = log(bT);
        grid = ta - 0.1 : dt/2 : ta + 0.1;
        chi = zeros(size(grid));
        for g = 1:numel(grid)
            q = p;  q(I.ta(j)) = grid(g);
            r = (y(win) - modelEval(q, t(win), I)) ./ e(win);
            chi(g) = r' * r;
        end
        [~, m] = min(chi);
        p(I.ta(j)) = grid(m);
    end

    % ---- the fit ------------------------------------------------------------
    h = zeros(size(p));
    for k = 1:numel(p)
        switch kinds{k}
            case 'lin';  h(k) = 1e-4 * lev0;
            case 'frac'; h(k) = 1e-5;
            case 'time'; h(k) = 2e-4;            % 0.2 ns
            case 'log';  h(k) = 1e-3;
            case 'rabi'; h(k) = 0.02;            % rad/us
        end
    end
    pLo = [P.lo]';  pHi = [P.hi]';
    clampFn = @(q) min(max(q, pLo), pHi);

    resid = @(q) (y - modelEval(q, t, I)) ./ e;
    [p, J, chi2, nIter, converged] = lmFit(resid, clampFn(p), h, clampFn);

    nP  = numel(p);
    dof = max(numel(y) - nP, 1);
    chi2red = chi2 / dof;
    Cov = pinv(J' * J) * max(1, chi2red);
    perr = sqrt(abs(diag(Cov)));
    curve = modelEval(p, t, I);

    % ---- edges --------------------------------------------------------------
    sg = exp(p(I.sg));
    lev = p(I.A);
    edges = setEdge(edges, L.onEdge,  I.on,  I.tr, +lev, p, Cov, I.sg, delayAt);
    edges = setEdge(edges, L.offEdge, I.off, I.tf, -lev, p, Cov, I.sg, delayAt);
    for j = 1:nd
        if ~I.seen(j); continue; end
        dA = p(I.C(j)) * lev;
        es = M.dips(j).startEdge;  ee = M.dips(j).endEdge;
        if isnan(edges(es).t0)      % a start shared by two windows keeps the first one's tau
            edges = setEdge(edges, es, I.ta(j), I.tda(j), -dA, p, Cov, I.sg, delayAt, I.W(j));
        end
        if isnan(edges(ee).t0)
            edges = setEdge(edges, ee, I.tb(j), I.tdb(j), +dA, p, Cov, I.sg, delayAt);
        end
    end

    % Sanity, on every fitted edge. These are the ways a fit can look like a number and
    % not be one.
    for i = 1:numel(edges)
        if ~edges(i).visible; continue; end
        del = edges(i).tSwitch - edges(i).commanded;
        if ~isfinite(edges(i).tSwitchSterr) || edges(i).tSwitchSterr > 0.5*cd
            edges(i).visible = false;
            edges(i).note = sprintf('REJECTED: not determined (+/- %.3g us)', edges(i).tSwitchSterr);
        elseif del < -0.02*cd
            edges(i).visible = false;
            edges(i).note = 'REJECTED: negative delay -- the fit caught the wrong edge';
        elseif del > 0.9*cd
            edges(i).visible = false;
            edges(i).note = sprintf('REJECTED: delay %.2f us runs into the next step', del);
        end
    end

    % ---- what the fit returns --------------------------------------------------
    R = struct();
    R.model = modelKind;
    R.names = pnames;  R.p = p;  R.perr = perr;  R.cov = Cov;
    R.chi2 = chi2;  R.dof = dof;  R.chi2red = chi2red;
    R.converged = converged;  R.iterations = nIter;
    R.t = t;  R.y = y;  R.e = e;  R.curve = curve;  R.residual = y - curve;
    R.dark = p(I.b);   R.darkErr = perr(I.b);
    R.level = lev;     R.levelErr = perr(I.A);
    R.sigma = sg;      R.sigmaErr = sg * perr(I.sg);
    R.tauRise = exp(p(I.tr));  R.tauRiseErr = R.tauRise * perr(I.tr);
    R.tauFall = exp(p(I.tf));  R.tauFallErr = R.tauFall * perr(I.tf);
    R.settle = NaN;  R.tauSettle = NaN;
    if I.K > 0
        R.settle = p(I.K);  R.tauSettle = exp(p(I.tK));
    end
    R.dips = struct('switchName', {}, 'seen', {}, 'significance', {}, 'contrast', {}, ...
                    'contrastErr', {}, 'tauDrop', {}, 'tauDropErr', {}, 'rabiMHz', {}, ...
                    'tauRecovery', {}, 'tauRecoveryErr', {}, ...
                    'slowFraction', {}, 'tauRecoverySlow', {}, 'endIsAwg', {});
    for j = 1:nd
        d = struct('switchName', M.dips(j).switchName, 'seen', I.seen(j), ...
                   'significance', dipInfo(j).z, 'contrast', NaN, 'contrastErr', NaN, ...
                   'tauDrop', NaN, 'tauDropErr', NaN, 'rabiMHz', NaN, ...
                   'tauRecovery', NaN, 'tauRecoveryErr', NaN, 'slowFraction', NaN, ...
                   'tauRecoverySlow', NaN, 'endIsAwg', M.dips(j).endIsAwg);
        if I.seen(j)
            d.contrast = p(I.C(j));             d.contrastErr = perr(I.C(j));
            d.tauDrop = exp(p(I.tda(j)));       d.tauDropErr = d.tauDrop * perr(I.tda(j));
            d.tauRecovery = exp(p(I.tdb(j)));   d.tauRecoveryErr = d.tauRecovery * perr(I.tdb(j));
            if I.W(j) > 0; d.rabiMHz = p(I.W(j)) / (2*pi); end
            if I.w(j) > 0
                d.slowFraction = p(I.w(j));     d.tauRecoverySlow = exp(p(I.tdb2(j)));
            end
        end
        R.dips(j) = d;
    end
    R.index = I;
end

function [P, k] = addPar(P, name, kind, v, lo, hi)
    if nargin < 5; lo = -Inf; hi = Inf; end
    P(end+1) = struct('name', name, 'kind', kind, 'v', min(max(v, lo), hi), 'lo', lo, 'hi', hi);
    k = numel(P);
end

function F = modelEval(p, t, I)
% The whole trace, see the header:
%   F = dark + light * L(t) * (1 - sum_j C_j * D_j(t))
    sg = exp(p(I.sg));
    Sr = stepResp(t, p(I.on), exp(p(I.tr)), sg);
    Lt = Sr - stepResp(t, p(I.off), exp(p(I.tf)), sg);
    if I.K > 0
        % the light starts K high and settles with tauSettle
        Lt = Lt .* (1 + p(I.K) * exp(-max(t - p(I.on), 0) / exp(p(I.tK))));
    end
    D = zeros(size(t));
    for j = 1:numel(I.seen)
        if ~I.seen(j); continue; end
        if I.W(j) > 0
            % full model: Rabi nutation, see fitTrace
            Dj = rabiResp(t, p(I.ta(j)), exp(p(I.tda(j))), p(I.W(j)));
        else
            Dj = stepResp(t, p(I.ta(j)), exp(p(I.tda(j))), sg);
        end
        Sb = stepResp(t, p(I.tb(j)), exp(p(I.tdb(j))), sg);
        if I.w(j) > 0
            % recovery: fast part and slow part, both starting at tEnd
            Sb = (1 - p(I.w(j))) * Sb + p(I.w(j)) * stepResp(t, p(I.tb(j)), exp(p(I.tdb2(j))), sg);
        end
        D = D + p(I.C(j)) * (Dj - Sb);
    end
    F = p(I.b) + p(I.A) * Lt .* (1 - D);
end

function s = stepResp(t, t0, tau, sg)
% H(t-t0)*(1 - exp(-(t-t0)/tau)), convolved with a Gaussian of width sg:
%
%   S = Phi(x/sg) - exp(-x/tau + sg^2/(2 tau^2)) * Phi(x/sg - sg/tau),   x = t - t0
%
% The second term is evaluated in whichever form cannot overflow: directly when
% x/sg - sg/tau >= 0 (the exponent is then <= 0), otherwise through erfcx, which turns
% it into exp(-x^2/(2 sg^2)) * erfcx(.)/2 with erfcx bounded.
    x = t - t0;
    if sg < 1e-6 * tau
        s = (x > 0) .* (1 - exp(-max(x, 0) / tau));
        return
    end
    z = x / sg - sg / tau;
    E = zeros(size(x));
    pos = z >= 0;
    E(pos)  = exp(-x(pos)/tau + sg^2/(2*tau^2)) .* 0.5 .* erfc(-z(pos)/sqrt(2));
    E(~pos) = 0.5 * exp(-x(~pos).^2 / (2*sg^2)) .* erfcx(-z(~pos)/sqrt(2));
    s = 0.5 * erfc(-x / (sg*sqrt(2))) - E;
end

function s = rabiResp(t, t0, tau, W)
% Step of a driven two-level system: H(t-t0) * (1 - exp(-x/tau) (cos Wx + sin Wx/(W tau))).
% Zero slope at the start, overshoot and ringing at W, settled after tau. W -> 0 is the
% critically damped step 1 - (1 + x/tau) exp(-x/tau).
    x = max(t - t0, 0);
    if abs(W) * tau < 1e-6
        s = 1 - (1 + x/tau) .* exp(-x/tau);
    else
        s = 1 - exp(-x/tau) .* (cos(W*x) + sin(W*x) / (W*tau));
    end
end

function q = completionOffset(frac, tau, sg, W)
% Time after t0 at which the edge's own step first reaches 'frac'. W = NaN: the smoothed
% exponential step (monotonic); otherwise the Rabi step, which can ring, so the first
% crossing is bracketed on a grid before bisecting.
    if nargin < 4; W = NaN; end
    if isnan(W)
        f = @(x) stepResp(x, 0, tau, sg);
        lo = -10*sg - 1e-9;
        hi = 10*sg + 60*tau;
    else
        f = @(x) rabiResp(x, 0, tau, W);
        per = tau;  if W > 0; per = min(tau, 2*pi/W); end
        xs = 0 : per/200 : 60*tau;
        k = find(f(xs) >= frac, 1);
        if isempty(k); q = NaN; return; end
        lo = xs(max(k-1, 1));  hi = xs(k);
    end
    for it = 1:80
        mid = 0.5*(lo + hi);
        v = f(mid);
        if v < frac
            lo = mid;
        else
            hi = mid;
        end
    end
    q = 0.5*(lo + hi);
end

function edges = setEdge(edges, i, kt, ktau, amp, p, Cov, ksg, delayAt, kW)
% kW: index of the Rabi frequency for a Rabi-shaped MW start, 0 for an exponential step.
    if nargin < 10; kW = 0; end
    t0  = p(kt);
    tau = exp(p(ktau));
    sg  = exp(p(ksg));
    if kW > 0
        W = p(kW);  idx = [kt, ktau, ksg, kW];
    else
        W = NaN;    idx = [kt, ktau, ksg];
    end
    Cs  = Cov(idx, idx);

    edges(i).t0        = t0;
    edges(i).t0Sterr   = sqrt(max(Cov(kt, kt), 0));
    edges(i).tau       = tau;
    edges(i).tauSterr  = tau * sqrt(max(Cov(ktau, ktau), 0));
    edges(i).sigma     = sg;
    edges(i).amplitude = amp;

    [edges(i).t90, edges(i).t90Sterr] = completionTime(0.9, t0, tau, sg, Cs, W);
    if ~ischar(delayAt) && tau > 0.2 && isnan(W)
        % A response this slow is the NV (re-polarisation), not the switch. Its 90% point
        % is a property of the NV, so even under a completion convention the switch is
        % read at the onset.
        delayAt = 'onset';
        edges(i).note = 'slow NV response: read at its onset, not its completion point';
    end
    if ischar(delayAt)
        edges(i).tSwitch      = edges(i).t0;
        edges(i).tSwitchSterr = edges(i).t0Sterr;
    else
        [edges(i).tSwitch, edges(i).tSwitchSterr] = completionTime(delayAt, t0, tau, sg, Cs, W);
    end
    edges(i).visible = true;
end

function [tc, err] = completionTime(frac, t0, tau, sg, Cs, W)
% tc = t0 + q(frac; tau, sigma[, W]), with its error from the covariance of
% (t0, ln tau, ln sigma[, W]).
    q  = completionOffset(frac, tau, sg, W);
    tc = t0 + q;
    d  = 1e-3;
    dq_dlntau = (completionOffset(frac, tau*exp(d), sg, W) - ...
                 completionOffset(frac, tau*exp(-d), sg, W)) / (2*d);
    dq_dlnsg  = (completionOffset(frac, tau, sg*exp(d), W) - ...
                 completionOffset(frac, tau, sg*exp(-d), W)) / (2*d);
    g = [1; dq_dlntau; dq_dlnsg];
    if ~isnan(W)
        dW = max(1e-3*W, 1e-3);
        g(4) = (completionOffset(frac, tau, sg, W + dW) - ...
                completionOffset(frac, tau, sg, max(W - dW, 0))) / (W + dW - max(W - dW, 0));
    end
    err = sqrt(max(g' * Cs * g, 0));
end

function [p, J, chi2, it, converged] = lmFit(resid, p, h, clampFn)
% Levenberg-Marquardt on weighted residuals. Central-difference Jacobian.
% A nearly singular step matrix is expected while the damping is being found, and is
% handled by raising the damping, so its warning is silenced here and restored on exit.
    ids = {'MATLAB:singularMatrix', 'MATLAB:nearlySingularMatrix'};
    if exist('OCTAVE_VERSION', 'builtin')
        ids = {'Octave:singular-matrix', 'Octave:nearly-singular-matrix'};
    end
    old = {};
    for k = 1:numel(ids)
        try
            old{end+1} = warning('query', ids{k}); %#ok<AGROW>
            warning('off', ids{k});
        catch
            % an identifier this version does not know: nothing to silence
        end
    end
    cleanup = onCleanup(@() restoreWarnings(old)); %#ok<NASGU>

    r = resid(p);
    chi2 = r' * r;
    lambda = 1e-3;
    converged = false;
    small = 0;
    J = jacobianOf(resid, p, h);
    for it = 1:300
        A = J' * J;
        g = J' * r;
        dA = diag(A);
        dA(dA <= 0) = max(dA) * 1e-12 + eps;
        stepped = false;
        while lambda <= 1e12
            dp = -(A + lambda * diag(dA)) \ g;
            pn = clampFn(p + dp);
            rn = resid(pn);
            chin = rn' * rn;
            if isfinite(chin) && chin < chi2
                stepped = true;
                break
            end
            lambda = lambda * 10;
        end
        if ~stepped
            converged = true;       % no downhill step at any damping: at the minimum
            break
        end
        rel = (chi2 - chin) / chi2;
        p = pn;  r = rn;  chi2 = chin;
        lambda = max(lambda / 10, 1e-12);
        J = jacobianOf(resid, p, h);
        if rel < 1e-7
            small = small + 1;
            if small >= 3; converged = true; break; end
        else
            small = 0;
        end
    end
end

function restoreWarnings(old)
    for k = 1:numel(old)
        warning(old{k}.state, old{k}.identifier);
    end
end

function J = jacobianOf(resid, p, h)
    r0 = resid(p);
    J = zeros(numel(r0), numel(p));
    for k = 1:numel(p)
        dp = zeros(size(p));  dp(k) = h(k);
        J(:, k) = (resid(p + dp) - resid(p - dp)) / (2*h(k));
    end
end


function tc = crossTime(t, y, level, after, direction)
% First time at or after 'after' where y crosses 'level' in 'direction' (+1 up, -1 down),
% linearly interpolated.
    tc = NaN;
    idx = find(t >= after);
    if numel(idx) < 2; return; end
    v = direction * (y - level);
    for k = idx(1)+1 : numel(t)
        if v(k-1) < 0 && v(k) >= 0
            frac = -v(k-1) / (v(k) - v(k-1));
            tc = t(k-1) + frac * (t(k) - t(k-1));
            return
        end
    end
end

function s = robustNoise(y)
    d = diff(y(:));
    s = 1.4826 * median(abs(d - median(d))) / sqrt(2);
    if ~(s > 0); s = std(y); end
end

% ======================================================================
% REPORTING
% ======================================================================

function printEdges(edges, R, delayAt)
% Columns: where the delay is read, its error, and the OTHER convention for reference.
    if ischar(delayAt)
        col = 'switch t0';  other = '90% point';
    else
        col = sprintf('%.0f%% point', 100*delayAt);  other = 'switch t0';
    end
    fprintf('\n%-11s %-4s %9s %10s %8s %10s %10s %9s   %s\n', 'channel', 'edge', 'commanded', ...
            col, '+/-', other, 'tau', 'DELAY', 'note');
    for i = 1:numel(edges)
        E = edges(i);
        if E.visible
            if ischar(delayAt); oth = E.t90; else; oth = E.t0; end
            fprintf('%-11s %-4s %9.4f %10.4f %8.4f %10.4f %10s %+9.4f   %s\n', E.channel, ...
                    E.kind, E.commanded, E.tSwitch, E.tSwitchSterr, oth, tauStr(E.tau), ...
                    E.tSwitch - E.commanded, E.note);
        else
            fprintf('%-11s %-4s %9.4f %10s %8s %10s %10s %9s   not visible: %s\n', E.channel, ...
                    E.kind, E.commanded, '-', '-', '-', '-', '-', E.note);
        end
    end

    fprintf('\nfit: chi2/dof = %.2f  (%d samples, %d parameters, %d iterations%s)\n', ...
            R.chi2red, numel(R.y), numel(R.p), R.iterations, ...
            ternary(R.converged, '', ', NOT CONVERGED'));
    fprintf('   dark %.4g, light %.4g (+/- %.1g)\n', R.dark, R.level, R.levelErr);
    fprintf('   laser: rise tau %s +/- %s, fall tau %s +/- %s, edge width sigma %s +/- %s\n', ...
            tauStr(R.tauRise), tauStr(R.tauRiseErr), tauStr(R.tauFall), tauStr(R.tauFallErr), ...
            tauStr(R.sigma), tauStr(R.sigmaErr));
    if isfinite(R.settle)
        fprintf('          after turn-on the light starts %.2f%% high and settles with tau %s\n', ...
                100*R.settle, tauStr(R.tauSettle));
    end
    for j = 1:numel(R.dips)
        d = R.dips(j);
        if ~d.seen
            fprintf('   MW window %d (%s): no dip seen (%.1f sigma)\n', j, d.switchName, abs(d.significance));
            continue
        end
        if isfinite(d.rabiMHz)
            fprintf(['   MW window %d (%s): drop %.2f +/- %.2f %% of the light; Rabi ringing at ' ...
                     '%.1f MHz, damped with tau %s\n'], j, d.switchName, 100*d.contrast, ...
                    100*d.contrastErr, d.rabiMHz, tauStr(d.tauDrop));
        else
            fprintf(['   MW window %d (%s): drop %.2f +/- %.2f %% of the light, falls with tau %s, ' ...
                     'recovers with tau %s\n'], j, d.switchName, 100*d.contrast, ...
                    100*d.contrastErr, tauStr(d.tauDrop), tauStr(d.tauRecovery));
        end
        if isfinite(d.slowFraction)
            fprintf('        recovery: %.0f%% fast (tau %s) + %.0f%% slow (tau %s)\n', ...
                    100*(1 - d.slowFraction), tauStr(d.tauRecovery), 100*d.slowFraction, ...
                    tauStr(d.tauRecoverySlow));
        end
    end
    if R.chi2red > 3
        fprintf(['   chi2/dof %.1f: the trace has structure the model does not (e.g. a Rabi ' ...
                 'transient at the MW start, a spike at a laser edge).\n   Look at the residual ' ...
                 'panel; the errors above are already inflated by sqrt(chi2/dof).\n'], R.chi2red);
    end
end

function compareStored(edges, stored, delayAt)
    fprintf('\n--- measured vs stored (setupInfo.json) ---\n');
    if ischar(delayAt)
        fprintf(['   NOTE the stored values were taken at the 90%% point of each edge; these ' ...
                 'are switch times (onset).\n   The laser differs by ~2.3*tau by definition. ' ...
                 'For a like-for-like check run with o.delayAt = 0.9.\n']);
    end
    for i = 1:numel(edges)
        E = edges(i);
        if ~E.visible; continue; end
        meas = E.tSwitch - E.commanded;
        if ~E.compare
            fprintf('   %-11s %-4s measured %+.4f +/- %.4f us   (%s)\n', E.channel, E.kind, ...
                    meas, E.tSwitchSterr, E.note);
            continue
        end
        key = matlab.lang.makeValidName(E.channel);
        if isfield(stored, key) && isfield(stored.(key), E.kind)
            sv = stored.(key).(E.kind);
            if ischar(delayAt)
                otherTxt = sprintf('90%% pt %+.4f', E.t90 - E.commanded);
            else
                otherTxt = sprintf('switch t0 %+.4f', E.t0 - E.commanded);
            end
            fprintf('   %-11s %-4s measured %+.4f +/- %.4f us (%s),  stored %+.4f', ...
                    E.channel, E.kind, meas, E.tSwitchSterr, otherTxt, sv);
            if sv == 0 && abs(meas) > 3*E.tSwitchSterr
                fprintf('   <-- stored is ZERO but a real delay is measured\n');
            elseif abs(meas - sv) > 3*E.tSwitchSterr
                fprintf('   <-- differ by %+.4f us (%.0f sigma)\n', meas - sv, ...
                        abs(meas - sv)/max(E.tSwitchSterr, eps));
            else
                fprintf('   agrees\n');
            end
        else
            fprintf('   %-11s %-4s measured %+.4f +/- %.4f us,  nothing stored\n', ...
                    E.channel, E.kind, meas, E.tSwitchSterr);
        end
    end
end

function plotFit(R, edges, names, delayAt, cd)
% Top: the whole trace with the fit over it. Under it: residuals in units of the error.
% Bottom row: a zoom on every fitted edge, each on its own vertical scale -- an MW dip is
% a few percent of the light and would be a flat line on the trace's scale.
    t = R.t;  y = R.y;  F = R.curve;  e = R.e;
    COL.data      = [0.15 0.42 0.69];      % blue
    COL.fit       = [0.85 0.33 0.10];      % orange -- distinguishable from the blue
    COL.commanded = [0.45 0.45 0.45];      % under deuteranopia, which red is not
    COL.switch    = [0.10 0.10 0.10];
    COL.band      = [0.98 0.88 0.80];

    vis = edges([edges.visible]);
    nz  = numel(vis);
    figure('Color', 'w', 'Name', 'calibChannelDelays');
    if nz > 0
        axMain = axes('Position', [0.06 0.62 0.91 0.29]);
        axRes  = axes('Position', [0.06 0.49 0.91 0.08]);
    else
        axMain = axes('Position', [0.06 0.30 0.91 0.60]);
        axRes  = axes('Position', [0.06 0.10 0.91 0.12]);
    end

    yl = [min([y; F]) max([y; F])];
    if ~(diff(yl) > 0); yl = yl(1) + [-1 1]*max(abs(yl(1)), 1)*0.1; end
    yl = yl + 0.08*diff(yl)*[-1 1];

    hold(axMain, 'on');
    for i = 1:nz
        x1 = vis(i).commanded;  x2 = vis(i).tSwitch;
        patch(axMain, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], COL.band, ...
              'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    hData = plot(axMain, t, y, '-', 'Color', COL.data, 'LineWidth', 1);
    hFit  = plot(axMain, t, F, '-', 'Color', COL.fit, 'LineWidth', 1.5);
    hCmd = [];  hSw = [];
    for i = 1:numel(edges)
        hc = plot(axMain, edges(i).commanded*[1 1], yl, ':', 'Color', COL.commanded, 'LineWidth', 1);
        if isempty(hCmd); hCmd = hc; else; set(hc, 'HandleVisibility', 'off'); end
    end
    for i = 1:nz
        hs = plot(axMain, vis(i).tSwitch*[1 1], yl, '-', 'Color', COL.switch, 'LineWidth', 1);
        if isempty(hSw); hSw = hs; else; set(hs, 'HandleVisibility', 'off'); end
        nearRight = vis(i).tSwitch > t(1) + 0.85*(t(end) - t(1));
        if nearRight
            txt = sprintf('%+.4f us ', vis(i).tSwitch - vis(i).commanded);  al = 'right';
        else
            txt = sprintf(' %+.4f us', vis(i).tSwitch - vis(i).commanded);  al = 'left';
        end
        text(axMain, vis(i).tSwitch, yl(2), txt, 'Color', COL.switch, 'FontWeight', 'bold', ...
             'VerticalAlignment', 'top', 'HorizontalAlignment', al, 'FontSize', 9);
    end
    xlim(axMain, [t(1) t(end)]);  ylim(axMain, yl);
    set(axMain, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', []);
    grid(axMain, 'on');
    ylabel(axMain, 'signal');
    if ischar(delayAt); swl = 'switch (fit t0)'; else; swl = sprintf('%.0f%% point', 100*delayAt); end
    hs = {hData, hFit, hCmd, hSw};
    ls = {'data', 'fit', 'commanded', swl};
    keep = cellfun(@(x) ~isempty(x), hs);
    legend(axMain, [hs{keep}], ls(keep), 'Location', 'best', 'Box', 'off');
    title(axMain, sprintf('channel delays  --  %s   (\\chi^2/dof %.2f)', strjoin(names, ', '), ...
                          R.chi2red), 'FontWeight', 'bold');
    try
        subtitle(axMain, summaryLine(vis));
    catch
        % subtitle is R2020b+; the numbers are on the edges and in the zooms anyway.
    end

    hold(axRes, 'on');
    plot(axRes, t, (y - F)./e, '.', 'Color', COL.data, 'MarkerSize', 4);
    plot(axRes, [t(1) t(end)], [0 0], '-', 'Color', COL.fit);
    xlim(axRes, [t(1) t(end)]);
    rr = (y - F)./e;  rl = max(4, min(30, 1.1*max(abs(rr))));
    ylim(axRes, [-rl rl]);
    set(axRes, 'Box', 'off', 'TickDir', 'out');
    grid(axRes, 'on');
    ylabel(axRes, 'resid/err');
    xlabel(axRes, 'time in detection window (us)');

    w = 0.91 / max(nz, 1);
    for k = 1:nz
        ax = axes('Position', [0.06 + (k-1)*w + 0.02, 0.07, w - 0.06, 0.30]);
        hold(ax, 'on');
        E = vis(k);
        span = min(max(6*E.tau, 0.15), 0.8*cd);
        lo = min(E.commanded, E.t0) - max([0.06, 5*E.sigma, 0.25*span]);
        hi = max(E.t0, E.tSwitch) + span;
        m = t >= lo & t <= hi;
        if nnz(m) < 2; m = abs(t - E.t0) <= 0.2; end
        ylz = [min([y(m); F(m)]) max([y(m); F(m)])];
        if ~(diff(ylz) > 0); ylz = ylz(1) + [-1 1]*1e-6; end
        ylz = ylz + 0.1*diff(ylz)*[-1 1];
        patch(ax, [E.commanded E.tSwitch E.tSwitch E.commanded], ylz([1 1 2 2]), COL.band, ...
              'EdgeColor', 'none');
        plot(ax, t(m), y(m), '.', 'Color', COL.data, 'MarkerSize', 7);
        plot(ax, t(m), F(m), '-', 'Color', COL.fit, 'LineWidth', 1.5);
        plot(ax, E.commanded*[1 1], ylz, ':', 'Color', COL.commanded, 'LineWidth', 1.2);
        plot(ax, E.tSwitch*[1 1], ylz, '-', 'Color', COL.switch, 'LineWidth', 1);
        if ~ischar(delayAt)
            plot(ax, E.t0*[1 1], ylz, '--', 'Color', COL.switch, 'LineWidth', 0.8);
        end
        xlim(ax, [lo hi]);  ylim(ax, ylz);
        set(ax, 'Box', 'off', 'TickDir', 'out');
        grid(ax, 'on');
        xlabel(ax, 'us');
        if k == 1; ylabel(ax, 'signal'); end
        title(ax, sprintf('%s %s\n%+.4f \\pm %.4f us   \\tau %s', E.channel, E.kind, ...
                          E.tSwitch - E.commanded, E.tSwitchSterr, tauStr(E.tau)), ...
              'FontWeight', 'normal', 'FontSize', 9);
    end
end

function s = summaryLine(vis)
    if isempty(vis); s = 'no edge detected'; return; end
    parts = cell(1, numel(vis));
    for i = 1:numel(vis)
        parts{i} = sprintf('%s %s %+.4f\\pm%.4f us', vis(i).channel, vis(i).kind, ...
                           vis(i).tSwitch - vis(i).commanded, vis(i).tSwitchSterr);
    end
    s = strjoin(parts, '   |   ');
end

function actions = delayVerdict(edges, stored, delayAt)
% What to change, where, and what to leave alone.
    fprintf('\n============================================================\n');
    fprintf('WHAT TO DO WITH THIS\n');
    fprintf('============================================================\n');
    actions = struct('channel', {}, 'kind', {}, 'path', {}, 'from', {}, 'to', {}, 'readAtOnset', {});

    nCompared = 0;
    for i = 1:numel(edges)
        E = edges(i);
        if ~E.visible || ~E.compare; continue; end
        key = matlab.lang.makeValidName(E.channel);
        if ~isfield(stored, key) || ~isfield(stored.(key), E.kind); continue; end
        nCompared = nCompared + 1;
        st = stored.(key);
        sv = st.(E.kind);
        meas = E.tSwitch - E.commanded;
        pth = '(unknown key)';
        if isfield(st, [E.kind 'Path']); pth = st.([E.kind 'Path']); end
        if abs(meas - sv) > 3*max(E.tSwitchSterr, eps)
            atOnset = ischar(delayAt) || ~isempty(strfind(E.note, 'slow NV response')); %#ok<STREMP>
            actions(end+1) = struct('channel', E.channel, 'kind', E.kind, ...
                                    'path', pth, 'from', sv, 'to', meas, ...
                                    'readAtOnset', atOnset);  %#ok<AGROW>
        end
    end

    if nCompared == 0
        fprintf(['\n1. NOTHING TO COMPARE. No stored delay was found for the edges measured ' ...
                 '(is JsonInfoReader on the\n   path?), so there is no verdict -- only the ' ...
                 'numbers above.\n']);
    elseif isempty(actions)
        fprintf(['\n1. CHANGE NOTHING. Every measured delay agrees with setupInfo.json within ' ...
                 'the error.\n   On a setup that has not been re-cabled that is the expected ' ...
                 'outcome -- this run was a CHECK, and it passed.\n']);
    else
        fprintf('\n1. These differ from setupInfo.json by more than 3 sigma:\n\n');
        for i = 1:numel(actions)
            a = actions(i);
            tag = '';
            if a.readAtOnset && ~ischar(delayAt); tag = '   (read at onset: slow NV recovery)'; end
            fprintf('       %-40s  %.4f  ->  %.4f%s\n', a.path, a.from, a.to, tag);
        end
        if ischar(delayAt)
            fprintf(['\n   CAREFUL: these are switch times (onset); the stored values are 90%% ' ...
                     'points. A difference\n   of ~2.3*tau on the laser is the change of ' ...
                     'definition, not a drift. Either rerun with\n   o.delayAt = 0.9 for a ' ...
                     'like-for-like check, or adopt the onset convention for EVERY channel at\n' ...
                     '   once and redo Level 2 afterwards (it reads the light at its 90%% point).\n']);
        end
        fprintf(['\n   The file is Control code\\Lab json\\Setup 7\\setupInfo.json, in MICROSECONDS. ' ...
                 'Restart main after\n   editing -- the json is read once at construction.\n']);
    end

    fprintf(['\n2. THESE ARE ABSOLUTE DELAYS. fixDelays was 0, so the number is the raw hardware ' ...
             'delay. Write it\n   in as it is; do not add it to the stored one.\n']);

    fprintf(['\n3. ON DELAYS MATTER MOST: the channel has to be on when the detection window ' ...
             'opens. An OFF\n   delay only moves dead time. The MW-off switch time is now the ' ...
             'START of the recovery, so it\n   is a real switch time, but the recovery itself ' ...
             '(tau ~ 0.5-1 us) is the NV being re-polarised.\n']);

    fprintf('\n4. WERE THE EDGES CLEAN?\n');
    for i = 1:numel(edges)
        E = edges(i);
        if ~E.visible
            if ~isempty(E.role)
                fprintf('   %-11s %-4s  NOT DETERMINED -- %s\n', E.channel, E.kind, E.note);
            end
            continue
        end
        if E.tau > 0.2
            note = 'slow response (NV physics); the switch time is the start of it';
        elseif E.tSwitchSterr > 0.01
            note = 'noisy -- raise averages before trusting it';
        else
            note = 'clean';
        end
        fprintf('   %-11s %-4s  tau %10s,  +/- %.4f us   %s\n', E.channel, E.kind, ...
                tauStr(E.tau), E.tSwitchSterr, note);
    end
    fprintf('\n============================================================\n\n');
end

% ======================================================================
% SMALL HELPERS
% ======================================================================

function d = checkDelayAt(d)
    if ischar(d) || (isstring(d) && isscalar(d))
        d = char(d);
        if ~strcmpi(d, 'onset')
            error('calibChannelDelays: delayAt must be ''onset'' or a fraction between 0 and 1.');
        end
        d = 'onset';
    elseif ~(isnumeric(d) && isscalar(d) && d > 0 && d < 1)
        error('calibChannelDelays: delayAt must be ''onset'' or a fraction between 0 and 1.');
    end
end

function s = tauStr(x)
    if ~isfinite(x)
        s = '-';
    elseif abs(x) < 0.2
        s = sprintf('%.1f ns', 1e3*x);
    else
        s = sprintf('%.3f us', x);
    end
end

function out = ternary(c, a, b)
    if c; out = a; else; out = b; end
end

function s = storedDelays()
    % Pull the on/off delays out of the setup json, wherever they live: lasers keep
    % theirs under 'switch', microwave switches under 'microwave', and the AWG trigger
    % under the frequency generator as onDelay1/offDelay1. Each entry also carries the
    % json key it came from, so a disagreement names the exact field to edit.
    s = struct();
    try
        j = JsonInfoReader.getJson();
    catch
        return
    end
    % jsondecode returns a 1x1 STRUCT for a single-element JSON array and a CELL for a
    % heterogeneous one. Handle both shapes.
    try
        k = 0;
        for L = asList(j.lasers)
            LL = L{1}; k = k + 1;
            if isstruct(LL) && isfield(LL, 'switch')
                s.(matlab.lang.makeValidName(LL.switch.switchChannelName)) = ...
                    struct('on', LL.switch.onDelay, 'off', LL.switch.offDelay, ...
                           'onPath',  sprintf('lasers(%d).switch.onDelay', k), ...
                           'offPath', sprintf('lasers(%d).switch.offDelay', k));
            end
        end
    catch
    end
    try
        f = fieldnames(j.microwave);
        for i = 1:numel(f)
            m = j.microwave.(f{i});
            s.(matlab.lang.makeValidName(m.switchChannelName)) = ...
                struct('on', m.onDelay, 'off', m.offDelay, ...
                       'onPath',  sprintf('microwave.%s.onDelay', f{i}), ...
                       'offPath', sprintf('microwave.%s.offDelay', f{i}));
        end
    catch
    end
    try
        k = 0;
        for G = asList(j.frequencyGenerators)
            g = G{1}; k = k + 1;
            if isfield(g, 'trigger1ChannelName')
                offd = 0;
                if isfield(g, 'offDelay1'); offd = g.offDelay1; end
                s.(matlab.lang.makeValidName(g.trigger1ChannelName)) = ...
                    struct('on', g.onDelay1, 'off', offd, ...
                           'onPath',  sprintf('frequencyGenerators(%d).onDelay1', k), ...
                           'offPath', sprintf('frequencyGenerators(%d).offDelay1', k));
            end
        end
    catch
    end
end

function e = blankEdge(channel, kind, commanded)
    e = struct('channel', channel, 'kind', kind, 'commanded', commanded, 'role', '', ...
               'dip', 0, 't0', NaN, 't0Sterr', NaN, 'tSwitch', NaN, 'tSwitchSterr', NaN, ...
               't90', NaN, 't90Sterr', NaN, 'amplitude', NaN, 'tau', NaN, 'tauSterr', NaN, ...
               'sigma', NaN, 'visible', false, 'compare', true, 'note', '');
end

function s = charOf(x)
    if iscell(x); s = [x{:}]; else; s = char(x); end
end

function o = withDefaults(o, d)
    if ~isstruct(o); error('calibChannelDelays: options must be a struct.'); end
    known = fieldnames(d);
    given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('calibChannelDelays: unknown option "%s". Valid: %s.', bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || (isempty(o.(known{i})) && ~strcmp(known{i}, 'data'))
            o.(known{i}) = d.(known{i});
        end
    end
end

function c = asList(x)
    % Whatever jsondecode produced, as a cell row.
    if iscell(x)
        c = reshape(x, 1, []);
    elseif isstruct(x)
        c = arrayfun(@(k) x(k), 1:numel(x), 'UniformOutput', false);
    else
        c = {};
    end
end
