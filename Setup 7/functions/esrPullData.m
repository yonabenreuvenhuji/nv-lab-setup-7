function [f, y, e, reps, nAvg, info] = esrPullData(src, opts)
%ESRPULLDATA Get (frequency, signal, error) out of an ESR experiment, struct or matrix.
%
%   [f, y, e] = esrPullData(a)                     % an AWGExpESR
%   [f, y, e] = esrPullData(prevSurveyOutput)      % a struct with the vectors in it
%   [f, y, e] = esrPullData([freq signal sterr])   % a matrix
%   [f, y, e, reps, nAvg, info] = esrPullData(a, struct('verbose', false))
%
% One place that knows how to get a spectrum out of whatever the caller has. It exists
% because two different functions needed it and got different answers: esrFitRobust could
% read an experiment, calibFieldSurvey's replay could only read a struct or a matrix, and
% `calibFieldSurvey(struct('replay', a))` therefore failed with
%
%     replay must be a struct or a matrix, got a AWGExpESR
%
% which is a documentation bug as much as a code one -- the run script and the level 3
% notes both said an experiment was accepted. Rather than copy twenty lines of extraction
% into the second caller, both now come here.
%
% ===================================================================
% WHY NOT extract1Ddata, AND NOT extractESRdata
% ===================================================================
%
% extractESRdata is a four-line wrapper whose body is a call to extract1Ddata (or
% extract2Ddata when mirrorSweepAround is set), so the objection is really to that.
%
% extract1Ddata applies three successive cuts to the ratio across averages, at 5, then 3,
% then ONE standard deviation:
%
%       s(abs(s - mean(s,2)) > 1*mean(std(s,1,2))) = NaN;
%
% A one-sigma cut discards about a third of perfectly good Gaussian data. It also runs
% after two harsher passes have already shrunk the spread it measures against, and it
% biases the survivors toward the mean -- so the scatter comes out smaller than it is and
% every error bar downstream is too small.
%
% Here nothing is discarded unless there are at least 5 averages to judge with, and then on
% a median-absolute-deviation rule applied once, to whole averages rather than to individual
% points.
%
% ===================================================================
% THE ERROR BAR IS THE LARGER OF TWO THINGS
% ===================================================================
%
% The propagated error answers "how many photons was this point built from". The scatter of
% the same point across averages answers "how much did everything drift between them". They
% are different questions and on this setup the second is routinely an order of magnitude
% bigger -- the same drift that made AWGExpESR's sterr read 42x the point-to-point scatter.
% Taking the larger is the only choice that does not understate the uncertainty.
%
% The scatter enhancement is OPTIONAL at runtime, not because it is doubtful but because
% reaching for the per-average data is the fragile part. An earlier version of the Rabi
% fitter had only that path and died with "Dot indexing is not supported for variables of
% this type" when the array came back a different shape. A fit that cannot run because of
% how it reached for the data is worse than one with slightly optimistic error bars, so a
% failure here is reported and the propagated numbers are used.
%
% OPTIONS
%   avgOutlierMad   drop an average whose ratio is this many MADs from the median, and only
%                   with 5+ averages to judge from (default 4; Inf keeps everything)
%   verbose         say which error bars were used (default true)
%   name            caller name for error messages (default 'esrPullData')
%
% OUTPUT
%   f, y, e         column vectors. e is NaN where no error was available.
%   reps, nAvg      repeats and averages, 1 if not readable
%   info.source     'experiment' | 'struct' | 'matrix'
%   info.dropped    indices of averages dropped by the MAD rule
%   info.enhanced   true if the scatter across averages was folded in
%
% Base MATLAB only.

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('avgOutlierMad', 4, 'verbose', true, ...
                                     'name', 'esrPullData'), 'esrPullData');

    reps = 1; nAvg = 1;
    info = struct('source', '', 'dropped', [], 'enhanced', false);

    % ---- a matrix -------------------------------------------------------
    if isnumeric(src)
        if size(src, 2) < 2
            error('%s: a replay matrix needs at least 2 columns (freq, signal).', opts.name);
        end
        f = src(:,1); y = src(:,2);
        if size(src, 2) >= 3; e = src(:,3); else; e = nan(size(f)); end
        info.source = 'matrix';
        return
    end

    % ---- a struct -------------------------------------------------------
    % isstruct FIRST, before the object branch. An experiment is not a struct, but a struct
    % can carry fields named like an object's properties, and testing in the other order
    % would send a plain replay struct down the experiment path.
    if isstruct(src)
        f = pickField(src, {'freqs', 'freq', 'frequency', 'frequencies', 'x'});
        y = pickField(src, {'signal', 'value', 'values', 'y'});
        e = pickField(src, {'sterr', 'err', 'error', 'errors', 'e'});
        if isempty(f) || isempty(y)
            error(['%s: this struct has no recognisable sweep in it. It needs freqs and ' ...
                   'signal (sterr optional) -- the fields a previous survey returns.'], ...
                  opts.name);
        end
        f = f(:); y = y(:);
        if isempty(e); e = nan(size(f)); else; e = e(:); end
        if numel(y) ~= numel(f) || numel(e) ~= numel(f)
            error('%s: the vectors have different lengths (%d, %d, %d).', ...
                  opts.name, numel(f), numel(y), numel(e));
        end
        info.source = 'struct';
        return
    end

    % ---- an experiment object -------------------------------------------
    if ~isobject(src)
        error(['%s: cannot read a spectrum out of a %s. Pass an ESR experiment, a struct ' ...
               'with freqs/signal/sterr, or an [f y e] matrix.'], opts.name, class(src));
    end
    info.source = 'experiment';

    f = [];
    try
        ax = src.mCurrentXAxisParam;
        if ~isempty(ax); f = ax.value(:); end
    catch
    end
    if isempty(f)
        try; f = src.frequency(:); catch; end
    end
    if isempty(f)
        error(['%s: could not read the frequency axis off this %s. Has the experiment ' ...
               'run?'], opts.name, class(src));
    end

    try
        y = src.signalParam.value(:);
        e = src.signalParam.sterr(:);
    catch err
        error(['%s: could not read signalParam off this %s (%s). Has the experiment ' ...
               'actually run?'], opts.name, class(src), err.message);
    end
    if isempty(y)
        error('%s: signalParam.value is empty -- the experiment has not run.', opts.name);
    end
    try; reps = src.repeats;  catch; end
    try; nAvg = src.averages; catch; end

    % The frequency axis and the signal can disagree in length when a sweep was interrupted
    % part way. Trim rather than error: a partial sweep is still worth fitting, and the
    % alternative is losing a run that took twenty minutes.
    n = min(numel(f), numel(y));
    if numel(f) ~= numel(y)
        if opts.verbose
            fprintf(['%s: %d frequencies but %d signal points -- using the first %d. A sweep ' ...
                     'that\n              stopped early still fits.\n'], ...
                    opts.name, numel(f), numel(y), n);
        end
        f = f(1:n); y = y(1:n);
        if numel(e) >= n; e = e(1:n); else; e = nan(n,1); end
    end
    if numel(e) ~= numel(y); e = nan(size(y)); end

    % ---- fold in the scatter across averages ----------------------------
    try
        s = ratiosPerAverage(src);
        if ~isempty(s) && size(s, 2) >= 2
            nA = size(s, 2);
            if nA >= 5 && isfinite(opts.avgOutlierMad)
                m   = median(s, 2, 'omitnan');
                dev = abs(s - m);
                mad = 1.4826 * median(dev(:), 'omitnan');
                badAvg = mean(dev > opts.avgOutlierMad*mad, 1, 'omitnan') > 0.5;
                if any(badAvg) && ~all(badAvg)
                    info.dropped = find(badAvg);
                    s(:, badAvg) = [];
                    nA = size(s, 2);
                end
            end
            scatterErr = std(s, 0, 2, 'omitnan') ./ sqrt(max(nA, 1));
            k = min(numel(scatterErr), numel(e));
            e(1:k) = max(e(1:k), scatterErr(1:k));
            info.enhanced = true;
            if opts.verbose
                fprintf(['%s: error bars are max(propagated, scatter across %d averages).\n'], ...
                        opts.name, nA);
                if ~isempty(info.dropped)
                    fprintf('              dropped average(s) %s on a %g-MAD rule\n', ...
                            mat2str(info.dropped), opts.avgOutlierMad);
                end
            end
        end
    catch err
        if opts.verbose
            fprintf(['%s: propagated error bars only -- the per-average ratios were not\n' ...
                     '              reachable (%s). Not fatal, but they may be optimistic.\n'], ...
                    opts.name, err.message);
        end
    end
end

% ======================================================================

function s = ratiosPerAverage(Exp)
% nPoints x nAverages of signal/reference. Tries the documented accessor first, then the
% raw array, and returns empty rather than throwing if neither is shaped as expected.
    s = [];
    try
        c = Exp.alternateSignal();
        if iscell(c) && numel(c) >= 2
            a = c{1}.value; b = c{2}.value;
            if ~isempty(a) && isequal(size(a), size(b)) && size(a, 2) >= 2
                s = a ./ b; return
            end
        end
    catch
    end
    try
        raw = Exp.signal_;                      % nPoints x nRepeats x nChannels
        if ndims(raw) == 3 && size(raw, 3) >= 2
            S1 = squeeze(raw(:, :, 1));
            S2 = squeeze(raw(:, :, 2));
            if isequal(size(S1), size(S2)) && size(S1, 2) >= 2
                s = S1 ./ S2;
            end
        end
    catch
    end
end

function val = pickField(s, names)
    val = [];
    for k = 1:numel(names)
        if isfield(s, names{k}) && ~isempty(s.(names{k})) && isnumeric(s.(names{k}))
            val = s.(names{k}); return
        end
    end
end

function o = withDefaults(o, d, who)
    if isempty(o); o = struct(); end
    if ~isstruct(o); error('%s: options must be a struct.', who); end
    known = fieldnames(d); given = fieldnames(o);
    bad = given(~ismember(given, known));
    if ~isempty(bad)
        error('%s: unknown option "%s". Valid: %s.', who, bad{1}, strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(o, known{i}) || isempty(o.(known{i}))
            o.(known{i}) = d.(known{i});
        end
    end
end
