function c = t1Contrast(args)
%T1CONTRAST The |0> minus |1> contrast of a T1 run, oriented from the DATA rather than the labels.
%
%   c = t1Contrast(struct('exp', t))
%   c = t1Contrast(struct('value', y0, 'sterr', e0, 'value2', y1, 'sterr2', e1, 'tau', tau))
%   c = t1Contrast(struct('grab', g))          % an expGrab struct
%
% Returns a struct: tau, contrast, cerr, flipped, n, and which projection it decided is which.
%
% ===================================================================
% WHY THE ORIENTATION CANNOT BE TAKEN FROM THE LABELS
% ===================================================================
%
% AWGExpT1's projection labels were the wrong way round. LoadAWG builds [zero, noPi, zero, pi]
% with triggers [1 0 1 0], so the FIRST pulse-generator half plays the zeroed pi -- the spin
% stays in |0> through the wait -- and the second plays the real pulse. Detection rows 1-2 are
% therefore the |0> relaxation and rows 3-4 the |1> one, which is the opposite of Echo and
% Ramsey; the labels had been copied from those and announced rows 1-2 as |1>.
%
% That is fixed in the class. But every file that forms a contrast has to give the same answer
% on data taken BEFORE the fix, on a machine where the fix has not been copied over yet, and on
% saved runs from either. Reading the sign off the field name cannot do that: it flips the sign
% of the contrast, and downstream that flips the sign of a leak rate or turns a T1 fit into a
% rising exponential. Silently, and with confidence.
%
% So the orientation is MEASURED. In a T1 the |1> trace moves with tau and the |0> trace is the
% flat control, so whichever trace varies more is |1> and the contrast is flat minus varying.
% Same rule the class's own checkProjectionLabels uses, in one place both can share.
%
% The variation is the 10-to-90 percentile range rather than max minus min, so one dropped point
% -- and the classic path does drop points, diluted toward no contrast -- cannot decide which
% trace is which.
%
% WHEN IT CANNOT TELL, IT SAYS SO. If the two spreads are within 30% of each other there is no
% honest basis for choosing: either both traces are flat because T1 is far longer than the tau
% range, or the run is too noisy, or something is wrong with one readout window. c.confident
% comes back false, the field order is used as the fallback, and the caller is expected to look
% at c.why rather than at the contrast.

    DEF = struct('exp', [], 'grab', [], 'value', [], 'sterr', [], ...
                 'value2', [], 'sterr2', [], 'tau', []);
    if ~exist('args', 'var') || isempty(args); args = struct(); end
    if ~isstruct(args); error('t1Contrast takes ONE struct.'); end
    o = DEF;
    fn = fieldnames(args);
    for i = 1:numel(fn)
        if ~isfield(DEF, fn{i})
            error('t1Contrast: unknown field "%s". Valid: %s.', fn{i}, strjoin(fieldnames(DEF)', ', '));
        end
        o.(fn{i}) = args.(fn{i});
    end

    c = struct('tau', [], 'contrast', [], 'cerr', [], 'flipped', false, 'n', 0, ...
               'confident', false, 'why', '', 'spreadA', NaN, 'spreadB', NaN);

    % ---- where the numbers come from -------------------------------------
    src = o;
    if ~isempty(o.grab)
        src.value  = o.grab.value;  src.sterr  = o.grab.sterr;
        src.value2 = o.grab.value2; src.sterr2 = o.grab.sterr2;
        if isempty(src.tau); src.tau = o.grab.x; end
    elseif ~isempty(o.exp)
        [src.value,  src.sterr]  = paramData(o.exp, 'signalParam');
        [src.value2, src.sterr2] = paramData(o.exp, 'signalParam2');
        if isempty(src.tau); src.tau = tryGet(o.exp, 'tau', []); end
    end

    a = src.value(:)';  ea = src.sterr(:)';
    b = src.value2(:)'; eb = src.sterr2(:)';
    tau = src.tau(:)';
    if isempty(a) || isempty(b)
        c.why = ['only one projection is present. A T1 contrast needs both -- and note that with ' ...
                 'doubleMeasurement off, AWGExpT1 runs only the pi half, so there is nothing to ' ...
                 'take a difference against.'];
        return
    end
    n = min([numel(a), numel(b), numel(tau)]);
    if n < 4
        c.why = sprintf('only %d overlapping points between the two projections and tau.', n);
        return
    end
    a = a(1:n); b = b(1:n); tau = tau(1:n);
    if numel(ea) < n; ea = zeros(1, n); else; ea = ea(1:n); end
    if numel(eb) < n; eb = zeros(1, n); else; eb = eb(1:n); end

    % ---- decide which trace is |1> ---------------------------------------
    va = spread(a); vb = spread(b);
    c.spreadA = va; c.spreadB = vb;
    ratio = max(va, vb) / max(min(va, vb), eps);
    c.confident = ratio > 1.3;
    if vb >= va
        flat = a; fe = ea; vary = b; ve = eb; c.flipped = false;
    else
        flat = b; fe = eb; vary = a; ve = ea; c.flipped = true;
    end
    if ~c.confident
        % No basis for choosing. Fall back to the field order -- value is |0> in the patched
        % class -- and say why, loudly enough that a caller checking c.confident will not miss it.
        flat = a; fe = ea; vary = b; ve = eb; c.flipped = false;
        c.why = sprintf(['the two traces vary by within %.0f%% of each other (%.2g vs %.2g), so ' ...
                         'which one is |1> cannot be told from the data. Either T1 is far longer ' ...
                         'than the tau range, or the run is too noisy. Using the field order; do ' ...
                         'not trust the sign.'], 100*(ratio-1), va, vb);
    end

    c.tau = tau;
    c.contrast = flat - vary;
    c.cerr = sqrt(fe.^2 + ve.^2);
    c.n = n;
    if isempty(c.why) && median(c.contrast) < 0
        c.confident = false;
        c.why = ['the contrast is negative on most points. That is not an orientation problem -- ' ...
                 'this function fixes those -- so the two traces are not a |0>/|1> pair at all. ' ...
                 'Check doubleMeasurement, and that both signalParams actually filled.'];
    end
end

% ======================================================================

function v = spread(y)
% 10-to-90 percentile range. Robust to the isolated points the classic path dilutes toward no
% contrast: with max-minus-min, one such point on the FLAT trace can make it look like the one
% that varies, and then every sign downstream is wrong.
    y = y(isfinite(y));
    if numel(y) < 3; v = 0; return; end
    s = sort(y);
    lo = s(max(1, round(0.10*numel(s))));
    hi = s(min(numel(s), max(1, round(0.90*numel(s)))));
    v = hi - lo;
end

function [v, s] = paramData(e, which)
    v = []; s = [];
    p = tryGet(e, which, []);
    if isempty(p); return; end
    v = tryGet(p, 'value', []); s = tryGet(p, 'sterr', []);
    v = v(:)'; s = s(:)';
end

function v = tryGet(o, name, dflt)
    v = dflt;
    try
        if isobject(o) || isstruct(o); v = o.(name); end
    catch
    end
    if isempty(v) && ~isempty(dflt); v = dflt; end
end
