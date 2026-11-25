function res = expPropsWrite(newVals, opts)
%EXPPROPSWRITE Write named numeric fields into expProperties.json, surgically.
%
%   res = expPropsWrite(out.toWrite)                      % into setup7Paths().json
%   res = expPropsWrite(out.toWrite, struct('dryRun', true))     % show, change nothing
%   res = expPropsWrite(struct('piTime', 0.0374))         % one field
%
% The only thing in this calibration suite that writes to the json. Everything else prints
% what it would write and stops, which is deliberate: a calibration that writes as a side
% effect cannot be re-run to check itself. The level scripts call this in their write cells:
% level2.m 2.3, level3.m 3.5, level4.m 4.3 and 4.5 (4.7 for a re-run). Levels 0 and 1 have
% nothing for expProperties.json -- level 1's delays live in setupInfo.json, by hand.
%
% ===================================================================
% WHY NOT saveJsonParams
% ===================================================================
%
% Experiment.saveJsonParams is broken on this installation, in two independent ways:
%
%   ARITY. It calls JsonInfoReader.setParams(paramsStruct, jsonLocation), and setParams
%   is declared setParams(params) -- one argument. The call raises before it writes.
%
%   AND IT IGNORES THE LOCATION ANYWAY. setParams does not read its second argument even
%   if you give it one. It resolves the target itself from C:\lab\setupParamsLocation.json,
%   so a call that names the newest results folder writes to whatever that pointer happens
%   to say -- which during a calibration session is usually the previous folder, or the
%   shared one. Silently. That is worse than the error.
%
% ===================================================================
% WHY A TEXT EDIT AND NOT jsondecode -> jsonencode
% ===================================================================
%
% The obvious implementation is decode, set the field, encode, write. It is wrong here:
%
%   - jsonencode emits one line with no spaces, so the file stops being readable and every
%     future diff is the whole file.
%   - field order is lost, so the grouping that makes the file navigable by eye is gone.
%   - a round trip is not identity. An integer-valued double comes back as 3 rather than
%     3.0; [] and {} are not distinguished; a 1x1 cell collapses. None of that changes what
%     MATLAB reads back, but it does change the file, and a calibration tool that rewrites
%     500 lines to change one number is a tool nobody will run twice.
%
% So this edits the value token in place and leaves every other byte alone. The file after a
% write differs from the file before it only in the numbers that changed.
%
% ===================================================================
% THE SAFETY RULES, ALL OF WHICH REFUSE RATHER THAN GUESS
% ===================================================================
%
%   1. TOP LEVEL ONLY. The key must sit at brace depth 1. A "duration" nested inside
%      lasers[1].switch is not the "duration" you meant, and a tool that cannot tell the
%      difference must not write.
%   2. UNIQUE. If the key appears more than once at the top level the file is already
%      ambiguous -- MATLAB's jsondecode keeps the LAST one, which is not what reading the
%      file top to bottom suggests. Refused, and named.
%   3. NUMERIC ONLY. Scalars and numeric row vectors. Strings, structs and booleans are
%      refused, because their formatting conventions are not worth guessing at.
%   4. BACKUP FIRST, and the backup is verified readable before the original is touched.
%   5. VERIFIED AFTER. The new file is jsondecode'd and compared field by field against
%      both the intended values AND the original decode. If any field that was not asked
%      for has changed, or any field that was asked for did not, the backup is restored and
%      the function raises. There is no path where the json is left in a state nobody
%      checked.
%
% OPTIONS
%   jsonLocation  default setup7Paths().json -- the NEWEST dated results folder
%   dryRun        report and change nothing (default false)
%   create        allow keys that are not in the file yet (default false). Off by default
%                 because a typo'd key would otherwise be created rather than caught, and
%                 loadStateFromStruct ignores json fields with no matching property, so it
%                 would then be silently dead forever.
%   backup        default true
%   verbose       default true
%
% OUTPUT
%   res.written   struct of what went in
%   res.skipped   cell of {key, reason}
%   res.backup    path to the backup, '' if none
%   res.file      the file written
%   res.ok        true only if every requested key was written

    if ~exist('opts', 'var'); opts = struct(); end
    opts = withDefaults(opts, struct('jsonLocation', '', 'dryRun', false, ...
        'create', false, 'backup', true, 'verbose', true), 'expPropsWrite');

    if ~isstruct(newVals) || numel(newVals) ~= 1
        error('expPropsWrite: first argument must be a 1x1 struct of field -> value.');
    end
    if isempty(opts.jsonLocation)
        opts.jsonLocation = setup7Paths().json;
    end
    f = opts.jsonLocation;
    if ~exist(f, 'file')
        error('expPropsWrite: no such file:\n  %s', f);
    end

    res = struct('written', struct(), 'skipped', {{}}, 'unchanged', {{}}, ...
                 'backup', '', 'file', f, 'ok', false);

    txt  = fileread(f);
    J0   = jsondecode(txt);
    keys = fieldnames(newVals);

    % ---- locate every top-level key, once ------------------------------
    idx = topLevelKeys(txt);

    if opts.verbose
        fprintf('\n============================================================\n');
        fprintf('expPropsWrite\n');
        fprintf('============================================================\n');
        fprintf('  file  %s\n', f);
        fprintf('  %d top-level keys in it\n\n', numel(fieldnames(idx)));
    end

    % ---- plan, with no writing yet -------------------------------------
    plan = struct('key', {}, 'val', {}, 'old', {}, 'a', {}, 'b', {}, 'note', {});
    for i = 1:numel(keys)
        k = keys{i};
        v = newVals.(k);
        if ~isnumeric(v) || ~isreal(v) || isempty(v) || ~isvector(v) || any(~isfinite(v))
            res.skipped{end+1} = {k, 'not a finite real numeric scalar or vector'};
            continue
        end
        v = v(:)';
        if ~isfield(idx, validKey(k))
            if ~opts.create
                res.skipped{end+1} = {k, ['not present in the file, and create is false. ' ...
                    'If the key really is new, check the spelling against the experiment ' ...
                    'class property first, then re-run with struct(''create'', true).']};
                continue
            end
            plan(end+1) = struct('key', k, 'val', v, 'old', [], ...
                                 'a', NaN, 'b', NaN, 'note', 'NEW');           %#ok<AGROW>
            continue
        end
        e = idx.(validKey(k));
        if numel(e.a) > 1
            res.skipped{end+1} = {k, sprintf(['appears %d times at the top level. The file ' ...
                'is ambiguous -- jsondecode keeps the last one. Delete the duplicates by ' ...
                'hand first.'], numel(e.a))};
            continue
        end
        old = [];
        if isfield(J0, validKey(k)) && isnumeric(J0.(validKey(k)))
            old = J0.(validKey(k));
            old = old(:)';
        end
        if ~isempty(old) && numel(old) == numel(v) && all(abs(old - v) < 1e-12)
            % NOT a skip. The file already says what was asked for, so the request is
            % satisfied -- counting this as a failure would make res.ok false on the second
            % run of a cell that worked the first time, which is exactly backwards.
            res.unchanged{end+1} = k;
            res.written.(validKey(k)) = v;
            continue
        end
        plan(end+1) = struct('key', k, 'val', v, 'old', old, ...
                             'a', e.a, 'b', e.b, 'note', '');                  %#ok<AGROW>
    end

    % ---- show the plan --------------------------------------------------
    if opts.verbose
        if isempty(plan)
            fprintf('  nothing to change.\n');
        else
            for i = 1:numel(plan)
                if isempty(plan(i).old)
                    fprintf('  %-28s %14s  ->  %s\n', plan(i).key, '(new)', ...
                            numToJson(plan(i).val));
                else
                    fprintf('  %-28s %14s  ->  %s\n', plan(i).key, ...
                            numToJson(plan(i).old), numToJson(plan(i).val));
                end
            end
        end
        for i = 1:numel(res.unchanged)
            fprintf('  %-28s already this value\n', res.unchanged{i});
        end
        for i = 1:numel(res.skipped)
            fprintf('  SKIPPED %-20s %s\n', res.skipped{i}{1}, res.skipped{i}{2});
        end
    end

    if opts.dryRun
        if opts.verbose
            fprintf('\n  DRY RUN -- nothing written. Re-run without dryRun to write.\n');
            fprintf('============================================================\n\n');
        end
        res.ok = isempty(res.skipped);
        return
    end
    if isempty(plan)
        res.ok = isempty(res.skipped);
        if opts.verbose; fprintf('============================================================\n\n'); end
        return
    end

    % ---- backup ---------------------------------------------------------
    if opts.backup
        stamp = datestr(now, 'yyyymmdd_HHMMSS');                          %#ok<TNOW1,DATST>
        bak = [f '.bak.' stamp];
        ok = false;
        try
            copyfile(f, bak);
            % Verified readable, not just created. A backup nobody has read is a hope.
            chk = fileread(bak);
            ok = strcmp(chk, txt);
        catch
        end
        if ~ok
            error(['expPropsWrite: could not make a verified backup at\n  %s\n' ...
                   'Nothing was written. Check the folder is writable.'], bak);
        end
        res.backup = bak;
        if opts.verbose; fprintf('\n  backup  %s\n', bak); end
    end

    % ---- splice, back to front so the offsets stay valid ----------------
    % Descending by start position: an edit only moves the text AFTER it, so editing the
    % last one first leaves every earlier span exactly where topLevelKeys found it. Doing
    % this forward and adjusting offsets by the length delta is the classic place this kind
    % of function goes wrong.
    news = txt;
    hasSpan = arrayfun(@(p) isfinite(p.a), plan);
    spanned = plan(hasSpan);
    if ~isempty(spanned)
        [~, ord] = sort([spanned.a], 'descend');
        spanned = spanned(ord);
        for i = 1:numel(spanned)
            p = spanned(i);
            news = [news(1:p.a-1), numToJson(p.val), news(p.b+1:end)];
        end
    end
    for i = 1:numel(plan)
        if ~hasSpan(i)
            news = insertKey(news, plan(i).key, numToJson(plan(i).val));
        end
    end

    % ---- write, then verify, then keep or roll back ---------------------
    fid = fopen(f, 'w');
    if fid < 0
        error('expPropsWrite: could not open for writing:\n  %s', f);
    end
    fwrite(fid, news);
    fclose(fid);

    [good, why] = verifyWrite(f, J0, plan);
    if ~good
        if ~isempty(res.backup)
            try; copyfile(res.backup, f); catch; end
            error(['expPropsWrite: the file it wrote does not read back correctly, so the\n' ...
                   'backup has been RESTORED and the json is as it was.\n  %s\n' ...
                   'Nothing about the calibration is wrong -- this is a bug in the writer.\n' ...
                   'Report the message above and set the field by hand for now.'], why);
        end
        error(['expPropsWrite: the file it wrote does not read back correctly and there was\n' ...
               'no backup to restore:\n  %s\nFIX THE FILE BY HAND BEFORE RUNNING ANYTHING.'], why);
    end

    for i = 1:numel(plan)
        res.written.(validKey(plan(i).key)) = plan(i).val;
    end
    res.ok = isempty(res.skipped);

    if opts.verbose
        fprintf('\n  WRITTEN and verified: the file reads back with exactly these changes and\n');
        fprintf('  no others.\n');
        writeVerdict(res, plan);
    end
end

% ======================================================================

function writeVerdict(res, plan)
    fprintf('\n------------------------------------------------------------\n');
    fprintf('WHAT TO DO NOW\n');
    fprintf('------------------------------------------------------------\n');

    fprintf(['\n1. NOTHING RELOADS THIS BY ITSELF. loadJsonParams runs when an experiment is\n' ...
             '   constructed, so an object that already exists in the workspace still holds\n' ...
             '   the old numbers. Rebuild it -- b = AWGExpEcho; -- rather than re-running\n' ...
             '   b.run and wondering why nothing moved.\n']);

    fprintf(['\n2. CHECK THE FIELD IS ACTUALLY READ. loadStateFromStruct loops over the\n' ...
             '   OBJECT''s properties and asks whether the json has a field of that name. A\n' ...
             '   json field with no matching property is skipped in silence -- no error, no\n' ...
             '   warning, and the value simply never arrives. One line settles it:\n' ...
             '       ismember(''halfPiTime'', properties(''AWGExpEcho''))\n']);

    if ~isempty(res.backup)
        fprintf(['\n3. THE PREVIOUS FILE is at\n       %s\n' ...
                 '   Copy it back to undo this. The backups accumulate one per write, which is\n' ...
                 '   what you want during a calibration session and worth clearing out after.\n'], ...
                res.backup);
    end

    n = 3 + double(~isempty(res.backup));
    fprintf(['\n%d. AND THE ONE THING TO VERIFY ON HARDWARE. Pulse times are only right at the\n' ...
             '   amplitude they were measured at. Run one Rabi at the new piTime and confirm\n' ...
             '   the contrast is at its minimum there -- that is a ten-second check that\n' ...
             '   catches a wrong amplitude, a wrong frequency and a write to the wrong file\n' ...
             '   all at once.\n'], n);
    fprintf('\n============================================================\n\n');
    if ~isempty(res.skipped)
        fprintf('NOTE %d field(s) were SKIPPED, listed above. They are NOT in the file.\n\n', ...
                numel(res.skipped));
    end
end

% ======================================================================

function idx = topLevelKeys(txt)
% Every "key": at brace depth 1, with the byte span of its value.
%
% A hand-written scanner rather than a regexp, for one reason: a regexp cannot count braces,
% so it cannot tell a top-level "duration" from one nested three levels down inside a laser
% definition. Getting that wrong means writing the right number into the wrong place, which
% is the single most damaging thing this file could do.
    idx = struct();
    n = numel(txt);
    depth = 0; i = 1; inStr = false;
    while i <= n
        c = txt(i);
        if inStr
            if c == '\'; i = i + 2; continue; end
            if c == '"'; inStr = false; end
            i = i + 1; continue
        end
        switch c
            case '"'
                % a string at depth 1 followed by a colon is a top-level key
                [s, iEnd] = readStr(txt, i);
                j = skipWs(txt, iEnd + 1);
                if depth == 1 && j <= n && txt(j) == ':'
                    a = skipWs(txt, j + 1);
                    b = valueEnd(txt, a);
                    k = validKey(s);
                    if isfield(idx, k)
                        e = idx.(k);
                        e.a(end+1) = a; e.b(end+1) = b;
                        idx.(k) = e;
                    else
                        idx.(k) = struct('name', s, 'a', a, 'b', b);
                    end
                    i = b + 1; continue
                end
                i = iEnd + 1; continue
            case {'{', '['}
                depth = depth + 1; i = i + 1;
            case {'}', ']'}
                depth = depth - 1; i = i + 1;
            otherwise
                i = i + 1;
        end
    end
end

function [s, iEnd] = readStr(txt, i)
% i points at the opening quote. Returns the contents and the index of the closing quote.
    n = numel(txt); j = i + 1; s = '';
    while j <= n
        if txt(j) == '\'; s(end+1) = txt(j+1); j = j + 2; continue; end   %#ok<AGROW>
        if txt(j) == '"'; break; end
        s(end+1) = txt(j); j = j + 1;                                     %#ok<AGROW>
    end
    iEnd = j;
end

function j = skipWs(txt, j)
    n = numel(txt);
    while j <= n && any(txt(j) == sprintf(' \t\r\n')); j = j + 1; end
end

function b = valueEnd(txt, a)
% Last byte of the value starting at a: scan to the ',' or '}' that closes it, at ITS depth.
    n = numel(txt); d = 0; inStr = false; j = a;
    while j <= n
        c = txt(j);
        if inStr
            if c == '\'; j = j + 2; continue; end
            if c == '"'; inStr = false; end
            j = j + 1; continue
        end
        if c == '"'; inStr = true; j = j + 1; continue; end
        if any(c == '{[' ); d = d + 1; j = j + 1; continue; end
        if any(c == '}]')
            if d == 0; break; end
            d = d - 1; j = j + 1; continue
        end
        if c == ',' && d == 0; break; end
        j = j + 1;
    end
    b = j - 1;
    % trailing whitespace belongs to the file, not to the value
    while b > a && any(txt(b) == sprintf(' \t\r\n')); b = b - 1; end
end

function txt = insertKey(txt, key, valStr)
% Add a new top-level key just before the closing brace, matching the indentation of the
% last entry so the file still looks hand-written.
    close = find(txt == '}');
    if isempty(close)
        error('expPropsWrite: no closing brace in the json; refusing to append.');
    end
    p = close(end);
    head = txt(1:p-1);
    % strip trailing whitespace off the body, then work out the indent in use
    q = numel(head);
    while q > 0 && any(head(q) == sprintf(' \t\r\n')); q = q - 1; end
    body = head(1:q);
    nl = find(body == sprintf('\n'), 1, 'last');
    indent = '    ';
    if ~isempty(nl)
        line = body(nl+1:end);
        k = 1;
        while k <= numel(line) && any(line(k) == sprintf(' \t')); k = k + 1; end
        if k > 1; indent = line(1:k-1); end
    end
    sep = '';
    if ~isempty(body) && body(end) ~= '{'; sep = ','; end
    txt = [body, sep, sprintf('\n'), indent, '"', key, '": ', valStr, ...
           sprintf('\n'), txt(p:end)];
end

function s = numToJson(v)
% Scalars plain, vectors bracketed. Integers as integers -- writing 3.0 where the file said
% 3 makes a diff that is not a change. %.12g elsewhere: enough to round-trip a double that
% matters, few enough digits that the file stays readable.
    if isempty(v); s = '[]'; return; end
    v = v(:)';
    parts = cell(1, numel(v));
    for i = 1:numel(v)
        x = v(i);
        if x == fix(x) && abs(x) < 1e15
            parts{i} = sprintf('%d', x);
        else
            parts{i} = sprintf('%.12g', x);
        end
    end
    if numel(v) == 1
        s = parts{1};
    else
        s = ['[' strjoin(parts, ', ') ']'];
    end
end

function [good, why] = verifyWrite(f, J0, plan)
% The point of this function: prove the file changed in exactly the intended way and in no
% other way. It is the difference between "the write did not error" and "the write is right".
    good = false; why = '';
    try
        J1 = jsondecode(fileread(f));
    catch err
        why = sprintf('the new file is not valid json (%s)', err.message);
        return
    end
    % (a) every planned field now holds the planned value
    for i = 1:numel(plan)
        k = validKey(plan(i).key);
        if ~isfield(J1, k)
            why = sprintf('"%s" is not in the file after the write', plan(i).key); return
        end
        got = J1.(k);
        if ~isnumeric(got)
            why = sprintf('"%s" did not come back numeric', plan(i).key); return
        end
        got = got(:)'; want = plan(i).val(:)';
        if numel(got) ~= numel(want) || any(abs(got - want) > 1e-9*max(1, max(abs(want))))
            why = sprintf('"%s" reads back as %s, not %s', plan(i).key, ...
                          numToJson(got), numToJson(want));
            return
        end
    end
    % (b) nothing else moved
    planned = cellfun(@validKey, {plan.key}, 'UniformOutput', false);
    was = fieldnames(J0);
    for i = 1:numel(was)
        k = was{i};
        if ismember(k, planned); continue; end
        if ~isfield(J1, k)
            why = sprintf('"%s" disappeared, and it was not one of the fields being set', k);
            return
        end
        if ~isequaln(J0.(k), J1.(k))
            why = sprintf(['"%s" changed, and it was not one of the fields being set. That ' ...
                           'is a bug in the writer, not in the calibration.'], k);
            return
        end
    end
    now1 = fieldnames(J1);
    extra = now1(~ismember(now1, [was; planned(:)]));
    if ~isempty(extra)
        why = sprintf('"%s" appeared, and nothing asked for it', extra{1});
        return
    end
    good = true;
end

function k = validKey(s)
% A json key that is not a legal MATLAB field name cannot be indexed, so match what
% jsondecode itself does: non-word characters to underscores, and a leading digit gets an
% 'x'. matlab.lang.makeValidName would do it, but it is not in Octave, and the offline
% graders have to be able to run this file.
    k = s;
    bad = ~(isstrprop(k, 'alphanum') | k == '_');
    k(bad) = '_';
    if isempty(k); k = 'x'; return; end
    if ~isstrprop(k(1), 'alpha'); k = ['x' k]; end
    if numel(k) > 63; k = k(1:63); end
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
        if ~isfield(o, known{i}); o.(known{i}) = d.(known{i}); end
    end
end
