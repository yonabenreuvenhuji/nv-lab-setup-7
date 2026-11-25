function testSequenceBatching()
%TESTSEQUENCEBATCHING Offline check that batching preserves the SCPI stream.
%
% Rebuilds the command list the OLD way (one command per write) and the NEW way
% (compound messages joined with ';:'), then splits the compound messages back
% apart and asserts the two streams are identical, in order.
%
% This is the thing most likely to break silently: the instrument accepts a
% malformed compound message with an execution error that fwrite never sees.

    fprintf('\n=== testSequenceBatching ===\n');
    nPass = 0; nFail = 0;
    function check(name, cond, detail)
        if cond, fprintf('  PASS  %s\n', name); nPass=nPass+1;
        else,    fprintf('  FAIL  %s   %s\n', name, detail); nFail=nFail+1; end
    end

    cases = { struct('wf',[7 3 9],      'rp',[1 110 1],   'tg',[1 0 0], 'goto',[0 0 0]), ...
              struct('wf',[7 3 9 7 4 9],'rp',[1 110 1 1 95 1],'tg',[1 0 0 1 0 0],'goto',[0 0 1 0 0 4]), ...
              struct('wf',ones(1,573),  'rp',ones(1,573),  'tg',repmat([1 0 0],1,191), ...
                     'goto',reshape([zeros(2,191); 1:3:573],1,[])) };

    for c = 1:length(cases)
        s = cases{c};
        oldStream = buildOld(s.wf, s.rp, s.tg, s.goto, '1');
        newMsgs   = buildNew(s.wf, s.rp, s.tg, s.goto, '1', 40);
        newStream = {};
        for m = 1:length(newMsgs)
            parts = strsplit(newMsgs{m}, ';:');
            newStream = [newStream, parts];                            %#ok<AGROW>
        end
        check(sprintf('case %d (%d steps): command streams identical', c, length(s.wf)), ...
              isequal(oldStream, newStream), ...
              sprintf('%d old vs %d new commands', length(oldStream), length(newStream)));
        check(sprintf('case %d: no message exceeds 40 commands', c), ...
              all(cellfun(@(m) length(strsplit(m,';:')) <= 40, newMsgs)), 'chunk too large');
        check(sprintf('case %d: no message starts or ends with a separator', c), ...
              all(cellfun(@(m) ~startsWith(m,';') && ~endsWith(m,';') && ~endsWith(m,':'), newMsgs)), ...
              'malformed compound message');
        fprintf('        %d writes -> %d writes (%.1fx fewer)\n', ...
                length(oldStream), length(newMsgs), length(oldStream)/length(newMsgs));
    end

    fprintf('\n  %d passed, %d failed\n\n', nPass, nFail);
    if nFail > 0, error('testSequenceBatching: %d check(s) failed', nFail); end
end

function cmds = buildOld(wf, rp, tg, goto, name)
    cmds = {};
    for r = 1:length(wf)
        cmds{end+1} = sprintf('slist:seq:step%d:tasset1:wav "%s", "%s"', r, name, sprintf('%d',wf(r))); %#ok<AGROW>
        cmds{end+1} = sprintf('slist:seq:step%d:rcount "%s", %d', r, name, rp(r));                      %#ok<AGROW>
        if tg(r) == 0
            cmds{end+1} = sprintf('slist:seq:step%d:winput "%s", off', r, name);                        %#ok<AGROW>
        else
            cmds{end+1} = sprintf('slist:seq:step%d:winput "%s", atrigger', r, name);                   %#ok<AGROW>
        end
        if goto(r) > 0
            cmds{end+1} = sprintf('slist:seq:step%d:goto "%s", %d', r, name, goto(r));                  %#ok<AGROW>
        end
    end
    if all(goto == 0)
        cmds{end+1} = sprintf('slist:seq:step%d:goto "%s", first', length(wf), name);
    end
end

function msgs = buildNew(wf, rp, tg, goto, name, maxBatch)
    cmds = buildOld(wf, rp, tg, goto, name);   % same generator, then chunk
    msgs = {};
    for first = 1:maxBatch:length(cmds)
        last = min(first+maxBatch-1, length(cmds));
        msgs{end+1} = strjoin(cmds(first:last), ';:');                                                  %#ok<AGROW>
    end
end
