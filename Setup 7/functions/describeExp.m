function describeExp(Exp, wanted)
%DESCRIBEEXP Say WHICH property access failed, when reading an experiment object throws.
%
%   describeExp(Exp)
%   describeExp(Exp, {'tau', 'signalParam'})
%
% "Dot indexing is not supported for variables of this type" names neither the variable nor
% the field. On an unfamiliar experiment class that is nearly useless -- it cost two runs on
% the setup before this diagnostic existed. So when a fitter cannot read an experiment, it
% walks the properties it needs and reports the class and size of each; one run then pins it.
%
% This lives in its own file because rabiFitRobust CALLED IT WITHOUT DEFINING IT. The call
% sat in the catch block that is only reached when data extraction fails, so the bug was
% invisible until the day something went wrong -- and on that day the helpful diagnostic was
% replaced by "undefined function 'describeExp'", which is worse than the error it was
% written to explain. Untested error paths are error paths that do not work.

    if ~exist('wanted', 'var') || isempty(wanted)
        wanted = {'tau', 'frequency', 'signalParam', 'signalParam2', 'signal', 'sterr', ...
                  'currIter', 'averages', 'repeats', 'halfPiTime', 'piTime', ...
                  'doubleMeasurement', 'mCurrentXAxisParam'};
    end
    fprintf('\nCould not read the experiment. What is actually there:\n');
    fprintf('  class            %s\n', class(Exp));
    for i = 1:numel(wanted)
        nm = wanted{i};
        try
            v = Exp.(nm);
            fprintf('  %-18s %s, size %s\n', nm, class(v), mat2str(size(v)));
        catch
            fprintf('  %-18s NOT READABLE\n', nm);
        end
    end
    fprintf(['\nIf the properties above look right, the fitter is asking for something else --\n' ...
             'read the error below for the field name.\n']);
end
