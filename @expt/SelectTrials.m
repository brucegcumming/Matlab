function [Expt, varargout] = SelectTrials(Expt, varargin)
%expt.SelectTrials(Expt, var, val, var, va...)
%restricts trials in Expt to those matching criteria
% [Expt.Trials.(var)] == val


if iscell(Expt)
    for j = 1:length(Expt)
        Expt{j} = expt.SelectTrials(Expt{j},varargin{:});
    end
    return;
end
j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'anything')
    else
        str = varargin{j};
        j = j+1;
        val = varargin{j};
        if isfield(Expt.Trials,str)
            id = find([Expt.Trials.(str)] == val);
            Expt.Trials = Expt.Trials(id);
        end
    end        
    j = j+1;
end