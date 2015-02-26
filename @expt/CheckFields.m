function E = CheckFields(E,varargin)
%Make sure widely used fields are preent in Expt/AllExpt struct


if isfield(E,'Spikes') && iscell(E.Spikes) %An AllExpt
    E.Expt = expt.CheckFields(E.Expt);
    return;
elseif iscell(E)
end

for j = 1:length(E.Trials)
    if ~isfield(E.Trials,'dur') || isempty(E.Trials(j).dur)
        E.Trials(j).dur = E.Trials(j).End(end)-E.Trials(j).Start(1);
    end
    if ~isfield(E.Trials,'id') || isempty(E.Trials(j).id)
        E.Trials(j).id = j;
    end
end