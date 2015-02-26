function Expt = InsertTrials(Expt, Trials, varargin)
% Expt = InsertTrials(Expt, Trials, inserts new trials into Expt
if isfield(Expt.Trials,'Trial')
    sortby = 'Trial';
elseif isfield(Expt.Trials,'id')
    sortby = 'id'l
else
    sortby = 'start';
end

for j = 1:length(Trials)
    Expt.Trials = CopySFields(AllExpt.Expt.Trials,'new',Trials(j));
end

X = expt.TrialVals(Expt,sortby);


