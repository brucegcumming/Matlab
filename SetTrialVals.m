function Trials = SetTrialVals(Expt,type)


Trials = Expt.Trials;
if(isfield(Expt.Trials,type))
     return;
elseif eval(['~isempty(Expt.Stimvals.' type ')']);
  val = eval(['Expt.Stimvals.' type]);
else
  fprintf('Can''t set %s for %s - no values or Stimvals set\n',Expt.Header.name,type);
     return;
end
eval(['[Trials.' type '] = deal(val);']);

