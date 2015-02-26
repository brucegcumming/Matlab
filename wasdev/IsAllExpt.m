function true = IsAllExpt(E)
  true = isfield(E,'Expt') && isfield(E,'Spikes');
  