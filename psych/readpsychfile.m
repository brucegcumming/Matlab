function [Data, Trials, blocklist] = readpsychfile(PSYCH, name)

blocklist = [];
Data = [];
Trials = [];
matname = sprintf('%s.mat',name);
if exist(matname,'file') & ~(PSYCH.opt.forceread > 0)
  load(matname);
  Data = CountPsychTrials(Trials,'sortby','sd','nmin', ...
			  PSYCH.opt.nmin,'skip',PSYCH.opt.first,'last',PSYCH.opt.last);
  return;
end

fid = fopen(name);
if isempty(fid) | fid < 0
  fprintf('Cannot open %s\n',name);
  return;
end

tline = fgetl(fid);
fclose(fid);
[a,b] = sscanf(tline,'%d %f');


if(strncmp(tline,'Reopened',8)) % this is a binoc output file
  Trials = readpsychtrials(name,'save');
  Data = CountPsychTrials(Trials,'sortby','sd','nmin', ...
			  PSYCH.opt.nmin,'skip',PSYCH.opt.first,'last',PSYCH.opt.last);
elseif length(a) & length(b)  %% a table of numbers, suggesting a online table
  tr = textread(name);
  id = find(tr(:,4) > 0 | tr(5,:) > 0)  %% made a response
  Trials.score = tr(id,4);
  Trials.sign = sign(tr(id,6));
  Trials.x = tr(id,6);
  Data = CountPsychTrials(Trials,'nmin', ...
			  PSYCH.opt.nmin,'skip',PSYCH.opt.first,'last',PSYCH.opt.last);
        
else
  Data = readpsychsum(name);
end
