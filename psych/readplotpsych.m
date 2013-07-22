function [Data, Trials, label, h] = readplotpsych(prefix, file, prop,varargin)

Trials = [];
PSYCH.opt.forceread = 0;
PSYCH.opt.first = 1;
PSYCH.opt.last = 0;
PSYCH.opt.nmin = 2;
PSYCH.opt.xmin = -100;
PSYCH.opt.xmax = 100;
color = 'r';
onesd = [];

if ~isempty(prefix)
  file = [prefix file];
end

j = 1;
if ~isempty(prop.skip)
  PSYCH.opt.first = prop.skip;
end
if ~isempty(prop.forceread)
  PSYCH.opt.forceread = prop.forceread;
end
if ~isempty(prop.last)
  PSYCH.opt.last = prop.last;
end
if ~isempty(prop.sd)
  onesd = prop.sd;
end
if ~isempty(prop.xmin)
  PSYCH.opt.xmin = prop.xmin;
end
if ~isempty(prop.nmin)
  PSYCH.opt.nmin = prop.nmin;
end
if ~isempty(prop.xmax)
  PSYCH.opt.xmax = prop.xmax;
end

if ~isempty(prop.color)
  color = prop.color;
end

j = 1;
while(j < nargin-2)
  if(strcmpi(varargin{j},'showinit'))
    showinit = 1;
  elseif(strncmpi(varargin{j},'twomeans',6))
    fitopt.twomeans = 1;
  elseif(strncmpi(varargin{j},'twosd',5))
    fitopt.twosd = 1;
  elseif(strncmpi(varargin{j},'twofits',6))
    fitopt.twofits = 1;
  end
  j = j+1;
end

[Data, Trials] = readpsychfile(PSYCH,file);
if isempty(Data)
    fprintf('Read Failed for %s\n',file);
    return;
end

if ~isempty(onesd)
     Data = Data(find([Data.sd] == onesd));
end
 
for n = unique([Data.expno])
    fit = fitpsf(Data,'xmin',PSYCH.opt.xmin,'xmax',PSYCH.opt.xmax,'expno',n);
if(~isnan(fit.fit(1)))
h = plotpsych(fit.data,fit.fit(1),fit.fit(2),'color',color);
  label = sprintf('%s %.3g',fit.data(1).name,fit.fit(2));
end
end


