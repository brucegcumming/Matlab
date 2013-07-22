function [Data, label, h] = readandplotpsych(file,varargin)

PSYCH.opt.forceread = 0;
PSYCH.opt.first = 1;
PSYCH.opt.last = 0;
PSYCH.opt.nmin = 2;
PSYCH.opt.xmin = -100;
PSYCH.opt.xmax = 100;
color = 'r';
onesd = [];

j = 1;
while(j < nargin)
  if(strcmpi(varargin{j},'showinit'))
    showinit = 1;
  elseif(strncmpi(varargin{j},'twomeans',6))
    fitopt.twomeans = 1;
  elseif(strncmpi(varargin{j},'twosd',5))
    fitopt.twosd = 1;
  elseif(strncmpi(varargin{j},'xmax',4))
    j = j+1;
    PSYCH.opt.xmax = varargin{j};
  elseif(strncmpi(varargin{j},'xmin',4))
    j = j+1;
    PSYCH.opt.xmin = varargin{j};
  elseif(strncmpi(varargin{j},'first',4))
    j = j+1;
    PSYCH.opt.first = varargin{j};
  elseif(strncmpi(varargin{j},'sd',2))
    j = j+1;
    onesd = varargin{j};
  elseif(strncmpi(varargin{j},'last',4))
    j = j+1;
    PSYCH.opt.last = varargin{j};
  elseif(strncmpi(varargin{j},'color',4))
    j = j+1;
    color = varargin{j};
  elseif(strncmpi(varargin{j},'twofits',6))
    fitopt.twofits = 1;
  end
  j = j+1;
end




[Data, Trials] = readpsychfile(PSYCH,file);
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


