function [sdf, nsac, nspikes, saccades] = mksacsdf(Trials, smooth, times, ...
					 varargin)

%
%[sdf, nsac, nspikes, saccades] = mksacsdf(Trials, smooth, times, type, varargin)
%
% makes saccade triggered spike density functions.
% 
%
maxsize = inf;
minsize = 0;
mincos = 0.5;
axis = NaN;
direction = NaN;
nvar = nargin - 4;
trigtype = 'start';
j = 1;
while j  <= length(varargin)
  str = varargin{j};
  if strncmpi(str,'MaxSize',4)
    maxsize = varargin{j+1};
    j = j+1;
  elseif strncmpi(str,'MinSize',4)
    minsize = varargin{j+1};
    j = j+1;
  elseif strncmpi(str,'Axis',3)
      axis = varargin{j+1};
    j = j+1;
  elseif strncmpi(str,'direction',6)
      direction = varargin{j+1};
    j = j+1;
  elseif strncmpi(str,'Trig',4)
      j = j+1;
    trigtype = varargin{j};
  elseif strncmpi(str,'DirWidth',4)
    dirwidth = varargin{j+1};
    j = j+1;
    mincos = cos(dirwidth);
  end
  j = j+1;
end

c.times = times;
if IsAllExpt(Trials)
    details.times = times;
    for j = 1:length(Trials.Header)
        Expt = All2Expt(Trials,Trials.Header(j).cellnumber);
        [sdf(:,j), nsac, nspikes, c] = mksacsdf(Expt.Trials,smooth, times, varargin{:});
        details.nsac = nsac;
        details.nspikes(j) = nspikes;
        details.cellnumber(j) = Trials.Header(j).cellnumber;
        details = CopyFields(details,c);
    end
    nsac = details;
    return;
end
sacamp = [];
sacdir = [];
sacdur = [];
sacpeak = [];
sacstarts = [];
for tr = 1:length(Trials);
  Saccades = Trials(tr).Saccades;
  if isempty(Saccades)
      sidx = [];
  else
      sidx = find([Saccades.start] > 0 & [Saccades.start] < Trials(tr).dur);
      Saccades = Saccades(sidx);
  end
  
  if isempty(Saccades)
      sidx = [];
  elseif ~isnan(axis) 
    dirdev = cos([Saccades.dir] - axis);
    sidx = find([Saccades.size] > minsize & [Saccades.size] < ...
		maxsize & abs(dirdev) > mincos);
  elseif ~isnan(direction) 
    dirdev = cos([Saccades.dir] - direction);
    sidx = find([Saccades.size] > minsize & [Saccades.size] < ...
		maxsize & dirdev > mincos);
  else
    sidx = find([Saccades.size] > minsize & [Saccades.size] < ...
		maxsize);
  end
  if isempty(sidx)
      Trials(tr).Trigger = [];
  else
      if strncmpi(trigtype,'start',4)
          Trials(tr).Trigger = [Trials(tr).Saccades(sidx).start];
      elseif strncmpi(trigtype,'peak',4)
          Trials(tr).Trigger = [Trials(tr).Saccades(sidx).peakt];
      elseif strncmpi(trigtype,'end',3)
          Trials(tr).Trigger = [Trials(tr).Saccades(sidx).end];
      else
          Trials(tr).Trigger = [Trials(tr).Saccades(sidx).start];
      end
      Trials(tr).SacTrigger = Trials(tr).Trigger;
      sacamp = [sacamp [Saccades(sidx).size]];
      sacdir = [sacdir [Saccades(sidx).dir]];
      sacdur = [sacdur [Saccades(sidx).end] - [Saccades(sidx).start]];
      sacpeak = [sacpeak [Saccades(sidx).peakt] - ...
          [Saccades(sidx).start]];
      sacstarts = [sacstarts [Saccades(sidx).start]];
  end
      spikes = Trials(tr).Spikes;
end
  saccades.amp = sacamp;
  saccades.dir = sacdir;
  saccades.peak = sacpeak;
  saccades.dur = sacdur;
  saccades.start = sacstarts;
  xtimes = [min(times)-smooth*2 times max(times)+smooth*2];
[sdf, nsac, nspikes] = trigsdfa(Trials, smooth, xtimes, varargin{:});
sdf = sdf(2:end-1);

