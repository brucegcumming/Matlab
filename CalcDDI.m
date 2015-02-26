function ddi = CalcDDI(Expt, varargin)

% ddi = CalcDDI(Expt, duration)
% Calculates a DDI over a set of disparities. Using the spikes in
% the first 'duration' timestamp units.

% returns a strucutre ddi
% ddi.ddi           the ddi
% ddi.max           Max value of mean(sqrt count).
% ddi.min           Min value of mean(sqrt count).
% ddi.sqerr         RSS of (sqrt(count))s 
% ddi.err           RSS of counts  (without sqrt transform).

latency = 500;
duration = min([Expt.Trials.End] - [Expt.Trials.Start]);

nmin = 2;
j = 1;
type = Expt.Stimvals.et;
if isnumeric(type)
    if type == 237
        if isfield(Expt.Trials,'dO')
            type = 'dO';
        elseif isfield(Expt.Trials,'dx')
            type = 'dx';
        end
    elseif type == 5
        type = 'dx';
    end
end

while j < length(varargin)
    if strncmpi(varargin{j},'Nmin',4)
        j = j+1;
        nmin = varargin{j};
    elseif strncmpi(varargin{j},'Duration',4)
        j = j+1;
        duration = varargin{j};
    elseif strncmpi(varargin{j},'type',4)
        j = j+1;
        type = varargin{j};
    end
    j = j+1;
end

scale = 10000/duration;

if strcmp(type,'dO') & ~isfield(Expt.Trials,'dy') & ~isfield(Expt.Trials,'dO') & isfield(Expt.Trials,'dx')
    type = 'dx';
end
for j = 1:length(Expt.Trials)
    dx = Expt.Trials(j).(type);
    dx = round(dx.*500)./500;
    Expt.Trials(j).(type) = dx;
end
dxs = unique([Expt.Trials.(type)]);

if(~isfield(Expt.Trials,'ce'))
  [Expt.Trials.ce] = deal(1);
end

if(~isfield(Expt.Trials,'me'))
  if(isempty(Expt.Stimvals.me))
    [Expt.Trials.me] = deal(0);
  else
    [Expt.Trials.me] = deal(Expt.Stimvals.me);
  end
end

if(~isfield(Expt.Trials,'st'))
  if(isempty(Expt.Stimvals.st))
    [Expt.Trials.st] = deal(0);
  else
    [Expt.Trials.st] = deal(Expt.Stimvals.st);
  end
end

sqmeans = [];
sqvar = [];
ns = [];
vars = [];
Expt = FillTrials(Expt,'ce');
Expt = FillTrials(Expt,'st');
ceval = max([Expt.Trials.ce]);
ndx = 1;
if strcmp(Expt.Stimvals.e2,'dp')
    id = find([Expt.Trials.dp] == 0);
    Expt.Trials = Expt.Trials(id);
end
for dx = dxs
    if strcmp(type,'dx') || strcmp(type,'dO')
        idx =find([Expt.Trials.(type)] == dx & [Expt.Trials.me] == 0 & [Expt.Trials.ce] ...
            == ceval & [Expt.Trials.st] > 0);
    else
        idx =find([Expt.Trials.(type)] == dx & [Expt.Trials.st] > 0);
    end
  k = 1;
  if length(idx)
  for j = idx
    counts(j) = length(find([Expt.Trials(j).Spikes] > latency & ...
			   [Expt.Trials(j).Spikes] < latency + duration));
    ddi.rates{ndx}(k) = 10000 * counts(j)/duration;
    k = k+1;
  end
  sqrates(idx) = deal(sqrt(counts(idx) .* scale));
  rates(idx) = deal(counts(idx) .* scale);
  sqmeans = [sqmeans mean(sqrates(idx))];
  sqvar = [sqvar var(sqrates(idx))];
  vars = [vars var(rates(idx))];
  ns = [ns length(idx)];
  ndx = ndx+1;
  end
end

idx = find(ns > nmin);
ddi.max = max(sqmeans(idx));
ddi.min = min(sqmeans(idx));
rmvar = sum(sqvar .* (ns-1))/sum(ns-1);
if isempty(idx)
    ddi.sqerr = NaN;
    ddi.ddi = NaN;
    ddi.err = 0;
else
    ddi.sqerr = sqrt(rmvar);
    ddi.err = sqrt(rmvar);
    ddi.ddi = (ddi.max - ddi.min)/(ddi.max - ddi.min + 2 * ddi.sqerr);
end
