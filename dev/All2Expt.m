function Expt = All2Expt(AllE, cell, varargin)
%Expt = All2Expt(AllExpt, cell, varargin)
%makes a single Expt struct from a multicell AllExpt struct
%
mu = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'all',2)
        mu= 2;
    elseif strncmpi(varargin{j},'mu',2)
        mu= 1;
    elseif strncmpi(varargin{j},'probe',2)
        findprobe =1;
    end
    j = j+1;
end

Expt = AllE.Expt;
cells = [AllE.Header.cellnumber];
probes = [AllE.Header.probe];

if mu == 0
    id = find(cells == cell);
elseif mu == 2
    id = cell;
elseif findprobe
    id = find(probes == cell);
    id = id(1); %is cell and MU, use cell
else
    id = find(probes == cell & cells == 0);
end

if isempty(id)
    return;
end
    
f = fields(AllE.Header);
for j = 1:length(f)
    Expt.Header.(f{j}) = AllE.Header(id).(f{j});
end

uset = ismember([Expt.Trials.id],AllE.Spikes{id}.trialid);
if length(uset) ~= length(AllE.Spikes{j}.trialid)
    fprintf('Trial Id length mismatch');
end

Expt.Trials = Expt.Trials(uset);
for j = 1:length(AllE.Spikes{id}.Spikes)
    Expt.Trials(j).Spikes = double(AllE.Spikes{id}.Spikes{j});
    Expt.Trials(j).OSpikes = double(AllE.Spikes{id}.OSpikes{j});
    Expt.Trials(j).Ocodes = double(AllE.Spikes{id}.Ocodes{j});
    count = sum(Expt.Trials(j).Spikes > 500 & Expt.Trials(j).Spikes < Expt.Trials(j).dur+500);
    Expt.Trials(j).count = count;
end
