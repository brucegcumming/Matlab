function  res = CheckIdx(name, varargin)
args = varargin;
if iscellstr(name)
    for j = 1:length(name)
        res{j} = CheckOneIdx(name{j}, args{:});
    end
elseif ischar(name)
    res = CheckOneIdx(name, args{:});
elseif iscell(name) && isfield(name{1},'Trials')
    res = CheckOneIdx(name, args{:});
end

function DATA = CheckOneIdx(name, varargin)

gaps = 4000;
plottrials = 0;
allgaps = [4000 8000 12000 14000 16000 18000 20000 22000 24000 26000 28000 30000 32000 64000];

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'gaps',4)
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            gaps = varargin{j};
        else
            gaps = allgaps;
        end
    elseif strncmpi(varargin{j},'plottrials',6)
        plottrials = 1;
    end
    j = j+1;
end

edur = [];
if ischar(name)
if exist(name,'file')
    load(name);
    Trials = Expt.Trials;
    for j = 1:length(Expts)
        if length(Expts(j).end)
        edur(j) = Expts(j).end-Expts(j).start;
        else
            edur(j) = 0;
        end
    end
else
    return;
end
elseif iscell(name)
    Expts = name;
    for j = 1:length(Expts)
        edur(j) = Expts{j}.Trials(end).End(end) - Expts{j}.Trials(1).Start(1);
    end
end

DATA.edur = edur;
DATA.gaps= gaps;
preperiod = 1000;
postperiod = 2000;
for j = 1:length(gaps)
    if iscell(name)
        DATA.blocks = CheckExptTime(name, preperiod, postperiod, gaps(j));
    else
        [DATA.blocks, edur] = CheckDataTime(Trials, preperiod, postperiod, gaps(j));
    end
    DATA.blockdur(j) = sum(diff(DATA.blocks))./10000;
    DATA.nblocks(j) = size(DATA.blocks,2);
    fprintf('%d blocks: %.1f sec\n',size(DATA.blocks,2),DATA.blockdur(j));
end
if length(gaps) > 1
    gaps = [gaps max(gaps)*2];
    if iscell(name)
        DATA.blockdur(j+1) = sum(edur)./10000;
        DATA.nexpts = length(edur);
        DATA.nblocks(j+1) = length(edur);
    elseif length(edur)
        DATA.blockdur(j+1) = sum(edur)./10000;
        DATA.nexpts = length(edur);
        DATA.nblocks(j+1) = length(edur);
    else
        DATA.blockdur(j+1) = (Trials.End(end)-Trials.Start(1))./10000;
        DATA.nexpts = 1;
        DATA.nblocks(j+1) = 1;
    end
    plot(gaps,DATA.blockdur./DATA.blockdur(1),'-o');
    hold on;
end
if plottrials && iscell(name)
    PlotExptBlocks(name, DATA);
elseif plottrials
    for j = 1:length(Trials.Start);
        if Trials.Result(j) > 0
            plot([Trials.Start(j) Trials.Start(j)],[0 1],'r-');
        else
            plot([Trials.Start(j) Trials.Start(j)],[0 1],'g-');
        end
        hold on;
    end
    for j = 1:length(Expts)
        if length(Expts(j).end)
            plot([Expts(j).start Expts(j).start Expts(j).end Expts(j).end Expts(j).start],[-0.2 1.2 1.2 -0.2 -0.2],'k');
        end
    end
    for j = 1:length(DATA.blocks)
        plot([DATA.blocks(1,j) DATA.blocks(2,j)],[-0.05 -0.05],'k','linewidth',2);
    end
    set(gca,'ylim', [-0.5 1.5]);
    title(name);
end


function blocks = CheckDataTime(Trials, preperiod, postperiod, gaplen)

nblk = 0;
id = find(Trials.Result ~= 0);
starts(1) = Trials.Start(id(1)) - preperiod;
for j = 2:length(id)
      dt = Trials.Start(id(j))-Trials.End(id(j-1));
      if dt > gaplen
          nblk = nblk+1;
          ends(nblk) = Trials.End(id(j-1))+postperiod;
          starts(nblk+1) = Trials.Start(id(j));
      end
  end
  ends(nblk+1) = Trials.End(end)+postperiod;
  blocks = [starts; ends];

function [blocks, gaps] = CheckExptTime(Expts, preperiod, postperiod, gaplen)
nblk = 0;

for j = 1:length(Expts)
    nblk = nblk+1;
    id = find([Expts{j}.Trials.RespDir] ~= 0);
    starts(nblk) = Expts{j}.Trials(1).Start(1) - preperiod;
    Trials = Expts{j}.Trials(id);
    for t = 2:length(Trials);
      gaps(t) = Trials(t).Start(1)-Trials(t-1).End(end);
      if gaps(t) > gaplen
          ends(nblk) = Trials(t-1).End(end)+postperiod;
          nblk = nblk+1;
          starts(nblk) = Trials(t).Start(1)-preperiod;
      end
  end
  ends(nblk) = Trials(end).End(end)+postperiod;
end
  blocks = [starts; ends];

  
  
function PlotExptBlocks(Expts, DATA)
for j = 1:length(Expts)
    id = find([Expts{j}.Trials.RespDir] ~= 0);
    Trials = Expts{j}.Trials(id);
    for t = 1:length(Trials);
        if Trials(t).RespDir ~= 0
            plot([Trials(t).Start Trials(t).Start],[0 1],'r-');
        else
            plot([Trials(t).Start Trials(t).Start],[0 1],'g-');
        end
    hold on;
    end
end
for j = 1:length(Expts)
    plot([Expts{j}.Header.Start Expts{j}.Header.Start Expts{j}.Header.End Expts{j}.Header.End Expts{j}.Header.Start],[-0.2 1.2 1.2 -0.2 -0.2],'k');
end
for j = 1:length(DATA.blocks)
    plot([DATA.blocks(1,j) DATA.blocks(2,j)],[-0.05 -0.05],'k','linewidth',2);
end
set(gca,'ylim', [-0.5 1.5]);
title(Expts{1}.Header.expname);
    
function  DATA = CheckExptsDur(Expts, varargin)
aps = 4000;
plottrials = 0;
allgaps = [1000 2000 4000 8000 12000 14000 16000 18000 20000 22000 24000 26000 28000 30000 32000 64000];

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'gaps',4)
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            gaps = varargin{j};
        else
            gaps = allgaps;
        end
    elseif strncmpi(varargin{j},'plottrials',6)
        plottrials = 1;
    end
    j = j+1;
end

edur = [];


DATA.gaps= gaps;
preperiod = 1000;
postperiod = 2000;

DATA.edur = edur;

