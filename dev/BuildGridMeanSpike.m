function result = BuildGridMeanSpike(dirname, varargin)
% result = BuildGridMeanSpike(dirname,expt, probe)
% given a source probe, build mean waveform on all probes
%
%BuildGridMeanSpike('X:/Utah/jbe/G040/',24,22)
%BuildGridMeanSpike('X:/Utah/jbe/G040/',24,22,'spts',-15000:1000,'distance',4)
j = 1;
cl = 1;
pw = 1; %distance to include
shuffletimes = 0;
spts = [];
filterA = [];
if ischar(dirname) && isdir(dirname)
    e = varargin{1};
    probes = varargin{2};
    j = 3;
end
while j <= length(varargin)
    if iscell(varargin{j}) && isfield(varargin{j}{1},'Trials')
        Expts = varargin{j};
    elseif strncmpi(varargin{j},'cutoff',6)
        j = j+1;
        fz = varargin{j};
        f = fz./15000;
        [filterA, filterB]  = butter(8,f);
    elseif strncmpi(varargin{j},'cluster',2)
        j = j+1;
        cl = varargin{j};
    elseif strncmpi(varargin{j},'distance',4)
        j = j+1;
        pw = varargin{j};
    elseif strncmpi(varargin{j},'shuffle',4)
        shuffletimes = 1;        
    elseif strncmpi(varargin{j},'spts',4)
        j = j+1;
        spts = varargin{j};
    end
    j = j+1;
end


if isstruct(dirname) %plot a previous result
    result = PlotGridMean(dirname,varargin{:});
    return;
end

cfile = [dirname '/Expt' num2str(e) 'ClusterTimes.mat'];
load(cfile);
if strcmp(probes,'cells')
    cellfile = [dirname '/CellList.mat'];
    if exist(cellfile,'file')
        load(cellfile);
    end
    probes = find(sum(CellList(e,:,:),3) > 0);
end

A = GetArrayConfig(dirname);
V = {};
for np = 1:length(probes)
    p = probes(np);
    x = A.X(p);
    y = A.Y(p);
    distances = abs(A.X-x + i *(A.Y-y));
    pid = find(abs(A.X-x) <= pw & abs(A.Y-y) <= pw);

    for  j = 1:length(pid)
        if pid(j) > length(V) || isempty(V{pid(j)})
            pfile = [dirname sprintf('/Expt%d.p%dFullV.mat',e,pid(j))];
            V{pid(j)} = LoadFullV(pfile,'nohighpass');
            if ~isempty(filterA)
                V{pid(j)}.V = filtfilt(filterA, filterB, V{pid(j)}.V);
            end
        end
    end
    for j = 1:length(V{p}.blkstart)
        bt{j} = V{p}.blkstart(j) + (1:V{p}.blklen(j)).* V{p}.samper;
    end
    t = cat(2,bt{:});
    if cl > 1
        [st,tid] = ismember(Clusters{p}.next{cl-1}.times,t);
    else
        if shuffletimes
            [st,tid] = ismember(Clusters{p}.times,t);
            ts = ShuffleTimes(Clusters{p}.times, Expts{e});
            [st,stid] = ismember(round(ts.*30000),round(t.*30000));
            stid = stid(stid>0);
            tid = tid(stid >0);
        else
            [st,tid] = ismember(Clusters{p}.times,t);
        end
    end
    if isempty(spts)
        spts = Clusters{p}.spts;
    end
    npts = length(V{p}.V);
    allid = repmat(tid,length(spts),1)' + repmat(spts',1,length(tid))';
    allid(allid<=0) = 1;
    allid(allid>=npts) = npts;
    Amps = [];
    MeanV = [];
    for  j = 1:length(pid)
        MeanV(j,:) = mean(V{pid(j)}.V(allid));
        result(np).X(j) = A.X(pid(j));
        result(np).Y(j) = A.Y(pid(j));
        result(np).distance(j) = distances(pid(j));
        result(np).mean = mean(V{pid(j)}.V);
    end


    pi = find(pid ==p);
    vsq = sum(MeanV(pi,:).^2);
    for  j = 1:length(pid)
        Amps(j) = sqrt(MeanV(j,:) * MeanV(pi,:)'/vsq);
    end
    result(np).MeanV = MeanV;
    result(np).pid = pid;
    result(np).Amps = Amps;
    plot(MeanV');
    refline(0);
    h = refline(0, result(np).mean);
    set(h,'color','r');
    if shuffletimes
        allid = repmat(stid,length(spts),1)' + repmat(spts',1,length(tid))';
        allid(allid<=0) = 1;
        allid(allid>=npts) = npts;
        for  j = 1:length(pid)
            result(np).ShuffleV(j,:) = mean(V{pid(j)}.V(allid));
        end
        if length(pid) ==1
            hold on;
            plot(result(np).ShuffleV,'r');
        end
    end
end

function PlotGridMovie(result, varargin);
nx = length(unique(result.X));
ny = length(unique(result.Y));
xp = unique(result.X);
yp = unique(result.Y);
first = 1;
last = size(result.MeanV,2);
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'range',4)
        j = j+1;
        first = varargin{j}(1);
        last = varargin{j}(2);
    end
    j = j+1;
end
X = zeros(nx,ny);
id = sub2ind(size(X),result.X,result.Y);
sm=50;
DATA.id = id;
DATA.smooth = sm;
subplot(1,1,1);
DATA.clim = [-800 400];
DATA.MeanV = result.MeanV;
DATA.image = X;
for j = first:sm:last
    DATA.tpoint = j;
    PlotImage(DATA);
end
DATA.tpoint = j;
set(gcf,'UserData',DATA);
set(gcf,'KeyPressFcn',@KeyPress);

function PlotImage(DATA)
    smid = DATA.tpoint:DATA.tpoint+DATA.smooth-1;
    if smid(end) <= size(DATA.MeanV,2)
        DATA.image(DATA.id) = mean(DATA.MeanV(:,smid),2);
        imagesc(DATA.image',DATA.clim);
        drawnow;
        title(sprintf('%d',DATA.tpoint))
    end
    set(gcf,'UserData',DATA);
        
function KeyPress(src, key)

DATA = get(src,'UserData');
if strcmp(key.Key,'rightarrow')
    DATA.tpoint = DATA.tpoint+DATA.smooth;
    PlotImage(DATA);
elseif strcmp(key.Key,'leftarrow')
    DATA.tpoint = DATA.tpoint-DATA.smooth;
    PlotImage(DATA);
end

function result = PlotGridMean(result, varargin)

nx = length(unique(result(1).X));
ny = length(unique(result(1).Y));
xp = unique(result(1).X);
yp = unique(result(1).Y);
pw = 100;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'movie',5);
        PlotGridMovie(result, varargin{:});
        return;
    elseif strncmpi(varargin{j},'distance',4)
        j = j+1;
        pw = varargin{j};
    end
    j = j+1;
end

if length(result) > 1
    [nx,ny] = Nsubplots(length(result)+1);
    for j = 1:length(result)
        mysubplot(ny, nx, j);
        hold off;
        id = find(result(j).distance <= pw);
        MMean(j,:) = mean(result(j).MeanV(id,:),1);
        plot(MMean(j,:));
        refline(0);
        if isfield(result,'ShuffleV')
            hold on;
            SMean(j,:) = mean(result(j).ShuffleV(id,:),1);
            plot(mean(result(j).ShuffleV(id,:),1),'r');
        end
        set(gca,'xtick',[],'ytick',[],'ylim',[-800 400],'xlim',[0 size(result(j).MeanV,2)]);
    end
    mysubplot(ny, nx, j+1);
    hold off;
    plot(mean(MMean,1));
    hold on;
    plot(mean(SMean,1),'r');
    return;
end
freeid = [];
for j = 1:nx;
    for k = 1:ny
        id = find(result.X == xp(j) & result.Y == yp(k));
        if length(id) == 1
            mysubplot(ny, nx, j + (k-1)*nx);
            plot(result.MeanV(id,:));
            if isfield(result,'ShuffleV')
                hold on;
                plot(result.ShuffleV(id,:),'r');
            end
            set(gca,'xtick',[],'ytick',[],'ylim',[-800 400],'xlim',[0 size(result.MeanV,2)]);
            text(0,400,sprintf('%d %d',result.X(id),result.Y(id)));
            refline(0);
            amps(j,k) = abs(result.Amps(id));
            if result.distance(id) == 0
                hold on;
                plot(result.MeanV(id,:)./10,'r');
            end
        else
            freeid = [freeid j + (k-1)*nx];
            amps(j,k) = NaN;
        end
    end
end
id = find(result.distance <= pw);
mysubplot(ny, nx, freeid(1));
imagesc(amps);
if length(freeid) > 1
mysubplot(ny, nx, freeid(2));
hold off;
plot(mean(result.MeanV(id,:),1));
refline(0);
if isfield(result,'ShuffleV')
    hold on;
    plot(mean(result.ShuffleV(id,:),1),'r');
end
set(gca,'xtick',[],'ytick',[],'ylim',[-800 400],'xlim',[0 size(result.MeanV,2)]);
end
