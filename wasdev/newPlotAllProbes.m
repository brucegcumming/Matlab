function result = PlotAllProbes(name, varargin)
% PlotAllProbes(name, varargin) reads in Expt files for each probe, and
% combines them

probes = 1:24;
cluster = 1;
plottype = 0;
delay = 0.05;
yvals = 1;
showblank = 0;
scales = ones(24,1);
timerange = [200 2000];
newscale = 0;
zvals = [];
isrc = [];
latecny = 500;
PLOTLFPDIFF=6;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'movie',4)
        plottype = 1;
    elseif strncmpi(varargin{j},'blank',4)
        showblank = 1;
    elseif strncmpi(varargin{j},'default',5)
        plottype = 0;
    elseif strncmpi(varargin{j},'LFPdiff',5)
        plottype = PLOTLFPDIFF;
    elseif strncmpi(varargin{j},'latency',4)
        j = j+1;
        latency = varargin{j};
    elseif strncmpi(varargin{j},'netspk',5)
        plottype = 2;
    elseif strncmpi(varargin{j},'probes',4)
        j = j+1;
        probes = varargin{j};
    elseif strncmpi(varargin{j},'pause',4)
        j = j+1;
        delay = varargin{j};
    elseif strncmpi(varargin{j},'scales',4)
        j = j+1;
        newscale = 1;
        scales(1:length(varargin{j})) = varargin{j};
    elseif strncmpi(varargin{j},'timerange',6)
        j = j+1;
        timerange = varargin{j};
    elseif strncmpi(varargin{j},'yval',4)
        j = j+1;
        yvals = varargin{j};
    elseif strncmpi(varargin{j},'zval',4)
        j = j+1;
        zvals = varargin{j};
    end
    j = j+1;
end



if isstruct(name)
    rc = name;
% if there is an 'lfp' matrix, uset that. Otherwise, if this is an rc file
% (output from PlotRevCOrAny, build the matrix
    if ~isfield(name,'lfp') & isfield(name,'sdfs') & isfield(name.sdfs,'lfp')
        [rc.lfp, rc.mlfp, rc.lfpblank] = CombineRCLFP(rc,scales);
        rc.lfptimes = rc.sdfs.lfptimes;
        isrc = 1;
    end
    rc.plot = 0;

    if isfield(rc,'tresps')
        if ndims(rc.tresps) == 3
            rc.spkresps = rc.tresps
        elseif ndims(rc.tresps) == 4
            rc.spkresps = squeeze(mean(rc.tresps(:,:,yvals,:),3));
            if isfield(rc,'blankresp')
                j = size(rc.spkresps,2)+1;
                rc.spkresps(:,j,:) = rc.blankresp;
            end
        end
        isrc = 1;
    end
    
        if isfield(rc,'lfp')
        if ndims(rc.lfp) == 3
            rc.lfpresp = rc.tresps
        elseif ndims(rc.lfp) == 4
            if isrc
                if plottype == PLOTLFPDIFF
                    yvals = [1 2];
                    zvals = [4 4];
                    rc.plot = PLOTLFPDIFF;
                end
                rc.lfpresp = squeeze(mean(rc.lfp(:,:,yvals,:),3));
                if length(rc.lfpblank) & showblank
                    j = size(rc.lfpresp,2)+1;
                    rc.lfpresp(:,j,:) = rc.lfpblank;
                end
                if ~isempty(zvals)
                    rc.lfprespa = squeeze(mean(rc.lfp(:,:,zvals,:),3));
                end
            else
            rc.lfpresp = squeeze(mean(rc.lfp(:,:,yvals,:),3));
            end
        end
        if newscale
            for j = 1:size(rc.lfpresp,3)
                rc.lfpresp(:,:,j) = rc.lfpresp(:,:,j) .* scales(j);
            end
        end
        end
    
        
    if plottype == 0
        if isrc
            plottype = 1;
        else
            plottype = 3; %image
        end
    end
    if plottype == 1
        if length(probes) > 1
            result = PlotMovie(rc,delay);
        else
            PlotTimeCourse(rc,probes);
        end
    elseif plottype == 2 %netspike plot for spikes
        tid = find(rc.times > timerange(1) & rc.times < timerange(2));
        meanresp = repmat(mean(rc.tresps(tid,:,1,:),2),[1 size(rc.tresps,2) 1 1]);
        netspk = rc.tresps(tid,:,:,:) - meanresp;
        netspk = squeeze(sum(netspk,1));
        imagesc(netspk);
        colorbar;
    elseif plottype == 3
        if size(rc.resps,2) > 1
            if isfield(rc,'lfp')
                tid = find(rc.lfptimes > latency)
                subplot(2,1,2);
                imagesc(rc.x(:,1), probes, squeeze(var(rc.lfpresp(tid,:,:))));
                subplot(2,1,1);
            else
                subplot(1,1,1);
            end
            
            imagesc(rc.x(:,1), probes, rc.resps);
        else
            if isfield(rc,'lfp')
                subplot(1,1,1);
                imagesc(squeeze(rc.lfp.lfp));
            else
            plot(rc.means);
            end
        end
    elseif plottype == PLOTLFPDIFF;
    end
    return;
elseif isdir(name)
    d = dir(name);
    nex = 1;
    for j = 1:length(d);
        if strfind(d(j).name,'.p1c1')
            fprintf('%s\n',d(j).name);
                        result.names{nex} = strrep(d(j).name,'.p1c1','.c1');
            result.exps{nex} = PlotAllProbes([name '/' result.names{nex}]);
            nex = nex+1;
        end
    end
    return;
end

if isempty(isrc)
if strfind(name,'RC')
    isrc = 1;
else
    isrc = 0;
end
end

cs = ['c' num2str(cluster) '.'];
ts = 34;
xvals = [];
rcargs = {};
if strfind(name,'DCORRC') 
    rcargs = {rcargs{:} 'psych' 'yvals' [0 0.06]};
end
for j = probes;
    ename = strrep(name,['.' cs],['.p' num2str(j) cs]);
    load(ename);
    if isrc
        res = PlotRevCorAny(Expt,'sdfw',166,'box',rcargs{:});
        ts = res.bestdelay;
        means(j) = mean(res.y(:));
        resps(j,:) = res.y(:,1,ts)./means(j);
        id = find(res.times > 166/2);
        ta = id(2);
        xvals = res.x;
        for iy = 1:size(res.sdfs.s,2)
            for ix = 1:size(res.sdfs.s,1)
            result.tresps(:,ix,iy,j) = res.sdfs.s{ix,iy}(ta:end)./means(j);
            result.vtimes(j,:) = res.delays;
            result.besttimes(j,:) = ts;
            end
        end
        if length(res.sdfs.extras) & strmatch(res.sdfs.extras{1}.label,'Blank')
            result.blankresp(:,j) = res.sdfs.extras{1}.sdf(ta:end)./means(j);
        end
        result.times = res.times(ta:end);
    else
        Expt = FillTrials(Expt,Expt.Stimvals.et); %in case
        res = PlotExpt(Expt);
        xvals = res.x;
        if isfield(res,'means')
        means(j) = mean(res.means(:));
        resps(j,:) = res.means(:)./means(j);
        end
    end
end

result.resps = resps;
result.means = means;
result.x = xvals;

ename = strrep(name,['.' cs],'.lfp.');
if exist(ename,'file')
    load(ename);
    if isrc
    rc = PlotRevCorAny(LFP,'lfp', rcargs{:})
    [result.lfp, mresp, result.lfpblank] = CombineRCLFP(rc, scales);
    result.mlfp = mresp;
    result.lfptimes = rc.sdfs.lfptimes;
    result.lfpwr = mean(abs(cat(3,LFP.Trials.FTlfp)),3);
    for j = find(scales ~= 1)
        result.lfpwr(:,j) = result.lfpwr(:,j) .* scales(j);
    end
    len = size(LFP.Trials(1).LFP,1)
    result.lfpfrq = (0:len-1)/(len * LFP.Header.LFPsamplerate .* 10);
    else
        LFP = FillTrials(LFP,Expt.Stimvals.et);
        lfp = PlotLFP(LFP,'lfp');
        result.lfp = lfp.lfp;
        result.lfptimes = lfp.lfptimes;
        result.lfpwr = lfp.lfppower;
        result.lfpfrq = lfp.ftfrq;
    end
end

function [sdfs, mresp, blank] = CombineRCLFP(rc, scales)

for j = 1:size(rc.sdfs.lfp,1)
    for k = 1:size(rc.sdfs.lfp,2)
    for ch = 1:size(rc.sdfs.lfp{j,k},2)
        sdfs(:,j,k,ch) = rc.sdfs.lfp{j,k}(:,ch) .* scales(ch);
    end
    end
end

mresp = squeeze(mean(sdfs,2));

if length(rc.sdfs.extras) & strmatch(rc.sdfs.extras{1}.label,'Blank');
        blank = rc.sdfs.extras{1}.lfp-mresp;
        for j = 1:size(blank,2)
            blank(:,j) = blank(:,j) .* scales(j);
        end
else
    blank = [];
end

for j = 1:size(rc.sdfs.lfp,1)
    sdfs(:,j,:,:) = squeeze(sdfs(:,j,:,:)) - mresp;
end


function NextFrame(a,b,step)

D = get(gcf,'UserData');
t = D.inow;
if step == 0
    t = t+1;
    if t > length(D.mtimes)
        D.mtimes = [D.mtimes D.mtimes(end)+mean(diff(D.mtimes))];
        t = length(D.mtimes);
    end
elseif step < 0
    t = t+step;
    if t < 1
        t = 1;
    end
end
PlotTimeSlice(D,D.mtimes(t));   
D.inow = t;
set(gcf,'UserData',D);

function PlayFrames(a,b,step)

D = get(gcf,'UserData');
t = D.inow;
set(findobj(gcf,'Tag','PlayStop'),'value',0);
if step > 0
    t = [D.inow:step:length(D.mtimes)];
elseif step < 0
    t = [D.inow:step:1];
end
for it = t
    stop = get(findobj(gcf,'Tag','PlayStop'),'value');
    if stop
        break;
    end
    PlotTimeSlice(D,D.mtimes(it));
    drawnow;
end
D.inow = it;
set(gcf,'UserData',D);

function res = PlotMovie(res, delay)

sustep = 1000;
lfpstep = 10000;
    

if isfield(res,'spkresps')
res.surange(1) = min(res.spkresps(:));
res.surange(2) = max(res.spkresps(:));
end
if isfield(res,'lfpresp')
    
res.lfprange(1) = min(res.lfpresp(:));
res.lfprange(2) = max(res.lfpresp(:));
end
res.mtimes = 2:10:1500;
PlotTimeSlice(res, res.mtimes(1));


fb = findobj(gcf,'Tag','NextFrame');
if isempty(fb)
    bp = [10 10 40 20];
    uicontrol(gcf, 'style','pushbutton','string','+','Position',bp,...
        'Tag','NextFrame','Callback',{@NextFrame, 0});
end
fb = findobj(gcf,'Tag','LastFrame');
if isempty(fb)
    bp = [50 10 40 20];
    uicontrol(gcf, 'style','pushbutton','string','-','Position',bp,...
        'Tag','LastFrame','Callback',{@NextFrame, -1});
end
fb = findobj(gcf,'Tag','PlayFwd');
if isempty(fb)
    bp = [90 10 40 20];
    uicontrol(gcf, 'style','pushbutton','string','>>','Position',bp,...
        'Tag','PlayFwd','Callback',{@PlayFrames,  1});
end
fb = findobj(gcf,'Tag','PlayBack');
if isempty(fb)
    bp = [130 10 40 20];
    uicontrol(gcf, 'style','pushbutton','string','<<','Position',bp,...
        'Tag','PlayBack','Callback',{@PlayFrames, -1});
end
fb = findobj(gcf,'Tag','PlayStop');
if isempty(fb)
    bp = [170 10 60 20];
    uicontrol(gcf, 'style','check','string','Stop','Position',bp,...
        'Tag','PlayStop');
end
set(gca,'UserData',res);
for ti = 1:length(res.mtimes)
    stop = get(findobj(gcf,'Tag','PlayStop'),'value');
    if stop
        ti = NaN;
    else
    PlotTimeSlice(res,res.mtimes(ti));
    res.inow = ti;
    set(gcf,'UserData',res);
    drawnow;
    if delay > 1
        pause
    else
       pause(delay);
    end
    end
end

function PlotTimeSlice(res,t);
    
if isfield(res,'lfprespa')
    subplot(2,1,1);
    [a,b] = min(abs(t-res.lfptimes));
    if res.plot
    imagesc(squeeze(res.lfprespa(b,:,:))' - squeeze(res.lfpresp(b,:,:))');
    else
    imagesc(squeeze(res.lfprespa(b,:,:))');
    end
    caxis(res.lfprange);
    title(sprintf('%.1f ms',res.lfptimes(b)./10));
elseif isfield(res,'spkresps')
    subplot(2,1,1);
    [a,b] = min(abs(t-mean(res.times,1)));
    imagesc(squeeze(res.spkresps(b,:,:))');
    caxis(res.surange);
    title(sprintf('%.1f ms',res.times(1,b)./10));
end

if isfield(res,'lfpresp')
    subplot(2,1,2);
    [a,b] = min(abs(t-res.lfptimes));
    if res.plot & isfield(res,'lfprespa')
        imagesc(squeeze(res.lfprespa(b,:,:)+res.lfpresp(b,:,:))./2);
    else
        imagesc(squeeze(res.lfpresp(b,:,:))');
    end
    caxis(res.lfprange);
    title(sprintf('%.1f ms',res.lfptimes(b)./10));
end

function PlotTimeCourse(rc, probe)

if isfield(rc,'lfp')
    subplot(2,1,2);
    imagesc(rc.x(:,1),rc.lfptimes./10,rc.lfpresp(:,:,probe));
    subplot(2,1,1);
else
    subplot(1,1,1);
end
imagesc(rc.spkresps(:,:,probe));

