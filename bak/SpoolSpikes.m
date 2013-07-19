function DATA = SpoolSpikes(DATA, varargin)
%
% GUI for playing back spikes from existing expts.
% Called from PlotExptSpikes
[DATA.spkvarnames, DATA.spkvarorder] = GetSpikeVals(DATA,NaN,NaN,NaN,[]);
DATA = PlaySpikes(DATA, DATA.Expts{DATA.currentexpt},varargin{:});

function DATA = PlaySpikes(DATA, Expt, varargin)

global mousept;
plotdprime = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'plotdprime',6)
        plotdprime = 1;
    end
    j = j+1;
end
expid = 0;
mode = 2;
cw = 10;
ch = 10;
rh = ch+10;
[xyfig, isnew] = GetFigure(DATA.tag.clusterxy);
DATA.xyfig = xyfig;
if isnew
    bp = [5 5 40 20];
    cp = bp;
    uicontrol(xyfig,'style','pop','string','1|2|3|4|5|6|7','Position',bp,'Tag','Clusterid',...
        'Callback',@Update);
    cp(3) = cw * 5;
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'Style', 'checkbox', 'Callback', @DensityPlot, ...
'String', 'Dens','Tag','Density', 'Position', cp);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', {@NextList, 'setfirst'}, ...
'String', 'Set+Next', 'Position', cp,'tag','Set+Next');
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', @NextList, ...
'String', 'Next', 'Position', cp);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', {@NextList, 'clearfirst'}, ...
'String', 'Clr+Next', 'Position', cp);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'style','pushbutton','string','spool','Position',cp,'Tag','SpoolSpikes',...
        'Callback', @SpoolExptSpikes);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'style','pushbutton','string','Optim','Position',cp,'Tag','OptimizeDprime',...
        'Callback', @OptimizeDprimeHit);
    cp(2) = cp(2)+rh;   
    if length(DATA.probelist) > 1
        cp(3) = cw*1.5;
    uicontrol(xyfig,'style','pushbutton','string','+','Position',cp, 'Callback',@AddOneCellToList);

    cp(1) = cp(1) + cp(3);
    cp(3) = cw*3;
    uicontrol(xyfig,'style','pop','string','1|2|3|4|5|6|7|8|9|10','Position',cp,'Tag','AddOneCellToList',...
        'Callback',@AddOneCellToList);
    end
    bp(1) = bp(1) + bp(3) + 10;
    uicontrol(xyfig,'style','pop','string','Not Set|Nothing|MU-|MU+|Poor|OK|Good|V Good|Excellent','Position',bp,'Tag','ClusterQuality',...
        'Callback',@Update);

    bp(1) = bp(1) + bp(3) + 10;
    uicontrol(xyfig,'style','pushbutton','string','Set','Position',bp,'Callback', @SetExptClusters);
    bp(1) = bp(1) + bp(3) + 10;
    bp(3)=cw;
    uicontrol(xyfig,'style','CheckBox','Position',bp,'Tag','ClusterIsSet');
    bp(1) = bp(1) + bp(3) + 10;
    bp(3)=cw*3;
    uicontrol(xyfig,'style','pushbutton','string','Del','Position',bp,'Callback', @DelClusterButton);
    bp(1) = bp(1) + bp(3) + 10;
    uicontrol(xyfig,'style','pushbutton','string','Clr','Position',bp,'Callback', @ClrSpkWin);
    bp(1) = bp(1) + bp(3) + 10;
    bp(3)=cw*5;
    uicontrol(xyfig,'style','text','string','Max: X','Position',bp);
    bp(1) = bp(1) + bp(3);
    uicontrol(xyfig,'style','edit','string',sprintf('%.2f',DATA.plot.clusterXrange(2)),'Position',bp,...
    'Callback', @RescaleClusterPlot,'Tag','ClusterXmax');

    bp(1) = bp(1) + bp(3) + 10;
    bp(3) = cw*1;
    uicontrol(xyfig,'style','text','string','Y','Position',bp);
    bp(1) = bp(1) + bp(3);
    bp(3) = cw*4;
    uicontrol(xyfig,'style','edit','string',sprintf('%.2f',DATA.plot.clusterYrange(2)),'Position',bp,...
    'Callback', @RescaleClusterPlot,'Tag','ClusterYmax');
    bp(1) = bp(1) + bp(3) + 10;
    bp(3) = cw*1;
    uicontrol(xyfig,'style','text','string','Z','Position',bp);
    bp(1) = bp(1) + bp(3);
    bp(3) = cw*4;
    uicontrol(xyfig,'style','edit','string',sprintf('%.2f',DATA.plot.clusterZrange(2)),'Position',bp,...
    'Callback', @RescaleClusterPlot,'Tag','ClusterZmax');
    bp(1) = bp(1) + bp(3);
    bp(3) = cw*5;
    uicontrol(xyfig,'style','CheckBox','string','auto','Position',bp,'Tag','AutoScale',...
        'value',(DATA.plot.autoscale > 0),'Callback',@Update);
    bp(1) = bp(1) + bp(3) + 10;
    hm = uimenu(xyfig,'Label','Plot3D','Tag','3DplotMenu');
     uimenu(hm,'Label','None','Callback',{@SetXYCluster, 3, 2, 0});
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@SetXYCluster, 3, 2,k});
    end
    hm = uimenu(xyfig,'Label','X','Tag','XClusterMenu');
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@SetXYCluster, 1, 2,k});
    end
    hm = uimenu(xyfig,'Label','Y','Tag','YClusterMenu');
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@SetXYCluster, 2, 2, k});
    end
    hm = uimenu(xyfig,'Label','Plot','Tag','PlotMenu');
    uimenu(hm,'Label','Track dprime','Callback',{@PlotDprime,1});
    uimenu(hm,'Label','Smoothed dprime','Callback',{@PlotDprime,2});
    uimenu(hm,'Label','Count','Callback',{@PlotDprime,3});
end
DATA.toplevel = xyfig;
set(xyfig,'UserData',DATA);
hold off;
%if isempty(findobj('Tag',DATA.tag.spikev
 %   sfig = figure('Renderer','painters','Tag','DATA.tag.spikev');
[sfig, isnew] = GetFigure(DATA.tag.spikev);
if isnew 

    x = 10;
    c = 10;
    bp = [10 10 40 20];
    uicontrol(sfig,'style','pushbutton','string','>>','Position',bp,'Tag','NextTrial',...
        'Callback', @PlayNextTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pushbutton','string','<<','Position',bp,'Tag','LastTrial',...
        'Callback', @PlayLastTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pushbutton','string','X','Position',bp,'Tag','CutTrial',...
        'Callback', @CutTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pop','string','1|2|3','Position',bp,'Tag','ChooseTrial',...
        'Callback', @SelectTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pushbutton','string','+','Position',bp,'Tag','NextSpike',...
        'Callback', @PlayNextSpike);
    bp(1) = bp(1)+bp(3);
    bp(3)=c*6;
    uicontrol(sfig,'style','pushbutton','string','spool','Position',bp,'Tag','SpoolSpikes',...
        'Callback', @SpoolExptSpikes);
    bp(1) = bp(1) + bp(3);
    
    bp(3)=c*9;
    uicontrol(gcf,'Style', 'checkbox',...
    'String', 'dVdt', 'Tag', 'dVdt', 'Position', bp,'value',DATA.plot.dvdt,...
    'Callback',@Update);
    bp(1) = bp(1) + bp(3);
    
    bp(3)=c*9;
    uicontrol(gcf,'Style', 'checkbox',...
    'String', 'Stop', 'Tag', 'StopSpool', 'Position', bp,'value',0,...
    'Callback',@Update);
    tmpdat.parentfigtag = DATA.tag.clusterxy;
    set(sfig,'UserData',tmpdat);
    set(sfig, 'WindowScrollWheelFcn',@ScrollTrial);

    DATA.hline = 0;    
end

    it = findobj(sfig,'Tag','ChooseTrial');
    set(it,'string',sprintf('%d|',[Expt.Trials.Trial]),'value',1);
    DATA.svfig = sfig;
ax = findobj(DATA.svfig,'type','ax');

if DATA.state.uselfp
    DATA.state.lfig = GetFigure('LFP');
    hold off;
    Expt = LoadSpike2LFP(Expt);
end
colors = mycolors;
Spks = DATA.AllData.Spikes;
tstart = now;
if isfield(Spks,'values')
nsamples = size(Spks.values,2);
else
    nsamples = 0;
end
if size(Spks.codes,2) > 1
nclusters = 1+max(Spks.codes(:,2));
else
nclusters = 1;
end
set(0,'CurrentFigure',xyfig);
if DATA.plot.autoscale == 0
    set(gca,'Xlim',DATA.plot.clusterXrange,'xlimmode','manual');
    set(gca,'Ylim',DATA.plot.clusterYrange,'ylimmode','manual');
end
set(0,'CurrentFigure',sfig);
xprobes = [];
probes = DATA.probe;
x = 32; %should calculate real size..
for k = 0:length(xprobes)
    ids{k+1} = FindSpikes(DATA, Expt.Header.trange,probes(k+1),[]);
    for j = 1:nclusters+1
        DATA.svh(j+k*(nclusters+1)) = plot([1:nsamples],[1:nsamples] * 0,'color',DATA.spkcolor{j});
        hold on;
    end
    if isfield(DATA,'AllSpikes')
        %nclusters needs to be max # used across all probes
        if length(DATA.plot.voffsets) > k & DATA.plot.prettyfigs
        h = text(25,DATA.plot.voffsets(k+1)+1,sprintf('Probe %d',probes(k+1)));
        set(h,'FontWeight','bold');
        else
        text(x,(k-(nprobes-1)/2) * DATA.plot.SpikeVsep,num2str(probes(k+1)));
        end
        ncut = sum(DATA.AllSpikes{probes(k+1)}.codes(:,2) > 0);
        if ncut
            id = find(ismember(DATA.AllSpikes{probes(k+1)}.codes(ids{k+1},2),spikelist));
%            ids{k+1} = ids{k+1}(id);
        end
        if ~isfield(DATA.AllSpikes{probes(k+1)},'dVdt')
            DATA.AllSpikes{probes(k+1)} = CleanSpikes(DATA.AllSpikes{probes(k+1)});
        end
    end
end
set(xyfig,'UserData',DATA);



reclassify = 0;
expspks = ids{1};
%DATA = SetSpkCodes(DATA, DATA.AllSpikes{probes(j)}.spklist, probes(j),0);
 
%FindSpikes might list some spikes just past the list recorded in
%Expts.gui.spkrange.
if isfield(DATA,'Spikes')
if max(ids{1}) > length(DATA.Spikes.cx)
    ispk = length(DATA.Spikes.cx):max(ids{1});
    DATA = CalcClusterVars(DATA,  ispk);
end
end
if isfield(DATA.AllData.Spikes,'codes')
id = find(DATA.AllData.Spikes.codes(ids{1},2) > nclusters);
if length(id)
    DATA.AllData.Spikes.codes(ids{1}(id),2) = 0;
end
end

DATA = CheckForPCA(DATA, expspks, 0);
        xprobes = [];
        
DATA.vstep = DATA.plot.SpikeMaxV * 0.8;
    if DATA.plot.SpikeVsep > 0
        DATA.vstep = DATA.plot.SpikeVsep;
    end
if DATA.plot.dvdt == 2
    set(gca,'Xlim',[-5 5],'Ylim',[-5 5]);
elseif length(xprobes)
    x = (length(xprobes)-1) .* DATA.plot.SpikeMaxV/2;
    if length(xprobes) > 3
        x = (length(xprobes)-1) .* DATA.plot.SpikeVsep;
    end
    x = (length(xprobes)) .* DATA.plot.SpikeVsep;
    if length(xprobes) 
    set(gca,'Xlim',[1 nsamples+1],'Ylim',[-(DATA.plot.SpikeMaxV+x/2) DATA.plot.SpikeMaxV+x/2]);
    end
elseif nsamples > 1        
    set(gca,'Xlim',[1 nsamples],'Ylim',[-DATA.plot.SpikeMaxV DATA.plot.SpikeMaxV]);
end
DATA.nclusters = nclusters;
nt = 1;
allspks = [];
firstspk = 0;
hline = 0;
if ~isfield(DATA,hline)
    DATA.hline = 0;
end

DATA.playingspk = 1;
DATA.nclusters = nclusters;


Aargs = {}; timemode = 0;
trialdur = Expt.Trials(1).End(end) - Expt.Trials(1).Start(1);
timemode = DATA.plot.showwave;

if isfield(DATA,'AllSpikes') & length(probes) > 1 & DATA.plot.syncoverlay == 0
    Aargs = {Aargs{:} 'timemode'};
    timemode = 1;
end



if timemode
        if  ~isfield(DATA,'timefig')
        DATA.timefig = GetFigure('SpikeTime');
    end
    set(0,'CurrentFigure',DATA.timefig);
    set(gca,'Xlim',[0 trialdur/10]);
    for k = 0:length(xprobes)
        for j = 1:nclusters+1
            DATA.tvh(j+k*(nclusters+1)) = plot([1:46],[1:46] * 0,'color',DATA.spkcolor{j});
            hold on;
        end
    end
end

probes = DATA.probe;
if length(probes) > 1 && isfield(DATA,'AllSpikes')
    if DATA.plot.synccluster > 0
        DATA.Spikes.cx(1:end) = 0;
        DATA.Spikes.cy(1:end) = 0;
    end
    step = DATA.plot.SpikeMaxV;
%syncsign = -1 show only sychronous spikes with negetive trigger
%syncsign = -1 show only sychronous spikes with negetive trigger
%syncsign = 11 show only sychronous spikes with positive trigger

    
% even if plotting non sync spikes, want to know which are synced
%    if DATA.syncsign < 2
        dt = 2;


        %find synchronous spikes

        for j = 2:length(probes)
            if isfield(DATA,'sids') && sum(ismember(DATA.sids{1},DATA.AllSpikes{probes(1)}.times(ids{1}))) > 100
                aid = DATA.sids{1};
                bid = DATA.sids{2};
            else
            [aid,bid] = FindSync(DATA.AllSpikes{probes(1)}.times(ids{1}),...
                DATA.AllSpikes{probes(j)}.times(ids{j}),dt);
            aid = ids{1}(aid);
            bid = ids{j}(bid);
            end
            dcs(1,1:length(aid)) = mean(DATA.AllSpikes{probes(1)}.values(aid,:),2);
            dcs(j,1:length(bid)) = mean(DATA.AllSpikes{probes(j)}.values(bid,:),2);
            if ismember(DATA.syncsign,[-1 3])
                id = find(DATA.AllSpikes{probes(1)}.values(aid,9) < 0 & ...
                    DATA.AllSpikes{probes(j)}.values(bid,9) < 0);
                aid = aid(id);
                bid = bid(id);
            elseif ismember(DATA.syncsign,[1 4])
                id = find(DATA.AllSpikes{probes(1)}.values(aid,9) > 0 & ...
                    DATA.AllSpikes{probes(j)}.values(bid,9) > 0);
                aid = aid(id);
                bid = bid(id);
            end
            DATA.sids{1} = aid;
            DATA.sids{j} = bid;
            aids{j} = aid;
            if ismember(DATA.plot.synccluster,[3 4 5 6]) % need PCA
                [E, V, DATA.pca] = SpikePCA(DATA.AllSpikes,[probes(1) probes(j)],{aid bid} );
            end
        end
        if length(probes) == 3 %intersection of all three
            [DATA.sids{1}, bi, ci] = intersect(aids{2},aids{3});
            DATA.sids{2} = DATA.sids{2}(bi);
            DATA.sids{3} = DATA.sids{3}(ci);        
        end

        if DATA.plot.prettyfigs ~= 1
        o = (length(probes)-1)/2;
        for j = 1:length(probes)
            if ~isempty(DATA.sids{j})
            mspk(j,:) = mean(DATA.AllSpikes{probes(j)}.values(DATA.sids{j},:),1);
            plot(mspk(j,:)+((j-1-o) * DATA.vstep),':');
            end
        end
        end

end
start = max([DATA.firsttrial 1]);
last = length(Expt.Trials);
% negative Trial numbers indicate manual exclusion
uset = start-1+find([Expt.Trials(start:last).Trial] > 0);
nsync = 0;nspks = 0;
tc = 1;
while tc <= length(uset)  %not for loop, so can change tc inside loop
    trial = Expt.Trials(uset(tc)).Trial;
    stop = get(findobj(DATA.svfig,'Tag','StopSpool'),'value');
    if stop
        set(findobj(DATA.svfig,'Tag','StopSpool'),'value',0);
        if 0
            DATA.playingspk = 0;
            return;
        else
            tc = length(uset);
            trial = Expt.Trials(uset(nt)).Trial;
        end
    end
    if DATA.state.uselfp
        GetFigure('TrialLFP');
    PlotLFPRaw(DATA.state,Expt.Trials(uset(nt)),Expt.Header.LFPsamplerate);
    end
    if timemode
    set(0,'CurrentFigure',DATA.timefig);
    else
    set(0,'CurrentFigure',sfig);
    end
    hold off;
    if mode == 1
        DATA = PlotTrialSpikes(DATA,itrial,colors, clusters);
    elseif mode == 2
        times(1) = Expt.Trials(uset(nt)).Start(1)-DATA.state.preperiod;
        times(2) = Expt.Trials(uset(nt)).End(end)+DATA.state.postperiod;
        Trial = Expt.Trials(uset(nt));
        if isfield(Trial,'uStim') && Trial.uStim > 0
            ustim = 1;
        else
            ustim = 0;
        end
        Trial.ed = GetEval(Expt,'ed',DATA.currenttrial);
        if DATA.syncsign < 2 & length(probes) > 1
            [spks, sspks, cx, cy] = PlotTrialSyncSpikes(DATA, times, [DATA.probe xprobes], colors,'Trial',Trial);
            if DATA.plot.synccluster > 0
            DATA.Spikes.cx(sspks(:,1)) = cx;
            DATA.Spikes.cy(sspks(:,1)) = cy;
            end
            nspks = nspks+length(spks);
            nsync = nsync+size(sspks,2);
            if DATA.plot.synccluster == 10
                    [DATA, spks] = APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probes', [probes(1:2)]);
            end

        else

            if length(probes) > 1
                if DATA.plot.syncoverlay
                for j= 1:length(xprobes)
                    APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe',[xprobes(j) 1+j length(probes)],'lineoff',j*(nclusters+1));
                end
                [DATA, spks] = APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe', [probes(1) 1 length(probes)]);
                end
                if timemode
                for j= 1:length(xprobes)
                    APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe',[xprobes(j) 1+j length(probes)],'lineoff',j*(nclusters+1),'timemode');
                end
                [DATA, spks] = APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe', [probes(1) 1 length(probes)],'timemode');
                end
            else %single electrode
                if isfield(Expt,'probes')
                DATA.probe = Expt.probes(uset(nt));
                probes(1) = DATA.probe;
                end
                if ustim == 0 || (probes(1) > 99 && ustim == 1)
                    [DATA, spks] = APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe', [probes(1) 1 length(probes)],Aargs{:});
                else
                    spks = [];
                end
                if timemode
                    APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe', [probes(1) 1 length(probes)],Aargs{:},'timemode');
                end
            end
        end
        allspks = [allspks spks'];
        nt = nt+1;
    end
    hold on;
    if tc == 1 && DATA.plot.autoscale == 0

    set(0,'CurrentFigure',xyfig);
    set(gca,'Xlim',DATA.plot.clusterXrange,'xlimmode','manual');
    set(gca,'Ylim',DATA.plot.clusterYrange,'ylimmode','manual');
end
    tc = tc+1;

end
if DATA.plot.showsync
    fprintf('%d/%d (%.3f) spikes Synchronized\n',nsync,nspks,nsync./nspks);
end
            

if DATA.probe == 100;
    Expt.MeanPulse = mean(Spks.values(allspks,:));
    plot(Expt.MeanPulse);
end
if DATA.state.uselfp  % look for calibration spikes

 GetFigure('LFP');
 hold off;
 CalcLFPPulse(Expt,DATA.AllData,'plot');
 GetFigure('SpikeV');
end
if isempty(allspks)
    DATA.spkrange(1) = 1;
    DATA.spkrange(2) = 1;
else
    DATA.spkrange(1) = min(allspks);
    DATA.spkrange(2) = max(allspks);
    DATA.spklist = allspks;
end

if isfield(Expt.Header,'Clusterid')
    cluster = Expt.Header.Clusterid;
else
    cluster = 1;
end
DATA.currentcluster = floor(median(cluster)); %in case its mixed

if plotdprime
[dp, details]  = PlotDprime(DATA,0,1);
[dp, a]  = CalcDprime(DATA,cluster,'track');
else
[dp, details]  = CalcDprime(DATA,cluster);
end
if isfield(DATA.AllData.Spikes,'dVdt')
[dps(1), a]  = CalcDprime(DATA,cluster,'xparam',1,'yparam',8);
[dps(2), a]  = CalcDprime(DATA,cluster,'xparam',1,'yparam',44);
[dps(3), a]  = CalcDprime(DATA,cluster,'xparam',1,'yparam',45);
[dps(4), a]  = CalcDprime(DATA,cluster,'xparam',46,'yparam',44);
[dps(5), a]  = CalcDprime(DATA,cluster,'xparam',46,'yparam',45);
else
    dps = [];
end
[dir, name] = fileparts(Expt.Header.Name);
title(sprintf('%scell%d: dprime %.2f',name,Expt.Header.cellnumber,dp));
Expt.gui.spkrange = DATA.spkrange;
Expt.gui.spks = allspks;
Expt.gui.s2clusters = 1+max(Spks.codes(allspks,1));
DATA.s2clusters = Expt.gui.s2clusters;
DATA.dprime = dp;
DATA.dprimes = dps;
if isfield(details,'dps')
DATA.dprimet = details.dps;
end
if isfield(details,'distances')
    DATA.distances = details.distances;
end
%FinishXYPlot(DATA);
%set(xyfig, 'KeyPressFcn',@KeyPressed);
%set(xyfig, 'WindowButtonDownFcn',@ButtonPressed);
%set(xyfig, 'WindowButtonMotionFcn',@ButtonDragged);
%set(xyfig, 'WindowButtonUpFcn',@ButtonReleased);
axis('manual');
%ClearMouse;
DATA.densityplot = 0;
xrange = get(gca,'Xlim');
yrange = get(gca,'Ylim');


if DATA.state.fixrange
    if length(allspks) > 1000
        prc = 99;
    else
        prc = 95;
    end
    xm = prctile(DATA.Spikes.energy(allspks),prc);
    ym = prctile(DATA.Spikes.vw(allspks),prc);
    if yrange(2) > 2 * ym
        set(gca,'Ylim',[0 ym * 1.5]);
    end
    if xrange(2) > 2 * xm
        set(gca,'Xlim',[0 xm * 1.5]);
    end
end
mousept.xrange = diff(xrange);
mousept.yrange = diff(yrange);
fprintf('Spool Took %.2f\n',(now-tstart) * 24 * 60 * 60);
if isfield(DATA,'AllClusters')
    set(0,'currentfig',DATA.xyfig);
    plot(DATA.Spikes.cx(allspks),DATA.Spikes.cy(allspks),'.','markersize',DATA.ptsize);
    DATA.allexp = DATA.currentexpt;
    DATA = PlotAllProbeXY(DATA,0);
end
if isfield(DATA,'AllSpikes') && nprobes == 2 && DATA.plot.xcorr
    GetFigure('Xcorr');
    xc = CalcXcorr(DATA,DATA.currentexpt,probes(1),probes(2));
end

DATA.playingspk = 0;
set(DATA.toplevel,'UserData',DATA);

function ClrSpkWin(caller,b)
%DATA = combine('getstate');
if isstruct(caller)
    DATA = caller;
else
    DATA = GetDataFromFig(caller);
end
GetFigure(DATA.tag.clusterxy);
ym = get(gca,'ylim');
xm = get(gca,'xlim');
hold off;
plot(0,0,'+');
set(gca,'ylim',ym);
set(gca,'xlim',xm);
hold on;
if ~isstruct(caller) % from a mouse button
end

function [dp,res] = PlotSpikeResult(a,b,type)
     DATA = GetDataFromFig(a);

function [dp,res] = PlotDprime(a,b,type)
     
     DATA = GetDataFromFig(a);
     if type == 2
         [dp, res] = CalcDprime(DATA,DATA.currentcluster,'track','smooth',2);
     else
         [dp, res] = CalcDprime(DATA,DATA.currentcluster,'track');
     end
     [f, isnew]= GetFigure('Dprime plot');
     if isnew
         tmp.parentfigtag = DATA.tag.clusterxy;
         set(f,'UserData',tmp);
     end
     hold off;
     plot(res.dps,'-');
     hold on;
     for j = 1:length(res.dps)
         plot(j,res.dps(j),'o','buttondownfcn',{@PlotOneTrial,DATA.Expts{1}.Trials(j).Trial});
     end
     plot(get(gca,'xlim'),[dp dp],'r');

     GetFigure('DprimeDist');
  hold off;
 Cx = DATA.Spikes.cx;
 Cy = DATA.Spikes.cy;
 Spks = DATA.AllData.Spikes;
 cspks = DATA.spklist;
 id = find(Spks.codes(cspks,2) == DATA.currentcluster);
 nid = find(Spks.codes(cspks,2) == 0);
 C.x(1) = mean(Cx(cspks(id)));
 C.x(3) = std(Cx(cspks(id)));
 C.y(1) = mean(Cy(cspks(id)));
 C.y(3) = std(Cy(cspks(id)));
 C.angle = 0;
 x = (Cx(cspks) - C.x(1))./C.x(3);
 y = (Cy(cspks) - C.y(1))./C.y(3);
 xr = x .* cos(C.angle) + y .* sin(C.angle);
 yr = y .* cos(C.angle) - x .* sin(C.angle);
 xc = mean(xr(id));
 yc = mean(yr(id));
 sy = ((yr-mean(yr(id)))./std(yr));
 sx = ((xr-mean(xr(id)))./std(xr));
 o = res.angle;
 d = sx .* cos(o) + sy.* sin(o);
 sd = std(d);
 [y,x] = smhist(d,'sd',sd/20);
 pid = find(y > max(y)/100);
 h = plot(sx,sy,'.','markersize',DATA.ptsize);
 axis('equal');
 refline(tan(o));
 hold on;
 yscale = diff(get(gca,'ylim'))./max(y);
 rx = (x-xc).*cos(-o) + (y-yc).*sin(-o)*yscale;
 ry = (y-yc).*cos(-o)*yscale - (x-xc).*sin(-o);
 plot(rx(pid),ry(pid));
 [y,x] = smhist(d(id),'sd',sd/20);
 pid = find(y > max(y)/100);
 rx = (x-xc).*cos(-o) + (y-yc).*sin(-o)*yscale;
 ry = (y-yc).*cos(-o)*yscale - (x-xc).*sin(-o);
 plot(rx(pid),ry(pid),'r');
 [y,x] = smhist(d(nid),'sd',sd/20);
 pid = find(y > max(y)/100);
 rx = (x-xc).*cos(-o) + (y-yc).*sin(-o)*yscale;
 ry = (y-yc).*cos(-o)*yscale - (x-xc).*sin(-o);
 plot(rx(pid),ry(pid),'g');

        
function PlotOneTrial(a,b, trial)
    DATA = GetDataFromFig(a);
    t = find([DATA.Expts{1}.Trials.Trial] == trial)
    PlayOneTrial(DATA,t, 0);
    
         
function [dprime, details] = CalcDprime(DATA,cluster, varargin)
    
    trackcl = 0;
    yparam = 0;
    xparam = 0;
    smw = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'track',4)
            trackcl = 1;
        elseif strncmpi(varargin{j},'smooth',4)
            j = j+1;
            smw = varargin{j};
        elseif strncmpi(varargin{j},'xparam',4)
            j = j+1;
            xparam = varargin{j};
        elseif strncmpi(varargin{j},'yparam',4)
            j = j+1;
            yparam = varargin{j};
        end
        j = j+1;
    end
    spkr = minmax(DATA.spklist);
        ispk = DATA.spklist;
    if yparam == 0
        yr = DATA.Spikes.cy(DATA.spklist);
    else
        yr = GetSpikeVals(DATA, ispk, DATA.AllData.Spikes.values(ispk,:), DATA.AllData.Spikes.dVdt(ispk,:),yparam, 1,[]);
    end
    if xparam == 0
    xr = DATA.Spikes.cx(DATA.spklist);
    else
        xr = GetSpikeVals(DATA, ispk, DATA.AllData.Spikes.values(ispk,:), DATA.AllData.Spikes.dVdt(ispk,:),xparam, 1,[]);
    end
    if length(cluster) > 1
        cluster = floor(median(cluster)); %Temporary. Need to make an array of length spklist with correct cluster in 
        fprintf('More than one clusterid in this cell\n');
    end
    
    id = find(DATA.AllData.Spikes.codes(DATA.spklist,2) == cluster);
    nid = find(DATA.AllData.Spikes.codes(DATA.spklist,2) == 0);
    sy = ((yr-mean(yr(id)))./std(yr));
    sx = ((xr-mean(xr(id)))./std(xr));
    sd = abs(sx+i*sy);
    o = [0:pi/40:pi];
    for j = 1:length(o)
        d = sx .* cos(o(j)) + sy.* sin(o(j));
        ds{j} = d;
%
%        dprimes(j) = abs(mean(d(nid(tid)))-mean(d(id)))./sqrt(mean([var(d(nid(tid))) var(d(id))]));
         dprimes(j) = abs(mean(d(nid))-mean(d(id)))./sqrt(mean([var(d(nid)) var(d(id))]));
%        dprimes(j) = abs(mean(d(nid))-mean(d(id)))./std(d(nid));
    end
    [dprime, maxi] = max(dprimes);
    details.distances(DATA.spklist) = ds{maxi};
    details.id = id;
    details.angle = o(maxi);
    if trackcl
        d = zeros(size(DATA.AllData.Spikes.codes,1),1);
        d(DATA.spklist) = ds{maxi};
        T = DATA.Expts{1}.Trials;
        for j = 1:length(T);
            ispk = [];
            codes = [];
            for k = max([1 j-smw]):min([length(T) j+smw])
                [a, ts, b] = FindSpikes(DATA,[T(k).Start(1) T(k).End(end)+500], DATA.probe, spkr);
                ispk = [ispk; a];
                codes = [codes; b];
            end
            nid = ispk(find(codes == 0));
            cid = ispk(find(codes == cluster));
            dp(j) = abs(mean(d(nid))-mean(d(cid)))./sqrt(mean([var(d(nid)) var(d(cid))]));
        end
        details.dps = dp;
    end

function [ispk, spktimes, codes] = FindSpikes(DATA, times, probe, range)

    spktimes = [];
    codes = [];
    cytpe = 1;
if isfield(DATA,'AllSpikes')
    if isempty(DATA.AllSpikes{probe}) | ~isfield(DATA.AllSpikes{probe},'times')
        ispk = [];
        return;
    end
    ispk = find(DATA.AllSpikes{probe}.times > times(1) &...
    DATA.AllSpikes{probe}.times < times(2));
    if nargout > 1
    spktimes = DATA.AllSpikes{probe}.times(ispk);
    codes = DATA.AllSpikes{probe}.codes(ispk,2);
    end
elseif isfield(DATA,'AllClusters')
    ispk = find(DATA.AllClusters(probe).times > times(1) &...
    DATA.AllClusters(probe).times < times(2));
elseif isempty(DATA.AllData.Spikes) || ~isfield(DATA.AllData.Spikes,'times')
    ispk = [];
else
    if size(DATA.AllData.Spikes.codes,2) == 1
        ctype = 1;
    end
    if ~isempty(range)
        ispk = find(DATA.AllData.Spikes.times > times(1) &...
            DATA.AllData.Spikes.times < times(2));
    else
        ispk = find(DATA.AllData.Spikes.times(range) > times(1) &...
            DATA.AllData.Spikes.times(range) < times(2));
        ispk = range(ispk);
    end
    spktimes = DATA.AllData.Spikes.times(ispk);
    codes = DATA.AllData.Spikes.codes(ispk,2);
end
        
function [DATA, ispk] = APlotTrialSpikes(DATA, times, colors, nclusters, classify, varargin)

j = 1;
Trial = [];
probe = DATA.probe;
lineoff = 0;
step = DATA.plot.SpikeMaxV;
step = DATA.vstep;
timemode = 0;
nprobes = 1;
voff = NaN;
ip = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'lineoff',6)
        j = j+1;
        lineoff = varargin{j};
    elseif strncmpi(varargin{j},'Probes',6)
        j = j+1;
        probe = varargin{j};
    elseif strncmpi(varargin{j},'Probe',4)
        j = j+1;
        probe = varargin{j}(1);
        if length(varargin{j}) > 1
            ip = varargin{j}(2); % # in list
        end
        if length(varargin{j}) > 2
            nprobes = varargin{j}(3);
        end
    elseif strncmpi(varargin{j},'timemode',4)
        timemode = 1;
    elseif strncmpi(varargin{j},'Trial',4)
        j = j+1;
        Trial = varargin{j};
    elseif strncmpi(varargin{j},'voff',4)
        j = j+1;
        voff = varargin{j};
    end
    j = j+1;
end
if isfield(DATA,'AllSpikes')
    Spks = DATA.AllSpikes{probe(1)};
    PCs = [];
elseif DATA.probe == 100
    Spks = DATA.AllData.UstimV;
else
    Spks = DATA.AllData.Spikes;
    if isempty(DATA.AllData.pcs)
        PCs = DATA.AllData.Spikes.codes;
    else
    PCs = DATA.AllData.pcs;
    end
end

if DATA.state.recut && size(Spks.codes,2) > 1
    ctype = 2;
else
    ctype = 1;
end

if ~isfield(Spks,'values')
    ispk = [];
    return;
end
splen = size(Spks.values,2);
% try calculating energy for all spikes in Expt in one step.
ispk = find(Spks.times > times(1) &...
        Spks.times < times(2));
if isfield(DATA,'AllClusters')
    return;
end
    if ispk
        
    if length(PCs) < max(ispk)
        PCs = Spks.codes;
    end
if isfield(DATA,'sids') 
    syncspk = find(ismember(ispk,DATA.sids{ip}));
    syncspklst = ismember(ispk,DATA.sids{ip});
    if length(probe) == 2
        Spks.values = (DATA.AllSpikes{probe(1)}.values(DATA.sids{1},:)+DATA.AllSpikes{probe(2)}.values(DATA.sids{2},:))/2;
        Spks.times =  DATA.AllSpikes{probe(1)}.times(DATA.sids{1});
        ispk = find(Spks.times > times(1) &...
            Spks.times < times(2));
    end
    if ismember(DATA.syncsign,[3 4 5])
    Spks.codes(ispk,3) = 0;
    Spks.codes(ispk(syncspk),3) = 1;
    ctype = 3;
    end
end

if DATA.probe < 100
    if isfield(Spks,'dVdt')
        [cx, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterX, classify,PCs(ispk,:));
        DATA.Spikes.cx(ispk) = cx;
  %      [cy, DATA] = GetSpikeVals(DATA,ispk, SPKVARE, classify);
        [cy, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterY, classify,PCs(ispk,:));
        DATA.Spikes.cy(ispk) = cy;
            dvdt = Spks.dVdt(ispk,:);
        if length(ispk)
        DATA.currentspike = ispk(1);
        end
    elseif isfield(Spks,'xy')
        if max(ispk) > size(Spks.xy,1)
            fprintf('Spike array mismatch, %d vs %d\n',max(ispk),size(Spks.xy,1));
            ispk = ispk(ispk < size(Spks.xy,1));
        end
        cx = Spks.xy(ispk,1);
        cy = Spks.xy(ispk,2);
    end
            adc = Spks.values(ispk,:);
% recut == 2 means that clusters are not set here, but clusters have
% been defined (previous list), so use those properites.
         if DATA.state.recut == 2 && DATA.probe == probe
             DATA = SetSpkCodes(DATA,ispk, DATA.probe,0);
             if isfield(DATA,'AllSpikes') %copy any changes into spks
                 Spks.codes = DATA.AllSpikes{DATA.probe}.codes;
             else
                 Spks.codes = DATA.AllData.Spikes.codes;
             end
         end
else
    Spks.codes = DATA.AllData.UstimV.codes;
            adc = Spks.values(ispk,:);
            cx = zeros(size(ispk));
            cy = zeros(size(ispk));
end

         
        for j = [1:nclusters+1]
            vs{j} = [];
            xs{j} = [];
        end
        for spk = 1:length(ispk);
            if DATA.syncsign < 2 & isfield(DATA,'sids')
                j = syncspklst(spk)+1;
            else
                j = Spks.codes(ispk(spk), ctype)+1;
            end
            if DATA.plot.dvdt ==2
                vs{j} = [vs{j} dvdt(spk,:) NaN];
                xs{j} = [xs{j} adc(spk,1:splen-1) NaN];
            elseif DATA.plot.dvdt
                vs{j} = [vs{j} dvdt(spk,:) NaN];
                xs{j} = [xs{j} [1:splen-1] NaN];
            elseif timemode
                vs{j} = [vs{j} adc(spk,:) NaN];
                xs{j} = [xs{j} [1:splen]+(Spks.times(ispk(spk))-times(1)).*timemode NaN];
            elseif DATA.plot.nodc
            vs{j} = [vs{j} adc(spk,:)-mean(adc(spk,:)) NaN];
            xs{j} = [xs{j} [1:splen] NaN];
            else
            vs{j} = [vs{j} adc(spk,:) NaN];
            xs{j} = [xs{j} [1:splen] NaN];
            end
        end
        if length(DATA.plot.voffsets) >= ip
            for j = [1:nclusters+1]
                vs{j} = vs{j}+DATA.plot.voffsets(ip);
            end
        else
        for j = [1:nclusters+1]
            vs{j} = vs{j}+(ip-(nprobes+1)/2).*step;
        end
        end
            
        if timemode && isfield(DATA,'timefig')
        set(0,'CurrentFigure',DATA.timefig);
        text(0,(ip-(nprobes+1)/2).*step,sprintf('%d',probe));
        vh = DATA.tvh;
        else
        set(0,'CurrentFigure',DATA.svfig);
        vh = DATA.svh;
        end
        nc = min([nclusters+1 length(DATA.svh)]);
  
    
        for j = 1:nc
            k = j+lineoff;
            if k > length(vh)
                    vh(k) = line('Xdata' , xs{j}, 'Ydata', vs{j});
            elseif ~isempty(xs{j}) 
                if ishandle(vh(k))
                    set(vh(k),'Xdata' , xs{j}, 'Ydata', vs{j});
                else
                    vh(k) = line('Xdata' , xs{j}, 'Ydata', vs{j});
                end
            elseif vh(k) && ishandle(vh(k)) %no spikes, wipe clean
                    set(vh(k),'Xdata' , 0, 'Ydata', 0);
            end
        end
        if ~isempty(Trial)
%            xc = ExtraLabels(Trial);
xc = '';
            nspk = sum(Spks.codes(ispk,2) == DATA.currentcluster);
            title(sprintf('Trial %d (id%d %.2f - %.2f) ed%.3f%s %d/%d spks P%d',abs(Trial.Trial),...
                Trial.id,Trial.Start(1)./10000,Trial.End(end)./10000,Trial.ed,xc,nspk,length(ispk),DATA.probe));
        end        

        if probe == DATA.probe
        drawnow;
        set(0,'CurrentFigure',DATA.xyfig);
        for j = 0:nclusters
            sp = find(Spks.codes(ispk, ctype) == j);
            plot(cx(sp),cy(sp),...
                '.','color',DATA.spkcolor{j+1},'markersize',DATA.ptsize);
            hold on; %% need this when called from PlotOneTrial
        end
        end
    end
    DATA.minplottime = 0.00;
    if DATA.minplottime > 0.001
        while toc < DATA.minplottime
        end
    end

    function [x,DATA] = GetSpikeVals(DATA, ispk, values, dVdt, type, recalc, pcs)

%Do NOT change the order of these definiitons. Cluster params are saved
%as integers....
SPKENERGY=1;
SPKVARE=2;
SPKMAXRATE = 3;
SPKPREMINRATE = 4;
SPKPEAKTIME = 5;
SPKPREMIN = 6;
SPKPEAK = 7;
SPKVAR = 8;
SPKSYMMETRY = 9;
SPKCENTROID = 10;

SPKMINRATE = 11;
SPKMAXRATEA = 12;
SPKMINRATEA = 13;
SPKMEANRATEA = 14;
SPKMAXRATEB = 15;
SPKMINRATEB = 16;
SPKMEANRATEB = 17;
SPKMEANA = 18;
SPKMEANB = 19;
SPKMAXA = 20;
SPKMAXB = 21;
SPKMIN =22;
SPKMINA = 23;
SPKMINB = 24;
SPKHEIGHT = 25;
SPKENERGYA=26;
SPKENERGYB=27;
SPKPREMAXRATE = 28;
TEMPLATEA = 29;
TEMPLATEB = 30;
TEMPLATEC = 31;
SPKSYMNEG = 32;
SPKPCA1= 33;
SPKPCA2= 34;
SPKVARTESTA= 35;
SPKVARTESTB= 36;
SPKPCA3= 37;

SPKPCA4= 38;
SPKPCA12= 39;
SPKMEAN =40;
SPKMAXACCEL = 41;
SPKMINACCEL = 42;
SPKMEANACCEL = 43;
SPKVAREA = 44;
SPKVAREB = 45;
SPKRMSV = 46;

if isnan(values) %% return names, order
    x = {'Energy' 'Var/Energy' 'MaxRate' 'PreMinRate' 'PeakTime' 'PreMin' 'Peak' 'Var' 'Symmetry' 'Centroid' 'Minrate' 'Maxrate(A)' ... 
    'Minrate(A)' 'Meanrate(A)' 'Maxrate(B)' 'Minrate(B)' 'Meanrate(B)' 'Mean(A)'...
    'Mean(B)' 'Max(A)' 'Max(B)' 'Min' 'Min(A)' 'Min(B)' 'Height' 'EnergyA' 'EnergyB' 'PreMaxRate' ...
    'Template 1' 'Template 2' 'Template 3' '-Symmetry' 'PCA1' 'PCA2' 'test1' 'test2' 'PCA3' 'PCA4' 'PCA1-PCA2' 'Mean' 'AccelMax' 'AccelMin' 'AccelMean'...
    'Var/sqrt(energy)' 'sqrt(Var/Energy)' 'sqrt(Energy)'};
   DATA = [SPKENERGY SPKVARE SPKMAXRATE SPKVAR SPKPEAK SPKMIN SPKPREMINRATE SPKPEAKTIME SPKPREMIN ...
       SPKSYMMETRY SPKCENTROID SPKMINRATE ...
       SPKMAXRATEA SPKMINRATEA SPKMEANRATEA SPKMAXRATEB SPKMINRATEB ...
    SPKMEANRATEB SPKMEAN SPKMEANA SPKMEANB SPKMAXA SPKMAXB SPKMIN SPKMINA SPKMINB ...
    SPKHEIGHT SPKENERGYA SPKENERGYB SPKPREMAXRATE TEMPLATEA TEMPLATEB ...
    TEMPLATEC SPKSYMNEG SPKPCA1 SPKPCA2 SPKPCA3 SPKPCA4 SPKPCA12  SPKMEAN SPKMAXACCEL SPKMINACCEL SPKMEANACCEL SPKVAREA SPKVAREB SPKRMSV];
return;
end


p = DATA.probe;
if p > size(DATA.cluster,2)
    DATA.cluster{DATA.currentcluster,p} = {};
end

arange = DATA.clusterArange;
brange = DATA.clusterBrange;
splen = size(values,2);
erange = intersect(DATA.clusterErange,1:splen-1);    
    

%dVdt and values are already just teh values for ispk
%ispk is only sent so that the correct entries in DATA.Spikes are 
%filled in
if type == SPKENERGY
    if recalc
        x  = sum(dVdt(:,erange).^2,2);
        x = x';
        DATA.Spikes.energy(ispk)= x;
    else
        x = DATA.Spikes.energy(ispk);
    end
elseif type == SPKRMSV
        x = sqrt(DATA.Spikes.energy(ispk));
elseif type == SPKENERGYA
    [mins, imins] = min(values'); 
    imins(find(imins<6)) = 6;
    imins(find(imins>26)) = 26;
    for j = 1:size(dVdt,1)
        x(j)  = sum(dVdt(j,imins(j)-5:imins(j)+5).^2,2);
    end
        DATA.Spikes.energy(ispk)= x;
elseif type == SPKVARE
    if recalc
        x = var(values(:,erange)');
        x = x./DATA.Spikes.energy(ispk);
        DATA.Spikes.vw(ispk) = x;
    else
        x = DATA.Spikes.vw(ispk);
    end
elseif type == SPKVAREA
    x = var(values(:,erange)');
    x = x./sqrt(DATA.Spikes.energy(ispk));
elseif type == SPKVAREB
    x = var(values(:,erange)');
    x = sqrt(x)./sqrt(DATA.Spikes.energy(ispk));
elseif type == SPKMAXRATE
    x = max(dVdt(:,:)');
elseif type == SPKPREMINRATE %min rate value preceding max rate value;
    [x,t] = max(dVdt(:,:)');
    for j = 1:length(t)
        x(j) = min(dVdt(j,1:t(j)));
    end
elseif type == SPKPREMAXRATE %MAX rate preceding peak V
    [x,t] = max(values(:,:)');
    t(find(t > size(dVdt,2))) = size(dVdt,2);
    for j = 1:length(t)
        x(j) = max(dVdt(j,1:t(j)));
    end
elseif type == SPKMAXRATEA
    x = max(dVdt(:,arange)');
elseif type == SPKMAXRATEB
    x = max(dVdt(:,brange)');
elseif type == SPKMINRATE
    x = min(dVdt(:,:)');
elseif type == SPKMINRATEA
    x = min(dVdt(:,arange)');
elseif type == SPKMINRATEB
    x = min(dVdt(:,brange)');
elseif type == SPKMEANA
    x = mean(values(:,arange)');
elseif type == SPKMEANB
    x = mean(values(:,brange)');
elseif type == SPKMEANRATEA
    x = mean(dVdt(:,arange)');
elseif type == SPKMEANRATEB
    x = mean(dVdt(:,brange)');
elseif type == SPKMAXA
    x = max(values(:,arange)');
elseif type == SPKPCA1
    x = pcs(:,1);
elseif type == SPKPCA2
    x = pcs(:,2);
elseif type == SPKPCA3
    x = pcs(:,3);
elseif type == SPKMAXB
    x = max(values(:,brange)');
elseif type == SPKMIN
    x = min(values(:,:)');
    if DATA.plot.nodc
        x = x - mean(values');
    end
elseif type == SPKMINA
    x = min(values(:,arange)');
elseif type == SPKMINB
    x = min(values(:,brange)');
elseif type == SPKPREMIN
    [x,t] = max(dVdt(:,:)');
    for j = 1:length(t)
        x(j) = min(values(j,1:t(j)));
    end
elseif type == SPKPEAKTIME
    [x,t] = max(dVdt(:,:)');
    zc = diff(sign(dVdt(:,:)')); 
    len = size(zc,1);
    for j = 1:length(t)
        tm = find(zc(:,j) < 0 & [1:len]' >= t(j));
        if isempty(tm)
            x(j) = 0;
        else
            x(j) = tm(1); % first zero crossing after peakrate
        end
    end
elseif type == SPKPEAK
    [x,t] = max(values(:,:)');
elseif type == SPKHEIGHT
    [x,t] = max(values(:,:)');
    [y,t] = min(values(:,:)');
    x = x-y;
elseif type == SPKVAR
    x = var(values(:,:)');
elseif type == SPKSYMMETRY || type == SPKSYMNEG
  len =  size(dVdt,2);
  cn  = round(sum((dVdt(:,:).^2)*[1:len]',2)./sum(dVdt(:,:).^2,2));
  for j = 1:length(ispk)
  z(j) = mean(values(j,max([cn(j)-5 1]):cn(j)),2);
  y(j) = mean(values(j,cn(j):min([cn(j)+5 len])),2);
  end
  if type == SPKSYMNEG
  x = atan2(-y,-z);
  else
  x = atan2(y,z);
  end
elseif type == SPKCENTROID
  len =  size(dVdt,2);
  x  = sum((dVdt(:,:).^2)*[1:len]',2)./sum(dVdt(:,:).^2,2);
elseif ismember(type, [TEMPLATEA TEMPLATEB TEMPLATEC])
    j = 1+type-TEMPLATEA;
    x = values * DATA.Templates(j,:)';
end

function [E,V, pc] = SpikePCA(spk, probes, ids)

if length(probes) == 2
    spks = [spk{probes(1)}.values(ids{1},:) spk{probes(2)}.values(ids{2},:)];
         [E,V] = eig(cov(spks));
        pc(1,:) = (E(:,end)'-mean(E(:,end))) * spks';
        pc(2,:) = (E(:,end-1)'-mean(E(:,end-1))) * spks';
        pc(3,:) = (E(1:32,end)'-mean(E(1:32,end-1))) * spks(:,1:32)';
        pc(4,:) = (E(33:end,end)'-mean(E(33:end,end-1))) * spks(:,33:end)';
        pc(5,:) = (E(1:32,end-1)'-mean(E(1:32,end-1))) * spks(:,1:32)';
        pc(6,:) = (E(33:end,end-1)'-mean(E(33:end,end-1))) * spks(:,33:end)';
end

function DATA = CheckForPCA(DATA, ispk, force)
   if sum(ismember([DATA.plot.clusterX DATA.plot.clusterY DATA.plot.clusterZ],[33 34 37])) %need PCA
       if isfield(DATA,'AllSpikes')
       elseif length(DATA.AllData.pcs) < max(ispk)
        [a,b,c] = OneProbePCA(DATA,ispk);
        DATA.AllData.pcs(ispk,:) = c';
       end
   end


function [E,V, pc] = OneProbePCA(DATA, ids)

        spks = DATA.AllData.Spikes.values(ids,:);
         [E,V] = eig(cov(spks));
        pc(1,:) = (E(:,end)'-mean(E(:,end))) * spks';
        pc(2,:) = (E(:,end-1)'-mean(E(:,end-1))) * spks';
        pc(3,:) = (E(:,end-2)'-mean(E(:,end-2))) * spks';
        pc(4,:) = (E(:,end-3)'-mean(E(:,end-3))) * spks';

        
        function PlayLastTrial(a, b)
DATA = GetDataFromFig(a);
PlayOneTrial(DATA, a,-1);

function PlayNextTrial(a, b)
DATA = GetDataFromFig(a);
PlayOneTrial(DATA,a,1);

function SelectTrial(a, b)
DATA = GetDataFromFig(a);
c = get(findobj(a,'Tag','ChooseTrial'),'value');
PlayOneTrial(DATA,c,0);

function PlayOneTrial(DATA, a, b)
%DATA = combine('getstate');
expid = DATA.currentexpt;

if b < 0  & DATA.currenttrial > 1 %step back
    DATA.currenttrial = DATA.currenttrial-1;
    while DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial < 0 & DATA.currenttrial > 1 
        DATA.currenttrial = DATA.currenttrial-1;
    end
elseif b > 0 & DATA.currenttrial < length(DATA.Expts{DATA.currentexpt}.Trials) %step back
    DATA.currenttrial = DATA.currenttrial+1;
    while DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial < 0 &  ...
            DATA.currenttrial < length(DATA.Expts{DATA.currentexpt}.Trials) 
        DATA.currenttrial = DATA.currenttrial+1;
    end
elseif b == 0 % step to a trial on list
    %could match Trial # in trial
    DATA.currenttrial = find([DATA.Expts{DATA.currentexpt}.Trials.Trial] == a);    
    %but in fact use # of trial in expt
    DATA.currenttrial = a;    
end
if DATA.currenttrial <= 0
    DATA.currenttrial = 1;
end
Trial = DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial);
if Trial.Trial < 0 && b == 0 %%manually select -Trial = include again
    DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial = abs(Trial.Trial);
    it = findobj(DATA.svfig,'Tag','ChooseTrial');
    set(it,'string',sprintf('%d|',[DATA.Expts{expid}.Trials.Trial]),'value',1);
end
%itrial = find(DATA.AllData.Trialids == abs(Trial.Trial));
set(0,'CurrentFigure',DATA.svfig);
hold off;
%DATA = PlotTrialSpikes(DATA, itrial, mycolors, DATA.clusters);
if DATA.state.recut
    nc = size(DATA.cluster,1)+1;
else
    nc = DATA.s2clusters;
end
Trial.ed = GetEval(DATA.Expts{DATA.currentexpt},'ed',DATA.currenttrial);
probes =  [DATA.probe DATA.xprobes];
if (length(probes) > 1 & DATA.plot.syncoverlay == 0) || DATA.plot.showwave
    Aargs = {'timemode'};
    if  ~isfield(DATA,'timefig')
        DATA.timefig = GetFigure('SpikeTime');
    end
    set(0,'CurrentFigure',DATA.timefig);
else
    Aargs = {};
end
 
if DATA.syncsign < 2
    PlotTrialSyncSpikes(DATA, [Trial.Start(1) Trial.End(end)], [DATA.probe DATA.xprobes], mycolors,'Trial',Trial);
else
for j = 2:length(probes)
DATA = APlotTrialSpikes(DATA, [Trial.Start(1) Trial.End(end)], mycolors, nc, 1,'Trial',Trial,'probe',[probes(j) j length(probes)],'lineoff',(j-1)*(1+DATA.nclusters), Aargs{:});
end
[DATA, ispk] = APlotTrialSpikes(DATA, [Trial.Start(1) Trial.End(end)], mycolors, nc, 1,'Trial',Trial,'probe',[probes(1) 1 length(probes)],Aargs{:});
end
if length(probes) > 1 & isfield(DATA,'timefig') %also show spikes over time
    set(0,'CurrentFigure',DATA.timefig);
    vh = DATA.tvh;
    if DATA.plot.timebyspikeprobe > 0
        p = DATA.probe;
        tw = 100;
        set(gca,'xlim',[0 tw*2]);
        [sid, stimes, codes] = FindSpikes(DATA, [Trial.Start(1) Trial.End(end)], p,[]);
        sid = sid(codes ==1);
        stimes = stimes(codes ==1);
        dt = [-110:10:110];
        dts = [];
        for j = 1:length(probes)
            [sids{j}, alltimes{j} allcodes{j}] = FindSpikes(DATA, [Trial.Start(1) Trial.End(end)], probes(j),[]);
            voff(j) = DATA.plot.SpikeVsep*(j-1);
            diffs{j} = [];
              for k = 1:length(stimes)
                  dts(k,:) = hist(stimes(k) - alltimes{j},dt);
              end
              corrs(j,:) = sum(dts);
        end
        nsmp = size(DATA.AllSpikes{probes(1)}.values,2);
        hscale = DATA.AllSpikes{probes(1)}.interval * 10000;
        hold off;
        for k = 1:length(sid)
            ts = stimes(k);
        for j = 1:length(probes)
            x = [];
            y = [];
            id = find(alltimes{j} > ts-tw & alltimes{j} < ts+tw);
            for t = 1:length(id)
                x = (alltimes{j}(id(t))-ts+tw) +[1:nsmp]*hscale;
                y = DATA.AllSpikes{probes(j)}.values(sids{j}(id(t)),:)+voff(j);
                plot(x,y,'color',DATA.spkcolor{allcodes{j}(id(t))+1});
                diffs{j} = [diffs{j} alltimes{j}(id(t))-ts];
                hold on;
            end
        end
%        hold off;
        set(gca,'xlim',[0 tw*2],'ylim',[-DATA.plot.SpikeMaxV max(voff)+DATA.plot.SpikeMaxV]);
        drawnow;
        pause(0.05);
        end
        DATA.currentspike = sid(1);
    else
    for j = 1:length(probes)
        DATA = APlotTrialSpikes(DATA, [Trial.Start(1) Trial.End(end)], mycolors, nc, 0,'Trial',Trial,'probe',[probes(j) j length(probes)],'lineoff',(j-1)*(1+DATA.nclusters),'timemode');
    end
    end
end
if DATA.state.uselfp
   GetFigure('TrialLFP');
    PlotLFPRaw(DATA.state,Trial,DATA.Expts{expid}.Header.LFPsamplerate);
end

if isfield(DATA,'distances')
cspk = find(DATA.AllData.Spikes.codes(ispk,2) == 1);
nspk = find(DATA.AllData.Spikes.codes(ispk,2) == 0);
d = DATA.distances;
dp = abs(mean(d(ispk(cspk)))-mean(d(ispk(nspk))))./sqrt(mean([var(d(ispk(cspk))) var(d(ispk(nspk)))]));
fprintf('Trial %d DPrime %.2f\n',Trial.Trial,dp);
end
%CheckSpoolButton(DATA);
set(DATA.toplevel,'UserData',DATA);
set(0,'CurrentFigure',DATA.svfig);
tid = DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).id;
it = findobj(DATA.svfig,'Tag','ChooseTrial');
set(it,'value',DATA.currenttrial);
    
function RescaleClusterPlot(a,b)
    DATA = GetDataFromFig(a);
it = findobj(DATA.xyfig,'Tag','ClusterXmax');
if it
    DATA.plot.clusterXrange(2) = str2num(get(it,'string'));
end

it = findobj(DATA.xyfig,'Tag','ClusterYmax');
if it
    DATA.plot.clusterYrange(2) = str2num(get(it,'string'));
end
it = findobj(DATA.xyfig,'Tag','ClusterZmax');
if it
    DATA.plot.clusterZrange(2) = str2num(get(it,'string'));
end
ax = findobj(DATA.xyfig,'type','ax');

%
%if the user has typed in a value, surely they want to manually scales
DATA.plot.autoscale = 0;
set(findobj(DATA.xyfig,'Tag','AutoScale'),'value',DATA.plot.autoscale);


    
    
if DATA.densityplot
    caxis([0 DATA.plot.clusterZrange(2)]);
end
set(DATA.toplevel,'UserData',DATA);

function SpoolExptSpikes(a,b)
DATA = GetDataFromFig(a);
DATA.spooling = 1;
if strmatch(get(a,'Tag'),'SpoolSpikes') & DATA.currenttrial > 1
    DATA.spooling = 2; %spool from current spike to end
end
DATA = PlaySpikes(DATA,DATA.Expts{DATA.currentexpt});
set(DATA.toplevel,'UserData',DATA)

function ScrollTrial(src, evnt)

DATA = GetDataFromFig(src);

if src ~= gcf
    return;
end

DATA = GetDataFromFig(src);
PlayOneTrial(DATA,src,sign(evnt.VerticalScrollCount));


function SetXYCluster(a,b, type, iplot, val)
DATA = GetDataFromFig(a);
recalcxy = 0;
if type == 1
DATA.plot.clusterX = val;
recalcxy = 1;
elseif type == 2
DATA.plot.clusterY = val;
recalcxy =1;
elseif type == 3
DATA.plot.clusterZ = val;
end
ispk = DATA.spklist;
if (iplot == 3 | type == 3) & DATA.plot.clusterZ > 0
    Plot3DClusters(DATA,recalcxy);
else
        DATA = CalcClusterVars(DATA, ispk);
        hold off;
    DATA = DrawXYPlot(DATA, ispk); %returned DATA has chnges to clusterX/Yrange
end
set(DATA.toplevel,'UserData',DATA);


function DATA = CheckPtSize(DATA, nspk)
    if DATA.plot.setptsize
        DATA.ptsize = DATA.plot.setptsize;
    elseif ~isfield(DATA,'ptsize')
        DATA.ptsize = 1;
    elseif nspk < 10000
        DATA.ptsize = 6;
    elseif nspk < 2000
        DATA.ptsize = 10;
    else
        DATA.ptsize = 4;
    end

function PlotSpikeXY(DATA, spkid, color)
CheckPtSize(DATA, length(spkid));
plot(DATA.Spikes.cx(spkid),DATA.Spikes.cy(spkid),'.',...
    'color',color,'markersize',DATA.ptsize);

function DATA = DrawXYPlot(DATA, expspks)
    
    ho = ishold;
    if ~isfield(DATA,'nclusters')
        DATA.nclusters = 0;
    end
    if DATA.state.recut
        ctype = 2;
    else
        ctype = 1;
    end
    DATA = CheckPtSize(DATA,length(expspks));

    expid = DATA.currentexpt;
    if isfield(DATA,'AllSpikes')
        expspks = DATA.Expts{DATA.currentexpt}.gui.spks;
        if ~isfield(DATA,'nclusters')
            DATA.nclusters = 0;
        end
        for j = 1:DATA.nclusters+1
            id = find(DATA.AllSpikes{DATA.probe}.codes(expspks,ctype) == j-1);
            PlotSpikeXY(DATA, expspks(id), DATA.spkcolor{j});
            hold on;
        end
    else
        if length(expspks) < 10 
            expspks = DATA.Expts{DATA.currentexpt}.gui.spks;
        end
        if DATA.densityplot
            PlotXYDensity(DATA.Spikes.cx(expspks),DATA.Spikes.cy(expspks));
        else
            for j = 1:DATA.nclusters+1
                id = find(DATA.AllData.Spikes.codes(expspks,ctype) == j-1);
                PlotSpikeXY(DATA, expspks(id), DATA.spkcolor{j});
                hold on;
            end
        end
    end
    hold on;
%    DATA = DrawClusters(DATA, DATA.cluster, 0);

%        CalcClusterdp(DATA,1);
    if ismember(DATA.plot.autoscale,[2 3]) & length(DATA.Expts{DATA.currentexpt}.gui.spks)
        expspks = DATA.Expts{DATA.currentexpt}.gui.spks;
        DATA.spklist = expspks;
        [xr, yr] = ClusterRange(DATA.cluster,DATA.probe);
        cx = DATA.Spikes.cx(expspks);
        cy = DATA.Spikes.cy(expspks);
        if DATA.plot.autoscale == 2
        y(2) = prctile(cy,99.5).*2;
        x(2) = prctile(cx,99.5).*2;
        else
        y(2) = prctile(cy,95).*2;
        x(2) = prctile(cx,95).*2;
        end
        y(2) = min([y(2) max(cy)]);
        y(2) = max([y(2) yr(2)]); %make sure cluster ellipse is visible
        x(2) = min([x(2) max(cx)]);
        x(2) = max([x(2) xr(2)]);
        if min(cy) > 0 & min(cy) < y(2)/10
            y(1) = 0;
        else
            y(1) = min(cy);
        end
        if min(cx) > 0 & min(cx) < x(2)/10
            x(1) = 0;
        else
            x(1) = min(cx);
        end
        if x(2) <= x(1)
            x(2) = x(1) + abs(x(1));
        end
        if y(2) <= y(1)
            y(2) = y(1) + abs(y(1));
        end
        DATA.plot.clusterYrange = y;
        DATA.plot.clusterXrange = x;
    end
    DATA = FinishXYPlot(DATA);
    if ~ho
        hold off;
    end

function DATA = SetXYRanges(DATA, cx, cy)

    
if nargin == 1
    cx = DATA.Spikes.cx(DATA.spklist);
    cy = DATA.Spikes.cy(DATA.spklist);
end

if DATA.plot.autoscale == 0
    return;
elseif DATA.plot.autoscale == 1
    DATA.plot.clusterXrange  = get(gca,'Xlim');
    DATA.plot.clusterYrange  = get(gca,'Ylim');
elseif ismember(DATA.plot.autoscale,[2 3 4])
    if min(cx) > 0
        minx = 0;
    else
        minx = prctile(cx,1 * 1.1);
    end
    if min(cy) > 0
        miny = 0;
    else
        miny = prctile(cy,1 * 1.1);
    end
    if DATA.plot.autoscale == 2
    DATA.plot.clusterYrange = [miny prctile(cy,99.9).*1.2];
    DATA.plot.clusterXrange = [minx prctile(cx,99.9).*1.2];
    elseif DATA.plot.autoscale == 3
    DATA.plot.clusterYrange = [miny prctile(cy,99.5).*1.2];
    DATA.plot.clusterXrange = [minx prctile(cx,99.5).*1.2];
    else
    DATA.plot.clusterYrange = [miny prctile(cy,95).*2];
    DATA.plot.clusterXrange = [minx prctile(cx,95).*2];
    end
    if min(cx) < 0
        DATA.plot.clusterXrange(1) = prctile(cx,1) * 1.1;
    end
    if min(cy) < 0
        DATA.plot.clusterYrange(1) = prctile(cy,1) * 1.1;
    end
elseif DATA.plot.autoscale == 4
    DATA.plot.clusterYrange = [min(cy) max(cy)];
    DATA.plot.clusterXrange = [min(cx) max(cx)];
end


function DATA = FinishXYPlot(DATA)
    if ismember(DATA.plot.autoscale,[0 2 3 4]) % 2 means I autoscale
        DATA = SetXYRanges(DATA);
        set(gca,'Xlim',DATA.plot.clusterXrange,'Ylim',DATA.plot.clusterYrange);

    end


        
        
    xlabel(DATA.spkvarnames{DATA.plot.clusterX});
    ylabel(DATA.spkvarnames{DATA.plot.clusterY});
    it = findobj(DATA.xyfig, 'Tag','Clusterid');
    if isempty(it)
        c = 1;
    else
        c = get(it,'value');
    end
    ei = DATA.currentexpt;
    if ~isempty(DATA.explabels{DATA.currentexpt})
        expname = DATA.explabels{DATA.currentexpt};
        expname = DATA.Expts{DATA.currentexpt}.Header.expname;
     else
        expname = DATA.Expts{DATA.currentexpt}.Header.expname;
    end
    if DATA.currenttrial > 1
        expname = [expname sprintf('from %d',DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial)];
    end
    p = DATA.probe;
    str = '';
    [dp, x] = CalcDprime(DATA,1);
    if isfield(DATA,'cluster') & iscluster(DATA.cluster,c,p)
        if isfield(DATA.cluster{c,p},'dprime')
            str = sprintf('d=%.2f',DATA.cluster{c,p}.dprime);
        end
        if isfield(DATA,'AllSpikes')
            id = find(DATA.AllSpikes{DATA.probe}.codes(DATA.spklist,2) == c);
        else
        id = find(DATA.AllData.Spikes.codes(DATA.spklist,2) == c);
        end
        if isfield(DATA.Expts{ei}.Stimvals,'du')
            rate = length(id)./(length(DATA.Expts{ei}.Trials) .* DATA.Expts{ei}.Stimvals.du);
        else
            rate = length(id)./length(DATA.Expts{ei}.Trials);
        end
       
    title(sprintf('%s: C%d %sX%s %s P%d %.1fHz',expname,c,DATA.spkvarnames{DATA.cluster{c,p}.params(1)},DATA.spkvarnames{DATA.cluster{c,p}.params(2)},str,DATA.probe,rate)); 
    else
    title(sprintf('%s: Dprime %.2f',expname,dp)); 
    end
    if DATA.plot.clusterZ > 0
        Plot3DClusters(DATA,0);
    end

    
function DATA = CalcClusterVars(DATA, ispk, varargin)
SPKENERGY=1;
SPKVARE = 2;
noforce  = 1;
j = 1;
probe = DATA.probe;
while j <= length(varargin)
    if isstruct(varargin{j}) & isfield(varargin{j},'values')
        j = j+1;
        Spikes = varargin{j};
    elseif strncmpi(varargin{j},'force',5)
        noforce = 0;
    elseif strncmpi(varargin{j},'probe',5)
        j =j+1;
        probe = varargin{j};
    end
    j = j+1;
end

if isfield(DATA,'AllClusters') && noforce % All calculated and stored, but No ADCs, so
  return;
end

if length(ispk) == 2
    ispk = ispk(1):ispk(2);
end
if ispk
    DATA = CheckForPCA(DATA,ispk, 1);
    if isfield(DATA,'AllSpikes')
        Spikes = DATA.AllSpikes{probe};
        PCs = DATA.AllSpikes{probe}.codes;
    else
        Spikes = DATA.AllData.Spikes;
        if isempty(DATA.AllData.pcs)
            PCs = DATA.AllData.Spikes.codes;
        else
        PCs = DATA.AllData.pcs;
        end
    end
    adc = Spikes.values(ispk,:);
    if isfield(DATA.AllData.Spikes,'dVdt')
        energy  = sum(Spikes.dVdt(ispk,:)'.^2);
    else
    energy  = sum(diff(adc').^2);
    end
    svar = var(adc');
    DATA.Spikes.energy(ispk)= energy;
    if DATA.plot.clusterX == SPKENERGY
        DATA.Spikes.cx(ispk)= energy;
    else
        DATA.Spikes.cx(ispk)= GetSpikeVals(DATA, ispk, Spikes.values(ispk,:), Spikes.dVdt(ispk,:),DATA.plot.clusterX, 1,PCs(ispk,:));
    end
    DATA.Spikes.vw(ispk) = svar./energy;
    if DATA.plot.clusterY == SPKVARE
        DATA.Spikes.cy(ispk)= DATA.Spikes.vw(ispk);
    else
        DATA.Spikes.cy(ispk)= GetSpikeVals(DATA, ispk, Spikes.values(ispk,:), Spikes.dVdt(ispk,:),DATA.plot.clusterY, 1,PCs(ispk,:));
    end
    if isfield(DATA,'AllSpikes')
        lastspk = min([length(DATA.AllSpikes{probe}.times) length(DATA.Spikes.cx)]);
        DATA.AllSpikes{probe}.cx = DATA.Spikes.cx(1:lastspk);
        DATA.AllSpikes{probe}.cy = DATA.Spikes.cy(1:lastspk);
    end
end

function [x,y, DATA] = GetClusterSpace(DATA, Expt)


    x = DATA.plot.clusterX;
    y = DATA.plot.clusterY;
    p = DATA.probe;
    if isfield(Expt,'Cluster') & ~isempty(Expt.Cluster)
        j = 1;
        while j < size(Expt.Cluster,1) && isempty(Expt.Cluster{j,p})
            j = j+1;
        end
        if j <= size(Expt.Cluster,1) & p <= size(Expt.Cluster,2) & ~isempty(Expt.Cluster{j,p}) & isfield(Expt.Cluster{j,p},'params')
            x = Expt.Cluster{j,p}.params(1);
            y = Expt.Cluster{j,p}.params(2);
            DATA.clusterArange = Expt.Cluster{j,p}.Arange;
            DATA.clusterBrange = Expt.Cluster{j,p}.Brange;
            if isfield(Expt.Cluster{j,p},'Erange')
                DATA.clusterErange = Expt.Cluster{j,p}.Erange;
            end
        end
    else
        x = DATA.plot.clusterX;
        y = DATA.plot.clusterY;
    end
       
function set = ClusterIsSet(Expt, probe)

if isfield(Expt,'Cluster') && size(Expt.Cluster,2) >= probe && ...
        ~isempty(Expt.Cluster{1,probe}) && isfield(Expt.Cluster{1,probe},'x')
    set = 1;
    if isfield(Expt.Cluster{1,probe},'autocut') & Expt.Cluster{1,probe}.autocut > 0
        set = 2;
    end
else
    set = 0;
end


function DATA = Plot3DClusters(a,b);
recalc = 0;
    if isfield(a,'AllData');
        DATA = a;
        recalc = b;
    else
    DATA = GetDataFromFig(a);
    end


[xyzfig, isnew] = GetFigure('Cluster3D');
if isnew
    if ~isfield(DATA.plot,'clusterZ')
        DATA.plot.clusterZ = 5;
    end
    z.parentfigtag = DATA.tag.top;
    z.parentfig = DATA.toplevel;
    set(xyzfig,'UserData',z);
%   uicontrol(gcf,'Style', 'pop','String',DATA.spkvarnames,'Position', bp,...
%      'Tag','SetClusterZ','Callback',@SetClusterZ,'value',DATA.plot.clusterZ);
    xyfig = xyzfig;
  hm = uimenu(xyfig,'Label','X','Tag','XClusterMenu');
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@SetXYCluster, 1, 3, k});
    end
    hm = uimenu(xyfig,'Label','Y','Tag','YClusterMenu');
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@SetXYCluster, 2,3, k});
    end
    hm = uimenu(xyfig,'Label','Z','Tag','3DplotMenu');
     uimenu(hm,'Label','None','Callback',{@SetXYCluster, 3, 0});
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@SetXYCluster, 3, 3, k});
    end

  
  DATA.figs.xyz = xyzfig;
    end
Spks = DATA.AllData.Spikes;
ispk = DATA.spklist;
DATA = CheckForPCA(DATA, ispk, 0);
if isempty(DATA.AllData.pcs)
    PCs = [];
else
PCs = DATA.AllData.pcs(ispk,:);
end
classify = 1;
[cz, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterZ, classify,PCs);
DATA.Spikes.cz(ispk) = cz;
if recalc
[cx, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterX, classify,PCs);
DATA.Spikes.cx(ispk) = cx;
[cy, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterY, classify,PCs);
DATA.Spikes.cy(ispk) = cy;
end
cl = unique(Spks.codes(ispk,2));
hold off;
for j = 1:length(cl)
    id = find(Spks.codes(ispk,2) == cl(j));
plot3(DATA.Spikes.cx(ispk(id)),DATA.Spikes.cy(ispk(id)),cz(id),'.',...
    'color',DATA.spkcolor{cl(j)+1},'markersize',DATA.ptsize);
     hold on;
end
xlabel(['X:' DATA.spkvarnames{DATA.plot.clusterX}]);
ylabel(['Y:' DATA.spkvarnames{DATA.plot.clusterY}]);
zlabel(['Z:' DATA.spkvarnames{DATA.plot.clusterZ}]);
if DATA.plot.autoscale == 0
    set(gca,'xlim',DATA.plot.clusterXrange,'ylim',DATA.plot.clusterYrange);
end

function res = isExptCluster(E,c,p)
    if ~isfield(E,'Cluster');
        res = 0;
    else
        res = iscluster(E.Cluster,c,p);
    end
    

    function res = iscluster(C,c,p)
    if isempty(C) | size(C,1) < c | size(C,2) < p | isempty(C{c,p})
        res = 0;
    elseif ~isfield(C{c,p},'params') | ~isfield(C{c,p},'x')
        res = 0;
    elseif ~isfield(C{c,p},'nspk') || C{c,p}.nspk == 0
        res = 2;
    else
        res = 1;
    end

    
function [x,y,z] = CalcDensity(DATA, expspks, mode)

if ismember(DATA.syncsign,[-1 1])
energy = DATA.Spikes.cx(DATA.sids{1});
vw = DATA.Spikes.cy(DATA.sids{1});
elseif isfield(DATA,'AllClusters')
    energy  = DATA.AllClusters(DATA.probe).cx(expspks);
    vw  = DATA.AllClusters(DATA.probe).cy(expspks);
else
energy = DATA.Spikes.cx(expspks);
vw = DATA.Spikes.cy(expspks);
end
DATA.plot.DensitySigma = [3 3];

if length(vw) > 10000
    lprc = 0.01;
    hprc = 99.99;
elseif length(vw) > 10000
    lprc = 0.1;
    hprc = 99.9;
    
elseif length(vw) > 1000
    lprc = 1;
    hprc = 99;
    DATA.plot.DensitySigma = [5 5];
else
    lprc = 5;
    hprc = 95;
    DATA.plot.DensitySigma = [10 10];
end    

if DATA.xyfig == gcf && mode < 3
erange = get(gca,'Xlim');
vrange = get(gca,'Ylim');
else
    if prctile(energy,hprc) > prctile(energy,hprc-1) * 3
        erange = [prctile(energy,lprc) prctile(energy,hprc-1)];
    else
        erange = [prctile(energy,lprc) prctile(energy,hprc)];
    end
vrange = [prctile(vw,lprc) prctile(vw,hprc)];
end
%GetFigure('DensityPlot');
hold off;
nbins = 200;
[x,y] = meshgrid(linspace(erange(1),erange(2),nbins),linspace(vrange(1),vrange(2),nbins));
tic;
if mode ==1 % add real gaussian to grid for each
    z = zeros(size(x));
    sx = (diff(erange)/100)^2;
    sy = (diff(vrange)/100)^2;
    for j=1:length(energy)
        z = z + exp(-(x-energy(j)).^2/sx - (y-vw(j)).^2/sy);
    end
elseif mode ==2 || mode == 3 %build fine 2-D histogram, then smooth
    [gx,gy] = meshgrid(-10:10,-10:10);
    sx= DATA.plot.DensitySigma(1);
    sy= DATA.plot.DensitySigma(2);
    G = exp(-(gx).^2/sx - (gy).^2/sy);
    G = G./sum(G(:));
    z = zeros(size(x));
% ignore spikes where both are set to 0
    idx = find(vw ~=0 | energy ~=0);
    
    vi = 1+floor(nbins * (vw(idx)-vrange(1))/diff(vrange));
    ei = 1+floor(nbins * (energy(idx)-erange(1))/diff(erange));
    idx = find(ei > 0 & ei <= nbins & vi > 0 & vi <= nbins);
    for j =idx
        z(vi(j),ei(j)) = z(vi(j),ei(j))+1;
    end
    z = conv2(z,G,'same');
end

function Update(a,b)
        
function DensityPlot(a,b)
global mousept;

%DATA = combine('getstate');
DATA = GetDataFromFig(a);

if isfield(DATA,'spklist') & length(DATA.spklist) > 10
    expspks = DATA.spklist;
else
    expspks = DATA.spkrange(1):DATA.spkrange(2);
end

    
    figure(DATA.xyfig);
if DATA.densityplot
    DATA.densityplot = 0;
    hold off;
    DATA = DrawXYPlot(DATA,expspks);
%    ClearMouse;
    set(DATA.toplevel,'UserData',DATA);
    return;
end

[x,y,z] = CalcDensity(DATA, expspks, 2);

erange = get(gca,'Xlim');
vrange = get(gca,'Ylim');
%GetFigure('DensityPlot');
hold off;
toc
pcolor(x,y,z);
shading('interp')
hold on;
if 0 %code for trackin peak ridge. May help track clusters/auto cut
[a,b] = max(z);
py = smooth(y(b,1),3,'gauss');
plot(x(1,:), py);
base = 0.1;
for j = length(py):-1:1
    dip(j) = a(j)./(base+max(a(j:end)).^2);
end
plot(dip);
end
%DATA = DrawClusters(DATA, DATA.cluster, 0);
DATA.densityplot = 1;
it = findobj(DATA.xyfig,'Tag','ClusterZmax');
x = caxis;
set(it,'string',sprintf('%.2f',x(2)));
set(DATA.toplevel,'UserData',DATA);


