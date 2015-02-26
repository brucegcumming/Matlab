function DATA = PlaySpikes(DATA, expid, varargin)

j = 1;
xyonly = 0;
plotxy = 1;
onetrial = 0;
setfig = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'onetrial',6)
        onetrial = 1;
    elseif strncmpi(varargin{j},'usefig',6)
        setfig = 0;
    elseif strncmpi(varargin{j},'quickspks',6)
        onetrial = 2;
        plotxy = 0;
    elseif strncmpi(varargin{j},'showxy',6)
        xyonly = 0;
    elseif strncmpi(varargin{j},'xyonly',6)
        xyonly = 1;
    end
    j = j+1;
end

if DATA.state.usexycache
    xyonly = 1;
end
mode = 2;
cw = DATA.plot.cw;
ch = DATA.plot.ch;
rh = ch+10;
if setfig
    [~, xyfig] = cmb.SetFigure(DATA.tag.clusterxy,DATA);
    isnew = 0;
%    [xyfig, isnew] = GetFigure(DATA.tag.clusterxy);
    DATA.xyfig = xyfig;
else
    isnew = 0;
    xyfig = gcf;
end

if isnew
    nf = strmatch(DATA.tag.clusterxy,{DATA.figpos.tag});
    if isempty(nf)
        bp = get(DATA.toplevel,'Position');
        fp = get(xyfig,'Position');
        nf = length(DATA.figpos)+1;
        DATA.figpos(nf).tag = 'xyfig';
        DATA.figpos(nf).pos =  [bp(1)+bp(3) bp(2)+bp(4)-fp(4) fp(3) fp(4)];
        if fp(4) == 0 %can happen with double clicks
            fp(4) = 420;
        end
        if fp(3) == 0
            fp(3) = 420;
        end
        
        if sum(DATA.figpos(nf).pos([2 4]))+100 > DATA.gui.scrsz(4)
            DATA.figpos(nf).pos(2) = DATA.gui.scrsz(4) - fp(4) - 100;
        end
        DATA.figpos(nf).tag = DATA.tag.clusterxy;
    end
    set(xyfig,'Position',DATA.figpos(nf).pos);
    setappdata(DATA.toplevel,'FigPos',DATA.figpos);
    bp = [5 5 40 20];
    cp = bp;
    uicontrol(xyfig,'style','pop','string','1|2|3|4|5|6|7|Artifact','Position',bp,'Tag','Clusterid',...
        'Callback',@cmb.Update);
    cp(3) = cw * 5;
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'Style', 'checkbox', 'Callback', @cmb.DensityPlot, ...
        'String', 'Dens','Tag','Density', 'Position', cp);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', {@cmb.NextList, 'setfirst'}, ...
        'String', 'Set+Next', 'Position', cp,'tag','Set+Next');
    cp(2) = cp(2)+rh;
    
    
    if length(DATA.probelist) > 1
        cp(3) = cw*1.5;
        uicontrol(xyfig,'style','pushbutton','string','+','Position',cp, 'Callback',@cmb.AddOneCellToList);
        
        cp(1) = cp(1) + cp(3);
        cp(3) = cw*3;
        uicontrol(xyfig,'style','pop','string','1|2|3|4|5|6|7|8|9|10','Position',cp,'Tag','AddOneCellToList',...
            'Callback',@cmb.AddOneCellToList);
        cp(2) = cp(2)+rh;
        cp(1) = 5;
        cp(3) = 40;
    end
    
    
    
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', @cmb.NextList, ...
        'String', 'Next', 'Position', cp);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', {@cmb.NextList, 'clearfirst'}, ...
        'String', 'Clr+Next', 'Position', cp);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'style','pushbutton','string','spool','Position',cp,'Tag','SpoolSpikes',...
        'Callback', @cmb.SpoolSpikes);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'style','pushbutton','string','Optim','Position',cp,'Tag','OptimizeDprime',...
        'Callback', @cmb.OptimizeDprimeHit);
    cp(2) = cp(2)+rh;
    
    bp(1) = bp(1) + bp(3) + 10;
    uicontrol(xyfig,'style','pop','string','Not Set|Nothing|MU-|MU+|Poor|OK|Good|V Good|Excellent|FromDprime','Position',bp,'Tag','ClusterQuality',...
        'Callback',@cmb.Update);
    
    bp(1) = bp(1) + bp(3) + 10;
    uicontrol(xyfig,'style','pushbutton','string','Set','Position',bp,'Callback', @cmb.SetExptClusters);
    bp(1) = bp(1) + bp(3) + 10;
    bp(3)=cw;
    uicontrol(xyfig,'style','CheckBox','Position',bp,'Tag','ClusterIsSet');
    bp(1) = bp(1) + bp(3) + 10;
    bp(3)=cw*3;
    uicontrol(xyfig,'style','pushbutton','string','Del','Position',bp,'Callback', @cmb.DelClusterButton);
    bp(1) = bp(1) + bp(3) + 10;
    uicontrol(xyfig,'style','pushbutton','string','Clr','Position',bp,'Callback', @cmb.ClrSpkWin);
    bp(1) = bp(1) + bp(3) + 10;
    bp(3)=cw*5;
    uicontrol(xyfig,'style','text','string','Max: X','Position',bp);
    bp(1) = bp(1) + bp(3);
    uicontrol(xyfig,'style','edit','string',sprintf('%.2f',DATA.plot.clusterXrange(2)),'Position',bp,...
        'Callback', @cmb.RescaleClusterPlot,'Tag','ClusterXmax');
    
    bp(1) = bp(1) + bp(3) + 10;
    bp(3) = cw*1;
    uicontrol(xyfig,'style','text','string','Y','Position',bp);
    bp(1) = bp(1) + bp(3);
    bp(3) = cw*4;
    uicontrol(xyfig,'style','edit','string',sprintf('%.2f',DATA.plot.clusterYrange(2)),'Position',bp,...
        'Callback', @cmb.RescaleClusterPlot,'Tag','ClusterYmax');
    bp(1) = bp(1) + bp(3) + 10;
    bp(3) = cw*1;
    uicontrol(xyfig,'style','text','string','Z','Position',bp);
    bp(1) = bp(1) + bp(3);
    bp(3) = cw*4;
    uicontrol(xyfig,'style','edit','string',sprintf('%.2f',DATA.plot.clusterZrange(2)),'Position',bp,...
        'Callback', @cmb.RescaleClusterPlot,'Tag','ClusterZmax');
    bp(1) = bp(1) + bp(3);
    bp(3) = cw*5;
    uicontrol(xyfig,'style','CheckBox','string','auto','Position',bp,'Tag','AutoScale',...
        'value',(DATA.plot.autoscale > 0),'Callback',@cmb.Update);
    bp(1) = bp(1) + bp(3) + 10;
    hm = uimenu(xyfig,'Label','Plot3D','Tag','3DplotMenu');
    uimenu(hm,'Label','None','Callback',{@cmb.SetXYCluster, 3, 2, 0});
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@cmb.SetXYCluster, 3, 2,k});
    end
    hm = uimenu(xyfig,'Label','X','Tag','XClusterMenu');
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@cmb.SetXYCluster, 1, 2,k});
    end
    hm = uimenu(xyfig,'Label','Y','Tag','YClusterMenu');
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@cmb.SetXYCluster, 2, 2, k});
    end
    hm = uimenu(xyfig,'Label','Ops','Tag','XYAddMenu');
    uimenu(hm,'Label','Force Cluster Ids','Callback',{@cmb.AddClusterIdButton});
    a = strmatch('Energy',DATA.spkvarnames,'exact');
    b = strmatch('Var/Energy',DATA.spkvarnames,'exact');
    uimenu(hm,'Label','Energy-Var/E','Callback',{@cmb.SetXYCluster, 4, 2, a(1), b(1)});
    b = strmatch('sqrt(Var/Energy)',DATA.spkvarnames,'exact');
    a = strmatch('sqrt(Energy)',DATA.spkvarnames,'exact');
    uimenu(hm,'Label','sqrt(Energy)-sqrt(Var/E)','Callback',{@cmb.SetXYCluster, 4, 2, a(1), b(1)});
    b = strmatch('Var/Energy',DATA.spkvarnames,'exact');
    a = strmatch('sqrt(Energy)',DATA.spkvarnames,'exact');
    uimenu(hm,'Label','sqrt(Energy)-Var/E','Callback',{@cmb.SetXYCluster, 4, 2, a(1), b(1)});
    b = strmatch('PCA1',DATA.spkvarnames,'exact');
    a = strmatch('PCA2',DATA.spkvarnames,'exact');
    uimenu(hm,'Label','PCA1 - 2','Callback',{@cmb.SetXYCluster, 4, 2, a(1), b(1)});
    b = strmatch('PCA2',DATA.spkvarnames,'exact');
    a = strmatch('PCA3',DATA.spkvarnames,'exact');
    uimenu(hm,'Label','PCA2 - 3','Callback',{@cmb.SetXYCluster, 4, 2, a(1), b(1)});
    if DATA.subprobes > 1
        b = strmatch('Energy 1',DATA.spkvarnames,'exact');
        a = strmatch('Energy 2',DATA.spkvarnames,'exact');
        c = strmatch('Energy 3',DATA.spkvarnames,'exact');
        uimenu(hm,'Label','Energy1-2-3','Callback',{@cmb.SetXYCluster, 5, 3, a(1), b(1), c(1)});
        uimenu(hm,'Label','Pts trode 1','Callback',{@cmb.SetXYCluster, 6,1});
        uimenu(hm,'Label','Pts trode 2','Callback',{@cmb.SetXYCluster, 6,2});
        uimenu(hm,'Label','Pts trode 3','Callback',{@cmb.SetXYCluster, 6,3});
        uimenu(hm,'Label','Pts trode 4','Callback',{@cmb.SetXYCluster, 6,4});
    end
    uimenu(hm,'Label','recalc PCA','Callback',{@cmb.SpkVMenu, 3});
    uimenu(hm,'Label','Gmix Distance','Callback',{@cmb.SpkVMenu, 2});
    tmpdat.parentfigtag = DATA.tag.top;
    set(xyfig,'UserData',tmpdat);
end

set(xyfig, 'KeyPressFcn',@cmb.KeyPressed);
set(xyfig, 'KeyReleaseFcn',@cmb.KeyReleased);
set(xyfig, 'WindowButtonDownFcn',@cmb.ButtonPressed);
set(xyfig, 'WindowButtonMotionFcn',@cmb.ButtonDragged);
set(xyfig, 'WindowButtonUpFcn',@cmb.ButtonReleased);
set(xyfig, 'WindowScrollWheelFcn',@cmb.ScrollWheel);
if xyonly
    cmb.ClearMouse;
    return;
end
hold off;
cmb.CheckSpoolButton(DATA);
%if isempty(findobj('Tag',DATA.tag.spikev
%   sfig = figure('Renderer','painters','Tag','DATA.tag.spikev');
[sfig, isnew] = GetFigure(DATA.tag.spikev);
if isnew
    nf = strmatch(DATA.tag.spikev,{DATA.figpos.tag});
    if  isempty(nf)
        nf = length(DATA.figpos)+1
        bp = get(DATA.toplevel,'Position');
        fp = get(sfig,'Position');
        DATA.figpos(nf).pos =  [bp(1) bp(2)-fp(4)-80 fp(3) fp(4)];
        DATA.figpos(nf).tag = DATA.tag.spikev;
    end
    set(sfig,'Position',DATA.figpos(nf).pos);
    setappdata(DATA.toplevel,'FigPos',DATA.figpos);
    x = 10;
    c = 10;
    bp = [10 10 40 20];
    uicontrol(sfig,'style','pushbutton','string','>>','Position',bp,'Tag','NextTrial',...
        'Callback', @cmb.PlayNextTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pushbutton','string','<<','Position',bp,'Tag','LastTrial',...
        'Callback', @cmb.PlayLastTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pushbutton','string','X','Position',bp,'Tag','CutTrial',...
        'Callback', @cmb.CutTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pop','string','1|2|3','Position',bp,'Tag','ChooseTrial',...
        'Callback', @cmb.SelectTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pushbutton','string','>>','Position',bp,'Tag','NextSpike',...
        'Callback', @cmb.PlayNextSpike);
    bp(1) = bp(1)+bp(3);
    bp(3)=c*6;
    uicontrol(sfig,'style','pushbutton','string','spool','Position',bp,'Tag','SpoolSpikes',...
        'Callback', @cmb.SpoolSpikes);
    bp(1) = bp(1) + bp(3);
    
    bp(3)=c*9;
    uicontrol(gcf,'Style', 'checkbox',...
        'String', 'dVdt', 'Tag', 'dVdt', 'Position', bp,'value',DATA.plot.dvdt,...
        'Callback',@cmb.Update);
    bp(1) = bp(1) + bp(3);
    
    bp(3)=c*9;
    uicontrol(gcf,'Style', 'checkbox',...
        'String', 'Stop', 'Tag', 'StopSpool', 'Position', bp,'value',0,...
        'Callback',@cmb.Update);
    tmpdat.parentfigtag = DATA.tag.top;
    set(sfig,'UserData',tmpdat);
    set(sfig, 'WindowScrollWheelFcn',@cmb.ScrollTrial);
    
    hm = uimenu(sfig,'Label','Set','Tag','SpikeVMenu');
    uimenu(hm,'Label','All Trials','Callback',{@cmb.SetTrialRange, 3});
    uimenu(hm,'Label','End Range','Callback',{@cmb.SetTrialRange, 1});
    uimenu(hm,'Label','Range->','Callback',{@cmb.SetTrialRange, 2});
    uimenu(hm,'Label','Vrange','Callback',{@cmb.SpkVMenu, 1});
    
    
    DATA.hline = 0;
end
it = findobj(sfig,'Tag','ChooseTrial');
set(it,'string',sprintf('%d|',[DATA.Expts{expid}.Trials.Trial]),'value',1);
DATA.svfig = sfig;
if DATA.spooling ~= 2 %2 == spool from current trial
    DATA.currenttrial = 0;
end
DATA.ISIpair = 0;
DATA.currentexpt(1) = expid;
if DATA.test.fastplot
    set(DATA.svfig,'renderer','painters');
    ax = findobj(DATA.svfig,'type','ax');
    set(ax,'DrawMode','fast','XTickMode','manual','YTickMode','manual');
    set(ax,'Xlim',[1 46],'Ylim',[-5 5]);
end
ax = findobj(DATA.svfig,'type','ax');



if DATA.state.uselfp
    DATA.state.lfig = GetFigure('LFP');
    hold off;
    DATA.Expts{expid} = LoadSpike2LFP(DATA.Expts{expid});
end
colors = mycolors;
if DATA.probe == 100
    Spks = DATA.AllData.UstimV;
elseif isfield(DATA,'AllClusters')
    if isfield(DATA.AllData,'Spikes') && DATA.AllData.Spikes.probe == DATA.probe
        Spks = DATA.AllData.Spikes;
    elseif (isdir(DATA.name) || DATA.state.somespikes)
        %remove monkey name then add again in case its not in prefix
        prefix = strrep(DATA.prefix,DATA.monkey,'');
        prefix = [DATA.monkey prefix];  %new style for spks files
        pid = GetProbe(DATA, DATA.currentexpt(1), DATA.probe);
        spkfile = sprintf('%s/Spikes/%s.p%dt%d.mat',DATA.datadir,prefix,pid,DATA.currentexpt(1));
        DATA.AllData.Spikes = cmb.GetProbeSpikes([],spkfile,'Spikes',pid);
        Spks = DATA.AllData.Spikes;
    end
elseif isfield(DATA,'AllSpikes')
    Spks = DATA.AllSpikes{DATA.probe};
else
    Spks = DATA.AllData.Spikes;
end
tstart = now;
if isempty(Spks)
    Spks.codes = 0;
    Spks.times = 0;
end
if size(Spks.codes,2) > 1
    ctype = 2;
else
    ctype = 1;
end
if isfield(Spks,'values')
    nsamples = size(Spks.values,2).* size(Spks.values,3);
else
    nsamples = 0;
end
if size(Spks.codes,ctype) > 1
    nclusters = 1+max(Spks.codes(:,ctype));
else
    nclusters = 1;
end
set(0,'CurrentFigure',xyfig);

if isfield(DATA.plot,'useprobe')
    xprobes = find(DATA.plot.useprobe);
    xprobes = setdiff(xprobes,DATA.probe);
else
    xprobes = [];
end
DATA.xprobes = xprobes;
probes = [DATA.probe xprobes];
nprobes = length(probes);



reclassify = 0;
cid = findobj('Tag','ClusterIsSet');
if DATA.spooling
    DATA.state.recut = 1;
elseif ClusterIsSet(DATA.Expts{expid},DATA.probe) & DATA.state.recut
    DATA.cluster = cmb.CopyClusters(DATA.cluster,DATA.Expts{expid}.Cluster);
    DATA.state.recut = 1;
    if ismember(DATA.Expts{expid}.gui.clustertype, [0 2 3]) %% The online cut or autocut
        set(cid,'value',0);
    else
        set(cid,'value',1);
    end
    if ~isfield(DATA.Expts{expid}.gui,'classified') | DATA.Expts{expid}.gui.classified ~= 1
        DATA = SetExptSpikes(DATA,expid,0);
    end
    if length(xprobes) & isfield(DATA,'AllSpikes')
        for j = 1:length(probes)
            if ~isfield(DATA.AllSpikes{probes(j)},'spklist')
                times(1) = DATA.Expts{expid}.Trials(1).Start(1)-500;
                times(2) = DATA.Expts{expid}.Trials(end).End(end)+500;
                DATA.AllSpikes{probes(j)}.spklist = FindSpikes(DATA, times, probes(j),[]);
            end
            DATA = SetSpkCodes(DATA, DATA.AllSpikes{probes(j)}.spklist, probes(j),0);
        end
    end
    
elseif ~isfield(DATA,'cluster')
    DATA.cluster = {};
elseif DATA.state.recut %No cluster yet defined
    DATA.state.recut = 2;
    % when inheriting a cluster, don't inherit any spike ranges
    nclusters = 1+max(Spks.codes(:,ctype));
    p = DATA.probe;
    for j = 1:nclusters
        if j <= size(DATA.cluster,1) & size(DATA.cluster,2) >= p & ~isempty(DATA.cluster{j,p})
            if isfield(DATA.Expts{expid}.gui,'spkrange') & ~isempty(DATA.Expts{expid}.gui.spkrange)
                DATA.cluster{j,p}.firstspk = DATA.Expts{expid}.gui.spkrange(1);
                DATA.cluster{j,p}.lastspk = DATA.Expts{expid}.gui.spkrange(2);
            else
                DATA.cluster{j,p}.firstspk = 1;
                DATA.cluster{j,p}.lastspk = length(Spks.times);
            end
        end
    end
    set(cid,'value',0);
end
nclusters = double(nclusters);
if plotxy ==0 %only replot waveforms
elseif DATA.plot.synccluster == 0
    cmb.ClrSpkWin(DATA);
    DATA = cmb.DrawClusters(DATA,DATA.cluster,0);
    if isfield(DATA.Expts{expid},'OnlineCluster') & DATA.state.showonlineclusters
        DATA = cmb.DrawClusters(DATA,DATA.Expts{expid}.OnlineCluster,0);
    end
else
    hold off;
    plot(0,0,'+');
    hold on;
    DATA = cmb.DrawClusters(DATA,DATA.cluster,0);
end
if DATA.plot.autoscale == 0
    set(gca,'Xlim',DATA.plot.clusterXrange,'Ylim',DATA.plot.clusterYrange);
end

set(0,'CurrentFigure',sfig);
hold off;

if DATA.spikelist == -1
    spikelist = [0 1 2 3 4];
else
    spikelist = DATA.spikelist;
end



if isfield(DATA,'AllSpikes')
    for k = 1:length(probes)
        nclusters = max([nclusters max(DATA.AllSpikes{probes(k)}.codes(:,2))]);
    end
end
n = cmb.CountClusters(DATA.cluster);
nclusters = max([n nclusters]);
x = 32; %should calculate real size..
for k = 0:length(xprobes)
    ids{k+1} = FindSpikes(DATA, DATA.Expts{expid}.Header.trange,probes(k+1),[]);
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
DATA.svhn = nclusters+1;

expspks = ids{1};
if ~isfield(DATA.Expts{expid}.gui,'spks') %haven't loaded this expt yet
    times = [DATA.Expts{expid}.Trials(1).Start(1) DATA.Expts{expid}.Trials(end).End(end)];
    if isfield(DATA,'AllSpikes')
        nspk = sum(DATA.AllSpikes{DATA.probe}.times > times(1) & DATA.AllSpikes{DATA.probe}.times < times(2));
    else
        nspk = sum(DATA.AllData.Spikes.times > times(1) & DATA.AllData.Spikes.times < times(2));
    end
    if nspk < 10000
        DATA.ptsize = 6;
    elseif nspk < 2000
        DATA.ptsize = 10;
    else
        DATA.ptsize = 4;
    end
    expspks = ids{1};
else
    %    expspks = DATA.Expts{expid}.gui.spks; % no good if > 1 probe loaded
end

if DATA.plot.autoVrange
    DATA.spklist = expspks;
    DATA = cmb.SpkVMenu(DATA,0,1);
end
DATA.ptsize = cmb.CheckPtSize(DATA, expspks);

%FindSpikes might list some spikes just past the list recorded in
%Expts.gui.spkrange.
if isfield(DATA,'Spikes') & isfield(DATA.Spikes,'cx')
    if max(ids{1}) > length(DATA.Spikes.cx)
        ispk = length(DATA.Spikes.cx):max(ids{1});
        DATA = CalcClusterVars(DATA,  ispk);
    end
end

if length(ids{1}) & isfield(DATA.AllData.Spikes,'codes') & size(DATA.AllData.Spikes.codes,1) >= ids{1}(end)
    id = find(DATA.AllData.Spikes.codes(ids{1},ctype) > nclusters);
    if length(id)
        DATA.AllData.Spikes.codes(ids{1}(id),ctype) = 0;
    end
end

DATA = CheckForPCA(DATA, expspks, 0);


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
    if isfield(Spks,'maxv') && Spks.maxv < DATA.plot.SpikeMaxV
        set(gca,'Xlim',[1 nsamples],'Ylim',[-Spks.maxv Spks.maxv]);
    else
        set(gca,'Xlim',[1 nsamples],'Ylim',[-DATA.plot.SpikeMaxV DATA.plot.SpikeMaxV]);
    end
end
DATA.nclusters = cmb.CountClusters(DATA.cluster);
nt = 1;
allspks = [];
firstspk = 0;
hline = 0;
if ~isfield(DATA,hline)
    DATA.hline = 0;
end

DATA.playingspk = 1;
DATA.nclusters = cmb.CountClusters(DATA.cluster);
if DATA.nclusters > nclusters
    nclusters = DATA.nclusters;
end
DATA.svhn = nclusters+1;

set(DATA.toplevel,'UserData',DATA);
if isempty(DATA.currentcluster) || length(DATA.cluster) < DATA.currentcluster
    DATA.currentcluster = 1;
end

Aargs = {}; timemode = 0;
trialdur = DATA.Expts{expid}.Trials(1).End(end) - DATA.Expts{expid}.Trials(1).Start(1);
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


nprobes = length(probes);
step = DATA.plot.SpikeMaxV;
DATA.syncprobes = [];
if length(probes) > 1 && isfield(DATA,'AllSpikes')
    if DATA.plot.synccluster > 0
        DATA.Spikes.cx(1:end) = 0;
        DATA.Spikes.cy(1:end) = 0;
        DATA.syncprobes = probes(1:2);
    else
    end
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
            [aid,bid] = cmb.FindSync(DATA.AllSpikes{probes(1)}.times(ids{1}),...
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
            [E, V, DATA.pca] = cmb.SpikePCA(DATA.AllSpikes,[probes(1) probes(j)],{aid bid} );
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
start = max([DATA.currenttrial 1]);
last = length(DATA.Expts{expid}.Trials);
% negative Trial numbers indicate manual exclusion
uset = start-1+find([DATA.Expts{expid}.Trials(start:last).Trial] > 0);
if onetrial
    uset = uset(1);
end
nsync = 0;nspks = 0;
tc = 1;
syncspikes = [];

while tc <= length(uset)  %not for loop, so can change tc inside loop
    trial = DATA.Expts{expid}.Trials(uset(tc)).Trial;
    stop = get(findobj(DATA.svfig,'Tag','StopSpool'),'value');
    if stop
        set(findobj(DATA.svfig,'Tag','StopSpool'),'value',0);
        if 0
            DATA.playingspk = 0;
            return;
        else
            tc = length(uset);
            trial = DATA.Expts{expid}.Trials(uset(nt)).Trial;
        end
    end
    if DATA.state.uselfp
        GetFigure('TrialLFP');
        cmb.PlotLFPRaw(DATA.state,DATA.Expts{expid}.Trials(uset(nt)),DATA.Expts{expid}.Header.LFPsamplerate);
    end
    if timemode
        set(0,'CurrentFigure',DATA.timefig);
    else
        set(0,'CurrentFigure',sfig);
    end
    hold off;
    itrial = find(DATA.AllData.Trialids == trial);
    if onetrial == 2
        step = round(length(DATA.spklist)/500);
        Aargs = {Aargs{:} 'spkid', DATA.spklist([1:step:length(DATA.spklist)]), 'noxy'};
    end
    if plotxy == 0;
        Aargs = {Aargs{:} 'noxy'};
    end
    if mode == 1
        DATA = cmb.PlotTrialSpikes(DATA,itrial,colors, clusters);
    elseif mode == 2
        times(1) = DATA.Expts{expid}.Trials(uset(nt)).Start(1)-DATA.state.preperiod;
        times(2) = DATA.Expts{expid}.Trials(uset(nt)).End(end)+DATA.state.preperiod;
        Trial = DATA.Expts{expid}.Trials(uset(nt));
        if isfield(Trial,'uStim') && Trial.uStim > 0
            ustim = 1;
        else
            ustim = 0;
        end
        Trial.ed = GetEval(DATA.Expts{expid},'ed',DATA.currenttrial);
        
        if DATA.syncsign < 2 & length(probes) > 1
            [spks, sspks, cx, cy] = cmb.PlotTrialSyncSpikes(DATA, times, [DATA.probe xprobes], colors,'Trial',Trial);
            if DATA.plot.synccluster > 0
                syncspikes = cat(1,syncspikes,sspks);;
                DATA.Spikes.cx(sspks(:,1)) = cx;
                DATA.Spikes.cy(sspks(:,1)) = cy;
                DATA.AllSpikes{DATA.probe}.cx(sspks(:,1)) = cx;
                DATA.AllSpikes{DATA.probe}.cy(sspks(:,1)) = cy;
                DATA.AllSpikes{xprobes(1)}.cx(sspks(:,2)) = cx;
                DATA.AllSpikes{xprobes(1)}.cy(sspks(:,2)) = cy;
            end
            nspks = nspks+length(spks);
            nsync = nsync+size(sspks,2);
            if DATA.plot.synccluster == 10
                [DATA, spks] = cmb.APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probes', [probes(1:2)]);
            end
            
        else
            
            if length(probes) > 1
                if DATA.plot.syncoverlay
                    for j= 1:length(xprobes)
                        cmb.APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe',[xprobes(j) 1+j length(probes)],'lineoff',j*(nclusters+1));
                    end
                    [DATA, spks] = cmb.APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe', [probes(1) 1 length(probes)]);
                end
                if timemode
                    for j= 1:length(xprobes)
                        cmb.APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe',[xprobes(j) 1+j length(probes)],'lineoff',j*(nclusters+1),'timemode');
                    end
                    [DATA, spks] = cmb.APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe', [probes(1) 1 length(probes)],'timemode');
                end
            else %single electrode
                if ustim == 0 || (probes(1) > 99 && ustim == 1)
                    [DATA, spks] = cmb.APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe', [probes(1) 1 length(probes)],Aargs{:});
                else
                    spks = [];
                end
                if timemode
                    cmb.APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe', [probes(1) 1 length(probes)],Aargs{:},'timemode');
                end
            end
        end
        allspks = [allspks spks'];
        nt = nt+1;
    end
    hold on;
    tc = tc+1;
end

if DATA.plot.synccluster > 0
    DATA.syncspikes = syncspikes;
end
if DATA.plot.showsync
    fprintf('%d/%d (%.3f) spikes Synchronized\n',nsync,nspks,nsync./nspks);
end


if DATA.probe == 100;
    DATA.Expts{expid}.MeanPulse = mean(Spks.values(allspks,:));
    DATA.Expts{expid}.PulseSD = std(Spks.values(allspks,:));
    plot(DATA.Expts{expid}.MeanPulse);
end
if DATA.state.uselfp  % look for calibration spikes
    
    GetFigure('LFP');
    hold off;
    CalcLFPPulse(DATA.Expts{expid},DATA.AllData,'plot');
    GetFigure('SpikeV');
end
if onetrial == 2 %Wehn calling quickspks, just paint the wwavefoems
    DATA.spklist = expspks;
    if plotxy
        DATA = cmb.DrawXYPlot(DATA,expspks);
    end
    DATA.spkrange = minmax(DATA.spklist);
elseif isempty(allspks)
    DATA.spkrange(1) = 1;
    DATA.spkrange(2) = 1;
else
    DATA.spkrange(1) = min(allspks);
    DATA.spkrange(2) = max(allspks);
    DATA.spklist = allspks;
end
DATA.Expts{expid}.gui.spkrange = DATA.spkrange;
DATA.Expts{expid}.gui.spks = allspks;
DATA.Expts{expid}.gui.s2clusters = 1+max(Spks.codes(allspks,1));
DATA.s2clusters = DATA.Expts{expid}.gui.s2clusters;
if plotxy
    DATA = cmb.FinishXYPlot(DATA);
end

axis('manual');
cmb.ClearMouse;
DATA.densityplot = 0;
xrange = get(gca,'Xlim');
yrange = get(gca,'Ylim');


if DATA.test.fastplot
    set(DATA.xyfig,'renderer','painters');
    ax = findobj(DATA.xyfig,'type','ax');
    set(ax,'DrawMode','fast','XTickMode','manual','YTickMode','manual');
    %set(,'DrawMode','fast');
end
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
if isfield(DATA,'AllClusters') && ~strncmp(DATA.filetype,'Grid',4)
    %    set(0,'currentfig',DATA.xyfig);
    %    plot(DATA.Spikes.cx(allspks),DATA.Spikes.cy(allspks),'.','markersize',DATA.ptsize);
    %    DATA.allexp = DATA.currentexpt(1);
    DATA = cmb.PlotAllProbeXY(DATA,0);
end
if isfield(DATA,'AllSpikes') && nprobes == 2 && DATA.plot.xcorr
    GetFigure(DATA.tag.xcorr);
    xc = cmb.CalcXcorr(DATA,DATA.currentexpt(1),probes(1),probes(2));
end

if DATA.plot.voltxy == 4
    cmb.PlotTrodeXcorr(DATA,0);
end

DATA.playingspk = 0;
set(DATA.toplevel,'UserData',DATA);
cmb.SetGui(DATA);


