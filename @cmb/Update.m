function Update(a,b)

guichange = 0;
DATA = GetDataFromFig(a);
DATA = cmb.combine('getstate');
caller = get(a,'Tag');

if strcmp(caller,'AllDensity')
    cmb.PlotAllProbeXY(DATA);
    return;
end

if strcmp(caller,'StopSpool') && get(a,'value') == 0
    DATA.playingspk = 0;
end


DATA.plot.quickspks = cmb.GetCheck('QuickSpks',DATA.toplevel);

it = findobj('Tag',DATA.tag.options);
if ~isempty(it)
    DATA.plot.autoclustermode = get(findobj(it,'Tag','AutoCutMode'),'value')-1;
    DATA.state.fixrange = cmb.GetCheck('FixRange');
    DATA.plot.showem = cmb.GetCheck('ShowEM',it);
    %DATA.plot.showcp = get(findobj(it,'Tag','ShowCP'),'value')-1;
    DATA.plot.condenseRC = cmb.GetCheck('Condense',it);
    DATA.plot.clusterX = get(findobj(it,'Tag','ClusterX'),'value');
    DATA.plot.clusterY = get(findobj(it,'Tag','ClusterY'),'value');
    DATA.plot.clusterXrange(1) = cmb.GetField('ClusterXmin',it);
    DATA.plot.clusterYrange(1) = cmb.GetField('ClusterYmin',it);
    DATA.clusterArange = cmb.GetField('ClusterArange',it);
    DATA.clusterBrange = cmb.GetField('ClusterBrange',it);
    DATA.clusterErange = cmb.GetField('ClusterErange',it);
    DATA.plot.nmin = cmb.GetField('Nmin',it);
    DATA.plot.nminrc = cmb.GetField('RCNmin',it);
    DATA.plot.sdfw = cmb.GetField('Sdfw',it);
    DATA.plot.SpikeMaxV = cmb.GetField('SpikeMaxV',it);
    DATA.plot.SpikeVsep = cmb.GetField('SpikeVsep',it);
    DATA.plot.DensitySigma = [3 3];
    DATA.plot.addhash = cmb.GetCheck('AddHash',it);
    DATA.plot.xcorr = cmb.GetCheck('ShowxCorr',it);
    DATA.plot.autoVrange = cmb.GetCheck('AutoVrange',it);
    DATA.plot.showartifacts = cmb.GetCheck('ShowArtifacts',it);
    if length(DATA.plot.clusterX) > 1 || length(DATA.plot.clusterY) > 1
        fprintf('ClusterX/Y is too big');
    end
    DATA.plot.acov = cmb.GetCheck('Acov',it);
    DATA.state.includeprobename = cmb.GetCheck('NameProbe',it);
    DATA.state.autofit = cmb.GetCheck('AutoFit',it);
    DATA.plot.flip = cmb.GetCheck('Flip',it);
    DATA.plot.collapse = cmb.GetCheck('Collapse1',it);
    DATA.plot.showISI = cmb.GetCheck('ISIH',it);
    [DATA.state.autolist, h] = cmb.GetCheck('AutoList');
    DATA.plot.setptsize = get(findobj(it,'Tag','SetPtSize'),'value')-1;
    DATA.plot.lfpplot = get(findobj(it,'Tag','LFPPlot'),'value')-1;
    DATA.plot.centerRFmeasures = cmb.GetCheck('CenterRFMeasures',it);
    DATA.state.autoplotnewprobe = cmb.GetCheck('AutoPlotNewProbe',it);
    DATA.state.autospool = cmb.GetCheck('AutoSpool',it);
    DATA.state.autoplotcells = cmb.GetCheck('AutoPlotCells',it);
    DATA.state.autoreplotgraph = cmb.GetCheck('AutoReplotGraph',it);
    DATA.state.autosetlit = cmb.GetCheck('AutoSetList',it);
    DATA.plot.autoscalemode = get(findobj(it,'Tag','AutoScaleMode'),'value');
    DATA.state.applylastcluster = cmb.GetCheck('ApplyLastCluster',it);
    
%    s = get(findobj(it,'Tag','SyncSign') ,'value')-2;
%    if isempty(DATA.xprobes) && ~isempty(s)
%        DATA.TriggerSign = s;
%    end
    if isfield(DATA,'AllSpikes')
        if ~isempty(s)
            DATA.syncsign = s;
            DATA.plot.synccluster = get(findobj(it,'Tag','SyncCluster') ,'value')-1;
            DATA.plot.timebyspikeprobe = cmb.GetCheck('SpikeBySpike',it);
        end
        DATA.plot.syncoverlay = cmb.GetCheck('SyncOverlay',it);
    end
end


if (double(DATA.xyfig) > 0  & get(a, 'Parent') == DATA.xyfig) || strcmp(caller,'AutoScaleMode')
    set(0,'currentfigure',DATA.xyfig);
    DATA.plot.autoscale = cmb.GetCheck('AutoScale',DATA.xyfig) * DATA.plot.autoscalemode;
    it = findobj(DATA.xyfig,'Tag','Clusterid');
    DATA.currentcluster = get(it,'value');
    ax = findobj(DATA.xyfig,'Type','axes');
    if DATA.plot.autoscale
        set(ax,'Ylimmode','auto','Xlimmode','auto');
        DATA.plot.clusterXrange  = get(ax,'Xlim');
        DATA.plot.clusterYrange  = get(ax,'Ylim');
        DATA = cmb.SetXYRanges(DATA);
        cmb.SetField(DATA.xyfig,'ClusterXmax',DATA.plot.clusterXrange(2));
        cmb.SetField(DATA.xyfig,'ClusterYmax',DATA.plot.clusterYrange(2));
    end
    set(ax,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
    DATA.state.currentclusterquality = get(findobj(DATA.xyfig,'Tag','ClusterQuality'),'value')-1;
    if strcmp(caller,'ClusterQuality')
        DATA.Expts{DATA.currentexpt(1)}.Cluster{DATA.currentcluster,DATA.probe}.quality = ...
            DATA.state.currentclusterquality;
    end
    it = findobj(DATA.xyfig,'Tag','ForceClusterid');
    if it
        if strcmp(caller,'Clusterid')
            set(it,'value',DATA.currentcluster);
        end
        DATA.forceclusterid = get(it,'value');
    end
elseif myhandle(DATA.xyfig) & strcmp(caller,'SetPtSize');
    figure(DATA.xyfig); hold off;
    cmb.DrawXYPlot(DATA,DATA.Expts{DATA.currentexpt(1)}.gui.spks);
end
DATA.plot.dvdt = cmb.GetCheck('dVdt');
DATA.plot.nodc = cmb.GetCheck('RemoveDC');
dvv = cmb.GetCheck('PhasePlot');
if DATA.plot.dvdt && dvv
    DATA.plot.dvdt = 2;
end

DATA.state.resetclusters = cmb.GetCheck('ResetClusters',DATA.toplevel);
DATA.state.showspkxy = cmb.GetCheck('SpkXY');
DATA.state.recount = cmb.GetCheck('Recount');
DATA.state.plotpsych = cmb.GetCheck('PlotPsych');
DATA.state.plotcombined = cmb.GetCheck('PlotCombined');
DATA.plot.plotmod = cmb.GetCheck('PlotMod');
DATA.plot.showsync = cmb.GetCheck('ShowSync');
DATA.plot.showwave = cmb.GetCheck('ShowWave');
DATA.plot.showN = cmb.GetCheck('ShowN');
DATA.state.showspikes = cmb.GetCheck('ShowSpikes');


[DATA.state.autoplot, h] = cmb.GetCheck('AutoPlot');
DATA.state.autonext = cmb.GetCheck('AutoNext');
if DATA.state.nospikes < 2
    DATA.state.nospikes = cmb.GetCheck('NoSpikes');
end
DATA.state.optimizeclusters = cmb.GetCheck('OptimizeClusters');
DATA.spikelist = cmb.WhichClusters(DATA.toplevel);
DATA.state.uselfp = get(findobj(DATA.toplevel,'Tag','UseLFP'),'value');
DATA.state.forcebuild = cmb.GetCheck('ForceBuild');
if isfigure(DATA.xyfig)
    DATA.state.currentclusterquality = get(findobj(DATA.xyfig,'Tag','ClusterQuality'),'value')-1;
end

DATA.state.plotseq = get(findobj(DATA.toplevel,'Tag','PlotSeq'),'value')-1;
strs = get(findobj(DATA.toplevel,'Tag','PlotSeq'),'string');
DATA.plot.type = deblank(strs(DATA.state.plotseq+1,:));
if isempty(DATA.plot.autoclustermode)
    DATA.plot.autoclustermode = 1;
end
for j = DATA.probelist
    DATA.plot.useprobe(j) = cmb.GetCheck(['UseProbe' num2str(j)]);
end

if DATA.state.online %%no LFP available
    DATA.state.uselfp = 0;
end
id = regexp(DATA.outname,'\.c[0-9]\.');
if id & DATA.spikelist >= 0
    DATA.outname(id+2) = num2str(DATA.spikelist(1));
    set(DATA.saveitem,'string',DATA.outname);
end

if DATA.state.plotseq == 5
    DATA.plot.showcp = 5;
end


needplot = 0;
if strmatch(caller,{'PlotSeq'})
    guichange = 1;
    s = get(a,'value');
    if ismember(s,[1 2 3])% no cp or psych
        DATA.state.plotpsych = 0;
        DATA.plot.showcp = 0;
    elseif s == 5 %CP - trials
        GetFigure(DATA.tag.dataplot);
        hold off;
        PlotRateSequence(DATA.Expts(DATA.expid));
        return;
    elseif s == 8 %CP - trials
        DATA.state.plotpsych = 1;
        DATA.plot.showcp = 2;
    elseif s == 9 %CP - Hist
        DATA.state.plotpsych = 1;
        DATA.plot.showcp = 3;
    elseif s == 6 %Psych Only
        DATA.state.plotpsych = 1;
        DATA.plot.showcp = 5;
    else
        DATA.plot.showcp = 0;
    end
    if isfield(DATA,'Expt') %have a combined expt
        cmb.PlotCombined(DATA,DATA.Expt);
    end
    needplot = 0;
elseif strmatch(caller,{'Psych'})
    cmb.PlotCombined(DATA,DATA.Expt);
    needplot = 0;
elseif strmatch(caller,{'SpkXY' 'ShowSpikes'})
    if (DATA.state.showspkxy || DATA.state.showspikes) && DATA.state.nospikes ~= 2
        DATA.state.nospikes = 0;
    end
elseif strmatch(caller,{'LFPPlot'}) & isappdata(DATA.toplevel,'LFPExpt')
    cmb.ShowLFPPlot(DATA);
    needplot = 1;
end
set(DATA.toplevel,'UserData',DATA);
if strmatch(caller,{'ClusterX' 'ClusterY'})
    if isfield(DATA,'spklist') && ~isempty(DATA.spklist)
        expspks = DATA.spklist;
    else
        expspks = DATA.spkrange(1):DATA.spkrange(2);
    end
    GetFigure(DATA.xyfig);
    hold off;
    DATA = CalcClusterVars(DATA, expspks);
    if DATA.plot.setptsize
        DATA.ptsize = DATA.plot.setptsize;
    end
    cmb.DrawXYPlot(DATA,expspks);
end

% used to read
%if DATA.state.autoplot & a ~= h & ~ DATA.state.showspikes
% but h is not set. What was this?

id = get(DATA.elst,'value');

if DATA.state.autoreplotgraph & needplot & (~DATA.state.showspikes | length(id) > 1)
    cmb.combine('setexp','flagchange');
end

if guichange
    cmb.SetGui(DATA);
end

