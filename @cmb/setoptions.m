function cntrl_box = setoptions(DATA, tag)

wsc = DATA.wsc;
SPACE = 3 * wsc(1);
VSPACE = 5 * wsc(2);
ch = DATA.plot.ch;
cw = DATA.plot.cw;
h = (ch+VSPACE)*18;
w = cw * 40;
scrsz = get(0,'Screensize');

cntrl_box = findobj('Tag',tag,'Name','Options');
if ~isempty(cntrl_box)
    figure(cntrl_box);
    return;
end

Figpos = getappdata(DATA.toplevel,'Figpos');
bp = get(DATA.toplevel,'Position');
wpos =  [bp(1)+bp(3) bp(2) w*wsc(1) h*wsc(2)];
if ~isempty(Figpos)
    nf = strmatch(tag,fields(Figpos));
    if isempty(nf)
        nf = length(DATA.figpos)+1;
        Figpos.(tag) =  [bp(1)+bp(3) bp(2) w*wsc(1) h*wsc(2)];
        setappdata(DATA.toplevel,'Figpos',Figpos);
    else
        wpos = Figpos.(tag);
    end
end
dat.parentfigtag = DATA.tag.top;
cntrl_box = figure('Position', wpos , 'Menubar', 'none',...
    'NumberTitle', 'off', 'Tag',tag,'Name','Options','UserData',dat);
set(cntrl_box,'DefaultUIControlFontSize',DATA.gui.FontSize);
setappdata(cntrl_box,'ParentFigure',DATA.toplevel);
hm=uimenu(cntrl_box,'label','Spikes');
uimenu(hm,'label','From NEV','tag', 'nev', 'callback',@cmb.SpikeUIMenu);
uimenu(hm,'label','From NS5','tag', 'ns5', 'callback',@cmb.SpikeUIMenu);
uimenu(hm,'label','From Spks.mat','tag', 'spk', 'callback',@cmb.SpikeUIMenu);
hm=uimenu(cntrl_box,'label','Trigger');
uimenu(hm,'label','Both','tag', 'TrigBoth', 'callback',@cmb.SpikeUIMenu);
uimenu(hm,'label','Negative','tag', 'Trig-', 'callback',@cmb.SpikeUIMenu);
uimenu(hm,'label','Positive','tag', 'Trig+', 'callback',@cmb.SpikeUIMenu);

top = num2str(double(DATA.toplevel));
HSPACE = 0.02;
VSPACE = 0.02;
bp(1) = HSPACE;
cols(1) = HSPACE;
cols(2) = 0.26;
cols(3) = 0.5;
cols(4) = 0.74;
if isfield(DATA,'AllSpikes') || length(DATA.probelist) > 2
    rows = VSPACE + [0:15]./17;
else
    rows = VSPACE + [0:14]./16;
end
%start from bottom
bp(3) = 0.24;
colw = bp(3);
bp(4) = mean(diff(rows));
bp(1) = cols(1,1);
nr = 1;
nc = 1;
bp(1) = cols(1);
bp(2) = rows(nr);
%row(1) is the bottom
a = uicontrol(gcf,'Style', 'CheckBox','String','AutoCut Empty','units','norm','Position', bp,...
    'Tag','AutoCutEmpt','Callback',@cmb.Update,'value',DATA.state.autofit);
bp(1) = cols(2);
a = uicontrol(gcf,'Style', 'CheckBox','String','Load Last Layout','units','norm','Position', bp,...
    'Tag','UseLastLayout','Callback',@cmb.UpdateItem,'value',DATA.options.uselastlayout);
bp(1) = cols(3);
a = uicontrol(gcf,'Style', 'CheckBox','String','Load Last Config','units','norm','Position', bp,...
    'Tag','UseLastConfig','Callback',@cmb.UpdateItem,'value',DATA.options.uselastconfig);

nr = nr+1;
bp(2) = rows(nr);
a = uicontrol(gcf,'Style', 'CheckBox','String','Fit Expts','units','norm','Position', bp,...
    'Tag','AutoFit','Callback',@cmb.Update,'value',DATA.state.autofit);
bp(1) = cols(2);
uicontrol(gcf,'Style', 'CheckBox','String','Spool if No Cluster','units','norm','Position', bp,...
    'Tag','AutoSpool','Callback',@cmb.Update,'value',DATA.state.autospool);
bp(1) = cols(3);
uicontrol(gcf,'Style', 'CheckBox','String','V range','units','norm','Position', bp,...
    'Tag','AutoVrange','Callback',@cmb.Update,'value',DATA.plot.autoVrange);
if length(DATA.probelist) > 1
    bp(1) = cols(4);
    uicontrol(gcf,'Style', 'CheckBox','String','Update CellList plot','units','norm','Position', bp,...
        'Tag','AutoPlotCells','Callback',@cmb.Update,'value',DATA.state.autoplotcells);
end

nr = nr+1;

bp(2) = rows(nr); bp(1) = cols(1);
uicontrol(gcf,'Style', 'CheckBox','String','Plot Graph if change','units','norm','Position', bp,...
    'Tag','AutoReplotGraph','Callback',@cmb.Update,'value',DATA.state.autoreplotgraph);
bp(1) = cols(2);
uicontrol(gcf,'Style', 'CheckBox','String','Set List latest','units','norm','Position', bp,...
    'Tag','AutoSetList','Callback',@cmb.Update,'value',DATA.state.autosetlist);

bp(1) = cols(3);
uicontrol(gcf,'Style', 'CheckBox','String','Use Last Cluster','units','norm','Position', bp,...
    'Tag','ApplyLastCluster','Callback',@cmb.Update,'value',DATA.state.applylastcluster);
bp(1) = cols(4);
if length(DATA.probelist) > 1
    uicontrol(gcf,'Style', 'CheckBox','String','Plot All For New Probe','units','norm','Position', bp,...
        'Tag','AutoPlotNewProbe','Callback',@cmb.Update,'value',DATA.state.autoplotnewprobe);
end





bp(1) = cols(1);
nr = nr+1;
bp(2) = rows(nr);
bp(3) = 0.5;
uicontrol(gcf,'Style','Text','String','Automatic actions','units','norm','Position',bp);


coff = 1;
nr = 5;
bp(3) = colw;
%First Column of Checkbodex
bp(2) = rows(nr);
bp(1) = cols(1);
uicontrol(gcf,'Style', 'CheckBox','String','Limit Range','units','norm','Position', bp,...
    'Tag','FixRange','Callback',@cmb.Update,'value',DATA.state.fixrange);

nr=nr+1;
bp(2) = rows(nr);
uicontrol(gcf,'Style', 'CheckBox','String','ShowEM','units','norm','Position', bp,...
    'Tag','ShowEM','Callback',@cmb.Update,'value',DATA.plot.showem);


nr = nr+1;
bp(2) = rows(nr);
uicontrol(gcf,'Style', 'CheckBox','String','Flip','units','norm','Position', bp,...
    'Tag','Flip','Callback',@cmb.Update,'value',DATA.plot.flip);


nr = nr+1;
bp(2) = rows(nr);
%   uicontrol(gcf,'Style', 'CheckBox','String','Collapse 1','units','norm','Position', bp,...
%      'Tag','Collapse1','Callback',@cmb.Update,'value',DATA.plot.collapse);
uicontrol(gcf,'Style', 'CheckBox','String','Acov','units','norm','Position', bp,...
    'Tag','Acov','Callback',@cmb.Update,'value',DATA.plot.acov);

nr = nr+1;
bp(2) = rows(nr);
uicontrol(gcf,'Style', 'CheckBox','String','Spool Spikes','units','norm','Position', bp,...
    'Tag','ShowSpikes','Callback',@cmb.Update,'value',DATA.state.showspikes);


%second column of checks


nr = 5;
bp(2) = rows(nr);
bp(1) = cols(2);
uicontrol(gcf,'Style', 'CheckBox','String','Force rebuild','units','norm','Position', bp,...
    'Tag','ForceBuild','Callback',@cmb.Update,'value',DATA.state.forcebuild);


bp(2) = rows(5+coff);
uicontrol(gcf,'Style', 'CheckBox','String','dV vs V','units','norm','Position', bp,...
    'Tag','PhasePlot','Callback',@cmb.Update,'value',(DATA.plot.dvdt == 2));
bp(2) = rows(6+coff);
uicontrol(gcf,'Style', 'CheckBox','String','Remove DC','units','norm','Position', bp,...
    'Tag','RemoveDC','Callback',@cmb.Update,'value',DATA.plot.nodc);

bp(2) = rows(7+coff);
uicontrol(gcf,'Style', 'CheckBox','String','Auto relist','units','norm','Position', bp,...
    'Tag','AutoList','Callback',@cmb.Update,'value',DATA.state.autolist);

bp(2) = rows(8+coff);
uicontrol(gcf,'Style', 'CheckBox','String','Condense RC','units','norm','Position', bp,...
    'Tag','Condense','Callback',@cmb.Update,'value',DATA.plot.condenseRC);


% third column of checkboxes
bp(2) = rows(4+coff);
bp(1) = cols(3);
uicontrol(gcf,'Style', 'CheckBox','String','.pN in name','units','norm','Position', bp,...
    'Tag','NameProbe','Callback',@cmb.Update,'value',DATA.state.includeprobename);

bp(2) = rows(5+coff);
uicontrol(gcf,'Style', 'CheckBox','String','Plot F1','units','norm','Position', bp,...
    'Tag','PlotMod','Callback',@cmb.Update,'value',DATA.plot.plotmod);

bp(2) = rows(6+coff);
uicontrol(gcf,'Style', 'CheckBox','String','Auto Advance','units','norm','Position', bp,...
    'Tag','AutoNext','Callback',@cmb.Update,'value',DATA.state.autonext);

bp(2) = rows(7+coff);
uicontrol(gcf,'Style', 'CheckBox','String','Hide Spikes','units','norm','Position', bp,...
    'Tag','NoSpikes','Callback',@cmb.Update,'value',DATA.state.nospikes>0);

bp(2) = rows(8+coff);
uicontrol(gcf,'Style', 'CheckBox','String','ISIH','units','norm','Position', bp,...
    'Tag','ISIH','Callback',@cmb.Update,'value',DATA.plot.showISI);

% fourth column of checkboxes
bp(1) = cols(4);
bp(2) = rows(4+coff);
if ~isfield(DATA.state, 'optimizeclusters')
    DATA.state.optimizeclusters = 0;
end
%
%  uicontrol(gcf,'Style', 'CheckBox','String','Optimize','units','norm','Position', bp,...
%     'Tag','OptimizeClusters','Callback',@cmb.Update,'value',DATA.state.optimizeclusters);
uicontrol(gcf,'Style', 'CheckBox','String','Show Artifacts','units','norm','Position', bp,...
    'Tag','ShowArtifacts','Callback',@cmb.Update,'value',DATA.plot.showartifacts);
bp(2) = rows(5+coff);
uicontrol(gcf,'Style', 'CheckBox','String','xCorr','units','norm','Position', bp,...
    'Tag','ShowxCorr','Callback',@cmb.Update,'value',DATA.plot.xcorr);
bp(2) = rows(6+coff);
uicontrol(gcf,'Style', 'CheckBox','String','+hash','units','norm','Position', bp,...
    'Tag','AddHash','Callback',@cmb.Update,'value',DATA.plot.addhash);
bp(2) = rows(7+coff);
uicontrol(gcf,'Style', 'CheckBox','String','Wave T','units','norm','Position', bp,...
    'Tag','ShowWave','Callback',@cmb.Update,'value',DATA.plot.showwave);
%back to left side
bp(2) = rows(8+coff);
uicontrol(gcf,'Style', 'CheckBox','String','RF relative','units','norm','Position', bp,...
    'Tag','CenterRFMeasures','Callback',@cmb.Update,'value',DATA.plot.centerRFmeasures);


bp(2) = rows(9+coff);
bp(1) = cols(1);
uicontrol(gcf,'Style', 'pop','String',DATA.spkvarnames,'units','norm','Position', bp,...
    'Tag','ClusterX','Callback',@cmb.Update,'value',DATA.plot.clusterX);

bp(1) = cols(2);
uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.clusterXrange(1)),'units','norm','Position', bp,...
    'Tag','ClusterXmin','Callback',@cmb.Update);

bp(1) = cols(3);
uicontrol(gcf,'Style', 'pushbutton','String','Plot ISIH','units','norm','Position', bp,...
    'Callback','cmb.combine(''PlotISI'')');


bp(1) = cols(1);
bp(2) = rows(10+coff);
uicontrol(gcf,'Style', 'pop','String',DATA.spkvarnames,'units','norm','Position', bp,...
    'Tag','ClusterY','Callback',@cmb.Update,'value',DATA.plot.clusterY);

bp(1) = cols(2);
uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.clusterYrange(1)),'units','norm','Position', bp,...
    'Tag','ClusterYmin','Callback',@cmb.Update);

bp(1) = cols(3); bp(3) = 0.1;
uicontrol(gcf,'Style', 'text','string','AutoScale','units','norm','Position',bp);

bp(1) = cols(3)+0.1; bp(3) = colw-0.1;
uicontrol(gcf,'Style', 'pop','String','Matlab|99.9|%99%|95%|Tight','units','norm','Position', bp,...
    'Tag','AutoScaleMode','Callback',@cmb.Update,'value',DATA.plot.autoscalemode);

if 0
bp(1) = cols(4); bp(3) = 0.1;
uicontrol(gcf,'Style', 'text','string','AutoCut','units','norm','Position',bp);

bp(1) = cols(4)+0.1; bp(3) = colw-0.1;
uicontrol(gcf,'Style', 'pop','String','Centiles|By Density|Test|None','units','norm','Position', bp,...
    'Tag','AutoCutMode','Callback',@cmb.Update,'value',DATA.plot.autoclustermode+1);
end
bp(1) = cols(1); bp(3) = 0.1;
bp(2) = rows(11+coff);
uicontrol(gcf,'Style', 'text','string','Pt Range A','units','norm','Position',bp);
bp(1) = cols(1)+0.1; bp(3) = colw-0.1;
uicontrol(gcf,'Style', 'edit','string',cmb.vec2str(DATA.clusterArange),'units','norm','Position', bp,...
    'Tag','ClusterArange','Callback',@cmb.Update);

bp(1) = cols(2); bp(3) = 0.1;
uicontrol(gcf,'Style', 'text','string','B','units','norm','Position',bp);
bp(1) = cols(2)+0.1; bp(3) = colw-0.1;
uicontrol(gcf,'Style', 'edit','string',cmb.vec2str(DATA.clusterBrange),'units','norm','Position', bp,...
    'Tag','ClusterBrange','Callback',@cmb.Update);
bp(1) = cols(3); bp(3) = 0.1;
uicontrol(gcf,'Style', 'text','string','E','units','norm','Position',bp);
bp(1) = cols(3)+0.1; bp(3) = colw-0.1;
uicontrol(gcf,'Style', 'edit','string',cmb.vec2str(DATA.clusterErange),'units','norm','Position', bp,...
    'Tag','ClusterErange','Callback',@cmb.Update);

bp(1) = cols(1); bp(3) = 0.1;
bp(2) = rows(12+coff);
uicontrol(gcf,'Style', 'text','string','N Min','units','norm','Position',bp);
bp(1) = cols(1)+0.1; bp(3) = colw-0.1;
uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.nmin),'units','norm','Position', bp,...
    'Tag','Nmin','Callback',@cmb.Update);

bp(1) = cols(2); bp(3) = colw;
uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.nminrc),'units','norm','Position', bp,...
    'Tag','RCNmin','Callback',@cmb.Update);

if 0  %all done from main menu now
    bp(1) = cols(3); bp(3) = 0.1;
    uicontrol(gcf,'Style', 'text','string','CP','units','norm','Position',bp);
    bp(1) = cols(3)+0.1; bp(3) = colw-0.1;
    uicontrol(gcf,'Style', 'pop','String','None|CP-time|CP trials|CP-hist|CP-EM|Psych Only|Psych Smooth','units','norm','Position', bp,...
        'Tag','ShowCP','Callback',@cmb.Update,'value',DATA.plot.showcp+1);
end
bp(1) = cols(4); bp(3) = 0.1;
uicontrol(gcf,'Style', 'text','string','PtSz','units','norm','Position',bp);
bp(1) = cols(4)+0.1; bp(3) = colw-0.1;
uicontrol(gcf,'Style', 'pop','String','Auto|1|2|3|4|5|6|7|8', 'units','norm','Position',bp,...
    'Tag','SetPtSize','Callback',@cmb.Update,'value',DATA.plot.setptsize+1);


bp(2) = rows(13+coff);
bp(1) = cols(1); bp(3) = 0.1;
uicontrol(gcf,'Style', 'text','string','sdfw','units','norm','Position',bp);
bp(1) = cols(1)+0.1; bp(3) = colw-0.1;
uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.sdfw),'units','norm','Position', bp,...
    'Tag','Sdfw','Callback',@cmb.Update);

bp(1) = cols(2); bp(3) = 0.1;
uicontrol(gcf,'Style', 'text','string','MLFP','units','norm','Position',bp);
if ~isfield(DATA.plot,'lfpplot');
    DATA.plot.lfpplot = 0;
end
bp(1) = cols(2)+0.1; bp(3) = colw-0.1;
uicontrol(gcf,'Style', 'pop','String','None|Default|Stack|Image|Blank|RC|Movie|OneStim|monocs|Eig|Var|BlankVar|Frameresp|Trial|CSD|TrialStart|FTpwr','units','norm', 'Position',bp,...
    'Tag','LFPPlot','Callback',@cmb.Update,'value',DATA.plot.lfpplot+1);


bp(1) = cols(1);
bp(2) = rows(14+coff); bp(3) = colw;
uicontrol(gcf,'Style', 'text','string','Spike Display MaxV','units','norm','Position',bp);
bp(1) = cols(2);
if ~isfield(DATA.plot,'SpikeMaxV')
    DATA.plot.SpikeMaxV = 5;
end
if ~isfield(DATA.plot,'SpikeMinV')
    DATA.plot.SpikeMinV = -5;
end
if ~isfield(DATA.plot,'SpikeVsep')
    DATA.plot.SpikeVsep = 3;
end
uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.SpikeMinV),'units','norm','Position', bp,...
    'Tag','SpikeMaxV','Callback',@cmb.Update);
bp(1) = cols(3);
uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.SpikeMaxV),'units','norm','Position', bp,...
    'Tag','SpikeMaxV','Callback',@cmb.Update);
bp(1) = cols(4);
uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.SpikeVsep),'units','norm','Position', bp,...
    'Tag','SpikeVsep','Callback',@cmb.Update);

bp(1) = cols(1);
bp(3) = 0.1;
%       bp(1) = SPACE;

if isfield(DATA,'AllSpikes')
    bp(2) = rows(15+coff);
    str = 'Sync';
uicontrol(gcf,'Style', 'text','string',str,'units','norm','Position',bp);
bp(3) = 0.2;
bp(1) = cols(1)+0.1;
uicontrol(gcf,'Style', 'pop','String','Neg|Either|Pos|None|PlotbyNeg|PlotByPos|PlotByAll', 'units','norm','Position',bp,...
    'Tag','SyncSign','Callback',@cmb.Update,'value',DATA.syncsign+2);
    bp(1) = cols(2);
    uicontrol(gcf,'Style', 'CheckBox','String','Overlay','units','norm','Position', bp,...
        'Tag','SyncOverlay','Callback',@cmb.Update,'value',DATA.plot.syncoverlay)
    
    bp(1) = cols(3);
    uicontrol(gcf,'Style', 'CheckBox','String','SpikeBySpike','units','norm','Position', bp,...
        'Tag','SpikeBySpike','Callback',@cmb.Update,'value',DATA.plot.timebyspikeprobe)
    bp(1) = cols(4);
    uicontrol(gcf,'Style', 'pop','String','None|Min|Max|Energy|PCA1-1|PCA2-2|PCA1-2|Xcorr-PCA|CX-CX|CY-CY|Sum|test', 'units','norm','Position',bp,...
        'Tag','SyncCluster','Callback',@cmb.Update,'value',DATA.plot.synccluster+1);
elseif length(DATA.probelist) > 2
    bp(2) = rows(15+coff);
    bp(1) = cols(1);
    bp(3) = colw;
    uicontrol(gcf,'Style', 'text','string','CCF:Probe','units','norm','Position',bp);
    bp(1) = cols(2); bp(3) = colw;
    str = num2str(DATA.probelist');
    uicontrol(gcf,'Style', 'pop','String',str, 'units','norm','Position',bp,...
        'Tag','CalcCCF','Callback',@cmb.UpdateItem,'value',DATA.probe);
end


