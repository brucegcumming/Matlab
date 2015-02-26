function [figpos, F] = SetFigure(tag, DATA, varargin)

figpos = getappdata(DATA.toplevel,'Figpos');
j = 1;
args = {'keepmenu'};
while j < length(varargin)
    args = {args{:} varargin{j}};
    j =  j+1;
end
cw = DATA.plot.cw;
ch = DATA.plot.ch;
rh = ch+10;

[F,isnew] = GetFigure(tag,args{:},'parent',DATA.toplevel);
SetFigPos(DATA, tag);
if isnew && strcmp(tag,DATA.tag.allexpts)
    hm = uimenu(F,'label','&Plots');
    uimenu(hm,'label','&Next','Callback',{@cmb.AllCellPlots, 'next'});
    uimenu(hm,'label','&Prev','Callback',{@cmb.AllCellPlots, 'prev'});
    uimenu(hm,'label','Choose','Tag','ExptCellChoose');
    sm = uimenu(hm,'label','&Sort');
    uimenu(sm,'label','by &Var','Callback',{@cmb.AllCellPlots, 'sortbyvar'});
    uimenu(sm,'label','by &Probe','Callback',{@cmb.AllCellPlots, 'sortbyprobe'});
    uimenu(hm,'label','&Save AllExpts','Callback',{@cmb.AllCellPlots, 'save'});
    uimenu(hm,'label','&Xcorrs','Callback',{@cmb.AllCellPlots, 'xcorr'});
    set(F,'KeyPressFcn',@cmb.AllCellKeyPressed);
    set(F,'UserData',DATA.toplevel);
    cmb.SetCellChooser(DATA);
    DATA.fig.CombinerAllExpts = F;
    set(DATA.toplevel,'UserData',DATA);
elseif isnew && strcmp(tag,DATA.tag.clusterxy)
    xyfig = F;
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
if strcmp(tag,DATA.tag.clusterxy)
    xyfig = F;
    set(xyfig, 'KeyPressFcn',@cmb.KeyPressed);
    set(xyfig, 'KeyReleaseFcn',@cmb.KeyReleased);
    set(xyfig, 'WindowButtonDownFcn',@cmb.ButtonPressed);
    set(xyfig, 'WindowButtonMotionFcn',@cmb.ButtonDragged);
    set(xyfig, 'WindowButtonUpFcn',@cmb.ButtonReleased);
    set(xyfig, 'WindowScrollWheelFcn',@cmb.ScrollWheel);
end

end

