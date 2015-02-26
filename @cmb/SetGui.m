function SetGui(DATA);
%cmb.SetCheck('Recut',DATA.state.recut > 0); %% in case it is 2
cmb.SetCheck('Recount',DATA.state.recount,DATA.toplevel);
cmb.SetCheck('AutoPlot',DATA.state.autoplot,DATA.toplevel);
cmb.SetCheck('ShowSpikes',DATA.state.showspikes,DATA.toplevel);
cmb.SetCheck('SpkXY',DATA.state.showspkxy,DATA.toplevel);
cmb.SetCheck('ShowN',DATA.plot.showN,DATA.toplevel);
cmb.SetCheck('PlotPsych',DATA.state.plotpsych,DATA.toplevel);
cmb.SetCheck('ResetClusters',DATA.state.resetclusters,DATA.toplevel);
cmb.SetCheck('QuickSpks',DATA.plot.quickspks,DATA.toplevel);

cmb.SetClusterCheck(DATA);
if isfigure(DATA.xyfig)
    it = findobj(DATA.xyfig, 'Tag','Density');
    if double(it)
        if DATA.densityplot
            set(it,'value',1);
        else
            set(it,'value',0);
        end
        
        ax = findobj(DATA.xyfig,'Type','axes');
        if DATA.plot.autoscale == 1
            set(ax,'Ylimmode','auto','Xlimmode','auto');
            DATA.plot.clusterXrange  = get(ax(1),'Xlim');
            DATA.plot.clusterYrange  = get(ax(1),'Ylim');
            set(ax,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
            cmb.SetField(DATA.xyfig,'ClusterXmax',DATA.plot.clusterXrange(2));
            cmb.SetField(DATA.xyfig,'ClusterYmax',DATA.plot.clusterYrange(2));
            set(ax,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
            cmb.SetField(DATA.xyfig,'ClusterXmax',DATA.plot.clusterXrange(2));
            cmb.SetField(DATA.xyfig,'ClusterYmax',DATA.plot.clusterYrange(2));
        elseif ismember(DATA.plot.autoscale,[2 3])
            x = get(ax(1),'Xlim');
            y = get(ax(1),'Ylim');
            set(ax,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
        else
            if diff(DATA.plot.clusterXrange) > 0
                set(ax,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
            end
        end
    end
end
if isfigure(DATA.optionfig)
    cmb.SetField(DATA.optionfig,'ClusterXmin',DATA.plot.clusterXrange(1));
    cmb.SetField(DATA.optionfig,'ClusterYmin',DATA.plot.clusterYrange(1));
    set(findobj(DATA.optionfig,'Tag','ClusterX'),'value',DATA.plot.clusterX);
    set(findobj(DATA.optionfig,'Tag','ClusterY'),'value',DATA.plot.clusterY);
end
if size(DATA.probenames,2) > 1
    DATA.probenames = DATA.probenames';
end
if DATA.listbycell == 0
    it = findobj(DATA.toplevel,'Tag','ProbeId');
    str = get(it,'string');
    if length(str) ~= length(DATA.probenames) || sum(strcmp(str,DATA.probenames)) ==0
        set(it,'string',DATA.probenames);
    end
end

