function PlotAllExptProbeMean(DATA, type, varargin)    if strcmp(type,'AllMeanIm')        type = 'imageonly';    elseif strcmp(type,'AllExptMeanIm')        type = 'smallimage';    elseif strcmp(type,'AllMean')        type = 'lineonly';    elseif strcmp(type,'AllExptIm')        type = 'exptim';    end            Clusters = getappdata(DATA.toplevel,'Clusters');    eid = DATA.currentpoint(1);    PC.SetFigure(DATA,DATA.tag.allxy);    cmenu = PC.AddCellContextMenu(DATA, 'subplot');    axdata.toplevel = DATA.toplevel;        if length(DATA.selectexpts)        expts = DATA.selectexpts;    else        expts = 1:length(Clusters);    end    if length(DATA.proberange)        probelist = DATA.proberange;    else       probelist = 1:DATA.nprobes;    end    nex = length(expts);    for eid = 1:length(expts)    for j = 1:length(Clusters{eid})        mins(eid,j) = min(Clusters{expts(eid)}{j}.MeanSpike.ms(:));        maxs(eid,j) = max(Clusters{expts(eid)}{j}.MeanSpike.ms(:));    end    end            clim = [min(mins(:)) max(maxs(:))];    ClearPlot;    for e = 1:length(expts)    for j = 1:length(probelist)        eid = expts(e);        p = probelist(j);        mysubplot(nex,length(probelist),j + (e-1)*length(probelist));        hold off;         PC.PlotMeanSpike(Clusters{eid}{p},0,1,type);        set(gca,'Xtick',[],'Ytick',[],'clim',clim);        h = get(gca,'title');        xl = get(gca,'xlim');        yl = get(gca,'ylim');        a = get(h,'position');        a(2) = yl(2);        a(1) = mean(xl);        set(h,'position',a,'VerticalAlignment','top');        delete(h);        axdata.probe = j;        axdata.eid = eid;        set(gca,'UIContextMenu',cmenu, 'UserData', axdata);    end    drawnow;    end    