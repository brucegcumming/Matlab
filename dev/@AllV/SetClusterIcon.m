function h = SetClusterIcon(DATA)h = NaN;        if DATA.elmousept.shape == 0        str = 'O';    elseif DATA.elmousept.shape < 0        str = '';    else        str = '/';    end    if ishandle(DATA.clustericon)        delete(DATA.clustericon);    end    if ishandle(DATA.maintitle)        mysubplot(2,4,1);    x = get(DATA.maintitle,'extent');    h = text(0,0,str,'units','Normalized','VerticalAlignment','Top','fontweight','bold','fontsize',DATA.gui.fontsize(1));    set(h,'color',DATA.colors{DATA.currentcluster+1});    else        mysubplot(2,4,1);        h = text(0,0,str,'units','Normalized','VerticalAlignment','Top','fontweight','bold','fontsize',DATA.gui.fontsize(1));        set(h,'color',DATA.colors{DATA.currentcluster+1});    end    DATA.clustericon = h;if DATA.plottype == 3 && DATA.currentcluster ~= DATA.templatecluster        mysubplot(2,4,4);    h = text(1,0,['Template From Cl ' num2str(DATA.templatecluster)],'units','Normalized',...        'VerticalAlignment','Top',...        'HorizontalAlignment','Right',...        'color',DATA.colors{1+DATA.templatecluster},...        'fontweight','bold','fontsize',DATA.gui.fontsize(1));    end