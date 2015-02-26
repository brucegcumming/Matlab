function PlayNextSpike(a,b)
%DATA = cmb.combine('getstate');
DATA = GetDataFromFig(a);
if DATA.ISIpair
    cmb.PlotISIPair(DATA,DATA.ISIpair);
    DATA.ISIpair = DATA.ISIpair+1;
    set(DATA.toplevel,'UserData',DATA);
    return;
end
set(0,'CurrentFigure',DATA.svfig);
if isfield(DATA,'AllSpikes')
    cmb.PlotSpike(DATA,DATA.currentspike,[DATA.probe DATA.xprobes]);
else
    cmb.PlotSpike(DATA,DATA.currentspike);
    if isfigure(DATA.xyfig)
        set(0,'CurrentFigure',DATA.xyfig);
        cmb.PlotSpikeXY(DATA,DATA.currentspike,DATA.spkcolor{DATA.AllData.Spikes.codes(DATA.currentspike,2)+1});
    end
end
DATA.currentspike = DATA.currentspike+1;

if length(DATA.probes) > 1 && DATA.plot.timebyspikeprobe > 0
    p = DATA.probe;
    probes = find(DATA.plot.useprobe);
    set(0,'CurrentFigure',DATA.timefig);
    hold off;
    tw = 100;
    set(gca,'xlim',[0 tw*2]);
    nsmp = size(DATA.AllSpikes{p}.values,2);
    hscale = DATA.AllSpikes{p}.interval * 10000;
    ts = DATA.AllSpikes{p}.times(DATA.currentspike);
    for j = 1:length(probes)
        [sids{j}, alltimes{j} allcodes{j}] = FindSpikes(DATA, [ts-tw ts+tw], probes(j),[]);
        voff(j) = DATA.plot.SpikeVsep*(j-1);
        x = [];
        y = [];
        id = find(alltimes{j} > ts-tw & alltimes{j} < ts+tw);
        for t = 1:length(id)
            x = (alltimes{j}(id(t))-ts+tw) +[1:nsmp]*hscale;
            y = DATA.AllSpikes{probes(j)}.values(sids{j}(id(t)),:)+voff(j);
            plot(x,y,'color',DATA.spkcolor{allcodes{j}(id(t))+1});
            hold on;
        end
    end
    set(gca,'xlim',[0 tw*2],'ylim',[-DATA.plot.SpikeMaxV max(voff)+DATA.plot.SpikeMaxV]);
    drawnow;
end


set(DATA.toplevel,'UserData',DATA);

