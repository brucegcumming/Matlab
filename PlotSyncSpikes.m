function PlotSyncSpikes(DATA, eid, probes, clnum, varargin)
j = 1;
if isfield(DATA.plot,'syncmax') && DATA.plot.synctmax > 0
    tlim = DATA.plot.synctmax .* 10;
else
    tlim = 20;
end
    
nstep = 10;
while j <= length(varargin)
    if strncmpi(varargin{j},'flip',4)
        probes = fliplr(probes);
    elseif strncmpi(varargin{j},'tmax',4)
        j = j+1;
        tlim = varargin{j}(1);
        if length(varargin{j}) > 1
            nstep = round(varargin{j}(2)/2);
        end
    end
    j = j+1;
end

Clusters = getappdata(DATA.toplevel,'Clusters');
if iscell(eid) %called from AllVPcs
    Spikes = eid;
elseif eid > 0 %called from PlotClusters
 AllSpikes = CheckAllSpikes(DATA, eid, probes);
 Spikes{1} = AllSpikes(eid,probes(1));
 Spikes{2} = AllSpikes(eid,probes(2));
 Clusters = getappdata(DATA.toplevel,'Clusters');
 Clusters = Clusters{eid};
 atid = find(Clusters{probes(1)}.clst == clnum(1)+1);
 btid = find(Clusters{probes(2)}.clst == clnum(2)+1);
 
 if isempty(AllSpikes)
     xct = -(tlim/10000):0.0002:(tlim/10000);
     [xc, b] = xcorrtimes(Clusters{probes(1)}.times(atid),Clusters{probes(2)}.times(btid),'times',xct);
     SetFigure(DATA,DATA.tag.spikes);
     plot(b.xpts, xc);
     return;
 end
 Spikes{1}.times = Spikes{1}.times(atid);
 Spikes{2}.times = Spikes{2}.times(btid);
 Spikes{1}.values = Spikes{1}.values(atid,:);
 Spikes{2}.values = Spikes{2}.values(btid,:);
end
 ta = Spikes{1}.times;
 tb = Spikes{2}.times;
 [a, aid, bid]= intersect(round(ta/10),round(tb/10));
% bid = find(ismember(round(ta/10),round(tb+10/10)));
 
 
aid = [];
bid = [];
dts = [];
for j = 1:length(ta)
    dt = ta(j)-tb;
%find all events witin tlime and show them    
    id = find(abs(dt) < tlim)';
    if length(id)
        aid = [aid j * ones(size(id))];
        bid = [bid id];
        dts = [dts dt(id)'];
    end
end
% id = union(aid,bid);
voff = Spikes{1}.maxv .* Spikes{1}.maxint/2;
x = [1:size(Spikes{1}.values,2)]./40;
hold off;
 for j = 1:length(aid)
     %+v toff = 1 after 2, so subtract from 2
     toff = Spikes{1}.times((aid(j)))-Spikes{2}.times((bid(j)));
     toff = toff ./10;
     plot(x,Spikes{1}.values((aid(j)),:),'k');
     hold on;
     plot(x-toff,double(Spikes{2}.values((bid(j)),:))+voff,'b');
 end
 text(x(end),voff,sprintf('%d/%d',probes(2),clnum(2)),'horizontalalignment','left','color','b','fontsize',20);
 text(x(end),0,sprintf('%d/%d',probes(1),clnum(1)),'horizontalalignment','left','color','k','fontsize',20);
 yl = get(gca,'ylim');
 xl = get(gca,'xlim');
 if tlim <= 50
     [y,x] = hist(dts,-tlim:2:tlim);
 else
     [y,x] = hist(dts,-tlim:5:tlim);
 end
 pmax = max(y)./length(Spikes{1}.times);
y = y .* yl(2)./max(y);
x  = -x./10;
title(sprintf('%d events (of %d) within +- %.2fms pmax %.3f',length(aid),length(Spikes{1}.times),tlim./10,pmax));
%x = (x-min(x)) .* xl(2)./max(x);
%x = x+xl(1);
plot(x,y,'ro-');
