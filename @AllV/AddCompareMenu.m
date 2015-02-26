function AddCompareMenu(DATA, varargin)

C = AllV.mygetappdata(DATA,'Clusters');
hm = findobj(DATA.toplevel,'type','uimenu','tag','CompareMenu');
if isempty(hm)
    a = findobj(DATA.toplevel,'type','uimenu','tag','ClusterMenu');
    hm = uimenu(a,'Label','Compare With','tag','CompareMenu');
else
    delete(get(hm,'Children'));
end

if isfield(DATA,'ArrayConfig') && isfield(DATA.ArrayConfig,'X')
    p = AllV.ProbeNumber(DATA);
    d = DATA.ArrayConfig.X - DATA.ArrayConfig.X(p) + i *  (DATA.ArrayConfig.Y - DATA.ArrayConfig.Y(p));
    d = d .* DATA.ArrayConfig.spacing;
else
    d = zeros(size(C));
end
    
for j = 1:length(C)
    if abs(d(j)) <= 200 && isfield(C{j},'next')
    str = sprintf('%d/1',j);
    pval = j + 0.1;
    sm = uimenu(hm,'Label',str,'callback', {@CompareClusters, pval},'foregroundcolor',DATA.colors{2});
    for k = 1:length(C{j}.next)
        if ~isempty(C{j}.next{k})
            str = sprintf('%d/%d',j,k+1);
            pval = j + (k+1)/10;
            sm = uimenu(hm,'Label',str,'callback', {@CompareClusters, pval},'foregroundcolor',DATA.colors{k+2});
        end
    end
    end
end


function CompareClusters(a,b, p)


DATA = GetDataFromFig(a);
if DATA.currentcluster > 1
    ms = DATA.cluster.next{DATA.currentcluster-1}.MeanSpike.ms;
    Ca = DATA.cluster.next{DATA.currentcluster-1};
else
    ms = DATA.cluster.MeanSpike.ms;
    Ca = DATA.cluster;
end
compcl = round(mod(p,1).*10);
probe = p;
DATA = AllV.CheckTrigTimes(DATA);
p = floor(p);
C = AllV.mygetappdata(DATA,'Clusters');
Cb = AllV.GetClusterDef(C,probe);
DATA = AllV.LoadTrigTimes(DATA,[p]);
ClusterDetails = AllV.mygetappdata(DATA, 'ClusterDetails')
msb = Cb.MeanSpike.ms;
V = AllV.mygetappdata(DATA,'AllVoltages');
V = reshape(V,size(V,1)*size(V,2),size(V,3));
x = ms(:)' * V;
y = msb(:)' * V;
F = AllV.SetFigure(DATA.tag.comparexy,DATA);
hold off;
plot(x,y,'.','color',DATA.colors{1},'markersize',1);
idx = find(DATA.clst == DATA.currentcluster+1);
ns(2) = length(idx);
hold on;
Vb = AllV.BuildAllV(DATA,DATA.trigtimes{floor(p)},Cb.spts,'noset');
[a,b,c] = CalcEfficacy(Ca.times,Cb.times);
ns(1) = length(c.id{1});
hold on;
Vb = reshape(Vb,size(Vb,1)*size(Vb,2),size(Vb,3));
xb = ms(:)' * Vb;
yb = msb(:)' * Vb;
plot(xb,yb,'.','color',[0.3 0.3 0.3],'markersize',1);


h = plot(x(idx),y(idx),'.','color',DATA.colors{2},'markersize',DATA.ptsz(1));
set(h,'UserData',Ca.times,'buttondownfcn',@HitPoint);
if length(idx) < 500
    set(h,'markersize',6);
end
set(DATA.toplevel,'UserData',DATA);
h = plot(x(idx(c.id{1})),y(idx(c.id{1})),'mo','buttondownfcn',@HitPoint);
set(h,'UserData',Ca.times(c.id{1}));

idx = find(ClusterDetails{p}.clst == compcl+1);
ns(3) = length(idx);
h = plot(xb(idx),yb(idx),'.','color',DATA.colors{3},'markersize',DATA.ptsz(1));
set(h,'UserData',Cb.times,'buttondownfcn',@HitPoint);
if length(idx) < 500
    set(h,'markersize',6);
end
h = plot(xb(idx(c.id{2})),yb(idx(c.id{2})),'ms');
set(h,'UserData',Cb.times(c.id{2}),'buttondownfcn',@HitPoint);
clusterpair.C{1} = Ca;
clusterpair.C{2} = Cb;
setappdata(F,'ClusterPair',clusterpair);

title(sprintf('Comparing %d.%d with %.1f: %d/(%d,%d) syncs',DATA.probe(1),DATA.currentcluster,probe,ns));
xlabel(sprintf('*MeanSpk %d.%d',DATA.probe(1),DATA.currentcluster));
ylabel(sprintf('*MeanSpk %d.%d',p,compcl));
%refline(1);
[xc, details] = xcorrtimes(Ca.times,Cb.times);
AllV.SetFigure(DATA.tag.xcorr, DATA);
plot(details.xpts, xc);


function HitPoint(a,b)
DATA = GetDataFromFig(a);
p = get(gca,'currentpoint');
x = get(gco,'xdata');
y = get(gco,'ydata');
xy= p(1,1) + i * p(1,2);
d = abs(xy - (x + i.*y));
[a,b] = min(d);
t = get(gco,'UserData');
tid = FindTrial(DATA.Expt,t(b));
fprintf('Spikes at %.3f Trial %d(%d) (id%d) \n',t(b),DATA.Expt.Trials(tid).Trial,tid,DATA.Expt.Trials(tid).id);
if ~isempty(t)
    C = getappdata(gcf,'ClusterPair');
    V = AllV.mygetappdata(DATA,'Vall');
    Vb = AllV.mygetappdata(DATA,'AllVoltages');
    AllV.SetFigure(DATA.tag.fullv,DATA);

    hold off;
    chspk = union(DATA.chspk, C.C{2}.chspk);
    chspk = union(chspk, C.C{1}.chspk);
    DATA.voffset = CalcVoffset(Vb, chspk,DATA.gui.spikeVoverlap);
    voff = DATA.voffset- DATA.voffset(AllV.ProbeNumber(DATA));
    
    [V,id] = AllV.PlotFullV(DATA,[t(b)-0.001 t(b)+0.001],'chspk',chspk,'markprobes');
    hold on;
    ts = t(b) + DATA.cluster.spts./40000;
    for j = chspk
        plot(ts, C.C{1}.MeanSpike.ms(j,:)+voff(j),'r');
        plot(ts, C.C{2}.MeanSpike.ms(j,:)+voff(j),'g');
    end
    [a,id] = min(abs(t(b)-DATA.t));
    if a < 0.001
        GetFigure(DATA.toplevel);
        if ismember(DATA.plottype, [3 4 7])
            p = DATA.tmplots;
            for j = 1:8
                mysubplot(2,4,j);
                hold on;
                h = plot(DATA.TemplateScores(id,p(j,1)),DATA.TemplateScores(id,p(j,2)),'kx','markersize',10);
            end
        end
    else
        fprintf('No Spike at %.4f for current probe\n',t(b));
    end
end
