function PlotISI(X, varargin)
% PlotISI(X, varargin)
%Plot ISI histograms for Cluster Structures

callback = {@HitIsi};
type = [];
j = 1;
while j <= length(varargin)
    if sum(strncmpi(varargin{j},{'scatter' 'bars'},4))
        type = varargin{j};
    end
    j = j+1;
end

c = iscluster(X);
if c > 1
    Array = GetArrayConfig(X);
    [nr,nc] = Nsubplots(Array);
    for j = 1:length(X)
        mysubplot(nr,nc,j);
        PlotISI(X{j}, varargin{:});
        axis('tight');
        set(gca,'xtick',[],'ytick',[]);
    end
    return;
end

AddISIMenu(gcf);
P = getappdata(gcf,'ISIproperties');
if isempty(type)
    if isfield(P,'plottype')
        type = P.plottype;
    else
        type = 'bars'; %default
    end
end
P.data = X;

if isfield(X,'times')
    t = X.times;
end

isi = diff(t);
if strcmp(type,'bars')
     [a,b] = hist(isi,[0:0.0005:0.1]);
     hold off;
     bar(b(1:end-1),a(1:end-1));
     id = find(isi < 0.002);
     if isfield(P,'xmax')
         set(gca,'xlim', [0 P.xmax]);
     end
elseif strcmp(type,'scatter')
    hold off;
    plot(t(1:end-1),isi,'o','buttondownfcn',callback);
%    for j = 1:length(isi)
%        plot(t(j+1),isi(j),'o','buttondownfcn',callback);
%        hold on;
%    end
    set(gca,'ylim',[0 0.005]);
end
setappdata(gcf,'ISIproperties',P);

function HitIsi(a,b)
p = get(gca,'currentpoint');
t = get(a,'Xdata');
isi = get(a,'Ydata');
xr = get(gca,'xlim');
yr = get(gca,'ylim');
d = abs((p(1,1)-t)./diff(xr) + i * (p(1,2)-isi)./diff(yr));
[a,b] = min(d);

fprintf('ISI %.1fms at %.3f\n',isi(b) .* 1000, t(b));

function AddISIMenu(F)

it = findobj(F,'tag','ISIMenu');
if ~isempty(it)
    return;
end

hm = uimenu(F,'label','ISI','Tag','ISIMenu');
sm = uimenu(hm,'label','Max Interval');
uimenu(sm,'label', '100ms', 'callback',{@SetISI, 'xmax', 0.1});
uimenu(sm,'label', '50ms', 'callback',{@SetISI, 'xmax', 0.05});
uimenu(sm,'label', '20ms', 'callback',{@SetISI, 'xmax', 0.02});
sm = uimenu(hm,'label','Type');
uimenu(sm,'label', 'Bars', 'callback',{@SetISI, 'type', 'bars'});
uimenu(sm,'label', 'Sequence', 'callback',{@SetISI, 'type', 'scatter'});

function SetISI(a, b, fcn, val)

F = GetFigure(a);
X = getappdata(F,'ISIproperties');

if strcmp(fcn,'type')
    X.plottype = val;
    PlotISI(X.data, X.plottype);
elseif strcmp(fcn,'xmax')
    X.xmax = val;
    set(gca,'xlim', [0 X.xmax]);
end
setappdata(F,'ISIproperties',X);



