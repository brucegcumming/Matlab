function cells = SpikeHist(list, varargin)
%
% [w, dp] = SpikeHist(list, ...)
% calculates a histogram of spike widths, for expt files name in 
% a cell stucture array list.
% SpikeHist(list, 'dp', th)
%     makes a histogram only for cells with isolation index > th
% if list is a single character string, the spike waveform for that
% file is shown.

j = 1;
dpcrit = 0;
plottype =1;
findarea = 2;
cellids = [];

while j <= nargin-1
    if strncmpi(varargin{j},'dprime',2)
        j = j+1;
        dpcrit = varargin{j};
    elseif strncmpi(varargin{j},'ids',2)
        j = j+1;
        cellids = varargin{j};
    elseif strncmpi(varargin{j},'plottype',4)
        j = j+1;
        if isnumeric(varargin{j})
            plottype = varargin{j};
        elseif strncmpi(varargin{j},'width',3)
            plottype = 1;
        elseif strncmpi(varargin{j},'zerocrossing',3)
            plottype = 2;
        end
    end
    j = j+1;
end

if isstruct(list)
    if ~isempty(cellids)
        PlotCells(list,cellids,'plot',plottype);
    else
        PlotCells(list,[1:length(list.w)]);
    end
    return;
end

if ischar(list)
    if strfind(list,'.mat')
        subplot(1,2,1);
        load(list);
        SpikeShape(Expt.Spike,'plot');
        title(splitpath(Expt.Header.Name));
        id = findobj('Tag','ReplaySpikes');
        if isempty(id)
            uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['ReplaySpikes(''' Expt.Header.LstName ''')'],...
                'String', 'Replay', 'Position', [10 10 50 20],'Tag','ReplaySpikes');
        else
            set(id, 'Callback', ['ReplaySpikes(''' Expt.Header.LstName ''')']);
        end
    return;
    else
        list = textread(list,'%s');
    end
end
k = 1;

for j = 1:length(list)
    if exist(list{j},'file')
        load(list{j});
        ds = SpikeShape(Expt.Spike);
        if isnan(ds.dprime)
            fprintf('dprime is NaN for %s\n',splitpath(Expt.Header.Name));
        end
        w(k) = ds.width;
        w(k) = ds.iw;
        vw(k) = ds.vw;
        dp(k) = ds.dprime;
        if isfield(Expt.Header,'Area')
            aid = strmatch(Expt.Header.Area,{'V1' 'V2' 'Vd'});
            if ~isempty(aid)
                area(k) = aid;
            else
                area(k) = 0;
            end
        else
            area(k) = 0;
        end
        [barea(k), depth(k)] = GetPenData(list{j});
        names{k} = list{j};
        k = k+1;
    end
end

cells.dp = dp;
cells.vw = vw;
cells.area = area;
cells.names = names;
cells.w = w;
cells.penarea = barea;
cells.depth = depth;

id = find(dp > dpcrit & ~isnan(w) & area == findarea);
PlotCells(cells, id);

function PlotCells(cells, id, varargin)
plottype = 1;
j = 1;
while j <= nargin-2
    if strncmpi(varargin{j},'plottype',4)
        j = j+1;
        plottype = varargin{j};
    end
    j = j+1;
end

if plottype == 1
    histvals = cells.vw(id);
elseif plottype == 2
    histvals = cells.w(id);
end    

subplot(1,2,1);
hold off;
hist(histvals,[1:0.5:max(histvals)]);
bins = hist(histvals,[1:0.5:max(histvals)]);
hold on;
[y,x] = smhist(histvals,'smoother',0.1);
plot(x,y .* max(bins)/max(y),'r');
[dips, pval] = hartigansDipSigniftest(sort(histvals),1000);
title(sprintf('Dip %.3f, p < %.3f, N = %d',dips,pval,length(id)));
subplot(1,2,2);
hold off;
for j = 1:length(id);
    if plottype == 1
        h = plot(cells.w(id(j)),cells.vw(id(j)), 'o', 'buttondownfcn',['SpikeHist(''' cells.names{id(j)} ''')']);
    else
        h = plot(cells.w(id(j)),cells.dp(id(j)), 'o', 'buttondownfcn',['SpikeHist(''' cells.names{id(j)} ''')']);
    end
    hold on;
end
