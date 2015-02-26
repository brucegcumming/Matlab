function h = myscatter(x,y,symbol,varargin)
%myscatter(x,y,symbol, ... makes a scatterplot where each datapoint has a callback
%function. Default is to print out index of data point in the vectors x,y
%
% myscatter(x,y,'.','buttondown',@fcn)
%               sets the function called after a press
% myscatter(x,y,'.','ids',idlist)
%               sets the id numbers that are passed to @fcn
% myscatter(x,y,'.','colors',colors)
%               attaches a label to each point in the callback
% myscatter(x,y,'.','labels',labels)
%               attaches a label to each point in the callback


fcn = @scatterhit;
idlist = [];
bidlist = [];
plotargs = {};
ptlabels = {};
colors = [];
colorids = [];
h = [];
if size(x,2) == 1
    x = x';
end
if size(y,2) == 1
    y = y';
end

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'buttonpress',10)
        j = j+1;
        fcn = varargin{j};
    elseif strncmpi(varargin{j},'ids',3)
        j = j+1;
        idlist = varargin{j};
    elseif strncmpi(varargin{j},'labels',5)
        j = j+1;
        ptlabels = varargin{j};
    elseif strncmpi(varargin{j},'colorids',8)
        j = j+1;
        colorids = varargin{j};
    elseif strncmpi(varargin{j},'colors',6)
        j = j+1;
        colors = varargin{j};
    else
        plotargs = {plotargs{:} varargin{j}};
    end
    j = j+1;
end

if ~isempty(colorids) && isempty(colors)
    colors = mycolors;
end
if isempty(idlist)
    idlist = 1:size(x,2);
elseif size(idlist,2) == 1
    idlist = idlist';
end
if isempty(bidlist)
   bidlist = 1:size(x,1);
end
np = 0;
for k = 1:size(x,1)
for j = 1:size(x,2)
    np = np+1;
    h(np) =  plot(x(k,j),y(k,j),symbol,'buttondownfcn',{fcn, idlist(:,j), bidlist(k)},plotargs{:});
    if ~isempty(colorids)
        set(h(np),'color',colors{colorids(k,j)});
    elseif ~isempty(colors)
        set(h(np),'color',colors{k,j});
    end
    hold on;
end
end
if ~isempty(ptlabels)
    setappdata(gcf,'Labels',ptlabels);
end


function scatterhit(a,b,ida, idb)

pos = get(gca,'currentpoint');
labels = getappdata(gcf,'Labels');
if isempty(labels)
    lbl = '';
elseif ida == 1
    lbl = labels{idb}
else
    lbl = labels{ida}
end
fprintf('Point %s, %d at %s %s %s\n',sprintf('%d,',ida), idb,num2str(pos(1,1)),num2str(pos(1,2)),lbl);

