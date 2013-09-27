function h = PlotND(X,plots, varargin)
%PlotND(X,plots, varargin)
%takes an MxN matrix and generates scatterplots X(:,1) vs X(:,2) etc
% plots is a 2xnp matrix where each row specifies 2 dimensions of X to
% plot.
%PlotND(X,plots, 'idlist', id) where id is a vector of length X(:,1)
%plots different colors according to the value of id
%    ...,'markersize',x)   sets markersize
%    ...,'marker',x)   sets marker type
idlist = [];
nid = 1;
colors = mycolors;
densityplot = 0;
marker = '.';
markersize = 1;
callback = [];
if length(X) < 1000
    markersize = 4;
end
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'colors',5)
        j = j+1;
        colors = varargin{j};
    elseif strncmpi(varargin{j},'callback',5)
        j = j+1;
        callback = varargin{j};
        j = j+1;
        callbackid = varargin{j};
    elseif strncmpi(varargin{j},'density',5)
        densityplot = 1;
    elseif strncmpi(varargin{j},'idlist',5)
        j = j+1;
        idlist = varargin{j};
        cids = unique(idlist);
        nid = length(cids);
    elseif strncmpi(varargin{j},'markersize',8)
        j = j+1;
        markersize = varargin{j};
    elseif strncmpi(varargin{j},'marker',5)
        j = j+1;
        marker = varargin{j};
    end
    j = j+1;
end

if nargin ==1 || isempty(plots)  %automatic
    nr=2;
    nc=4;
    if min(size(X)) == 2
        plots = [1 2];
        nplots = 1;
        nr=1;
        nc=1;
    elseif min(size(X)) == 3
        plots = [1 2; 1 3; 2 3;]
        nplots = 3;
        nr=2;
        nc=2;
    elseif min(size(X)) == 4
        plots = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4;]
        nplots = 6;
        nr=2;
        nc=2;
    elseif min(size(X)) == 5
        plots = [1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4;];
        nplots = 8;
    elseif min(size(X)) > 5
        plots = [1 2; 1 3; 1 4; 1 5; 1 6; 2 2; 2 5; 2 4;];
        nplots = 8;
    end
    
else
    nplots = size(plots,1);
    [nr,nc] = Nsubplots(nplots);
end

if densityplot
    for j = 1:nplots
        subplot(nr,nc,j);
        hold off;
        DensityPlot(X(:,plots(j,1)),X(:,plots(j,2)));
        set(gca,'ydir','normal');
    end
else
    for k = 1:nid
        if isempty(idlist)
            id = 1:size(X,1);
        else
            id = find(idlist == cids(k));
        end
        for j = 1:nplots
            subplot(nr,nc,j);
            if k == 1
                hold off;
            end
            if ~isempty(callback)
                h = myscatter(X(id,plots(j,1)),X(id,plots(j,2)),marker,...
                    'ids',callbackid,...
                    'buttonpress', callback,...
                    'markersize',markersize,'color',colors{k});
            else
                h(j) = plot(X(id,plots(j,1)),X(id,plots(j,2)),marker,'markersize',markersize,'color',colors{k});
            end
            hold on;
            title(sprintf('%d Vs %d',plots(j,1),plots(j,2)));
        end
    end
end
