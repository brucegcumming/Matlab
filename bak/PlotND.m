function PlotND(X,plots, varargin)
%PlotND(X,plots, varargin)
%takes an MxN matrix and generates scatterplots X(:,1) vs X(:,2) etc
% plots is a 2xnp matrix where each row specifies 2 dimensions of X to
% plot.
%PlotND(X,plots, 'idlist', id) where id is a vector of length X(:,1)
%plots different colors according to the value of id
idlist = [];
nid = 1;
colors = mycolors;
densityplot = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'colors',5)
        j = j+1;
        colors = varargin{j};
    elseif strncmpi(varargin{j},'density',5)
        densityplot = 1;
    elseif strncmpi(varargin{j},'idlist',5)
        j = j+1;
        idlist = varargin{j};
        nid = length(unique(idlist));
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
    elseif min(size(X)) == 4
        plots = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4;]
        nplots = 6;
    elseif min(size(X)) == 5
        plots = [1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4;];
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
    end
else
    for k = 1:nid
        if isempty(idlist)
            id = 1:size(X,1);
        else
            id = find(idlist == k);
        end
        for j = 1:nplots
            subplot(nr,nc,j);
            if k == 1
                hold off;
            end
            plot(X(id,plots(j,1)),X(id,plots(j,2)),'.','markersize',1,'color',colors{k});
            hold on;
        end
    end
end
