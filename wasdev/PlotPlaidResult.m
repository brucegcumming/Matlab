function result = PlotPlaidResult(result, varargin)
%PlotPlaidResult(result) replots a result returned PlotExpt for Plaid expts, and
%plots up special verions
% 'xcorr'  plots R(C) vs R(AC) as a scatterplot. If expt type is mixac, will do 
% scatterplots for 0.5 and 0.

plottype = 'counts';
coloroffset = 0;
showcounts = 1;
colors = mycolors;
showors = [];
j = 1;
cvals= [];
bvals= [];
allh = [];
while j <=length(varargin)
    if strncmpi(varargin{j},'coloroff',8)
        j = j+1;
        coloroffset = varargin{j};
    elseif strncmpi(varargin{j},'bvals',4)
        j = j+1;
        bvals = varargin{j};
    elseif strncmpi(varargin{j},'cvals',4)
        j = j+1;
        cvals = varargin{j};
    elseif strncmpi(varargin{j},'ors',3)
        j = j+1;
        showors = varargin{j};
    elseif strncmpi(varargin{j},'plot',4)
        j = j+1;
        plottype = varargin{j};
    end
    j = j+1;
end



res = result;

allh = [];
labels = {};

if strncmp(plottype,'subplot',6)
    nor = unique(res.y(:));
    [nr,nc] = Nsubplots(length(nor));
    for j = 1:length(nor)
        subplot(nr,nc, j);
        hold off;
        PlotPlaidResult(result, 'ors', nor(j));
        title(sprintf('or %.0f',nor(j)));
        axs(j) = gca;
    end
    yl = minmax(cell2mat(get(axs,'Ylim')));
    set(axs,'Ylim',yl);
        
elseif ~isempty(showors)
    for j = 1:length(showors)
        id = find(res.y(1,:,1) == showors(j));
        if ~isempty(id)
            for k = 1:size(res.means,3)
                h = PlotDataLine(res.x(:,id,k),res.means(:,id,k),res.sd(:,id,k),res.n(:,id,k),[],'color',colors{k});
                allh = [allh h(1)];
                zv = res.z(:,:,k);
                labels{length(allh)} = sprintf('%s=%.1f',res.type{3},mean(zv(:)));
            end
        end
    end
    mylegend(allh, labels);
end

xv = unique(res.x(:));

if strcmp(res.type{3},'mixac')
    acs = unique(res.z(:));
    yv = unique(res.y(:));
end


if strcmp(plottype,'xcorr')
    hold off;
    colors = mycolors;
    for j = 1:length(yv)
%corr vs 50% mix        
        plot(res.means(:,j,1),res.means(:,j,2),'o','color',colors{j}, 'markerfacecolor', colors{j});
        hold on;
%corr vs AC        
        plot(res.means(:,j,1),res.means(:,j,3),'o','color',colors{j});
    end
end