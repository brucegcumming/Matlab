function result = PlotACResult(result, varargin)
%PlotACResult(result) replots a result returned PlotExpt for AC expts, and
%plots up special verions
% 'xcorr'  plots R(C) vs R(AC) as a scatterplot. If expt type is mixac, will do 
% scatterplots for 0.5 and 0.

plottype = 'counts'
coloroffset = 0;
showcounts = 1;
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
    elseif strncmpi(varargin{j},'xcorr',4)
        plottype = 'xcorr';
    end
    j = j+1;
end



res = result;

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