function X = PlotLFPpwr(X, varargin);

smw = 0;
colors = mycolors;
tight = 1;
printlabels = 0;
calcresp = 0;
plotargs = {};
pid = [];
explot = [];
j =1;
while j <= length(varargin)
    if strncmpi(varargin{j},'calcresp',4)
        calcresp = 1;
        if length(varargin) > j && ischar(varargin{j+1}) && isfield(X,varargin{j+1})
            j = j+1;
            explot = varargin{j};
        elseif iscell(varargin{j+1});
            j = j+1;
            expts = varargin{j};
            explot = expts{1};
        end
    elseif strncmpi(varargin{j},'ids',3)
        j = j+1;
        pid = varargin{j};
    elseif strncmpi(varargin{j},'print',3)
        printlabels = 1;
    elseif strncmpi(varargin{j},'smooth',3)
        j = j+1;
        smw = varargin{j};
    else
        plotargs = {plotargs{:} varargin{j}};
    end
    j = j+1;
end

if iscell(X)
    All = X{1};
    for j = 2:length(X)
        All.labels = {All.labels{:} X{j}.labels{:}};
        All.n = [All.n  X{j}.n];
        All.resps = cat(3,All.resps, X{j}.resps);
    end
    PlotLFPpwr(All,varargin{:});
    return;
end
if isempty(pid)
    pid = 1:size(X.resps,3);
end
if calcresp
    for j = 1:size(X.resps,3)
        smoothed(:,j) = smooth(X.resps(:,1,j),5,'gauss');
    end
    rvar = std(smoothed,[],2)./mean(smoothed,2);
    [rmax, tmax] = max(rvar);
    x = prctile(rvar,[5 95]);
    crit = x(1) + diff(x) *0.1; %10% of max from baseline
    fid = find(rvar < crit);
    [rmax, tmax] = max(rvar(fid));
    a = find(fid(1:tmax) > crit);
    b = find(fid(tmax:end) > crit);
    fid = fid(a(end))+1:fid(tmax+b(1));
    baseline = min(smoothed(fid,:),[],2);
    for j = 1:size(X.resps,3)
        X.gamma(j) = sum(X.resps(fid,:,j)-baseline);
    end
end
if smw
    for j = 1:size(X.resps,3)
        X.resps(:,1,j) = smooth(X.resps(:,1,j),smw,'gauss');
    end
end
hold off;
baseline = min(smooth(squeeze(X.resps)',4,'gauss')',[],2);
for j  = 1:length(pid)
    c = colors{j};
    plot(X.resps(:,1,pid(j)),'-','color',c,plotargs{:});
    hold on;
end
legend(X.labels(pid));
%plot(baseline,'k');
set(gca,'yscale','log');
if tight
    axis('tight');
end

if printlabels
    for j = 1:length(X.labels)
        fprintf('%d: %s\n',j,X.labels{j});
    end
end


if ~isempty(explot)
    GetFigure('GammaResp');
    hold off;
    if length(expts) > 1
        for j = 2:length(expts)
            ux = unique(X.(expts{j}));
            uvals{j} = ux;
            for k = 1:length(ux);
                ids{j-1}{k} = find(X.(expts{j}) == ux(k));
            end
            nv(j-1) = length(ux);
            cx{j-1} = 0;
        end
        nx = 0;
        allidx = [];
        for j = 1:prod(nv)
            [cx{:}] = ind2sub(nv,j);
            id = ids{1}{cx{1}};
            for k = 2:length(cx);
                id = intersect(id, ids{k}{cx{k}});
            end
            id = setdiff(id, allidx);
            if ~isempty(id)
                nx = nx+1;
                allids{nx} = id;
                allidx = [allidx id];
                labels{nx} = '';
                for k = 1:length(cx)
                labels{nx} = sprintf('%s%s=%.2f ',labels{nx},expts{k+1},uvals{k+1}(cx{k}));
                end
            end
        end
    end
    colors = mycolors;
    for j = 1:length(allids)
        plot(X.(explot)(allids{j}),X.gamma(allids{j}),'o-','color',colors{j},...
            'markerfacecolor',colors{j});
        hold on;
    end
    legend(labels);
end
