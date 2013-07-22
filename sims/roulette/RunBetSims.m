function [stopt, result, res] = RunBetSims(edges, nruns, varargin)
runl = 10000;
maxls = [1000 10000];
li = 2;
maxl = 1000;
type = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'len',3)
        runl = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'ids',3)
        j= j+1;
        
        ids = varargin{j};
        ei = ids(1);
        li = ids(2);
    elseif strncmpi(varargin{j},'maxl',4)
        maxls = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'anal',4)
       result = edges;
       res = nruns;
       if length(varargin) > j & isnumeric(varargin{j+1})
           j = j+1;
           type = varargin{j};
       else
           type = 1;
       end
       [stopt,result] = AnalyseRuns(res, result, maxls,ei,li, type);
       return;
    end
    j = j+1;
end


for j = 1:length(edges)
    [a, res{j}] = PropBet(edges(j), nruns,'len',runl,'maxl',maxl);
    result.xids{j} = res{j}.xids;
    res{j}.pots = a;
end
result.edges = edges;
result.runl = runl;
result.nruns = nruns;
ei = length(edges);
li = length(maxls);
[stopt,result] = AnalyseRuns(res, result, maxls,ei, li, type);

function [stopt, result] = AnalyseRuns(res, result, maxls,ei,li, type)

plottype = type;
edges = result.edges;
runl = result.runl;
nruns = result.nruns;
totallens = [50000 100000 500000 1000000];
stopt = zeros(size(edges,2), size(maxls,2), size(totallens,2));
pots = zeros(size(edges,2),size(maxls,2),runl);
for l = 1:length(totallens)
for j = 1:length(edges)
    for k = 1:nruns  %% find gap of > runl, with positive sign at end
        for m = 1:length(maxls)
            if length(res{j}.xids{k}) > 2
                id = find(diff(abs(res{j}.xids{k})) > maxls(m) & res{j}.xids{k}(2:end) > 0);
%negative xid means going negative. So if only one is negative, would have
%stopped
            elseif length(res{j}.xids{k}) == 1 & res{j}.xids{k}(1) < 0
                id = 1;
            else
                id = [];
            end
            id= id(find(abs(res{j}.xids{k}(id)) <= totallens(l)));
            if length(id)
                stoptimes(j,m,l,k) = abs(res{j}.xids{k}(id(1)));
                stopt(j,m,l) = stopt(j,m,l)+1;
                result.pots(j,m,k) = res{j}.allstoppots{k}(id(1));
                result.nbets(j,m,k) = res{j}.nbets{k}(id(1));
%                result.nbets(j,m,k) = res{j}.xids{k}(id(1));
            else
                stoptimes(j,m,l,k) = totallens(l);
                if isempty(res{j}.nbets{k})
                result.nbets(j,m,k) = 0;
                else
                    result.nbets(j,m,k) = res{j}.nbets{k}(end);
                end
                result.pots(j,m,k) = res{j}.pots(k);
            end
        end
    end
end
end
j = 1;

if plottype == 3
    for m = 1:length(totallens)
        subplot(2,2,m);
        hold off;
        plot(result.edges,squeeze(stopt(:,:,m)));
        title(sprintf('Len %d',totallens(m)));
    end
    fid = fopen('stops.csv','w');
    fprintf(fid,'"Egde"')
    for j = 1:length(maxls)
        fprintf(fid,',"Stop %d"',maxls(j));
    end
    fprintf(fid,'\n');
    for j = 1:length(result.edges)
        fprintf(fid,'%.2f',result.edges(j));
        fprintf(fid,', %.4f',squeeze(stopt(j,:,end))./1000);        
        fprintf(fid,'\n');
    end
    fclose(fid);
elseif plottype == 4
    colors = mycolors;
    for m = 1:length(totallens)
        subplot(2,2,m);
        hold off;
        for j = 1:length(maxls)
            [y,x] = smhist(squeeze(stoptimes(1,j,stm,:)));
            plot(x,y,'color',colors{j});
            hold on;
        end
        title(sprintf('Len %d',totallens(m)));
    end
elseif plottype == 2
    hist(result.pots(1,1,:));
    title(sprintf('Edge %.2f',result.edges(1)));
else
subplot(2,2,1);
%imagesc(stopt);
plot(squeeze(stopt(:,:,1)));
subplot(2,2,2);
plot(squeeze(stopt(:,:,2)));
title(sprintf('Run Length %d',totallens(2)));

subplot(2,2,3);
hist(result.pots(end,1,:));
title(sprintf('Edge %.2f',result.edges(end))); 
subplot(2,2,4);
hold off;
plot(squeeze(result.nbets(ei,li,:)),squeeze(result.pots(ei,li,:)),'o');
title(sprintf('Edge %.2f, Ncrit %.2f',result.edges(ei),maxls(li))); 

end