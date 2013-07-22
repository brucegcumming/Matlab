function [cps, means, sds, xcs, all] = cpsim(pools, varargin)                       
%[cps, means, sds] = cpsim(pools, varargin)                       
%
%cpsim quick simulation of CP only, with no stim noise.
% pools is a vector of pool sizes

bothcorr = 0;
ntrials = 100;
ntotal = 2000;
r = 0.15;
rb = 0.12;
all = [];
CMatrix = [];
rs = ones(size(pools)) .* r;
makefigure = 0;
Wi = [];
j = 1;
while j <= nargin -1
    if strncmpi(varargin{j},'corrs',5)
        j = j+1;
        rs = varargin{j};
    elseif strncmpi(varargin{j},'corr',3)
        j = j+1;
        r = varargin{j}(1);
        rs = ones(size(pools)) .* r;
        if size(varargin{j},1) > 1 && size(varargin{j},2) == size(varargin{j},1)
            bothcorr = 3;
            CMatrix = varargin{j};
        elseif length(varargin{j}) > 1
            rb = varargin{j}(2);
            bothcorr = 2;
            if rb == r
                bothcorr = 1;
            end
            if rb == 0
                bothcorr = 0;
            end
        end
    elseif strncmpi(varargin{j},'ncp',3)
        j = j+1;
        ntotal = varargin{j};
        if(ntotal < max(pools))
            ntotal = max(pools);
        end
    elseif strncmpi(varargin{j},'weights',5)
        j = j+1;
        Wi = varargin{j};
    elseif strncmpi(varargin{j},'ntrials',3)
        j = j+1;
        ntrials = varargin{j};
    elseif strncmpi(varargin{j},'plotresps',6)
        makefigure = 1;
    end
    j = j+1;
end


for pj = 1:length(pools)
    poolsize = pools(pj);
    if length(Wi) ~= poolsize
        Wi = ones(1, poolsize);
    end
    nloops = round(ntotal/(poolsize.*2));
    presps = [];
    nresps = [];
    for k = 0:nloops-1
    if (bothcorr == 1)
        allresps = corr_counts(r,ntrials,2*poolsize);
        presps = allresps(1:poolsize,:);
        nresps = allresps(poolsize+1:2*poolsize,:);
    elseif (bothcorr == 2)
        [presps, nresps] = twopoolcounts(r,rb,ntrials,poolsize);
    elseif (bothcorr == 3)
        np = size(CMatrix,1); %number of subgroups
        ps = poolsize * 2./np; %number of cells in each subgroup
        for j = 1:np
            for m = 1:np
                starta = 1 + ps * (j-1);
                startb = 1 + ps * (m-1);
                CM(starta:starta+ps-1,startb:startb+ps-1) = CMatrix(j,m);
            end
        end
        for j = 1:size(CM,1)
            CM(j,j) = 1;
        end
        [presps, nresps] = twopoolcounts(r,rb,ntrials,poolsize,'CMatrix',CM);
    else
        presps = corr_counts(rs(pj),ntrials,poolsize);
        nresps = corr_counts(rs(pj),ntrials,poolsize);
    end
    
    psum = (Wi * presps)./poolsize;
    nsum = (Wi * nresps)./poolsize;
    dvar = psum-nsum;
    choices =  dvar > 0;
    for j = 1:poolsize
        cps(pj,j+k*poolsize*2) = CalcCP(presps(j,find(choices)),presps(j,find(~choices)));
        cps(pj,j+poolsize+(k*poolsize*2)) = CalcCP(nresps(j,find(~choices)),nresps(j,find(choices)));
        xc = corrcoef(presps(j,:),dvar);
        xcs(pj,j+k*poolsize*2) = xc(1,2);
        xc = corrcoef(nresps(j,:),-dvar);
        xcs(pj,j+poolsize+(k*poolsize*2)) = xc(1,2);
    end
    end
    means(pj) = mean(cps(pj,:));
    sds(pj) = std(cps(pj,:));
end

all.presps =  presps;
all.nresps = nresps;
all.choices = choices;
colors = mycolors;

if makefigure ==1;
    colors{1} = [0.8 0.2 0];
    colors{2} = [0 0.5 0 ];
    colors{3} = [0.2 0.2 1];
    colors{4} = [1 0 1];
    colors{5} = [0 1 1];
    base = 50;
   o = [1:poolsize] .* 3.*std(presps(:));
   ax = subplot(2,1,1);
   hold off;
   c = 1+mod([1:poolsize]-1,length(colors));
   for j = 1:poolsize
       plot(presps(j,:)+o(j),'o','color',colors{c(j)},'markerfacecolor',colors{c(j)});
       hold on;
       plot([1 ntrials],[base+o(j) base+o(j)],'-','color',colors{c(j)});
       text(0.5,base+o(j),sprintf('Cell %d',j),'horizontalalignment','right','color',colors{c(j)},'fontsize',24);
   end
   axis('tight');
   title('Pool ''A''');
   set(gca,'Yticklabel',[]);
   set(gca,'Xticklabel',[]);
   yl = get(gca,'ylim');
   y = yl(1) - diff(yl)/5;
   for j = 1:size(presps,2)
       if choices(j)
           text(j,y,'A','fontsize',24);
       else
           text(j,y,'B','fontsize',24);
       end
   end
   for j = 1:size(presps,1)
   end
   text(21,y,'Choice','horizontalalignment','left','fontsize',24);
   subplot(2,1,2);
   hold off;

   for j = 1:poolsize
       plot(nresps(j,:)+o(j),'o','color',colors{c(j)},'markerfacecolor',colors{c(j)});
       hold on;
       plot([1 ntrials],[base+o(j) base+o(j)],'-','color',colors{c(j)});
       text(0.5,base+o(j),sprintf('Cell %d',j),'horizontalalignment','right','color',colors{c(j)},'fontsize',24);
   end
   axis('tight');
   title('Pool ''B''');
   set(gca,'Yticklabel',[]);
   xlabel('Trial Number');
   figify(gcf,gca);
   figify(gcf,ax);
   
end
   
