function res = Nsinesdf(Expt, im, varargin)

%res = Nsinesdf(Expt, im, varargin)
%
fz = 10000/Expt.Stimvals.fz;
fz = fz/2;
ntrig = 0;
sdfw = 104;
plottype = 0;
times = [0:10:1000];
sfid = [1:size(im.sfs,2)];
legendpos = 1;
triggers = [];
rtrig = [];
ltrig = [];
res = [];
smean = [];
if isfield(im,'sfvals')
    sfvals = im.sfvals;
else
    sfvals = [];
end
labels = {};
sdfargs = {'clip'};
holdon = 0;
showim = 0;
MONOCS = 7;
CORRPLOT = 8;
corrs = [];
cth = [];
seedstep = 2;

j = 1;
while j < nargin - 1;
    if strncmpi(varargin{j},'all',3)
        plottype = 1;
    elseif strncmpi(varargin{j},'box',3)
        j = j+1;
        sdfargs = {sdfargs{:} 'box'};
        sdfw = varargin{j};
    elseif strncmpi(varargin{j},'corrth',6)
        j = j+1;
        for k = 1:size(corrs,2)
            cth(:,k) = varargin{j};
        end
        res.corrth = cth;
    elseif strncmpi(varargin{j},'corrs',3)
        j = j+1;
        corrs = varargin{j};
        
        if isempty(cth)
            for k = 1:size(corrs,2)
                cth(:,k) = prctile(corrs(:,k),[25 75]);
            end
            res.corrth = cth;
        end
        plottype = CORRPLOT;
    elseif strncmpi(varargin{j},'sfdp',4)
        plottype = 6;
        if length(varargin) > j  & isnumeric(varargin{j+1})
            j = j+1;
            dpval = varargin{j};
        else
            dpval = im.dpvals(1);
        end
    elseif strncmpi(varargin{j},'dp',2)
        if strncmpi(varargin{j},'dpa',3)
            plottype = 4;
        else
            plottype = 3;
        end
        if length(varargin) > j  & isnumeric(varargin{j+1})
            j = j+1;
            sfid = varargin{j};
        else
            sfid = 4;
        end
    elseif strncmpi(varargin{j},'diffs',3)
        plottype = 2;
    elseif strncmpi(varargin{j},'hold',3)
        holdon = 1;
    elseif strncmpi(varargin{j},'legend',3)
        j = j+1;
        legendpos = varargin{j};
    elseif strncmpi(varargin{j},'image',3)
        showim = 1;
    elseif strncmpi(varargin{j},'monocs',3)
        plottype = MONOCS;
    elseif strncmpi(varargin{j},'plot',3)
        j = j+1;
        plottype = varargin{j};
    elseif strncmpi(varargin{j},'rc',2)
        j = j+1;
        nsres = varargin{j};
    elseif strncmpi(varargin{j},'sdfw',3)
        j = j+1;
        sdfw = varargin{j};
    elseif strncmpi(varargin{j},'sfid',3)
        j = j+1;
        sfid = varargin{j};
    elseif strncmpi(varargin{j},'times',3)
        j = j+1;
        times = varargin{j};
    elseif strncmpi(varargin{j},'triggers',3)
        j = j+1;
        triggers = varargin{j};
    elseif strncmpi(varargin{j},'vals',3)
        j = j+1;
        sfvals = varargin{j};
    end
    j = j+1;
end


if isempty(triggers)
    for j = length(Expt.Trials):-1:1
        nseeds = Expt.Trials(j).Nf * seedstep;
        ds = Expt.Trials(j).ls - im.seeds; 
        id = find(ds >= 0 & ds <= nseeds);
        ds = ds(id); %% seeds that were used in this trial)
        if ~isempty(corrs)
            for k = 1:size(corrs,2)
                ida  = find(corrs(id,k) > cth(2,k));
                if isempty(ida)
                    triggers(k).Trial(j).trigger = [];
                else
                    triggers(k).Trial(j).trigger = fz * (nseeds - ds(ida))';
                end
                ida = find(corrs(id,k) < cth(1,k));
                if isempty(ida)
                    atriggers(k).Trial(j).trigger = [];
                else
                    atriggers(k).Trial(j).trigger = fz * (nseeds - ds(ida))';
                end
                ida = find(corrs(id,k) >= cth(1,k) & corrs(id,k) <= cth(2,k));
                if isempty(ida)
                    ztriggers(k).Trial(j).trigger = [];
                else
                    ztriggers(k).Trial(j).trigger = fz * (nseeds - ds(ida))';
                end
            end
        elseif ismember(plottype,[3 4])
            for k = 1:length(im.dpvals)
                ida  = find(im.sfs(id,sfid) > 0 & im.dp(id,sfid) == im.dpvals(k));
                triggers(k).Trial(j).trigger = fz * (nseeds - ds(ida))';
            end
        elseif ismember(plottype,[9 10 11 12])            
            if plottype == 9
                lc = 0.5;
                rc = 0.5;
            elseif plottype == 10
                lc = 1;
                rc = 1;
            elseif plottype == 11
                lc = 0.5;
                rc = 1;
            elseif plottype == 12
                lc = 1;
                rc = 0.5;
            end
            for k = 1:size(im.sfs,2)
                ida  = find(im.rsfs(id,k) == rc & im.lsfs(id,k) == lc & im.dp(id,k) == dpval);
                triggers(k).Trial(j).trigger = fz * (nseeds - ds(ida))';
            end
        elseif ismember(plottype,[5 6])
            for k = 1:size(im.sfs,2)
                if isnan(dpval)
                    ida  = find(im.sfs(id,k) > 0);
                else
                    ida  = find(im.sfs(id,k) > 0 & im.dp(id,k) == dpval);
                end
                triggers(k).Trial(j).trigger = fz * (nseeds - ds(ida))';
            end
        else
            for k = 1:size(im.sfs,2)
                ida  = find(im.sfs(id,k) > 0);
                triggers(k).Trial(j).trigger = fz * (nseeds - ds(ida))';
            end
        end
        if isfield(im,'rsfs')
            for k = 1:size(im.sfs,2)
                ida  = find(im.rsfs(id,k) > 0);
                rtrig(k).Trial(j).trigger = fz * (nseeds - ds(ida))';
            end
        end
        if isfield(im,'lsfs')
            for k = 1:size(im.sfs,2)
                ida  = find(im.lsfs(id,k) > 0);
                ltrig(k).Trial(j).trigger = fz * (nseeds - ds(ida))';
            end
        end
        alltriggers.Trial(j).trigger = fz * 2 * [0:200];
    end
    if plottype == 6
        res.triggers = triggers;
    else
        res.triggers = alltriggers;
    end
end

if plottype == MONOCS
    for j = 1:length(ltrig)
        [res.lsdf(j,:), nt, nsp, spk, ptime] = trigsdfa(Expt.Trials,sdfw,times,'trig',ltrig(j).Trial,sdfargs{:});
        res.n(j,1) = nt;
    end
    for j = 1:length(rtrig)
        [res.rsdf(j,:), nt, nsp, spk, ptime] = trigsdfa(Expt.Trials,sdfw,times,'trig',rtrig(j).Trial,sdfargs{:});
        res.n(j,2) = nt;
        res.times = ptime;
    end
    return;
end

if ~holdon
    hold off;
end
if plottype == 1 %plot sdfs
    colors = mycolors;
    for j = 1:length(triggers)
        [sdf(j,:), nt, nsp, spk, ptime] = trigsdfa(Expt.Trials,sdfw,times,'trig',triggers(j).Trial,sdfargs{:});
        plot(ptime./10, sdf(j,:),'color', colors{j});
        if ~isempty(sfvals)
            labels{j} = sprintf('%.1f',sfvals(j));
        end
        hold on;
    end
    smean = mean(sdf);
    sdf(j+1,:) = smean;
    plot(ptime./10, smean,'k', 'linewidth',2);
elseif ~isempty(corrs)
    colors = mycolors;
    for j = 1:length(triggers)
        [sdf(j,:), nt, nsp, spk, ptime] = trigsdfa(Expt.Trials,sdfw,times,'trig',triggers(j).Trial,sdfargs{:});
        [asdf(j,:), nt, nsp, spk, ptime] = trigsdfa(Expt.Trials,sdfw,times,'trig',atriggers(j).Trial,sdfargs{:});
        [zsdf(j,:), nt, nsp, spk, ptime] = trigsdfa(Expt.Trials,sdfw,times,'trig',ztriggers(j).Trial,sdfargs{:});
    end    
    imagesc(sdf(40:end-40,:));
    res.asdfs = asdf;
    res.zsdfs = zsdf;
elseif plottype == 2 %%plot sdfs - mean;
    colors = mycolors;
    hold off;
    for j = 1:length(triggers)
        [sdf(j,:), nt, nsp, spk, ptime] = trigsdfa(Expt.Trials,sdfw,times,'trig',triggers(j).Trial,sdfargs{:});
        res.n(j) = nt;
        nspikes(j) = nsp;
    end
    smean = mean(sdf);
    for j = 1:length(sfid)
        plot(ptime./10, sdf(sfid(j),:)-smean,'color', colors{j});
        hold on;
        if ~isempty(sfvals)
            labels{j} = sprintf('%.1f',sfvals(sfid(j)));
        end
    end
elseif ismember(plottype,[3 4 5 6 9 10 11 12]) %%line for each dp, single SF
    tic;
    colors = mycolors;
    hold off;
    for j = length(triggers):-1:1
        [sdf(j,:), nt, nsp, spk, ptime] = trigsdfa(Expt.Trials,sdfw,times,'trig',triggers(j).Trial,sdfargs{:});
        res.n(j) = nt;
    end
    smean = mean(sdf);
    if ismember(plottype,[4 6])
        [smean, nt, nsp, spk, ptime] = trigsdfa(Expt.Trials,sdfw,times,'trig',alltriggers.Trial,sdfargs{:});
        smean = smean';
    else
        smean = mean(sdf);
    end
    res.buildtime = toc;
    for j = 1:length(triggers)
        plot(ptime./10, sdf(j,:)-smean,'color', colors{j});
        hold on;
        if j <= length(sfvals)  %% canhapppen if n dp val s> sfs
            labels{j} = sprintf('%.1f',sfvals(j));
        end
    end
    plot(ptime./10,smean-mean(smean),'k','linewidth',2);
end
figure(gcf);

res.sdfs = sdf;
res.times = ptime;
res.mean = smean;

if ~isempty(labels) & legendpos
    legend(labels);
end

if showim
    hold off;
    imagesc(ptime./10,sfvals,sdf);
%    colorbar;
end

if ismember(plottype,[9 10 11 12])
    title(sprintf('L%.1f,R%.1f',lc,rc));
elseif ismember(plottype,[3 4]) %%line for each dp, single SF
    title(sprintf('SF=%.2f',sfvals(sfid)));
elseif ismember(plottype,[5 6]) %%line for each dp, single SF
    title(sprintf('dp=%.2f',dpval));
end