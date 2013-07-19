function result = PlotResult(result, varargin)
%PlotResult(result) replots a result returned by PlotExpt()
getcp = 1;
nresample = 0;
npsych = 5;
getem = 0;
holdon = 0;
showseq = 0;
newids = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'newids',6)
        j = j+1;
        newids = varargin{j};
    end
    j = j+1;
end

if isfield(result,'Data')
    Expt = result(1).Data;
    Trials = result(1).Data.Trials;
    [Trials.e0] = deal(0);
elseif iscell(result)
   if isfield(result{1},'plotres')
       PlotExptRes(result{1}.plotres,varargin{:});
   end
   return;
elseif isfield(result,'cp')
else
    if ~isempty(newids)
        result = RecalcResult(result, newids,varargin{:});
    end
    PlotExptRes(result,varargin{:});
    return;
end

%
% need to determine which value of exp is pref, and which value of
% exp2 is ambigous. If maxob is set, this is an expt where large values
% of exp2 are the ambiguous ones. If Maxob is Nan, the smallest value of
% expt2 (usually 0) is ambiguous
maxob = NaN;
minob = NaN;
sigsign = 1;
if strcmp(result(1).type{1}, 'cvsign')
    maxob = 100;
    ctype = 'e0';
    ptype = result(1).type{1};
    if max([Trials.ob]) < maxob
        maxob = max([Trials.ob]) - 10;
    end
    e2vals = [Trials.or];
    e1vals = [Trials.cvsign];
elseif strcmp(result(1).type{2}, 'ob')
    maxob = 100;
    ctype = 'or';
    ptype = result(1).type{2};
    if max([Trials.ob]) < maxob
        maxob = max([Trials.ob]) - 10;
    end
    e2vals = [Trials.ob];
    e1vals = [Trials.or];
    sigsign = -1; %large ob = small signal
elseif strcmp(result(1).type{2},'Tw')
    ctype = 'or';
    maxob = 18;
    if max([Trials.Tw]) < maxob
        maxob = max([Trials.Tw]);
    end
    e2vals = [Trials.Tw];
    e1vals = [Trials.or];
elseif strcmp(result(1).type{2},'Dc')
    ctype = 'dx';
    e2vals = [Trials.Dc];
    e1vals = [Trials.dx];
elseif strcmp(result(1).type{2},'Id')
    ctype = Expt.Stimvals.et;
    e1vals = [Trials.Id];
    e2vals = [Trials.(ctype)];
    minob = 0;
elseif strcmp(result(1).type{2},'dx')
    ctype = Expt.Stimvals.et;
    e1vals = [Trials.e0];
    e2vals = [Trials.(ctype)];
    minob = 0;
    ptype = ctype;
elseif strcmp(result(1).type{1},'Dc')
    ctype = Expt.Stimvals.e2;
    e2vals = [Trials.Dc];
    e1vals = [Trials.(ctype)];
    ptype = ctype;
    minob = min(e2vals);
    sigsign = 1;
else
    e1vals = [Trials.(result(1).type{1})];
    e2vals = [Trials.(result(1).type{2})];
    ctype = Expt.Stimvals.e2;
    ptype = Expt.Stimvals.et;
end

if isnan(minob)
    minob = min(e2vals);
end

if length(result(1).type) < 2 || strcmp(result(1).type{2},'e0')
    pidx = find(e1vals == 0 & [Trials.RespDir] < 0 & [Trials.Trial] > 0);
    nidx = find(e1vals == 0 & [Trials.RespDir] > 0 & [Trials.Trial] > 0);
    goodidx = find(e2vals >0 & [Trials.RespDir] .*  [Trials.rwdir] == 1 & [Trials.Trial] > 0);
    erridx = find(e2vals > 0 & [Trials.RespDir] .* [Trials.rwdir] == -1 & [Trials.Trial] > 0);
elseif isnan(maxob)
    pidx = find(e2vals == minob & [Trials.RespDir] < 0 & [Trials.Trial] > 0);
    nidx = find(e2vals == minob & [Trials.RespDir] > 0 & [Trials.Trial] > 0);
    goodidx = find(e2vals >0 & [Trials.RespDir] .*  [Trials.rwdir] == 1 & [Trials.Trial] > 0);
    erridx = find(e2vals > 0 & [Trials.RespDir] .* [Trials.rwdir] == -1 & [Trials.Trial] > 0);
else
    pidx = find(e2vals >= maxob & [Trials.RespDir] < 0 & [Trials.Trial] > 0);
    nidx = find(e2vals >= maxob & [Trials.RespDir] > 0 & [Trials.Trial] > 0);
    goodidx = find(e2vals < 51 & [Trials.RespDir] .*  [Trials.rwdir] == 1 & [Trials.Trial] > 0);
    erridx = find(e2vals < 51 & [Trials.RespDir] .* [Trials.rwdir] == -1 & [Trials.Trial] > 0);
end
xs = unique([Trials.(ctype)]);
ps = unique([Trials.(ptype)]);


if strcmp(ctype,'e0')
    for j = 1:2;
        if j ==1
        idx = find(e1vals > 0);
        else
        idx = find(e1vals > 0);
        end
        bestrates(j) = mean([Trials(idx).count]);
        rates(j) = mean([Trials(idx).count]);
        choices(j) = mean([Trials(idx).RespDir]);
        vars(j) = var([Trials(idx).count]);
        bestn(j) = length(idx);
        if isfield(Trials,'vs') % hn 12/4/07
            dirs(j) = mean([Trials(idx).vs]);
        elseif isfield(Trials,'sq')
            dirs(j) = mean([Trials(idx).sq]);
        end
    end
else

for j = 1:length(xs)
    if sigsign < 0
        idx = find(e2vals < mean(e2vals) & e1vals == xs(j));
    else
        idx = find(e2vals > mean(e2vals) & e1vals == xs(j));
    end
    rates(j) = mean([Trials(idx).count]);
    choices(j) = mean([Trials(idx).RespDir]);
    if isfield(Trials,'vs') % hn 12/4/07
        dirs(j) = mean([Trials(idx).vs]);
    elseif isfield(Trials,'sq')
        dirs(j) = mean([Trials(idx).sq]);
    end
    idx = find(e2vals == minob & e1vals == xs(j));
    bestrates(j) = mean([Trials(idx).count]);
    vars(j) = var([Trials(idx).count]);
    bestn(j) = length(idx);
end
end
%if rates(1) is highest, pref stim is #1
%if dirs(1) < dirs(end) means this stim was the negative saccade,
%which is Positive respdir.  So swap pidx, nidx. same if both signs
% other way.
dprime = abs(bestrates(end)-bestrates(1))/sqrt((vars(1) * bestn(1) + vars(end) * bestn(end))/(bestn(1)+bestn(end)));
if((rates(1) - rates(end)) * (dirs(1) - dirs(end)) < 0)
    idx = pidx;
    pidx = nidx;
    nidx = idx;
end
if rates(1) > rates(end)
    prefstim = xs(1);
    nullstim = xs(end);
    prefdir = choices(1)-choices(2);
    prefchoice = 1;
    nullchoice = -1;
else
    prefstim = xs(end);
    nullstim = xs(1);
    prefdir = choices(2)-choices(1);
    prefchoice = -1;
    nullchoice = 1;
end

k = 1;
args = {};
if nresample
    args{k} = 'resample';
    args{k+1} = nresample;
    k = k + 2;
end



[cp, details] = CalcCP([Trials(pidx).count],[Trials(nidx).count],args{:});
result(1).cp.cp = cp;
result(1).cp.prefdir = prefdir;
result(1).cp.p = details.pval;
result(1).cp.counts = details.n;
result(1).cp.dprime = dprime;
result(1).cp.prefstim = prefstim;
result(1).cp.nullstim = nullstim;
if ~isnan(details.pval)
    str = sprintf('Pref %.2f CP %.3f p < %.3f, N=',prefstim,cp,details.pval);
else
    str = sprintf('Pref %.2f CP %.3f, N=',prefstim,cp);
end

counts = [Trials.count];
pTrials.sval = [Trials.(ctype)];
pTrials.mag = [Trials.(ptype)];
for j = 1:length(ps)
    for k = 1:length(xs)
        id = find([Trials.(ctype)] == xs(k) & [Trials.(ptype)] == ps(j));
        pTrials.zscores(id) = (counts(id)  - mean(counts(id)))./std(counts(id));
        pp(j,k) = sum([Trials(id).RespDir] > 0)/ sum([Trials(id).RespDir] ~= 0);
        np(j,k) = length(id);
        ip{j,k} = id;
    end
end
ratio = 0.5 - abs(0.5-pp);
gcp.id = find(ratio .* np >= npsych);
tid = [ip{gcp.id}];
pids = find([Trials(tid).RespDir] == prefchoice);
nids = find([Trials(tid).RespDir] == nullchoice);
[gcp.cp, gcp.details] = CalcCP(pTrials.zscores(tid(pids)),pTrials.zscores(tid(nids)),args);
result(1).gcp = gcp;


args = {};
if getem == 3
    args = {args{:} 'emdiff' em 'emskip' emskip};
elseif getem
    args = {args{:} 'em' em 'emskip' emskip};
end
if holdon
    args = {args{:}, 'holdon'};
end
if getcp == 1
    PlotSequence(result,pidx,nidx,'str',str,'correct',goodidx,'wrong',erridx,args{:});
    %        PlotSequence(result,pidx,nidx,'str',str,args{:});
elseif getcp == 3
    PlotSequence(result,pidx,nidx,'str',str,'hist');
else
    [a,b,c,d] = PlotTimeCourse(result,pidx,nidx,'str',str,'sdfw',sdfw,args{:});
    result(1).cp.psdf = a;
    result(1).cp.nsdf = b;
    if getem
        result(1).cp.em = c;
        if ~isempty(d)
            result(1).cp.sample_rate = d.sample_rate;
            result(1).cp.prefsac = d.prefsac;
            result(1).cp.nullsac = d.nullsac;
        end
    end
    result(1).fig = gcf;

end

if ~isempty(result)
    if showseq == 2
        PlotSequence(result,find([Expt.Trials.Trial] > 0),[],'showtimes');
    elseif showseq
        PlotSequence(result,find([Expt.Trials.Trial] > 0),[]);
    end
end
result(1).fig = gcf;

%replot a result from PlotExpt
function PlotSequence(result,fidx, nidx, varargin)

showhist = 0;
showtrials = 1;
showtimes = 0;
holdon = 0;
erridx = [];
goodidx = [];
str = ' ';
    j = 1;
    while j < nargin -2
        if strncmpi(varargin{j},'cp',2)
            j = j+1;
            cp = varargin{j};
            j = j+1;
            cpval = vararargin{j};
        elseif strncmpi(varargin{j},'correct',5)
            goodidx = varargin{j+1};
            j = j+1;
         elseif strncmpi(varargin{j},'showtime',8)
            showtimes = 1;
        elseif strncmpi(varargin{j},'hist',3)
            showhist = 1;
        elseif strncmpi(varargin{j},'holdon',3)
            holdon = 1;
        elseif strncmpi(varargin{j},'str',3)
            j = j+1;
            str = varargin{j};
        elseif strncmpi(varargin{j},'wrong',5)
            erridx = varargin{j+1};
            j = j+1;
        end
        j = j+ 1;
    end
Trials = result(1).Data.Trials;
if isfield(result(1).Data.Header,'BlockStart')
    blocks = result(1).Data.Header.BlockStart;
else
    blocks = [];
end
    
if ~holdon
    hold off;
    
end
    for j = 1:length(Trials)
        starts(j) = Trials(j).Start(1);
    end
    fcounts = [Trials(fidx).count];
    tstart = min(starts)/10000;
    if isfield(Trials,'vs') & mean([Trials(fidx).vs]) < 0
        prefdir = 'down';
        nulldir = 'up';
    else
        prefdir = 'up';
        nulldir = 'down';
    end
    ncounts = [Trials(nidx).count];
    if showhist
        [y,hx] = hist([fcounts ncounts]);
        hist(fcounts,hx);
        hold on;
        if length(fcounts) > 1
        [y,x] = smhist(fcounts,'smoother',0.1);
        pl(1) = plot(x,y,'b','LineWidth',2);
        end
        hold on;
        [y,x] = hist(ncounts,hx);
        h = bar(x,y,1);
        set(h,'FaceColor','none','EdgeColor','r','LineWidth',2);
        if length(ncounts) > 1
        [y,x] = smhist(ncounts,'smoother',0.1);
        pl(2) = plot(x,y,'r','LineWidth',2);
        end
        ylabel('N');
        xlabel('Spike Count');
    else
        if showtimes
            x{1} = starts(fidx)/10000;
            x{2} = starts(nidx)/10000;
            x{3} = [Trials(goodidx).Start]./10000;
            x{4} = [Trials(erridx).Start]./10000;
        else
            x{1} = [Trials(fidx).Trial];
            x{2} = [Trials(nidx).Trial];
            x{3} = [Trials(goodidx).Trial];
            x{4} = [Trials(erridx).Trial];
        end
        if showtrials
            if length(fcounts)
                pl(1) = plot(x{1}(1,:),fcounts,'ob');
                hold on;
            end
            if length(ncounts)
                pl(2) = plot(x{2}(1,:),ncounts,'or');
            end
            ylim = get(gca,'ylim');
            if length(goodidx)
                pl(2) = plot(x{3},ones(size(x{3})).*ylim(2),'ok');
            end
            if length(erridx)
                pl(2) = plot(x{4},zeros(size(x{4})),'ok');
            end
            xlim = get(gca,'xlim');
            for j = 1:length(blocks)
                if blocks(j) < length(Trials) && blocks(j) > 0
                    trial = Trials(blocks(j)).Trial;
                    plot([trial trial],[0 xlim(1)/4],':');
                end
            end
        else
            pl(1) = plot([Trials(fidx).Start]/10000 - tstart,fcounts,'ob');
            hold on;
            pl(2) = plot([Trials(nidx).Start]/10000 - tstart,ncounts,'or');
        end
        if isempty(ncounts)
        else
            [a, isort] = sort(ncounts);
            xlim = get(gca,'xlim');
            plot(xlim(1)+[1:length(ncounts)].* diff(xlim)/length(ncounts),ncounts(isort),'r-');
            [a, isort] = sort(fcounts);
            plot(xlim(1)+[1:length(fcounts)].* diff(xlim)/length(fcounts),fcounts(isort),'b-');
        end
    end
    if sum(ishandle(pl)) == 2 & ~ishold %in case one did not plot
        legend(pl,{sprintf('Pref (%s) Choice',prefdir) sprintf('Null (%s) Choice',nulldir)});
    end
    title(sprintf('%s, %s %d,%d',splitpath(result.name),str,length(fcounts),length(ncounts)));    
    if showtimes
        xlabel(sprintf('Time (sec) - %d',tstart));
    else
         xlabel(sprintf('Trial'));
    end
    ylabel(sprintf('Spike Count (%.2f sec)',result(1).duration./10000));
    
    
   
    
function PlotExptRes(res,varargin)
        
coloroffset = 0;
showcounts = 1;
plottype = 'default';
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
    end
    j = j+1;
end
        
     if iscell(res)
         for j = 1:length(res)
             PlotExptRes(res{j});
         end
         return;
     end
     
     if strcmp(plottype,'acxc')
         PlotACResult(res, 'xcorr');
         return;
     end
     colors = mycolors;
     symbols = 'ososooooooooooooooooooooooooooooo';
     linestyle = {'-',':','-.','--'};
     fillsymbols = 1;
     
      labels = {};  
     if isfield(res,'bestdelay')
         a = PlotRC(res, 'sdf', varargin{:});
         title(sprintf('Cell %d P%.1f',res.cellid,res.Header.probe));
     elseif isfield(res(1),'x2') && res(1).x2 == 0  
        
           for nc = 1:size(res(1).x,3)
       nx = sum(res.n(:,:,nc),1);
       ny = sum(res.n(:,:,nc),2);
       [a,xi] = max(nx);
       [a,yi] = max(ny);
       c = colors{1};
       er = errorbar(res.x(:,xi,nc),res.means(:,xi,nc),res.sd(:,xi,nc) ./ sqrt(max(res.n(:,xi,nc),1)),'o-');
       hold on;
       set(er,'color',colors{1},'linestyle',linestyle{nc});
       allh(nc) = plot(res.x(:,xi,nc),res.means(:,xi,nc),symbols(nc),'color',c);
       res.allhandles = [res.allhandles er];
       if fillsymbols & nc == 1
           set(allh,'MarkerFaceColor',colors{1});
       end
       
       c = colors{2};
       er = errorbar(res.y(yi,:,nc),res.means(yi,:,nc),res.sd(yi,:,nc) ./ sqrt(max(res.n(yi,:,nc),1)),'o-');
       set(er,'color',colors{2},'linestyle',linestyle{nc});
       res.allhandles = [res.allhandles er];
       allh(nc) = plot(res.y(yi,:,nc),res.means(yi,:,nc),symbols(nc),'color',c);
       cpoint = mean(res.x(:));
       labels{nc} = num2str(res.z(1,1,nc));
           set(allh,'MarkerFaceColor',colors{2});
       
           end
                   for j = 1:length(res(1).extras.means)
            if res(1).extras.n(j) > 0
                c =colors{j};
            errorbar(cpoint,res(1).extras.means(j),res(1).extras.sd(j)./sqrt(res(1).extras.n(j)),'color',colors{j});
            hold on;
            h = plot(cpoint,res(1).extras.means(j),'o','color',c);
            set(h,'MarkerFaceColor',c);
            end
                   end
        mylegend(allh,labels);
      else
       cpoint = mean(res(1).x(:));
       if isempty(cvals)
           cvals = 1:size(res(1).x,3);
       end
       if isempty(bvals)
           bvals = 1:size(res(1).x,2);
       end
        for k= cvals;
        for j= bvals;
            if coloroffset
                c = colors{j+coloroffset};
            elseif isfield(res,'colors')
                c = res(1).colors{j,k};
            else
                c = colors{j};
            end
            l = errorbar(res(1).x(:,j,k),res(1).means(:,j,k),res(1).sd(:,j,k)./sqrt(res(1).n(:,j,k)),'color',c);
            if k > 1
                set(l,'linestyle','--');
            end
            hold on;
            dx = range(res(1).x(:,j,k))./50;
            allh(j,k) = plot(res(1).x(:,j,k),res(1).means(:,j,k),'o','color',c,'MarkerFaceColor',c);
            if length(res(1).type) > 1
                labels{j} = sprintf('%s=%.1f',res(1).type{2},res(1).y(1,j,1));
            end
            if showcounts
               h = text(res(1).x(:,j,k)+dx,res(1).means(:,j,k),num2str(res(1).n(:,j,k)),...
                    'horizontalalignment','left','verticalalignment','bottom');
            end
        end
        end
        for j = 1:length(res(1).extras.means)
            c =colors{j};
            if res(1).extras.n(j) > 0
            errorbar(cpoint,res(1).extras.means(j),res(1).extras.sd(j)./sqrt(res(1).extras.n(j)),'color',colors{j});
            h = plot(cpoint,res(1).extras.means(j),'o','color',c);
            set(h,'MarkerFaceColor',c);
            end
        end
     end
     if ~isempty(allh)
      mylegend(allh(:,1),labels);
      if isfield(res,'title')
          title(res(1).title);
      end
     end
      
      
function E = RecalcResult(E, ids, varargin)
% E = RecalcResult(E, ids)
%recalculates a PlotExpt Result using only a subset of Trials
%

dosqrt = 0;
j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'sqrt',4)
        dosqrt = 1;
    end
    j = j+1;
end

for j = 1:length(E.x(:))
    [a, tid] = intersect(E.ids{j},ids);
    if dosqrt == 1
        E.means(j) = mean(sqrt(E.counts{j}(tid)));
        E.sd(j) = std(sqrt(E.counts{j}(tid)));
    else
        E.means(j) = mean(E.counts{j}(tid));
        E.sd(j) = std(E.counts{j}(tid));
    end
    E.n(j) = length(tid);
end
    
for j = 1:length(E.extras.id)
    [a, tid] = intersect(E.extras.id{j},ids);
    if dosqrt == 1
        E.extras.means(j) = mean(sqrt(E.extras.counts{j}(tid)));
        E.extras.sd(j) = std(sqrt(E.extras.counts{j}(tid)));
    else
        E.extras.means(j) = mean(E.extras.counts{j}(tid));
        E.extras.sd(j) = std(E.extras.counts{j}(tid));
    end
    E.extras.n(j) = length(tid);
end
    
    



      
      