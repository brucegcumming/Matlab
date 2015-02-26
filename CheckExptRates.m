function [err, counts] = CheckExptRates(Expt, varargin)
% [err, counts] = CheckExptRates(Expt, varargin)
% Checks an Expt, Set of Expots, or AllExpts struct for suspicious changes
% in rate across blocks. Looks for large blockwise Fano factor (> 5) Coefficient of variaion (> 2) or very
% skewed distribution (-1.5)
% This only checks rates for each cell/cluster separately. To look for
% trials where all cells lack spkies, use PlotRateSequence(AllExpt,'check');
% Given a cell string array, will load each named file ad check that.
%
% to chagne thresholds for reporting errors
% CheckExptRates(Expt,'ffcrit',x) for fano facto
% CheckExptRates(Expt,'skewcrit',x) for skeneww (NB negative number)
% CheckExptRates(Expt,'cvcrit',x) for coefficient of variation.
%
% Also responds to its own ouput
% CheckExptRates(err,'print') prints out suspicous cells, 
% CheckExptRates(err,'plot') plots results for suspicious sequences
% CheckExptRates(err,'plot', 'all') plots results for all
% CheckExptRates(err,'plot','label') adds a label by each point identifying
%                                          file/cell
%   each point is a single expt. Click on the point to bring up the rate
%   sequence plot for that cell. Labels for Expt blocks also provide a
%   context menu for inspecting data/adding comment
%
% For the ORBW Project use
% lst = scanlines('/b/bgc/anal/orbw/allcelllist')
% errs = CheckExptRates(lst)
% CheckExptRates(errs,'plot')


printwarn = 0;
warncolor = 'red';
useall = 0;
crit = [10 -1.5 5 2];

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'print',5)
        printwarn = 1;
    elseif strncmpi(varargin{j},'cvcrit',4)
        j = j+1;
        crit(4) = varargin{j};
    elseif strncmpi(varargin{j},'ffcrit',4)
        j = j+1;
        crit(1) = varargin{j};
        crit(3) = crit(1);
    elseif strncmpi(varargin{j},'skewcrit',6)
        j = j+1;
        crit(2) = varargin{j};
    elseif strncmpi(varargin{j},'plot',4)
        PlotRateCheck(Expt, varargin{j:end});
        return;
    end
    j = j+1;
end

    if ischar(Expt)
        name = Expt;
        Expt = LoadExpt(name);
        [err, counts] = CheckExptRates(Expt, varargin{:});
        return;
    end
    if iscellstr(Expt)
        names = Expt;        
        for j = 1:length(names)
            [err{j}, counts{j}]  = CheckExptRates(names{j}, varargin{:}); 
        end
        return;
    end
    if iscell(Expt)
        if isfield(Expt{1},'err')
            PrintErrors(Expt);
        else
            [err, counts] = CheckExptsRateSequence(Expt, crit, varargin{:});
        end
        return;
    end
    if isfield(Expt,'Spikes') %An AllExpt Structure
        AllExpt = Expt;
        for j = 1:length(AllExpt.Spikes)
            if useall 
                Expt = All2Expt(AllExpt,j,'all');
            elseif AllExpt.Header(j).cellnumber > 0
                Expt = All2Expt(AllExpt,AllExpt.Header(j).cellnumber);
            else
                Expt = [];
            end
            if ~isempty(Expt)
                [X.err{j}, X.counts{j}] = CheckExptRates(Expt,varargin{:});
                [~, eid] = GetExptNumber(Expt);
                nx = 0;
                lastb = 0;
                for b = 1:length(Expt.Header.BlockStartid)
                    if b == length(Expt.Header.BlockStartid)
                        nextb = max([Expt.Trials.id]);
                    else
                        nextb = Expt.Header.BlockStartid(b+1);
                    end
                    bid = find([Expt.Trials.id] >= Expt.Header.BlockStartid(b) & [Expt.Trials.id] < nextb);
                    if ~isempty(bid)
                        nx = nx+1;
                        X.err{j}.blockstart(nx) = bid(1);
                        lastb = Expt.Trials(bid(1)).id;
                        
                        X.err{j}.exptno(nx)  = eid(b);
                        if length(Expt.Header.Clusters) >= nx
                            X.err{j}.probe(nx)  = Expt.Header.Clusters{nx}.probe;
                            X.err{j}.cl(nx)  = Expt.Header.Clusters{nx}.cluster;
                            if isfield(Expt.Header.Clusters(nx),'savetime')
                                X.err{j}.savetime(nx)  = Expt.Header.Clusters{nx}.savetime;
                            end
                        end
                       X.err{j}.CombineDate = Expt.Header.CombineDate;     
                    end
                end
                X.err{j}.cellnumber = AllExpt.Header(j).cellnumber;
                X.name = GetName(Expt);
            end
        end
        err = X;
        counts = X.counts;
        return;
    end
    err.warning = zeros(1,5);
    counts = [Expt.Trials.count];
    if length(counts) > 100
        smw = 10;
    elseif length(counts) > 20
        smw = 5;
    elseif length(counts) == 0
        smw = 2;
    else
        smw = 2;
    end
    idstr = [GetName(Expt) Expt2Name(Expt,'addprobe')];
    

    smc = smooth(counts, smw);
    err.ff = var(smc)./mean(smc);
    b = polyfit(1:length(counts),sqrt(counts),1);
    err.mean = mean(counts);
    err.slope = length(counts).*b(1)./mean(sqrt(counts));
    err.smw = smw;
    if isfield(Expt.Header,'exptno')
        err.exptno = Expt.Header.exptno;
    end
    if sum(counts) == 0
        err.warning(6) = 1;
        mycprintf('blue','E%dCell%d no spikes (%s)\n',err.exptno,Expt.Header.cellnumber,idstr);
        err.slope = 0;
    end
    if err.ff > crit(3)
        err.warning(4) = 1;
        err.msg = sprintf('%s:Fano Factor %d pt avg %.1f\n',idstr,smw,err.ff);
        if printwarn
            mycprintf(warncolor,err.msg);
        end
    end
    if isfield(Expt.Header,'BlockStart')
        trls = [Expt.Trials.Trial];
        for b = 1:length(Expt.Header.BlockStart)
            if b < length(Expt.Header.BlockStart)
                last = Expt.Header.BlockStart(b+1);
            else
                last = max(trls)+1;
            end
            id = find( trls >= Expt.Header.BlockStart(b) & trls < last);
            blkmean(b) = mean([Expt.Trials(id).count]);
            blkn(b) = length(id);
        end
        err.blkcv = sqrt(nanvar(blkmean))./nanmean(blkmean);
        err.blkff = nanvar(blkmean)./nanmean(blkmean);
        err.blkmean = blkmean;
        err.blkn = blkn;
        err.blkskew = skewness(blkmean(~isnan(blkmean)));
        if err.blkff > crit(1)
            err.warning(5) = 1;
            if printwarn
                mycprintf(warncolor,'%s:Block (N%d-%d) Fano Factor %.1f\n',idstr,min(err.blkn),max(err.blkn),err.blkff);
            end
        end
        if err.blkskew < crit(2)
            err.warning(2) = 1;
            if printwarn
                mycprintf(warncolor,'%s:Block (N%d-%d) Skewness %.1f\n',idstr,min(err.blkn),max(err.blkn), err.blkskew);
            end
        end
        if  err.blkcv > crit(4)
             err.warning(5) = 1;
            if printwarn
                mycprintf(warncolor,'%s:Block (N%d-%d) Fano Factor %.1f CV %.2f\n',idstr,min(err.blkn),max(err.blkn), err.blkff,err.blkcv);
            end
        end

    end
    
function [err, counts] = CheckExptsRateSequence(Expts,crit, varargin)
printwarn = 0;
warncolor = 'red';
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'print',5)
        printwarn = 1;
    end
    j = j+1;
end
err.warning = zeros(1,5);
counts = [];
    if iscell(Expts)
        nx = 0;
        for j = 1:length(Expts)
            if isfield(Expts{j},'Trials')
                nx = nx+1;
                [ err.errs{nx}, counts{nx}] = CheckExptRates(Expts{j});
                m(nx) = mean(counts{nx});
                err.blkn(nx) = length(counts{nx});
                idstr = Expt2Name(Expts{j},'addprobe');
            end
        end
        if nx
            err.errs = CellToStruct(err.errs);
        err.blkmean = m;
        id = find(m > 0);
        mrate = mean(cat(2,counts{id}));
        medrate = median(m(m>0));
        sd = std(cat(2,counts{id}));
        err.blkcv = std(m(id))./mean(m(id));
        err.blkff = var(m(id))./mean(m(id));
        err.blkskew = skewness(m); %N.B include all for this
        for j = 1:length(m)
            if m(j) > 0;
             err.diffs(j) = ( err.errs(j).mean)./medrate;
            else
             err.diffs(j) = 1;
            end
        end
        if  err.blkff > crit(1) || err.blkcv > crit(4)
             err.warning(5) = 1;
            if printwarn
                mycprintf(warncolor,'%s:Block (N%d-%d) Fano Factor %.1f CV %.2f\n',idstr,min(err.blkn),max(err.blkn), err.blkff,err.blkcv);
            end
        end
       
        if err.blkskew < crit(2)
            err.warning(2) = 1;
            if printwarn
                mycprintf(warncolor,'%s:Block (N%d-%d) Skewness %.1f\n',idstr,min(err.blkn),max(err.blkn), err.blkskew);
            end
        end
        if sum(m ==0) % should not happen for defined cell
            err.warning(3) = 1;
        end
        end
    end
    
function PrintErrors(errs)



for j = 1:length(errs)
    if isfield(errs{j},'err')
        dpath = BuildFileName(errs{j}.name,'datadir');
        CM = PlotComments(dpath,'loadonly');
        for k = 1:length(errs{j}.err)
            E = errs{j}.err{k};
            if sum(E.warning)
                fprintf('%s Cell %d %d %.2f %.2f   ',errs{j}.name,E.cellnumber,sum(find(E.warning)),E.blkff,E.blkskew);
                ShowComments(0, CM, 'expts', E.exptno,'probes', E.probe);
                fprintf('\n');
            end
        end
    end
end
    
    
    
function PlotRateCheck(errs, varargin)

showall = 0;
showlabel = 0;
marked = [];
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'marked')
        X = varargin{j};
        marked = X.marked;
    elseif strncmpi(varargin{j},'all',3)
        showall = 1;
    elseif strncmpi(varargin{j},'label',5)
        showlabel = 1;
    end
    j = j+1;
end
       

n = 0;
if length(errs) ==1 && isstruct(errs) && isfield(errs,'err')
    x{1} = errs;
    errs = x;
end
X = [];
for j = 1:length(errs)
    if iscell(errs{j})
    for k = 1:length(errs{j})
        E = errs{j}{k};
        if sum(E.warning) || showall
            n = n+1;
            X.meanrate(n) = nanmean(E.blkmean);
            X.blkskew(n) = E.blkskew;
            X.blkff(n) = E.blkff;
            X.ff(n) = E.ff;
            X.name{n} = errs{j}{k}.name;
            X.id(n,:) = [j k];
        end
    end
    elseif isfield(errs{j},'err')
    for k = 1:length(errs{j}.err)
        E = errs{j}.err{k};
        if sum(E.warning) || showall
            n = n+1;
            X.meanrate(n) = nanmean(E.blkmean);
            X.blkskew(n) = E.blkskew;
            X.blkff(n) = E.blkff;
            X.ff(n) = E.ff;
            X.id(n,:) = [j k];
            X.errid = j;
            X.name{n} = errs{j}.name;
            X.cellno(n) = E.cellnumber;
        end
    end
    end
end
if isempty(X)
    fprintf('No rate warnings. Try ''showall''\n');
else
X.x = X.blkff;
X.y = X.blkskew;
if length(marked) < length(X.x)
    marked(length(X.x)) = 0;
end
X.marked = marked;

GetFigure('RateCheck');
hold off; 
colors = mycolors('white');
h = plot(X.x,X.y,'o','buttondownfcn',{@HitPlot},'color',colors{1});
if showlabel
    for j = 1:length(X.x)
        t = text(X.x(j),X.y(j),[' ' X.name{j} '.C' num2str(X.cellno(j))]);
        set(t, 'buttondownfcn',{@HitPlot});
    end
end
cmenu = uicontextmenu;
uimenu(cmenu,'label','Mark Fixed','Callback',{@HitPlot, 'fixed'});
uimenu(cmenu,'label','Mark Come back','Callback',{@HitPlot, 'return'});
uimenu(cmenu,'label','UnMark','Callback',{@HitPlot, 'unmark'});
set(h,'uicontextmenu',cmenu);
if sum(marked)
    hold on;
    id = find(marked ==1)
    h = plot(X.x(id),X.y(id),'o','buttondownfcn',{@HitPlot},'color',colors{2});    
    set(h,'uicontextmenu',cmenu);
    id = find(marked ==2)
    h = plot(X.x(id),X.y(id),'o','buttondownfcn',{@HitPlot},'color',colors{3});    
    set(h,'uicontextmenu',cmenu);
end

setappdata(gcf,'ScatterData',X);
setappdata(gcf,'RateCheckData',errs);
end

function HitPlot(a,b, varargin)
pos = get(gca,'currentpoint');

alt = get(gcf,'SelectionType');
X = getappdata(gcf,'ScatterData');
errs = getappdata(gcf,'RateCheckData');
d = (pos(1,1) + i .* pos(1,2)) - (X.x + X.y .* i);
[a,b] = min(abs(d));


markstrs = {'unmark' 'fixed' 'return'};
j = 1;
while j <= length(varargin)
    if sum(strcmp(varargin{j},markstrs))
        X.marked(b) = find(strcmp(varargin{j},markstrs))-1;
        PlotRateCheck(errs, X);
        return;
    end
    j = j+1;
end
if strcmp(alt,'normal');
F = GetFigure('BlockMeans','parent',gcf);
AddMenus(F,'BlockMeans');
hold off;

id = X.id(b,:);
if iscell(errs{id(1)});
E = errs{id(1)}{id(2)};
else
    X = errs{id(1)};
E = X.err{id(2)};
end

dpath = BuildFileName(X.name,'datadir');
CM = PlotComments(dpath,'loadonly');

plot(X.counts{id(2)});
str = [];
if isfield(E,'CombineDate')
    str = datestr(E.CombineDate);
end
title(sprintf('%sCell%d%s',X.name,E.cellnumber,str));
    yl = get(gca,'ylim');
    if ~isfield(E,'cl')
        E.cl = zeros(size(E.probe));
    end
for j = 1:length(E.blockstart)
    h = line([E.blockstart(j) E.blockstart(j)],yl);
    set(h,'linestyle',':','color','k');
    h = text(E.blockstart(j),yl(2),sprintf('E%dP%d.%d',E.exptno(j),E.probe(j),E.cl(j)),'verticalalignment','top');
    id(3) = j;
    ShowComments(h, CM, 'expts', E.exptno(j),'probes', E.probe(j));    
    cmenu = uicontextmenu;
    uimenu(cmenu,'label','AllVPcs (quick)','Callback',{@ShowData, 'allvpcs', id});
    uimenu(cmenu,'label','AllVPcs (FullV)','Callback',{@ShowData, 'fullv', id});
    uimenu(cmenu,'label','XYPlot','Callback',{@ShowData, 'ClusterXY', id});
    uimenu(cmenu,'label','Check Backups','Callback',{@ShowData, 'listbackup', id});
    uimenu(cmenu,'label','Clear From CellList','Callback',{@ShowData, 'clearcell', id});
    uimenu(cmenu,'label','Add Comment','Callback',{@ShowData, 'comment', id});
    set(h,'uicontextmenu',cmenu');
end
setappdata(gcf,'CurrentErr',id);
end

function ShowData(a,b,type, xid);


F = GetParentFigure(a);
errs = getappdata(F,'RateCheckData');

if nargin < 4
    xid = getappdata(GetFigure(a),'CurrentErr');
end
X = errs{xid(1)};
E = X.err{xid(2)};
exptno = E.exptno(xid(3));
p = E.probe(xid(3));
dpath = BuildFileName(X.name,'datadir');
if strcmp(type,'ClusterXY')
    C = LoadCluster(dpath,exptno,'getxy');
    GetFigure('XYPlot');
    hold off;
    PC.PlotClusterPoints(C{p},[],E.cl(xid(3)));
elseif strcmp(type,'comment')
    PlotComments(dpath,'expt',exptno,'probe',p,'program','CheckExptRates');
elseif strcmp(type,'allvpcs')
    name = sprintf('%s/Expt%dSpikes.mat',dpath,exptno);
    AllVPcs(name,'tchan',p,'reapply');
elseif strcmp(type,'fullv')
    name = sprintf('%s/Expt%dFullV.mat',dpath,exptno);
    AllVPcs(name,'tchan',p,'reapply');
elseif strcmp(type,'clearcell')
    celllist.clearcell(dpath, exptno, E.cellnumber);
elseif strcmp(type,'listbackup')
    ListClusterBackup(dpath,'expts',exptno);
elseif strcmp(type,'meanspikes')
    [C, details] = LoadCluster(dpath,E.exptno,'nodetails');
    GetFigure('MeanSpike');
    colors = mycolors;
    hold off;
    for j = 1:length(C)
        p = E.probe(j);
        if size(C{j}{p}.MeanSpike.ms,1) >= p
            plot(C{j}{p}.MeanSpike.ms(p,:),'color',colors{j});
        else
            plot(C{j}{p}.MeanSpike.ms,'color',colors{j});
        end
        hold on;
        labels{j} = sprintf('E%dP%d',E.exptno(j),E.probe(j));
    end
    legend(labels);
end

function ShowComments(h, CM, varargin);

j = 1; 
expts = [];
probes = [];
for j = 1:length(varargin)
    if strncmpi(varargin{j},'expts',5)
        j = j+1;
        expts = varargin{j};
    elseif strncmpi(varargin{j},'probes',5)
        j = j+1;
        probes = varargin{j};
    end
    j = j+1;
end

if ~isfield(CM,'Comments') || isempty(CM.Comments)
    return;
end

for j = 1:length(CM.Comments)
    cexpts(j) = CM.Comments(j).ex(1);
    cprobes(j) = CM.Comments(j).p;
end

if isempty(expts)
    expts = unique(cexpts);
end

if isempty(probes)
    probes = unique(cprobes);
end

id = find(ismember(cprobes,probes) & ismember(cexpts,expts));
for j = 1:length(id)
    C = CM.Comments(id(1));
    fprintf('E%dP%d %s\n',C.ex,C.p,C.comment); 
    if h > 0
        x = get(h,'extent');
        h = text(x(1),x(2)-x(4),C.comment);
    end
end

function AddMenus(F, type)
if strcmp(type,'BlockMeans')
    id = findobj(F,'type','uimenu','tag','PlotMenu');
    if isempty(id)
        hm = uimenu(F,'label','Plots','tag','PlotMenu');
        sm = uimenu(hm, 'label', 'MeanSpikes', 'callback', {@ShowData, 'meanspikes'});
    end
end
