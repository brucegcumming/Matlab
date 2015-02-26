function [res, details] = ExptListCorrs(list, varargin)
%
%ExptListCorrs(list....
%takes a cell array of strings in liest, reads in those files, and
%calculates correlations between any simultaneously recorded cells. 
res = [];
verbose = 0;
args = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'verbose',4)
        verbose = 1;
        args = {args{:} varargin{j}};
    end
    j = j+1;
end


if iscell(list) && ~iscellstr(list)
    [res, details] = CalcCorrs(list, [], varargin);
    return;
end
if isfield(list,'Expt')
    for j = 1:length(list.Spikes)
        Expts{j} = All2Expt(list, j, 'all');
    end
    [res, details] = CalcCorrs(Expts, [], varargin);
    return;
end
for j = 1:length(list)
    suffs{j} = regexprep(list{j}, '.*cell[0-9]*\.','');
    t = regexprep(list{j}, '.*cell([0-9]*)\..*','$1');
    cellst(j) = str2num(t);
end

types = unique(suffs);

details.types = types;
k = 0;
for j = 1:length(types)
    id = strmatch(types{j},suffs);
    cells = cellst(id);
    details.good(j) = 0;
    if length(cells) > 1
        if verbose
            fprintf('Calculating %s,%s....\n',list{id(1)},list{id(2)});
        end
        [a, b] = CalcCorrs(list(id), cells,args{:});
        if ~isempty(a)
            k = k+1;
            res{k}.c = a;
            res{k}.type = types{j};
            res{k}.ztype = b.ztype;
            details.good(j) = 1;
        end
    end
end




function [res, details] = CalcCorrs(list, cells, varargin)

nmin = 3;
j = 1;
verbose = 0;

while j <= length(varargin)
    if strncmpi(varargin{j},'nmin',4)
        j = j+1;
        nmin = varargin{j};
    elseif strncmpi(varargin{j},'verbose',4)
            verbose = 1;
    end
    j = j+1;
end

details.ztype = 0;

ng = 0;
if iscell(list) && ~iscellstr(list)
    Expts = list;
    for j = 1:length(list)
        cells(j) = j;
    end
    ng = length(list);
else
    for j = 1:length(list)
        load(list{j});
        if exist('cExpt','var')
            Expts{j} = cExpt;
            Expts{j}.Header.cellnumber = cells(j);
            ng = ng+1;
        else
            Expts{j} = [];
        end
        Expts{j}.Header.filename = list{j};
    end
end
if ng == 0
    res = [];
    return;
end

et = Expts{1}.Stimvals.et;
e2 = Expts{1}.Stimvals.e2;
isrc = Expts{1}.Header.rc;

condition = {};
if strcmp(et,'or') && strcmp(e2,'ob')
    condition{1} = 'ob] > 120';
elseif strcmp(et,'Dc') || strcmp(e2,'Dc')
    condition{1} = 'Dc] == 0';
elseif isfield(Expts{1}.Trials,'co')
    condition{1} = 'co] > 0';
end


nc = 0;
if length(condition)
    details.ztype = 2;
   for j = 1:length(Expts)
       for k = j:length(Expts)
           [a,b] = ExptCorr(Expts{j},Expts{k},condition);
           if ~isempty(b) & ~isnan(b.rawxc)
           b.xc = a;
           b.cells = [cells(j) cells(k)];
           b.probes = [mean(Expts{j}.Header.probe) mean(Expts{k}.Header.probe)];
           b.state.dprime = [GetExptQuality(Expts{j}) GetExptQuality(Expts{k})];
           nc = nc+1;
           res(nc) = b;
           end
       end
   end
else
    details.ztype = 1;
    times = -3010:10:3010;
   for j = 1:length(Expts)
       seeds = [];
       if isfield(Expts{j}.Trials,'seed')
           for k = 1:length(Expts{j}.Trials)
               seeds(k) = Expts{j}.Trials(k).seed(end);
           end
       end
       for k = j:length(Expts)
           b.ntrials = 0;
           b.rsc = NaN;
           b.cells = [cells(j) cells(k)];
           aid = find(ismember([Expts{j}.Trials.id],[Expts{k}.Trials.id]));
           if length(aid) > nmin
           ares = PlotExpt(Expts{j},'noplot','condense','Trials',aid);
           nx = 0;
           for c = 1:length(ares.counts(:));
               if length(ares.counts{c}) > nmin
                   nx = nx+1;
                   az{nx} = (ares.counts{c} - mean(ares.counts{c}))/std(ares.counts{c});
                   for d = 1:length(ares.ids{c})
                   aspks{nx}{d} = Expts{j}.Trials(ares.tidx{c}(d)).Spikes;
                   end
                   eb(nx,1) = ares.y(nx);
               end
           end
           bid = find(ismember([Expts{k}.Trials.id],[Expts{j}.Trials.id]));
           bres = PlotExpt(Expts{k},'noplot','condense','Trials',bid);
           nx = 0;
           if isfield(bres.Data.Trials,'se')
               [secounts, ses] = Counts([bres.Data.Trials.se]);
               nrpt = sum(secounts > 1);
               if nrpt > 10
                   serpt = 1;
               end
           else
               serpt = 0;
           end
           for c = 1:length(bres.counts(:));
               if length(bres.counts{c}) > nmin
                   nx = nx+1;
                   bz{nx} = (bres.counts{c} - mean(bres.counts{c}))/std(bres.counts{c});
                   for d = 1:length(bres.ids{c})
                       bspks{nx}{d} = Expts{k}.Trials(bres.tidx{c}(d)).Spikes;
                   end
                   eb(nx,2) = ares.y(nx);
               end
           end
           if nx
               xc = zeros(size(times));
               sxc = zeros(size(times));
               [d, zbin] = min(abs(times));
               for c = 1:length(aspks)
                   sid = randperm(length(aspks{c}));
               for a = 1:length(aspks{c})
                   for s = 1:length(aspks{c}{a})
                   cc = hist(aspks{c}{a}(s) - bspks{c}{a},times);
                   scounts(1,a) = cc(zbin);
                   xc = xc + cc;
                   cc = hist(aspks{c}{a}(s) - bspks{c}{sid(a)},times);
                   sxc = sxc + cc;
                   scounts(2,a) = cc(zbin);
                   end
               end
               end
               b.xc = xc(2:end-1);
               b.sxc=sxc(2:end-1);
               a = cat(2,az{:});
               x = cat(2,bz{:});
               id = (~isnan(a) & ~isnan(x));
               xc = corrcoef(a(id),x(id));
               if isfield(Expts{j},'probes')
               b.probesep = mean((Expts{j}.probes(aid)-Expts{k}.probes(bid)));
               b.probes = [mean(Expts{j}.probes(aid)) mean(Expts{k}.probes(bid))];
               else
               b.probesep = mean((Expts{j}.Header.probe-Expts{k}.Header.probe));
               b.probes = [mean(Expts{j}.Header.probe) mean(Expts{k}.Header.probe)];
               end
               ps = ReadPenSep(Expts{j}.Header);
               if isfield(Expts{j}.Header,'probesep') & Expts{j}.Header.probesep > 0
                    if isnan(ps) || Expts{j}.Header.probesep == ps  %% false if ps is NaN. SO use ps if its not NaN and is different
                        b.probesep = b.probesep * Expts{j}.Header.probesep;
                    else
                        fprintf('%s Distance mismatch: %d vs %d\n',list{j},ps,Expts{j}.Header.probesep);
                        b.probesep = b.probesep * ps;
                    end
               elseif ps > 0
                    b.probesep = b.probesep * ps;
               else
                   [dir, name] = fileparts(list{j});
    id = strmatch(name(1:7),{'lemM010' 'lemM011'  'lemM019' 'lemM020' 'lemM023' 'lemM027' 'lemM029' 'lemM040' 'lemM043' 'lemM051' 'lemM055'...
        'lemM062' 'lemM067' 'lemM072' 'lemM098' 'lemM106' 'lemM107' 'lemM110'});
    tmpps = [100 100 50 50 50 50 50 150 150 150 150 150 150 150 75 75 75 75];
                   if id
                       ps = tmpps(id(1));
                       b.probesep = b.probesep * ps;
                   else
                       ps = NaN;
                   end
                   fprintf('%s has no probesep set to %d\n',list{j},ps);
               end
               b.ntrials = length(id);
               b.rsc(1) = xc(1,2);
               ebvals = unique(eb(~isnan(eb)));
               b.counts = [a; x];
               b.cells = [cells(j) cells(k)];
               b.state.means = cat(2,ares.means(:),bres.means(:));
               b.state.ebvals = ebvals;
               for e = 1:length(ebvals)
                   id = find(eb(:,1) == ebvals(e));
                   a = cat(2,az{id});
                   x = cat(2,bz{id});
                   id = (~isnan(a) & ~isnan(x));
                   xc = corrcoef(a(id),x(id));
                   b.ysc(e) = xc(1,2);
                   b.ym(e,:) = [mean(a) mean(x)];
               end
               xc = corr(b.state.means);
               b.state.rsignal = xc(1,2);
               b.state.dprime = [GetExptQuality(Expts{j}) GetExptQuality(Expts{k})];
               nc = nc+1;
               res(nc) = b;
           end
           end
       end
   end
end
if nc == 0
    res = [];
end


function dprime = GetExptQuality(Expt)

if ~isfield(Expt,'probes')
    if isfield(Expt.Header,'cellnumber')
        mahal = CellToMat(Expt.Header.Clusters,'mahal');
        dprime = mean(max(mahal));
    elseif isfield(Expt.Header,'dips') && length(Expt.Header.dips{1}) > 2
        for j = 1:length(Expt.Header.dips)
            if length(Expt.Header.dips{j}) > 2
            dprimes(j) = abs(Expt.Header.dips{j}(3));
            end
        end
        dprime = mean(dprimes(dprimes > 0));
    else
        dprime = NaN;
    end
    return;
end
ps = Expt.probes;
blks = [Expt.Header.BlockStart Expt.Trials(end).Trial];
for j = 1:length(blks)-1;
    Trials = blks(j):blks(j+1);
    id = find(ismember([Expt.Trials.Trial],Trials));
    if length(id)
        p = median(Expt.probes(id));
        C = Expt.Header.Clusters{j}{1,p};
        if isfield(C,'dprime')
            dp(j) = C.dprime;
        else  %Shouldn't happen
            dp(j) = NaN;
        end
        n(j) = length(id);
    else
        dp(j) = 0;
        n(j) = 0;
    end
end
dprime = WeightedSum(dp,n);

if dprime == 0
    details.dp = dp;
    details.n = n;
    if sum(isnan(dp))
        fprintf('%s Missing Clusters\n',Expt.Header.filename);
    else
        fprintf('%s No Trials\n',Expt.Header.filename);
    end
end


function d = ReadPenSep(Header)

d = NaN;             
if isfield(Header,'Peninfo') && isfield(Header.Peninfo,'trode')
    id = strfind(Header.Peninfo.trode,'Contact');
    if length(id)
        x = id(1);
        id = strfind(Header.Peninfo.trode(x:end),' ');
        if strncmp(Header.Peninfo.trode(id(1)+x:end),'CNT',3)
            x = x+id(2);
        else
            x = x+id(1);
        end
        x = sscanf(Header.Peninfo.trode(x:end),'%d');
        d = x;
    end
end
