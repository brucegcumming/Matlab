function [result, details] = CalcCellSize(C, varargin) 
%[result, details] = CalcCellSize(C, varargin) 
%
%Reads Clusters (AllVPcs) and calculates size, and vertical spread
%If C is a directory name, will read all ClusterTimes files
%If C is a cell array of strings, will read each file/directory in list
%
%CalcCellSize(C, 'nofit') read data but does not fit Gaussians = much quicker
%
%CalcCellSize(result,'plot',plotttype) plots results
%                 default is to plot mahal distance(1D) vs SD of spread
%                           '2dgauss' plots fitted dprime gauss vs SD
%                            'sdspkw' plots spike width vs SD.
% result{}.dip is [1-Dmahal fitdprime 2-Dmahal];

%get voltage amplitude too
%make plotting routine so that can interact.

CellList = [];
CellDetails = [];
result = {};
details = [];
prefix = [];
onepercell = 0;
savecells = 0; 
dofit = 1;
useauto = 0;
ClusterDetails = [];
loadxy = 1;
parallel = 0;

plottype = 'mahal1';
figlabel = 'CellShapes';
j = 1;
while j <= length(varargin)
    if IsClusterDetails(varargin{j})
        ClusterDetails = varargin{j};
    elseif strncmpi(varargin{j},'celllist',7)
        j = j+1;
        CellList = varargin{j};
        j = j+1;
        CellDetails = varargin{j};
    elseif strncmpi(varargin{j},'cellsonly',6)
    elseif strncmpi(varargin{j},'loadxy',6)
        loadxy = 1;
    elseif strncmpi(varargin{j},'onepercell',6)
        onepercell = 1;
    elseif strncmpi(varargin{j},'plot',4)
        j = j+1;
        plottype = varargin{j};
    elseif strncmpi(varargin{j},'nofit',5)
        dofit = 0;
    elseif strncmpi(varargin{j},'parallel',5)
        parallel = 1;
    elseif strncmpi(varargin{j},'prefix',5)
        j = j+1;
        prefix = varargin{j};
    elseif strncmpi(varargin{j},'save',4)
        savecells = 1;
    end
    j = j+1;
end    


if iscellstr(C)
    details.dir = C;
if parallel
    parfor j = 1:length(C)
        results{j} = CalcCellSize(C{j},varargin{:});
    end
else
    for (j = 1:length(C))
        results{j} = CalcCellSize(C{j},varargin{:});
    end
end
    for j = 1:length(results)
        if ~isempty(results{j})
        details.nres(1,j) = length(result);
        result = {result{:} results{j}{:}};
        details.nres(2,j) = length(result);
        end
    end
    return;
end

if ischar(C) & isdir(C)
    details.dir = C;
    details.starttime = now;
    cellfile = [C '/CellList.mat']; 
    if exist(cellfile)
        load(cellfile);
    end
    d = dir([C '/*ClusterTimes.mat']);
    nc = 0;
    for j = 1:length(d)
        if useauto || isempty(regexp(d(j).name,'Auto'))
            cname = [C '/' d(j).name];
            load(cname);
            if loadxy
                ClusterDetails = LoadClusterDetails(cname);
            end
        nc = nc+1;
        res = CalcCellSize(Clusters,ClusterDetails,'CellList',CellList, CellDetails, varargin{:});
        if onepercell
            res = OncePerCell(res);
        end
        result = {result{:} res{:}};
        details.nres(1,nc) = length(res);
        details.nres(2,nc) = length(result);
        end
    end
    if onepercell
        result = OncePerCell(result);
    end
    if savecells
        save([details.dir '/cellshapes.mat'],'result');
    end
    details.endtime = now;
    return;
end

if iscell(C) && isfield(C{1},'amp')
    if onepercell
        C = OncePerCell(C,varargin{:});
    end
    PlotCellSizeResult(C, varargin{:});
    result = C;
    return;
end

if iscell(C)
    x = 0;
    nres = 1;
    [C, errs] = FixCluster(C);
    if ~isempty(errs)
        result{nres}.errs = errs;
    end
    for j = 1:length(C)
        res = CalcCellSize(C{j},'CellList',CellList, CellDetails, varargin{:});
        if length(ClusterDetails) >= j && isfield(ClusterDetails{j},'xy')
            res.xy = ClusterDetails{j}.xy;
            res.clst = ClusterDetails{j}.clst;
        elseif loadxy
             fprintf('Error loading XY %s E%dP%dcl1\n',C{j}.spkfile,C{j}.exptno,C{j}.probe(1));
        end
        if length(res) == 1
            result{nres} = res;
        else
            result(nres:nres+length(res)-1) = res(:);
        end
        nres = nres+length(res);
        if isfield(C{j},'next')
        for k = 1:length(C{j}.next)
            if isfield(C{j}.next{k},'MeanSpike')
                result{nres} = CalcCellSize(C{j}.next{k},'CellList',CellList, varargin{:});
                if length(ClusterDetails) >= j && isfield(ClusterDetails{j},'next') && length(ClusterDetails{j}.next) >= k
                    result{nres}.xy = ClusterDetails{j}.next{k}.xy;
                    result{nres}.clst = ClusterDetails{j}.clst;
                elseif loadxy
                    fprintf('Error loading XY %s E%dP%dcl%d\n',C{j}.spkfile,C{j}.exptno,C{j}.probe(1),k+1);
                end
                nres = nres+1;
            end
        end
        end
    end
    if onepercell
        result = OncePerCell(result);
    end
    return;
    end

 if isfield(C,'bytes') % a directory result
     for j = 1:length(C);
         dirs{j} = [prefix '/' C(j).name];
     end
     [result, details] = CalcCellSize(dirs,varargin{:});
 end
    

if isfield(C, 'MeanSpike')
    [C, errs] = FixCluster(C);
    nerr = 0;
    if ~isfield(C,'cluster')
        nerr = nerr+1;
       result.errs{nerr} = 'clustermissing';
       C.cluster = 1;
    end
    sds = std(C.MeanSpike.ms,[],2);
    result.vmax = max(abs(C.MeanSpike.ms(:)));
    result.triggerpt = find(C.spts ==0);
    result.cluster = C.cluster;
    [fit, maxi] = FitSDs(sds, dofit);
    result.V = C.MeanSpike.ms(maxi,:);
    if maxi == 1
        result.nextV = C.MeanSpike.ms([maxi+1 maxi+2],:);
    elseif maxi == size(C.MeanSpike.ms,1)
        result.nextV = C.MeanSpike.ms([maxi-1 maxi-2],:);
    else
        result.nextV = C.MeanSpike.ms([maxi-1 maxi+1],:);
    end
    if C.probe == 1
        chspk = [1 2 3];
    elseif C.probe == size(C.MeanSpike.ms,1)
        chspk = C.probe + [-2 -1 0];
    else
        chspk = C.probe + [-1 0 1];
    end
    result.chspk = chspk;
    vdiff = C.MeanSpike.ms-C.MeanSpike.mu;
    if isfield(C.MeanSpike,'vdprime')
        if size(C.MeanSpike.vdprime,1) >= max(chspk)
            result.vdprime = C.MeanSpike.vdprime(chspk,:);
        else
            result.vdprime = C.MeanSpike.vdprime;
        end
    end
    result.vdiff = vdiff(result.chspk,:);
    result.amp = sds;
    result.sd = abs(fit.sd);
    if isfield(fit,'guess')
        resuld.isd = abs(fit.guess(2));
    end
    usds = std(C.MeanSpike.mu,[],2);
    [ufit, maxi] = FitSDs(usds,dofit);
    result.usd = abs(ufit.sd);

    if isfield(C,'mahal')
        result.dip(1) = C.mahal(4); %1d
        result.dip(3) = C.mahal(1); %2d
    end
    if isfield(C,'fitdprime')
        result.dip(2) = C.fitdprime(1);
    else
        result.dip(2) = NaN;
    end
    result.probe = C.probe;
    result.eid = C.exptno;
    result.cell = isacell(C,CellList, CellDetails);
end


function good = IsClusterDetails(C)
good = 0;
if iscell(C) 
 for j = 1:length(C)
     if isfield(C{j},'xy') && isfield(C{j},'Evec')
         good = good+1;
     end
 end
end

function [fit, b] = FitSDs(sds, dofit)
    [a,b] = max(sds);
    trueb = b;
    if b ==1 
        sds(3:end) = sds(1:end-2);
        sds(2) = sds(4);
        sds(1) = sds(5);
        [a,b] = max(sds);
    elseif b == length(sds)
        sds(1:end-2) = sds(3:end);
        sds(end-1) = sds(end-3);
        sds(end) = sds(end-4);
        [a,b] = max(sds);
    end
    x = [1:length(sds)]';
    guess(1) = b;
    guess(2) = std(((x-mean(x)) .* sds))/mean(sds);
    guess(2) = 1;
    guess(3) = a;
    guess(4) = prctile(sds,40);
    c = sort(sds-guess(4),'descend');
    sdr = mean(c(4:5))./c(1);
    guesssd(1) = sqrt(-8/log(sdr));
    sdr = mean(c(2:3))./c(1);
    guesssd(2) = sqrt(-2/log(sdr));
    guess(2) = min(guesssd);
    if dofit
        fit = FitGauss(1:length(sds),sds','freebase','guess',guess);
        if trueb == length(sds)
            fit.fitted(3:end) = fit.fitted(1:end-2);
            fit.mean = fit.mean+2;
        elseif trueb == 1
            fit.fitted(1:end-2) = fit.fitted(3:end);
            fit.mean = fit.mean-2;
        end
    else
        fit.sd = guess(2);
        fit.amp = guess(3);
        fit.params = guess;
    end


function cell = isacell(C, CellList, CellDetails)

eid =[];
if isfield(CellDetails,'exptids')
    eid = find(CellDetails.exptids == C.exptno);
end
if isempty(eid)
    eid = floor(C.exptno);
end
if isempty(CellList)  || eid < 1 || eid > size(CellList,1) || C.cluster > size(CellList,3)
    cell = NaN;
else
    cell = CellList(eid, C.probe, C.cluster);
end

function PlotCellSizeResult(R, varargin)

figlabel = 'CellShapes';
plottype = 'mahal1';
cellsonly = 0;
details = [];

callback = @HitScatter;
selectcrit = 0;
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'dir')
        details = varargin{j};
    elseif strncmpi(varargin{j},'callback',8)
        j = j+1;
        callback = varargin{j};
    elseif strncmpi(varargin{j},'cellsonly',6)
        cellsonly = 1;
    elseif strncmpi(varargin{j},'plot',4)
        j = j+1;
        plottype = varargin{j};
    elseif strncmpi(varargin{j},'select',4)
        j = j+1;
        selectcrit = varargin{j};
    end
    j = j+1;
end
sds = CellToMat(R,'sd');
usds = CellToMat(R,'usd');
dips = CellToMat(R,'dip');
cells = CellToMat(R,'cell');
cls = CellToMat(R,'cluster');
cid = find(cells > 0);
nid = find(~(cells > 0)); %include NaNs
c2id = find(cls > 1);
if cellsonly
    nid = [];
end
[a, isnew] = GetFigure(figlabel);
DATA.plot.showxy = 1;
DATA.plot.density = 0;
DATA.plot.refit = 0;
DATA.plot.nrefit = 0;
DATA.plot.plottype = plottype;
DATA.cells = R;
DATA.details = details;
DATA.prefix = '';

if isnew
    hm = uimenu(a,'Tag','OptionMenu','label','Options');
    sm = uimenu(hm,'label','Refit','tag','refit','callback',{@OptionMenu, 'setbytag'});
    sm = uimenu(hm,'label','DensityPlot','tag','density','callback',{@OptionMenu, 'setbytag'});
    sm = uimenu(hm,'label','Show XY','tag','showxy','callback',{@OptionMenu, 'setbytag'});
    sm = uimenu(hm,'label','Refit several mixtures','tag','nrefit','callback',{@OptionMenu, 'setbytag'});
    sm = uimenu(hm,'label','Refit 3','tag','refit3','callback',{@OptionMenu, 'refit3'});
    sm = uimenu(hm,'label','Refit 4','tag','refit4','callback',{@OptionMenu, 'refit4'});
    sm = uimenu(hm,'label','Refit 5','tag','refit5','callback',{@OptionMenu, 'refit5'});
    sm = uimenu(hm,'label','Refit All','tag','refit5','callback',{@OptionMenu, 'refitall'});
    SetMenuChecks(hm,DATA.plot);
end
hold off;
if strcmp(plottype,'2dgauss')
    myscatter(dips(nid,2),sds(nid),'o','ids',nid,'buttonpress',callback);
    myscatter(dips(cid,2),sds(cid),'ro','ids',cid,'buttonpress',callback);
elseif strcmp(plottype,'mahal') %dips 2 vs 1 = my fit vs 1D GM
    myscatter(dips(nid,2),dips(nid,1),'o','ids',nid,'buttonpress',callback);
    myscatter(dips(cid,2),dips(cid,1),'ro','ids',cid,'buttonpress',callback);
    myscatter(dips(c2id,2),dips(c2id,1),'go','ids',c2id,'buttonpress',callback);
elseif strcmp(plottype,'mahal2')%dips 2 vs 3 = my fit vs 2D GM
    myscatter(dips(nid,2),dips(nid,3),'o','ids',nid,'buttonpress',callback,'color',[0.5 0.5 0.5]);
    colors = mycolors('spkcolors');
    for j = 1:length(cid)
        cc{j} = colors{cls(cid(j))+1};
    end
    myscatter(dips(cid,2),dips(cid,3),'ro','ids',cid,'colors', cc,'buttonpress',callback);
 %   myscatter(dips(c2id,2),dips(c2id,3),'go','ids',c2id,'buttonpress',callback);
    xlabel('dprime from indep 1-D fits');
    ylabel('dprime 2D GM fit');
elseif strcmp(plottype,'mahal3')%dips 1 vs 3 = 1DGM vs 2D GM
    myscatter(dips(nid,1),dips(nid,3),'o','ids',nid,'buttonpress',callback);
    myscatter(dips(cid,1),dips(cid,3),'ro','ids',cid,'buttonpress',callback);
    myscatter(dips(c2id,1),dips(c2id,3),'go','ids',c2id,'buttonpress',callback);
    xlabel('dprime 1D GM fit');
    ylabel('dprime 2D GM fit');
elseif strcmp(plottype,'muamp')
    myscatter(sds(nid),usds(nid),'o','ids',nid,'buttonpress',callback);
    myscatter(sds(cid),usds(cid),'ro','ids',cid,'buttonpress',callback);
    xlabel('Sigma for SU');
    xlabel('Sigma for MU');
elseif strcmp(plottype,'mahalmind')
    [a,b] = max(CellToMat(R,'mind')');
    myscatter(dips(nid,2),a(nid),'o','ids',nid,'buttonpress',callback);
    myscatter(dips(cid,2),a(cid),'ro','ids',cid,'buttonpress',callback);
    myscatter(dips(c2id,1),a(c2id),'go','ids',c2id,'buttonpress',callback);
    xlabel('2-D GM');
    xlabel('max(mind)');
elseif strcmp(plottype,'mind')
    colors = mycolors;;
    hold off;
    for j = 1:length(R) 
         if isfield(R{j},'mind')
             ci = 1+mod(j-1,length(colors));
             plot(R{j}.mind, R{j}.maxd,'o-','color',colors{ci},'buttondownfcn',{@HitScatter, j,0});
             hold on;
         else
             fprintf('Missing mind %s (%d)\n',IDstr(DATA,j),j);
         end
     end;
elseif strcmp(plottype,'sdhist')
    hist(sds(dips(:,2) <-3),100);
elseif strcmp(plottype,'dprimeb') %dips 2 vs 1 = my fit vs 1D GM
    PlotDprimes(DATA.cells,cid,'oprobe');
elseif strcmp(plottype,'dprime') %dips 2 vs 1 = my fit vs 1D GM
    PlotDprimes(DATA.cells,cid);
elseif strcmp(plottype,'pcs') %dips 2 vs 1 = my fit vs 1D GM
    PlotPCs(DATA.cells,cid);
elseif strcmp(plottype,'shapes') %dips 2 vs 1 = my fit vs 1D GM
    PlotShapes(DATA.cells,cid);
elseif sum(strcmp(plottype,{'sdspkw' 'spkw'}))
    for j = 1:length(R)
        [a,b] = max(R{j}.V);
        [c,d] = min(R{j}.V);
        a = max([a -c]);
        id = find(abs(R{j}.V) > a/5);
        if isempty(id)
            w(j,:) = NaN;
        else
        w(j,1) = abs(b-d);
        w(j,2) = id(end)-id(1);
        end
    end
    if selectcrit == 0
        id = 1:size(dips,1);
    else
        id = find(dips(:,2) < -2.5);
    end
    if strcmp(plottype,'spkw')
        myscatter(w(id,1),w(id,2),'o','buttonpress',callback);
    else
        myscatter(w(id,2),sds(id),'o','buttonpress',callback);
        id = find(dips(:,2) < -2.5);
    end
else
    myscatter(dips(nid,1),sds(nid),'o','ids',nid,'buttonpress',callback);
    myscatter(dips(cid,1),sds(cid),'ro','ids',cid,'buttonpress',callback);
end
DATA.toplevel = gcf;
set(gcf,'UserData',DATA);


function OptionMenu(a,b, fcn)
DATA = GetDataFromFig(a);
onoff = {'off' 'on'};
if strcmp(fcn,'setbytag')
    f = get(a,'Tag');
    if ~isfield(DATA.plot,f)
        DATA.plot.(f) = 1;
    else
        DATA.plot.(f) = ~DATA.plot.(f);
    end
    set(a,'checked',onoff{DATA.plot.(f)+1});
elseif sum(strcmp(fcn,{'refit3' 'refit4' 'refit5'}))
    DATA.plot.refit = 2+find(strcmp(fcn,{'refit3' 'refit4' 'refit5'}));
elseif strcmp(fcn,'refitall')
    DATA = RefitAll(DATA);
end
set(DATA.toplevel,'UserData',DATA);
HitScatter(DATA,  [], DATA.selected);

function SetMenuChecks(hm, S)
sms = findobj(hm, 'type', 'uimenu');
onoff = {'off' 'on'};
for j = 1:length(sms)
    t = get(sms(j),'tag');
    if isfield(S,t)
        set(sms(j),'checked',onoff{1+S.(t)});
    end
end

function PlotDprimes(R, cid, varargin)

otherprobe = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'oprobe',5)
        otherprobe = 1;
    end
    j = j+1;
end

for k = 1:length(cid)
    C = R{cid(k)};
    p = find(C.chspk == C.probe);
    vsd = C.vdiff./C.vdprime;
    vdiff = C.vdiff - repmat(C.vdiff(:,C.triggerpt),1, size(C.vdiff,2));
    vdprime = vdiff./vsd;
    s = abs(smooth([vdprime(p,1)/2 vdprime(p,:) vdprime(p,end)/2],1,'gauss'));
    s = s(2:end-1);
    id = 1+LocalMaxima(abs(s));
    if isempty(id)
        [a,id] = max(s);
    end
    [a,b] = sort(s(id),'descend');
    dpi(k,1:length(id)) = id(b)-C.triggerpt;
    dp(k,1:length(id)) = a;
    if size(R{cid(k)}.vdprime,1) > p       
%        vdiff = R{cid(k)}.MeanSpike.ms-R{cid(k)}.MeanSpike.mu;
        s = abs(smooth(R{cid(k)}.vdprime(p+1,:),2,'gauss'));
        id = 1+LocalMaxima(abs(s));
        [a,b] = sort(s(id),'descend');
        adpi(k,1:length(id)) = id(b)-C.triggerpt;
        adp(k,1:length(id)) = a;
    end
    if p > 1
        s = abs(smooth(R{cid(k)}.vdprime(p-1,:),2,'gauss'));
        id = 1+LocalMaxima(abs(s));
        [a,b] = sort(s(id),'descend');
        bdpi(k,1:length(id)) = id(b)-C.triggerpt;
        bdp(k,1:length(id)) = a;
    end
end
if otherprobe == 0
    myscatter(dpi(:,1),dp(:,1),'o','ids',cid,'buttonpress',@HitScatter);
    gid = find(dpi(:,2) ~= 0);
    myscatter(dpi(gid,2),dp(gid,2),'ro','ids',cid,'buttonpress',@HitScatter);
else
    dpi = cat(1, adpi,bdpi);
    dp = cat(1, adp,bdp);
    id = find(dp(:,1) > 0);
    dp = dp(id,:);
    dpi = dpi(id,:);
    myscatter(adpi(:,1),dp(:,1),'o','ids',cid,'buttonpress',@HitScatter);
    myscatter(bdpi(:,2),dp(:,2),'ro','ids',cid,'buttonpress',@HitScatter);
    myscatter(bdpi(:,1),dp(:,1),'o','ids',cid,'buttonpress',@HitScatter);
    myscatter(bdpi(:,2),dp(:,2),'ro','ids',cid,'buttonpress',@HitScatter);
    gid = find(dpi(:,2) ~= 0);
end
[y,x] = smhist(dpi(:,1));
yl = get(gca,'ylim');
plot(x,y .* yl(2)./max(y));
[y,x] = smhist(dpi(gid,2));
yl = get(gca,'ylim');
plot(x,y .* yl(2)./max(y),'r-');
x = unique(dpi(:));
for j = 1:length(x);
    id = find(dpi(:,1)==x(j));
    sums(j,1) = sum(dp(id,1));
    sums(j,2) = sum(dp(id,2));
end
plot(x,sums(:,1).*yl(2)./max(sums(:,1)),'r');
hold on;
plot(x,sums(:,2).*yl(2)./max(sums(:,2)),'b');
%GetFigure('CellHistogram');

function id = LocalMaxima(x)

sgn = diff(sign(diff(x)));
id = find(sgn < 0);

function PlotPCs(R, id)
for j = 1:len(id)
    V(j,:) = R{id(j)}.V;
end
[a,b] =eig(cov(V));
pcs = V * a;
PlotND(pcs(:,37:40),[],'marker','o','callback',@HitScatter, id);

function PlotShapes(R, id)

for j = 1:length(id)
    V = R{id(j)}.V;
    V = V./std(V);
    [a,b] = min(V);
    minv(j) = a;
    [c,d] = max(V(1:b));
    premaxpt(j) = d-b;
    premax(j) = c;
    [c,d] = max(V(b:end));
    postmaxpt(j) = d;
    postmax(j) = c;
    endpts = [b+1:length(V)];
    prepts = [1:b-1];
    postshape(j) = V(endpts)*endpts'./sum(endpts);
    preshape(j) = V(prepts)*prepts'./sum(prepts);
end
nr=2;
nc=3;
subplot(nr,nc,1);
hold off;
myscatter(premaxpt,premax,'o','ids',id,'buttonpress',@HitScatter);
xlabel('Pre max pt');
ylabel('Pre max');

subplot(nr,nc,2);
hold off;
myscatter(postmaxpt,postmax,'o','ids',id,'buttonpress',@HitScatter);
xlabel('Post max pt');
ylabel('Post max');

subplot(nr,nc,3);
hold off;
myscatter(minv,postmax,'o','ids',id,'buttonpress',@HitScatter);
xlabel('Min V');
ylabel('Post max');


subplot(nr,nc,4);
hold off;
myscatter(minv,premax,'o','ids',id,'buttonpress',@HitScatter);
xlabel('Min V');
ylabel('pre max');

subplot(nr,nc,5);
hold off;
myscatter(postmax,premax,'o','ids',id,'buttonpress',@HitScatter);
xlabel('Post max');
ylabel('Pres max');

subplot(nr,nc,6);
hold off;
myscatter(preshape,postshape,'o','ids',id,'buttonpress',@HitScatter);
xlabel('PreSum');
ylabel('PostSum');

function HitScatter(a,b,id, idb)
DATA = GetDataFromFig(a);
C = DATA.cells{id};
fit = [];
if isfield(C,'tag')
fprintf('%s:',C.tag)
elseif isfield(DATA.details,'dir') && ischar(DATA.details.dir) %just one directory
    aid = 0;
    cname = DATA.details.dir; 
elseif isfield(DATA.details,'nres')
    aid = find(DATA.details.nres(1,:) <= id & id <= DATA.details.nres(2,:));
    fprintf('%s:',DATA.details.dir{aid});
    cname = DATA.details.dir{aid(1)}; 
else
    aid = 0;
end
if strcmp(DATA.plot.plottype,'mahal1') && DATA.plot.refit
    [fit, maxi] = FitSDs(DATA.cells{id}.amp,1);
    DATA.cells{id}.sd = abs(fit.sd);
    set(DATA.toplevel,'UserData',DATA);
end
DATA.selected = id;

fprintf('%d cell %d E%dP%dcl%d dip %s\n',id,C.cell,C.eid,C.probe,C.cluster,sprintf('%.1f ',C.dip));
GetFigure('Spike');
subplot(2,1,1);
hold off;
plot(C.amp);
if isfield(fit,'fitted')
    hold on;
    plot(fit.fitted,'r');
end
title(sprintf('E%dP%d',C.eid,C.probe));
subplot(2,1,2);
plot(C.V);

title(sprintf('E%dP%dcl%d',C.eid,C.probe,C.cluster));


if DATA.plot.showxy && length(aid) == 1
    p = C.probe;
    if length(DATA.prefix) > 1  && aid > 0
        cname = strrep(DATA.details.dir{aid},'/data/',DATA.prefix);
    end
    
    if isfield(C,'xy')
        X = C;
    else
        [Clusters,a,b] = LoadCluster(cname, C.eid,'getxy');
        X = Clusters{p};
        fprintf('Load  took %.2f\n',b.loadtime);
    end
    if isfield(X,'xy')
        colors = mycolors('spkcolors');
        GetFigure('XYplot');
        if aid > 0
            xid = DATA.details.nres(1,aid):DATA.details.nres(2,aid);
            xid = xid(xid <= length(DATA.cells)); %in case removed some
        else
            xid = 1:length(DATA.cells);
        end
        expts = CellToMat(DATA.cells(xid),'eid');
        probes = CellToMat(DATA.cells(xid),'probe');
        oid = find(expts == X.eid & probes == X.probe);
        if length(oid) > 1 %may be more spaces
            xy = [];
            for j = 1:length(oid)
                xy = cat(2,xy,DATA.cells{xid(oid(j))}.xy(:,1));
                xy = cat(2,xy,DATA.cells{xid(oid(j))}.xy(:,2));
            end    
            c = corrcoef(xy);
            need = ones(1,size(c,1));
            for j = 2:size(c,1)
                for k = 1:j-1
                    if c(j,k) > 0.98
                        need(k) = 0;
                    end
                end
            end
            xy = xy(:,find(need));
        else
            xy = X.xy;
        end
        if DATA.plot.density
        PlotND(xy,[],'idlist',X.clst,'density');
        else
        PlotND(xy,[],'idlist',X.clst,'colors',colors);
        end
        
        if DATA.plot.nrefit 
            [X,G] = Refit(DATA, X);
            GetFigure('FitDistance');
            newid = cluster(G, X.xy);
            plot(X.mind,X.maxd,'ro-');
            GetFigure('XYplot');
            PlotND(X.xy,[],'idlist',newid,'colors',colors);
        end
        
        if DATA.plot.refit
            G = GMfit(X.xy,2,1,'idlist',X.clst);
            [X,Gn] = Refit(DATA, X);
            nc = size(Gn.mu,1);
            newid = cluster(Gn, X.xy);
            PlotND(X.xy,[],'idlist',newid,'colors',colors);
            title(sprintf('New Fit %.2f(2) %.2f(%d)',GMdprime(G),GMdprime(Gn),nc));
        end
    end
end
set(DATA.toplevel,'UserData',DATA);


function str = IDstr(DATA,id)

C = DATA.cells{id};
cname = [];
if isempty(C)
    str= 'Empty';
    return;
end
if isfield(DATA.details,'dir') && ischar(DATA.details.dir) %just one directory
    aid = 0;
    cname = DATA.details.dir; 
elseif isfield(DATA.details,'nres')
    aid = find(DATA.details.nres(1,:) <= id & id <= DATA.details.nres(2,:));
    cname = DATA.details.dir{aid}
end
str = sprintf('%s E%dP%dcl%d (cell%d)',cname,C.eid,C.probe,C.cluster,C.cell);

function DATA = RefitAll(DATA)
for j = 1:length(DATA.cells)
    if isfield(DATA.cells{j},'xy')
    DATA.cells{j} = Refit(DATA, DATA.cells{j});
    else
        fprintf('%d No XY\n',j);
    end
end



function [X, Gn] = Refit(DATA, X)

    
    
    if DATA.plot.nrefit
        for nd = 1:4
            G{nd} = GMfit(X.xy,nd+1,1);
            [a,b] = gmdprime(G{nd});
            mind(nd) = min(b.d(:));
            maxd(nd) = max(b.d(:));
        end
        id = find(mind > 2);
        if ~isempty(id)
            best = id(end);
        else
            best = 1;
        end
        [a,b] = max(mind);
        if b > best
            best = b;
        end
        X.mind = mind;
        X.maxd = maxd;
        Gn = G{best};
    end

    if DATA.plot.refit
        if DATA.plot.refit > 1
            nc = DATA.plot.refit;
        else
            nc = length(unique(X.clst));
        end
        [Gn, allg] = GMfit(X.xy,nc,1,'idlist',X.clst);
    end

function res = OncePerCell(res, varargin)
details = [];
tag = [];
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'dir')
        newres = {};
        details = varargin{j};
        for k = 1:size(details.nres,2)
            a = details.nres(1,k);
            b = details.nres(2,k);
            [c, dirname] = fileparts(details.dir{k});
            if b > a
                c = OncePerCell(res(1+a:b),'tag',dirname);
                if isempty(c)
                    fprintf('No cells in %s\n',details.dir{k});
                else
                    newres = {newres{:} c{:}};
                end
            end
        end
        res = newres;
        return;
    elseif strcmpi(varargin{j},'tag')
        j = j+1;
        tag = varargin{j};
    end
    j = j+1;
end

cellids = CellToMat(res,'cell');
dips = CellToMat(res,'dip');
cells = unique(cellids);
cells = cells(cells > 0);
if isempty(cells)
    res = {};
else
for j = 1:length(cells)
    id = find(cellids == cells(j));
    [a,b] = min(dips(id,1));
    goodid(j) = id(b);
    if ~isempty(tag)
        res{goodid(j)}.tag = tag;
    end
end
res = res(goodid);
end