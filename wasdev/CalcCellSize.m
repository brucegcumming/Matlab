function [result, details] = CalcCellSize(C, varargin) 
%Reads Clusters (AllVPcs) and calculates size, and vertical spread
%[result, details] = CalcCellSize(C, varargin) 
%
%If C is a directory name, will read all ClusterTimes files
%If C is a cell array of strings, will read each file/directory in list
%
%CalcCellSize(C, 'nofit') read data but does not fit Gaussians = much quicker
%
%CalcCellSize(result,'onepercell') returns a result containing only marked cells, and only 1 file per cell
%
%CalcCellSize(result,'loadmean') reads in meanspike data for each element
%
%
%CalcCellSize(result,'plot',plotttype) plots results
%                 default is to plot mahal distance(1D) vs SD of spread
%                    '2dgauss' plots fitted dprime gauss vs SD
%                    'sdspkw' plots spike width vs SD.
%                    'mahal' fitted dprime vs 1D mahal distances
%                    'mahal2' fitted dprime vs 2D mahal distances
%                    'mahal3' 2D vs 2D mahal
%                    'muamp' spread SD for mu vs sd for SU
%                    'dprime' histogram of usefulness of differen voltge samples, using dprime
%                    'dprimeb' histogram of usefulness voltage samples from adjacent probes
%                    'dprimec' scatterplot most useful vs 2nd most useful
%                    'pcs' PCA analysis of all meanpike V.  (Need to limit this to same trigger point/sign)
%                    'shapes' plots Pre/Post max V and sample # shape metrics
%                    'shapeim' image plot of all shapes
%                    'spkw' compares 2 width measures
%                    'sdspkw' width vs spread
%
%   Additional args
%           'dotprod' uses dot-product amplitude estimates rather than SD, if available
%
%            'csdonly', 'withcsd', 'withdy', 'imageonly' are passed to PlotMeanSpike
%
%
%  To build a set from scratch
%        [shapes, details] = CalcCellSize('/b/data','build','lem')  seaches lem dirs, builds for each
%        cellshapes = CalcCellSize(shapes, details, 'onepercell','loadmean'); 
%                   selects only identified cells, and one shape for each.
%                   Then loads Full MeanSpike  and xy plot for each into the struct.
%        cellshpaes = CalcCellSize(cellshapes, details, 'loadmean'); 
%
%  If a refit has been run to define multiple clusters, then max and min of mahal
%  across differet N cluster are also recoded
%                    'mind'  plot min(D) vs max(D)
%                    'mahalmind' 2D mahal vs min(D)
%
% result{}.dip is [1-Dmahal fitdprime 2-Dmahal];
%
%  To load existing cellshapes files for many folders:
% d = TreeFind(rootpath,'name','cellshapes.mat');  
% CalcCellSize(d)


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
loadmean = 0;
recalc = 0;
removemean = 0;

plottype = 'mahal1';
figlabel = 'CellShapes';
j = 1;
while j <= length(varargin)
    if IsClusterDetails(varargin{j})
        ClusterDetails = varargin{j};
    elseif strncmpi(varargin{j},'build',5)
        j = j+1;
        [result, details] = BuildMeanData(C, varargin{j:end});
        return;
    elseif strncmpi(varargin{j},'celllist',7)
        j = j+1;
        CellList = varargin{j};
        j = j+1;
        CellDetails = varargin{j};
    elseif strncmpi(varargin{j},'cellsonly',6)
    elseif strncmpi(varargin{j},'loadmean',6)
        loadmean = 1;
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
    elseif strncmpi(varargin{j},'recalc',5)
        recalc = 1;
    elseif strncmpi(varargin{j},'save',4)
        savecells = 1;
    end
    j = j+1;
end    


if iscellstr(C)
    details.dir = C;
    if parallel
        parfor j = 1:length(C)
            fprintf('Worker %s Calculating Size/Shape for %s\n',WorkerString(),C{j});
            results{j} = CalcCellSize(C{j},varargin{:});
            results{j}.dirid = j;
        end
    else
        for (j = 1:length(C))
            fprintf('Calculating Size/Shape for %s\n',C{j});
            results{j} = CalcCellSize(C{j},varargin{:});
            for k = 1:length(results{j});
                results{j}{k}.dirid = j;
            end
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
    savefile = [details.dir '/cellshapes.mat'];
    if exist(savefile) && ~recalc
        ts = now;
        fprintf('Loading %s',savefile);
        load(savefile);
        fprintf(' took %.2f\n',mytoc(ts));
        if onepercell
            result = OncePerCell(result);
        end
        details.endtime = now;
        details.fromdist = 1;
        return;
    end
    details.fromdisk = 0;
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

if ischar(C) && exist(C)  %file
    details.dir = fileparts(C);
    details.starttime = now;
    load(C);
        details.endtime = now;
        return;   
end

if iscell(C) && isfield(C{1},'amp')
    if onepercell
        [C, details]= OncePerCell(C,varargin{:});
    end
    if loadmean
        C = AddMeanSpike(C, varargin{:});
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
    meanamp = mean(C.MeanSpike.ms); %will be zero if mean subtraction not reversed
    if removemean
        for j = 1:size(C.MeanSpike.ms,1)
            C.MeanSpike.ms(j,:) = C.MeanSpike.ms(j,:) - meanamp; 
        end
    end
    result.meanamp = std(meanamp);
    sds = std(C.MeanSpike.ms,[],2);
    result.vmax = max(abs(C.MeanSpike.ms(:)));
    result.triggerpt = find(C.spts ==0);
    result.cluster = C.cluster;
    [fit, maxi] = FitSDs(sds, dofit);
    result.V = C.MeanSpike.ms(maxi,:);
    result.muV = C.MeanSpike.mu(maxi,:);
    if maxi == 1
        result.nextV = C.MeanSpike.ms([maxi+1 maxi+2],:);
    elseif maxi == size(C.MeanSpike.ms,1)
        result.nextV = C.MeanSpike.ms([maxi-1 maxi-2],:);
    else
        result.nextV = C.MeanSpike.ms([maxi-1 maxi+1],:);
    end
    result.dotA =  C.MeanSpike.ms * result.V';
    A = result.dotA./max(result.dotA);
    
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
if trueb == 1
    if b > 2
        b = b-2;
    else
        b = 4-b;
    end
elseif trueb == length(sds)
    if b > length(sds)
        b = trueb - (b-trueb);
    else
        b = b+2;
    end
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
argon = {};
plotargs = {};
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'dir')
        details = varargin{j};
    elseif strncmpi(varargin{j},'callback',8)
        j = j+1;
        callback = varargin{j};
    elseif strncmpi(varargin{j},'cellsonly',6)
        cellsonly = 1;
    elseif sum(strncmpi(varargin{j},{'csdonly' 'smooth' 'withcsd' 'withdy'},6))
        plotargs = {plotargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'plot',4)
        j = j+1;
        plottype = varargin{j};
    elseif strncmpi(varargin{j},'select',4)
        j = j+1;
        selectcrit = varargin{j};
    else 
        argon = {argon{:} varargin{j}};
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
DATA.plotargs = plotargs;

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
    xlabel('Fit Dprime')
    ylabel('Amp SD');
elseif strcmp(plottype,'mahal') %dips 2 vs 1 = my fit vs 1D GM
    myscatter(dips(nid,2),dips(nid,1),'o','ids',nid,'buttonpress',callback);
    myscatter(dips(cid,2),dips(cid,1),'ro','ids',cid,'buttonpress',callback);
    myscatter(dips(c2id,2),dips(c2id,1),'go','ids',c2id,'buttonpress',callback);
    xlabel('Fit Dprime')
    ylabel('Mahel 1-D');
elseif strcmp(plottype,'mahal2')%dips 2 vs 3 = my fit vs 2D GM
    myscatter(dips(nid,2),dips(nid,3),'o','ids',nid,'buttonpress',callback,'color',[0.5 0.5 0.5]);
    colors = mycolors('spkcolors');
    for j = 1:length(cid)
        cc{j} = colors{cls(cid(j))+1};
    end
    myscatter(dips(cid,2),dips(cid,3),'ro','ids',cid,'colors', cc,'buttonpress',callback);
 %   myscatter(dips(c2id,2),dips(c2id,3),'go','ids',c2id,'buttonpress',callback);
    xlabel('dprime from indep 1-D fits');
    ylabel('Mahal 2D');
elseif strcmp(plottype,'mahal3')%dips 1 vs 3 = 1DGM vs 2D GM
    myscatter(dips(nid,1),dips(nid,3),'o','ids',nid,'buttonpress',callback);
    myscatter(dips(cid,1),dips(cid,3),'ro','ids',cid,'buttonpress',callback);
    myscatter(dips(c2id,1),dips(c2id,3),'go','ids',c2id,'buttonpress',callback);
    xlabel('Mahal 1D GM fit');
    ylabel('Mahal 2D GM fit');
elseif strcmp(plottype,'muamp')
    myscatter(sds(nid),usds(nid),'o','ids',nid,'buttonpress',callback);
    myscatter(sds(cid),usds(cid),'ro','ids',cid,'buttonpress',callback);
    xlabel('Sigma for SU');
    ylabel('Sigma for MU');
    set(gca,'xlim',[0 4],'ylim',[0 4]);
    refline(1);
elseif strcmp(plottype,'mahalmind')
    [a,b] = max(CellToMat(R,'mind')');
    if length(a) >= size(dips,1)
    myscatter(dips(nid,2),a(nid),'o','ids',nid,'buttonpress',callback);
    myscatter(dips(cid,2),a(cid),'ro','ids',cid,'buttonpress',callback);
    myscatter(dips(c2id,1),a(c2id),'go','ids',c2id,'buttonpress',callback);
    xlabel('2-D GM');
    xlabel('max(mind)');
    end
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
elseif strcmp(plottype,'distance') %dips 2 vs 1 = my fit vs 1D GM
    PlotDistance(DATA.cells,cid,nid,argon{:});
elseif strcmp(plottype,'dprime') %dips 2 vs 1 = my fit vs 1D GM
    PlotDprimes(DATA.cells,cid);
elseif strcmp(plottype,'dprimec') %dprime best sampel
    PlotDprimes(DATA.cells,cid,'scatter');
elseif strcmp(plottype,'pcs') %dips 2 vs 1 = my fit vs 1D GM
    PlotPCs(DATA.cells,cid);
elseif strncmp(plottype,'csdpc',5) %dips 2 vs 1 = my fit vs 1D GM
    PlotCSDPCs(DATA.cells,cid);
elseif strcmp(plottype,'shapes') %dips 2 vs 1 = my fit vs 1D GM
    PlotShapes(DATA.cells,cid);
elseif strcmp(plottype,'shapeim') %dips 2 vs 1 = my fit vs 1D GM
    PlotShapes(DATA.cells,cid,'image');
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
        xlabel('Width (max-min)');
        ylabel('Width (20%%max');
    else
        myscatter(w(id,2),sds(id),'o','buttonpress',callback);
        xlabel('Width (max-min)');
        ylabel('Spread (sigma)');
        id = find(dips(:,2) < -2.5);
    end
else
    myscatter(dips(nid,1),sds(nid),'o','ids',nid,'buttonpress',callback);
    myscatter(dips(cid,1),sds(cid),'ro','ids',cid,'buttonpress',callback);
    xlabel('Mahal distance');
    ylabel('Size(SD)')
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

function PlotDistance(R, cid, nid, varargin)
scatter = 0;
colors = mycolors;
nbin = 0:0.01:1;
usedotproduct = 0;
signed = 0;
someprobes = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'dotpro',5)
        usedotproduct = 1;
    elseif strncmpi(varargin{j},'oprobe',5)
        otherprobe = 1;
    elseif strncmpi(varargin{j},'someprobes',5)
        j = j+1;
        someprobes = varargin{j};
    elseif strncmpi(varargin{j},'signed',5)
        signed = 1;
    elseif strncmpi(varargin{j},'scatter',5)
        scatter = 1;
    end
    j = j+1;
end


subplot(1,2,1);
hold off;
ds = [];
amps = [];
for k = 1:length(cid)
    C = R{cid(k)};
    if isempty(someprobes) || ismember(C.probe,someprobes)
    if isfield(C,'dotA') && usedotproduct %use dot product if available
        [a,p] = max(C.dotA);
        probes = 1:length(C.dotA);
        A = C.dotA ./a;
    else
        [a,p] = max(C.amp);
        probes = 1:length(C.amp);
        A = C.amp ./a;
    end
    if signed
        d = probes(:) - p;
    else
        d = abs(probes(:) - p);
    end
    ds = [ds d];
    amps = [amps A(:)];
    plot(d, A(:),'ro','buttondownfcn',{@HitDistance,cid(k)});
    hold on;
    end
end
subplot(1,2,2);
hold off;

alld = unique(ds(:));
for j = 1:length(alld)
    id = find(ds == alld(j));
    a(j) = nanmean(amps(id));
    if abs(alld(j)) > 0
    subplot(1,2,2);
    [y,x] = hist(amps(id),nbin);
    plot(x,y,'color',colors{j});
    hold on;
    end
end
subplot(1,2,1);
plot(alld, a,'ko','markerfacecolor','k');
dx = 0:max(alld)/100:max(alld);
%
%Exponential fit is approxiamte. Even a true exponential decay 
%gets distorted by discretized sampling, since the true peak may fall 
%between two probes. see matlab/sims/DistanceVoltage
baseline = max([min(a) 0]);
A = min(baseline) + (1-baseline) .* exp(-dx * 1.1);
plot(dx, A,'k-','linewidth',2);
refline(0);

if length(amps) < 20
ds = [];
amps = [];
hold on;
for k = 1:length(nid)
    C = R{nid(k)};
    if ~isempty(C)
    probes = 1:length(C.amp);
    [a,p] = max(C.amp);
    A = C.amp ./a;
    ds = [ds  abs(probes(:) - p)];
    amps = [amps A(:)];
    end
end
plot(ds, amps,'bo');
alld = unique(ds(:));
for j = 1:length(alld)
    id = find(ds == alld(j));
    a(j) = mean(amps(id));
end
hold on;
plot(alld, a,'ko','markerfacecolor','k');
end

function HitDistance(a, b, cell)

DATA = GetDataFromFig(a);
X = DATA.cells{cell};
GetFigure('Amps');
hold off; 
plot(X.amp./max(X.amp));
if isfield(X,'dotA')
    hold on;
    plot(X.dotA./max(X.dotA),'r');
end
refline(0)
if isfield(DATA.details,'dir')
    if iscell(DATA.details.dir)
        d = fileparts(DATA.details.dir{X.dirid});
    elseif isdir(DATA.details.dir) %single session of data
        d = DATA.details.dir;
    else    %single session of data
        d = fileparts(DATA.details.dir);
    end
    cfile = sprintf('%s/Expt%dClusterTimes.mat',d,X.eid);
    if exist(cfile)
        load(cfile);
        GetFigure('MeanSpike');
        PlotMeanSpike(Clusters{X.probe},DATA.plotargs{:});
    end
end

function PlotDprimes(R, cid, varargin)

otherprobe = 0;
scatter = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'oprobe',5)
        otherprobe = 1;
    elseif strncmpi(varargin{j},'scatter',5)
        scatter = 1;
    end
    j = j+1;
end

for k = 1:length(cid)
    C = R{cid(k)};
    if isfield(C,'chspk')
        p = find(C.chspk == C.probe);
    else
        p = C.probe;
    end
    if size(C.vdprime,1) ~= size(C.vdiff,1)
    else
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
    end
    if size(R{cid(k)}.vdprime,1) > p       
%        vdiff = R{cid(k)}.MeanSpike.ms-R{cid(k)}.MeanSpike.mu;
        s = abs(smooth(R{cid(k)}.vdprime(p+1,:),2,'gauss'));
        id = 1+LocalMaxima(abs(s));
        [a,b] = sort(s(id),'descend');
        adpi(k,1:length(id)) = id(b)-C.triggerpt;
        adp(k,1:length(id)) = a;
    else
        adpi(k,:) = 0;
        adp(k,:) = 0;
    end
    if p > 1
        s = abs(smooth(R{cid(k)}.vdprime(p-1,:),2,'gauss'));
        id = 1+LocalMaxima(abs(s));
        [a,b] = sort(s(id),'descend');
        bdpi(k,1:length(id)) = id(b)-C.triggerpt;
        bdp(k,1:length(id)) = a;
    else
        bdpi(k,:) = 0;
        bdp(k,:) = 0;
    end
end
if scatter
    subplot(1,1,1);
    gid = find(dpi(:,2) ~= 0);
    myscatter(dpi(:,1),dpi(:,2),'o','ids',cid,'buttonpress',@HitScatter);
    xlabel('best sample');
    ylabel('second best sample');
    return;
elseif otherprobe == 0
    myscatter(dpi(:,1),dp(:,1),'o','ids',cid,'buttonpress',@HitScatter);
    gid = find(dpi(:,2) ~= 0);
    myscatter(dpi(gid,2),dp(gid,2),'ro','ids',cid,'buttonpress',@HitScatter);
else
    n = min([size(adpi,2) size(bdpi,2)]);
    dpi = cat(1, adpi(:,1:n),bdpi(:,1:n));
    dp = cat(1, adp(:,1:n),bdp(:,1:n));
    id = find(dp(:,1) > 0);
    dp = dp(id,:);
    dpi = dpi(id,:);
    myscatter(adpi(:,1),dp(:,1),'o','ids',cid,'buttonpress',@HitScatter);
    myscatter(adpi(:,2),dp(:,2),'ro','ids',cid,'buttonpress',@HitScatter);
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
if otherprobe
    xlabel('Voltage Sample (adjacent probe)');
else
    xlabel('Voltage Sample');
end
ylabel('Dprime SU vs MU');
title('Spikes aligned at trigger point');

function id = LocalMaxima(x)

sgn = diff(sign(diff(x)));
id = find(sgn < 0);

function PlotPCs(R, id)
for j = 1:length(id)
    sV = R{id(j)}.V;
    V(j,1:length(sV)) = sV;
end
V(isnan(V)) = 0;
[a,b] =eig(cov(V));
pcs = V * a;
PlotND(pcs(:,37:40),[],'marker','o','callback',@HitScatter, id);

function PlotCSDPCs(R, id)
plotpeak = 2;
order = 1;
for j = 1:length(id)
    X = R{id(j)};
    smoothw = 1;
    [~,~,Z] = gauss2d(smoothw,-5:5);

    if isfield(X,'MeanSpike')
        [a, maxi] = max(X.amp);
        if maxi > 2 && maxi < length(X.amp)-1 %can calc csd
            csd = diff(squeeze(X.MeanSpike(1,:,:)));
            if order == 2
                csd = diff(csd);
            end
            csd = conv2(csd,Z,'same');
            [~, maxt(j)] = max(std(csd));
            [~, maxp(j)] = max(std(csd,[],2));
            [xm, hsd(j)] = Hist2Gauss(std(csd));
            [xm, vsd(j)] = Hist2Gauss(std(csd,[],2));

            peakcsd(j) = csd(maxp(j),maxt(j));
        sV = csd(maxi-2:maxi,:);
        V(j,1:length(sV(:))) = sV(:);
        end
    end
end
V(isnan(V)) = 0;
[a,b] =eig(cov(V));
pcs = V * a;
PlotND(pcs(:,37:40),[],'marker','o','callback',@HitScatter, id);

GetFigure('CSDparams');
subplot(1,2,1);
plot(hsd,vsd,'o', 'ButtonDownFcn', {@HitScatter});
xlabel('Time SD');
ylabel('Spread');
subplot(1,2,2);
plot(maxp,peakcsd,'o', 'ButtonDownFcn', {@HitScatter});
xlabel('probe');
ylabel('CSD peak');

function PlotShapes(R, id, varargin)

plotimage = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'image',4)
        plotimage = 1;
    end
    j = j+1;
end

for j = 1:length(id)
    V = R{id(j)}.V;
    V = V./std(V);
    [a,b] = min(V);
    minv(j) = a;
    minpt(j) = b;
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

if plotimage
    [a, sid] = sort(postshape);
    for j = 1:length(id);
        V = R{id(sid(j))}.V;
        Im(j,1:length(V)) = V;
    end
    subplot(1,1,1);
    imagesc(Im);
    return;
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
x = get(a,'xdata');
y = get(a,'ydata');
if nargin == 2 %not plotted one at at time
    datapts = x + i * y;
    p = get(gca,'currentpoint');
    [~,id] = min(abs(p(1,1)+i*p(1,2)-datapts));
    C = DATA.cells{id};
    x = x(id);
    y = y(id);
else
    C = DATA.cells{id};
end
fit = [];

if isfield(C,'monkey')
    fprintf('%s',C.monkey)
end
if isfield(C,'tag')
fprintf('%s',C.tag)
cname = C.tag;
   for j = 1:length(DATA.details.dir)
       if ~isempty(strfind(DATA.details.dir{j},C.tag))
           cname = DATA.details.dir{j};
       end
   end
aid = 0;
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
fprintf('at %.2f,%.2f: ',x,y);
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
else
    [xm, sd] = Hist2Gauss(C.amp,[1:length(C.amp)],max(C.amp)/8);
    xi = 1:0.1:length(C.amp);
    fitted = Gauss([xm sd max(C.amp)],xi);
    hold on;
    plot(xi,fitted,'g');
end
title(sprintf('E%dP%d',C.eid,C.probe));
subplot(2,1,2);
hold off;
plot(C.V,'r');
if isfield(C,'muV')
    hold on;
    plot(C.muV,'color',[0.5 0.5 0.5]);
end
if isfield(C,'MeanSpike')
    cl = C.cluster;
    hold on;
    plot(squeeze(C.MeanSpike(:,C.probe,:))');
    X.MeanSpike.ms = squeeze(C.MeanSpike(cl,:,:));
    X = CopyFields(X,C,{'probe' 'cluster'});
    X.exptno = C.eid;
    X.mahal = C.dip;
    
    GetFigure('MeanSpike');
    PlotMeanSpike(X,DATA.plotargs{:});
    GetFigure('Spike');
end

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
        if length(Clusters) >= p
            X = Clusters{p};
            DATA.cells{id}.xy = X.xy;
            DATA.cells{id}.clst = X.clst;
            fprintf('Load  took %.2f\n',b.loaddur);
        end
    end
    if isfield(X,'xy')
        colors = mycolors('spkcolors');
        GetFigure('XYplot');
        if isfield(X,'eid')
        if aid > 0
            xid = DATA.details.nres(1,aid):DATA.details.nres(2,aid);
            xid = xid(xid <= length(DATA.cells)); %in case removed some
        else
            xid = 1:length(DATA.cells);
        end
        expts = CellToMat(DATA.cells(xid),'eid');
        probes = CellToMat(DATA.cells(xid),'probe');
        oid = find(expts == X.eid & probes == X.probe);
        if isfield(X,'dirid')
            dirs = CellToMat(DATA.cells(xid),'dirid');
            oid = intersect(oid, find(dirs == X.dirid));
        end
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
            eid = expts(find(need));
            probe = probes(find(need));
        else
            xy = X.xy;
            eid = expts(oid);
            probe = probes(oid);
        end
        else
            xy = X.xy;
        end
        if DATA.plot.density
        PlotND(xy,[],'idlist',X.clst,'density');
        else
        PlotND(xy,[],'idlist',X.clst,'colors',colors);
        end
        title(sprintf('E%dP%d',eid,probe));
        if DATA.plot.nrefit 
            [X,G] = Refit(DATA, X);
            GetFigure('FitDistance');
            newid = cluster(G, X.xy);
            plot(X.mind,X.maxd,'ro-');
            xlabel('smallest Mahal distance');
            ylabel('largest distance');
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

function [xm, sd] = Hist2Gauss(y,x, crit)

if nargin ==1
    x = 1:length(y(:));
end
if nargin < 3
    crit = max(y)/10;
end
id = find(y >crit);
x = x(id);
y = y(id);
xm = sum(x(:) .* y(:))./sum(y(:));
sd = sqrt(sum((x(:)-xm).^2 .* y(:))./sum(y(:)));


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
            [a,b] = GMdprime(G{nd});
            mind(nd) = min(b.d(:));
            maxd(nd) = max(b.d(:));
            ll(nd) = G{nd}.NlogL;
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
        X.ll = ll;
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
    
    
function res = AddMeanSpike(C, varargin)    
details = [];
tag = [];
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'dir')
        details = varargin{j};
    end
    j = j+1;
end

for j = 1:length(C)
    dirtag{j} = C{j}.tag;
    tags{j} = sprintf('%s.%d.%d',C{j}.tag,C{j}.eid,C{j}.probe);
end
[a,b] = Counts(tags);
dirs = unique(dirtag);
for j = 1:length(dirs)
    id = find(strcmp(dirs{j},dirtag));
    if isfield(C{id(1)},'dirid')
        dirid = C{id(1)}.dirid;
        name = sprintf('%s/CellList.mat',details.dir{dirid});
    else
        name = sprintf('\\data/lem/%s/CellList.mat',dirs{j});
        dirid = 0;
    end
    if exist(name,'file');
    load(name);
    for k = 1:length(id)
        eid = C{id(k)}.eid;
        if size(CellList,1) >= eid
            xid = find(CellList(eid,C{id(k)}.probe,:) == C{id(k)}.cell);
        else
            xid = [];
        end
        if ~isempty(xid)
              cid(id(k)) = xid(1);
              if ~isfield(C{id(k)},'cluster')
                  C{id(k)}.cluster = xid(1);
              end
        elseif ~isfield(C{id(k)},'cluster')
                  C{id(k)}.cluster = 0;
        end
    end
    else
        fprintf('No CellList %s\n',name);
    end
end
for j = 1:length(C)
    if ~isempty(details)
        name = sprintf('%s/Expt%dClusterTimes.mat',details.dir{C{j}.dirid},C{j}.eid);
    else
        name = sprintf('\\data/lem/%s/Expt%dClusterTimes.mat',C{j}.tag,C{j}.eid);
    end
    if exist(name)
        load(name);
        Cl = Clusters{C{j}.probe};
        nc = 1;
        clear M;
        M(nc,:,:) = Cl.MeanSpike.ms;
        if isfield(Cl,'next')
        for k = 1:length(Cl.next)
            if isfield(Cl.next{k},'MeanSpike')
                nc = nc+1;
                M(1+k,:,:) = Cl.next{k}.MeanSpike.ms;
            end
        end
        end
        C{j}.MeanSpike = M;
    else
        fprintf('Cant file Clusters for %s%s\n',C{j}.monkey,C{j}.tag);
    end
end
res = C;
    
function [res, details] = OncePerCell(res, varargin)
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
            monkey = GetMonkeyName(details.dir{k});
            [c, dirname] = fileparts(details.dir{k});
            if b > a
                c = OncePerCell(res(1+a:b),'tag',dirname);
                if isempty(c)
                    fprintf('No cells in %s\n',details.dir{k});
                else
                    for ci = 1:length(c)
                        c{ci}.dirid = k;
                        c{ci}.monkey = monkey;
                    end
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
    if min(dips(id,1) < -0.5)
        [a,b] = min(dips(id,1));
    else
        [a,b] = max(dips(id,1));
    end
    goodid(j) = id(b);
    if ~isempty(tag)
        res{goodid(j)}.tag = tag;
    end
end
res = res(goodid);
end

function [result, details] = BuildMeanData(basedir,monkey, varargin)

if strcmp(monkey,'lem')
    d = mydir({[basedir '/lem/M18*'] [basedir '/lem/M18*'] [basedir '/lem/M2*']});
end
[result, details]  = CalcCellSize({d.name},varargin{:});

