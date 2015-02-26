function AllExpt = MakeAllExpts(name, varargin)
% MakeAllExpts(name, varargin) reads in Expt files for each probe, and
% combines them into an AllExpts struct. For older files make from swatches
%
% e.g. 
% PlotAllProbes('lemM025.c1.image.DCORRC.mat')
%
%
probes = 1:24;
cluster = 1;
plottype = 0;
delay = 0.05;
yvals = 1;
xvals = 1;
showblank = 0;
scales = ones(24,1);
timerange = [200 2000 500 1300];
newscale = 0;
zvals = [];
isrc = [];
latency = 500;
mintrials = 10;
PLOTLFPDIFF=6;
BLANKRESP = 7;
FRAMERESP = 8;
LFPFREQ = 9;
LINEPLOT = 10;
PLOTLFPEIG = 11;
PLOTMONOC = 12;
PLOTONESTIM = 13;
PLOTONESTIMMU = 14;
XTPROBELFP = 15;
XTPROBELFPMU = 16;
FINDINVERSION = 17;
LFPPWRDIFF = 18;
STIMVAR = 19; 
STIMLATENCY = 20;
BLANKVAR = 21; %Stimvar and Blank plots.
LFPBANDPWR = 22;
MUVAR = 23;
CSDSUM = 24;
LFPSIGPWR = 25;
RESPXBLANK = 26;
PLOTLFPEIGVECTORS = 27;
LFPTRIAL = 28;
OCULARITY = 29;

setlfpsign = 2; %0 use blank resp, 1 use continuity, 3 use blank * resp
saveresult = 0;
fixscale = 1;
calcsptrig = 0;
nr=2;
nc=1;
pn = [1 2];
freqs = 2:50;
freqbands = [2 10 30 60 90];
freqrange = [1 100];
latency = 500;
lineplot = 0;
sumy = 0;
spacing = 0.2;
argon = {};
rcnmin = 10;  %get all the data first, worry about n when collapsing.
figures = {'PlotAllProbes', 'PlotAllb'};
sptrig = [];
plotargs = {};
csdsk = [2 2]; %smoothing kernel for CSD
interpolate = 0;
sumx = 0;
zoomstart = [];


j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'save',4)
        saveresult = 1;
    end
    j = j+1;
end

if isstruct(name)
    rc = name;
% if there is an 'lfp' matrix, uset that. Otherwise, if this is an rc file
% (output from PlotRevCOrAny, build the matrix
  if ~isfield(rc.Header,'exptype') & isfield(rc.Header,'expname')
      if strfind(rc.Header.expname,'OpRC')
          rc.Header.exptype = 'OPRC';
      end
  end
      
elseif ischar(name) && isdir(name)
    result = ProcessDir(name, saveresult, argon{:});
    return;
end

if isempty(isrc) && ischar(name)
    if strfind(name,'RC')
        isrc = 1;
    else
        isrc = 0;
    end
end

cs = ['c' num2str(cluster) '.'];
ts = 34;
xvals = [];
rcargs = {'nmin' rcnmin};
ispsych = 0;
if strfind(name,'DCORRC') 
    rcargs = {rcargs{:} 'psych' 'yvals' [0 0.06] 'nmin', 10};
    ispsych =1;
end
if strfind(name,'ORBW') 
    ispsych =1;
end
if strfind(name,'BOPRC') 
    rcargs = {rcargs{:} 'twoslice'};
end
clear LFP;

%Load the LFP file first, so we can build avgs
lfpname = strrep(name,['.' cs],'.lfp.');
if exist(lfpname,'file')
     load(lfpname);
else
     fprintf('No file %s\n',lfpname);
end


Expts = {};
for j = probes;
    ename = strrep(name,['.' cs],['.p' num2str(j) cs]);
    if exist(ename,'file')
        load(ename);
        Expt.Header.probe = j;
        Expts{j} = Expt;
        clear C;
        [ed, eds] = GetEval(Expt,'ed');
        if isrc | Expt.Header.rc == 1
            [res, bExpt] = PlotRevCorAny(Expt,'sdfw',166,'box',rcargs{:});
            result.type = res.type;
            if isfield(res,'sdfs')
                if ~isnan(res.bestdelay)
                    ts = res.bestdelay;
                else
                    ts = 50;
                end
                means(j) = mean(res.y(:));
                resps(j,:) = res.y(:,1,ts)./means(j);
                result.spkn = res.sdfs.n;
                id = find(res.times > 166/2);
                ta = id(2);
                xvals = res.x;
                for iy = 1:size(res.sdfs.s,2)
                    for ix = 1:size(res.sdfs.s,1)
                        if res.sdfs.n(ix,iy) > 0 & ~isempty(res.sdfs.s{ix,iy})
                            result.tresps(:,ix,iy,j) = res.sdfs.s{ix,iy}(ta:end)./means(j);
                            result.vtimes(j,:) = res.delays;
                            result.besttimes(j,:) = ts;
                        else
                            result.tresps(:,ix,iy,j) = 0;
                        end
                    end
                    yvals(iy) = res.sdfs.y(1,iy);
                end
                id = find(res.sdfs.extraval == -1009);
                if length(id) == 2
                    result.blankresp(:,j) = (res.sdfs.extras{id(1)}.sdf(ta:end)+res.sdfs.extras{id(2)}.sdf(ta:end))./( 2.*means(j));
                elseif length(id) == 1
                    result.blankresp(:,j) = res.sdfs.extras{id}.sdf(ta:end)./means(j);
                end
                result.times = res.times(ta:end);
                if isfield(bExpt.Trials,'Dc');
                    result.Dc = [bExpt.Trials.Dc];
                end
                if isfield(bExpt.Trials,'ori');
                    result.or = [bExpt.Trials.ori];
                end
                if ispsych
                    result.RespDir = [bExpt.Trials.RespDir];
                    bres = PlotExpt(Expt,'condense');
                    id = find(bres.x(:,1) >= 0.1);
                    count = sum(bres.means(id,:) .* bres.n(id,:))./sum(bres.n(id,:));
                    res.cp.sigcount = count; %mean count for prefs and nulls
                    res.cp.sigcounts = bres.means;
                    res.cp.sigsd = bres.sd;
                    if count(1) > count(end)
                        res.cp.prefbycount = mean(bres.y(:,1));
                    else
                        res.cp.prefbycount = mean(bres.y(:,end));
                    end
                    result.cp(j) = res.cp;
                end
                if strmatch(res.type{1},{'Op' 'Pp'})
                    result.fit{j} = FitExpt(res,'plotfit');
                end

                for k = 1:length(bExpt.Header.Clusters)
                    C(k) = GetClusterInfo(bExpt.Header.Clusters{k},j);
                    if isempty(C(k).autocut)
                        if isfield(bExpt.Header,'Combineids')
                            fprintf('No autocut in Expt %d, Probe %d\n',bExpt.Header.Combineids(k),j);
                        else
                            fprintf('No autocut in Expt Sub%d Probe %d\n',bExpt.Header.Combined(k),j);
                        end
                    end
                end
                result.Cluster(j).autocut = [C.autocut];
                result.Cluster(j).dprime = [C.dprime];
 % sptrig dimesnions are time, LFP probe yval, spike probe
                if calcsptrig
                    if strmatch(Expt.Stimvals.e3, 'me')
                        [avg, a] = ExptTrigLFP(Expt,LFP,'split','me',[-1 0 1],'nohist');
                        sptrig.lfp(:,:,:,j) = avg;
                    else
                        [sptrig.lfp(:,:,:,j), a] = ExptTrigLFP(Expt,LFP,'nohist');
                    end
                    sptrig.lfpn(j,:) = [a.nspk];
                    sptrig.lfptimes = a(1).times./10;
                end
            end
        else
            Expt = FillTrials(Expt,Expt.Stimvals.et); %in case
            res = PlotExpt(Expt);
            if ispsych
                pres = PlotExpt(Expt,'psych','cp','noplot');
                res.cp = pres(1).cp;
                res.cp.gcp = pres(1).gcp.cp;
            end
            if isfield(res,'means')
                xvals = res.x;
                yvals = res.y(1,:);
                means(j) = mean(res.means(:));
                resps(:,:,j) = res.means./means(j);
                if calcsptrig
                    if strmatch(Expt.Stimvals.e2, 'ce')
                        [avg, a] = ExptTrigLFP(Expt,LFP,'split','ce',[-1 1],'nohist');
                    elseif strmatch(Expt.Stimvals.e2, 'me')
                        [avg, a] = ExptTrigLFP(Expt,LFP,'split','me',[-1 0 1 NaN],'nohist');
                    else
                        [avg, a] = ExptTrigLFP(Expt,LFP,'nohist');
                    end
                    if ~isempty(avg)
%the size of a CAN vary with probe. If an Expt combines different Expts, and one probe is only 
% cut for some expts, can end up with different # stimuli
                    for k = 1:size(avg,3)
                        sptrig.lfp(:,:,k,j) = avg(:,:,k);
                        sptrig.lfpn(j,k) = [a(k).nspk];
                    end
                    sptrig.lfptimes = a(1).times./10;
                    end
                end
                if strmatch(res.type,{'Op' 'Pp'})
                    result.fit{j} = FitExpt(res);
                end
                result.type = res.type; 
            end
        end
    for k = 1:length(Expt.Header.Clusters)
        C(k) = GetClusterInfo(Expt.Header.Clusters{k},j);
        if j <= length(Expt.Header.Clusters{k})
            Clusters{k}{j} = Expt.Header.Clusters{k}{j};
        else
            Clusters{k}{j} = Expt.Header.Clusters{k}{j};
        end
    end
    result.isrc = isrc;

    result.Cluster(j).autocut = [C.autocut];
    result.Cluster(j).dprime = [C.dprime];
    if isfield(Expt.Header,'SpkStats')
        result.Cluster(j).SpkStats = Expt.Header.SpkStats;
    end
    end
end
if isempty(Expts)
fprintf('No Matching Expts\n');
AllExpt = {};
return;
end

result.Stimvals.fz = Expt.Stimvals.fz;
result.Stimvals.et = Expt.Stimvals.et;
result.Stimvals.e2 = Expt.Stimvals.e2;
result.Stimvals.Ro = Expt.Stimvals.Ro;


clst = regexprep(name,'.c1.*','.cells.mat');
isolation = 5; %%quality must exceed this to count
if exist(clst,'file')
nsux =0;
nmux = 0;
    Expt = FillTrials(Expt,'dur');
    load(clst);
    for j = 1:size(CellList,1)
        probe = CellList(j,:);
        quality = CellQuality(j,:);
        if exist('CellListCluster','var')
          cluster =  CellListCluster(j,:);
        else
            cluster = ones(size(CellQuality(j,:)));
        end
        cExpt = Expt;
        muExpt = Expt;
        trials = [Expt.Trials.Trial];
        res = [];
        ures = [];
%
% CellList can be shorter that # trials, if some are undefined.
        trials = trials(find(trials < size(CellList,2)));
        nsu = 0;
        nmu = 0;
        dur = mean([Expt.Trials.dur]);
        if sum(probe(trials)>0) > mintrials
        for t = 1:length(trials)
            if probe(trials(t)) > 0 && quality(trials(t)) > isolation
                if cluster(trials(t)) <= 1
                    spks = Expts{probe(trials(t))}.Trials(t).Spikes;
                    cExpt.Trials(t).Spikes = spks;
                    cExpt.Trials(t).count = sum(spks > 500 & spks < dur+500);
                else
                    cid = find(Expts{probe(trials(t))}.Trials(t).Ocodes == cluster(trials(t)));
                    spks = Expts{probe(trials(t))}.Trials(t).OSpikes(cid);
                    cExpt.Trials(t).Spikes = spks;
                    cExpt.Trials(t).count = sum(spks > 500 & spks < dur+500);
                end
                nsu = nsu+1;
            end
            if probe(trials(t)) > 0 && quality(trials(t)) <= isolation
                muExpt.Trials(t).Spikes = Expts{probe(trials(t))}.Trials(t).Spikes;
                nmu = nmu+1;
            end
        end
        cid = find(probe(trials) > 0 & quality(trials) > isolation);
        cExpt.Trials = cExpt.Trials(cid);
        cExpt.probes = probe(trials(cid));
        cl = unique(cluster(trials(cid)));
        if length(cl) == 1
            cExpt.Header.Clusterid = cl;
        else
            cExpt.Header.Clusterid = cluster(trials(cid));
        end
        bk = Expt.Header.BlockStart;
        bk = [bk Expt.Trials(end).Trial+1];
        for b = 1:length(bk)-1
            tid = bk(b):bk(b+1)-1;
            id = find(ismember([cExpt.Trials.Trial],tid));
            if isempty(id)
                cExpt.Header.BlockCount(b) = NaN;
                cExpt.Header.blockprobe(b) = NaN;
            else
                cExpt.Header.BlockCount(b) = mean([cExpt.Trials(id).count]);
                p = mode(probe([cExpt.Trials(id).Trial]));
                cExpt.Header.Clusters{b} = Expts{p}.Header.Clusters{b};
                cExpt.Header.blockprobe(b) = p;
            end
        end
        cExpt.Header.cellnumber = j;
        cExpt.Header.probe = mean(probe(trials(cid)));
        cExpt.Header.quality = mean(quality(trials(cid)));

        x = GetProbeSep(cExpt.Header);
        if ~isnan(x)
            cExpt.Header.probesep = x;
        end
        muid = find(probe(trials) > 0 & quality(trials) <= isolation);
        muExpt.Trials = muExpt.Trials(find(probe(trials) > 0 & quality(trials) <= isolation));
        muExpt.probes = probe(trials(find(probe(trials) > 0 & quality(trials) <= isolation)));
        muExpt.Header.probe = median(probe(trials(muid)));
        muExpt.Header.cellnumber = j;
        muExpt.Header.quality = mean(quality(trials(muid)));
        if isrc
            if length(cExpt.Trials)
            [res, bExpt] = PlotRevCorAny(cExpt,'sdfw',166,'box',rcargs{:});
            if calcsptrig
 % lfptrig dimensions are time, LFP probe, yval    
                [T.lfptrig, a] = ExptTrigLFP(cExpt,LFP,'nohist');
                T.triglfpnspk = [a.nspk];
                T.triglfptimes = a(1).times./10;
                T.probe = mean(cExpt.probes);
%                cExpt.triglfp = T;
            end
            if ispsych
                bres = PlotExpt(cExpt,'condense');
                id = find(bres.x(:,1) >= 0.1);
                count = sum(bres.means(id,:) .* bres.n(id,:))./sum(bres.n(id,:));
                res.cp.sigcount = count;
                if count(1) > count(end)
                    res.cp.prefbycount = mean(bres.y(:,1));
                else
                    res.cp.prefbycount = mean(bres.y(:,end));
                end
            end
            if strmatch(res.type{1},{'Op' 'Pp'})
                res.fit = FitExpt(res);
            end
            end

            if nmu > mintrials
            [ures, uExpt] = PlotRevCorAny(muExpt,'sdfw',166,'box', rcargs{:});
            end
        else
            if length(cExpt.Trials) 
                cExpt = FillTrials(cExpt,Expt.Stimvals.et); %in case
                res = PlotExpt(cExpt);
                if ispsych
                    pres = PlotExpt(cExpt,'psych','cp','noplot');
                    res.cp = pres(1).cp;
                    res.cp.gcp = pres(1).gcp.cp;
                end
                if calcsptrig
                    if strmatch(Expt.Stimvals.e2, 'ce')
                        [avg, a] = ExptTrigLFP(cExpt,LFP,'split','ce',[-1 1],'nohist');
                    else
                        [avg, a] = ExptTrigLFP(cExpt,LFP,'nohist');
                    end
                    if ~isempty(avg)
                        T.lfptrig = avg;
                        T.triglfpnspk = [a.nspk];
                        T.triglfptimes = a(1).times./10;
                        T.probe = mean(cExpt.probes);
                        T.probes = cExpt.probes;
                        if strmatch(res.type{1},{'Op' 'Pp'})
                            res.fit = FitExpt(res);
                            T.type = res.type;
                        end
                    end
                end
                res = rmfield(res,'Data');
            end
              
        end
        if calcsptrig && exist('T','var')
            trigs(j) = T;
            sptrig.cells = trigs;
        end
        if length(cid)
        pt = find(diff(probe(trials(cid))) ~= 0);  % trials where probe changes
        res.probes = probe(trials(cid([1 pt])));
        res.probestep = trials(cid(pt));
        end
        res.cellnumber = j;
        outname = strrep(name,'.c1.',['.cell' num2str(j) '.']);
        res.Header = cExpt.Header;
        res.Header.filename = outname;
        res.Trials = [cExpt.Trials.Trial];
        if length([res.Trials])
            result.spkres{j} = res;
        else
            result.spkres{j} = [];
        end
        if ~isempty(ures)
            ures.Header = Expt.Header;
            ures.Trials = [uExpt.Trials.Trial];
            pt = find(diff(probe(trials(muid))) ~= 0);  % trials where probe changes
            ures.probes = probe(trials(muid([1 pt])));
            ures.probestep = trials(muid(pt));
            result.mures{j} = ures;
        end
        if nsu > mintrials
            nsux = nsux +1;
            suExpts{nsux} = cExpt; 
        end
        if nmu > mintrials
            nmux = nmux +1;
            muExpts{nmux} = muExpt; 
            muprobes(nmux) = muExpt.Header.probe;
        end
        if nsu < mintrials && nmu < mintrials
            nmux = nmux +1;
            muExpts{nmux} = muExpt; 
            muprobes(nmux) = muExpt.Header.probe;
        end
        end
    end
end

for j = 1:nsux
    E = suExpts{j};
    T = E.Trials;
    Headers(j).cellnumber = E.Header.cellnumber;
    Headers(j).probe = E.Header.probe;
    Headers(j).quality = E.Header.quality;
    Spikes{j}.trialid = [E.Trials.id];
    Spikes{j}.Trial = [E.Trials.Trial];
    for t = 1:length(T)
        Spikes{j}.Spikes{t} = T(t).Spikes;
        Spikes{j}.OSpikes{t} = T(t).OSpikes;
        Spikes{j}.Ocodes{t} = T(t).Ocodes;
    end
end
if nmux > 0
[a,b] = sort(muprobes);
muExpts = muExpts(b);
for j = 1:nmux
    k = j+nsux;
    E = muExpts{j};
    T = E.Trials;
    Headers(k).cellnumber = E.Header.cellnumber;
    Headers(k).probe = E.Header.probe;
    Headers(k).quality = E.Header.quality;
    Spikes{k}.trialid = [E.Trials.id];
    Spikes{k}.Trial = [E.Trials.Trial];
    for t = 1:length(T)
        Spikes{k}.Spikes{t} = T(t).Spikes;
        Spikes{k}.OSpikes{t} = T(t).OSpikes;
        Spikes{k}.Ocodes{t} = T(t).Ocodes;
    end
end
end
AllExpt.Spikes = Spikes
AllExpt.Header = Headers;
AllExpt.Expt = Expt;
if saveresult
    outname = strrep(name,'.c1.','');
    outname = strrep(outname, '.mat','.Cells.mat');
    save(outname,'AllExpt');
end

function C = GetClusterInfo(Cluster, probe)

C.autocut = [];
C.dprime = [];
if probe > size(Cluster,2)
    return;
end

if isfield(Cluster{1,probe}, 'autocut')
    C.autocut = Cluster{1,probe}.autocut;
else 
    C.autocut = 0; 
end

if isfield(Cluster{1,probe},'dprime')
    C.dprime = Cluster{1,probe}.dprime;
end


