function D = FullV2LFP(name,varargin)
%D = FullV2LFP(name,....)
%    if name is a directory, Builds LFP files from each FullV file
%D = FullV2LFP(FullV,....)
%    where FullV is a loaded FullV structure, Builds LFP files from that
%    struct
%...,'expts',exlist)   only makes LFP for expts in exlist (ints)
%...,'parallel') uses parfor to loop over expts
%...,'kernel',K) use K as the smoothing kernel to lowpass the FullV record
%                default is a Gauassian with SD of 24 samples = 200Hz sigma
%                for 30KHz sampling
%...,'decimate',n) take every nth sample after smoothing (default 60);


mkexpts = [];
G = [];
savelist = 0;
rebuild = 0;
LFPversion = 1.0;
parallel = 0;
LFP.sd = 24; %24 at 30KHz = SD of 0.8ms, = SD of 200Hz in Freq
LFP.decimate = 60;
checkversion = 0;
V = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'checkversion',10)
        checkersion = 1;
    elseif strncmpi(varargin{j},'decimate',5)
        j = j+1;
        LFP.decimate = varargin{j};
    elseif strncmpi(varargin{j},'expts',5)
        j = j+1;
        mkexpts = varargin{j};
    elseif strncmpi(varargin{j},'kernel',5)
        j = j+1;
        G = varargin{j};
    elseif strncmpi(varargin{j},'parallel',5)
        parallel = 1;
    elseif strncmpi(varargin{j},'rebuild',5)
        rebuild = 1;
    end
    j = j+1;
end


if iscellstr(name) %list of dirs
    for j = 1:length(name)
        D{j} = FullV2LFP(name{j});
    end
    return;
end

if iscell(name)
    lfp = name;
    goodch = zeros(1,length(lfp));
    LFP.rawlfp = [];
    for p = 1:length(lfp)
        if isfield(lfp{p},'lfp')
            rawlfp(1:length(lfp{p}.lfp),p) = lfp{p}.lfp;
            if isfield(lfp{p},'goodch')
                goodch(p) = lfp{p}.goodch;
            else
                goodch{p} = 1;
            end
            if isfield(lfp{p},'t')
                LFP.t = lfp{p}.t; %should always be the same
            end
        end
        LFP = CopyFields(LFP,lfp{p},'samper','version');
        if isfield(lfp{p},'name')
            LFP.name = lfp{p}.name;
        end
    end
    intscale(3) = max(max(abs(rawlfp(:,find(goodch)))));
    intscale(2) = 32000;
    intscale(1) = intscale(2) ./ intscale(3);
    LFP.rawlfp = int16(round(rawlfp.*intscale(1)));
    LFP.intscale = intscale;
    D.trange = minmax(LFP.t);
    outname = sprintf('%s/Expt%d.lfp.mat',LFP.name,mkexpts(1));
    fprintf('Saving %dx%d samples (max %d), to %s\n',size(rawlfp,1),size(rawlfp,2),max(LFP.rawlfp(:)),outname);
    ts = now;
    save(outname,'LFP','-v7.3');
    fprintf('took %.2f sec\n',mytoc(ts));
    return;
end


if isempty(G)
    G = Gauss(LFP.sd,-100:100);
end
LFP.kernel = G;

if isstruct(name) && isfield(name,'V')
        lfp = conv(name.V,G,'same');
        LFP.lfp = downsample(lfp,LFP.decimate); %down to 500Hz;
        if ~isfield(name,'t')
            name.t = BuildFullVt(name);
        end
        LFP.t = downsample(name.t,LFP.decimate);
        LFP.samper = name.samper.*LFP.decimate; 
        LFP.name = fileparts(name.name);
        LFP.version = LFPversion;
        if name.skew > 5 && name.kurtosis > 1000
            LFP.goodch = 0;
        else
            LFP.goodch =1;
        end
        D = LFP;
    return;
end

ts = now;
LFP.name = name;
LFP.version = LFPversion;
lfpage = [];
d = dir([name '/Expt*.lfp.mat']);
for j = 1:length(d)
    ex = sscanf(d(j).name,'Expt%d');
    lfpage(ex) = d(j).datenum;
    D.lfpnames{ex} = d(j).name;
    lfpindex(ex) = j;
end

d = mydir([name '/Expt*FullV.mat']);
if isempty(d)
    D.loadsum = 0;
    D = AddError(D,'No FullV  files in %s\n',name);
    return;
end

neednew = [];
needbuild = [];
ts = now;
for j = 1:length(d)
    
    if isempty(regexp(d(j).name,'Expt[0-9]*.p[0-9]*FullV.mat'))
        good(j) = 0;
    else
        good(j) = 1;
    end
    expt(j) = sscanf(d(j).filename,'Expt%d');
    if length(lfpage) >= expt(j)
        if lfpage(expt(j)) < d(j).datenum
            neednew(expt(j)) = 1;
            fprintf('Need to Rebuild %s\n',d(j).name);
        elseif checkversion %check to see if old version of FullV2LPF was used.  This is slow
            x = matfile(d(j).name);
            if sum(strcmp('LFPHeader',fields(x)))
                H = x.LFPHeader;
                if isfield(H,'version')
                    fprintf('%s made by version %.2f\n',H.version);
                else
                    fprintf('%s made No Version record\n',H.version);
                end
            end
        end
    else
        needbuild(expt(j)) = 1;
    end
end
if sum(needbuild)
    fprintf('%s: No LFP Files for Expts%s\n',name,sprintf(' %d',find(needbuild)));
end
neednew = union(find(neednew),find(needbuild));
fprintf('Checking %s took %.2f\n',name,mytoc(ts))
expt = expt(good > 0);
d = d(good > 0);

G = G./sum(G);
expts = unique(expt);
if isempty(mkexpts)
    mkexpts = expts;
    savelist = 1;
else
    mkexpts = intersect(mkexpts,expts);    
end

if rebuild == 0
    mkexpts = intersect(neednew,expts);
end

if isempty(mkexpts)
    fprintf('All LFP files are up to date\n');
    D.duration = mytoc(ts);
    D.loadsum = 0;
    return;
end
ts = now;
sumload = 0;
if parallel
    parfor j = (1:length(mkexpts))
        Di{j} = MakeLFP(mkexpts(j),expt, d, LFP);
    end
else
    for j = 1:length(mkexpts)
        Di{j} = MakeLFP(mkexpts(j),expt, d, LFP);
    end
end
D = Di;


loadsum = 0;
duration = 0;
for j = 1:length(Di)
    if isfield(Di{j},'sumload') && isfield(Di{j},'duration')
        loadsum = loadsum + Di{j}.sumload;
        duration = duration + Di{j}.duration;
    end
end

fprintf('%.2f/%.2f (%.1f%%) of total time was load\n',loadsum,duration, loadsum./duration);
if savelist
    outname = sprintf('%s/LFPlist.mat',name);
    if exist(outname)
        X = load(outname);
    else
        X.D = {};
    end
    if isfield(X,'D') && isstruct(X.D) %old format
        X.D = D;
    else
        for j = 1:length(D)
            X.D{D{j}.exptno} = D{j};
        end
    end
    save(outname,'-struct','X');
end

function D = MakeLFP(e, expt, d, LFP)
ts = now;
id = find(expt == e);
    clear rawlfp;
    G = LFP.kernel;
    D.sumload = 0;
    for p = 1:length(id)
        fprintf('Loading %s ',d(id(p)).name);
        V = LoadFullV(d(id(p)).name,'nohighpass','addtime');
        D.sumload = D.sumload + V.loadtime;
        probe = Name2Probe(d(id(p)).name);
        if probe ~= V.chspk
            cprintf('red','Probe name mismatch\n');
        end
        if ~isfield(V,'skew')
            D.skew(probe) = skewness(V.V);
            D.kurt(probe) = kurtosis(V.V);
        else
            D.skew(probe) = V.skew;
            D.kurt(probe) = V.kurtosis;
        end
        D.max(probe) = max(V.V).*1000000;
        if D.skew(probe) > 5 && D.kurt(probe) > 1000
            goodch(probe) = 0;
        else
            goodch(probe) = 1;
        end
        fprintf('skew %.2f Kurtosis %.1f Max %.1fuV %d pts\n',D.skew(probe),D.kurt(probe),D.max(probe),length(V.V));
        if ~isfield(V,'highpass')
            D = AddError(D,'-show','%s highpass not recorded ? filtered\n',d(id(p)).name);
            goodch(probe) = -1;
        elseif V.highpass == 0 %data does not need offline high pass for spike channels...
            D = AddError(D,'-show','%s Looks like low F were filtered out\n',d(id(p)).name);
            goodch(probe) = -2;
        end
        lfp = conv(V.V,G,'same');
        lfplen(probe) = length(lfp);
        if p > 1
            gid = setdiff(find(goodch>0),probe); %current probe may bee in goodch at this point
            if isempty(gid)
                gid = setdiff(find(goodch>0 | goodch == -2 | goodch==-1),probe); %current probe may bee in goodch at this point
            end
                
            if length(unique(lfplen(gid))) == 1
                npts = unique(lfplen(gid));
            else
                D = AddError(D,'%s P%d Lengths dont match',d(id(p)).name,probe);
                npts = round(mean(lfplen(gid)));
            end
        else
            npts = lfplen(probe);
        end
        if lfplen(probe)  < npts
            D = AddError(D,'%s P%d missing data %d points, array is %d',d(id(p)).name,probe,lfplen(p),npts);
            goodch(probe) = -3;
        elseif lfplen(probe)  > npts
            D = AddError(D,'%s P%d %d points, array is %d',d(id(p)).name,probe,lfplen(p),npts);
             goodch(probe) = -4;
        else
            rawlfp(:,probe) = decimate(lfp,LFP.decimate); %down to 500Hz;
        end
    end
    if sum(goodch > 0)
        intscale(3) = max(max(abs(rawlfp(:,find(goodch>0)))));
    else
        intscale(3) = max(max(abs(rawlfp(:,find(ismember(goodch,[-2 -1 1]))))));
    end
    intscale(2) = 32000;
    intscale(1) = intscale(2) ./ intscale(3);
    LPF.badchannels = find(goodch <=0);
    LFP.rawlfp = int16(round(rawlfp.*intscale(1)));
    LFP.intscale = intscale;
    LFP.t = downsample(V.t,LFP.decimate);
    LFP.samper = V.samper.*LFP.decimate; % = in sec, not silly units
    LFPHeader.goodch = goodch;
    LFPHeader.version = LFP.version;
    LFPHeader.nsamples = lfplen;
    
    if isfield(D,'errs')
        LFP.errs = D.errs;
        LFPHeader.errs = D.errs;
    end
    LFPHeader = CopyFields(LFPHeader,LFP,'errdata','badchannels','kernel','samper','decimate', 'nsamples');
    D.Header = LFPHeader;
    D.trange = minmax(LFP.t);
    D.duration = mytoc(ts);
    D.exptno = e;

    ts = now;
    outname = sprintf('%s/Expt%d.lfp.mat',LFP.name,e);
    fprintf('Saving %d samples (max %d), to %s ',length(rawlfp),max(LFP.rawlfp(:)),outname);
    save(outname,'LFP','LFPHeader','-v7.3');
    D.savedur = mytoc(ts);
    fprintf('took %.2f sec\n',D.savedur);    

