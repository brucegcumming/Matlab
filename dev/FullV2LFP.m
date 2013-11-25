function D = FullV2LFP(name,varargin)
%D = FullV2LFP(name,varargin)


mkexpts = [];
G = [];
savelist = 0;
version = 1.0;
parallel = 0;
V = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'expts',5)
        j = j+1;
        mkexpts = varargin{j};
    elseif strncmpi(varargin{j},'kernel',5)
        j = j+1;
        G = varargin{j};
    elseif strncmpi(varargin{j},'parallel',5)
        parallel = 1;
    end
    j = j+1;
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
LFP.sd = 24; %24 at 30KHz = SD of 0.8ms, = SD of 200Hz in Freq
LFP.decimate = 60;

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
        LFP.version = version;
        if name.skew > 5 && name.kurtosis > 1000
            LFP.goodch = 0;
        else
            LFP.goodch =1;
        end
        D = LFP;
    return;
end

LFP.name = name;
d = mydir([name '/Expt*FullV.mat']);
for j = 1:length(d)
    expt(j) = sscanf(d(j).filename,'Expt%d');
end

G = G./sum(G);
expts = unique(expt);
if isempty(mkexpts)
    mkexpts = expts;
    savelist = 1;
else
    mkexpts = intersect(mkexpts,expts);
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
D.duration = mytoc(ts);
D.loadsum = sumload;
fprintf('%.2f of total time was load\n',D.loadsum./D.duration);
if savelist
    outname = sprintf('%s/LFPlist.mat',name);
    save(outname,'D');
end

function D = MakeLFP(e, expt, d, LFP)
    id = find(expt == e);
    clear rawlfp;
    G = LFP.kernel;
    D.sumload = 0;
    for p = 1:length(id)
        fprintf('Loading %s',d(id(p)).name);
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
        if D.skew(probe) > 5 && D.kurt(e,probe) > 1000
            goodch(probe) = 0;
        else
            goodch(probe) = 1;
        end
        fprintf('skew %.2f Kurtosis %.1f Max %.1fuV\n',D.skew(probe),D.kurt(probe),D.max(probe));
        if V.highpass == 0
            cprintf('red','%s Looks like low F were filtered out\n',d(id(p)).name);
        end
        lfp = conv(V.V,G,'same');
        rawlfp(:,probe) = downsample(lfp,LFP.decimate); %down to 500Hz;
    end
    intscale(3) = max(max(abs(rawlfp(:,find(goodch)))));
    intscale(2) = 32000;
    intscale(1) = intscale(2) ./ intscale(3);
    LPF.badchannels = find(goodch ==0);
    LFP.rawlfp = int16(round(rawlfp.*intscale(1)));
    LFP.intscale = intscale;
    LFP.t = downsample(V.t,LFP.decimate);
    LFP.samper = V.samper.*LFP.decimate; % = in sec, not silly units
    LFP.version = version;
    D.trange = minmax(LFP.t);
    outname = sprintf('%s/Expt%d.lfp.mat',LFP.name,e);
    fprintf('Saving %d samples (max %d), to %s\n',length(rawlfp),max(LFP.rawlfp(:)),outname);
    save(outname,'LFP','-v7.3');

