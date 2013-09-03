function D = FullV2LFP(name,varargin)
mkexpts = [];
G = [];
savelist = 0;
V = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'expts',5)
        j = j+1;
        mkexpts = varargin{j};
    elseif strncmpi(varargin{j},'kernel',5)
        j = j+1;
        G = varargin{j};
    end
    j = j+1;
end

d = mydir([name '/Expt*FullV.mat']);
for j = 1:length(d)
    expt(j) = sscanf(d(j).filename,'Expt%d');
end

LFP.sd = 24; %24 at 30KHz = SD of 0.8ms, = SD of 200Hz in Freq
if isempty(G)
    G = Gauss(LFP.sd,-100:100);
end
G = G./sum(G);
expts = unique(expt);
if isempty(mkexpts)
    mkexpts = expts;
    savelist = 1;
else
    mkexpts = intersect(mkexpts,expts);
end

LFP.decimate = 60;
LFP.kernel = G;
ts = now;
sumload = 0;
for j = 1:length(mkexpts)
    e = mkexpts(j);
    id = find(expt == mkexpts(j));
    clear rawlfp;
    for p = 1:length(id)
        fprintf('Loading %s',d(id(p)).name);
        V = LoadFullV(d(id(p)).name,'nohighpass','addtime');
        sumload = sumload + V.loadtime;
        probe = Name2Probe(d(id(p)).name);
        if probe ~= V.chspk
            cprintf('red','Probe name mismatch\n');
        end
        if ~isfield(V,'skew')
            D.skew(e,probe) = skewness(V.V);
            D.kurt(e,probe) = kurtosis(V.V);
        else
            D.skew(e,probe) = V.skew;
            D.kurt(e,probe) = V.kurtosis;
        end
        D.max(e,probe) = max(V.V).*1000000;
        if D.skew(e,probe) > 5 && D.kurt(e,probe) > 1000
            goodch(probe) = 0;
        else
            goodch(probe) = 1;
        end
        fprintf('skew %.2f Kurtosis %.1f Max %.1fuV\n',D.skew(e,probe),D.kurt(e,probe),D.max(e,probe));
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
    LFP.samper = V.samper.*LFP.decimate;
    D.trange(e,:) = minmax(LFP.t);
    outname = sprintf('%s/Expt%d.lfp.mat',name,mkexpts(j));
    fprintf('Saving %d samples (max %d), to %s\n',length(rawlfp),max(LFP.rawlfp(:)),outname);
    save(outname,'LFP','-v7.3');
end
D.duration = mytoc(ts);
D.loadsum = sumload;
fprintf('%.2f of total time was load\n',D.loadsum./D.duration);
if savelist
    outname = sprintf('%s/LFPlist.mat',name);
    save(outname,'D');
end
    
