function [result, sptrig] = PlotAllProbes(name, varargin)
% PlotAllProbes(name, varargin) reads in Expt files for each probe, and
% combines them. name shoudl be the file name with the probe number removed
%
% e.g. 
% PlotAllProbes('lemM025.c1.image.DCORRC.mat')
%
% if name is a directory, it does this for each file *.p1c1.* in the directory.
% if name is a stucture (previosly returned), the results are just plotted
%
% ..,'sptrig'      includes calculation of spike-triggered LFPs (takes longer).

% TrackBlank needs to track the frameresp reversal. 
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
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'movie',4)
        plottype = 1;
    elseif strncmpi(varargin{j},'blank',4)
        plottype = BLANKRESP;
    elseif strncmpi(varargin{j},'+blank',4)
        showblank = 1;
    elseif strncmpi(varargin{j},'default',5)
        plottype = 0;
    elseif strncmpi(varargin{j},'csdsum',6)
        plottype = CSDSUM;
    elseif strncmpi(varargin{j},'csdplot',3)
        lineplot = 3;
        plotargs = {plotargs{:} varargin{j}};
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            csdsk = varargin{j};
        end
            plotargs = {plotargs{:} csdsk};            
    elseif strncmpi(varargin{j},'dvdt',4)
        plotargs = {plotargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'figs',4)
        j = j+1;
        figures = varargin{j};
    elseif strncmpi(varargin{j},'frame',4)
        plottype = FRAMERESP;
    elseif strncmpi(varargin{j},'freqs',4)
        j = j+1;
        freqrange = varargin{j};
    elseif strncmpi(varargin{j},'findinv',6)
        plottype = FINDINVERSION;
    elseif strncmpi(varargin{j},'interpolate',5)
        interpolate = 1;
    elseif strncmpi(varargin{j},'LFPdiff',5)
        plottype = PLOTLFPDIFF;
    elseif strncmpi(varargin{j},'LFPeigv',7)
        plottype = PLOTLFPEIGVECTORS;
    elseif strncmpi(varargin{j},'LFPeig',5)
        plottype = PLOTLFPEIG;
%        plottype = RESPXBLANK;
    elseif strncmpi(varargin{j},'rexpXblank',5)
        plottype = RESPXBLANK;
    elseif strncmpi(varargin{j},'LFPcmp',5)
        plottype = 0;
        yvals = [1 2];
        zvals = [3 4];
    elseif strncmpi(varargin{j},'LFPpwrdiff',9)
        plottype = LFPPWRDIFF;
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            freqs = varargin{j};
        end
    elseif strncmpi(varargin{j},'LFPsigpwr',9) | strncmpi(varargin{j},'LPFsigdiff',9)
        plottype = LFPSIGPWR;
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            freqs = varargin{j};
        end
    elseif strncmpi(varargin{j},'LFPpwr',5)
        plottype = LFPFREQ;
    elseif strncmpi(varargin{j},'LFPband',5)
        plottype = LFPBANDPWR;
    elseif strncmpi(varargin{j},'LFPtrial',5)
        plottype = LFPTRIAL;
    elseif strncmpi(varargin{j},'Lineplot',5)
        lineplot = 1;
    elseif strncmpi(varargin{j},'muvar',5)
        plottype = MUVAR;
    elseif strncmpi(varargin{j},'scatterplot',5)
        lineplot = 2;
    elseif strncmpi(varargin{j},'Lines',5)
        plottype = LINEPLOT;
    elseif strncmpi(varargin{j},'latency',4)
        plottype = STIMLATENCY;
    elseif strncmpi(varargin{j},'latency',4)
        j = j+1;
        latency = varargin{j};
    elseif strncmpi(varargin{j},'monoc',5)
        plottype = PLOTMONOC;
    elseif strncmpi(varargin{j},'netspk',5)
        plottype = 2;
    elseif strncmpi(varargin{j},'ocularity',4)
        plottype = OCULARITY;
    elseif strncmpi(varargin{j},'onestimMU',9)
        plottype = PLOTONESTIMMU;
    elseif strncmpi(varargin{j},'onestim',4)
        plottype = PLOTONESTIM;
    elseif strncmpi(varargin{j},'probes',4)
        j = j+1;
        probes = varargin{j};
    elseif strncmpi(varargin{j},'pause',4)
        j = j+1;
        delay = varargin{j};
    elseif strncmpi(varargin{j},'runlist',4)
        fid = fopen(name,'r');
        names = textscan(fid,'%s');
        for k = 1:length(names{1})
            [dir, file] = fileparts(names{1}{k});
            PlotAllProbes(dir,varargin{1:j-1});
        end
    elseif strncmp(varargin{j},'subplot',4)
        j = j+1;
        subplots = varargin{j};
        nr = subplots(1);
        nc = subplots(2);
        pn = subplots(3:end);
    elseif strncmpi(varargin{j},'save',4)
        saveresult = 1;
        argon = {argon{:} varargin{j}};
    elseif strncmpi(varargin{j},'scales',4)
        j = j+1;
        newscale = 1;
        scales(1:length(varargin{j})) = varargin{j};
    elseif strncmpi(varargin{j},'setsign',4)
        j = j+1;
        setlfpsign = varargin{j};
    elseif strncmpi(varargin{j},'sptrig',4)
        calcsptrig = 1;
        argon = {argon{:} varargin{j}};
    elseif strncmpi(varargin{j},'sumy',4)
        sumy = 1;
    elseif strncmpi(varargin{j},'diffy',5)
        sumy = 2;
    elseif strncmpi(varargin{j},'timerange',6)
        j = j+1;
        timerange = varargin{j};
    elseif strncmpi(varargin{j},'stimvar',6)
        plottype = STIMVAR; %stim induced variance over time
    elseif strncmpi(varargin{j},'varblank',6)
        plottype = BLANKVAR; %stim induced variance over time
    elseif strncmpi(varargin{j},'writerf',7)
        WriteRFLoc(name);
        return;
    elseif strncmpi(varargin{j},'xtprobesmu',10)
        plottype = XTPROBELFPMU;
    elseif strncmpi(varargin{j},'xtprobes',4)
        plottype = XTPROBELFP;
    elseif strncmpi(varargin{j},'xval',4)
        j = j+1;
        xvals = varargin{j};
    elseif strncmpi(varargin{j},'yval',4)
        j = j+1;
        yvals = varargin{j};
           
    elseif strncmpi(varargin{j},'zval',4)
        j = j+1;
        zvals = varargin{j};
    elseif isdir(varargin{j})
        result = ProcessDir(varargin{j},saveresult, argon{:});
    end
    j = j+1;
end

if isstruct(name)
    rc = name;
% if there is an 'lfp' matrix, uset that. Otherwise, if this is an rc file
% (output from PlotRevCOrAny, build the matrix
if sum(yvals) == 0
    if isfield(rc,'tresps')
        yvals = 1:size(rc.tresps,3);
    end
end
    if ~isfield(name,'lfp') & isfield(name,'sdfs') & isfield(name.sdfs,'lfp')
        [rc.lfp, rc.mlfp, rc.lfpblank,a ] = CombineRCLFP(rc,scales);
        rc.lfpblankn = a.blankn(1);
        rc.lfptimes = rc.sdfs.lfptimes;
        if ~isfield(rc,'x')
            rc.x = a.x';
        end
        if ~isfield(rc,'y')
            rc.y = a.y;
        end
        if isfield(a,'lm')
            rc.lfplm = a.lm;
        end
        if isfield(a,'rm')
            rc.lfprm = a.rm;
        end
        if isfield(a,'uc')
            rc.lfpuc = a.uc;
        end
        isrc = 1;
        if ~isfield(rc.Stimvals,'xtype')
            rc.Stimvals.xtype = rc.Stimvals.et;
            rc.Stimvals.ytype = rc.Stimvals.e2;
        end
    end
    rc.plot = 0;

    if isfield(rc,'pk');
        ofig = gcf;
        if isfield(rc,'lfp')
            rc.tid = find(rc.lfptimes > timerange(1) & rc.lfptimes < timerange(2));
            [rc.eigs, rc.eigresp, rc.eigblank,x] = PlotLFPEig(rc,NaN, sumy, setlfpsign);
        [a,b] = NetSpkPref(rc.eigresp,rc.x(:,1),'zerobase');
        lfppref = angle(b.vsum);
        rc.lfpprefs = x.prefs;
        rc.lfpprefdir = lfppref;
        end

        GetFigure('Psych');
        id = find(rc.pk.xv > -1000);
        if isfield(rc,'pp')
            subplot(3,1,2);
            hold off;
            nmin = max([rc.pp.n])./5;
            fitpsf(rc.pp,'nmin',nmin, 'showfit','shown');
            subplot(3,1,1);
        else
            subplot(1,1,1);
        end
        hold off;
        plot(rc.pk.xv(id),rc.pk.kernel(id));
        rc.pk.xv = rc.pk.xv(id);
        rc.pk.kernel = rc.pk.kernel(id);
        rc.pk.n = rc.pk.n(id,:);
        sigs = unique(rc.or);
        sigs= mod(sigs,180);
        for j = 1:length(sigs)
            hold on;
            plot([sigs(j) sigs(j)],[0 max(rc.pk.kernel)],'r:');
        end
        title(sprintf('mean %.1f frames',mean(rc.pk.n(:))));
        [rcprefs, x] = GetNetSpkPref(rc,[timerange(3) timerange(4)]);
        prefo = rc.x(rcprefs,1);
        rc.cv = abs(x.vsum);
        cproj(:,1) = cos((prefo-mean([rc.cp.upstim])) * pi/180);
        cproj(:,2) = cos((prefo-mean([rc.cp.dnstim])) * pi/180);
        [a, prefchoice] = max(abs(cproj),[],2);
        pid = find(prefchoice == 1);
        nid = find(prefchoice == 2);
        if isfield(rc,'cp')
            for j = 1:length(rc.cp)
                if isempty(rc.cp(j).cp)
                    rc.cp(j).cp = NaN;
                end
            end
            subplot(3,1,3);
            hold off;
            plot([rc.cp.gcp],'ro');
            hold on;
            plot(pid,[rc.cp(pid).cp],'o','MarkerFaceColor','b');
            [rc.cp.pcp] = deal(rc.cp.cp);
%pcp is calculated accoring to pref ori from netspk tuning
            for j = 1:length(nid)
                rc.cp(nid(j)).pcp = 1-rc.cp(nid(j)).cp;
            end
            plot(nid,[rc.cp(nid).cp],'o');
            cp = mean([rc.cp.cp]);
            plot([1 24],[cp cp]);
            cp = mean([rc.cp.gcp]);
            plot([1 24],[cp cp],'r');
            plot([1 24],[0.5 0.5],'k:');
            title(sprintf('CP for Pref = %.1f',mean([rc.cp.upstim])));
        end
        if isfield(rc,'spkres')
        for j = 1:length(rc.spkres)
            if isfield(rc.spkres{j},'cp')
                [netspk, details] = PlotRC(rc.spkres{j},'netspk','zcheck');
                [maxr, pref] = max(mean(netspk,2));
                [a,x] = NetSpkPref(netspk,details.xv);
                rc.cellcv(j) = abs(x.vsum);
                cellproj(j,1) = cos((rc.x(pref,1)-mean([rc.cp.upstim])) * pi/180);
                cellproj(j,2) = cos((rc.x(pref,1)-mean([rc.cp.dnstim])) * pi/180);
                h = plot(mean(rc.spkres{j}.probes),rc.spkres{j}.cp.cp,'ko');
%cellcp is the CP according to pref by count for signal trials.
                if rc.spkres{j}.cp.prefbycount == rc.spkres{j}.cp.upstim
                    rc.cellcp(j) = rc.spkres{j}.cp.cp;
                else
                    rc.cellcp(j) = 1-rc.spkres{j}.cp.cp;
                end

                rc.cellprobe(j) = mean(rc.spkres{j}.probes);
%CP according to cells own pref, from subspace map
                if diff(abs(cellproj(j,:))) < 0 %prefers upstim
                    set(h,'markerfacecolor','k');
                    rc.cellpcp(j) = rc.spkres{j}.cp.cp;
                else
                    rc.cellpcp(j) = 1-rc.spkres{j}.cp.cp;
                end
%CP according to LFP subspace pref of this probe
                if rc.cellprobe(j) > 0
                    cellproj(j,1) = cos(lfppref(round(rc.cellprobe(j)))-mean([rc.cp.upstim] * pi/180));
                    cellproj(j,2) = cos(lfppref(round(rc.cellprobe(j)))-mean([rc.cp.dnstim] * pi/180));
                    if diff(abs(cellproj(j,:))) < 0 %LFP prefers upstim
                        rc.celllcp(j) = rc.spkres{j}.cp.cp;
                    else
                        rc.celllcp(j) = 1-rc.spkres{j}.cp.cp;
                    end
                else
                    rc.celllcp(j) = NaN;
                end
                
            end
        end
        end
        if isfield(rc,'lfpblank')
            [a, bmax] = max(max(rc.lfpblank));
            ylim = get(gca,'Ylim');
            plot([bmax bmax],ylim,'k');
            rc.blankmax = bmax;
        end
        set(0,'CurrentFigure',ofig);
    end
        
        
    if isfield(rc,'tresps')
        if ndims(rc.tresps) == 3
            rc.spkresps = rc.tresps
        elseif ndims(rc.tresps) == 4
            if max(yvals) > size(rc.tresps,3)
                yvals = yvals(find(yvals <= size(rc.tresps,3)));
            end
            rc.spkresps = squeeze(mean(rc.tresps(:,:,yvals,:),3));
            if isfield(rc,'blankresp') & showblank;
                j = size(rc.spkresps,2)+1;
                rc.spkresps(:,j,:) = rc.blankresp;
            end
            nsum = sum(rc.spkn(:));
            rc.iscell(1:size(rc.tresps,4)) = 0;
            if isfield(rc,'spkres')
                
                for j = 1:length(rc.spkres)
                    if length(rc.spkres{j});
                        rate = mean(rc.spkres{j}.y(:));
                        ssum = sum(rc.spkres{j}.sdfs.n(:));
                        p = round(mean(rc.spkres{j}.probes));
                        tid = find(rc.spkres{j}.times >= rc.times(1) & rc.spkres{j}.times <= rc.times(end)); 
                        if ssum > nsum/2
                            rc.iscell(p) = 1;
                            for n = 1:size(rc.tresps,3);
                            for m = 1:size(rc.tresps,2);
                                if ~isempty(rc.spkres{j}.sdfs.s{m,n})
                                rc.tresps(:,m,n,p) = rc.spkres{j}.sdfs.s{m,n}(tid)./rate;
                                end
                            end
                            end
                        end
                    end
                end
            end
        end
        isrc = 1;
    elseif isfield(name,'sdfs') & isfield(name.sdfs,'lfp')
        isrc = 1;
    end
    if isfield(rc,'isrc')
        isrc = rc.isrc;
    end
    rc.plottype = plottype;
    
    if isfield(rc,'lfp')
        tid = find(rc.lfptimes > timerange(1) & rc.lfptimes < timerange(2));
        rc.tid = tid;
        if ndims(rc.lfp) == 3
            rc.lfpresp = rc.tresps
        elseif ndims(rc.lfp) == 4
            if size(rc.lfp,4) < max(probes)
                probes = [1:size(rc.lfp,4)];
            end
            if sum(yvals) == 0
                yvals = 1:size(rc.lfp,3);
            end
            if max(yvals) > size(rc.lfp,3)
                yvals = yvals(find(yvals <= size(rc.lfp,3)));
            end
            if max(zvals) > size(rc.lfp,3)
                zvals = zvals(find(zvals <= size(rc.lfp,3)));
            end
            if max(xvals) > size(rc.lfp,2)
                xvals = xvals(find(xvals <= size(rc.lfp,2)));
            end

            if ~isfield(rc,'lfpn') & isfield(rc,'sdfs')
                rc.lfpn = rc.sdfs.n;
            end
            nrep = rc.lfpn(:,yvals);
            nrep = mean(sum(nrep,2));
            rc.probes = probes;
            rc.csdk = csdsk;
            if ismember(plottype,[PLOTLFPEIG RESPXBLANK PLOTLFPEIGVECTORS])
                [rc.eigs, rc.eigresp, rc.eigblank,x] = PlotLFPEig(rc,lineplot, sumy, setlfpsign);
                [a,b] = NetSpkPref(rc.eigresp, rc.x(:,1),'zerobase');
                rc.lfpprefs = x.prefs;
                rc.lfpprefdir = angle(b.vsum);
                if plottype == RESPXBLANK
                    subplot(1,1,1);
                    imagesc(x.blankresp');
                    hold on;
                elseif plottype == PLOTLFPEIGVECTORS
                    subplot(3,1,1);
                    hold off;
                    imagesc(rc.lfptimes(tid)./10,probes,rc.eigs(:,probes)');
                    subplot(3,1,2); hold off;
                    imagesc(rc.lfptimes(tid)./10,probes,x.beigs(:,probes)');
                    subplot(3,1,3); hold off;
                    imagesc(rc.lfptimes(tid)./10,probes,x.ceigs(:,probes)');
                else
                if lineplot == 3
                subplot(3,1,2);
                plot(rc.lfpprefdir.*180,probes(2:end-1),'b:');
                elseif lineplot == 0
                    if strcmp(rc.Stimvals.et,'or')
                        subplot(3,1,2);
                        plot(rc.lfpprefdir(probes).*180,probes,'b:');
                    end
                end
                end
            elseif plottype == STIMVAR %var across stim, wrt t, probe
                if isfield(rc,'tresps')
                    subplot(2,1,1);
                    if sumy
                        resp = squeeze(sum(rc.tresps,3)); %% pool all Y to get eigenvectors
                    else
                        resp = reshape(rc.tresps,[size(rc.tresps,1) size(rc.tresps,2)*size(rc.tresps,3) size(rc.tresps,4)]);
                    end
                    rc.spkim.z = squeeze(var(resp,[],2));
                    hold off;
                    [sdfl, details] = SDFlatency(rc,sumy);
                    rc.latencies = sdfl;
                    rc.latsig = details.maxsdr;
                    
                    if lineplot == 1
                    h = plot(rc.times./10,rc.spkim.z);
                    legend(h,num2str(rc.probes'));
                    else
                    imagesc(rc.times./10,probes,rc.spkim.z');
                    hold on; 
                    plot(sdfl(probes)./10,probes,'w:');
                    id = find(rc.iscell);
                    for j = 1:length(id)
                        plot(rc.times([1 10])./10,[id(j) id(j)],'r-');
                    end
                    end
                    
                    subplot(2,1,2);
                end
                [ll, details] = LFPLatency(rc.lfp,rc.lfptimes);
                rc.lfplatencies = ll;
                rc.lfplatsig = details.maxsdr;
                v = squeeze(mean(var(rc.lfp,[],2),3));
                if isfield(rc,'times')
                tid = find(rc.lfptimes > timerange(1) & rc.lfptimes < max(rc.times) & rc.lfptimes < timerange(2));
                else
                    tid = find(rc.lfptimes > timerange(1) & rc.lfptimes < timerange(2));
                end
                imagesc(rc.lfptimes(tid)./10,probes,v(tid,:)');
                [peakv, pmax] = max(max(v));
                hold on
                plot([1 size(v,1)],[pmax pmax],'r');
                plot(ll(probes,4)./10,probes,'w:');
                if lineplot == 2
                    GetFigure('Scatter');
                    plot(ll(probes,4)./10,sdfl(probes)./10,'o-');
                    refline(1);
                end
            elseif plottype == BLANKVAR %var across stim, wrt t, probe
                preid = find(rc.lfptimes > 100 & rc.lfptimes < 500);
                if size(rc.lfp,2) > 1
                subplot(nr,nc,pn(1));
                [a,b] = LFPLatency(rc.lfp,rc.lfptimes);
                rc.lfplatencies  = a;
                hold off;
                rc.lfpim.z = b.var(tid,:)';
                rc.lfpim.x = b.times(tid)./10;
                rc.lfpim.y = probes;
                imagesc(b.times(tid)./10,probes,b.var(tid,:)');
                hold on; plot(probes,a(:,3)./10,'r:');
                [t, probe] = min(a(:,3));
                [mv, mvt] = max(b.var(tid,:));
                [a, mp] = max(mv);
                rc.maxvarprobe = mp;
                rc.minlatprobe = probe;
                title(sprintf('Min Lat %.1f Probe %d, Max Var %d',t./10,probe,mp));
                end
                if isfield(rc,'mlfp')
                    x = GetFrameReversal(rc.mlfp, rc.lfptimes, probes);
%                    [rc.framereverseprobe, rc.frameminprobe,x] = GetFrameReverse(rc, probes);
                    rc.framemaxprobe = x.maxprobe;
                    rc.framereverseprobe = x.revp;
                    rc.frameminprobe = x.minprobe;
                end
                if isfield(rc,'lfpblank')
                subplot(nr,nc,pn(2));
                hold off;
                prev = rc.lfpblank(preid,:);
                rc.blanksdr = max((rc.lfpblank(tid,:)))./std(prev);
                rc.lfpim.zb = rc.lfpblank(tid,:)';
                imagesc(rc.lfptimes(tid)./10,probes, rc.lfpblank(tid,:)');
                hold on;
                [v, t] = max(abs(rc.lfpblank));
                [a,probe] = max(v);
                plot(probes,rc.lfptimes(t(probes))./10,'r:');
                z = conv2(rc.lfpblank,csdsk,'same');
                [a,b] = LFPLatency(z,rc.lfptimes);
                plot(probes,a(probes,3)./10,'w:');
                title(sprintf('Max Resp at %.1f (%.1fV) Probe %d',rc.lfptimes(t(probe))./10,v(probe),probe));
                end
            elseif plottype == STIMLATENCY %var across stim, wrt t, probe
                [a,b] = LFPLatency(rc.lfp,rc.lfptimes);
                hold off;
                imagesc(probes,b.times./10,b.var);
                hold on; plot([1:24],a(:,3)./10,'r:');
            elseif ismember(plottype, [LFPTRIAL])
                lfp = zeros(size(rc.lfp,1),size(rc.lfp,4));
                for k = 1:size(rc.lfp,3);
                for j = 1:size(rc.lfp,2);
                    lfp = lfp + squeeze(rc.lfp(:,j,k,:)) * rc.lfpn(j,k);
                end
                end
                hold off;
                imagesc(lfp);
            elseif ismember(plottype, [LFPPWRDIFF LFPSIGPWR]) %power spec for pref vs null choices
                % sigpwr is power diff by stim summed over freqs, for each probe
                % upstim (lower value) is first
                if plottype == LFPSIGPWR
                    bychoice = 0; %show pwr diff by signal
                else
                    bychoice = 1; %show pwr diff by signal
                end
                freqs = find(rc.lfpfrq > freqrange(1) & rc.lfpfrq < freqrange(2));
                if isfield(rc,'lfpblank')
                [rc.eigs, rc.eigresp, rc.eigblank,x] = PlotLFPEig(rc,NaN, sumy, setlfpsign);
                [a, b] = NetSpkPref(rc.eigresp, rc.x(:,1),'zerobase');
                x.prefdir = angle(b.vsum) .* 180/pi;
                cliprange = [-65 115];
                prefo = rc.x(x.prefs,1);
                id = find(prefo > cliprange(2));
                while length(id) > 0
                    prefo(id) = prefo(id) - 180;
                    id = find(prefo > cliprange(2));
                end
                id = find(prefo < cliprange(1));
                while length(id) > 0
                    prefo(id) = prefo(id) + 180;
                    id = find(prefo < cliprange(1));
                end
                id = find(x.prefdir > cliprange(2));
                while length(id) > 0
                    x.prefdir(id) = x.prefdir(id) - 180;
                    id = find(x.prefdir > cliprange(2));
                end
                id = find(x.prefdir < cliprange(1));
                while length(id) > 0
                    x.prefdir(id) = x.prefdir(id) + 180;
                    id = find(x.prefdir < cliprange(1));
                end

                [a, bmax] = max(max(rc.lfpblank));
                end
                if isfield(rc,'lfpsigpwr')
                    upstim = mean([rc.cp.upstim]);
                    theta = (rc.x(x.prefs(probes),1) -upstim) * pi/180;
                    reverse = find(abs(sin(theta)) > sin(pi/4));
                    for j = 1:length(freqbands)-1
                        ifq = find(rc.lfpfrq >= freqbands(j) & rc.lfpfrq < freqbands(j+1));
                        bandpwr(:,j) =  squeeze(sum(diff(abs(rc.lfpsigpwr(:,ifq,:))),2));
                    end
                        

                    sigpwr =  squeeze(sum(diff(abs(rc.lfpsigpwr(:,freqs,:))),2));
                    pwrdiff = squeeze(diff(abs(rc.lfpsigpwr)));
                    pwrsum = squeeze(sum(abs(rc.lfpsigpwr)))./2;
                    if bychoice
                        pwrsum = squeeze(sum(abs(rc.lfpwr),1))./2;
                        pwrdiff = squeeze(diff(abs(rc.lfpwr)));
                    end
                    upstim = mean([rc.cp.upstim]); 
                    dnstim = mean([rc.cp.dnstim]); 

                    prefstim = ones(size(sigpwr)) .* upstim; 
                    nullstim = ones(size(sigpwr)) .* dnstim; 
                    id = find(sigpwr) > 0;
                    if length(id)
                        prefstim(id) = dnstim;
                        nullstim(id) = upstim;
                    end    
                    theta = (rc.x(x.prefs(probes),1) -prefstim) * pi/180;
                    reverse = find(abs(sin(theta)) > sin(pi/4));
                    rc.prefs = x.prefs;
                    rc.prefdir = x.prefdir;
                    rc.pwrpref = prefstim;
                    rc.bandsigpwr = bandpwr;
                else
                    prefstim = NaN;
                    reverse = [];
                    if ndims(rc.lfpwr) > 2
                        pwrsum = squeeze(sum(rc.lfpwr,1));
                        pwrdiff = squeeze(diff(rc.lfpwr));
                    else
                        pwrsum = squeeze(rc.lfpwr,1);
                        pwrdiff = zeros(size(pwrsum));
                    end
                end
                
                
                if isfield(rc,'spkres')
                for j = 1:length(rc.spkres)
                    if isfield(rc.spkres{j},'cp')
                    probe = round(mean(rc.spkres{j}.probes));
                        if rc.pwrpref(probe) == dnstim
                            rc.cellpwrcp(j) = 1 -rc.spkres{j}.cp.cp;
                        else
                            rc.cellpwrcp(j) = rc.spkres{j}.cp.cp;
                        end
                    end
                end
                end
                        
                        
                normalizediff = 0;
                pwrmax = max(max(pwrsum(freqs,:)));
                subplot(4,1,1);
                if normalizediff
                    diffresp = pwrdiff(freqs,:)./pwrsum(freqs,:);
                else
                    diffresp = pwrdiff(freqs,:);
                end
                if lineplot == 1
                    
                plot(rc.lfpfrq(freqs),diffresp);
                else
                imagesc(probes, rc.lfpfrq(freqs),diffresp);
                colorbar;
                end
                rc.pwrdiff = pwrdiff./pwrsum;
                if ~isnan([rc.cp.upstim])
                title(sprintf('%.2f - %.2f',rc.cp(1).dnstim, rc.cp(1).upstim));
                end
%                caxis([0 pwrmax]);
                subplot(4,1,2);
                if lineplot == 1
                plot(rc.lfpfrq(freqs),abs(pwrsum(freqs,:))');
                else
                imagesc(probes, rc.lfpfrq(freqs),abs(pwrsum(freqs,:)));
                caxis([0 pwrmax]);
                colorbar;
                end
                if ~isnan(prefstim)
                    title(sprintf('%.0f,',mean(rc.x(x.prefs,1))));
                end
                subplot(4,1,3);
                pwrpref(find(sigpwr < 0)) = mean([rc.cp.upstim]);
                pwrpref(find(sigpwr > 0)) = mean([rc.cp.dnstim]);
                hold off;
                plot(probes,prefo,'o-');
                hold on;
                plot(probes,pwrpref(probes),'ko-');
                plot(probes,x.prefdir(probes),'mo-');
                plot([bmax bmax],get(gca,'ylim'),'k:');
                rc.bmaxprobe = bmax;
                subplot(4,1,4);
                rc.freqbands = freqbands;
                for j = 1:length(freqbands)-1
                    id = find(rc.lfpfrq >= freqbands(j) & rc.lfpfrq <= freqbands(j+1));
                    rc.pwrbands(j,:) = mean(pwrdiff(id,:))./mean(pwrsum(id,:));
                    labels{j} = sprintf('%.1f-%.1f',freqbands(j),freqbands(j+1));
                end
                if length(reverse) %
                plot(probes(reverse),rc.pwrbands(:,reverse)','o');
                hold on;
                end
                plot(probes,rc.pwrbands(:,probes)');
                legend(labels);
            elseif isrc
                if plottype == PLOTLFPDIFF & rc.ispsych
                    yvals = [1 2];
                    zvals = [3 4];
                    rc.plot = PLOTLFPDIFF;
                end
                rc.lfpresp = squeeze(mean(rc.lfp(:,:,yvals,:),3));

                if isfield(rc,'lfpblank') & length(rc.lfpblank) & showblank
                    j = size(rc.lfpresp,2)+1;
                    rc.lfpresp(:,j,:) = rc.lfpblank;
                end
                if ~isempty(zvals)
                    rc.lfprespa = squeeze(mean(rc.lfp(:,:,zvals,:),3));
                end

            else
                rc.lfpresp = squeeze(mean(rc.lfp(:,:,yvals,:),3));
            end
            rc.titlestr{1} = sprintf('N%.1f',nrep);
            if isfield(rc,'cp')
            rc.titlestr{2} = sprintf('Choice %.1f - %.1f',rc.cp(1).dnstim,rc.cp(1).upstim);
            else
                rc.titlestr{2} = '';
            end
            rc.nrep = nrep;
        end
        if newscale
            for j = 1:size(rc.lfpresp,3)
                rc.lfpresp(:,:,j) = rc.lfpresp(:,:,j) .* scales(j);
            end
        end
    end

    if plottype == FINDINVERSION
        FindSdfCorr(rc);
        return;
    end
    
    mtimes = timerange(1):10:timerange(2);
    if plottype == 0
        if isrc
            plottype = 1;
        else
            plottype = 3; %image
        end
    end
    if ~isfield(rc,'x') | length(rc.x) == 1
        plottype = 3;
    end
    
    if isfield(rc,'means') && length(probes) > length(rc.means)
    rc.probes = probes(1:length(rc.means));
    probes = rc.probes;
    end
    if isfield(rc,'lfp')
        tid = find(rc.lfptimes > timerange(1) & rc.lfptimes < timerange(2));
        if size(rc.lfp,4) < max(probes)
            probes = [1:size(rc.lfp,4)];
        end
    end
    if isfield(rc,'tresps')
        stid = find(rc.times > timerange(1) & rc.times < timerange(2));
    end
    rc.yvals = yvals;
    
    if isfield(rc,'Stimvals') & ~isfield(rc.Stimvals,'xtype')
        rc.Stimvals.xtype = rc.Stimvals.et;
        rc.Stimvals.ytype = rc.Stimvals.e2;
    end
    rc.lineplot = lineplot;
    rc.csdsk = csdsk;
    if plottype == 1
        if length(probes) > 1
            if size(rc.x,1) == 1
                rc.x = rc.x';
            end
            result = PlotMovie(rc,delay,mtimes);
        else
            PlotTimeCourse(rc,probes, lineplot, sumy);
        end
    elseif plottype == STIMVAR & isfield(rc,'tresps') & ~isfield(rc,'lfp')
        if showblank
            np = 3;
            subplot(np,1,3);
            imagesc(probes,rc.lfptimes,rc.lfpblank);
        else
            np = 2;
        end
        
        subplot(np,1,1);
        v = squeeze(mean(var(rc.tresps,[],2),3));
        imagesc(probes,rc.times./10,v);
        subplot(np,1,2);
        v = squeeze(mean(var(rc.lfp,[],2),3));
        imagesc(probes,rc.lfptimes./10,v);
        
    elseif plottype == OCULARITY %netspike plot for spikes
        l = squeeze(rc.tresps(:,1,1,:));
        r = squeeze(rc.tresps(:,3,1,:));
        z = rc.blankresp;
        odi = (l-r)/((l-z)+(r-z));
        imagesc(odi);
    elseif plottype == 2 %netspike plot for spikes
        [prefs, x] = GetNetSpkPref(rc, [timerange(3) timerange(4)]);
        netspk = x.netspk;
        if isfield(x,'blankspk')
            netspk(size(netspk,1)+1,:,:) = x.blankspk;
        end
        if lineplot == 2
            ns = squeeze(sum(netspk,1));
            ya = round(size(ns,2)/2);
            yb = size(ns,2);
            for j = 1:length(probes)
                subplot(6,4,j);
                scatter(ns(:,ya,j),ns(:,yb,j));
                [a,b] = fit_bothsubj2error(ns(:,ya,j),ns(:,yb,j));
                title(sprintf('Slope %.2f',b));
                refline(b,a);
            end
        elseif sumy
            subplot(1,1,1);
            netspk = squeeze(mean(netspk,2));
            if lineplot
                subplot(2,1,1);
                plot(rc.x(:,1),netspk);
                subplot(2,1,2);
                plotyy(probes,abs(x.vsum(probes)),probes,angle(x.vsum(probes)).*90/pi);
            else
            imagesc(rc.x(:,1),probes,netspk');
            end
        else
            ns = netspk;
            cr = [min(ns(:)) max(ns(:))];
            ny = size(netspk,2);
            for j = 1:ny
                subplot(ny,1,j);
                if lineplot 
                plot(rc.x(:,1),squeeze(netspk(:,j,:)));
                set(gca,'ylim',cr);
                else
                imagesc(rc.x(:,1),probes,squeeze(netspk(:,j,:))');
                caxis(cr);
                end
            end
        end
        if isfield(rc,'spkres')
            colors = mycolors;
            GetFigure('Cells');
            subplot(1,1,1);
            hold off;
            nc = 1;
            for j = 1:length(rc.spkres)
                if isfield(rc.spkres{j},'sdfs')
                    netspk = PlotRC(rc.spkres{j},'netspk');
                    plot(rc.x(:,1),sum(netspk,2),'color',colors{j});
                    if isfield(rc.spkres{j},'cp')
                    labels{nc} = sprintf('%.1f (%.0f)',mean(rc.spkres{j}.probes),rc.spkres{j}.cp.prefbycount);
                    else
                    labels{nc} = sprintf('%.1f ',mean(rc.spkres{j}.probes));
                    end
                    nc = nc+1;
                    hold on;
                end
            end
            if nc > 1
            legend(labels);
            end
        end
%        colorbar;
    elseif plottype == 3
        if isfield(rc,'resps') & size(rc.resps,2) > 1
            if isfield(rc,'lfp')
                tid = find(rc.lfptimes > latency);
                subplot(2,1,2);
                hold off;
                if length(probes) == 1
                    if lineplot
                        plot(rc.lfptimes./10, rc.lfpresp(:,:,probes));
                    else
                        imagesc(rc.x(:,1), rc.lfptimes(tid)./10, rc.lfpresp(tid,:,probes));
                    end
                else
                    imagesc(rc.x(:,1), probes, squeeze(var(rc.lfpresp(tid,:,:)))');
                end
                subplot(2,1,1);
            else
                subplot(1,1,1);
                hold off;
            end
            set(gca,'xscale','lin')
            if sumy
                imagesc(rc.x(:,1), probes, squeeze(sum(rc.resps(:,yvals,:),2)));
            else
                imagesc(rc.x(:,1), probes, squeeze(rc.resps(:,yvals(1),:)));
            end
        else 
            if isfield(rc,'lfp')
                if max(probes) > size(rc.lfp,4)
                    probes
                end
                if size(rc.lfp,4) < max(probes)
                    probes = [1:size(rc.lfp,4)];
                end
                tid = find(rc.lfptimes > timerange(1) & rc.lfptimes < timerange(2));
                if lineplot == 1
                    if spacing > 0
                        hold off;
                        colors = mycolors;
                        for j = 1:length(probes)
                            plot(rc.lfptimes./10,rc.lfpresp(:,j) + spacing.*j,'color',colors{1+mod(j-1,5)});
                            hold on;
                        end
                        hold off;
                    else
                        plot(rc.lfptimes./10,rc.lfpresp);
                    end
                elseif size(rc.lfp,2) == 1 & size(rc.lfp,3) == 1 | strmatch(rc.Header.exptype,'Flash') % one condition - plot like Flash
                    rc = PlotOneLFP(rc, timerange, probes, plotargs{:});
                    %                subplot(1,1,1);
                else
                    if isfield(rc,'lfpresp')
% Non rc plot. plot Variance as a function of stim, probe
                    rc.lfpim.z = squeeze(var(rc.lfpresp))';
                    rc.lfpim.x = rc.lfptimes(tid)./10;
                    rc.lfpim.y = 1:size(rc.lfpim.z,1);
                    else
                    rc.lfpim.z = squeeze(rc.lfp);
                    rc.lfpim.x = rc.lfptimes(tid)./10;
                    rc.lfpim.y = 1:size(rc.lfpim.z,1);
                    end
%                    rc.lfpim.z = rc.lfpim.z(:,tid);
                    hold off;
                    imagesc(rc.lfpim.x,rc.lfpim.y,rc.lfpim.z);
                end
            else
                plot(rc.means);
            end
        end
    elseif plottype == PLOTLFPDIFF;
        if length(probes) > 1
            result = PlotMovie(rc,delay,mtimes);
        else
            PlotTimeCourse(rc,probes,lineplot, sumy);
        end
    elseif plottype == BLANKRESP;
        if ~isfield(rc,'lfpblank')
            return;
        end
        if isfield(rc,'blankresp')
        subplot(2,1,1);
        if lineplot == 1
        plot(rc.times(stid)./10,rc.blankresp(stid,probes)');
        else
        imagesc(rc.times(stid)./10,probes,rc.blankresp(stid,probes)');
        end
        subplot(2,1,2);
        end
        hold off;
        if lineplot == 3
            imagesc(CalcCSD(rc.lfpblank(tid,probes)','smooth',csdsk));
        elseif lineplot == 1
            plot(rc.lfptimes(tid)./10,rc.lfpblank(tid,:)');
        else
        imagesc(rc.lfptimes(tid)./10,probes, rc.lfpblank(tid,:)');
        end
        rc.lfpim.x = rc.lfptimes(tid)./10;
        rc.lfpim.y = probes;
        rc.lfpim.z = rc.lfpblank(tid,:)';
        [rc.lfplatencies, a] = LFPLatency(rc.lfpblank, rc.lfptimes);
        [rc.blank.minlat(1), rc.blank.minlat(2)]  = min(rc.lfplatencies(:,3));
        rc.blank.maxprobe = a.maxprobe;
        hold on;
        if ismember(lineplot,[0 3])
        plot(rc.lfplatencies(probes,3)./10,probes,'r:');
        end
        result = rc;
    elseif plottype == FRAMERESP;
        if isfield(rc,'tresps');
        subplot(2,1,1);
        sumresp = [];
            for j = 1:size(rc.tresps,2)
                for k = 1:size(rc.tresps,3)
                    sumresp = cat(3,sumresp,squeeze(rc.tresps(:,j,k,:)).*rc.spkn(j,k));
                end
            end
        sumresp = sum(sumresp,3)./sum(rc.spkn(:));
        imagesc(sumresp');
        set(gca,'xlim',[75 125]);
        subplot(2,1,2);
        hold off;
        end
        
        if ndims(rc.mlfp) == 3
            mlfp = squeeze(mean(rc.mlfp,2));
        else
            mlfp = rc.mlfp;
        end
        tid = find(rc.lfptimes > 750);
        [rc.framereverseprobe, rc.frameminprobe, x] = GetFrameReverse(rc, probes);
        rc.framemaxprobe = x.maxprobe;
        mlfp(tid(1),probes) = x.amps;
        if lineplot == 3
            csd = CalcCSD(mlfp(probes,:),'flip','smooth',csdsk);
            imagesc(rc.lfptimes./10,probes(2:end-1),csd);
        else
            imagesc(rc.lfptimes./10,probes,mlfp');
        end
        set(gca,'xlim',[75 125]);
        hold on;
        plot(16.6 * 6 + x.phases .* 16.6/( 2 * pi),probes,'r:');
        plot([75 125],[rc.framereverseprobe rc.framereverseprobe],'r')
        
    elseif plottype == LFPBANDPWR; %Stimulus related Power in Band
        subplot(1,1,1);
        freqs = find(rc.lfpfrq > freqrange(1) & rc.lfpfrq < freqrange(2));
        hold off;

        if lineplot == 1
            if ndims(rc.lfpwr) == 2
                plot(rc.lfpfrq(freqs),rc.lfpwr(freqs,probes));
            elseif ndims(rc.lfpwr) == 3
             plot(rc.lfpfrq(freqs),abs(squeeze(sum(rc.lfpwr(:,freqs,probes))))');
            elseif ndims(rc.lfpwr) == 4
             plot(rc.lfpfrq(freqs),abs(squeeze(sum(sum(rc.lfpwr(freqs,:,yvals,probes),3),2)))');
            title(sprintf('%.1f-%.1fHz',freqrange(1),freqrange(2)));
            end
        elseif ndims(rc.lfpwr) == 2
            imagesc(rc.lfpfrq(freqs),probes,rc.lfpwr(freqs,probes)');
        elseif ndims(rc.lfpwr) == 3
            imagesc(rc.lfpfrq(freqs),probes,abs(squeeze(sum(rc.lfpwr(:,freqs,probes))))');
        elseif ndims(rc.lfpwr) == 4
            subplot(3,1,1);
            lfpwr = reshape(rc.lfpwr,[size(rc.lfpwr,1) size(rc.lfpwr,2)*size(rc.lfpwr,3)*size(rc.lfpwr,4)]);
            imagesc(rc.x(:,1),probes,squeeze(sum(sum(rc.lfpwr(freqs,:,yvals,probes),3),1))');
            title(sprintf('%.1f-%.1fHz',freqrange(1),freqrange(2)));
            [E,V] = eig(cov(lfpwr(freqs,:)'));
            for j = 1:size(rc.lfpwr,2)
            eva(j,:) =  squeeze(rc.lfpwr(freqs,j,1,:))' * E(:,end);
            evb(j,:) =  squeeze(rc.lfpwr(freqs,j,1,:))' * E(:,end-1);
            end
            subplot(3,1,2);
            imagesc(probes, rc.x(:,1),eva(:,probes));
            subplot(3,1,3);
            imagesc(probes,rc.x(:,1),evb(:,probes));
        end        
    elseif plottype == LFPFREQ;  %Look at Frequency spectrum, summed over stim
        subplot(1,1,1);
        freqs = find(rc.lfpfrq > freqrange(1) & rc.lfpfrq < freqrange(2));
        hold off;
        if ndims(rc.lfpwr) == 2
            lfpwr = rc.lfpwr(freqs,:)';
        elseif ndims(rc.lfpwr) == 3
            lfpwr = abs(squeeze(sum(rc.lfpwr(:,freqs,probes))))';
        elseif ndims(rc.lfpwr) == 4
            lfpwr = abs(squeeze(sum(sum(rc.lfpwr(freqs,:,yvals,probes),3),2)))';
        end
        [E,V] = eig(cov(lfpwr));
        V = diag(V);

        ev = lfpwr * E(:,end-1:end);
        if lineplot == 1
            subplot(2,1,1);
            hold off;
            plot(rc.lfpfrq(freqs),lfpwr);
            hold on;
            plot(rc.lfpfrq(freqs),E(:,end) .* mean(ev(:,2)),'k');            
            plot(rc.lfpfrq(freqs),E(:,end-1) .* mean(ev(:,2)),'k:');            
            title(sprintf('%.1f-%.1fHz',freqrange(1),freqrange(2)));
            subplot(2,1,2);
            plot(ev(:,2),ev(:,1),'o-');
        elseif ndims(rc.lfpwr) == 2
            imagesc(rc.lfpfrq(freqs),probes,rc.lfpwr(freqs,probes)');
        elseif ndims(rc.lfpwr) == 3
            subplot(3,1,1);
            imagesc(rc.lfpfrq(freqs),probes,abs(squeeze(sum(rc.lfpwr(:,freqs,probes))))');
            subplot(3,1,2);
            ev = lfpwr * E(:,end-1:end);
            if sum(ev(:,2)) < 0
                ev(:,2) = ev(:,2) .* -1;
            end
            plot(ev);
        elseif ndims(rc.lfpwr) == 4
            imagesc(rc.lfpfrq(freqs),probes,squeeze(sum(sum(rc.lfpwr(freqs,:,yvals,probes),3),2))');
            title(sprintf('%.1f-%.1fHz',freqrange(1),freqrange(2)));
        end
%        set(gca,'xlim',[75 125]);
    elseif plottype == LINEPLOT
        plot(rc.lfpresp);
    elseif plottype == MUVAR %% now with stimvar
    elseif plottype == PLOTONESTIMMU
        if isfield(rc,'tresps')
        if sum(xvals) == 0
            xvals = 1:size(rc.tresps,2);
        end
        nx = length(xvals);
        ny = length(yvals);
        tid = find(rc.times > timerange(1) & rc.times < timerange(2));
        allr = rc.tresps(tid,xvals,yvals,probes);
        if sumy
            ny = 1;
            resps = sum(rc.tresps,3);
            yvals = 1;
        else
            resps = rc.tresps;
        end
        cr = [min(resps(:)) max(resps(:))];
        for k = 1:ny
            for j = 1:length(xvals)
                subplot(nx,ny,(j) * ny + k-ny);
                rc.spkim.z = squeeze(resps(tid,xvals(j),yvals(k),probes))';
                rc.spkim.x = rc.times(tid)./10;
                rc.spkim.y = probes;
                imagesc(rc.spkim.x,rc.spkim.y,rc.spkim.z);
                caxis(cr);
                title(sprintf('%s=%.2f,%s=%.1fN%d',rc.Stimvals.xtype,rc.x(xvals(j),1),rc.Stimvals.ytype,rc.y(yvals(k)),rc.lfpn(xvals(j),yvals(k))))
            end
        end
        else
        end
    elseif ismember(plottype,[PLOTONESTIM, CSDSUM])
        nplots = 1;
        if isfield(rc,'tresps');
            if sum(xvals) == 0
                xvals = 1:size(rc.tresps,2);
            end
            subplot(2,1,1);
            nplots = 2;
            rc.spkim.z = squeeze(mean(rc.tresps(:,xvals(1),yvals,:),3))';
            rc.spkim.x = rc.times./10;
            rc.spkim.y = 1:size(rc.tresps,4);
            imagesc(rc.times./10,1:size(rc.tresps,4),rc.spkim.z);
        end
        if sum(xvals) == 0
            xvals = 1:size(rc.lfp,2);
        end

        nx = length(xvals)+showblank;
        if nx > 1
            xo = 0;
        else
            xo = nplots-1;
        end
        ny = length(yvals);
        allr = rc.lfp(tid,:,:,:);
        cr = [min(allr(:)) max(allr(:))];
        if sumy == 1
            ny = 1;
            allr = sum(rc.lfp(tid,:,yvals,:),3);
            yvals = 1;
        end
        if sumy == 2 & size(rc.lfp,3) > 1
            ny = 2;
            allr(:,:,1,:) = sum(rc.lfp(tid,:,yvals,:),3)./2;
            if length(yvals) == 4
            allr(:,:,2,:) = rc.lfp(tid,:,yvals(1),:) + rc.lfp(tid,:,yvals(2),:) - ...
                (rc.lfp(tid,:,yvals(3),:) + rc.lfp(tid,:,yvals(4),:));
            else
            allr(:,:,2,:) = rc.lfp(tid,:,yvals(1),:) - rc.lfp(tid,:,yvals(2),:);
            yvals = [yvals(1) yvals(2)];
            end
        end
        for k = 1:ny
            for j = 1:length(xvals)
                latency = LFPLatency(squeeze(rc.lfp(:,xvals(j),yvals(k),probes)),rc.lfptimes);
                rc.lfplatencies(j,k,:) = latency(:,3);
                subplot(nx+xo,ny,(j+xo) * ny + k-ny);
                rc.lfpim.z = squeeze(allr(:,xvals(j),yvals(k),probes))';
                if lineplot == 3
                    csd = CalcCSD(rc.lfpim.z,'smooth',csdsk);
                    imagesc(csd);
                    [maxc, maxt] = max(abs(csd'));
                    [a, maxp(j,k)] = max(maxc);
                    [a, fastp] = min(maxt);
                    range(j,k,:) = minmax(csd(:));
                    csds(j+(k-1)*nx,:,:) = csd';
                elseif lineplot == 1
                    plot(rc.lfptimes(tid)./10,rc.lfpim.z);
                else
                imagesc(rc.lfptimes(tid)./10,probes,rc.lfpim.z);
                caxis(cr);
                csds(j+(k-1)*nx,:,:) = rc.lfpim.z';
                end
                rc.lfpim.x = rc.lfptimes(tid)./10;
                rc.lfpim.y = probes;
                if diff(size(rc.x)) > 0
                    rc.x = rc.x';
                end
                if isfield(rc,'Stimvals')
                title(sprintf('%s=%.2f,%s=%.2f N%d',rc.Stimvals.xtype,rc.x(xvals(j),1),rc.Stimvals.ytype,rc.y(yvals(k)),rc.lfpn(xvals(j),yvals(k))))
                end
                if j < length(xvals)
                set(gca,'xticklabel',[]);
                end
                if k == 1
                    xlabel('Time (ms)');
                end
            end
        end
        if showblank
            subplot(nx+xo,ny,(nx+xo) * (ny));
            if lineplot == 3
                csd = CalcCSD(rc.lfpblank(tid,probes)','smooth',csdsk);
                imagesc(csd);
            else
                imagesc(rc.lfptimes(tid)./10, probes,rc.lfpblank(tid,probes)');
            end
            title(sprintf('Blank'));
        end
        if lineplot == 3 || lineplot == 0%CSD
            
            if lineplot ==3 % CSD, need to set color scale again
            cr = [min(range(:)) max(range(:))];
        for k = 1:ny
            for j = 1:length(xvals)
                subplot(nx+xo,ny,(j+xo) * ny + k-ny);
                caxis(cr);
            end
        end
            end
        if sumy == 2
            csdsum = squeeze(sum(abs(csds(1:nx,:,:))));
        else
            csdsum = squeeze(sum(abs(csds),1));
        end
            [a,b] = max(max(csdsum));
            [a, tmax] = max(csdsum(:,b));
            signs = sign(csds(:,tmax,b));
            csum = csdsum';
            csdsum = zeros(size(csdsum));
            if sumy == 2
                csddiff = zeros(size(csdsum));
                for j = 1:nx
                    csdsum = csdsum + squeeze(csds(j,:,:)) .* signs(j);
                    csddiff = csddiff + squeeze(csds(j+nx,:,:)) .* signs(j);
                end
            else
            for j = 1:size(csds,1)
                csdsum = csdsum + squeeze(csds(j,:,:)) .* signs(j);
            end
            end
            csdsum = csdsum';
            pmax = csdsum(:,tmax);
            id = find(pmax < 0); %probes with -ve csd
            upid = id(find(id < b));
            downid = id(find(id > b));
            if lineplot == 3
                p = probes(2):probes(end-1);
                typelabel = 'CSD';
            else
                p = probes;
                typelabel = 'LFP';
            end
            ip = p(1):0.1:p(end);
            if isempty(downid)
                rc.csdzc(1) = NaN;
            else
            id = [downid(1) downid(1)-1];
            rc.csdzc(1) = interp1(pmax(id),p(id),0);
            end
            if isempty(upid)
                rc.csdzc(2) = NaN;
            else
            id = [upid(end) upid(end)+1];
            rc.csdzc(2) = interp1(pmax(id),p(id),0);
            end
            rc.csdmax = p(b);
            if plottype == CSDSUM
                subplot(1,1,1);
            else
                GetFigure(figures{2});
            end
            if sumy ==2
                subplot(2,1,1);
            imagesc(csdsum);
                subplot(2,1,2);
            imagesc(csddiff');
            else
            imagesc(csdsum);
            end

            if interpolate
            isum =smooth2(rc.lfptimes(tid),probes(2):probes(end-1),csdsum,rc.lfptimes(tid),ip,[10 1]);
            imagesc(isum);
            [c,d] = max(max(isum'));
            title(sprintf('%s max at %.0f, %.1f\n',typelabel,p(b),ip(d)));
            rc.csdimax = ip(d);
            else
                title(sprintf('%s max at %.0f',typelabel,p(b)));
            end
        end
    elseif plottype == PLOTMONOC
        diffscale = 1;
        if isfield(rc,'lfplm') & isfield(rc,'lfprm')
       cm(1) = min([rc.lfpresp(:); rc.lfpblank(:); rc.lfplm(:); rc.lfprm(:)]);
       cm(2) = max([rc.lfpresp(:); rc.lfpblank(:); rc.lfplm(:); rc.lfprm(:)]);
        subplot(2,2,1);
        imagesc(rc.lfpresp);
        caxis(cm);
        title('binoc');
        subplot(2,2,2);
        imagesc(rc.lfpblank);
        caxis(cm);
        title('blank');
        subplot(2,2,3);
        imagesc(rc.lfplm);
        caxis(cm);
        title('Left');
        subplot(2,2,4);
        imagesc(rc.lfprm);
        caxis(cm);
        title('right');
        else
       cm(1) = min([rc.lfpresp(:); rc.lfpblank(:)]);
       cm(2) = max([rc.lfpresp(:); rc.lfpblank(:)]);
        subplot(2,2,1);
        imagesc(squeeze(rc.lfpresp(:,2,:)));
        caxis(cm);
        title('binoc');
        subplot(2,2,2);
        imagesc(rc.lfpblank);
        caxis(cm);
        title('blank');
        
        subplot(2,2,3);
        imagesc(squeeze(rc.lfpresp(:,1,:)));
        if diffscale
            title(sprintf('Left x %.1f',range(cm)./range(caxis)));
        else
            caxis(cm);
        title('Left');
        end

        
        subplot(2,2,4);
        imagesc(squeeze(rc.lfpresp(:,3,:)));
        if diffscale
            title(sprintf('Right x %.1f',range(cm)./range(caxis)));
        else
            caxis(cm);
        title('right');
        end
        end
    elseif plottype == XTPROBELFP
        PlotTimeCourseProbes(rc,probes,lineplot,sumy);
    elseif plottype == XTPROBELFPMU
        PlotTimeCourseProbesMU(rc,probes,lineplot,sumy, fixscale);
    end
    result = rc;
    return;
elseif isdir(name)
    result = ProcessDir(name, saveresult, argon{:});
    return;
end

if isempty(isrc)
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

if isempty(strfind(name,'.lfp.')) %if name is .lfp, only do .lfp files

    Clusters = {};

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
            result.types = res.types;
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
                if strmatch(res.types{1},{'Op' 'Pp'})
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
                result.types = res.type; 
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
    end
    result.isrc = isrc;

    result.Cluster(j).autocut = [C.autocut];
    result.Cluster(j).dprime = [C.dprime];
    if isfield(Expt.Header,'SpkStats')
        result.Cluster(j).SpkStats = Expt.Header.SpkStats;
    end
end

result.Stimvals.fz = Expt.Stimvals.fz;
result.Stimvals.et = Expt.Stimvals.et;
result.Stimvals.e2 = Expt.Stimvals.e2;
result.Stimvals.Ro = Expt.Stimvals.Ro;



clst = regexprep(name,'.c1.*','.cells.mat');
isolation = 5; %%quality must exceed this to count
if exist(clst,'file')
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
        cExpt.Header.probe = median(probe(trials(cid)));
        x = GetProbeSep(cExpt.Header);
        if ~isnan(x)
            cExpt.Header.probesep = x;
        end
        muid = find(probe(trials) > 0 & quality(trials) <= isolation);
        muExpt.Trials = muExpt.Trials(find(probe(trials) > 0 & quality(trials) <= isolation));
        muExpt.probes = probe(trials(find(probe(trials) > 0 & quality(trials) <= isolation)));
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
            if strmatch(res.types{1},{'Op' 'Pp'})
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
                    end
                    T.probe = mean(cExpt.probes);
                    T.probes = cExpt.probes;
                end
                if strmatch(res.type{1},{'Op' 'Pp'})
                    res.fit = FitExpt(res);
                    T.types = res.type;
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
        if saveresult
            if nsu > mintrials
            outname = strrep(name,'.c1.',['.cell' num2str(j) '.']);
            fprintf('Saving %s\n',outname);
            save(outname,'cExpt');
            end
            if nmu > mintrials
            outname = strrep(name,'.c1.',['.mu' num2str(j) '.']);
            fprintf('Saving %s\n',outname);
            save(outname,'muExpt');
            end
        end
        end
    end
end

if ~isempty(sptrig)
    sptrig.name = name;
end
result.resps = resps;
result.means = means;
result.x = xvals;
result.y = yvals;
if strmatch(Expt.Stimvals.et,'Dc');
result.Stimvals.xtype = Expt.Stimvals.e2;
result.Stimvals.ytype = Expt.Stimvals.et;
else
result.Stimvals.xtype = Expt.Stimvals.et;
result.Stimvals.ytype = Expt.Stimvals.e2;
end
result.Stimvals.et = Expt.Stimvals.et;
result.Stimvals.e2 = Expt.Stimvals.e2;
if isrc
    [result.rcprefs, x] = GetNetSpkPref(result,[timerange(3) timerange(4)]);
    result.netspksd = x.netsd;
end
id = find(ismember([Expt.Trials.Trial],[Expt.Header.BlockStart]));

for j = 1:length(id)
    result.Header.BlockStart(j) = Expt.Trials(id(j)).Start(1);
    result.Header.BlockTrial(j) = Expt.Trials(id(j)).Trial;
end
result.Header.eds = eds;
if isfield(LFP.Header,'LFPamps')
result.Header.LFPamps = LFP.Header.LFPamps;
end
if isfield(Expt.Header,'depths')
    result.Header.depths = Expt.Header.depths;
end

if ischar(name)
result.Header.filename = name;
end
if ispsych
    result.pp = ExptPsych(Expt);
end

if ispsych & isrc
    [result.pk.kernel, details]= CalcPk(bExpt);
    result.pk.xv = details.xv;
    result.pk.n = details.n;
    result.pk.signal = details.signal;
    result.pk.signaln = details.signal;
    result.RespDir = [bExpt.Trials.RespDir];
end
result.ispsych = ispsych';
end



ename = strrep(name,['.' cs],'.lfp.');
if ~exist('LFP','var') & exist(ename,'file')
    load(ename);
end
if exist('LFP','var')
    rescale = find(scales ~=1);
    if isrc
        LFP = CheckLFP(LFP,'fix');
        [rc, blfp] = PlotRevCorAny(LFP,'lfp', rcargs{:});
        [result.lfp, mresp, result.lfpblank, a] = CombineRCLFP(rc, scales);
        if isfield(rc.sdfs,'triallfp')
            result.lfptrial = rc.sdfs.triallfp;
        end
        if ~isfield(rc,'x')
            rc.x = a.x;
            result.x = rc.x;
        end
        if isfield(a,'lm')
            result.lfplm = a.lm;
        end
        if isfield(a,'rm')
            result.lfprm = a.rm;
        end
        if isfield(a,'uc')
            result.lfpuc = a.uc;
        end
        result.mlfp = mresp;
        result.lfpmean = rc.sdfs.alllfp;
        result.lfptimes = rc.sdfs.lfptimes;
        result.lfpn = rc.sdfs.n;
        if ~isfield(result,'sdfs')
        result.sdfs = rc.sdfs;
        end
        ck = CheckLFP(blfp);
        len = mean(ck.lens);
        if isfield(rc,'cp')
                result.lfpwr = rc.cp.lfpwr;
                result.lfpfreq = (0:len-1)/(len * blfp.Header.LFPsamplerate);
                result.lfpsigpwr = rc.cp.siglfpwr; %lfp power spectrum by signal
                result.lfpsigt = rc.cp.siglfpt; %mean trial timecourse, by siggal
        elseif isfield(blfp.Trials,'FTlfp') %sometimes this is omitted to save memory
            if length(unique(ck.lens)) == 1
                result.lfpwr = mean(abs(cat(3,blfp.Trials.FTlfp)),3);
            end
        else
            len = mean(ck.lens);
            result.lfpwr = mean(abs(fft(cat(3,blfp.Trials.LFP))),3);
            result.lfpfreq = (0:len-1)/(len * blfp.Header.LFPsamplerate);
        end
        if length(rescale)
            for j = find(scales ~= 1)
                result.lfpwr(:,j) = result.lfpwr(:,j) .* scales(j);
            end
        end
        len = size(LFP.Trials(1).LFP,1);
        result.lfpfrq = (0:len-1)/(len * LFP.Header.LFPsamplerate);
    else
        LFP = CheckLFP(LFP,'fix');
        if exist('Expt','var') %may not exist if doing LFP only
        LFP = FillTrials(LFP,Expt.Stimvals.et);
        end
        lfp = PlotLFP(LFP,'lfp');
        if length(rescale)
            for j = rescale
                lfp.lfp(:,:,:,j) = lfp.lfp(:,:,:,j) .* scales(j);
                lfp.lfppower(:,:,:,j) = lfp.lfppower(:,:,:,j) .* scales(j);
            end
        end
        result.lfp = lfp.lfp;
        result.lfptimes = lfp.lfptimes;
        result.lfpwr = lfp.lfppower;
        result.lfpfrq = lfp.ftfrq;
        result.lfpn = lfp.n;
        if isfield(lfp,'extra')
            id = strmatch('Blank',lfp.extra.label);
            if length(id) && lfp.extra.n(id) > 0
                result.lfpblank = lfp.extra.lfp{id};
                result.lfpblankn = lfp.extra.n(id);
                result.lfpblankpwr = lfp.extra.lfppower{id};
            end
        end
    end
end

function rc = PlotOneLFP(rc, timerange, probes, varargin)

dvdt = 0;
csd = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'csd',3)
        csd = 1;
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            csdsk = varargin{j};
        end
    elseif strncmpi(varargin{j},'dvdt',4)
        dvdt = 1;
    end
    j = j+1;
end
tid = find(rc.lfptimes > timerange(1) & rc.lfptimes < timerange(2));
rc.lfpim.z = squeeze(rc.lfp(tid,1,1,probes))';

rc.lfpim.x = rc.lfptimes(tid)./10;
if dvdt
    rc.lfpim.z = diff(rc.lfpim.z,[],2);
    rc.lfpim.x = rc.lfptimes(tid(2:end))./10;        
end
if csd
    rc.lfpim.z= CalcCSD(rc.lfpim.z,'smooth',csdsk);
    csdv = CalcCSDPeaks(rc.lfpim.z);
    rc.lfpim.y = rc.probes(1:end-2);
    rc.csdzc = csdv.zc;
else
    rc.lfpim.y = rc.probes;
end

hold off;
imagesc(rc.lfpim.x,rc.lfpim.y,rc.lfpim.z);
[v,tmax] = max(rc.lfpim.z');
[a,probe] = max(v);

if dvdt == 0 && csd == 0
pid = find(rc.lfptimes < 500 & rc.lfptimes > 0);
base = squeeze(mean(rc.lfp(pid,1,1,probes)));
zid = tmax(probe)-30:tmax(probe)+40;
zid = zid(find(zid > 0 & zid <= size(rc.lfpim.z,2)));
for j = zid
    uvp(:,j) = smooth(rc.lfpim.z(:,j),1,'gauss');
    vp(:,j) = uvp(:,j) - base;
    id = find(vp(:,j) <0);
    if length(id)
        zc(j) = id(end);
    else
        zc(j) = NaN;
    end
    id = find(uvp(:,j) <0);
    if length(id)
        uzc(j) = id(end);
    else
        uzc(j) = NaN;
    end
end
hold on;
plot(rc.lfptimes(tid(zid))./10,zc(zid),'r:');
plot(rc.lfptimes(tid(zid))./10,uzc(zid),'m:');
rc.lfpzc = zc(zid);
sumresp = sum(rc.lfpim.z);
respvar = var(rc.lfpim.z);
hvid = find(respvar > max(respvar)/2);
plot([hvid(1) hvid(1)],[1 24],'b:');
plot([hvid(end) hvid(end)],[1 24],'b:');
vid = find(respvar > max(respvar)/10);
plot([vid(1)+5 vid(1)+5],[1 24],'b:');
plot([vid(end) vid(end)],[1 24],'b:');
rc.lfpztimes= rc.lfptimes(tid(zid));
rc.maxprobe = probe;
if hvid(end) > length(zc)
else
title(sprintf('ZC at %d,%d',zc(vid(1)+5),zc(hvid(end))));
end
elseif csd
        hold on; 
    plot(rc.lfpim.x(csdv.t),csdv.zc,'r:');
end
%plot(vp);

function [revp, minp, details] = GetFrameReverse(rc, probes)

mlfp = rc.mlfp;
if ndims(mlfp) == 3
    mlfp = squeeze(sum(mlfp,2));
end

for j = probes
    [a,b] = famp(rc.lfptimes, mlfp(:,j),1/166.66);
    phases(j) = angle(b);
    amps(j) = abs(b);
end
[a,b] = max(abs(diff(phases)));
if a >  pi
    id = find(phases < phases(b-1));
    id = find(phases < 0);
    phases(id) = phases(id) + 2 * pi;
end
details.amps = smooth(amps,2,'gauss');
[a,b] = max(details.amps);
details.maxprobe = probes(b);
[a,b] = min(details.amps(1:b));
minp = probes(b);


[y,x] =smhist(phases);
[a,b] = max(y);
details.phase = x(b); %dominant phase
sp = details.phase +pi/2; %90 deg ahead - halfway to 180;
details.phases = smooth(phases,2,'gauss');
id = find(details.phases > details.phase);
lastid = find(diff(id) > 1); %find break;
if length(lastid)
id = id(1:lastid(1));
end
[a,b] = max(abs(diff(details.phases(id))));
revp = probes(b);

if sp > max(details.phases)
    revp = interp1(details.phases(id),probes(id),sp,'pchip','extrap');
    if revp < 0 || revp > max(probes)
        revp = NaN;
    end
else
    revp = interp1(details.phases(id),probes(id),sp);
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

function WriteRFLoc(All)

for j = 1:length(All.exps)
    if isfield(All.exps{j},'types')
        types{j} = All.exps{j}.types{1};
    elseif isfield(All.exps{j},'type')
        types{j} = All.exps{j}.type{1};
    else
        types{j} = 'e0';
    end
    if isempty(types{j})
        types{j} = 'e0';
    end
    if isfield(All.exps{j},'spkres')
        for k = 1:length(All.exps{j}.spkres)
            cells(j,k) = ~isempty(All.exps{j}.spkres{k});
        end
    end
end 

probes = 8;

ip = strmatch('Pp',types);
io = strmatch('Op',types);
for p = 1:probes
for j = 1:length(ip)
    Pp(j,p) = All.exps{ip(j)}.fit{p}.mean;
    Ppsd(j,p) = All.exps{ip(j)}.fit{p}.sd;
    Por(j) = All.exps{ip(j)}.Stimvals.Ro;
end
for j = 1:length(io)
    fit = All.exps{io(j)}.fit{p};
    Op(j,p) = fit.mean;
    Opsd(j,p) = fit.sd;
    Oor(j) = All.exps{io(j)}.Stimvals.Ro;
%    Opve(j) = fit.ve;
end

Ro = median([Por Oor]);
x(1) = mean(Op);
x(2) = mean(Pp); 
rf(:,p) = op2xy(x, Ro);
end
clear Pp;
clear Op;
for j = 1:length(ip)
for p = 1:length(All.exps{ip(j)}.spkres)
    S = All.exps{ip(j)}.spkres{p};
    if ~isempty(S)
        Pp(j,p) = S.fit.mean;
    end
end
end
for j = 1:length(io)
for p = 1:length(All.exps{io(j)}.spkres)
    S = All.exps{io(j)}.spkres{p};
    if ~isempty(S)
        Op(j,p) = S.fit.mean;
    end
end
end

id = find(Pp ~= 0 & Op ~= 0);

function [result, sptrig] = ProcessDir(name, saveresult, varargin)
    d = dir(name);
    nex = 1;
    monkey = 'lem';
    ntrig = 0;
    result = [];
    for j = 1:length(d);
        if strfind(d(j).name,'.p1c1')
            fprintf('%s\n',d(j).name);
                        result.names{nex} = strrep(d(j).name,'.p1c1','.c1');
            [result.exps{nex}, a] = PlotAllProbes([name '/' result.names{nex}],varargin{:});
            if ~isempty(a)
                sptrig{nex} = a;
                ntrig = ntrig+1;
            end
            nex = nex+1;
        end
    end
    if isempty(result)
        fprintf('No Expts for %s\n',name);
    end
    if saveresult && exist('result','var')
        [a,b] = fileparts(name);
        AllExpts = result;
        savename = [name '/' monkey b '.all.mat'];
        fprintf('Saving %s\n',savename);
        save(savename,'AllExpts')
        if ntrig
            savename = [name '/' monkey b '.sptrig.mat'];
            fprintf('Saving %s\n',savename);
            save(savename,'sptrig')
        end
    end

function [netspk, meanspk, details] = GetNetSpk(rc,timerange)

tid = find(rc.times > timerange(1) & rc.times < timerange(2));
if size(rc.spkn,2) > 1
xid = find(prod(rc.spkn') > 0);
else
xid = find(rc.spkn > 0);
end
details.xv = rc.x(xid,1);
meanresp = repmat(mean(rc.tresps(tid,xid,:,:),2),[1 length(xid) 1 1]);
netspk = (rc.tresps(tid,xid,:,:) - meanresp)./meanresp;
netspk = shiftdim(mean(netspk,1)); % sum over t to get mean rate diff
netspk = netspk .* diff(timerange)./10000;
meanrate = squeeze(mean(mean(mean(meanresp,1),2),3));
if isfield(rc.Stimvals,'fz')
    meanspk = meanrate./rc.Stimvals.fz;
else
    meanspk = meanrate/60;
end
if isfield(rc,'blankresp')
    m = squeeze(meanresp(:,1,1,:));
    bk = (rc.blankresp(tid,:) - m)./m;
    bk = mean(bk) .* diff(timerange)./10000;
    details.blankspk = bk+meanspk';
end
for j = 1:size(netspk,1)
    for k = 1:size(netspk,2)
        netspk(j,k,:) = squeeze(netspk(j,k,:)) + meanspk;
    end
end
        
function [prefs, x] = GetNetSpkPref(rc, timerange)

[netspk, meanspk, details] = GetNetSpk(rc,timerange);
[prefs, x] = NetSpkPref(netspk, details.xv);
if isfield(details,'blankspk')
    x.blankspk = details.blankspk;
end

function [prefs,x] = NetSpkPref(netspk, xv, varargin)
%
% Calculate pref resp from Max value, and from direction
% of mean vector (angles doubled - for ori, not dir)
noneg = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'zerobase',5)
        noneg = 1;
    end
    j = j+1;
end

if noneg
    for j = 1:size(netspk,3)
        netspk(:,:,j) = netspk(:,:,j) - min(min(netspk(:,:,j)));
    end
end
x.netspk = netspk;
sumy = squeeze(mean(x.netspk,2));
        [x.maxr, prefs] = max(sumy);
        x.netsd = std(sumy);
        oris = 2 * xv(:,1) .* pi/180; % double angles for CV
        oris = repmat(oris,1,size(sumy,2));
        x.vsum = (sum(sumy .* cos(oris)) + i * sum(sumy .* sin(oris)))./sum(abs(sumy));
        x.vsum = abs(x.vsum) .* (cos(angle(x.vsum)/2) + i * sin(angle(x.vsum)/2));
        
function [sdfs, mresp, blank, details] = CombineRCLFP(rc, scales)

for j = 1:size(rc.sdfs.lfp,1)
    for k = 1:size(rc.sdfs.lfp,2)
    for ch = 1:size(rc.sdfs.lfp{j,k},2)
        sdfs(:,j,k,ch) = rc.sdfs.lfp{j,k}(:,ch) .* scales(ch);
    end
    details.y(k) = rc.sdfs.y(1,k);
    end
    details.x(j) = rc.sdfs.x(j,1);
end

mresp = squeeze(mean(sdfs,2));
nsdf = prod(size(rc.sdfs.lfp));

blanks = [];
blank = [];
blankn = 0;
no=0;
for k = 1:length(rc.sdfs.extras) 
    if isempty(rc.sdfs.extras{k})
        no = no+1;
    else
    blanks(:,k-no,:) = rc.sdfs.extras{k}.lfp;
    blankn(:,k-no,:) = rc.sdfs.extras{k}.lfpn;
    id = strmatch(rc.sdfs.extras{k}.label,{'Blank' 'Right' 'Left' 'Uncorr'});
    if length(id)
        types(k-no) = id(1);
    else
        types(k-no) = 0;
    end
    end
end
if ~isempty(blanks)
for j = 1:size(blanks,3)
    blanks(:,j) = blanks(:,j) .* scales(j);
end


if ndims(mresp) > 2
% including blanks is a tad tricky for this
    %    mresp = cat(2,mresp .* nsdf, blanks);
%    mresp = squeeze(mean(mresp,2));
    bmresp = squeeze(mean(mresp,2));
else
    if ~isempty(blanks)
    mresp = (mresp .* nsdf + squeeze(sum(blanks,2)))./(nsdf+length(types));
    end
    bmresp = mresp;
%    blank = squeeze(mean(blanks(:,find(types== 1),:),2));
end
if sum(types == 1)
blank = squeeze(mean(blanks(:,find(types== 1),:),2)) - bmresp;
end

for j = 1:size(blank,2)
 %   blank(:,j) = blank(:,j) - mean(blank(:,j));
end
if sum(types == 2)
details.rm = squeeze(mean(blanks(:,find(types== 2),:),2)) - bmresp;
end
if sum(types == 3)
details.lm = squeeze(mean(blanks(:,find(types== 3),:),2)) - bmresp;
end
if sum(types == 4)
details.uc = squeeze(mean(blanks(:,find(types== 4),:),2)) - bmresp;
end
end
details.blankn = blankn;

for j = 1:size(rc.sdfs.lfp,1)
    sdfs(:,j,:,:) = squeeze(sdfs(:,j,:,:)) - mresp;
end


function NextFrame(a,b,step)

D = get(gcf,'UserData');
t = D.inow;
if step == 0
    t = t+1;
    if t > length(D.mtimes)
        D.mtimes = [D.mtimes D.mtimes(end)+mean(diff(D.mtimes))];
        t = length(D.mtimes);
    end
elseif step < 0
    t = t+step;
    if t < 1
        t = 1;
    end
end
PlotTimeSlice(D,D.mtimes(t));   
D.inow = t;
set(gcf,'UserData',D);

function FrameSlider(a,b,step)

D = get(gcf,'UserData');
t = get(a,'value') .* 10;
[dt, it] = min(abs(t-D.mtimes)); 
PlotTimeSlice(D,D.mtimes(it));
drawnow;
D.inow = it;
set(gcf,'UserData',D);

function PlayFrames(a,b,step)

D = get(gcf,'UserData');
t = D.inow;
set(findobj(gcf,'Tag','PlayStop'),'value',0);
if step > 0
    if D.inow >= length(D.mtimes)
        D.inow = 1;
    end
    t = [D.inow:step:length(D.mtimes)];
elseif step < 0
    t = [D.inow:step:1];
end
for it = t
    stop = get(findobj(gcf,'Tag','PlayStop'),'value');
    if stop
        break;
    end
    PlotTimeSlice(D,D.mtimes(it));
    drawnow;
end
D.inow = it;
set(gcf,'UserData',D);

function res = PlotMovie(res, delay, mtimes)

sustep = 1000;
lfpstep = 10000;
    
if res.lineplot == 3
    for j = 1:size(res.lfpresp,2)
      lfpresp(:,j,:) = CalcCSD(squeeze(res.lfpresp(:,j,res.probes))','smooth',res.csdsk)';
    end
    res.lfpresp = lfpresp;
    if isfield(res,'lfprespa')
        for j = 1:size(res.lfprespa,2)
            lfpresp(:,j,:) = CalcCSD(squeeze(res.lfprespa(:,j,res.probes))','smooth',res.csdsk)';
        end
        res.lfprespa = lfpresp;
    end
        
end

if isfield(res,'spkresps')
res.surange(1) = min(res.spkresps(:));
res.surange(2) = max(res.spkresps(:));
end
if isfield(res,'lfpresp')
    
res.lfprange(1) = min(res.lfpresp(:));
res.lfprange(2) = max(res.lfpresp(:));
end
res.mtimes = mtimes;
PlotTimeSlice(res, res.mtimes(1));


fb = findobj(gcf,'Tag','NextFrame');
if isempty(fb)
    bp = [10 10 40 20];
    uicontrol(gcf, 'style','pushbutton','string','+','Position',bp,...
        'Tag','NextFrame','Callback',{@NextFrame, 0});
end
fb = findobj(gcf,'Tag','LastFrame');
if isempty(fb)
    bp = [50 10 40 20];
    uicontrol(gcf, 'style','pushbutton','string','-','Position',bp,...
        'Tag','LastFrame','Callback',{@NextFrame, -1});
end
fb = findobj(gcf,'Tag','PlayFwd');
if isempty(fb)
    bp = [90 10 40 20];
    uicontrol(gcf, 'style','pushbutton','string','>>','Position',bp,...
        'Tag','PlayFwd','Callback',{@PlayFrames,  1});
end
fb = findobj(gcf,'Tag','PlayBack');
if isempty(fb)
    bp = [130 10 40 20];
    uicontrol(gcf, 'style','pushbutton','string','<<','Position',bp,...
        'Tag','PlayBack','Callback',{@PlayFrames, -1});
end
fb = findobj(gcf,'Tag','PlayStop');
if isempty(fb)
    bp = [170 10 60 20];
    uicontrol(gcf, 'style','check','string','Stop','Position',bp,...
        'Tag','PlayStop');
    bp(1) = bp(1)+bp(3);
    bp(3) = 200;
    res.slider = uicontrol(gcf, 'style','slider','string','t','Position',bp,...
        'Tag','PlaySlider','Min',0,'max',200,'value',0,'sliderstep',[0.005 0.025],'CallBack',@FrameSlider);
else
    set(fb(1),'value',0);
end
set(gcf,'UserData',res);
ofig = gcf;
subplot(2,1,1); hold off;
subplot(2,1,2); hold off;
for ti = 1:length(res.mtimes)
    stop = get(findobj(gcf,'Tag','PlayStop'),'value');
    if stop
        ti = NaN;
    else


        set(0,'currentfigure',ofig);
    PlotTimeSlice(res,res.mtimes(ti));
    res.inow = ti;
    set(gcf,'UserData',res);
    drawnow;
    if delay > 1
        pause
    else
       pause(delay);
    end
    end
end

function PlotTimeSlice(res,t);

nplots = 1;
if isfield(res,'lfprespa')
    subplot(2,1,1);
    nplots = 2;
    [a,b] = min(abs(t-res.lfptimes));
    if res.plot

    imagesc(res.x(:,1),res.probes,squeeze(res.lfprespa(b,:,:))' - squeeze(res.lfpresp(b,:,:))');
    else
        imagesc(res.x(:,1),res.probes,squeeze(res.lfprespa(b,:,:))');
    end
    caxis(res.lfprange);
    if isfield(res,'titlestr')
        title(sprintf('%.1f ms %s',res.lfptimes(b)./10,res.titlestr{2}));
    else
        title(sprintf('%.1f ms',res.lfptimes(b)./10));
    end
elseif isfield(res,'spkresps')
    subplot(2,1,1);
    nplots = 2;
    [a,b] = min(abs(t-mean(res.times,1)));
    imagesc(squeeze(res.spkresps(b,:,:))');
    caxis(res.surange);
    title(sprintf('%.1f ms',res.times(1,b)./10));
end

if isfield(res,'lfpresp')
    if nplots == 1
        subplot(1,1,1);
    else
        subplot(2,1,2);
    end
    [a,b] = min(abs(t-res.lfptimes));
    if res.plot & isfield(res,'lfprespa')
        title(sprintf('%.1f ms',res.lfptimes(b)./10));
        imagesc(squeeze(res.lfprespa(b,:,:)+res.lfpresp(b,:,:))'./2);
    else
        imagesc(res.x(:,1),res.probes,squeeze(res.lfpresp(b,:,:))');
    end
    caxis(res.lfprange);
    if isfield(res,'titlestr')
        title(sprintf('%.1f ms %s',res.lfptimes(b)./10,res.titlestr{1}));
    else
        title(sprintf('%.1f ms',res.lfptimes(b)./10));
    end
    if isfield(res,'slider')
        t = res.lfptimes(b)./10;
        if t <= get(res.slider,'Max') & t >= get(res.slider,'Min')
            set(res.slider,'value',res.lfptimes(b)./10);
        end
    end
end


function PlotTimeCourseProbes(rc, probes, style, sumy)

showvar = 1;

[nr,nc] = Nsubplots(length(probes));
if sumy == 2
    ny = ceil(length(rc.yvals)/2);
    yv = rc.yvals(1:ny);
    zv = rc.yvals((ny+2):length(rc.yvals));
end

cr = [min(rc.lfp(:)) max(rc.lfp(:))];
for probe = probes
    subplot(nr,nc,probe);
    if sumy == 1
        if style == 1
            plot(rc.lfptimes./10,rc.lfpresp(:,:,probe)');
            if showvar
                v = var(rc.lfpresp(:,:,probe),[],2);
                scale = max(max(rc.lfpresp(:,:,probe)))./max(v);
                hold on;
                plot(rc.lfptimes./10,v.*scale,'k');
                hold off;
                end
        else
            imagesc(rc.x(:,1),rc.lfptimes./10,rc.lfpresp(:,:,probe));
            caxis(cr);
        end
    elseif sumy == 2  %diff
        if style == 1
            plot(rc.lfptimes./10,sum(rc.lfp(:,:,yv,probe),3)-sum(rc.lfp(:,:,zv,probe),3));
        else
            imagesc(rc.x(:,1),rc.lfptimes./10,sum(rc.lfp(:,:,yv,probe),3)-sum(rc.lfp(:,:,zv,probe),3));
            caxis(cr);
        end
    end
end

function PlotTimeCourseProbesMU(rc, probes, style, sumy, fixscale)

[nr,nc] = Nsubplots(length(probes));
if sumy == 2
    ny = ceil(length(rc.yvals)/2);
    yv = rc.yvals(1:ny);
    zv = rc.yvals((ny+2):length(rc.yvals));
end
cr = [min(rc.tresps(:)) max(rc.tresps(:))];
for probe = probes
    subplot(nr,nc,probe);
    if sumy == 1
    if style == 1
    plot(rc.times./10,rc.spkresps(:,:,probe)');
    else
    imagesc(rc.x(:,1),rc.times./10,rc.spkresps(:,:,probe));
    if fixscale
            caxis(cr);
    end
    end
    elseif sumy == 2
        if style == 1
            plot(rc.times./10,sum(rc.tresps(:,:,yv,probe),3)-sum(rc.tresps(:,:,zv,probe),3));
        else
            imagesc(rc.x(:,1),rc.times./10,sum(rc.tresps(:,:,yv,probe),3)-sum(rc.tresps(:,:,zv,probe),3));
    if fixscale
            caxis(cr);
    end
        end
    end
end

function PlotTimeCourse(rc, probe, style, sumy)

if sumy == 0 
    nr = size(rc.lfp,3);
    else
    nr = 1;
end
for j = 1:nr
if isfield(rc,'lfp')
    subplot(2,nr,j*2-1);
    if style == 1
        if sumy
            plot(rc.lfptimes./10,rc.lfpresp(:,:,probe));
        else
            plot(rc.lfptimes./10,rc.lfp(:,:,j,probe));
        end
    else
        if sumy
        imagesc(rc.x(:,1),rc.lfptimes./10,rc.lfpresp(:,:,probe));
        else
        imagesc(rc.x(:,1),rc.lfptimes./10,rc.lfp(:,:,j,probe));
        end
    end
    subplot(2,nr,j*2);
else
    subplot(1,1,1);
end
if style ==1 
    if sumy
        plot(rc.times./10,rc.spkresps(:,:,probe).*rc.means(probe));
    else
        plot(rc.times./10,rc.tresps(:,:,j,probe).*rc.means(probe));
    end
else
    if sumy
        imagesc(rc.spkresps(:,:,probe));
    else
        imagesc(rc.x(:,1),rc.times./10,rc.tresps(:,:,j,probe));
    end
end
end

function [eigs, resps, blr, details] = PlotLFPEig(rc,lineplot, sumy, setsign)

plottype = 0;
eigs = [];
resps = [];
blr = [];
usey = 0;
details = [];
mkfit = 1;

if size(rc.lfp,3) == 2 && strcmp(rc.Stimvals.ytype,'ce')
    usey = 1;
end
tid = rc.tid;
blankr = rc.lfpblank(tid,:);

if lineplot == 3
    for j = 1:size(rc.lfp,2)
        lfpresp(:,j,:) = CalcCSD(squeeze(rc.lfp(tid,j,1,rc.probes)),'flip','smooth',rc.csdk);
    end
    blankr = CalcCSD(rc.lfpblank(tid,rc.probes),'flip','smooth',rc.csdk);
elseif usey == 1
    lfpresp = reshape(rc.lfp(tid,:,:,:),[length(tid) size(rc.lfp,2)*size(rc.lfp,3) size(rc.lfp,4)]);
else
    lfpresp = squeeze(sum(rc.lfp(tid,:,:,:),3)); %% pool all Y to get eigenvectors
end
if ndims(lfpresp) == 3 && lineplot ~= 3
lfpresp(:,end+1,:) = rc.lfpblank(tid,:);
end

for j = 1:size(lfpresp,2);
    details.blankresp(j,:) = sum(squeeze(lfpresp(:,j,:)) .* blankr);
end

if size(lfpresp,2) < 2
    fprintf('Need at least 2 resps to calculated Eigenvectors\n');
    return;
end

for probe = 1:size(lfpresp,3)
    [A,B] = eig(cov(lfpresp(:,:,probe)'));
    eigs(:,probe) = A(:,end);
    beigs(:,probe) = A(:,end-1);
    ceigs(:,probe) = A(:,end-2);
end

[a,b] = max(var(eigs));

%
% need to determine sign.
% could use dissimilarity from blank

if 0
for probe = 1:size(rc.lfp,4)
    sgn = sign(sum(eigs(:,probe) .* eigs(:,b)));
    eigs(:,probe) = eigs(:,probe) .* sgn;
end
end

if ~isempty(rc.lfpblank) & setsign == 0 %sign set by eig * blank
    for k = 1:size(rc.lfpblank,3)
        xresps(k,:) = sum(rc.lfpblank(tid,:,k) .* eigs);
    end
    for probe = 1:size(rc.lfp,4)
        eigs(:,probe) = eigs(:,probe) .* -sign(xresps(1,probe));
    end
elseif setsign == 1 % sign set by similarity to previous probe
    for probe = 1:size(lfpresp,3)
        if probe > 1
        ds = sum(eigs(:,probe) .* eigs(:,probe-1));
        if ds < 0
            eigs(:,probe) = eigs(:,probe) .*  -1;
        end
        end

    end
end

nx = size(rc.lfp,2);
for probe = 1:size(rc.lfp,4)
        for k = 1:nx
            for m = 1:size(rc.lfp,3)
                resps(k, m,probe) = sum(rc.lfp(tid,k,m,probe) .* eigs(:,probe));
                bresps(k, m,probe) = sum(rc.lfp(tid,k,m,probe) .* beigs(:,probe));
                cresps(k, m,probe) = sum(rc.lfp(tid,k,m,probe) .* ceigs(:,probe));
            end
        end
        xc = corrcoef(squeeze(sum(resps(1:nx,:,probe),2)),details.blankresp(1:nx,probe));
        probesign(probe) = xc(1,2);
        if ~isempty(rc.lfpblank)
            blr(1,probe) = sum(rc.lfpblank(tid,probe) .* eigs(:,probe));
            blr(2,probe) = sum(rc.lfpblank(tid,probe) .* beigs(:,probe));
            blr(3,probe) = sum(rc.lfpblank(tid,probe) .* ceigs(:,probe));
        end
end

if setsign == 3 % sign set by similarity of tuning curve to blankresp tuning curve
    id = find(probesign > 0);
    eigs(:,id) = eigs(:,id) * -1;
    resps(:,:,id) = resps(:,:,id) * -1;
end

if size(resps,2) == 4
    if sumy > 0 & lineplot == 2 %show scatterplot of eigs
        aresp = squeeze(sum(resps,2))';
        bresp = squeeze(sum(bresps,2))';
        cresp = squeeze(sum(cresps,2))';
    else
        aresp = squeeze(sum(resps(:,[1 2],:),2))';
        bresp = squeeze(sum(resps(:,[3 4],:),2))';
    end
    arange = [min([aresp(:); bresp(:)]) max([aresp(:); bresp(:)])];
    brange = arange;
    plottype = 1;
elseif size(resps,2) == 2 && sumy %COmpare 2 sets
    aresp = squeeze(mean(resps(:,:,:),2))';
    bresp = eigs';
    arange = [min(aresp(:)) max(aresp(:))];
    brange = [min(bresp(:)) max(bresp(:))];
elseif size(resps,2) == 2
    aresp = squeeze(resps(:,1,:))';
    bresp = squeeze(resps(:,2,:))';
    arange = [min([aresp(:); bresp(:)]) max([aresp(:); bresp(:)])];
    brange = arange;
    if strmatch(rc.Header.exptype,'OPRC') & mkfit
        for j = 1:size(aresp,1)
        x = [rc.x(:,1)' Inf];
        y = [aresp(j,:) xresps(j)];
        bfits(j) = FitGauss(x,y,'freebase');
        z = [bresp(j,:) xresps(j)];
        wfits(j) = FitGauss(x,z,'freebase');
        end
        details.fitratios = [wfits.amp]./[bfits.amp];
        details.bwvar = std(aresp')./std(bresp');
    end
    aresp = squeeze(sum(resps(:,1,:),2))';
    bresp = eigs';
    arange = [min(aresp(:)) max(aresp(:))];
    brange = [min(bresp(:)) max(bresp(:))];
elseif size(resps,2) == 1 
    aresp = squeeze(mean(resps(:,:,:),2))';
    bresp = squeeze(mean(bresps(:,:,:),2))';
    arange = [min(aresp(:)) max(aresp(:))];
    brange = [min(bresp(:)) max(bresp(:))];
end
[maxr, prefs] = max(squeeze(sum(resps,2)));
details.prefs = prefs;
details.beigs = beigs;
details.ceigs = ceigs;
details.bresps = bresps;
details.cresps = cresps;
if lineplot == 1
subplot(2,1,1);
plot(rc.x(:,1),aresp);
subplot(2,1,2);
if size(resps,2) == 2
    plot(rc.x(:,1),bresp);
else
    plot(rc.lfptimes./10,bresp);
end
elseif lineplot == 4 %scatterplot, all on one axis
    subplot(1,1,1);
    hold off;
    colors = 'rgbcymkrgbcymkrgbcymk';
    symbols = '.ovsn';
    for j = 1:size(rc.lfp,4)
        k = ceil(j/4);
        sym = 1+mod(j-1,4);
        plot(aresp(j,:),bresp(j,:),symbols(sym),'color',colors(k));
        hold on;
    end
elseif lineplot == 2 %scatterplot, one for each probe
    [nr,nc] = Nsubplots(size(rc.lfp,4));
    for j = 1:size(rc.lfp,4)
        subplot(nr,nc,j);
        hold off;
        [a,b] = fit_bothsubj2error(aresp(j,:),bresp(j,:));
        plot(aresp(j,:),bresp(j,:),'o');
        title(sprintf('Slope %.2f',b));
        refline(b,a);
        hold on;
        plot(blr(1,j),blr(2,j),'o','markerfacecolor','b');
        plot(aresp(j,:),cresp(j,:),'ro');
    end
elseif plottype == 1  && lineplot >= 0
subplot(3,1,1);
imagesc(rc.x(:,1),1:size(rc.lfp,4),aresp);
caxis(arange);
subplot(3,1,2);
imagesc(rc.x(:,1),1:size(bresp,1),bresp);
hold on;
plot(rc.x(prefs,:),1:size(bresp,1),'w:');
caxis(brange);
subplot(3,1,3);
imagesc(rc.x(:,1),1:size(rc.lfp,4),aresp-bresp);
%caxis(brange);
else
subplot(2,1,1);
imagesc(rc.x(:,1),1:size(rc.lfp,4),cat(2, aresp, blr(1,:)'));
text(rc.x(end,1),0,'blank');
caxis(arange);
subplot(2,1,2);
imagesc(rc.x(:,1),1:size(rc.lfp,4),cat(2, bresp, blr(1,:)'));
caxis(brange);
end

function FindSdfCorr(rc)
%look for inverted responses (AC, Luminance specific RF)

for p = 1:size(rc.tresps,4);
    subplot(6,4,p);
    hold off;
    [a,b] = max(max(rc.tresps(:,:,1,p)));
    [c,d] = max(max(rc.tresps(:,:,2,p)));
    x = [rc.tresps(:,b,1,p) rc.tresps(:,d,1,p)];
    y = [rc.tresps(:,b,2,p) rc.tresps(:,d,2,p)];
    scatter(x(:,1),y(:,1));
    hold on;
    scatter(x(:,2),y(:,2),'r');
    r = corrcoef(x(:),y(:));
    title(sprintf('Corr %.2f',r(1,2)));
end

function [latencies, details] = SDFlatency(rc, sumy, varargin)

tid = [];
details = [];
j = 1;
separatey = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'times',5)
        j = j+1;
        times = varargin{j}
        tid = find(rc.times >= times(1) & rc.times <= times(2));
    end
    j = j+1;
end


if sumy
    resp = mean(rc.tresps,3);
elseif separatey
    resp = rc.tresps;
else
    resp = reshape(rc.tresps,[size(rc.tresps,1) size(rc.tresps,2).*size(rc.tresps,3) size(rc.tresps,4)]);
end
if isempty(tid)
    tid = 1:length(rc.times);
end
stimvar = reshape(var(resp,[],2),[size(resp,1) size(resp,3) size(resp,4)]);
details.stimvar = stimvar;
pretid = find(rc.times > 0 & rc.times < 400);
for k = 1:size(stimvar,2)
[maxr(:,k), peakt(:,k)] = max(stimvar(tid,k,:));
base(:,k) = mean(stimvar(pretid,k,:));
basesd(:,k) = std(stimvar(pretid,k,:));
for j = 1:size(stimvar,3);
    latencies(j,k) = NaN;
    halfmax = base(j,k)/2 + maxr(j,k)/2;
    id = find(stimvar(1:peakt(j,k),k,j) < halfmax);
    if ~isempty(id)
    halft = id(end);  %time at halfmax
    id = find(stimvar(1:halft,k,j) < base(j,k)+ 2* basesd(j,k));
    ts = id(end);  %time when first exceed 2D of baseline
    if ts < halft
    lr = polyfit(rc.times(ts:halft),stimvar(ts:halft,1,j)'-base(j,k),1);
    latencies(j,k) = -lr(2)/lr(1);
    end
    end
end
end
details.alllatencies = latencies;
if size(latencies,2) == 2
id = find(abs(diff(latencies,[],2)) < 100);
else
    id = 1:size(latencies,1);
end
[latencies, mid] = min(latencies,[],2);
for j = 1:length(latencies)
    details.maxr(j) = maxr(j,mid(j));
    details.basesd(j) = basesd(j,mid(j));
    details.peakt(j) = peakt(j,mid(j));
    details.type = mid(j);
end
for j = 1:length(latencies)
    nid = find(~isnan(details.alllatencies(j,:)));
latencies(j) = mean(details.alllatencies(j,nid),2);
end
if length(id)
details.peakt(id) = mean(peakt(id,:),2);
details.maxr(id) = mean(maxr(id,:),2);
details.basesd(id) = mean(basesd(id,:),2);
detials.type(id) = 0;
end
details.base = mean(base,2)';
details.maxsdr = (details.maxr-details.base)./details.basesd;
details.allsdr = (maxr-base)./basesd;
    
function CheckSigns(LFP, probe)

uid = find([LFP.Trials.Dc] == 0 & [LFP.Trials.RespDir] > 0);
did = find([LFP.Trials.Dc] == 0 & [LFP.Trials.RespDir] < 0);
a = abs(cat(3,LFP.Trials(uid(1:end)).FTlfp));
upchoicepwr = mean(a(1:250,7,:),3);
a = abs(cat(3,LFP.Trials(did(1:end)).FTlfp));
dnchoicepwr = mean(a(1:250,7,:),3);

mo = mean([LFP.Trails.ori]);
uid = find([LFP.Trials.Dc] > 0.1 & [LFP.Trials.ori] > mo);
did = find([LFP.Trials.Dc] > 0.1 & [LFP.Trials.ori] < mo);
a = abs(cat(3,LFP.Trials(uid(1:end)).FTlfp));
upsigpwr = mean(a(1:250,7,:),3);
a = abs(cat(3,LFP.Trials(did(1:end)).FTlfp));
dnsigpwr = mean(a(1:250,7,:),3);
plot(upsigpwr-dnsigpwr);
hold on;
plot(upchoicepwr-dnchoicepwr,'r');




