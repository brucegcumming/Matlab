function [DATA, dprime, details] = SetSpkCodes(DATA, expspks, probe, show, varargin)

onecluster = 0;
nexp = DATA.currentexpt;
     j = 1;
while j <= length(varargin)
    j = j+1;
end
if DATA.plot.synccluster & DATA.syncsign ~= 2
    onecluster = DATA.currentcluster;
end
dprime = 0;
details.nc = 0;
if isempty(expspks)
    return;
end


if ~isfield(DATA,'cluster')
 nclusters = 0;
 DATA.cluster = {};
end

clusterdefined = cmb.iscluster(DATA.cluster,1,DATA.probe);

if isfield(DATA,'AllClusters')
    if iscell(DATA.AllClusters)
        DATA.AllClusters{nexp}(probe).codes(expspks) = 0;
        if isfield(DATA.AllData.Spikes,'codes')
            DATA.AllData.Spikes.codes(expspks,1) = 0;
        end
        Cx = DATA.AllClusters{nexp}(probe).cx;
        Cy = DATA.AllClusters{nexp}(probe).cy;
    else
  DATA.AllClusters(probe).codes(expspks) = 0;
  Cx = DATA.AllClusters(probe).cx;
  Cy = DATA.AllClusters(probe).cy;
    end
elseif isfield(DATA,'AllSpikes')
    if onecluster
        id = find(DATA.AllSpikes{probe}.codes(expspks,2) == onecluster);
        DATA.AllSpikes{probe}.codes(expspks(id),2) = 0;
    else
        DATA.AllSpikes{probe}.codes(expspks,2) = 0;
    end
  if iscluster(DATA.cluster,1,probe) & sum(ismember(DATA.cluster{1,probe}.params, [DATA.plot.clusterX DATA.plot.clusterY]) < 2)
      if ~isfield(DATA.AllSpikes{probe},'dVdt')
      DATA.AllSpikes{probe}.dVdt = diff(DATA.AllSpikes{probe}.values,1,2);
      end
      Spikes = DATA.AllSpikes{probe};
      if isfield(DATA.AllSpikes{probe},'pcs')
      PCs = DATA.AllSpikes{probe}.pcs;
      else
      PCs = DATA.AllSpikes{probe}.codes;
      end
  end
  if isfield(DATA.AllSpikes{probe},'cx')
%used to check 
%     && DATA.plot.synccluster == 0
%but now Allspikes has cx. AND must use this - DATA.Spikes.cx  indices only
%match one probe
      Cx = DATA.AllSpikes{probe}.cx;
      Cy = DATA.AllSpikes{probe}.cy;
  elseif isfield(DATA,'Spikes')
      Cx = DATA.Spikes.cx;
      Cy = DATA.Spikes.cy;
  else
      return;
  end
elseif isfield(DATA,'Spikes')
  Cx = DATA.Spikes.cx;
  Cy = DATA.Spikes.cy;
  %if no cluster defined, don't wipe out codes - might have been set by
  %autocut
  if DATA.state.online == 0 || clusterdefined
    DATA.AllData.Spikes.codes(expspks,2) = 0;
  end
    Spikes = DATA.AllData.Spikes;
    if isempty(DATA.AllData.pcs) || length(DATA.AllData.pcs) < max(expspks)
        PCs = DATA.AllData.Spikes.codes;
    else
        PCs = DATA.AllData.pcs;
    end
else
    return; % No spikes
end

%
%really want to limit this to spike in the current scope;id = find(DATA.Spks.cluster == 0);
nclusters = size(DATA.cluster,1);
p = probe;
if nclusters >1
%    fprintf('%s: %d Clusters for Probe %d\n',DATA.explabels{DATA.currentexpt},nclusters,p);
end
sumplot = 0;
details.maxx = max(Cx(expspks));
details.maxy = max(Cy(expspks));
if onecluster
    cllist = onecluster;
elseif nclusters > 7 %includes and artifact. Do this last
    cllist = [nclusters:-1:1 8];
else
    cllist = nclusters:-1:1;
end

for cl = cllist
    cspks = expspks;
    nspk = 0;
    if p > size(DATA.cluster,2) %no cluster defined in combine for this yet. Stick with exisitng
        C = [];
    else
        C = DATA.cluster{cl,p};
    end
    if isfield(C,'Cluster') && ~isempty(C.Cluster) 
        splitlist = 1;
    else
        splitlist = 0;
    end
    if isfield(C,'forceid') & C.forceid
        clid = C.forceid;
    else
        clid = cl;
    end
    while ~isempty(C) & isfield(C,'x')
        if C.params(1) ~= DATA.plot.clusterX || C.params(2) ~= DATA.plot.clusterY
            if isfield(Spikes,'cache') %will need spike values again
                DATA = cmb.LoadSpikes(DATA, DATA.currentexpt(1), probe,'force');
                Spikes = DATA.AllData.Spikes;
            end
        end
        if ~isfield(C,'lastspk')  & ~isempty(DATA.spklist)
            C.lastspk = DATA.spklist(end);
        end
        if ~isfield(C,'firstspk') & ~isempty(DATA.spklist)
            C.firstspk = DATA.spklist(1);
        elseif isfield(C,'firstspk') & C.firstspk > 0 & splitlist
            cspks = expspks(find(expspks >= C.firstspk & expspks <= C.lastspk));
        end
    if C.params(1) == DATA.plot.clusterX  %%otherwise need to calc values again)     
        x = (Cx(cspks) - C.x(1))./C.x(3);
    else
%Jan 2015 what was here before can't work unless max(csspks) = length(cspks)        
%        x = GetSpikeVals(DATA, cspks, Spikes.values(cspks,:), Spikes.dVdt(cspks,:),C.params(1), 1,PCs(cspks,:));
        x = GetSpikeVals(DATA, 1:length(cspks), Spikes.values(cspks,:), Spikes.dVdt(cspks,:),C.params(1), 1,PCs(cspks,:));
        x = (x- C.x(1))./C.x(3);
    end
    if C.params(2) == DATA.plot.clusterY  %%otherwise need to calc values again)     
        y = (Cy(cspks) - C.y(1))./C.y(3);
    else
%        y = GetSpikeVals(DATA, cspks, Spikes.values(cspks,:), Spikes.dVdt(cspks,:),C.params(2), 1,PCs(cspks,:));        
        y = GetSpikeVals(DATA, 1:length(cspks), Spikes.values(cspks,:), Spikes.dVdt(cspks,:),C.params(2), 1,PCs(cspks,:));
        y = (y- C.y(1))./C.y(3);
    end
    xr = x .* cos(C.angle) + y .* sin(C.angle);
    yr = y .* cos(C.angle) - x .* sin(C.angle);
    d = (yr./C.y(2)*C.y(3)).^2 + (xr./C.x(2)*C.x(3)).^2;
    
    
%####################################    
%Here is where cluser is actually set
%###################################
    id = find(d <= 1);
    nid = find(d > 1);
    nspk = nspk + length(id);
%####################################    
%Here is where cluser is actually set
%###################################

    
    
    if length(id) & length(nid)
    dprime(cl) = (mean(d(nid))-mean(d(id)))./sqrt(mean([var(d(nid)) var(d(id))]));
    else
        dprime(cl) = NaN;
    end
    tic;
    if length(nid) > 10000
        pcrit = 1;
    elseif length(nid) > 1000
        pcrit = 10;
    elseif length(nid) > 200
        pcrit = 20;
    else
        pcrit = 50;
    end
    %    sy = ((yr-mean(yr(id)))./std(yr(id)));
%    sx = ((xr-mean(xr(id)))./std(xr(id)));   
    sy = ((yr-mean(yr(id)))./std(yr));
    sx = ((xr-mean(xr(id)))./std(xr));
    sd = abs(sx+i*sy);
    o = [0:pi/40:pi];
    for j = 1:length(o)
        d = sx .* cos(o(j)) + sy.* sin(o(j));
        a = prctile(abs(d(nid)),pcrit);
        tid = find(abs(d(nid))< a);
        ds{j} = d;
%
%        dprimes(j) = abs(mean(d(nid(tid)))-mean(d(id)))./sqrt(mean([var(d(nid(tid))) var(d(id))]));
         dprimes(j) = abs(mean(d(nid))-mean(d(id)))./sqrt(mean([var(d(nid)) var(d(id))]));
         bii(j) = (1+skewness(d).^2)./(kurtosis(d)+3);
%        dprimes(j) = abs(mean(d(nid))-mean(d(id)))./std(d(nid));
    end
    [dprime(cl), maxi] = max(dprimes);
    DATA.cluster{cl,p}.dprimepar = [o(maxi)];
    details.truedprime = abs(mean(d(nid))-mean(d(id)))./sqrt(mean([var(d(nid)) var(d(id))]));
    details.bii = max(bii);
    bii = max(bii);
    debug = 0;
    if debug
        [a,b] = smhist(ds{maxi}(id),'sd',0.1);
        [c,d] = smhist(ds{maxi},'sd',0.1);
        [e,f] = smhist(ds{maxi}(tid),'sd',0.1);
        of = gcf;
        GetFigure('DDF');
        hold off;
        plot(b,a);
        hold on;
        plot(d,c,'r');
        plot(f,e,'g');
        title(sprintf('Spks %d/%d/%d',length(ds{maxi}),length(id),length(tid)));
        figure(of);
    end
    details.nc(cl) = length(id);

    
    if isfield(DATA,'AllClusters')
        if iscell(DATA.AllClusters)
            DATA.AllClusters{nexp}(probe).codes(cspks(id)) = clid;
            DATA.AllData.Spikes.codes(cspks(id))=clid;
        else
            DATA.AllClusters(probe).codes(cspks(id)) = clid;
        end
    elseif isfield(DATA,'AllSpikes')
        DATA.AllSpikes{probe}.codes(cspks(id),2) = clid;
    else
        DATA.AllData.Spikes.codes(cspks(id),2) = clid;
    end
    DATA = cmb.SpkCache(DATA,DATA.currentexpt(1),DATA.probe,'set');
    DATA.cluster{cl,p}.dprime = dprime(cl);
    DATA.cluster{cl,p}.bii = max(bii);
    if isfield(C,'Cluster')
        C = C.Cluster;
    else
        C = {};
    end
    if show &&  (~DATA.densityplot && ~DATA.alldensityplot)
        DATA.ptsize = CheckPtSize(DATA, length(expspks));
        plot(Cx(cspks(id)),Cy(cspks(id)),'.',...
            'color',DATA.spkcolor{clid+1},'markersize',DATA.ptsize);
        sumplot = sumplot + length(id);
        if sumplot
            hold on;
        end
        ei = DATA.currentexpt;
        if isfield(DATA.Expts{ei}.Stimvals,'du')
            rate = length(id)./(length(DATA.Expts{ei}.Trials) .* DATA.Expts{ei}.Stimvals.du);
        else
            rate = length(id)./length(DATA.Expts{ei}.Trials);
        end
        if isfield(DATA,'explabels') & cl == DATA.currentcluster
    title(sprintf('%s dprime %.2f(%.2f) %.1fHz:%d/%dspks',DATA.explabels{DATA.currentexpt},dprime(cl),bii,rate,length(id),length(cspks))); 
        end
    else
        if isfield(DATA,'explabels') && length(DATA.explabels) >= DATA.currentexpt(1) && cl == DATA.currentcluster
    title(sprintf('%s dprime %.2f %d/%dspks',DATA.explabels{DATA.currentexpt(1)},dprime(cl),length(id),length(cspks))); 
        end
    end
    end
    DATA.cluster{cl,p}.nspk = nspk;
end %for cl = cllist
if DATA.state.somespikes ==2
    return;
end
if length(dprime) >= DATA.currentcluster
    dprime = dprime(DATA.currentcluster);
else
    dprime = dprime(1);
end
    if ~DATA.densityplot && show > 1 && ~DATA.alldensityplot %show unclassified also
    if isfield(DATA,'AllClusters')
        if iscell(DATA.AllClusters)
           id = DATA.AllClusters{nexp}(probe).codes(cspks) == 0;
           DATA.AllData.Spikes.codes(cspks(id))=0;
        else
            id = DATA.AllClusters(probe).codes(cspks) == 0;
        end
    elseif isfield(DATA,'AllSpikes')
        id = DATA.AllSpikes{probe}.codes(cspks,2) == 0;
    else
        id = DATA.AllData.Spikes.codes(cspks,2) == 0;
    end
        id = cspks((find(id)));
        plot(Cx(id),Cy(id),'.','color',DATA.spkcolor{1},'markersize',DATA.ptsize);
        hold on;
    end
