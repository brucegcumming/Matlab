function DATA = CalcClusterVars(DATA, ispk, varargin)
%DATA = CalcClusterVars(DATA, ispk, varargin)


SPKENERGY=1;
SPKVARE = 2;
noforce  = 1;
j = 1;
probe = DATA.probe;
while j <= length(varargin)
    if isstruct(varargin{j}) & isfield(varargin{j},'values')
        j = j+1;
        Spikes = varargin{j};
    elseif strncmpi(varargin{j},'force',5)
        noforce = 0;
    elseif strncmpi(varargin{j},'probe',5)
        j =j+1;
        probe = varargin{j};
    end
    j = j+1;
end

if isfield(DATA,'AllClusters') && noforce % All calculated and stored, but No ADCs, so
    DATA.Spikes.cx = DATA.AllClusters{DATA.currentexpt}(DATA.probe).cx;
    DATA.Spikes.cy = DATA.AllClusters{DATA.currentexpt}(DATA.probe).cy;
  return;
end

if length(ispk) == 2
    ispk = ispk(1):ispk(2);
end
if ispk
    DATA = CheckForPCA(DATA,ispk, 1);
    if isfield(DATA,'AllSpikes')
        Spikes = DATA.AllSpikes{probe};
        if isfield(DATA.AllSpikes{probe},'pcs')
        PCs = DATA.AllSpikes{probe}.pcs;
        else
        PCs = DATA.AllSpikes{probe}.codes;
        end
        adc = Spikes.values(ispk,:);
        if ~isfield(Spikes,'dVdt')
            Spikes.dVdt(ispk,:) = diff(adc,[],2);
        end
    else
        Spikes = DATA.AllData.Spikes;
        adc = Spikes.values(ispk,:);
        if isempty(DATA.AllData.pcs) || length(DATA.AllData.pcs) < max(ispk) %might have pcs computed for earlier spikes
            PCs = DATA.AllData.Spikes.codes;
        else
        PCs = DATA.AllData.pcs;
        end
    end
    adc = double(Spikes.values(ispk,:));
    if isfield(Spikes,'dVdt')
        energy  = sum(Spikes.dVdt(ispk,:)'.^2);
    else
    energy  = sum(diff(adc').^2);
    end
    svar = var(adc');
    DATA.Spikes.energy(ispk)= energy;
    if DATA.plot.clusterX == SPKENERGY && noforce
        DATA.Spikes.cx(ispk)= energy;
    else
        DATA.Spikes.cx(ispk)= GetSpikeVals(DATA, ispk, Spikes.values(ispk,:), Spikes.dVdt(ispk,:),DATA.plot.clusterX, 1,PCs(ispk,:));
    end
    DATA.Spikes.vw(ispk) = svar./energy;
    if DATA.plot.clusterY == SPKVARE & noforce
        DATA.Spikes.cy(ispk)= DATA.Spikes.vw(ispk);
    else
        DATA.Spikes.cy(ispk)= GetSpikeVals(DATA, ispk, Spikes.values(ispk,:), Spikes.dVdt(ispk,:),DATA.plot.clusterY, 1,PCs(ispk,:));
    end
    if isfield(DATA,'AllSpikes')
        lastspk = min([length(DATA.AllSpikes{probe}.times) length(DATA.Spikes.cx)]);
        DATA.AllSpikes{probe}.cx = DATA.Spikes.cx(1:lastspk);
        DATA.AllSpikes{probe}.cy = DATA.Spikes.cy(1:lastspk);
    end
end
