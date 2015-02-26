function nc = MaxSpikeCode(DATA, expspks)

ctype = 2;

if isfield(DATA,'AllClusters')
    probe = GetProbe(DATA,DATA.currentexpt(1), DATA.probe);
    if iscell(DATA.AllClusters)
        ctype = 1;
        if isempty(expspks)
            nc = [];
        else
            nc = max(DATA.AllClusters{DATA.currentexpt(1)}(probe).codes(expspks,ctype));
        end
    else
        nc = max(DATA.AllClusters(DATA.probe).codes(expspks,ctype));
    end
else
    nc = max(DATA.AllData.Spikes.codes(expspks,ctype));
end


