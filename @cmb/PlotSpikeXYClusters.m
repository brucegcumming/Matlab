function PlotSpikeXYClusters(DATA, expspks, clusters, varargin)
%do the xy plot for each defined ckuster.,
ctype = 2;
cmb.PlaySpikes(DATA, DATA.currentexpt(1), varargin{:}, 'xyonly');

probe = GetProbe(DATA,DATA.currentexpt(1), DATA.probe);
for j = clusters
    if isfield(DATA,'AllClusters')
        if iscell(DATA.AllClusters)
            ctype = 1;
            id = find(DATA.AllClusters{DATA.currentexpt(1)}(probe).codes(expspks,ctype) == j-1);
        else
            id = find(DATA.AllClusters(probe).codes(expspks,ctype) == j-1);
        end
    elseif isfield(DATA,'AllSpikes')
        id = find(DATA.AllSpikes{DATA.probe}.codes(expspks,ctype) == j-1);
    else
        id = find(DATA.AllData.Spikes.codes(expspks,ctype) == j-1);
    end
    if length(id)
        cmb.PlotSpikeXY(DATA, expspks(id), DATA.spkcolor{j});
        hold on;
    end
end
cmb.ShowCellLabels(DATA);
if DATA.plot.showartifacts
    cmb.PlotSpikeXYClusters(DATA,expspks,9);
end


