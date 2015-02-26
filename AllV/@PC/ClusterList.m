function cells = ClusterList(DATA, Clusters, varargin)
%cells = ClusterList(Clusters) build a list of all clusters cut 
nc = 1;
for k = 1:length(Clusters)
    if DATA.plot.dprimemin == 0 || (Clusters{k}.fitdprime(1) < DATA.plot.dprimemin)
        cells(nc).p = k;
        cells(nc).cl = 1;
        cells(nc).dp = PC.DistanceMeasure(Clusters{k},1,DATA.mahaltype,'zeropad');
        cells(nc).dropi = Clusters{k}.dropi(3);
        if isfield(Clusters{k},'clst')
            t = find(Clusters{k}.clst == 2);
            if isfield(Clusters{k},'t')
                Clusters{k}.t = Clusters{k}.t(t);
            else
                Clusters{k}.times = Clusters{k}.times(t);
            end
        end
        nc = nc+1;
    end
    for j = 1:length(Clusters{k}.next)
        if isfield(Clusters{k}.next{j},'times') && ...
                (DATA.plot.dprimemin == 0 || Clusters{k}.next{j}.fitdprime(1) < DATA.plot.dprimemin)
            cells(nc).p = k;
            cells(nc).cl = 1+j;
            cells(nc).dp = PC.DistanceMeasure(Clusters{k},1+j,DATA.mahaltype,'zeropad');
            cells(nc).dropi = Clusters{k}.next{j}.dropi(3);
            nc = nc+1;
        end
    end
end
