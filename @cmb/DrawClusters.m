function DATA = DrawClusters(DATA, cluster, setfig)
if setfig
    set(0,'CurrentFigure',DATA.xyfig);
end

if ~iscell(cluster)
    return;
end
p = DATA.probe;
for j = 1:min([size(cluster,1) 7])
    if p > size(cluster,2)
        C = [];
    else
        C = cluster{j,p};
    end
    while ~isempty(C)
        if isfield(C,'params')
            color = DATA.spkcolor{j+1};
            if isfield(C,'touched')
                if C.touched == 0
                    color = color./2;
                elseif C.touched == -1 %online cluster
                    color = (1- color)/2;
                end
            end
            if C.params(1) == DATA.plot.clusterX & ...
                    C.params(2) == DATA.plot.clusterY
                h = cmb.DrawCluster(C,color);
                if myhandle(h)
                    DATA.cluster{j,p}.h = h;
                    if isfield(C,'autocut') & C.autocut
                        set(h,'LineStyle','--');
                    end
                end
                hold on;
            end
        end
        if isfield(C,'Cluster')
            C = C.Cluster;
        else
            C ={};
        end
    end
end

