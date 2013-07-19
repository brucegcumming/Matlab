function set = ClusterIsSet(Expt, probe)

if isfield(Expt,'Cluster') && size(Expt.Cluster,2) >= probe && ...
        ~isempty(Expt.Cluster{1,probe}) && isfield(Expt.Cluster{1,probe},'x')
    set = 1;
    if isfield(Expt.Cluster{1,probe},'autocut') & Expt.Cluster{1,probe}.autocut > 0
        set = 2;
    end
else
    set = 0;
end
