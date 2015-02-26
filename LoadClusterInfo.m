function Expt = LoadClusterInfo(Expt, varargin)
% Expt = LoadClusterInfo(Expt, varargin)  Read additional
%    information from Cluster Files for Expt, including
%  Mean Spike Waveform



loaddir = fileparts(Expt.Header.loadname);

if isfield(Expt.Header,'Combineids')
    for j = 1:length(Expt.Header.Combineids)
        p = Expt.Header.Clusters{j}.probe;
        name = sprintf('%s/Expt%dClusterTimes.mat',loaddir,Expt.Header.Combineids(j));
        if exist(name)
        load(name);
        Expt.Header.Clusters{j}.MeanSpike = Clusters{p}.MeanSpike;
        Expt.Header.Clusters{j}.trigt = find(Clusters{p}.spts ==0);
        end
    end
end