function MergeCluster(Clusts,source,destination)

%takes a matrix of exptno/probe row vectors and copies the relevant
%ClusterTimes and Details chunks from one directory to another
% 'Clusters' input should be of the form:
%  [exptno1 probe1
%   exptno2 probe2
%     ...    ...
%   exptnoN probeN];

tic; 
[a b] = size(Clusts);
for i=1:a
    exptno=Clusts(i,1);
    probe=Clusts(i,2);
    SourceCluster=load([source,'/Expt',num2str(exptno),'ClusterTimes.mat']);
    SourceClusterDetails=load([source,'/Expt',num2str(exptno),'ClusterTimesDetails.mat']);
    load([destination,'/Expt',num2str(exptno),'ClusterTimes.mat']);
    load([destination,'/Expt',num2str(exptno),'ClusterTimesDetails.mat']);
    Clusters{probe}=SourceCluster.Clusters{probe};
    ClusterDetails{probe}=SourceClusterDetails.ClusterDetails{probe};
    save([destination,'/Expt',num2str(exptno),'ClusterTimes.mat'],'Clusters');
    save([destination,'/Expt',num2str(exptno),'ClusterTimesDetails.mat'],'ClusterDetails');
end
F=toc;
fprintf('\n Took %s seconds to copy %s chunks. \n',num2str(round(F)),num2str(a));

