%List of Matlab functions/structure for Clustering
%AllVPcs  Reads FullV files, makes ExptNClusterTimes
%PlotClusters Plots results in ExptNClusterTimes fiels
%ListClusterBackup  Looks at ExptNClusterTimes files in backup
%LoadCluster  Load Cluster +/- Details given folder, expt #
%ListTree Searches Listing of backup disks
%
%Fields in Cluster
% mahal(1) 2D mahal distance
% mahal(2) dprime of 2D mahal distance
% mahal(3) ND mahal distance
% mahal(4) 1D mahal distance
% fitdprime(1) is dprime of Gaussians Fit to each half of decision
%                  histogram. Negative values are correct direction
%          (2)  1/0 whether fitted mean 1 is right side of boundary
%          (3)  1/0 whether fitted mean 2 is right side of boundary
%          (4) bin width for histogram
%          (5) location of the minimum in the fit

