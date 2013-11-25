function ClusterParam(C)
%Cluster Structures created from FullV Files
%
% dropi(3) is the distance (in SDs) from the Trigger level to the fitter
%          mean of spike trigger points. so > 2 = >95% of spikes caught
%
% mahal(1) mahalanobis distance for 2-D decision space
% mahal(4) mahalanobis distance for 1-D projection of criterion value
%
% fitdprime(1) dprime distance between Gaussian fit only to data either side
%            of criterion (negative)