function [Lratio,iso_distance] = CalcLratio(X,gmm_fit,cluster_labels, varargin)

if ~isobject(gmm_fit)
    Lratio = nan;
    iso_distance = nan;
    return;
end

N_sus = nanmax(cluster_labels)-1;
[N_spks,df] = size(X);
if N_sus == 0
    Lratio = nan;
    iso_distance = nan;
    return;
end

Lratio = nan(N_sus,1);
iso_distance = nan(N_sus,1);
for ii = 1:N_sus
    cluster_comps = ii+1;
    clust_spikes = find(cluster_labels == ii+1);
    non_clust_spikes = setdiff(1:N_spks,clust_spikes);
    
    %if single-component cluster, use the model stats for that comp
    if length(cluster_comps) == 1
        clust_mean = gmm_fit.mu(cluster_comps,:);
        clust_Sigma = squeeze(gmm_fit.Sigma(:,:,cluster_comps));
    else %use empirical mean if more than one gaussians used for cluster
        clust_mean = mean(X(clust_spikes,:));
        clust_Sigma = cov(X(clust_spikes,:));
    end
    mean_seps = bsxfun(@minus,X,clust_mean);
    D = sum((mean_seps/clust_Sigma) .* mean_seps,2); %mah D
    
    Lratio(ii) = sum(1 - chi2cdf(D(non_clust_spikes),df));
    
    
    if isempty(clust_spikes) || isempty(non_clust_spikes)
        iso_distance(ii) = nan;
    else
        D = sqrt(sort(D(non_clust_spikes)));
        if length(clust_spikes) <= length(non_clust_spikes)
            iso_distance(ii) = D(length(clust_spikes));
        else
            iso_distance(ii) = D(end);
        end
    end
end