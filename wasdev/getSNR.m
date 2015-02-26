function stats = getSNR(Spikes,varargin)

% calculates the mean, std, and SNR for a 
% cluster, given a Spikes file

%SNR defined as in Kelly et al (J. Neurosci., 2007).

j=1;
cluster=0;
while j<=length(varargin)
    if strncmpi(varargin{j},'cluster',2)
        j=j+1;
        cluster=varargin{j};
    else
    end
    j=j+1;
end

if cluster
    clusters=cluster;
else
    clusters=1:6;
end

for c=clusters
    cind = find(Spikes.codes==c);
    if isempty(cind)
        continue
    end
    W{c}=Spikes.values(cind,:);
    W_hat{c}=mean(W{c},1);
    for i=1:size(W,1)
        e(i,:)=W{c}(i,:)-W_hat{c};
    end
    SNR(c)=( max(W_hat{c}) - min(W_hat{c}) ) / (2.*std(e(:)));
    stats.ms{c}=W_hat{c};
    stats.sd{c}=std(W{c},1);
end
stats.SNR=SNR;
