function [d, details] = GMdprime(G, varargin)
%calcualte drpime between two Gaussians in gmdistribution fit        
%calcualte drpime between two Gaussians in gmdistribution fit        

if ~isobject(G)
        d = NaN;
        details.d = [NaN NaN NaN];
        return;
    end
    nc =size(G.mu,1);
if size(G.mu,1) == 3
    %find three distances. one is allowed to be smaltt (splitting hash into
    %two. But one cluster must be distant from both of these. So take
    %lowest of top two = middle value
distance = mahal(G,G.mu);
    d(1) = sqrt(2./((1./distance(1,2))+(1./distance(2,1))));
    d(2) = sqrt(2./((1./distance(1,3))+(1./distance(3,1))));
    d(3) = sqrt(2./((1./distance(2,3))+(1./distance(3,2))));
    details.d = d;
    ds = sort(d);
    d  = d(2);  
elseif nc > 3
    distance = mahal(G,G.mu);
    for j = 1:nc
        for k = 1:j-1
            d(j,k) = sqrt(2./((1./distance(j,k))+(1./distance(k,j))));
            d(k,j) = d(j,k);
        end
        d(j,j) = NaN;
    end
%need to return a scalar for typical use.
%put full matix in details
    details.d = d;
    d  = mean(d(d>0));  
%    d = sqrt(2./((1./distance(2,1))+(1./distance(1,2))));
else
    distance = mahal(G,G.mu);
    d = sqrt(2./((1./distance(2,1))+(1./distance(1,2))));
    details.d = d;
end
