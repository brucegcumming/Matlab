function [xc, details] = xcorr(E, varargin)
%calculate cross correlation between spikes 
%if length(Expts) > 2, returns a matrix of
%cross correlations


nc = 0;
spkt = expt.spktimes(E);
for j = 2:length(spkt)
    for k = 1:(j-1)
        nc = nc+1;
        [xc(nc,:), details(nc)] = xcorrtimes(spkt{j},spkt{k});        
    end
end

