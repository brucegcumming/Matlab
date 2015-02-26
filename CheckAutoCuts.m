function CheckAutoCuts(X,varargin)

j = 1;
while j <= length(varargin)
    j = j+1;
end

if isfield(X,'cuts')
    for j = 1:length(X.cuts)
        nsu(j) = sum(X.cuts{j}.cluster.eckercluster.SU);
    end
end