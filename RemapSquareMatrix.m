function [remapedDM] = remapSquareMatrix(dm, map)
% remapSquareMatrix(dm, map)  re-order a square matix using order defined in map
%if length(map) < size(dm,1), will return a samaller matrix
remapedDM = zeros(length(map));
for i = 1:length(map)
    for j  = 1:length(map)
        remapedDM(i,j) = dm(map(i), map(j));
    end
end

end
