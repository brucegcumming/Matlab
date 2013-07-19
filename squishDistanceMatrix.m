function [squishedMatrix, score, cellMap] = squishDistanceMatrix( originalDistanceMatrix )
%squishDistanceMatrix rearranges the input distance matrix to maximize the 
% sum of the product of distances of cells to diagonal and their values 
%   this will take a very long time for matrices larger than 10 by 10
% 
% output:
% squishedMatrix is the rearranged matrix with maximum score
%
% score is the score (sum of product of distances to diagonal and cell
% values)
%
% cellMap is vecor of length n (where the distance matrix is n by n) and
% contains the mapping of each row (or column) in the original distance
% matrix to the squished matrix

p = perms(1:size(originalDistanceMatrix,1));
scores = [];

for ip = 1:size(p,1)
    remapedDM = remapDistanceMatrix(originalDistanceMatrix, p(ip,:));
    scores(ip) = squeezeScore(remapedDM);
end

permutations = p;
maxidx = find(scores==max(scores), 1);
cellMap = permutations(maxidx, :);
squishedMatrix = remapDistanceMatrix(originalDistanceMatrix, cellMap);
score(1) = scores(maxidx);
score(2) = max(squishedMatrix(sub2ind([6 6],[2 3 4 5 6],[1 2 3 4 5])));
end

function [remapedDM] = remapDistanceMatrix(dm, map)

remapedDM = zeros(size(dm));
for i = 1:size(dm,1)
    for j  = 1:size(dm,2)
        remapedDM(i,j) = dm(map(i), map(j));
    end
end

end

function [ score ] = squeezeScore( d )
%squeezeScore gets a square matrix, which is a distance matrix
%   and calculates sum of the product of values and distances from the
%   diagonal

s = 0;
rawsum = 0;
for i = 1:size(d,1)
    for j = 1:size(d,1)
      s = s + d(i,j) * abs(i-j);  
      rawsum = rawsum + d(i, j);
    end
end

score = s / rawsum;
end