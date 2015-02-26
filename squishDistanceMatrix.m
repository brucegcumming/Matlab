function [squishedMatrix, score, cellMap, lapscores] = squishDistanceMatrix( originalDistanceMatrix, varargin )
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

findmin = 0;
nrnd = 100000;
method = 6;
laps = 10;
trackplot = 0;
markwait = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'min',3)
        findmin = 1;
    elseif strncmpi(varargin{j},'nrnd',3)
        j = j+1;
        nrnd = varargin{j};
    elseif strncmpi(varargin{j},'laps',3)
        j = j+1;
        laps = varargin{j};
    elseif strncmpi(varargin{j},'track',3)
        trackplot = 1;
    elseif strncmpi(varargin{j},'waitbar',3)
        markwait = 1;
        trackplot = 1;
    end
    j = j+1;
end
scores = [];
score(3) = NaN;
len = size(originalDistanceMatrix,1);
dd = abs(bsxfun(@minus,[1:len],[1:len]'));
if markwait
    markwait = waitbar(0, 'Finding Matrix');
end
if len < 9
    p = perms(1:len);
    for ip = 1:size(p,1)
        remapedDM = remapDistanceMatrix(originalDistanceMatrix, p(ip,:));
        scores(ip) = squeezeScore(remapedDM,dd);
    end
else
    p = zeros(nrnd,len);
    if method == 0
        for ip = 1:nrnd
            p(id,:) = randperm(len);
            remapedDM = remapDistanceMatrix(originalDistanceMatrix, p(ip,:));
            scores(ip) = squeezeScore(remapedDM,dd);
        end
    elseif method ==1 
    oldorder = 1:len;
    scores(1) = squeezeScore(originalDistanceMatrix,dd);
    bestscore = scores(1);
    p(1,:) = oldorder;
    keep(1) = 0;
    lastswap = 2;
    for ip = 2:nrnd
        if ip - length(keep) > 100
            oldorder = randperm(len);
%            if swplen < len && lastswap == swplen;
%                swplen = swplen+1;
%            end
            swplen = 0;
        else
            swplen = 2;
        end
        swplens(ip) = swplen;
        x = randperm(len,swplen);
        neworder = oldorder;
        neworder(x) = oldorder(circshift(x,[0 1]));
        remapedDM = remapDistanceMatrix(originalDistanceMatrix, neworder);
        scores(ip) = squeezeScore(remapedDM,dd);
        p(ip,:) = neworder;
        if swplen == 0 %new permutation
            bestscore = scores(ip);
        end
        if scores(ip) <= bestscore
            oldorder = neworder;
            keep(ip) = swplen;
            bestscore = scores(ip);
            score(3) = ip;
            lastswap = swplen;
        end
    end
    elseif method ==3
        oldorder = 1:len;
        scores(1) = squeezeScore(originalDistanceMatrix,dd);
        bestscore = scores(1);
        p(1,:) = oldorder;
        keep(1) = 0;
        lastswap = 2;
        stop = 0;
        ip = 0;
        swaps = nchoosek(1:len,4);
        while stop == 0
            j = 1;
            while j < size(swaps,1)
                x = swaps(j,:);
                neworder = oldorder;
                neworder(x) = oldorder(circshift(x,[0 1]));
                remapedDM = remapDistanceMatrix(originalDistanceMatrix, neworder);
                ip = ip+1;
                scores(ip) = squeezeScore(remapedDM,dd);
                p(ip,:) = neworder;
                if scores(ip) < bestscore %when score improves, check all permutations again
                    oldorder = neworder;
                    keep(ip) = length(x);
                    bestscore = scores(ip);
                    score(3) = ip;
                    lastswap = keep(ip);
                    j = NaN;
                else
                    j = j+1;
                end
            end
            if j == size(swaps,1) || ip > nrnd
                if size(swaps,2) > 2
                swaps = nchoosek(1:len,2);
                else
                stop = 1;
                end
            end
        end
        
    elseif method ==4 %just find best across nchoose(len,k).
        %no good. need perms for each choose to make this comprehensive
        oldorder = 1:len;
        scores(1) = squeezeScore(originalDistanceMatrix,dd);
        bestscore = scores(1);
        p(1,:) = oldorder;
        keep(1) = 0;
        lastswap = 2;
        swaps = nchoosek(1:len,5);
        stop = 0;
        ip = 0;
        p = zeros(size(swaps,1),len);
        scores = zeros(1,size(swaps,1));
        for ip = 1:size(swaps,1)
            x = swaps(ip,:);
            neworder = oldorder;
            neworder(x) = oldorder(circshift(x,[0 1]));
            remapedDM = remapDistanceMatrix(originalDistanceMatrix, neworder);
            scores(ip) = squeezeScore(remapedDM,dd);
            p(ip,:) = neworder;
        end
    elseif method == 5 %just find best across prems(nchoose(len,k))
        oldorder = 1:len;
        scores(1) = squeezeScore(originalDistanceMatrix,dd);
        bestscore = scores(1);
        p(1,:) = oldorder;
        keep(1) = 0;
        lastswap = 2;
        swaps = nchoosek(1:len,3);
        subs = perms(1:3);
        stop = 0;
        ip = 0;
        p = zeros(size(swaps,1).*size(subs,1),len);
        scores = zeros(1,size(swaps,1).*size(subs,1));
        ns = 1;
        for ip = 1:size(swaps,1)
            for j = 1:size(subs,1)
                x = swaps(ip,subs(j,:));
                neworder = oldorder;
                neworder(x) = oldorder(circshift(x,[0 1]));
                remapedDM = remapDistanceMatrix(originalDistanceMatrix, neworder);
                scores(ns) = squeezeScore(remapedDM,dd);
                p(ns,:) = neworder;
                ns = ns+1;
            end
        end
    elseif method ==6
%for well conditioned distance matrices composed N samples from M positions,
%so that the asnwere consists of squares along the diagnonal.
%pairwise swapping will find the correct answer.
        oldorder = 1:len;
        rng(101);
%        oldorder = randperm(len);
        scores(1) = squeezeScore(originalDistanceMatrix,dd);
        bestscore = scores(1);
        p(1,:) = oldorder;
        keep(1) = 0;
        lastswap = 2;
        stop = 0;
        ip = 0;
        swaps = nchoosek(1:len,2);
        n = size(swaps,1) * 2 * laps;
%        p = zeros(n,len);
%        scores = zeros(1,n);
        ts = now;
        while stop == 0
            j = 1;
            while j < size(swaps,1)
                x = swaps(j,:);
                neworder = oldorder;
                neworder(x) = oldorder(circshift(x,[0 1]));
                remapedDM = remapDistanceMatrix(originalDistanceMatrix, neworder);
                d = squeezeScore(remapedDM,dd);
                if d < bestscore %when score improves, check all permutations again
                    ip = ip+1;
                    scores(ip) = squeezeScore(remapedDM,dd);
                    p(ip,:) = neworder;
                    oldorder = neworder;
                    keep(ip) = length(x);
                    bestscore = scores(ip);
                    score(3) = ip;
                    lastswap = ip;
                    if trackplot && ip > trackplot+100
                        if markwait == 0
                            plot(scores);
                            drawnow;
                            fprintf('%d(%d/%d):%.4f (%.2f sec)\n',ip, j,size(swaps,1),scores(ip),mytoc(ts));
                        end
                        trackplot = ip;
                        waitbar(j./size(swaps,1));
                    end
                    j = NaN;
                else
                    j = j+1;
                end
            end
            if j == size(swaps,1)
                ip = ip+1;
               keep(ip) = lastswap;
                if laps > 1
                    laps = laps -1;
                    oldorder = randperm(len);
                    remapedDM = remapDistanceMatrix(originalDistanceMatrix, oldorder);
                    bestscore = squeezeScore(remapedDM,dd);
                else
                    stop = 1;
                end
            end
        end
        lapscores = scores(keep(keep > 2));
    end
end

permutations = p;
if findmin
    [a, maxidx] = min(scores);
else
    [a, maxidx] = max(scores);
end
cellMap = permutations(maxidx, :);
squishedMatrix = remapDistanceMatrix(originalDistanceMatrix, cellMap);
score(1) = scores(maxidx);
if len > 5
score(2) = max(squishedMatrix(sub2ind([6 6],[2 3 4 5 6],[1 2 3 4 5])));
end

if markwait
    delete(markwait);
end


function [remapedDM] = remapDistanceMatrix(dm, map)

remapedDM = zeros(size(dm));
for i = 1:size(dm,1)
    for j  = 1:size(dm,2)
        remapedDM(i,j) = dm(map(i), map(j));
    end
end


function [ score ] = squeezeScore( d, dd )
%squeezeScore gets a square matrix, which is a distance matrix
%   and calculates sum of the product of values and distances from the
%   diagonal

%s = dd .*d;
%score = sum(s(:));
%return;

s = 0;
rawsum = 0;
for i = 1:size(d,1)
    for j = 1:size(d,1)
      s = s + d(i,j) * abs(i-j);  
      rawsum = rawsum + d(i, j);
    end
end

score = s / rawsum;
