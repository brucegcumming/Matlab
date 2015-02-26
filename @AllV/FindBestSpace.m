function [bestspace, DATA] = FindBestSpace(DATA, C, varargin)
%Given a Cluster classification, find an ND space that optimizes the
%separation measure using mahal dprime
% by default ccompares ! vs 0
%
cl{1} = DATA.currentcluster+1;
cl{2} = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'clusters',5)
        j = j+1;
        cl{1} = varargin{j};
        j = j+1;
        cl{2} = varargin{j};        
    end
    j = j+1;
end

if DATA.usegmcid
    nid = find(ismember(DATA.gmcid,cl{1}));
    gid = find(ismember(DATA.gmcid,cl{2}));
else
    nid = find(ismember(DATA.clst,cl{1}));
    gid = find(ismember(DATA.clst,cl{2}));
end
    pts = [nid(:)' gid(:)'];
    clst = DATA.clst;
    clst(nid) = 1;
    clst(gid) = 2;

bestspace.space = C.space;    
method =2;
T = [];
if C.space(1) == 3
    T = DATA.TemplateScores;
elseif C.space(1) == 1
    T = DATA.pcs;
    if size(T,2) > 10
        T = T(:,1:6);
    end
elseif C.space == 2
    T = [];
else
    T = [];
end
if ~isempty(T)
if method ==1 %just look at this cluster vs everyone else
    dims = C.space(2:3);
    newdims = setdiff(1:size(T,2),dims);
    d(1,:) = CalcIsolation(T(pts,dims),clst(pts),2);
    for j = 1:length(newdims)        
        d(j+1,:) = CalcIsolation(T(pts,[dims newdims(j)]),clst(pts),2);
    end
    max(d);
elseif method ==2 %try all pairings
    n = 1;
    for j = 2:size(T,2)
        for k = 1:j-1;
            score(n,:) = CalcIsolation(T(pts,[j k]),clst(pts),2);
            space(n,1) = j;
            space(n,2) = k;
            n = n+1;
        end
    end
    [a,b] = max(score);
    dims = unique(space(b,:))';
    if length(dims) > 2
        d = CalcIsolation(T(pts,dims),clst(pts),2);
        if sum(d-a) < 0.5  %no real improvement
            dims = space(b(1),:);
        else
            a = d;
        end            
    end
    newdims = setdiff(1:size(T,2),dims);
    ds(1,:) = CalcIsolation(T(pts,dims),clst(pts),2);
    for j = 1:length(newdims)        
        ds(j+1,:) = CalcIsolation(T(pts,[dims newdims(j)]),clst(pts),2);
    end
    [c,d] = max(ds);
    if sum(c-a) > 0.5 %3D space really better
        if c(1) - a(1) > 0.1 %in case c() ~ a(1) and c(2) >> a(2) 
            if d(1) > 1
                dims = [dims newdims(d(1)-1)];
            else
                cprintf('red','BestSpace: 3D > 2D but no third dimension used!\n');
            end
            bestspace.isolation = ds(d(1),:);
        else
            dims = [dims newdims(d(2)-1)];
            bestspace.isolation = ds(d(2),:);
        end
    else
        bestspace.isolation = ds(1,:);
        if sum(d > 1)  %a 3d space was a bit better
            dm =[ 0 0 ];
            id = find(d > 1);
            dm(id) = newdims(d(id)-1);
            fprintf('Best 3D space isolation %.3f,%.3f including space %d,%d\n',ds(d(1),1),ds(d(2),2),dm(1),dm(2));
        end
    end    
end
bestspace.space = [C.space(1) dims];
fprintf('Best Mahal distance cl%d vs cl%d  %.3f,%.3f in space%s\n',cl{1}-1,cl{2}-1,bestspace.isolation,sprintf(' %d', bestspace.space));
    
end





