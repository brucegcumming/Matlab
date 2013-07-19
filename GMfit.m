function [G, D, all] = GMfit(X, nd, nr, varargin)
%[G, D, all] = GMfit(X, nd, nr, ..)
%Fit Gaussian mixture model wiht nd componentc to values in matrix X
% X is an NxD matrix
% if nr >1, nr different fits are attempted.
% G returns the fit with teh best separation between components
%

idlist = [];
    if nr < 1
        nr = 1;
    end
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'idlist',5)
            j =j+1;
            idlist = varargin{j};
        end
        j = j+1;
    end
    ts = now;
%include my own starting point
      C = cov(X);
      pvar = diag(C);
        [E,V] = eig(C);
        pc = E(:,end);
        pcb = E(:,end-1);
        pc = pc./max(abs(pc));
        pcb = pcb./max(abs(pcb));
        for j = 1:size(X,2)
            S.mu(1,j) = mean(X(:,j)) + pc(j);
            S.mu(2,j) = mean(X(:,j)) - pc(j);
            if nd == 3
                S.mu(3,j) = mean(X(:,j)) + pcb(j);
            end
        end
        for j = 1:nd
        S.Sigma(:,:,j) = C./sqrt(2);
        end
        all.guess = S;

    
    for j = 1:nr
        try
       all.Gs{j} = gmdistribution.fit(X,nd,'Options',statset('MaxIter',1000));
       [all.ds(j), all.alld{j}] = gmdprime(all.Gs{j});
        catch
            fprintf('GM Fit fail\n');
            all.Gs{j} = S;
            all.Gs{j}.Converged = -1;
            all.Gs{j}.NlogL = NaN;
        end
    end
    j = j+1;
    try
    all.Gs{j} = gmdistribution.fit(X,nd,'Options',statset('MaxIter',1000),'Start',S);
    [all.ds(j), all.alld{j}] = gmdprime(all.Gs{j});
    catch
        fprintf('GM Fit (My Start) fail\n');
        all.Gs{j} = S;
        all.Gs{j}.Converged = -1;
        all.Gs{j}.NlogL = NaN;
        all.ds(j) = 0;
    end
    if length(idlist) == length(X)
        j = j+1;
        try
            id = find(idlist == 2);
            nid = find(idlist==1);
            S.mu(1,:) = mean(X(nid,:));
            S.mu(2,:) = mean(X(id,:));
            S.Sigma(:,:,1) = cov(X(nid,:));
            S.Sigma(:,:,2) = cov(X(id,:));
            S.PComponents(1) = length(nid)./size(X,1);
            S.PComponents(2) = length(id)./size(X,1);
            all.Gs{j} = gmdistribution.fit(X,nd,'Options',statset('MaxIter',1000),'Start',S);
            [all.ds(j), all.alld{j}] = gmdprime(all.Gs{j});
        catch
            fprintf('GM Fit (Start with Classification) fail\n');
            all.Gs{j} = S;
            all.Gs{j}.Converged = -1;
            all.Gs{j}.NlogL = NaN;
            all.ds(j) = 0;
        end
    end
    if all.ds(j) > max(all.ds(1:j-1)) * 1.1
        [a,b] = max(all.ds(1:j-1));
           fprintf('Best manual start %.2f vs %.2f (%d of %d)\n',all.ds(j),a,b,j-1);
    end
    [D,b] = max(all.ds);
    G = all.Gs{b};
    all.took = mytoc(ts);
    
function [d, ds] = gmdprime(G, varargin)
%calcualte drpime between two Gaussians in gmdistribution fit        
if size(G.mu,1) == 3
    %find three distances. one is allowed to be smaltt (splitting hash into
    %two. But one cluster must be distant from both of these. So take
    %lowest of top two = middle value
distance = mahal(G,G.mu);
    d(1) = sqrt(2./((1./distance(1,2))+(1./distance(2,1))));
    d(2) = sqrt(2./((1./distance(1,3))+(1./distance(3,1))));
    d(3) = sqrt(2./((1./distance(2,3))+(1./distance(3,2))));
    ds = d;
    d  = d(2);  
else
distance = mahal(G,G.mu);
    d = sqrt(2./((1./distance(2,1))+(1./distance(1,2))));
    ds = d;
end
    
