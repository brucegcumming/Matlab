
function DATA = CheckForPCA(DATA, ispk, force)
    needpca = 0;
    p = DATA.probe;
    if isempty(DATA.cluster)
        needpca = 0;
    elseif sum(ismember([DATA.plot.clusterX DATA.plot.clusterY DATA.plot.clusterZ],[33 34 37 38 39])) %need PCA
       needpca = 1;
   elseif p <= size(DATA.cluster,2) & isfield(DATA.cluster{1,p},'params') & sum(ismember(DATA.cluster{1,p}.params,[33 34 37 38]))
       needpca = 1;
   end
   if needpca
       if isfield(DATA,'AllSpikes')
           for j = p
               if ~isempty(DATA.AllSpikes{j})
                   [a,b,c] = OneProbePCA(DATA,ispk);
                   DATA.AllSpikes{j}.pcs(ispk,:) = c';
                   DATA.AllSpikes{j}.EigVal = diag(b);
               end
           end
       elseif force || length(DATA.AllData.pcs) < max(ispk) || std(DATA.AllData.pcs(ispk,1)) < 0.01
           id = find(DATA.AllData.Spikes.codes(ispk,2) < 1);
        [a,b,c] = OneProbePCA(DATA,ispk,id);
        DATA.AllData.pcs(ispk,:) = c';
        DATA.AllData.EigVal = diag(b);
        DATA.AllData.EigVec = a;
       end
   end


function [E,V, pc] = OneProbePCA(DATA, ids, subid)

    if isfield(DATA,'AllSpikes')
        spks = DATA.AllSpikes{DATA.probe}.values(ids,:);
    else
        spks = DATA.AllData.Spikes.values(ids,:);
    end
    if nargin  == 3
         [E,V] = eig(cov(spks(subid,:)));
    else
         [E,V] = eig(cov(spks));
    end
        pc(1,:) = (E(:,end)'-mean(E(:,end))) * spks';
        pc(2,:) = (E(:,end-1)'-mean(E(:,end-1))) * spks';
        pc(3,:) = (E(:,end-2)'-mean(E(:,end-2))) * spks';
        pc(4,:) = (E(:,end-3)'-mean(E(:,end-3))) * spks';
