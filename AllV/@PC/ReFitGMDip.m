 function DATA = ReFitGMDip(DATA, varargin)      ts = now;     Clusters = getappdata(DATA.toplevel,'Clusters');     cid = 1:length(Clusters);     j = 1;     while j <= length(varargin)         if strncmpi(varargin{j},'expts',5)             j = j+1;             cid = varargin{j};         end         j = j+1;     end     warning('off','stats:gmdistribution:FailedToConverge');     warning('off','stats:gmdistribution:MaxIterations');     for j = cid         for k = 1:length(Clusters{j})         C = Clusters{j}{k};         [a,b] = GMDip(C.xy(:,1),[],'idlist',C.clst,'noplot');         DATA.GaussFitdp(j,k,1) = b.gmdprime;         [DATA.GaussFitdp(j,k,2), b, c] = PC.Fit2Gauss(C);         DATA.GaussFitdp(j,k,3) = now; %record fit time so can update         DATA.gmfitpos(j,k,:) = c.fitpos;         fprintf('E%dP%d %.2f %.2f %d %d\n',j,k,DATA.GaussFitdp(j,k,1),DATA.GaussFitdp(j,k,2),DATA.gmfitpos(j,k,1),DATA.gmfitpos(j,k,2));         end         fprintf('Expt %d took %.0f sec\n',j,mytoc(ts));     end     mytoc(ts);     PC.SetFigure(DATA,DATA.tag.clusters);     hold off;     for j = 1:length(Clusters)         for k = 1:length(Clusters{j})             if sum(DATA.gmfitpos(j,k,:)) == 2         plot(DATA.GaussFitdp(j,k,1),DATA.GaussFitdp(j,k,2),'o','buttondownfcn',{@PC.HitPopPoint, j,k});             else         plot(DATA.GaussFitdp(j,k,1),DATA.GaussFitdp(j,k,2),'ro','buttondownfcn',{@PC.HitPopPoint, j,k});             end         hold on;         end     end     refline(1);     PC.SaveExtras(DATA);     set(DATA.toplevel,'UserData',DATA);