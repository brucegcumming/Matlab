function DATA = CalcDistances(DATA)    Clusters = getappdata(DATA.toplevel,'Clusters');    nc = max(DATA.nclusters(:));    tic;     for e = 1:length(Clusters)         for p = 1:length(Clusters{e})             for c = 1:nc                 d(e,p,c) = PC.DistanceMeasure(Clusters{e}{p},c,DATA.mahaltype);             end         end     end    toc    DATA.mdistance = d;    