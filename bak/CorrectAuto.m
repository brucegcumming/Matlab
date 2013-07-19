function res = CorrectAuto(dir,varargin)

j=1;
save=0;
while j<=length(varargin)
    if strncmpi(varargin{1},'save')
        save=1;
    else
    end
    j=j+1;
end
a=TreeFind(dir,'name','AutoClusterTimes.mat')
for i=1:length(a)
    k{i}=regexprep(a{i},'.*Expt([0-9])*.*','$1');
end
for i=1:length(k)
    T(i)=str2num(k{i});
end
N=sort(T);
res=NaN(length(k),24);
for k=N
    G=[dir,'/Expt',int2str(k),'AutoClusterTimes.mat'];
    load(G);
    AutoClusters=Clusters;
    try
    F=[dir,'/Expt',int2str(k),'ClusterTimes.mat'];
    load(F);
    catch
        continue
    end
    for j=1:24
        if (AutoClusters{j}.savetime - Clusters{j}.savetime(1:2)) == 0
            if save
                Clusters{j}.auto = 1;
                save(F,'Clusters')
            end
            res(k,j+1)=1;
        else
            res(k,j+1)=0;
        end
    end
end
return
