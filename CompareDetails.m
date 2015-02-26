function CompareDetails(dpath, varargin)
% CompareDetails(dpath, varargin) Check folder dpath for erros in
% ClusterDetails files. Checks that AutoClusterDetails is not newer.
% dpath should end with '/' if its a folder. If it expts ExptN, then
%only files matching ExptN*AutoClusterTimesDetails.mat are checked

d = mydir([dpath '*AutoClusterTimesDetails.mat']);
fixdate = 0;

j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'fixdate',6)
        fixdate = 1;
    end
    j = j+1;
end

for j = 1:length(d)
    cfile = strrep(d(j).name,'AutoCluster','Cluster');
    if exist(cfile)
        fixed = 0;
        fprintf('Checking %s\n',cfile);
        X = load(d(j).name);
        A = X.ClusterDetails;
        D = load(cfile);
        C = D.ClusterDetails;
        X = load(strrep(cfile,'TimesDetails','Times'));
        for c = 1:length(C)
            if C{c}.ctime < A{c}.ctime
                if X.Clusters{c}.ctime > C{c}.ctime
                    str = 'And ClusterTimes';
                    if fixdate && X.Clusters{c}.auto == 1
                        C{c} = A{c};
                        fixed = fixed+1;
                    end
                else
                    str = 'Just Details';
                end
                fprintf('AutoDetails %d,%d is newer (%s) %s vs %s\n',j,c,str,datestr(C{c}.ctime),datestr(A{c}.ctime));
            elseif X.Clusters{c}.ctime > C{c}.ctime                
                if X.Clusters{c}.manual == 2
                    fprintf('Unquantified PlotClusters Cut in %d,%d\n',j,c);
                elseif X.Clusters{c}.quick
                    fprintf('Unquantified Quick Cut in %d,%d\n',j,c);
                else
                fprintf('CLusters Newer than Details %d,%d: %s vs %s\n',j,c,datestr(X.Clusters{c}.ctime),datestr(C{c}.ctime));               
                end
            end
        end
        if fixed > 0
            ClusterDetails = C;
            FullVData = D.FullVData;
            save(cfile,'ClusterDetails','FullVData');
        end
    end
end
