function X = ListClusterBackup(prefix,varargin)
%ListClusterBackup(prefix)
%Finds cluster backup files and compares them with current definitions
j = 1;
while j <= length(varargin)
    j = j+1;
end
if isfield(prefix,'savetime') %plot previous result
    PlotBackup(prefix,varargin{:});
    return;
end

d = mydir([prefix '/*ClusterTimes*.mat']);
path = fileparts(prefix);
[dates, did] = sort([d.datenum],'descend');
nb = 0;
for j = 1:length(d)
    clear Clusters;
    load(d(did(j)).name);
    if exist('Clusters','var')
        nb = nb+1;
        saves = CellToMat(Clusters,'savetime');
        [a,b] = max(saves(:,end));
        X(nb).savetime =a;
        X(nb).p = Clusters{b}.probe(1);
    end
%    fprintf('%s %d %d %s\n',d(did(j)).filename,b,Clusters{b}.manual,datestr(a));
end
PlotBackup(X);

function PlotBackup(X, varargin)
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'maxage',6)
        j = j+1;
        id = find(now - [X.savetime] < varargin{j});
        X = X(id);
    end
    j = j+1;
end


plot([X.savetime],[X.p],'o');
datetick;


