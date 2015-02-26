function [DATA, D] = CheckClusterLoaded(DATA, eid, pid, varargin)
%[DATA, D] = CheckClusterLoaded(DATA, eid, pid, varargin)
% eid is index in DATA.expntos
needcx = 0;
reload = 0;
D.loaddur = [0 0 0];
j = 1;

if nargin == 2
    if isfield(DATA,'probe')
        pid = DATA.probe;
    else
        pid = 1;
    end
end
if DATA.listbycell == 2
    needcx = 1;
end

if DATA.state.somespikes == 2 || DATA.state.usensx == 2
    %if usenex == 2 it means reading spikes from NEV files. For now, can read these every time
    if (length(DATA.AllClusters) < eid || length(DATA.AllClusters{eid}) < pid)...
            || DATA.state.usensx == 1
        if strncmp(DATA.filetype,'Grid',4)
            DATA =  cmb.GetNS5Spikes(DATA, DATA.currentexpt(1),  pid);
        end
    end
elseif length(DATA.AllClusters) < eid || isempty(DATA.AllClusters{eid}) ...
        || isempty(DATA.AllClusters{eid}(pid).times) ...
        || (isempty(DATA.AllClusters{eid}(pid).cx) && needcx)
    [DATA, D]  = cmb.ReadCluster(DATA, eid, pid);
elseif isfield(DATA,'Clusterfile') && length(DATA.Clusterfile) >= eid
    X = DATA.Clusterfile{eid};
    if isfield(X,'loadname')
        d = dir(X.loadname);
        if d.datenum > X.loadtime || reload
            fprintf('%s Modified since %s. Reloading\n',X.loadname,datestr(X.loadtime));
            [DATA, D] = cmb.ReadCluster(DATA, eid, pid);
        end
    end
end
if DATA.state.showspikes && DATA.state.somespikes == 1
    DATA = cmb.LoadSpikes(DATA, eid,'nocheck');
end
if isempty(D)
    D.loaddur = [0 0 0];
else
    DATA.newload = 1;
end


