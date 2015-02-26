function Cd = UpdateClusterDetails(DATA, varargin)
%Not a good idea. Do't want to change details if the cluster is not saved
Cd = {};

p = AllV.ProbeNumber(DATA);
Cd = getappdata(DATA.toplevel,'ClusterDetails')
if length(Cd) >= p && isfield(Cd{p},'clst')
    Cd{p}.clst = DATA.clst;
    Cd{p}.t = DATA.t;
    Cd{p}.triggerV = DATA.rVAllV;
    AllV.mysetappdata(DATA,'ClusterDetails',Cd);
end