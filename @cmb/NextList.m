function NextList(a,b, varargin)
%DATA = cmb.combine('getstate');

j = 1;
saved = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'clearfirst',6)
cmb.DeleteCluster(0,a,'nodraw');
cmb.SetExptClusters(a,b,'nosave');
elseif strncmpi(varargin{j},'setfirst',6)
saved = cmb.SetExptClusters(a,b,'nosave');
if saved == 0
set(a,'backgroundcolor','r');  %as its not saved
end
end
j = j+1;
end

DATA.savedclusters = saved;
DATA = GetDataFromFig(a);
playspk = get(findobj(DATA.toplevel,'Tag','ShowSpikes'),'value');

if DATA.playingspk
set(findobj(DATA.svfig,'Tag','StopSpool'),'value',1);
return;
end
strs = get(DATA.elst,'string');
val = get(DATA.elst,'value');
if val < length(strs)
set(DATA.elst,'value',val+1)
if isfield(DATA,'AllClusters') & playspk == 0
cmb.combine('setexpt',DATA);
if 0 % need to use setexp to clear cluster.touched, etc
id = get(DATA.elst,'value');
DATA.currentexpt = DATA.expid(id);
DATA.allexp = DATA.currentexpt(1);
cmb.PlotAllProbeXY(DATA);
end
else
cmb.combine('setexpt',DATA);
end
end


