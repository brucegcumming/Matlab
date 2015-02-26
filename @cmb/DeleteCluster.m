function DeleteCluster(cl,callfig, varargin)
DATA = GetDataFromFig(callfig);
redraw = 1;

j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'nodraw',4)
redraw = 0;
end
j = j+1;
end
p = get(gca,'UserData');
if isempty(p)
p = DATA.probe;
end
if cl == 0
for j = 1:size(DATA.cluster,1)
cmb.DeleteCluster(j, callfig, varargin{:});
end
return;
end
if ~isempty(DATA.cluster{cl,p}) & isfield(DATA.cluster{cl,p},'h') &  ishandle(DATA.cluster{cl,p}.h)
delete(DATA.cluster{cl,p}.h);
end
DATA.cluster{cl,p} = [];
DATA.cluster{cl,p}.touched = 1; %records active deletion. Will be copied to Expt
DATA.cluster{cl,p}.quality = 0; %records active deletion. Will be copied to Expt
DATA.cluster{cl,p}.deleted = 1; %records active deletion. Will be copied to Expt
if isfield(DATA.Expts{DATA.currentexpt(1)}.Cluster{cl,p},'forceid')
DATA.Expts{DATA.currentexpt(1)}.Cluster{cl,p} = rmfield(DATA.Expts{DATA.currentexpt(1)}.Cluster{cl,p},'forceid');
end
DATA.forceclusterid = 0;
DATA.Expts{DATA.currentexpt(1)}.Cluster{cl,p}.quality = 0;
if isfield(DATA.Expts{DATA.currentexpt(1)}.gui,'spkrange')
ispk = DATA.Expts{DATA.currentexpt(1)}.gui.spkrange;
ispk = [ispk(1):ispk(2)];
if strncmp(DATA.filetype,'Grid',4)
ctype = 1;
else
ctype = 2;
end
if isfield(DATA,'AllSpikes')
spks = find(DATA.AllSpikes{DATA.probe}.codes(ispk,ctype) == cl);
DATA.AllSpikes{DATA.probe}.codes(ispk(spks),ctype) = 0;
else
spks = find(DATA.AllData.Spikes.codes(ispk,ctype) == cl);
DATA.AllData.Spikes.codes(ispk(spks),ctype) = 0;
end
if redraw
cmb.PlotSpikeXY(DATA,ispk(spks),DATA.spkcolor{1});
end
% Shouldn't need this. not deleting all clusters. And if we were would
% imply that we are setting not clusters
% DATA.Expts{DATA.currentexpt(1)}.gui.classified = 0;
end
set(DATA.toplevel,'UserData',DATA);
if cl > 1
newc = cl-1;
elseif cmb.iscluster(DATA.cluster,cl+1,DATA.probe) == 1
newc = cl+1;
else
newc = 1;
end
it = findobj('Tag','Clusterid');
set(it,'value',newc);

