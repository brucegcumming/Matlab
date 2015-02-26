function MakeTemplates(a,b, varargin)
DATA = GetDataFromFig(a);
eid = DATA.currentexpt(1);
if isfield(DATA.Expts{eid}.gui,'spks')
ispk = DATA.Expts{eid}.gui.spks;
else
end
if isfield(DATA,'AllSpikes')
Spks = DATA.AllSpikes{DATA.probe};
else
Spks = DATA.AllData.Spikes;
end
nc = length(unique(Spks.codes(ispk,2)));
for j = 1:nc
id = find(Spks.codes(ispk,2) == j-1);
Template(j,:) = mean(Spks.values(ispk(id),:));
DATA.TemplateInfo(j).exid = eid;
DATA.TemplateInfo(j).probe = DATA.probe;
DATA.TemplateInfo(j).cluster = DATA.currentcluster-1;
DATA.TemplateInfo(j).cellid = 1;
end
DATA.Templates = Template;
cmb.PlotTemplates(DATA);
cmb.SaveCellList(DATA);

set(DATA.toplevel,'UserData',DATA);
cmb.NotBusy(DATA);

