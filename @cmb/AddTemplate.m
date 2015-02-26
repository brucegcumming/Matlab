function AddTemplate(a,b, varargin)
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
id = find(Spks.codes(ispk,2) == DATA.currentcluster);
Template(1,:) = mean(Spks.values(ispk(id),:));
DATA.Templates = cat(1, DATA.Templates, Template);
nt = size(DATA.Templates,1);
DATA.TemplateInfo(nt).exid = eid;
DATA.TemplateInfo(nt).probe = DATA.probe;
DATA.TemplateInfo(nt).cluster = DATA.currentcluster;
cmb.GetCellNumber(DATA, probe);

DATA.TemplateInfo(nt).cellid = nt;
cmb.PlotTemplates(DATA);
cmb.SaveCellList(DATA);
set(DATA.toplevel,'UserData',DATA);


