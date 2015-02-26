function filename = GetProbeFilename(DATA, eid, probe)

id = find(DATA.probelist == probe);
if DATA.state.online == 0
if length(id) > 1
else
if DATA.probelist(id) > 16
filename = strrep(DATA.datafilename,'.mat',sprintf('A.p%s.mat',DATA.probevars{id}(3:end)));
else
filename = strrep(DATA.datafilename,'.mat',sprintf('.p%s.mat',DATA.probevars{id}(3:end)));
end
end
else
filename = ['C:' DATA.Expts{eid}.Header.Name];
if DATA.probelist(id) > 16
filename = strrep(filename,'/Expt','A/Expt');
end
end




