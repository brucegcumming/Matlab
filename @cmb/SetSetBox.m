function SetSetBox(DATA, ei)

cid = findobj('Tag','ClusterIsSet');
if ismember(DATA.Expts{ei}.gui.clustertype, [0 2 3]) %% The online cut or autocut
set(cid,'value',0);
else
set(cid,'value',1);
end

