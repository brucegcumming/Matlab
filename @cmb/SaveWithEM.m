function SaveWithEM(a,b)
DATA = GetDataFromFig(a);
outname = get(DATA.saveitem,'string');
Expt = DATA.Expt;
Expt.Header.CombineDate = now;
if DATA.state.online
for j = 1:length(DATA.combineids)
end
else
Expt = LoadEmData(Expt);
end
fprintf('Saving %s with Eye position Data\n',outname);
save(outname,'Expt');


