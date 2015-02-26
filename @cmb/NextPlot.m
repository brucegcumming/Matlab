function P =  NextPlot(P, AllCellRes, step)
R = [];
if P.sortbyvar
var = cat(1,AllCellRes.var);
id = find(var(:,2) < 0.4);
var(id,2) = 0.4;
[a, id] = sort(var(:,1)./var(:,2),'descend');
else
id = 1:P.nplots;
end
if step == 0 && P.currentcell < P.nplots
hold off;
P.currentcell = P.currentcell+1;
elseif step == -1 && P.currentcell > 1
hold off;
P.currentcell = P.currentcell-1;
else
P.currentcell = step;
end
R = PlotResult(AllCellRes(id(P.currentcell)));
if P.sortbyvar && ~isempty(R)
t = get(gca,'title');
set(t,'string',[get(t,'string') sprintf('  VR %.3f (%.3f/%.3f)',R.var(1)./R.var(2),R.var(1),R.var(2))]);
end

