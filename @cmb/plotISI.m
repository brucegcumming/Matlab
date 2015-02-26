function plotISI(DATA)

if DATA.plot.showISI
GetFigure('ISI');
isis = CalcISI(DATA.Expts{DATA.currentexpt(1)}.Trials);
id = find(isis < 1000)
hist(isis(id),100);
end

