function NewCluster(a)
DATA = GetDataFromFig(a);
j = size(DATA.cluster,1);
while isempty(DATA.cluster{j,DATA.probe}) & j > 1
j = j-1;
end
newc = j+1;
it = findobj('Tag','Clusterid');
set(it,'value',newc);

