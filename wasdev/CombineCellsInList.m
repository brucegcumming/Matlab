function CellList = CombineCellsInList(name, cells, varargin)
%CombineCellsInList takes a CellList data file for combine and
%adds two cells toghether. Cells is a 2xn list of combinations
% cells(1,:) is moved to cells(2,:)


saving = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'save',4)
        saving = 1;
    end
    j = j+1;
end
load(name);

for j = 1:size(cells,1)
    to = cells(j,2);
    from = cells(j,1);
    id = find(CellList(from,:) > 0);
    if length(id) && sum(CellList(to,id)) == 0 %not already defined
        CellList(to,id) = CellList(from,id);
        CellQuality(to,id) = CellQuality(from,id);
        CellListCluster(to,id) = CellListCluster(from,id);
        CellList(from,id) = 0;
        CellQuality(from,id) = 0;
        CellListCluster(from,id) = 0;
    end
end

if saving
    save(name,'CellList','CellQuality','Templates','CellListCluster');
end