function map = loadmap(file)

%
% map = loadmap(file)
%  reads a text file with penetration data in columns
% RFx RFy RFw RFh anem pen# penX penY VisualArea depth date
% and returns a matlab structure with corresponding elememts
%
MapDefs;
file
[rx, ry, rw, rh, ro, names, ids, px, py, area, depth, dates, hemi, fitx, fity] = textread(file,'%f%f%f%f%f%s%d%f%f%s%d%s%s%f%f%*[^\n]');


for j = 1:length(rx)
    map.rf(j,1) = rx(j);
    map.rf(j,2) = ry(j);
    map.rf(j,3) = rw(j);
    map.rf(j,4) = rh(j);
    map.rf(j,5) = ro(j);
    map.rf(j,6) = fitx(j);
    map.rf(j,7) = fity(j);
    map.cellname{j} = regexprep(names{j},'\.[0-9]*$','');
    if regexp(names{j},'M[0-9][0-9][0-9]')
        map.types(j) = MULTICONTACT;
    else
        map.types(j) = NORMAL;
    end
    map.pen(j,1) = ids(j);
    map.pen(j,2) = px(j);
    map.pen(j,3) = py(j);
    map.areastr{j} = area{j};
    aid = strmatch(area{j},areanames,'exact');
    if isempty(aid) 
        map.area(j) = VUNKNOWN;
    elseif  aid == VDEFAULT
        map.area(j) = V1;
    else
        map.area(j) = aid;
    end
    if hemi{j}(1) == 'R'
        map.hemisphere(j) = 1;
    else
        map.hemisphere(j) = 0;
    end
    if length(hemi{j}) > 1 && hemi{j}(2) == 'M'
        map.types(j) = MULTICONTACT;
    end
        
    map.depth(j) = depth(j);
    map.datestr{j} = dates{j};
    map.datenum(j) = datenum(dates{j});
    map.age(j) = now - map.datenum(j);
end    

cellnames = unique(map.cellname);
skipids = [];
for j = 1:length(cellnames)
    id = strmatch(cellnames{j},map.cellname,'exact');
    if length(id) > 1 && ~strcmp(cellnames{j},'Nocell')
        skipids = [skipids id(2:end)'];
    end
end
if length(skipids)
    id = setdiff(1:length(map.cellname),skipids);
    map.rf = map.rf(id,:);
    map.cellname = map.cellname(id);
    map.types = map.types(id);
    map.pen = map.pen(id,:);
    map.area = map.area(id);
    map.hemisphere = map.hemisphere(id);
    map.depth = map.depth(id);
    map.datestr = map.datestr(id);
    map.datenum = map.datenum(id);
    map.age = map.age(id);
end
map.filename = file;

if 0 
    map(j).rf(1) = rx(j);
    map(j).rf(2) = ry(j);
    map(j).rf(3) = rw(j);
    map(j).rf(4) = rh(j);
    map(j).rf(5) = ro(j);
    map(j).name = names{j};
    map(j).pen(1) = ids(j);
    map(j).pen(2) = px(j);
    map(j).pen(3) = py(j);
    map(j).area = area{j};
    map(j).depth = depth(j);
    map(j).datestr = dates{j};
    map(j).datenum = datenum(dates{j});
end    