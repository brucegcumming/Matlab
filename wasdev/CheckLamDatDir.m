function CheckLamDatDir(path, varargin)

d = dir(path);
root = fileparts(path);
nlist = 0;
ncell = 0;
nc = 0;
for j = 1:length(d)
    if d(j).isdir && d(j).name(1) ~= '.'
        CheckLamDatDir([root '/' d(j).name], varargin{:});
    elseif regexp(d(j).name, 'p[0-9]*cl.mat') %cluster defs
        nc = nc+1;
        cid(nc) = j;
    elseif regexp(d(j).name, '.cell[0-9].*.mat')
        ncell = ncell+1;
        cellid(ncell) = j;
    elseif regexp(d(j).name, '.cells.mat')
        nlist = nlist+1;
        listid(nlist) = j;
    end
end
if nc || nlist
    if nc
        lastclusterset = max([d(cid).datenum]);
    end
    if nlist 
        listdate = min([d(listid).datenum]);
        if nc && listdate < lastclusterset
            fprintf('Clusters changed since list\n');
        end
    end
    if ncell
        oldcell = min([d([cellid]).datenum]);
        fprintf('Clusters last set %s Oldeset Cell %s, List %s\n',...
            datestr(lastclusterset),datestr(oldcell),datestr(listdate));
        fprintf('%d cells, %d clusters\n',ncell,nc);
        if ncell && oldcell < max([lastclusterset listdate])
            id = find([d(cellid).datenum] < max([lastclusterset listdate]));
            for j = id
                fprintf('%s (%s) out of date\n',[d(cellid(id(j))).name],datestr(d(cellid(id(j))).datenum));
            end
        end
    end
end
    