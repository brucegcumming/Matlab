function pid = GetProbe(DATA, eid, probe)
%  pid = GetProbe(DATA, eid, probe)
%Finds currently selected probe from GUI
%If what is selected is a cell, returns correct probe from DATA.CellList

    if  DATA.listbycell == 1
        if strncmp(DATA.filetype,'Grid',4)
            rowid = DATA.cellexid(eid);
            if rowid > 0
                [pid, cid] = find(squeeze(DATA.CellList(rowid(1),:,:)) == DATA.probe);
                if size(DATA.CellList,3) == 1
                    pid = cid;
                    cid = 1;
                end
                if length(pid) ~= 1 %> 1 if cell defined twice
                    pid = 0;
                end
            else
                pid = 0;
            end
        else
        cid = regexp(DATA.Expts{eid}.Header.Name,'\.[0-9]*\.mat');
        en = sscanf(DATA.Expts{eid}.Header.Name(cid+1:end),'%d');
        rowid = find(DATA.CellDetails.exptids == en);
        if isempty(rowid)
            pid = 0;
        else
            [pid, cid] = find(squeeze(DATA.CellList(rowid(1),:,:)) == DATA.probe);
            if length(pid) ~= 1 %> 1 if cell defined twice
                pid = 0;
            end
        end
        end
    else
        pid = probe;
    end
