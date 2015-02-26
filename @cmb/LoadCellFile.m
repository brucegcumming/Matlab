function DATA = LoadCellFile(DATA)
DATA.CellDetails = []; 
DATA.muCellList = [];

if exist(DATA.cellfile,'file')
    load(DATA.cellfile);
    DATA.CellList = CellList;
    if exist('CellDetails','var')
        DATA.CellDetails = CellDetails;
        if isfield(CellDetails,'fitjumps')
            DATA.CellDetails.drift = [0 cumsum(CellDetails.fitjumps)];
        end
    end
    if exist('CellQuality','var')
        if size(CellQuality,2) < size(CellList,2)
            nq = size(CellQuality,2);
            nl = size(CellList,2);
            CellQuality(:,nq+1:nl) = 0;
        end
        DATA.CellQuality = CellQuality;
    else
        DATA.CellQuality = zeros(size(DATA.CellList));
    end
    if length(DATA.CellQuality) < length(DATA.CellList)
        lq = length(DATA.CellQuality);
        ll = length(DATA.CellList);
        DATA.CellQuality(:,lq+1:ll) = 0;
    end
    if exist('CellListCluster','var')
        DATA.CellListCluster = CellListCluster;
    else
        DATA.CellListCluster = zeros(size(CellList));
        id = find(CellList > 0);
        DATA.CellListCluster(id) = 1;
    end
    if exist('Templates','var')
        DATA.Templates = Templates;
    end
    if exist('TemplateInfo','var') && isempty(DATA.TemplateInfo)
        DATA.TemplateInfo = TemplateInfo;
    end
    if exist('muCellList','var')
        DATA.muCellList = muCellList;
    end
        
end
DATA.exabsid = 1:length(DATA.Expts);
DATA = cmb.CheckCellExptList(DATA);
if DATA.state.fixdrift
    if ~isfield(DATA.CellDetails,'probedrift')
        DATA.CellDetails.probedrift = celllist.drift(DATA);
    end
end

