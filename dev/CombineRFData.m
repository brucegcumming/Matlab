function allrfs = CombineRFData(varargin)
%rflist = CombineRFData(rflist, fits)
%Combines together RF date from ulf files (rflist) with fitted RFs
%(for multiple probes in array data) to produce one rflist for PlotMat

fits = []
rfs = [];
j = 1;
while j <= length(varargin)
    if iscell(varargin{j})
        if length(CellToMat(varargin{j},'depth')) > 0
            rfs = varargin{j};
        else
        %if length(CellToMat(varargin{j},'proberf')) > 0
            fits = varargin{j};
        end
    end
    j = j+1;
end

names = CellToMat(rfs,'name');
dates = CellToMat(rfs,'date');
[a,b] = Counts(names);
id = find(a > 1);
xid = [];
for j = 1:length(id)
    nid = find(strcmp(b{id(j)},names));
    if length(unique(dates(nid))) > 1
        xid = [xid nid(1:end-1)];
    else
        xid = [xid nid(1:end-1)];
    end
end
gid = setdiff(1:length(rfs),xid);
rfs = rfs(gid);
for j = 1:length(rfs)
    [a, rfnames{j}] = fileparts(rfs{j}.name);
end
allrfs = rfs;
for j = 1:length(fits)
    if isfield(fits{j},'dirname');
        [a, fitnames{j}] = fileparts(fits{j}.dirname);
        id = find(strcmp(fitnames{j},rfnames));
        if isempty(id)
            nfit(j,1) = 0;
        else
            nfit(j,1) = size(fits{j}.proberf,1);
            for k = 1:nfit(j,1)
                newrf.rf(1:10) = fits{j}.proberf(k,1:10); %x,y
                newrf.depth = fits{j}.proberf(k,11);
                newrf.probe = fits{j}.proberf(k,12);
                newrf.name = fits{j}.dirname;
                newrf.electrode = fits{j}.electrode;
                allrfs = {allrfs{:} newrf};
            end
        end
        nfit(j,2) = length(fits{j}.fits);
    else
        nfit(j,1) =0;
    end
end

sum(nfit);
id = find(nfit(:,2) == 0);
for j = 1:length(id)
    fprintf('%s No fit Data\n',rfs{id(j)}.name);
end
id = find(nfit(:,1) == 0 & nfit(:,2) > 0);
for j = 1:length(id)
    fprintf('%s Too few fits\n',rfs{id(j)}.name);
end

