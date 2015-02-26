function [allrfs, details] = CombineRFData(varargin)
%rflist = CombineRFData(rflist, fits)
%Combines together RF date from ulf files (rflist) with fitted RFs
%(for multiple probes in array data) to produce one rflist for PlotMat

fits = [];
rfs = [];
d = [];
j = 1;
while j <= length(varargin)
    if iscell(varargin{j})
        if sum(CellToMat(varargin{j},'nfits')) > 0
            if isempty(fits)
                fits = varargin{j};
            else %two sets of fits
                rfs = varargin{j};
                [allrfs, details] = CombineFits(fits, varargin{j});
                return;
            end
        elseif length(CellToMat(varargin{j},'readtime')) > 0
            rfs = varargin{j};
        end
    elseif isstruct(varargin{j}) && isfield(varargin{j},'name')
        d = varargin{j};
    end
    j = j+1;
end

%first get names/dates. Don't use CellToMat
% it may be that some rfs are structs and some are cells
for j = 1:length(rfs)
    if isfield(rfs{j},'name')
        names{j} = rfs{j}.name;
    elseif iscell(rfs{j})
        names{j} = 'CelStruct';
        xid(j) = 1;
    else
        names{j} = '??';
    end
    if isfield(rfs{j},'date')
        dates{j} = rfs{j}.date;
    end
end


%remove duplicate names. DOn't actually have a good way
% yet, Just takes last. But could use date - started below
[a,b] = Counts(names);
id = find(a > 1);
for j = 1:length(id)
    nid = find(strcmp(b{id(j)},names));
    if length(unique(dates(nid))) > 1
        xid(nid(1:end-1)) = 1;
    else
        xid(nid(1:end-1)) = 1;
    end
end


gid = setdiff(1:length(rfs),find(xid>0));
rfs = rfs(gid);
for j = 1:length(rfs)
    [a, rfnames{j}] = fileparts(rfs{j}.name);
end
allrfs = rfs;
fits = RemoveDuplicates(fits);
for j = 1:length(fits)
    if isempty(fits{j}) && length(d) == length(fits)
        fits{j}.dirname = d(j).name;
    end
    if isfield(fits{j},'fits');
        if ~isfield(fits{j},'spacing')
            cprintf('red','Missing spacing info for %s\n',fits{j}.dirname);
            fits{j}.spacing = 1;
        end
        [a, fitnames{j}] = fileparts(fits{j}.dirname);
        id = find(strcmp(fitnames{j},rfnames));
        if isempty(id) || isempty(fits{j}.proberf)
            nfit(j,1) = 0;
        else
            if isfield(fits{j},'area')
                area = fits{j}.area;
            else
                area = 'unknown';
            end
            nfit(j,1) = size(fits{j}.proberf,1);
            rfgroups = [0 nfit(j,1)];
            if std(fits{j}.proberf(:,1)) + std(fits{j}.proberf(:,2)) > 1
                depths = fits{j}.proberf(:,11)+fits{j}.proberf(:,12) .* fits{j}.spacing./1000;
                [depths,b] = sort(depths);
                fits{j}.proberf = fits{j}.proberf(b,:);
                rfs = (fits{j}.proberf(:,1)-fits{j}.proberf(:,9)) +i*(fits{j}.proberf(:,2)-fits{j}.proberf(:,10));
                rfd = [];
                for k = 1:size(fits{j}.proberf,1)-1
                    rfd(k) = abs(rfs(k)-rfs(k+1));
                end
                dd = diff(depths)';
                id = find(dd > 0.5 & rfd > 1); %jumps in RF with gap
                if ~isempty(id)
                    rfgroups = [0 id(1) length(dd)];
                end
            end
            for na = 1:length(rfgroups)-1
                fitid = (1+rfgroups(na)):rfgroups(na+1);
            for k = fitid(:)'
                newrf.rf(1:10) = fits{j}.proberf(k,1:10); %x1,y
                newrf.depth = fits{j}.proberf(k,11)+fits{j}.proberf(k,12) .* fits{j}.spacing./1000;
                newrf.probe = fits{j}.proberf(k,12);
                newrf.name = fits{j}.dirname;
                newrf.electrode = fits{j}.electrode;
                if na > 2
                    newrf.area = [area '*']; 
                elseif na > 1
                    newrf.area = [area '*']; 
                else
                    newrf.area = area; 
                end
                allrfs = {allrfs{:} newrf};
            end
            end
        end
        nfit(j,2) = length(fits{j}.fits);
    else
        nfit(j,1) =0;
    end
end


    
sum(nfit);
id = find(nfit(:,2) == 0);
nerr = 0;
for j = 1:length(id)
    nerr = nerr+1;
    fprintf('%s No fit Data\n',fits{id(j)}.dirname);
    details.missing(nerr).name = GetMonkeyName(fits{id(j)}.dirname,'expname');
    details.missing(nerr).id = id(j);
    details.missing(nerr).err = 0;
    details.missing(nerr).nfit = 0;
end
id = find(nfit(:,1) == 0 & nfit(:,2) > 0);
for j = 1:length(id)
    nerr = nerr+1;
    fprintf('%s Too few fits\n',fits{id(j)}.dirname);
    details.missing(nerr).name = GetMonkeyName(fits{id(j)}.dirname,'expname');
    details.missing(nerr).id = id(j);
    details.missing(nerr).err = 1;
    details.missing(nerr).nfit = nfit(id(j),2);
end

function [fits, details] = CombineFits(fits, newfits)
details.modified = 0;
nfit = CellToMat(fits,'nfits','pad');
rf = CellFields(fits,'rf');
fits = fits(find(nfit(:) > 0 | rf(:) > 0));
for j = 1:length(fits)
    [a,b,c,oldnames{j}] = GetMonkeyName(fits{j}.dirname);
    if isfield(fits{j},'date')
        olddates(j) = fits{j}.date;
    end
    if isfield(fits{j},'savedate')
        oldfitdates(j) = fits{j}.savedate;
    elseif isfield(fits{j},'fitdate')
        oldfitdates(j) = fits{j}.fitdate;
    end
end
nfit = CellToMat(newfits,'nfits','pad');
rf = CellFields(newfits,'rf');
newfits = newfits(find(nfit(:) > 0 | rf(:) >0));
for j = 1:length(newfits)
    [a,b,c,newnames{j}] = GetMonkeyName(newfits{j}.dirname);
    if isfield(newfits{j},'date')
        newdates(j) = newfits{j}.date;
    end
    if isfield(newfits{j},'savedate')
        newfitdates(j) = newfits{j}.savedate;
    elseif isfield(newfits{j},'fitdate')
        newfitdates(j) = newfits{j}.fitdate;
    end
end
[a,b,c] = intersect(newnames,oldnames);
newer = find(newfitdates(b) > oldfitdates(c));
if ~isempty(newer)
    fits(c(newer)) = newfits(b(newer));
    details.modified = 1;
    details.updated = b(newer);
end
[a,b] = setdiff(newnames,oldnames);
if ~isempty(a)
    fits = {fits{:} newfits{b}};
    details.modified = 1;
    details.newfiles =1;
end

function allrfs = RemoveDuplicates(rfs)
allrfs = {};

names = GetMonkeyName(CellToMat(rfs,'name'),'expname');
[n, unames] = Counts(names);
id = find(n > 1);
for j = 1:length(n)
    rid = find(strcmp(unames{j},names));
    allrfs(j) = rfs(rid(1));
    for k = 2:length(rid)
        if isfield(rfs{rid(k)},'proberf')
            if isfield(allrfs{j},'proberf')
                allrfs{j}.proberf = cat(1,allrfs{j}.proberf,rfs{rid(k)}.proberf);
            else
                allrfs(j) = rfs(rid(k));
            end
        end
    end
end


