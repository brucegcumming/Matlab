function gid = SelectFits(fits,varargin)
%gid = SelectFits(fits,varargin)
%Find fits in cell array of fits

depthrange = [];
probes = [];
varcrit = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'depth',4)
        j = j+1;
        depthrange = varargin{j};
    elseif strncmpi(varargin{j},'probe',4)
        j = j+1;
        probes = varargin{j};
    elseif strncmpi(varargin{j},'pvar',4)
        j = j+1;
        varcrit = varargin{j};
    end
    j = j+1;
end

pv = CellToMat(fits,'pvar');
gid = find(pv >= varcrit);
if ~isempty(depthrange)
    d = CellToMat(fits,'depth');
    id = find(d >= depthrange(1) & d <= depthrange(2));
    gid = intersect(gid, id);
end
if ~isempty(probes)
    p = CellToMat(fits,'probe');
    id = find(ismember(p,probes));
    gid = intersect(gid, id);
end
