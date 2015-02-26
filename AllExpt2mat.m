function [M, details] = AllExpt2Mat(E, ptype,varargin)
%[M, details] = AllExpt2Mat(E, ptype,varargin)
%Extract a matrix of data from an allexpt
% ptype can be an index to a particular SDF from a RC expt
%  'Blank' to get the blank response
%
%the allexpt has to have the SDFs built first. See PlotAllCellFiles
M = [];
smsd = 0;
Array = [];
splitdistance = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'smooth',5)
        j = j+1;
        smsd = varargin{j};
    elseif strncmpi(varargin{j},'splitdistance',8)
        splitdistance = 1;
    elseif isfield(varargin{j},'X')
        Array = varargin{j};
    end
    j = j+1;
end

if isfield(E{1}.Header,'loadname')
    details.loadname = E{1}.Header.loadname;
end
details.celllist = CellToMat(E,'cellid');
details.probes = CellToMat(E,'Header.probe');
[probes, pid] = sort(details.probes);
details.probes = probes;
details.celllist = details.celllist(pid);
if strcmp(ptype,'blank')
    for j = 1:length(E)
        p = pid(j);
        if isfield(E{p}.plotres,'sdfs')
            id = find(E{p}.plotres(1).sdfs.extraval == -1009);
            M(j,:) = E{p}.plotres(1).sdfs.extras{id}.sdf;
            if smsd > 0
                M(j,:) = smooth(M(j,:),smsd,'gauss');
            end
            details.times = E{p}.plotres(1).times./10;
            details.n(p) = E{p}.plotres(1).sdfs.extras{id}.n;
        end
    end
elseif isnumeric(ptype)
    xvals = ptype;
    yvals = 1;
    for j = 1:length(E)
        p = pid(j);
        if isfield(E{p}.plotres,'sdfs')
            M(j,:) = E{p}.plotres(1).sdfs.s{xvals(1),yvals(1)};
            if smsd > 0
                M(j,:) = smooth(M(j,:),smsd,'gauss');
            end
            details.times = E{p}.plotres(1).times./10;
        end
    end    
end
        

if splitdistance
    M = myNormalize(M,'mean');
    x = mean(Array.X);
    y = mean(Array.Y);
    d = abs((Array.X + i .* Array.Y) - (x + i*y));
    x = prctile(d,50);
    a = find(d < x);
    b = find(d >= x);
    details.split(1,:) = mean(M(a,:));
    details.split(2,:) = mean(M(b,:));
end

