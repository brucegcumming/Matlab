function [pos, details] = drift(C, varargin)
%pos = celllist.drift(CellList) Estimates electrode drfit from cell numbers
%returns a vecotor of estimated electrode positions, as floats.
%an offset has been added to minimize errors from round(pos);

if isfield(C,'CellList') && isfield(C,'CellDetails')
    DATA = C;
    CD = DATA.CellDetails;
    C = DATA.CellList;
end


maxiter = 10000;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'save',4);
        saveres = 1;
    end
    j = j+1;
end

if ischar(C)
    if exist(C)
        DATA = load(C);
        DATA.loadname = C;
    end
    CD = DATA.CellDetails;
    C = DATA.CellList;
end
cells = unique(C(:));
cells = cells(cells>0);
nprobes = size(C,2);
T = zeros(length(cells),nprobes);
for j = 1:length(cells)
    [a,b] = find(C == cells(j)); %a is expts
    [c,d] = ind2sub([size(C,2) size(C,3)],b); %c is probes, disregargind cluster
    cellim(j, a) = c;
    im(j,:) = hist(c,[1:nprobes]);
    [a,b] = max(im(j,:));
    T(j,b) = a;
end
imagesc(im);
[a,p] = max(im');

%Each row of T is porbe estimates for cell. Start with 
%maxmima 

AllT(:,:,1) = T;
for j = 1:size(cellim,2)
    shiftim(j,:) = FindShift(T,cellim(:,j)); 
end
    shifts = im2shift(shiftim);
    AllShifts(1,:) = shifts;
for n = 1:3
    [a,b] = max(shiftim');
    T = MakeTemplate(cellim,b-8);
    AllT(:,:,1+n) = T;
    for j = 1:size(cellim,2)
        shiftim(j,:) = FindShift(T,cellim(:,j));
    end
    shifts = im2shift(shiftim);
    AllShifts(1+n,:) = shifts;
end
[a,b] = max(shiftim');
T = MakeTemplate(cellim,b-8);
%p is starting guess for cell probe location. 


function shift = im2shift(im)

shifts = 1:size(im,2);
for j = 1:size(im,1)
    shift(j) = sum(im(j,:).*shifts)./sum(im(j,:));
end


function T = MakeTemplate(cellim, shift)

T = zeros(size(cellim,1),24);
for j = 1:size(cellim,2)
    id = find(cellim(:,j)>0);
    for k = 1:length(id)
        p = id(k)+shift(j); %probe ins shifted mean
        if(p > 0)
            T(p,cellim(id(k),j)) = T(p,cellim(id(k),j)) +1;
        end
    end
end
size(T);
function score = FindShift(T, im)
%im is a list for each cell of which probe it is on.
%im(j) = 0 means that cell not present in this row.
for k = -7:7
    si = k+8;
    score(si) = 0;
    for j = 1:length(im)
        m = im(j)+k;
        if m > 0 && m < size(T,2)
            score(si) = score(si) + (im(j)>0) .* T(j,m);
        end
    end
end

return;

%think this is leftover stuff
[fitd, fitdetails] = FitDriftMatrix(d,'maxiter',maxiter);
pos = [0 cumsum(fitdetails.jumps(:)')];
dxs = 0:0.01:1;
for j = 1:length(dxs);
    x = pos+dxs(j);
    err(j) = sum((round(x)-x).^2);.1
end
[a,b] = min(err);
pos = pos+dxs(b);
[a,b] = Counts(round(pos));
pos = pos - b(1); %make 0 most common value
details.distances = d;
try
    details.mdpos = mdscale(d,1);
catch ME
    str = CheckExceptions(ME);
    details = AddError(details,str);
end

if saveres && ~isempty(loadname);
    DATA.CellDetails.probedrift = pos;
    save(loadneam,'-struct',DATA);
end
 