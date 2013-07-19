function [shift, xcs] = AlignMatrices(A,B, dim, varargin)

interp = 0;
j = 1;
sd = 0;

while j <= length(varargin)
    if strncmpi(varargin{j},'interp',4)
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            sd = varargin{j};
        else
            sd = 1;
        end
    end
    j = j+1;
end
if length(dim) == 2 %use variance over dim(1), summed across dim(2)
    A = squeeze(sum(var(A,1,dim(1)),dim(2)));
    B = squeeze(sum(var(B,1,dim(1)),dim(2)));
    dim = 2;
end

if dim == 2 && ndims(A) == 2;
    A = A';
    B = B';
    dim = 1;
end

nch = size(A,dim);
if nch == 24 %kludge. Need to figure this out properly...
    nch = 23;
end
ns = nch-2;

for j = -ns:0
    if ndims(A) == 2
    a = A(1:nch+j,:);
    b = B(1-j:nch,:);
    elseif ndims(A) == 4 && dim == 4
        a = A(:,:,:,1:nch+j);
        b = B(:,:,:,1-j:nch);
    end
    xc = corrcoef(a(:),b(:));
    xcs(1+j+ns) = xc(1,2);
    diffs(1+j+ns) = mean(abs(a(:)-b(:)));
    shifts(1+j+ns) = j;
end
for j = 1:ns
    if ndims(A) == 2
        a = A(1+j:nch,:);
        b = B(1:nch-j,:);
    elseif ndims(A) == 4 && dim == 4
        a = A(:,:,:,1+j:nch);
        b = B(:,:,:,1:nch-j);
    end
    xc = corrcoef(a(:),b(:));
    xcs(j+ns+1) = xc(1,2); 
    diffs(j+ns+1) = mean(abs(a(:)-b(:)));
    shifts(1+j+ns) = j;
end
if sd > 0
    ip = shifts(1):0.1:shifts(end);
    for j = 1:length(ip)
        ixcs(j) = sum(exp(-(shifts-ip(j)).^2./(2.*sd^2)) .* xcs);
    end
    [maxc,shift] = max(ixcs);
    shift = ip(shift);
    xcs = ixcs;
else
[maxc,shift] = max(xcs);
shift = shift- (ns+1);
end
