function C = CondenseDetails(C, varargin)

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'reverse',4)
        C = UnCondenseDetails(C);
        return;
    end
    j = j+1;
end
if isfield(C,'xy')
    if isinteger(C.xy)
        return;
    end
    C.maxv = max(abs(C.xy));
    C.maxint = 32000;
    for j = 1:size(C.xy,2)
        xy(:,j) = int16(round(C.xy(:,j) .* C.maxint./C.maxv(j)));
    end
    C.xy = xy;
end
if isfield(C,'triggerV')
    C.maxint = 32000;
    C.triggermax = max(abs(C.triggerV));
    t =  int16(round(C.triggerV .* C.maxint./C.triggermax));
    C.triggerV = t;
end
if isfield(C,'clst')
    C.clst = int16(C.clst);
end
if isfield(C,'next')
    for j = 1:length(C.next)
        C.next{j} = CondenseDetails(C.next{j});
    end
end

function C = UnCondenseDetails(C)

if ~isfield(C,'maxint')
    return;
end
if isfield(C,'xy')
    for j = 1:size(C.xy,2)
        xy(:,j) = double(C.xy(:,j)) .* C.maxv(j)./C.maxint;
    end
    C.xy = xy;
end    
if isfield(C,'triggerV')
    t = double(C.triggerV) .* C.triggermax./C.maxint;
    C.triggerV = t;
end
if isfield(C,'next')
    for j = 1:length(C.next)
        C.next{j} = UnCondenseDetails(C.next{j});
    end
end
C = rmfields(C,{'maxint' 'maxv' 'triggermax'});