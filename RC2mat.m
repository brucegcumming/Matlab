function x = RC2mat(rc, varargin)
%x = RC2mat(rc, varargin)
% convert an RC structure (from PlotRevCorAny)
% to a 3-D matrix (R(t), xval, yval)
% and a 4-D matcix (R(t), xval, yval, probe) for LFP

addblank = 0;
adduc = 0;
addextras = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'addextra',4)
        addextras = 1;
    end
    j = j+1;
end

if isfield(rc.sdfs,'lfp')
    for j = size(rc.sdfs.s,1):-1:1
        for k = size(rc.sdfs.s,2):-1:1
            x(:,j,k,:) = rc.sdfs.lfp{j,k};
        end
    end
else
    for j = size(rc.sdfs.s,1):-1:1
        for k = size(rc.sdfs.s,2):-1:1
            x(:,j,k) = rc.sdfs.s{j,k};
        end
    end
end
