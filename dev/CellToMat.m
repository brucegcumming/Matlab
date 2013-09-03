function M = CellToMat(C, f, varargin)
%CellToMat M = CellToMat(C, f,...
%takes the field f from each element of C, and puts into
% a matrix. Only works for cell arrays with 1 or 2 dimensions
%
% if the field is omitted, then CellToMat returns a matrix of size(C)
% elements that contain the length of C{i,j};
M = [];
if nargin < 2
    f = '';
end
dot = strfind(f,'.');
if ~isempty(dot)
    fa = f(dot+1:end);
    f = f(1:dot-1);
else 
    fa = [];
end
if size(C,1) == 1 && size(C,2) > 1
    for k = 1:length(C)
        if iscell(C{k})
            for j= 1:length(C{k})
                if isempty(f)
                    M(j,k) = size(C{k}{j});
                elseif isfield(C{k}{j},f)
                    M(k,j,1:length(C{k}{j}.(f))) = C{k}{j}.(f);
                else
                    M(k,j,:) = NaN;
                end
            end
        elseif isempty(f)
            M(k) = length(C{k});
        elseif isfield(C{k},f)
            if isempty(fa)
                x = C{k}.(f);
            elseif isfield(C{k}.(f),fa)
                x = C{k}.(f).(fa);
            end
           if iscell(x) || isstruct(x) 
            M{k,1:length(x)} = x;
           else
            M(k,1:length(x)) = x;
           end
        end
    end
else
    if iscell(C{1})
        C = C{1};
    end
for j = 1:size(C,1)
    for k = 1:size(C,2)
        if isfield(C{j,k},f)
            if isempty(fa)
                x = C{j,k}.(f);
            elseif isfield(C{j,k}.(f),fa)
                x = C{j,k}.(f).(fa);
            end
            M(j,k,1:length(x)) = x;
        end
    end
end
end