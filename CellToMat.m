function [M, details] = CellToMat(C, f, varargin)
%CellToMat M = CellToMat(C, f,...
%takes the field f from each element of C, and puts into
% a matrix. Only works for cell arrays with 1 or 2 dimensions
%
% CellToMat(C, f, 'pad') forces the length of M to match the length of C
% if the field is omitted, then CellToMat returns a matrix of size(C)
% elements that contain the length of C{i,j};
M = [];
padlength = 0;
if nargin < 2
    f = '';
end

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'pad',3)
        padlength = 1;
    end
    j = j+1;
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
                    M(j,k) = length(C{k}{j});
                elseif isfield(C{k}{j},f)
                    M(k,j,1:length(C{k}{j}.(f))) = C{k}{j}.(f);
                else
                    M(k,j,:) = NaN;
                end
            end
           details.found(k) = 1;
        elseif isempty(f)
            M(k) = length(C{k});
           details.found(k) = 1;
        elseif isfield(C{k},f)
            if isempty(fa)
                x = C{k}.(f);
            elseif isfield(C{k}.(f),fa)
                x = C{k}.(f).(fa);
            else
                x = NaN;
            end
           if iscell(x) || isstruct(x) 
               M{k,1:length(x)} = x;
            elseif ischar(x)
                M{k} = x;
           else
               if iscell(M)
                   M{k} = x;
               elseif sum(size(x)>1) ==2
                   M(k,:,:) = x;
               else
                   M(k,1:length(x)) = x;
               end
           end
           details.found(k) = 1;
        else
           details.found(k) = 0;
        end
    end
    if padlength && details.found(end) == 0
        M(k) = NaN;
    end
else
    if ~isempty(C) && iscell(C{1})
        C = C{1};
    end
for j = 1:size(C,1)
    for k = 1:size(C,2)
        if isempty(f)
            M(j,k) = length(C{j,k});
        elseif isfield(C{j,k},f)
            if isempty(fa)
                x = C{j,k}.(f);
            elseif isfield(C{j,k}.(f),fa)
                x = C{j,k}.(f).(fa);
            else
                x = NaN; %subfield absent
            end
            M(j,k,1:length(x)) = x;
        end
    end
end
end