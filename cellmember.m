function x = cellmember(a,b, varargin)
%cellmember(a,b) ismember for cell arrays of doubles.
%if a is a scalar, b a cell array, returns which cellarrays contain a

if isnumeric(a) && iscell(b)
    for j = 1:length(b)
        if isnumeric(b{j})
            x(j) = sum(ismember(a,b{j}));
        end
    end
end