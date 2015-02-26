function S = CellToStruct(C, varargin)
% S = CellToStruct(C, varargin)
%converts a cell array to a structure. N.B quite different from cell2struct
%This is Inteded for cell arrays that have similar elements. cell2struct
%takes all the fields in each element of C, and then makes a structure with
%these fields.

allfields = {};
for j = 1:length(C);
    if isstruct(C{j})
    f = fields(C{j});
    allfields = unique({f{:} allfields{:}});
    for k = 1:length(f)
        S(j).(f{k}) = C{j}.(f{k});
    end
    end
end

