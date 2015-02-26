function fc = CellFields(C, f, varargin)
%CellFields(C, field, varargin) finds elements of C contianing field f

for j = 1:length(C)
    if isfield(C{j},f)
        fc(j) = 1;
    else
        fc(j) = 0;
    end
end
    