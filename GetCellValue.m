function x = GetCellValue(C, n, f)
%GetCellValue(C, n, field) returns C{n}.field IF
%it exists
x = NaN;
if n <=0
    return;
end
if length(C) >= n && isfield(C{n},f)
    x = C{n}.(f);
end