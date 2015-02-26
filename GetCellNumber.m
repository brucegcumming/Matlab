function cell = GetCellNumber(name)
%Get a cell number from a string,Expt or Cluster
%not all implemented yet!


if isfield(name, 'Header') && isfield(name.Header,'cellnumber')
    cell = name.Header.cellnumber;
elseif ischar(name)
    cell = NumberFromString(name, 'cell');
end


