function [j, value] = GetOptionalArg(args,j, default)
% [j, value] = GetOptionalArg(args,j, default) get number from next arg in list,
% if available

if length(args) > j && isnumeric(args{j+1})
    j = j+1;
    value = args{j};
else
    value = default;
end