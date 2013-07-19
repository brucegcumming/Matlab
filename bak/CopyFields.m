function to = CopyFields(to, from, varargin)
%to = CopyFields(to, from)
%Copy fields in structre from to struture to
%to = CopyFields(to, from, fields)
%      where fields is a cell array of strings, copies the field s
%named in fields

j = 1;
f = {};

while j <= length(varargin)
    if iscellstr(varargin{j})
        f = varargin{j};
    end
    j = j+1;
end

if isempty(from)
    return;
end

if isempty(f) %copy all fields
    f = fields(from);
end

for j = 1:length(f)
    if isfield(from,f{j})
        to.(f{j}) = from.(f{j});
    end
end