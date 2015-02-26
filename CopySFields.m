function to = CopyFields(to, en, from, varargin)
%to = CopySFields(to, element, from) fopy fileds to structure element
%copyies files in from yo to(element N.B. TO is the struct, not one element)
%Copy fields in structre from to struture to
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

if en == 0 | strcmp(en,'new') %adda ne element
    en = length(to)+1;
end
if isempty(f) %copy all fields
    f = fields(from);
end

for j = 1:length(f)
    if isfield(from,f{j})
        to(en).(f{j}) = from.(f{j});
    end
end