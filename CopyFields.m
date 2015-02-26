function to = CopyFields(to, from, varargin)
%to = CopyFields(to, from)
%Copy fields in structre from to struture to
%to = CopyFields(to, from, fields)
%      where fields is a cell array of strings, copies the field s
%named in fields
%to = CopyFields(to, from, field1, field2, field3) also works
%
%CopyFields(to, from, ..., '-noempty'
%    only copies fields in from that are not empty
%see also CopySFields to copy safely into a strcture element
j = 1;
useall = 1;
newonly = 0;
f = {};

while j <= length(varargin)
    if iscellstr(varargin{j})
        f = varargin{j};
    elseif strcmp(varargin{j},'-noempty')
        useall = 0;
    elseif strcmp(varargin{j},'-ifnew')
        newonly = 1;
    else 
        f = {f{:} varargin{j}};
    end
    j = j+1;
end

if isempty(from) || ~isstruct(from)
    return;
end

if isempty(f) %copy all fields
    f = fields(from);
end

for j = 1:length(f)
    if isfield(from,f{j}) && (useall || ~isempty(from.(f{j}))) && (newonly == 0 || ~isfield(to,f{j}))
        to.(f{j}) = from.(f{j});
    end
end
