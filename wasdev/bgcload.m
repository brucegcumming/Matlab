function varargout= bgcload(name)
%load matlab data, and include pathname
%also allows workaround for 2013b '/' bug on PC


if ispc && name(1) == '/'
    name(1) = '\';
end

a = load(name);
f = fields(a);

%now add some loadname depending on struct
if sum(strcmp('Clusters',f))
    for j = 1:length(a.Clusters)
        a.Clusters{j}.loadname = name;
    end    
end

for j = 1:length(f)
    varargout{j} = a.(f{j});
end
