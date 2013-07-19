function mylegend(h, labels, varargin)
% mylegend(h, labels, varargin)
% calls legend with a list of handles and labels, but checks these
% first to avoid annoying legend crashes.
ilb = []; %in case labels is empty.
if isempty(labels)
    return;
end
for j = 1:length(labels)
    if ischar(labels{j})
        ilb(j) = 1;
    else
        ilb(j) = 0;
    end
end
id = intersect(find(h & ishandle(h)),find(ilb));
if length(id)
legend(h(id),{labels{id}}, varargin{:});
end

