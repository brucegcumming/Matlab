function WinFront(tag)
% WinFront(tag)
% bring to the front all figures with tags listed in structre tag.
f = fields(tag);
for j = 1:length(f)
    it = findobj('Tag',tag.(f{j}));
    if ~isempty(it)
        figure(it);
    end
end
