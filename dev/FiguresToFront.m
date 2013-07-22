function FiguresToFront(tag, varargin)
%for all the tags conainted in the structure
%bring the corresponding figure window forward

f = fields(tag);
for j = 1:length(f)
    it = findobj('type','figure','Tag', tag.(f{j}));
    if length(it) == 1
        figure(it);
    end
end
        