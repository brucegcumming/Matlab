function Figpos = GetFigPos(parent, varargin)
% GetFigPos(parent, varargin) find current locations of all child figures
% and then records these in appdata 'Figpos' in the parent
Figpos = getappdata(parent,'Figpos');
f = findobj(get(0,'children'),'flat','type', 'figure');
for j = 1:length(f)
    p = getappdata(f(j), 'ParentFigure');
    if p == parent
        tag = get(f(j),'tag');
        tag = regexprep(tag,'/','');
        if isvarname(tag)
            Figpos.(tag) = get(f(j),'position');
        else
            Figpos.(genvarname(tag)) = get(f(j),'position');
        end
    end
end
tag = get(parent,'Tag');
Figpos.(tag) = get(parent,'position');
setappdata(parent,'Figpos', Figpos);