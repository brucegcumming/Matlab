function F = GetParentFigure(a, varargin)

F = GetFigure(a);
if isappdata(F,'ParentFigure')
    F = getappdata(F,'ParentFigure');
else
    D = get(F,'UserData');
    if isfield(D,'toplevel');
        F = D.toplevel;
    end
end