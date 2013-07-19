function SetFigPos(Figpos, tag)
%SetFigPos(X, tag)
%Sets figure size an location for figure whose tag is tag
%Data used to determine figure pos depends on X
%If X is a figure, getappdata(X,'Figpos');
%If X is strucutre with a field name matching tag, then 
%thie field of the structure is used
%if X is a typical figure data structure withe the field 'toplevel,
%then getappdata(X.toplevel,'Figpos');
if isfigure(Figpos)
    Figpos = getappdata(Figpos,'Figpos');
elseif isfield(Figpos,'toplevel')
    Figpos = getappdata(Figpos.toplevel,'Figpos');
end
if isfield(Figpos,tag)
    it = findobj('type','figure','Tag',tag);
    if length(it) == 1
        set(it,'position',Figpos.(tag));
    end
end
