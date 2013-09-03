function SetFigPos(Figpos, tag)
%SetFigPos(X, tag)
%Sets figure size an location for figure whose tag is tag
%Data used to determine figure pos depends on X
%If X is a figure, getappdata(X,'Figpos');
%If X is strucutre with a field name matching tag, then 
%thie field of the structure is used
%if X is a typical figure data structure withe the field 'toplevel,
%then getappdata(X.toplevel,'Figpos');
toplevel = 0;
forcepos = 0;
if isfigure(Figpos)
    toplevel = Figpos;
    Figpos = getappdata(Figpos,'Figpos');
elseif isfield(Figpos,'toplevel')
    toplevel = Figpos.toplevel;
    Figpos = getappdata(Figpos.toplevel,'Figpos');
end
    it = findobj('type','figure','Tag',tag);
    if isempty(it)
        return;
    end
if isfield(Figpos,tag)
    go = 1;
    if length(Figpos.(tag)) > 4 && Figpos.(tag)(5) == 1 %already set once - user mmay have moved
        go = forcepos;
    end
    if length(it) == 1 && go
        set(it,'position',Figpos.(tag)(1:4));
        Figpos.(tag)(5) = 1;
    end
else
   Figpos.(tag)  = get(it(1),'position');
end
if isfigure(toplevel)
    setappdata(toplevel,'Figpos',Figpos);
end
