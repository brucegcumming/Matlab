function SetFigureName(F, name)


if ~isfigure(F)
    return;
end

set(F,'Name', name);
drawnow;
