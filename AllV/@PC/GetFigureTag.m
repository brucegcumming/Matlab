function tag = GetFigureTag(src)    a = src;    while ~isfigure(a) && a ~= 0        a = get(a,'Parent');    end    tag = get(a,'Tag');                     