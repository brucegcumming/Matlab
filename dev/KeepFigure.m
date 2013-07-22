function KeepFigure(a,b)
%Changes a figure Tag so that subsequent calls to GetFigure do not
%return this figure.  

if isfigure(a)
    F = a;
else
    F= GetFigure(a);
end

tag = get(F,'Tag');
set(F,'tag', ['Keep' tag]);
