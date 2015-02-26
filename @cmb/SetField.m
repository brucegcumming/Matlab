function SetField(parent, tag, value, varargin)
if isfigure(parent)
it = findobj(parent,'Tag', tag);
else
it = findobj('Tag', tag);
end
if double(it) > 0
set(it,'string',num2str(value));
end

