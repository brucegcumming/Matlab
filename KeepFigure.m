function F = KeepFigure(a, varargin)
%KeepFigure(figure,  ...) Changes a figure Tag so that subsequent calls to GetFigure do not
%return this figure.  
%figure can be a figure handle or a Tag.  
%KeepFigure(figure,  'detach'...) reomves any association with a parent, so
%that Closing the the parent will not close the kept figure
if isfigure(a)
    F = a;
else
    F= GetFigure(a);
end

detatch = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'detatch',4)
        detatch = 1;
    end
    j = j+1;
end

tag = get(F,'Tag');
set(F,'tag', ['Keep' tag]);
name = get(F,'name');
set(F,'name', ['Keep' name]);
setappdata(F,'KeepFigure',tag); %record original tag
if detatch
    if isappdata(F,'ParentFigure')
        rmappdata(F,'ParentFigure');
    end
end