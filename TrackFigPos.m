function TrackFigPos(F, b, varargin)
%called when fiugres are closed to store current position


parent = getappdata(F,'ParentFigure');
if ~isempty(parent) && isfigure(parent)
    Figpos = getappdata(parent,'Figpos');
    tag = get(F,'tag');
    Figpos.(tag) = get(F,'position');
    setappdata(parent,'Figpos',Figpos);
end

j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'close')
        delete(F);
    end
    j = j+1;
end