function [success] = SetCheck(tag, value,  varargin)

success = 0;
if nargin == 3 & isfigure(varargin{1})
    it = findobj(allchild(varargin{1}),'flat','Tag',tag);
elseif nargin == 3 & double(varargin{1}) == 0
    return;
else
    it = findobj('Tag',tag);
end
if ~isempty(it)
    set(it(1),'value',value);
    success = 1;
end

