function [success] = SetCheck(tag, value, varargin)
%SetCheck(tag, value)
%sets a checkbox with Tag = tag to value 
%if varargin includes a figure handle, search for Tag restricted to that
%figure.

if nargin == 3 & isfigure(varargin{1})
    it = findobj(varargin{1},'Tag',tag,'style','checkbox');
else    
    it = findobj('Tag',tag,'style','checkbox');
end

if length(it) == 1
    set(it,'value',value);
    success = 1;
elseif length(it) > 1
    fprintf('Tag %s has %d matches\n',tag, length(it));
    success = 0;
else
    success = 0;
end

