function value = GetField(tag, varargin)

if nargin == 2 & isfigure(varargin{1})
it = findobj(varargin{1},'Tag',tag);
else    
it = findobj('Tag',tag);
end
if ~isempty(it) 
value = str2num(get(it(1),'string'));
else
value = NaN;
end

