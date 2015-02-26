function [value, it] = GetShow(DATA,tag)
it = findobj(DATA.showid, 'Tag',tag);
if ~isempty(it) 
value = get(it(1),'value');
else
value = 0;
end

