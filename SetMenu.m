function SetMenuCheck(F, tag, value, varargin)

exclusive = 0;
j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'exclusive')
        exclusive = 1;
    end
    j = j+1;
end

% SetMenu(F, tag, value)
% Find a menu belonging to figure F, set it checked. 

onoff = {'off' 'on'};
if ischar(F) %tag for a figure  
    F = findobj('type','figure','tag',F);
    if isempty(F);
        return;
    end
end

if ishandle(tag) && isnumeric(value)
    if value > 0 && exclusive
        c = get(get(tag,parent),'children')
        set(c,'Checked','off')
    end
    set(it,'Checked',onoff{1+value});
else
    it = findobj(F,'Tag',tag);
    if length(it) == 1
        if exclusive
            c = get(it,'children')
            set(c,'Checked','off');
        end
        mit = findobj(it,'tag',value);
        set(mit,'Checked',onoff{2});
    end
end
