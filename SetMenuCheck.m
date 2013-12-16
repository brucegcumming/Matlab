function SetMenuCheck(F, tag, value, varargin)
% Turns on/off checked proerty in menu itesm
% SetMenuCheck(F, tag, value) finds by tag
% SetMenuCheck(handle, value) Just sets one item
% SetMenuCheck(handle, 'exclusive' value) or
% SetMenuCheck(F, tag, value,'exclusive') 
%           unsets other members of the same menu
% Find a menu belonging to figure F, set it checked/unchecked according to value.
% SetMenu(F, tag, value,'exclusive') turns off other items in the same
% menu, if value > 0
%Set Check marks on menus
%if F is a handle to a menu object, means this has just been selected
% SetMenucheck(a, 'exclusive') turns off checks for other members of menu
%so turn this on
%Otherwised
exclusive = 0;
j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'exclusive')
        exclusive = 1;
    end
    j = j+1;
end


onoff = {'off' 'on'}; %in case value > 1
%would be nice to have second check type
%or could have value = value > 0?? 
if ischar(F) %tag for a figure  
    F = findobj('type','figure','tag',F);
    if isempty(F);
        return;
    end
end

if ishandle(F) & ~isfigure(F)
    
    if nargin > 1 && strncmp(tag,'exclusive',5)
        exclusive = 1;
    end
    if exclusive
        c = get(get(F,'parent'),'children');
        set(c,'Checked','off')
    end
    if nargin ==2 && length(tag) == 0 && tag == 0
        set(F,'Checked','off');
    else
        set(F,'Checked','on');
    end
elseif ischar(tag)
    it = findobj(F,'Tag',tag);
    if length(it) == 1
        if exclusive
            c = get(it,'children');
            set(c,'Checked','off');
        end
        if ischar(value)
            mit = findobj(it,'tag',value);
            set(mit,'Checked',onoff{2});
        elseif isnumeric(value)
            value = value > 0;
            set(it,'Checked',onoff{1+value});
        end
    end
elseif ishandle(tag) && isnumeric(value)
    if value > 0 && exclusive
        c = get(get(tag,parent),'children')
        set(c,'Checked','off')
    end
    set(it,'Checked',onoff{1+value});
end
