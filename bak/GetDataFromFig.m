function [DATA, fn] = GetDataFromFig(a, varargin)
%if a subrouting has been called from a uicontrol, this find
% UserData from the parent figure
fn = 0;
    if isstruct(a) %can call a callback routine manually, passing DATA
        DATA = a;
        return;
    end
        
    DATA = get(a,'UserData');

    b = a;
    while isempty(DATA) & ~isempty(b)
        b = get(b,'parent');
        if isfigure(b)
            fn = b;
        end
        DATA = get(b,'UserData');
    end
    if isfigure(DATA)
        DATA = get(DATA,'UserData');
    elseif isfield(DATA,'parentfig')
        DATA = get(DATA.parentfig,'UserData');
    elseif isfield(DATA,'toplevel')
        DATA = get(DATA.toplevel,'UserData');
    elseif isfield(DATA,'parentfigtag')
        DATA = get(findobj('Tag',DATA.parentfigtag,'Type','Figure'),'UserData');
    end

    