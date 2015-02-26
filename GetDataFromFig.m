function [DATA, fn] = GetDataFromFig(a, varargin)
%if a subrouting has been called from a uicontrol, this find
% UserData from the parent figure
fn = 0;
javaobj = 0;
    if isstruct(a) %can call a callback routine manually, passing DATA
        DATA = a;
        return;
    end

    if strncmp(class(a),'javahandle',8)
        a = get(a,'uipanel');
        javaobj = 1;
    end
    b = GetFigure(a,'noforce');
    if isappdata(b,'ParentFigure')
        b = getappdata(b,'ParentFigure');
        DATA = get(b,'UserData');
        return;
    elseif javaobj
        DATA = get(b, 'UserData');
        return;
    end
    DATA = get(a,'UserData');
    b = a;
    while isempty(DATA) & ~isempty(b)
        if isfigure(b)
            fn = b;
            if isempty(DATA)
                p = getappdata(b,'ParentFigure');
                if ~isempty(p)
                    DATA = get(p,'UserData');
                    break;
                end
            end
            b = get(b,'parent');
        else
            b = GetFigure(b,'noforce');
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

    