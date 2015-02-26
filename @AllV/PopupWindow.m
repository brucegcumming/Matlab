function PopupWindow(toplevel, tag, label, callback, varargin)

DATA = get(toplevel,'UserData');
if isfield(DATA,'tag') && isfield(DATA.tag,'popup')
    wtag = DATA.tag.popup;
else
    wtag = 'popup';
end

if nargin == 1
    tag = '';
    label = '';
    callback = [];
end

truetag = tag;
truelabel = label;
if isempty(tag)
    tag = label;
end
tag = strrep(tag,'&','');
label = strrep(label,'&','');
X = getappdata(toplevel,'PopupWindowData');
if ~iscell(X) 
    X = {};
    if nargin ==1 % dont popup if empty
    return;
    end
end
[F, isnew] = GetFigure(wtag,'parent',toplevel,'trackpos','setpos');
if isnew
    set(F,'menubar','none');
    sm = uimenu(F,'label','Remove','tag','RemoveMenu');
    sm = uimenu(F,'label','Add','tag','AddMenu');
    uimenu(sm,'label','Last Operation',...
        'callback',@AddCallback,...
        'tag','AddCallbackButton');
    if size(X,2) > 2
        a = CellToStruct(X(:,3));
        nr = 1+sum([a.on]);
        id = find([a.on]);
        if sum(strcmp(truelabel,X(id,2))) == 0
            nr = nr-1;
        end
    else
        nr = 1;
    end
    for j = 1:size(X,1)
        if X{j,3}.on
            it =FindChild(DATA.toplevel,'tag',X{j,1},'label',X{j,2});
            if ~isempty(it)
                X{j,3}.buttonid = uicontrol(F,'style','pushbutton','string',X{j,2},...
                    'callback',get(it(1),'Callback'),...
                    'tag',X{j,1},...
                    'units','normalized','position',[0.02 (j-1)/nr 0.96 0.96/nr]);
            end
            
        end
    end
end

it = findobj(F,'type','uicontrol','style','pushbutton','tag',tag);
if isempty(it)
    it = findobj(F,'type','uicontrol','style','pushbutton');
    nr = length(it)+1;
    for j = 1:length(it)
        set(it(j), 'units','normalized','position',[0.02 (j-1)/nr, 0.96, 0.96/nr]);
    end
    
    j = nr;
    if ~isempty(callback)
        state.on = 1;
    X{end+1,1} = truetag;
    X{end,2} = truelabel;
    state.buttonid = uicontrol(F,'style','pushbutton','string',tag,...
        'callback',callback,...
        'tag',tag,...
        'units','normalized','position',[0.02 (j-1)/nr 0.96 0.96/nr]);
    X{end,3} = state;
    end

    sm  = findobj(F,'type','uimenu','tag','RemoveMenu');
    uimenu(sm,'label',label,'tag',tag,'callback',{@RemoveButton, tag});
    setappdata(toplevel,'PopupWindowData',X);
end

function RemoveButton(a,b,tag)

F = GetFigure(a);
it = findobj(F,'type','uicontrol','style','pushbutton','tag',tag);
delete(it);
sm  = findobj(F,'type','uimenu','tag','RemoveMenu');
it = findobj(sm,'type','uimenu','tag',tag);
sm  = findobj(F,'type','uimenu','tag','AddMenu');
ait = findobj(sm,'type','uimenu','tag',tag);
if isempty(ait)
    args = CopyUiProperties(it(1),{'tag' 'label'});
    hm = uimenu(sm,args{:});
    tag = get(it,'tag');
    set(hm, 'callback', {@AddButton, tag},'userdata',get(it,'callback'));
end
delete(it);
parent = getappdata(F,'ParentFigure');
X = getappdata(parent,'PopupWindowData');
state = CellToStruct(X(:,3));
id = find([state.buttonid] == it);
if ~isempty(id)
    X{id,3}.on = 0;
end
setappdata(parent,'PopupWindowData',X);

function AddButton(a,b,tag)

F = GetFigure(a);
toplevel = getappdata(F,'ParentFigure');
sm = get(a,'parent');
%original callback stored in userdataeta
AllV.PopupWindow(toplevel, get(a,'tag'),get(a,'label'),get(a,'userdata'));


function AddCallback(a,b)

F = GetFigure(a);
DATA = GetDataFromFig(F);
a = DATA.guistate.lastguihandle;
if ishandle(a)
    AllV.PopupWindow(DATA.toplevel, get(a,'tag'),get(a,'label'),get(a,'callback'));
end

