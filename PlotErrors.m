function D = PlotErrors(varargin)
%PlotErrors(name...   Show errors and warnings for and expt directory

parent = 0;
D.program = 'PlotErrors';
D.interactive = 0;
j = 1;
if ischar(varargin{1})
    if isdir(varargin{1})
        path = varargin{1};
    else
        [path, filename] = fileparts(varargin{1});
    end
    D.filename = dir2name(varargin{1},'filename');
    D.ErrorFile = [path '/Errors.mat'];
    D.penlabel = D.filename;
    Expts = ReadExptDir(path);
    D.Errors = ShowExptErrs(Expts);
    j = 2;
end
fontsize = 18;

while j <= length(varargin)
    if strncmpi(varargin{j},'checkandsave',10)
        if ~exist(D.ErrorFile) && ~isempty(D.Errors)
            Errors = D.Errors;
            save(D.ErrorFile,'Errors');
        end
        return;
    end
    j = j+1;
end

[F, isnew] = GetFigure('ErrorList');

if isnew
    cmenu = uicontextmenu;
    uimenu(cmenu,'Label','Edit','Callback',{@EditError, 'edit'});
    uimenu(cmenu,'Label','Delete (Del)','Callback',@DeleteError);
    uimenu(cmenu,'Label','Set Current Expt','Callback',{@EditError, 'setexp'});
    set(F,'DefaultUIControlFontSize',fontsize,'UIcontextmenu',cmenu);
    lst = uicontrol(F, 'Style','listbox','String', 'No Errors',...
        'Callback', {@HitErrorList},...
        'Tag', 'ErrorDisplay',...
        'KeyPressFcn', @KeyPressed,...
        'fontsize',fontsize,...
        'UIcontextmenu',cmenu,...
        'max', 2,...
        'units','norm', 'Position',[0 0.1 1 0.9]);
    D.lstgui = lst;
    a = uicontrol(F, 'Style','edit','String', '', ...
        'tag','AcceptError',...
        'fontsize',fontsize,...
        'units','norm', 'Position',[0 0 0.85 0.1],'callback',@ErrorText);
    setappdata(F,'ErrorListData',D);
    uicontrol(F, 'Style','Pushbutton','String', 'Close',...
        'Callback', {@Finished},...
        'fontsize',fontsize,...
        'units','norm', 'Position',[0.85 0 0.15 0.1]);
end
if isfield(D,'filename')
    set(F,'Name','Errors');
end
DisplayErrors(F);

function ErrorText(a,b)
   f = GetFigure(a);
   if get(f,'currentcharacter') %Hit return
       AddErrorList(a,b,'refresh');
   end
   
function DisplayErrors(F, varargin)
D = getappdata(F,'ErrorListData');
setpos = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'setpos',6)
        setpos = 1;
    end
    j = j+1;
end

    nc = 0;
    ne = 0;
    gothemisphere = 0;
    str = {};
    ts = [];
    strid = [];
    str{1} = D.filename;
    for j = 1:length(D.Errors)
        C = D.Errors{j};
        if isfield(C,'accepted') && ~isempty(C.accepted)
            chr = sprintf('(%s)*',C.user);
        else
            chr = '';
        end
        if isfield(C,'p')
            pstr = sprintf('P%d',C.p);
        else
            pstr = '';
        end
        str{j+1} = sprintf('%sE%d%s %s',chr,C.eid,pstr,C.s);
        strid(j+1) = j;
        if strcmp(D.Errors{j}.type,'cluster')
            nc = nc +1;
            linelist(j,:) = [0 nc];
        else
            ne = ne+1;
            linelist(j,:) = [ne 0];
        end        
    end
    setappdata(F,'linelist',linelist);
    lst = findobj(F,'Tag','ErrorDisplay');
    oldstrid = get(lst,'userdata');
    if ~isempty(str)
        set(lst,'String',str,'userdata',strid);
        if isempty(oldstrid) || setpos || get(lst,'value') > length(str) %new window or Errors
            set(lst,'value',length(str));
        end
    end
    

 function Finished(a,b)
   F = GetFigure(a);
   close(F);
        
function LoadErrors(name)

function DoOptions(a,b)
F = GetFigure(a);
D = getappdata(F,'ErrorListData');
f = get(a,'Tag');
if strcmp(f,'plotpen')
    ReadPen(FindPenLog(fileparts(D.ErrorFile)));
end
if strcmp(f,'errors')
    PlotErrors(fileparts(D.ErrorFile));    
end


function SetOptions(a,b)
F = GetFigure(a);
D = getappdata(F,'ErrorListData');
f = get(a,'Tag');
onoff = {'off' 'on'};
if isfield(D.show,f)
D.show.(f) = ~D.show.(f); 
set(a,'checked',onoff{1+D.show.(f)});
end
setappdata(F,'ErrorListData',D);
DisplayErrors(F);

function AddErrorList(a,b, type)
F = get(a,'parent');
D = getappdata(F,'ErrorListData');

%reload from disk then modify
[E, details] = ShowExptErrs(D.ErrorFile,'silent');
if ~isempty(E)
    D.Errors = E;
else %use app data
%    Errors = [];
%    Tagged = [];
end

a = findobj(F,'Tag','AcceptError');
txt = get(a,'string');

if length(txt) > 0

DATA = GetDataFromFig(F);
D = AddError(D, txt, DATA);
set(a,'string','');
SaveErrorList(D);
setappdata(F,'ErrorListData',D);
end

DisplayErrors(F);
%set(a,'value',n);


function SaveErrorList(D)
    Errors = {};
for j = 1:length(D.Errors)
    if strcmp(D.Errors{j}.type,'matfile')
        Errors{end+1} = D.Errors{j};
    end
end
if ~isempty(Errors)
    save(D.ErrorFile,'Errors');
end


function C = AddError(D, txt, DATA,  varargin)
Errors = D.Errors;
line = get(D.lstgui,'value');
strid = get(D.lstgui,'UserData');
n = strid(line);
if n > 0 && n <= length(Errors)
    Errors{n}.user = GetUserName;
    Errors{n}.usertime = now;
    Errors{n}.accepted = txt;
end

C = D;
C.Errors = Errors;

function KeyPressed(src, ks)
if sum(strcmp(ks.Key,{'delete' 'backspace'})) 
    DeleteError(src);
end
    
function EditError(a,b, fcn)
F = GetFigure(a);
a = findobj(F,'tag','ErrorDisplay');
str = get(a,'string');
line = get(a,'value');
strid = get(a,'userdata');
D = getappdata(F,'ErrorListData');

if strcmp(fcn,'edit')
    l = strid(line);
    l = l(l > 0);
    if ~isempty(l)
        load(D.ErrorFile);
        a = findobj(F,'tag','AcceptError');
        set(a,'string',D.Errors(l(1)).Error);
    end
elseif strcmp(fcn,'setexp')
    l = strid(line);
    ex = abs(l(l < 0));
    b = l(l>0);
    if ~isempty(b)
        b = D.Errors(b).exptno;
    end
    DATA = get(F,'UserData');
    DATA.currentexpt = [ex b];
    if sum(DATA.currentexpt)
        str = sprintf('%s:E%s',D.filename,sprintf('%d,',DATA.currentexpt));
    else
        str = D.filenamel;
    end
    set(F,'UserData',DATA,'Name',str);
end


function DeleteError(a,b)
%Deletes Error currently selected in text list
F = get(a,'parent');
a = findobj(F,'tag','ErrorDisplay');
line = get(a,'value');
strid = get(a,'userdata');
D = getappdata(F,'ErrorListData');
load(D.ErrorFile);
id = setdiff(1:length(Errors),strid(line));
Errors = Errors(id);
BackupFile(D.ErrorFile);
save(D.ErrorFile,'Errors','Tagged');
D.Errors = Errors;
setappdata(F,'ErrorListData',D);
if line > 1
    set(a,'value',line-1);
end
DisplayErrors(F);

function HitErrorList(a,b)
line = get(a,'value');
strid = get(a,'userdata');
F = get(a,'parent');
txtui = findobj(F,'tag','AcceptError');
D = getappdata(F,'ErrorListData');
C = [];
for j = 1:length(line)
    if strid(line(j)) >0 && strid(line(j)) < 1000
        C = D.Errors{strid(line(j))};
        if isfield(C,'accepted') && ~isempty(C.accepted)
            set(txtui,'string',C.accepted);
        else
            set(txtui,'string','');            
        end
    end
end
if isfield(C,'program')
    fprintf(' from %s',C.program);
end
