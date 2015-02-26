function D = PlotComments(varargin)
%Plotcomments(name...   Shows comments for an experimen/fiel/directory, and
%allows additions
%PlotCommnets(D, 'addhidden', text)  adds comment text to the file, with no
%GUI
%
%Uses two files  ExptComments.mat, which is generated from the .smr files
%and Comments.mat which are manually added comments
parent = 0;
D.program = 'PlotComments';
D.closeonadd = 0;
D.includeexptlist = 1;
D.interactive = 1;
D.show.user = 0;
D.show.program  = 0;
D.show.time = 0;
D.show.exptprobe = 0;
D.show.depth = 0;
DATA = [];
addline = [];
filename = [];

if ischar(varargin{1})
    if isdir(varargin{1})
        path = varargin{1};
    else
        [path, filename] = fileparts(varargin{1});
    end
    D.filename = dir2name(varargin{1},'filename');
    [D.ExptComments, D.ExptList, D.CellInfo] = BuildExptComments(varargin{:});
    cname = [varargin{1} '/Comments.mat'];
    D.CommentFile = cname;
    if exist(cname)
        load(cname);
        D.Comments = Comments;
        D.Tagged = Tagged;
    else
        D.Comments = [];
        D.Tagged = [];
    end
    for j = 1:length(D.Comments)
        if iscell(D.Comments(j).comment)
            D.Comments(j).comment = D.Comments(j).comment{1};
        end
    end
    j = 2;
else
    if isfield(varargin{1},'Comments') && isfield(varargin{1},'CommentFile')
        D = varargin{1};
    end
    j = 1;
end
while j <= length(varargin)
    if isfield(varargin{j},'text') && isfield(varargin{j},'time')
        D.Comments = varargin{j};
    elseif isfield(varargin{j},'name') && isfield(varargin{j},'start')
        D.ExptList = varargin{j};
    elseif isfield(varargin{j},'progname') && isfield(varargin{j},'toplevel')
        DATA = varargin{j};
    elseif strncmpi(varargin{j},'exptno',4)
        j = j+1;
        DATA.exptno = varargin{j};
    elseif strncmpi(varargin{j},'probe',4)
        j = j+1;
        DATA.currentprobe = varargin{j};
    elseif strncmpi(varargin{j},'addone',5)
        D.closeonadd = 1;
    elseif strncmpi(varargin{j},'loadonly',5)
        return;
    elseif strncmpi(varargin{j},'noninteractive',5)
        D.interactive = 0;
    elseif strncmpi(varargin{j},'addhidden',5)
        D.interactive = 0;
        j = j+1;
        addline = varargin{j};
    elseif strncmpi(varargin{j},'parent',5)
        j = j+1;
        parent = varargin{j};
    end        
    j = j+1;
end
fontsize = 18;

if ~isempty(addline) %just add a comment without GUI
    D = AddComment(D,addline, DATA);
    SaveComments(D);
    return;
end
    if parent > 0
        if D.interactive
            [F, isnew] = GetFigure('Comments',parent,'front');
        end
        X = get(parent,'UserData');
        if ~isempty(X)
            DATA = X;
        end
%        D.closeonadd = 0;
        if isfield(DATA,'gui')
            if isfield(DATA.gui,'fontsize')
                fontsize = DATA.gui.fontsize;
            elseif isfield(DATA.gui,'FontSize')
                fontsize = DATA.gui.FontSize;
            end
        end
        if isfield(DATA,'progname')
            D.program = DATA.progname;
        end
    else
        [F, isnew] = GetFigure('Comments');
        D.program = 'PlotComments';
        if ~isempty(DATA)
            if ~isfield(DATA,'toplevel')
                DATA.toplevel = F;
            end
            set(F, 'UserData',DATA);
        end
    end    
    if D.interactive 
        if ~isempty(filename)
            D.filename = filename;
        end
        str = dir2name(D.filename,'filename');
        if isfield(D,'CellInfo') && isfield(D.CellInfo,'pen')
            str = [str sprintf('   Pen%d: %.1f,%.1f',D.CellInfo.pen)];
        end
        D.penlabel = str;
        if isnew

            cmenu = uicontextmenu;
            uimenu(cmenu,'Label','Edit','Callback',{@EditComment, 'edit'});
            uimenu(cmenu,'Label','Delete (Del)','Callback',@DeleteComment);
            uimenu(cmenu,'Label','Set Current Expt','Callback',{@EditComment, 'setexp'});
            set(F,'DefaultUIControlFontSize',fontsize,'UIcontextmenu',cmenu);
            lst = uicontrol(F, 'Style','listbox','String', 'No Comments',...
                'Callback', {@HitCommentList},...
                'Tag', 'CommentDisplay',...
                'KeyPressFcn', @KeyPressed,...
                'fontsize',fontsize,...
                'UIcontextmenu',cmenu,...
                'max', 2,...
                'units','norm', 'Position',[0 0.1 1 0.9]);
            D.lstgui = lst;
            a = uicontrol(F, 'Style','edit','String', '', ...
                'tag','NewComment',...
                'fontsize',fontsize,...
                'units','norm', 'Position',[0 0 0.85 0.1],'callback',@CommentText);
            if D.closeonadd == 0
                set(a,'position',[0 0 0.7 0.1]);
                a = uicontrol(F, 'Style','Pushbutton','String', 'Add',...
                    'Callback', {@AddCommentList, 'refresh'},...
                    'fontsize',fontsize,...
                    'units','norm', 'Position',[0.7 0 0.15 0.1]);
                D.addbutton = a;                
            end
            setappdata(F,'CommentData',D);
            uicontrol(F, 'Style','Pushbutton','String', 'Done',...
                'Callback', {@AddCommentList, 'close'},...
                'fontsize',fontsize,...
                'units','norm', 'Position',[0.85 0 0.15 0.1]);
            hm = uimenu(F,'Label','Show');
            uimenu(hm,'Label','User','tag','user','Callback',@SetOptions);
            uimenu(hm,'Label','Time','tag','time','Callback',@SetOptions);
            uimenu(hm,'Label','Expt/Probe','tag','exptprobe','Callback',@SetOptions);
            uimenu(hm,'Label','Electrode Depth','tag','depth','Callback',@SetOptions);
            uimenu(hm,'Label','Program','tag','program','Callback',@SetOptions);
            uimenu(hm,'Label','Plot Penetration','tag','plotpen','Callback',@DoOptions);
            uimenu(hm,'Label','Errors','tag','errors','Callback',@DoOptions);
        else %check if file changed while window up
            X = getappdata(F,'CommentData');
            if ~strcmp(X.filename,D.filename)
                setappdata(F,'CommentData',D);
            end
        end
            if isfield(D,'filename')
                set(F,'Name',D.penlabel);
            end
        DisplayComments(F);
    end
 
function CommentText(a,b)
   f = GetFigure(a);
   if get(f,'currentcharacter') %Hit return
       AddCommentList(a,b,'refresh');
   end
   
function DisplayComments(F, varargin)
D = getappdata(F,'CommentData');
setpos = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'setpos',6)
        setpos = 1;
    end
    j = j+1;
end

    nc = 0;
    gothemisphere = 0;
    str = {};
    ts = [];
    str{1} = D.penlabel;
    if D.includeexptlist
        for j = 1:length(D.ExptList)
            epstr = [];
            dstr = [' at ' datestr(D.ExptList(j).start,'HH:MM')];
            if D.show.depth && isfield(D.ExptList,'depth')
                epstr = [epstr ' ' sprintf('ed=%.3f',mean(D.ExptList(j).depth))];
            end
            str{j+1} = sprintf('E%d %s: %d trials%s %s',j,D.ExptList(j).name,D.ExptList(j).ntrials,dstr,epstr);
            ts(j+1) = D.ExptList(j).start;
            if isfield(D.ExptList,'suffix')
                ex = D.ExptList(j).suffix;
            else
                ex = j;
            end
            strid(j+1) = -ex;
        end
    end
    nc = length(str);
    %Cell array or struct? WEhn called on a grid dir, seems to be a strcut
    nl = 0;
    for j = 1:length(D.ExptComments)
        epstr = [];
        if D.show.depth && isfield(D.ExptList,'depth')
            epstr = [epstr ' ' sprintf('ed=%.3f',mean(D.ExptComments(j).depth))];
        end
        dstr = [' at ' datestr(D.ExptComments(j).time,'HH:MM')];
        go = 1;
        if strcmp(D.ExptComments(j).text,'RightHemisphere') && gothemisphere
            go = 0;
        end
        if findstr(D.ExptComments(j).text,'RightHemisphere')
            gothemisphere = 1;
        end
        if go
            nl = nl+1;
            str{nc+nl} = [D.ExptComments(j).text dstr epstr];
            ts(nc+nl) = D.ExptComments(j).time;
            strid(nc+nl) = 1000+j;
        end
    end
    [a,b] = sort(ts);
    str = str(b);
    strid = strid(b);
    nc = length(str);
    if isfield(D,'Comments')
    for j = 1:length(D.Comments)
        if D.Comments(j).exptno > 0 & D.show.exptprobe
            epstr = sprintf('E%dP%d ',D.Comments(j).ex(1),D.Comments(j).p);
        else
            epstr = '';
        end
        if D.show.user
            epstr = [D.Comments(j).user ' ' epstr];
        end
        if D.show.time
            epstr = [datestr(D.Comments(j).time) ' ' epstr];
        end
        if D.show.program && isfield(D.Comments,'program')
            epstr = [D.Comments(j).program ' ' epstr];
        end
        if D.show.depth && isfield(D.Comments,'depth')
            epstr = [sprintf('%.3f') D.Comments(j).depth ' ' epstr];
        end
        str{j+nc} = [ epstr D.Comments(j).comment];
        strid(j+nc) = j;
    end
    end
    lst = findobj(F,'Tag','CommentDisplay');
    oldstrid = get(lst,'userdata');
    if ~isempty(str)
        set(lst,'String',str,'userdata',strid);
        if isempty(oldstrid) || setpos || get(lst,'value') > length(str) %new window or comments
            set(lst,'value',length(str));
        end
    end
    
    
function LoadComments(name)

function DoOptions(a,b)
F = GetFigure(a);
D = getappdata(F,'CommentData');
f = get(a,'Tag');
if strcmp(f,'plotpen')
    ReadPen(FindPenLog(fileparts(D.CommentFile)));
end
if strcmp(f,'errors')
    PlotErrors(fileparts(D.CommentFile));    
end


function SetOptions(a,b)
F = GetFigure(a);
D = getappdata(F,'CommentData');
f = get(a,'Tag');
onoff = {'off' 'on'};
if isfield(D.show,f)
D.show.(f) = ~D.show.(f); 
set(a,'checked',onoff{1+D.show.(f)});
end
setappdata(F,'CommentData',D);
DisplayComments(F);

function AddCommentList(a,b, type)
F = get(a,'parent');
D = getappdata(F,'CommentData');
if exist(D.CommentFile,'file')
    load(D.CommentFile);
D.Comments = Comments;
D.Tagged = Tagged;
else %use app data
%    Comments = [];
%    Tagged = [];
end

a = findobj(F,'Tag','NewComment');
txt = get(a,'string');

if length(txt) > 0

DATA = GetDataFromFig(F);
D = AddComment(D, txt, DATA);
set(a,'string','');
SaveComments(D);
setappdata(F,'CommentData',D);
end

if D.closeonadd || strcmp(type,'close')
    close(F);
else
DisplayComments(F,'setpos');
%set(a,'value',n);
end

function SaveComments(D)
Comments = D.Comments;
Tagged = D.Tagged;
save(D.CommentFile,'Comments', 'Tagged');


function C = AddComment(D, txt, DATA,  varargin)
Comments = D.Comments;
n = length(Comments) + 1;
Comments(n).user = GetUserName;
Comments(n).time = now;
cx = sscanf(txt,'E%fP%f');
if isfield(DATA,'exptno')
    ex = DATA.exptno;
elseif isfield(DATA,'currentexpt')
    ex = DATA.currentexpt;
else
    line = get(D.lstgui,'value');
    strid = get(D.lstgui,'UserData');
    ex = strid(line);
    ex = ex(ex<0);
    if isempty(ex) 
        ex = 0;
        if ~isempty(cx) && cx(1) > 0 %Read from line
            ex = cx(1);
        end
    else
        ex = -ex;
    end
end
Comments(n).ex = ex;
Comments(n).exptno = ex;

if isfield(DATA,'currentprobe')
    Comments(n).p = DATA.currentprobe;
else
    if length(cx) > 1 && cx(2) > 0 %Read from line
        Comments(n).p = cx(2);
    else
        Comments(n).p = 0;
    end
end
if iscell(txt) %return in text win makes string a cell with blank
    Comments(n).comment = txt{1};
else
    Comments(n).comment = txt;
end
if isfield(D,'program')
    Comments(n).program = D.program;
elseif isfield(DATA,'program')
    Comments(n).program = DATA.program;
end
C = D;
C.Comments = Comments;

function KeyPressed(src, ks)
if sum(strcmp(ks.Key,{'delete' 'backspace'})) 
    DeleteComment(src);
end
    
function EditComment(a,b, fcn)
F = GetFigure(a);
a = findobj(F,'tag','CommentDisplay');
str = get(a,'string');
line = get(a,'value');
strid = get(a,'userdata');
D = getappdata(F,'CommentData');
if strcmp(fcn,'edit')
    l = strid(line);
    l = l(l > 0);
    if ~isempty(l)
        load(D.CommentFile);
        a = findobj(F,'tag','NewComment');
        set(a,'string',D.Comments(l(1)).comment);
    end
elseif strcmp(fcn,'setexp')
    l = strid(line);
    ex = abs(l(l < 0));
    b = l(l>0);
    if ~isempty(b)
        b = D.Comments(b).exptno;
    end
    DATA = get(F,'UserData');
    DATA.currentexpt = [ex b];
    if sum(DATA.currentexpt)
        str = sprintf('%s:E%s',D.penlabel,sprintf('%d,',DATA.currentexpt));
    else
        str = D.penlabel;
    end
    set(F,'UserData',DATA,'Name',str);
end


function DeleteComment(a,b)
%Deletes comment currently selected in text list
F = get(a,'parent');
a = findobj(F,'tag','CommentDisplay');
line = get(a,'value');
strid = get(a,'userdata');
D = getappdata(F,'CommentData');
load(D.CommentFile);
id = setdiff(1:length(Comments),strid(line));
Comments = Comments(id);
BackupFile(D.CommentFile);
save(D.CommentFile,'Comments','Tagged');
D.Comments = Comments;
setappdata(F,'CommentData',D);
if line > 1
    set(a,'value',line-1);
end
DisplayComments(F);

function HitCommentList(a,b)
line = get(a,'value');
strid = get(a,'userdata');
F = get(a,'parent');
D = getappdata(F,'CommentData');
C = [];
for j = 1:length(line)
    if strid(line(j)) >0 && strid(line(j)) < 1000
        C = D.Comments(strid(line(j)));
        fprintf('%s: %s at %s',C.comment, C.user, datestr(C.time));
    end
end
if isfield(C,'program')
    fprintf(' from %s',C.program);
end
fprintf('\n');