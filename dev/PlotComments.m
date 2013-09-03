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
D.show.program  =0;
D.show.time = 0;
D.show.exptprobe = 0;
DATA = [];
addline = [];

if ischar(varargin{1})
    if isdir(varargin{1})
        path = varargin{1};
    else
        [path, filename] = fileparts(varargin{1});
    end
    [D.ExptComments, D.ExptList] = BuildExptComments(varargin{1});
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
    elseif strncmpi(varargin{j},'addone',5)
        D.closeonadd = 1;
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
        F = GetFigure('Comments',parent);
        end
        DATA = get(parent,'UserData');
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
        F = GetFigure('Comments');
    end
    if D.interactive
        set(F,'DefaultUIControlFontSize',fontsize);

        setappdata(F,'CommentData',D);
        lst = uicontrol(F, 'Style','listbox','String', 'No Comments',...
            'Callback', {@HitCommentList},...
            'Tag', 'CommentDisplay',...
            'KeyPressFcn', @KeyPressed,...
            'fontsize',fontsize,...
            'units','norm', 'Position',[0 0.1 1 0.9]);  
        a = uicontrol(F, 'Style','edit','String', '', ...
            'tag','NewComment',...
            'fontsize',fontsize,...
            'units','norm', 'Position',[0 0 0.85 0.1]);  
        if D.closeonadd == 0
            set(a,'position',[0 0 0.7 0.1]);
            uicontrol(F, 'Style','Pushbutton','String', 'Add',...
                'Callback', {@AddCommentList, 'refresh'},...
                'fontsize',fontsize,...
                'units','norm', 'Position',[0.7 0 0.15 0.1]);
        end
        uicontrol(F, 'Style','Pushbutton','String', 'Done',...
            'Callback', {@AddCommentList, 'close'},...
            'fontsize',fontsize,...
            'units','norm', 'Position',[0.85 0 0.15 0.1]);  
        hm = uimenu(F,'Label','Show');
        uimenu(hm,'Label','User','tag','user','Callback',@SetOptions);
        uimenu(hm,'Label','Time','tag','time','Callback',@SetOptions);
        uimenu(hm,'Label','Expt/Probe','tag','exptprobe','Callback',@SetOptions);
        uimenu(hm,'Label','Program','tag','program','Callback',@SetOptions);
        DisplayComments(F);
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
    str = {};
    if D.includeexptlist
        for j = 1:length(D.ExptList)
            str{j} = sprintf('E%d %s: %d trials',j,D.ExptList(j).name,D.ExptList(j).ntrials);
            strid(j) = 0;
        end
    end
    nc = length(str);
    %Cell array or struct? WEhn called on a grid dir, seems to be a strcut
    for j = 1:length(D.ExptComments)
        str{nc+j} = D.ExptComments(j).text;
        strid(nc+j) = 1000+j;
    end
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
        if D.show.program
            epstr = [D.Comments(j).program ' ' epstr];
        end
        str{j+nc} = [ epstr D.Comments(j).comment];
        strid(j+nc) = j;
    end
    end
    lst = findobj(F,'Tag','CommentDisplay');
    oldstrid = get(lst,'userdata');
    if ~isempty(str)
        set(lst,'String',str,'userdata',strid);
        if isempty(oldstrid) || setpos %new window
            set(lst,'value',length(str));
        end
    end
    
    
function LoadComments(name)

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
if isfield(DATA,'exptno')
    ex = DATA.exptno;
elseif isfield(DATA,'currentexpt')
    ex = DATA.currentexpt;
else
    ex = 0;
end
Comments(n).ex = ex;
Comments(n).exptno = ex;

if isfield(DATA,'currentprobe')
    Comments(n).p = DATA.currentprobe;
else
    Comments(n).p = 0;
end
if iscell(txt) %return in text win makes string a cell with blank
    Comments(n).comment = txt{1};
else
    Comments(n).comment = txt;
end
if isfield(DATA,'program')
Comments(n).program = DATA.program;
end
C = D;
C.Comments = Comments;

function KeyPressed(src, ks)
if strmatch(ks.Key,'delete') 
    DeleteComment(src);
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
if strid(line) >0 && strid(line) < 1000
C = D.Comments(strid(line));
fprintf('%s at %s: %s\n',C.user, datestr(C.time), C.comment);
end