function PlotComments(varargin)
%Plotcomments(name...   Shows comments for an experimen/fiel/directory, and
%allows additions
%
%Uses two files  ExptComments.mat, which is generated from the .smr files
%and Comments.mat which are manually added comments
parent = 0;
D.program = 'PlotComments';
D.closeonadd = 0;
D.includeexptlist = 1;

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
    j = 1;
end
while j <= length(varargin)
    if isfield(varargin{j},'text') && isfield(varargin{j},'time')
        D.Comments = varargin{j};
    elseif isfield(varargin{j},'name') && isfield(varargin{j},'start')
        D.ExptList = varargin{j};
    elseif strncmpi(varargin{j},'addone',5)
        D.closeonadd = 1;
    elseif strncmpi(varargin{j},'parent',5)
        j = j+1;
        parent = varargin{j};
    end        
    j = j+1;
end
fontsize = 18;

    if parent > 0
        F = GetFigure('Comments',parent);
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
    DisplayComments(F);

    
function DisplayComments(F)
D = getappdata(F,'CommentData');
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
        if D.Comments(j).exptno > 0
            epstr = sprintf('E%dP%d ',D.Comments(j).ex,D.Comments(j).p);
        else
            epstr = '';
        end
        str{j+nc} = [ epstr D.Comments(j).comment];
        strid(j+nc) = j;
    end
    end
    lst = findobj(F,'Tag','CommentDisplay');
    oldstrid = get(lst,'userdata');
    if ~isempty(str)
        set(lst,'String',str,'userdata',strid);
        if isempty(oldstrid) %new window
            set(lst,'value',length(str));
        end
    end
    
    
function LoadComments(name)

function AddCommentList(a,b, type)
F = get(a,'parent');
D = getappdata(F,'CommentData');
if exist(D.CommentFile,'file')
    load(D.CommentFile);
else
    Comments = [];
    Tagged = [];
end
a = findobj(F,'Tag','NewComment');
txt = get(a,'string');

if length(txt) > 0
n = length(D.Comments) + 1;

DATA = GetDataFromFig(F);

Comments(n).user = GetUserName;
Comments(n).time = now;
if isfield(DATA,'exptno')
    Comments(n).ex = DATA.exptno;
    Comments(n).exptno = DATA.exptno;
else
    Comments(n).ex = 0;
    Comments(n).exptno = 0;
end
if isfield(DATA,'currentprobe')
Comments(n).p = DATA.currentprobe;
else
    Comments(n).p = 0;
end
a = findobj(F,'Tag','NewComment');
Comments(n).comment = get(a,'string');
Comments(n).program = D.program;
set(a,'string','');

save(D.CommentFile,'Comments', 'Tagged');
D.Comments = Comments;
D.Tagged = Tagged;
setappdata(F,'CommentData',D);
end
if D.closeonadd || strcmp(type,'close')
    close(F);
else
DisplayComments(F);
end

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