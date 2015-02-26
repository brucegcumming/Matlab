function DATA = listadd(varargin)
%
%listadd(varagin)
%starter .m file for adding extra features to runlist
%listadd.m
fname = 'listadd';

listtag = 'RevCorrTop';
j = 1;
if isempty(findobj('Tag',listtag))
    if nargin 
        list = varargin{1};
    else
        list = 'rclist.mall';
    end
    j = 2;
    DATA = runlist(list,'Tag',listtag,'function',fname,varargin{2:end});
    DATA.fname = fname;
    DATA.plot.rcsmooth = 20;
    DATA.plot.simplify = 0;
    if ~isfield(DATA,'statefile')
        DATA.statefile = 'rundtDATA.mat';
    end

%
% Add Plots to the Main Plot Type Menu
    AddStrings('plottype','Xplot');
    AddStrings('ListType','Xlist');
    AddStrings('PopPlot','Xpop');
    AddStrings('extraplot','Xextra');
    AddStrings('xop','Xextra');

    uimenu(DATA.menus.mark,'Label','LoadDATA','Callback',[fname '(''load'');']);
    uimenu(DATA.menus.options,'Label','Relist','Callback',[fname '(''relist'',0);']);


    runlist('store',DATA);
end

DATA = runlist('getstate');

while j < nargin+1
    if strncmpi(varargin{j},'setid',5)
        DATA = doentry(DATA,DATA.id);
    elseif strncmpi(varargin{j},'setselect',5)
        AddSelect(DATA.cntrlbox,DATA);
    elseif strncmpi(varargin{j},'update',5)
        DATA.plot.rcsmooth = str2num(get(findobj('Tag','rcsmooth'),'string'));
        DATA.plot.simplify = get(findobj('Tag','Simplify'),'value');
        runlist('store',DATA);
        if DATA.plot.autoplot
            doentry(DATA,DATA.id);
        end
    elseif strncmpi(varargin{j},'DATA',4)
        j = j+1;
        DATA = varargin{j};
    elseif strncmpi(varargin{j},'popplot',4)
        DATA = PopPlot(DATA);
    elseif strncmpi(varargin{j},'touchpoint',6)
        j = j+2;
        DATA.id = varargin{j};
        GetFigure(DATA.tag.pop);
        PlotPoint(DATA,DATA.plot.allpts(:,DATA.id),DATA.id,2);
        fprintf('%s\n',DATA.fstrings{DATA.id});
        doentry(DATA,DATA.id);
    elseif strncmpi(varargin{j},'All',3)
        for j = 1:length(DATA.fstrings)
            res = PlotExpt(DATA.fstrings{j});
            DATA.results(j) = res;
            nu = strmatch('Uncorr',res.extras.label);
            lm = strmatch('Left',res.extras.label);
            rm = strmatch('Right',res.extras.label);
            blank = strmatch('Blank',res.extras.label);
            if ~iempty(lm) & ~isempty(rm) & ~isempty(blank)
                
            end
        end
    else
        fprintf('%s\n',varargin{j});
    end
    j = j+1;
end


function AddSelect(f, DATA)
%
% AddSelect pops up an additional set of options in a window
%

SPACE = 5;
VSPACE = 5;
cw = 16;
ch = 18;
bh = 18*DATA.wsc;
bp(1) = 0;
bp(2) = DATA.selectlast + bh;
bp(3) = 40;   
bp(1) = SPACE;
bp(4) = bh;
fname = DATA.fname;

  uicontrol(gcf,'Style', 'text','String','rcsdf','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
uicontrol(f,'Style', 'edit', 'Callback', [fname '(''update'')'],...
    'String', sprintf('%d',DATA.plot.rcsmooth), 'Tag', 'rcsmooth','Position', bp);

bp(1) = bp(1) +bp(3) + SPACE;
bp(3) = cw *6;
uicontrol(gcf,'Style', 'CheckBox','String','Simplify','Position', bp,...
      'Tag','Simplify','Callback',[fname '(''update'')'],'Value',DATA.plot.simplify);


function DATA = doentry(DATA, id)


argon{1} = 'legendpos';
argon{2} = DATA.plot.legendpos;

if strmatch(DATA.plot.onetype,'ACtuning')


end


function AddStrings(tag, varargin)

j = 1;
k = 1;
while j < nargin
    j = j+1;
end

mnu = findobj('Tag',tag);
if ~isempty(mnu)
    strs = cellstr(get(mnu,'String'));
    newstr = {strs{:}, varargin{k:end}};
    set(mnu,'String',newstr);
end
    

function PlotPoint(DATA, pos, id, type)
symbol = 'o';
if type == 1
    color = 'b';
    fillcolor = 'b';
else
    color = 'r';
    fillcolor = 'r';
end
plot(pos(1),pos(2),symbol,'buttondownfcn',[DATA.fname '(''touchpoint'',gcf,' num2str(id) ');'],'color',color,...
    'MarkerFaceColor',fillcolor);
hold on;

   
function DATA = PopPlot(DATA)

GetFigure(DATA.tag.pop);
idx = 1:length(DATA.fstrings);
%idx = find([cps.nmin] > 10 & [cps.proportion] > 0.2 & [cps.dprime] > DATA.plot.mindprime);
%limit plot to those things currently on list
idx = intersect(idx, DATA.listids);
val =  strmatch(DATA.plot.poptype,{'Xpop'});

for j = 1:length(idx)
    pos(1) = j;
    pos(2) = j;
    type = 1; % set this for sig/ns or whatever want symbols for
    PlotPoint(DATA,pos,idx(j),type);
    DATA.plot.allpts(:,idx(j)) = pos;
end
