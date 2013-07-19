function [out, n] = PlotExpts(list, varargin)

% PlotExpts(dir, ...)
% list all .mat files for a cell in GUI , and plot results if click. 
% dir should be a directorys.
%
% PlotExpts(dir, 'pattern', expr) includes only files matching regular
% expression expr

global TOPTAG;
if isempty(TOPTAG)
    TOPTAG = 'PlotExptsTOPlevel';
end

listargs = {};

j = 1;
while j <= length(varargin)
    if ischar(varargin{j})
    if strncmpi(varargin{j},'cell',3)
        listargs = {listargs{:} varargin{j} varargin{j+1}};
        j = j+1;
    elseif strncmpi(varargin{j},'pattern',3)
        listargs = {listargs{:} varargin{j} varargin{j+1}};
        j = j+1;
    elseif strncmpi(varargin{j},'Tag',3)
        TOPTAG = varargin{j+1};
        j = j+1;
    end
end
    j = j+1;
end

if(isempty(findobj('Tag',TOPTAG)))
  OTTF.tag.top = TOPTAG;
  OTTF.tag.fig = 'OTTFDataPlot';
  OTTF.tag.extra = 'OTTFExtraPlot';
  OTTF.fstrings = FindExpts(list,listargs{:});
  OTTF.listtag = 'OTTFlist';
  OTTF.tag.data = 'OTTFdata';
  OTTF.show.or = 0;
  OTTF.path = list;
  out = InitOTTF(list, OTTF, varargin{:});
  return;
else
    top = findobj('Tag',TOPTAG);
    if strmatch('store',list)
        set(top,'UserData',varargin{1});
        return;
    end
    OTTF = get(top,'UserData');
end


if strcmpi(list,'getstate')
     out = OTTF;
     it = findobj('Tag',OTTF.listtag);
     n = get(it, 'value');
     return;    
elseif strmatch(list,'close')
  exit_ottf(OTTF);
elseif strmatch(list,'setentry')
  [OTTF, n] = PlotExpts('getstate');
  OTTF = doentry(OTTF.listids(n),OTTF,'setid');
  set(OTTF.toplevel,'UserData',OTTF);
elseif strmatch(list,'nexti')
     it = findobj(gcf, 'Tag',OTTF.listtag);
     n = get(it, 'value');
     if n < length(OTTF.listids);
       n = n+1;
       set(it, 'value',n);
       doentry(OTTF.listids(n),OTTF,'setid');
     end
     figure(OTTF.toplevel);
elseif strmatch(list,'prev')
     it = findobj(gcf, 'Tag',OTTF.listtag);
     n = get(it, 'value');
     if n > 1
       n = n-1;
       set(it, 'value',n);
       doentry(OTTF.listids(n),OTTF,'setid');
     end
     figure(OTTF.toplevel);
elseif strmatch('show',list)
    type = varargin{1};
    if isfield(OTTF.show,type)
        OTTF.show.(type) = ~OTTF.show.(type);
    else
        OTTF.show.(type) = 0;
    end
    set(OTTF.toplevel,'UserData',OTTF);
elseif strmatch('store',list)
     top = findobj('Tag',TOPTAG);
     set(top,'UserData',varargin{1});
elseif strmatch(list,'update')
     top = num2str(OTTF.toplevel);
     OTTF.prefix = get(findobj('Tag','Prefix'),'String');
     OTTF.args = get(findobj('Tag','Args'),'String');
     OTTF.plot.sdfw = str2num(get(findobj('Tag','sdfw'),'String'));
     OTTF.plot.legendpos = get(findobj('Tag','LegendPos'),'value') -1;
     PlotExpts('store',OTTF);
     PlotCurrent;
elseif list(1) == '/'
    NewList(list,OTTF,listargs{:});
end



function PlotCurrent(varargin)
[OTTF,  n] = PlotExpts('getstate');
doentry(OTTF.listids(n),OTTF,'setid');

function OTTF = NewList(path, OTTF, varargin)

j =1;
findcell = 0;
while j <= length(varargin)
    if strncmpi(varargin{j},'cell',4)
    end
    j = j+1;
end
OTTF.path = path;
OTTF.fstrings = FindExpts(OTTF.path, varargin{:});
OTTF.id = 1;
lst = findobj('Tag',OTTF.listtag,'Parent',OTTF.toplevel);
set(lst,'string',OTTF.fstrings,'value',1);
OTTF.listids = 1:length(OTTF.fstrings);
set(OTTF.toplevel,'UserData',OTTF);

function fstrings = FindExpts(epath, varargin)

j =1;
findcell = 0;
rexp = [];
while j <= length(varargin)
    if strncmpi(varargin{j},'cell',4)
        j = j+1;
        findcell = varargin{j};
    elseif strncmpi(varargin{j},'pattern',4)
        j = j+1;
        rexp = varargin{j};
    end
    j = j+1;
end

if length(rexp)
    fstrings = TreeFind(epath,'name',rexp);
elseif findcell > 99
    fstrings = TreeFind(epath,'name','cell.*.mat');
elseif findcell > 0
    fstrings = TreeFind(epath,'name',['cell' num2str(findcell) '.*.mat']);
else
    fstrings = TreeFind(epath,'name','.mat');
end

function OTTF = doentry(id, OTTF,varargin)
   
j = 1;
while(j < nargin - 1)
    if strncmpi(varargin{j},'setid',5)
        OTTF.id = id;
    end
    j = j+1;
end

funcalled = 0;
argon = {'sdfw' OTTF.plot.sdfw, 'rcnmin', 10};
if get(findobj('Tag','LFP'),'value')
       argon = {argon{:} {'lfp'}};
end
   if get(findobj(OTTF.toplevel,'Tag','ShowN'),'value')
       argon = {argon{:} {'ShowN'}};
   end
   if get(findobj(OTTF.toplevel,'Tag','Sdf'),'value')
       argon = {argon{:} {'sdfall'}};
   end
   if get(findobj(OTTF.toplevel,'Tag','PlotPsych'),'value')
       argon = {argon{:} {'psych'}};
   end
   if get(findobj(OTTF.toplevel,'Tag','Pcolor'),'value')
       argon = {argon{:} {'pcolor'}};
   end
   if OTTF.plot.plotall
       argon = {argon{:} {'plotall'}};
   end
   if length(OTTF.args) > 1
       commas = strfind(OTTF.args,',');
       if isempty(commas)
           argon = {argon{:} OTTF.args};
       else
           last = 1;
           for j = 1:length(commas)
               argon = {argon{:} OTTF.args(last:commas(j)-1)};
               last = commas(j)+1;
           end
           argon = {argon{:} OTTF.args(commas(j)+1:end)};
       end
   end
   plottype = get(findobj(OTTF.toplevel,'Tag','PlotType'),'value');
   if plottype == 2 %CP hist
       argon = {argon{:} 'psych' 'CPhist'};
   elseif plottype == 3  %trials and choices
    argon = {argon{:} 'psych'  'cpt'};
   elseif plottype == 4  %trials and choices
    argon = {argon{:} 'psych' 'cp'};
   elseif plottype == 5  %cp and EM timecourse
    argon = {argon{:} 'cpt' 'emdiff' 1 'emskip' 100};
   elseif plottype == 6  %trials and choices
    ExptPsych(Expt);
    return;
   elseif plottype == 7  %trials and choices
    argon = {argon{:} 'psych' 'cpgrand'};
   elseif plottype == 8  %trials and choices
    argon = {argon{:} 'psych' 'cp' 'detrend'};
   elseif plottype == 9  %trials and choices
    argon = {argon{:} 'collapse' 2};
   elseif plottype == 10  %sequence
    argon = {argon{:} 'seq'};
   end
   
   GetFigure(OTTF.tag.fig);
   hold off; 
   colors = mycolors;
   res = PlotExpt(OTTF.fstrings{id(1)},'fbox',argon{:});
   for j = 2:length(id)
       res = PlotExpt(OTTF.fstrings{id(j)},'fbox',argon{:},'hold','forcecolor',colors{j});
   end
   OTTF.Expt = res(1).Data;
   tstr = get(get(gca,'title'),'string');
   f = fields(OTTF.show);
   for j = 1:length(f)
       if OTTF.show.(f{j})
           tstr = sprintf('%s %s=%.2f',tstr,f{j},GetEval(res.Data,f{j}));
       end
   end
   title(tstr);
   

   


function OTTF = InitOTTF(list, OTTF, varargin)

wsize = [40 35]; %in columns
j = 1;
while j < nargin -2
    if strncmpi(varargin{j},'sizestrings',7)
        OTTF.sizestrings = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'size',2)
        j = j+1;
        wsize = varargin{j};
    end
    j = j+1;
end

scrsz = get(0,'Screensize');

if scrsz(3) == 1920 & scrsz(4) == 1200
    wsc = 1.1;
elseif scrsz(3) == 1024 & scrsz(4) == 768
    wsc = 0.9;
else
    wsc = 1;
end

OTTF.wsc = wsc;
OTTF.plot.plotall = 0;
OTTF.plot.sdfw = 20;
OTTF.args = '';


cw = 10 * wsc;
bh = 9 * wsc;
listlen = 15;
wsiz = [wsize(1)*bh cw * wsize(2)];


OTTF.listids = 1:length(OTTF.fstrings);


  cntrl_box = figure('Position', [100 scrsz(4)-(wsiz(2)+bh*4) wsiz(1) wsiz(2)],...
       'NumberTitle', 'off', 'Tag',OTTF.tag.top,'Name',sprintf('OTTF %s (%d)',OTTF.path,length(OTTF.fstrings)));
  lst = uicontrol(gcf, 'Style','listbox','String',OTTF.fstrings,...
		'Callback', ' PlotExpts(''setentry'')','Tag',OTTF.listtag,...
        'Max',3,'Min',1,...
		'Position',[10 10 wsiz(1)-20 cw*listlen]);
  OTTF.toplevel = gcf;

  SPACE = 2 * wsc;
  VSPACE = 5 * wsc;

  OTTF.toplevel = gcf;
  bp(1) = SPACE; bp(2) = cw*(listlen+1); bp(3) = 40; bp(4) = 22;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' PlotExpts(''nexti'')',...
'String', '>>', 'Tag',OTTF.tag.data,'Position', bp);

  bp(1) = bp(1) + bp(3)+SPACE;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' PlotExpts(''prev'')',...
'String', '<<', 'Position', bp);

  bp(1) = bp(1) + bp(3)+SPACE;
  bp(3) = 40;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' PlotExpts(''All'')',...
'String', 'All', 'Position', bp);


  bp(1) = bp(1) + bp(3)+SPACE;
  bp(3) = 40;

  uicontrol(gcf,'Style', 'text','String','sdfw','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw*3;
  uicontrol(gcf,'Style', 'edit', 'Callback', ' PlotExpts(''update'')',...
      'String', OTTF.plot.sdfw, 'Tag', 'sdfw','Position', bp);
  

    bp(1) = SPACE; bp(3) = 25; bp(2) = bp(2) + bp(4)+VSPACE; bp(4) = 22;
bp(3) = 50;
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'LFP', 'Tag', 'LFP', 'Position', bp);

  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw*6;
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'ShowN', 'Tag', 'ShowN', 'Position', bp);

bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = 6*cw;
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'Pcolor', 'Tag', 'Pcolor', 'Position', bp);

bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = 6*cw;
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'Psych', 'Tag', 'PlotPsych', 'Position', bp);

bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = 6*cw;
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'sdf', 'Tag', 'Sdf', 'Position', bp);

  bp(1) = SPACE;
  bp(2) = bp(2)+bp(4)+VSPACE;
  bp(3) = 40;
  uicontrol(gcf,'Style', 'text','String','Plot','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw*30;
  uicontrol(gcf,'Style', 'pop',...
'String', 'Default|CP Hist|CP Timecourse| CP Trials| CP| Psych only|Grand CP|CP detrend|Collapse2|Seq', 'Tag', 'PlotType', ...
'Callback', ' PlotExpts(''update'')', 'Position', bp);

  bp(1) = SPACE;
  bp(2) = bp(2)+bp(4)+VSPACE;
  bp(3) = 40;
  uicontrol(gcf,'Style', 'text','String','Args','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw*30;
  uicontrol(gcf,'Style', 'edit', 'Callback', ' PlotExpts(''update'')',...
      'String', OTTF.args, 'Tag', 'Args','Position', bp);

  
  hm = uimenu(gcf,'Label','Mark');
  uimenu(hm,'Label','&Mark','Callback','PlotExpts(''Mark'')');
  uimenu(hm,'Label','UnMark','Callback','PlotExpts(''UnMark'')');
  uimenu(hm,'Label','Close','Callback',' PlotExpts(''close'')');
  set(gcf,'Menubar','none');
  set(gcf,'Menubar','none');

  
  OTTF.id = 1;
set(cntrl_box,'UserData',OTTF);


function CloseTag(tag)
it = findobj('Tag',tag);
if ~isempty(it)
  close(it);
end


function exit_ottf(OTTF)
  names = fieldnames(OTTF.tag);
  for j = 1:length(names)
      eval(['CloseTag(OTTF.tag.' names{j} ');']);
  end
  return;
