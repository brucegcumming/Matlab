function [out, n] = PlotFiles(list, varargin)

% runonline(list)
%
% plots online data files via matlab structures, and allows 
% these matlab structures to be saved.


global TOPTAG;
if isempty(TOPTAG)
    TOPTAG = 'plotfilesTOPlevel';
end

octc = 0;

j = 1;
while j < nargin
    if ischar(varargin{j})
    if strncmpi(varargin{j},'Tag',3)
        TOPTAG = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'function',7)
        OTTF.xfunc = varargin{j+1};
        j = j+1;
    end
end
    j = j+1;
end


if(isempty(findobj('Tag',TOPTAG)))
  if ~exist('list','var') | isempty(list)
    fprintf('Must name a file\n');
    return;
else
   % TOPTAG = list;
  end
    OTTF.tag.top = TOPTAG;
     OTTF.tag.fig = 'OTTFDataPlot';
     OTTF.tag.extra = 'OTTFExtraPlot';
     OTTF.tag.figb = 'OTTFDataBPlot';
     OTTF.tag.pop = 'OTTFPopPlot';
     OTTF.tag.select = 'OTTFSelections';
     OTTF.listtag = 'OTTFlist';
     OTTF.savefile = [list '/test.mat'];
     OTTF.listname = list;
     [path, dir]  = name2path(OTTF.listname);
     matlist = TreeFind(dir,'name','.mat');
     out = InitOTTF(matlist, OTTF, varargin{:});
else
     top = findobj('Tag',TOPTAG);
     OTTF = get(top,'UserData');
 end
 
 if strcmpi(list,'getstate')
     out = OTTF;
     it = findobj('Tag',OTTF.listtag);
     n = get(it, 'value');
     return;    
elseif strmatch('store',list)
     top = findobj('Tag',TOPTAG);
     set(top,'UserData',varargin{1});
 elseif strmatch('Save',list)
SaveMatFile(OTTF);
elseif strmatch('Refresh',list)
     OTTF.listname = get(findobj('Tag','Prefix'),'String');

     OTTF = relist(OTTF);
     PlotFiles('store',OTTF);
 elseif strmatch('loadone',list)
  LoadOTTFData(OTTF.fstrings{OTTF.id});
elseif strmatch(list,'load')
  LoadAllData;
elseif strmatch(list,'Mark')
    OTTF.data.marked(OTTF.id) = 1;
    set(findobj('Tag',TOPTAG),'UserData',OTTF);
elseif strmatch(list,'UnMark')
    OTTF.data.marked(OTTF.id) = 0;
    set(findobj('Tag',TOPTAG),'UserData',OTTF);
    RePlot(OTTF);
elseif strmatch(list,'PrintMarkList')
    idx = find(OTTF.data.marked > 0);
    fprintf('%s\n',OTTF.fstrings{idx});

elseif strmatch(list,'close')
  exit_ottf(OTTF);
elseif strmatch(list,'popselect')
  OTTF.cntrlbox = setselect(OTTF,OTTF.tag.select);
  if ~isempty(OTTF.cntrlbox) && isfield(OTTF,'xfunc')
     feval(OTTF.xfunc,'setselect');
 end
elseif strmatch(list,'TouchPoint')
    j = varargin{2};
    GetFigure(OTTF.tag.extra);
    fprintf('%s\n',OTTF.fstrings{j});
    if OTTF.plot.labels(j) > 0
        delete(OTTF.plot.labels(j));
        OTTF.plot.labels(j) = 0;
    else
        OTTF.plot.labels(j) = text(OTTF.plot.allx(j),OTTF.plot.ally(j),splitpath(OTTF.fstrings{j}));
    end
    OTTF.id = j;
    set(findobj('Tag',TOPTAG),'UserData',OTTF);
    doentry(j,OTTF);
elseif strmatch(list,'popplot')
    GetFigure(OTTF.tag.extra);
    hold off;
 for j = 1:length(OTTF.fstrings)
     if exist(OTTF.fstrings{j},'file')
         load(OTTF.fstrings{j});
         quality = CheckSpike(Expt,'noplot');
         plot(quality(1),quality(3),'o','buttondownfcn',['PlotFiles(''TouchPoint'',gcf,' num2str(j) ')']);
         OTTF.plot.allx(j) = quality(1);
         OTTF.plot.ally(j) = quality(3);
         OTTF.plot.labels(j) = 0;
         plot(OTTF.plot.allx(j),OTTF.plot.ally(j),'o','buttondownfcn',['PlotFiles(''TouchPoint'',gcf,' num2str(j) ')']);
         if(OTTF.state.verbose)
             fprintf('%s: %.2f %.2f\n',OTTF.fstrings{j},quality(1),quality(3));
         end
         hold on;
     end
 end
elseif strmatch(list,'next')
     it = findobj(gcf, 'Tag',OTTF.listtag);
     n = get(it, 'value');
     if n < length(OTTF.listids);
       n = n+1;
       set(it, 'value',n);
       doentry(OTTF.listids(n),OTTF,'setid');
     end
elseif strmatch(list,'prev')
     it = findobj(gcf, 'Tag',OTTF.listtag);
     n = get(it, 'value');
     if n > 1
       n = n-1;
       set(it, 'value',n);
       doentry(OTTF.listids(n),OTTF,'setid');
     end
 elseif strmatch(list,'replot')
     OTTF.plot.type = get(findobj('Tag','PlotType'),'value') -1;
     set(findobj('Tag',TOPTAG),'UserData',OTTF);

     RePlot(OTTF);
     
 elseif strmatch(list,'update')
     OTTF.prefix = get(findobj('Tag','Prefix'),'String');
     OTTF.plot.legendpos = get(findobj('Tag','LegendPos'),'value') -1;
     OTTF.plot.type = get(findobj('Tag','PlotType'),'value') -1;
     OTTF.state.verbose = get(findobj('Tag','Verbose'),'value');
     OTTF.plot.combine = get(findobj('Tag','Combine'),'value');
     OTTF.plot.reverse = get(findobj('Tag','Reverse'),'value');
     OTTF.plot.autoplot = get(findobj('Tag','Autoplot'),'value');
     if ~isempty(findobj('Tag','ForceZero'))
         OTTF.plot.forcezero = get(findobj('Tag','ForceZero'),'value');
         OTTF.plot.smoothing = get(findobj('Tag','Smooth'),'value');
         OTTF.plot.sdfw = str2num(get(findobj('Tag','sdfw'),'string'));
         OTTF.plot.rcsdfw = str2num(get(findobj('Tag','rcsdfw'),'string'));
         OTTF.plot.nmin = str2num(get(findobj('Tag','Nmin'),'string'));
         OTTF.savefile = get(findobj('Tag','SaveFile'),'string');
     end
     PlotFiles('store',OTTF);
     if OTTF.plot.autoplot
         PlotCurrent;
     end

elseif strmatch(list,'setplot')
     SetPlotTypes;
     if OTTF.plot.autoplot
         PlotCurrent;
     end
elseif strmatch(list,'setentry')
    if OTTF.plot.autoplot
        PlotCurrent;
    end
elseif strmatch(list,'Plot')
    PlotCurrent;
elseif strmatch(list,'relist')
     listtype= get(findobj('Tag','ListType'),'value');
     relist(listtype);
end


function OTTF = RePlot(OTTF)

argon = {};
if get(findobj('Tag','LFP'),'value')
    argon = {argon{:} {'lfp'}};
end
if get(findobj('Tag','ShowN'),'value')
    argon = {argon{:} {'ShowN'}};
end

if get(findobj('Tag','Pcolor'),'value')
    argon = {argon{:} {'pcolor'}};
end

if strncmpi('rates',OTTF.plot.onetype,5)
    argon = {argon{:} {'rates'}};
end
if strncmp('cycSDF',OTTF.plot.onetype,6)
    argon = {argon{:} {'sdf'} {'periodic'} {'sdfw'} OTTF.plot.sdfw};
end
if strncmp('SDF',OTTF.plot.onetype,3)
    argon = {argon{:} {'sdf'} {'showsum'} {'sdfw'} OTTF.plot.sdfw};
end

GetFigure(OTTF.tag.fig);
Expt = OTTF.ExptData;
if ~isempty(OTTF.Txtdata)
    if OTTF.plot.type == 1
        NsineRC(OTTF.ExptData,OTTF.Txtdata,'monoc');
    else
        NsineRC(Expt,OTTF.Txtdata,[],{});
    end
elseif(Expt.isrc)
%             if strmatch(Expt.Stimvals.et,'dO')
%                 PlotRevCor(Expt,'sdfw',100,argon{:});
 %                PlotRevCorAny(Expt,'sdfw',OTTF.plot.rcsdfw,argon{:});
 %           else
 if OTTF.plot.type == 2
     PlotRevCorAny(Expt,'sdfw',OTTF.plot.rcsdfw,'nmin',OTTF.plot.nmin,argon{:},'yval',0,'psych');
 elseif OTTF.plot.type == 3
     PlotRevCorAny(Expt,'sdfw',OTTF.plot.rcsdfw,'nmin',OTTF.plot.nmin,argon{:},'yval',0);
 else
     PlotRevCorAny(Expt,'sdfw',OTTF.plot.rcsdfw,'nmin',OTTF.plot.nmin,argon{:});
 end
 %end
elseif(OTTF.plot.reverse)
    OTTF.data.plotdata{OTTF.id} = PlotRates(Expt, Expt.Stimvals.e2, 'type2',Expt.Stimvals.et, 'legendpos',OTTF.plot.legendpos,argon{:});
else
    OTTF.data.plotdata{OTTF.id} = PlotRates(Expt, Expt.Stimvals.et, 'type2',Expt.Stimvals.e2, 'legendpos',OTTF.plot.legendpos,argon{:});
end

function doentry(id, OTTF,varargin)
   
j = 1;
while(j < nargin - 1)
    if strncmpi(varargin{j},'setid',5)
        OTTF.id = 1;
    end
    j = j+1;
end

argon = {};
if get(findobj('Tag','LFP'),'value')
    argon = {argon{:} {'lfp'}};
end
if get(findobj('Tag','ShowN'),'value')
    argon = {argon{:} {'ShowN'}};
end

if get(findobj('Tag','Pcolor'),'value')
    argon = {argon{:} {'pcolor'}};
end

if strncmpi('rates',OTTF.plot.onetype,5)
    argon = {argon{:} {'rates'}};
end
if strncmp('cycSDF',OTTF.plot.onetype,6)
    argon = {argon{:} {'sdf'} {'periodic'} {'sdfw'} OTTF.plot.sdfw};
end
if strncmp('SDF',OTTF.plot.onetype,3)
    argon = {argon{:} {'sdf'} {'showsum'} {'sdfw'} OTTF.plot.sdfw};
end

GetFigure(OTTF.tag.fig);
Txt = [];
np = 1;
[nr,nc] = Nsubplots(length(id));
for j = 1:length(id)
    subplot(nr,nc,j);
    PlotExpt(OTTF.fstrings{id(j)},argon{:});
end
PlotFiles('store',OTTF);

 function SaveMatFile(OTTF)
 
  Expt = OTTF.ExptData;
  save(OTTF.savefile,'Expt');

 function args = CheckExpt(Expt,args)

 if Expt.Stimvals.st == 2 || Expt.Stimvals.st == 15  %% RDS or RLS
     args = {args{:} 'Uncorr'};
 end
 if strcmp(Expt.Stimvals.et,'sf') | strcmp(Expt.Stimvals.et,'tf')
     args = {args{:} 'LogX'};
 end
 if strmatch(Expt.Stimvals.et,{'dp', 'dO', 'dx'}) & ~ Expt.isrc
     args = {args{:} 'AddMonoc'};
 end

function PlotCurrent()
    [OTTF, n] = PlotFiles('getstate');

    doentry(OTTF.listids(n),OTTF,'setid');

function SetPlotTypes()
%SetPlotTypes looks at the values of the menus for selecting
%plot types, and records them in OTTF for use elsewhere

OTTF = PlotFiles('getstate');

it = findobj('Tag','plottype');
type = get(it, 'value');
str = get(it, 'String');
OTTF.plot.onetype = strrep(str(type,:),' ',''); %The string shown on the menu.

it = findobj('Tag','extraplot');
if ~isempty(it)
    type = get(it, 'value');
    str = get(it, 'String');
    OTTF.plot.extraplot = strrep(str(type,:),' ',''); %The string shown on the menu.
end
doentry(OTTF.id,OTTF);
PlotFiles('store',OTTF);


function LoadAllData()
    global OTTF;

    for j = 1:length(OTTF.fstrings)
        fprintf('%s',OTTF.fstrings{j});
        LoadOTTFData(OTTF.fstrings{j});
    end

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


function OTTF = relist(OTTF)


    [path, dir]  = name2path(OTTF.listname);
    matlist = TreeFind(dir,'name','.mat');
    OTTF.strings = matlist;
    OTTF.fstrings = matlist;
    OTTF.listids = 1:length(OTTF.fstrings);
    lst = findobj('Tag',OTTF.listtag)
    set(lst, 'String',OTTF.strings,'value',1);

function OTTF = InitOTTF(mfiles, OTTF, varargin)
global bgcfileprefix;

OTTF.fstrings = mfiles;
OTTF.strings = mfiles;

if ~isempty(bgcfileprefix)
    OTTF.prefix = bgcfileprefix;
else
    OTTF.prefix = '';
end

OTTF.listids = 1:length(OTTF.fstrings);
OTTF.state.showonset= 0;
OTTF.state.minsacs = 1;
OTTF.state.minspikes = 1;
OTTF.state.mintrials = 1;
OTTF.state.minsacsize = 0.1;
OTTF.state.maxsacsize = 1.0;
OTTF.state.addedspikes = 0;
OTTF.state.invertspikes = 0;
OTTF.state.rebuild = 0;
OTTF.state.verbose = 0;
OTTF.state.includeratio = 0.5;
OTTF.state.dirw = 0.5;
OTTF.plot.onetype = 'Saccade';
OTTF.plot.type = 0;
OTTF.plot.nmin = 100;
OTTF.plot.tftype = 1;
OTTF.plot.sizetype = 1;
OTTF.plot.poptype = 1;
OTTF.plot.sdfw = 150;
OTTF.plot.rcsdfw = 100;
OTTF.plot.sdftype = 'exp';
OTTF.plot.legendtype = 1;
OTTF.plot.sdfstart = -5000;
OTTF.plot.sdfend = 5000;
OTTF.plot.sdfstep = 20;
OTTF.plot.legendpos = 2;
OTTF.plot.extraplot = 'None';
OTTF.plot.forcezero = 0;
OTTF.plot.autoplot = 0;
OTTF.plot.smoothing = 0;
OTTF.plot.combine = 1;
OTTF.plot.reverse = 0;

j = 1;
while(j < nargin-1)
    if strncmpi(varargin{j},'prefix',5)
        j = j+1
        OTTF.prefix = varargin{j};
    elseif strncmpi(varargin{j},'verbose',5)
        j = j+1
        OTTF.state.verbose = 1;
    end
    j = j+1;
end



for j = 1:length(OTTF.fstrings)
    OTTF.data.marked(j) = 0;
end


  scrsz = get(0,'Screensize');
  boxw = 300;
  cntrl_box = figure('Position', [100 scrsz(4)-320 boxw 280],...
       'NumberTitle', 'off', 'Tag',OTTF.tag.top,'Name',sprintf('Online Data %s (%d)',OTTF.listname,length(OTTF.fstrings)));
  lst = uicontrol(gcf, 'Style','listbox','String',OTTF.strings,...
      'Max',10,...
		'Callback', ' PlotFiles(''setentry'')','Tag',OTTF.listtag,...
		'Position',[10 10 boxw-20 150]);

  SPACE = 3;
  VSPACE = 10;
  cw = 10;

  bp(1) = SPACE; bp(2) = 160; bp(3) = 40; bp(4) = 22;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' PlotFiles(''next'')',...
'String', '>>', 'Position', bp);

  bp(1) = bp(1) + bp(3)+SPACE;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' PlotFiles(''prev'')',...
'String', '<<', 'Position', bp);



bp(1) = bp(1) + bp(3)+SPACE;
bp(3) = 40;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' PlotFiles(''Refresh'')',...
'String', 'Refresh', 'Position', bp);
bp(1) = bp(1) + bp(3)+SPACE;
bp(3) = 40;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' PlotFiles(''Plot'')',...
'String', 'Plot', 'Position', bp);
bp(1) = bp(1) + bp(3)+SPACE;
bp(3) = 40;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' PlotFiles(''Save'')',...
'String', 'Save', 'Position', bp);
bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = 70;
uicontrol(gcf,'style','pop','string','Counts|Rates|SDF|cycSDF', ...
    'Callback', ' PlotFiles(''setplot'')', 'Tag','plottype',...
    'position',bp,'value',1);


bp(1) = SPACE;
bp(2) = bp(2)+bp(4)+VSPACE;
bp(3) = boxw-SPACE;
uicontrol(gcf,'Style', 'edit', 'Callback', 'PlotFiles(''Refresh'')',...
'String', OTTF.listname, 'Tag', 'Prefix','Position', bp);


    bp(1) = SPACE; bp(3) = 25; bp(2) = bp(2) + bp(4)+VSPACE; bp(4) = 22;
bp(3) = 50;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' PlotFiles(''update'')',...
'String', 'Auto', 'Tag', 'Autoplot', 'Position', bp,'value',OTTF.plot.autoplot);


  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' PlotFiles(''update'')',...
'String', 'ShowN', 'Tag', 'ShowN', 'Position', bp);

bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = 6*cw;
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'Pcolor', 'Tag', 'Pcolor', 'Position', bp);

bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' PlotFiles(''update'')',...
'String', 'Verbose', 'Tag', 'Verbose', 'Position', bp);

bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' PlotFiles(''update'')',...
'String', 'Combine', 'Tag', 'Combine', 'Position', bp);
bp(1) = SPACE;
bp(2) = bp(2)+bp(4)+SPACE;
  uicontrol(gcf,'Style', 'checkbox','Callback', ' PlotFiles(''update'')',...
'String', 'Reverse', 'Tag', 'Reverse', 'Position', bp);

  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw * 6;
  uicontrol(gcf,'Style', 'text','String','Legend','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw * 5;
  uicontrol(gcf,'style','pop','string','Auto|UR|UL|LL|LR|Off|None', ...
		    'Callback', ' PlotFiles(''update'')', 'Tag','LegendPos',...
		    'position',bp);
  
  bp(1) = bp(1) +bp(3) + SPACE;
  uicontrol(gcf,'Style', 'text','String','Type','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw *7;
  uicontrol(gcf,'style','pop','string','Default|Monocs|PsyRC|RC0|New', ...
		    'Callback', ' PlotFiles(''replot'')', 'Tag','PlotType',...
		    'position',bp);

        x = 10; y = 180;

  hm = uimenu(gcf,'Label','Mark');
  uimenu(hm,'Label','&Mark','Callback','PlotFiles(''Mark'')');
  uimenu(hm,'Label','UnMark','Callback','PlotFiles(''UnMark'')');
  uimenu(hm,'Label','List','Callback','PlotFiles(''PrintMark'')');
  uimenu(hm,'Label','UnMark All','Callback','PlotFiles(''UnMarkAll'')');
  uimenu(hm,'Label','Load One','Callback',' PlotFiles(''loadone'')');
  uimenu(hm,'Label','Load All','Callback',' PlotFiles(''load'')');
  uimenu(hm,'Label','Close','Callback',' PlotFiles(''close'')');
  hm = uimenu(gcf,'Label','Options');
  uimenu(hm,'Label','Options','Callback',' PlotFiles(''popselect'')');
  set(gcf,'Menubar','none');
  hm = uimenu(gcf,'Label','&Next','Callback','PlotFiles(''next'');');

OTTF.figid.mainplot = figure('Tag',OTTF.tag.fig,'Position',[300 scrsz(4)-470 512 ...
		    400]);
OTTF.id = 1;
relist(OTTF);
set(cntrl_box,'UserData',OTTF);

function data = Condense(res)


data.sf = res.x;
data.n = res.n;
data.eye = res.y;
data.count = res.means;
data.sd = res.sd;
data.name = res.name;
data.extras = res.extras;
data.width = GetEVal(res.Data,'wi');
data.height = GetEVal(res.Data,'hi');
data.ori = GetEVal(res.Data,'or');
data.x = GetEVal(res.Data,'xo');
data.y = GetEVal(res.Data,'yo');


fz = GetEval(res.Data,'fz');
if isnan(fz)
    fz = 72.24;
end
dur = mean([res.Data.Trials.End] - [res.Data.Trials.Start]);
data.duration = round(dur * fz/10000) * 1000/fz;


function OTTF = GetSimple(OTTF, varargin)


[a, b, names] = textread('../tf/simplelist','%f %f %s');

for j = 1:length(OTTF.fstrings)
    name = splitpath(OTTF.fstrings{j});
    OTTF.data.sxcx(j) = NaN;
    for k = 1:length(names)
        if strncmp(names{k},name,11)
            OTTF.data.sxcx(j) = b(k);
        end
    end
end
PlotFiles('store',OTTF);


function cntrl_box = setselect(OTTF, tag)

cntrl_box = [];

if isempty(findobj('Tag',tag))
    
    if isempty(OTTF.plot.forcezero)
        OTTF.plot.forcezero = 0;
    end
    if isempty(OTTF.plot.autoplot)
        OTTF.plot.autoplot = 0;
    end
    if isempty(OTTF.plot.smoothing)
        OTTF.plot.smoothing = 0;
    end

  SPACE = 5;
  VSPACE = 2;
  h = 220;
  scrsz = get(0,'Screensize');
 wsc = scrsz(3) /1000;
  cntrl_box = figure('Position', [200 scrsz(4)-(h+30)*wsc 280*wsc h*wsc], 'Menubar', 'none',...
       'NumberTitle', 'off', 'Tag',tag,'Name','Section Criteria');

   cw = 10*wsc;
   ch = 18;
   bh = 18*wsc;
   bp(1) = SPACE;
   bp(2) = ch+VSPACE;
   bp(3) = cw*7;
   bp(4) = ch;
  
  uicontrol(gcf,'Style', 'CheckBox','String','Force Zero','Position', bp,...
      'Tag','ForceZero','Callback','PlotFiles(''update'')','value',OTTF.plot.forcezero);
  bp(1) = bp(1) + bp(3) + SPACE;
  bp(3) = cw*5;
    uicontrol(gcf,'Style', 'CheckBox','String','Smooth','Position', bp,...
      'Tag','Smooth','Callback','PlotFiles(''update'')','value',OTTF.plot.smoothing);

  bp(1) = SPACE;
  bp(2) = bp(2) + ch + VSPACE;
  bp(3) = cw*4;
  uicontrol(gcf,'Style', 'text','String','sdfw','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  uicontrol(gcf,'Style', 'edit', 'Callback', ' PlotFiles(''update'')',...
    'String', sprintf('%d',OTTF.plot.sdfw), 'Tag', 'sdfw','Position', bp);

  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw*6;
  uicontrol(gcf,'Style', 'text','String','RCsdfw','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 3 * cw;
  uicontrol(gcf,'Style', 'edit', 'Callback', ' PlotFiles(''update'')',...
    'String', sprintf('%d',OTTF.plot.rcsdfw), 'Tag', 'rcsdfw','Position', bp);

  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = cw*6;
  uicontrol(gcf,'Style', 'text','String','Nmin','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 3 * cw;
  uicontrol(gcf,'Style', 'edit', 'Callback', ' PlotFiles(''update'')',...
    'String', sprintf('%d',OTTF.plot.nmin), 'Tag', 'Nmin','Position', bp);

bp(1) = SPACE;
  bp(2) = bp(2) + ch + VSPACE;
  bp(3) = cw*6;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ' PlotFiles(''Save'')',...
'String', 'Save As', 'Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 20 * cw;
  uicontrol(gcf,'Style', 'edit', 'Callback', ' PlotFiles(''update'')',...
    'String', OTTF.savefile, 'Tag', 'SaveFile','Position', bp);

OTTF.wsc = wsc;
  OTTF.cntrlbox = cntrl_box;
  OTTF.selectlast = bp(2);
  PlotFiles('store',OTTF);
end
