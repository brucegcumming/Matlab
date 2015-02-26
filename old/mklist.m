function runottf(list, varargin)

global OTTF;

OTTF.tag = 'OTTF-TopLevel';
OTTF.figtag = 'OTTFDataPlot';
OTTF.xfigtag = 'OTTFExtraPlot';
OTTF.figbtag = 'OTTFDataBPlot';
OTTF.poptag = 'OTTFPopPlot';
OTTF.seltag = 'OTTFSelections';


if(isempty(findobj('Tag',OTTF.tag)))
  if ~exist('list','var')
    list = './ottf.malist';
  end
  InitOTTF(list, varargin{:});
elseif strmatch('loadone',list)
  LoadOTTFData(OTTF.fstrings{OTTF.id});
elseif strmatch(list,'load')
  LoadAllData;
elseif strmatch(list,'All')
  OTTF.starttime = now;
  tic;
  it = findobj(gcf, 'Tag','TheList');
  for n = 1:length(OTTF.listids);
    set(it, 'value',n);
    doentry(OTTF.listids(n));
  end
  OTTF.endtime = now;
  fprintf('Done in %.2f sec\n',(OTTF.endtime - OTTF.starttime) * 24 * 60 * 60);
elseif strmatch(list,'close')
  exit_ottf(OTTF);
elseif strmatch(list,'popselect')
  setselect;
elseif strmatch(list,'next')
     it = findobj(gcf, 'Tag','TheList');
     n = get(it, 'value');
     if n < length(OTTF.listids);
       n = n+1;
       set(it, 'value',n);
       doentry(OTTF.listids(n));
     end
elseif strmatch(list,'prev')
     it = findobj(gcf, 'Tag','TheList');
     n = get(it, 'value');
     if n > 1
       n = n-1;
       set(it, 'value',n);
       doentry(OTTF.listids(n));
     end
elseif strmatch(list,'update')
     OTTF.prefix = get(findobj('Tag','Prefix'),'String');
     OTTF.plot.legendpos = get(findobj('Tag','LegendPos'),'value') -1;
     PlotCurrent(OTTF);
 elseif strmatch(list,'setplot')
  SetPlotTypes;
elseif strmatch(list,'setentry')
  it = findobj('Tag','TheList');
  n = get(it, 'value');
  doentry(OTTF.listids(n));
 elseif strmatch(list,'relist')
     listtype= get(findobj('Tag','ListType'),'value');
     if listtype == 3
        OTTF = relist(0.5,OTTF);
    elseif listtype == 1
        OTTF =relist(0,OTTF);
     elseif listtype == 2
        OTTF =relist(1,OTTF);
    end
end


function PlotCurrent(OTTF)
  it = findobj('Tag','TheList');
  n = get(it, 'value');
  doentry(OTTF.listids(n));

function SetPlotTypes()
%SetPlotTypes looks at the values of the menus for selecting
%plot types, and records them in OTTF for use elsewhere
global OTTF;

  it = findobj('Tag','plottype');
  type = get(it, 'value');
  str = get(it, 'String');
  OTTF.plot.onetype = strrep(str(type,:),' ',''); %The string shown on the menu.
 
  it = findobj('Tag','extraplot');
  type = get(it, 'value');
  str = get(it, 'String');
  OTTF.plot.extraplot = strrep(str(type,:),' ',''); %The string shown on the menu.
 
  it = findobj('Tag','tftype');
  type = get(it, 'value');
  OTTF.plot.tftype = type; 

  it = findobj('Tag','sizetype');
  type = get(it, 'value');
  OTTF.plot.sizetype = type; 
  doentry(OTTF.id);

  
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
  CloseTag(OTTF.tag);
  CloseTag(OTTF.figtag);
  CloseTag(OTTF.figbtag);
  CloseTag(OTTF.poptag);
  CloseTag(OTTF.seltag);
  CloseTag(OTTF.xfigtag);
  return;


function OTTF = relist(step,OTTF)

n = 1;
OTTF.listids = [];
for j = 1:length(OTTF.fstrings)
    load(OTTF.fstrings{j});
    sf = GetEval(Expt,'sf','mode');
    fs = GetEval(Expt,'Fs','mode');
    OTTF.sfsteps(j) = fs * sf;
    x = mod(OTTF.sfsteps(j),1.0);
    if((x < 0.9 & x > 0.1 & step == 0.5) | step == 0 | ( x < 0.1 & step ==1))
        fstrings{n} = OTTF.fstrings{j};
        OTTF.listids(n) = j;
    n = n+1;
    end
end

  lst = findobj('Tag','TheList');
  set(lst, 'String',fstrings);
  
function InitOTTF(list, varargin)
global OTTF;

OTTF.fstrings = textread(list,'%s');
OTTF.prefix = '';  
OTTF.list = list;  
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
OTTF.state.verbose = 1;
OTTF.state.includeratio = 0.5;
OTTF.state.dirw = 0.5;
OTTF.plot.onetype = 'Saccade';
OTTF.plot.tftype = 1;
OTTF.plot.sizetype = 1;
OTTF.plot.poptype = 1;
OTTF.plot.sdfsmooth = 150;
OTTF.plot.sdftype = 'exp';
OTTF.plot.legendtype = 1;
OTTF.plot.sdfstart = -5000;
OTTF.plot.sdfend = 5000;
OTTF.plot.sdfstep = 20;
OTTF.plot.legendpos = 2;
OTTF.plot.extraplot = 'None';

j = 1;
while(j < nargin)
     if strncmpi(varargin{j},'prefix',5)
     j = j+1
     OTTF.prefix = varargin{j};
     end
j = j+1;
end

idx = findstr(list,'/');
if isempty(idx)
  idx = findstr(list,'\\');
end
if isempty(idx)
  listname = list;
else
  listname = list(idx(end)+1:end);
end

  scrsz = get(0,'Screensize');
  cntrl_box = figure('Position', [100 scrsz(4)-290 270 200],...
       'NumberTitle', 'off', 'Tag',OTTF.tag,'Name',sprintf('OTTF %s',listname));
  lst = uicontrol(gcf, 'Style','listbox','String',OTTF.fstrings,...
		'Callback', 'runottf(''setentry'')','Tag','TheList',...
		'Position',[10 10 250 90]);

  SPACE = 3;
  bp(1) = SPACE; bp(2) = 130; bp(3) = 40; bp(4) = 20;
  uicontrol(gcf,'Style', 'checkbox',...
'String', 'ReFit', 'Tag', 'ReFit', 'Position', bp);

  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = 70;
  uicontrol(gcf,'style','pop','string','None|OTuning|sacdir|sacortho|Jump|SacJump|spike|TF|SZ|OP|OT|SF|SZcmp', ...
		    'Callback', 'runottf(''setplot'')', 'Tag','extraplot',...
		    'position',bp,'value',1);
        
   bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = 70;
  uicontrol(gcf,'style','pop','string','Saccade|OnSac|Onset|Compare TF|OTuning|sacdir|sacortho|Jump|SacJump|spike|TF|SZ|OP|OT|SF', ...
		    'Callback', 'runottf(''setplot'')', 'Tag','plottype',...
		    'position',bp,'value',7);
  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'style','pop','string','TF0|TFOpt|Both|Combine|PJump|OJump', ...
		    'Callback', 'runottf(''setplot'')', 'Tag','tftype',...
		    'position',bp);
  bp(1) = SPACE;
  bp(2) = bp(2) + bp(4) + SPACE;
  bp(3) = 40;
  uicontrol(gcf,'Style', 'text','String','Size','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 70;
  uicontrol(gcf,'style','pop','string','All|First|Last|Both|Combine|2|3|4|5|6', ...
		    'Callback', 'runottf(''setplot'')', 'Tag','sizetype',...
		    'position',bp);

  bp(1) = bp(1) +bp(3) + SPACE;
  uicontrol(gcf,'Style', 'text','String','Legend','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 70;
  uicontrol(gcf,'style','pop','string','Auto|UR|UL|LL|LR|Off|None', ...
		    'Callback', 'runottf(''update'')', 'Tag','LegendPos',...
		    'position',bp);
  
  bp(1) = SPACE;
  bp(2) = bp(2) + bp(4) + SPACE;
  uicontrol(gcf,'Style', 'text','String','List','Position', bp);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 70;
  uicontrol(gcf,'style','pop','string','All|FullPeriod|HalfPeriod', ...
		    'Callback', 'runottf(''relist'')', 'Tag','ListType',...
		    'position',bp);
        
  bp(1) = SPACE; bp(3) = 25; bp(2) = 110; bp(4) = 22;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'runottf(''next'')',...
'String', '>>', 'Position', bp);

  bp(1) = bp(1) + bp(3)+SPACE;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'runottf(''prev'')',...
'String', '<<', 'Position', bp);

  bp(1) = bp(1) + bp(3)+SPACE;
  bp(3) = 40;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'runottf(''All'')',...
'String', 'All', 'Position', bp);

  bp(1) = bp(1) + bp(3)+SPACE;
  bp(3) = 40;
  uicontrol(gcf,'Style', 'edit', 'Callback', 'runottf(''update'')',...
'String', OTTF.prefix, 'Tag', 'Prefix','Position', bp);
  x = 10; y = 180;

  hm = uimenu(gcf,'Label','Mark');
  uimenu(hm,'Label','Mark','Callback','listcallback(''Mark'')');
  uimenu(hm,'Label','UnMark','Callback','listcallback(''UnMark'')');
  uimenu(hm,'Label','UnMark All','Callback','listcallback(''UnMarkAll'')');
  uimenu(hm,'Label','Load One','Callback','runottf(''loadone'')');
  uimenu(hm,'Label','Load All','Callback','runottf(''load'')');
  uimenu(hm,'Label','Close','Callback','runottf(''close'')');
  hm = uimenu(gcf,'Label','Options');
  uimenu(hm,'Label','Options','Callback','runottf(''popselect'')');
  set(gcf,'Menubar','none');

fign = figure('Tag',OTTF.figtag,'Position',[300 scrsz(4)-470 512 ...
		    400]);
OTTF.id = 1;

