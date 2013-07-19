function Result = psfgui(varargin)

global PSYCH;
PSYCH.tag = 'PSYCH-TopLevel';

if(isempty(findobj('Tag',PSYCH.tag)))

    clear global PSYCH;
    global PSYCH;
    PSYCH.tag = 'PSYCH-TopLevel';
    PSYCH.figtag = 'PSYCH-PsfPlot';
    PSYCH.figbtag = 'PSYCH-DatPlot';
PSYCH.seltag = 'PSYCH-Select';
PSYCH.plotsds = 0;
PSYCH.plotmeans = 0;
PSYCH.colors = mycolors;
PSYCH.Result = [];
PSYCH.doblocks = [];
PSYCH.opt.nmin = 10;
PSYCH.opt.xmin = -100;
PSYCH.opt.xmax = 100;
PSYCH.opt.showinit = 0;
PSYCH.opt.meanlabel = 0;
PSYCH.opt.counterrs = 0;
  PSYCH.opt.conditions = {};
PSYCH.opt.vars{1} = '';
PSYCH.opt.vars{2} = '';
PSYCH.opt.vars{3} = '';
PSYCH.opt.sdlabel = 0;
PSYCH.opt.first = 0;
PSYCH.opt.last = 0;
PSYCH.opt.forceread = 1;
PSYCH.opt.collapsedata = 0;
PSYCH.opt.holdon = 0;
PSYCH.opt.blockexpts = 0;
PSYCH.opt.blocksize = 0;
PSYCH.opt.pcorrect = 0;
PSYCH.opt.prefix = '';
  j = 1;
  nfiles = 1;
while j <= nargin
  if(strncmpi(varargin{j},'sortby',4))
    j = j+1;
  elseif(strncmpi(varargin{j},'prefix',4))
    j = j+1;
    PSYCH.opt.prefix = varargin{j};
  else
    PSYCH.files{nfiles} = varargin{j};
    PSYCH.data{nfiles} = [];
    PSYCH.listids(nfiles) = nfiles;
    nfiles = nfiles+1;
  end
  j = j+1;
end
  InitPSYCH(PSYCH);
  return;
else
  flag = varargin{1};
end

if strmatch(flag,'close')
  exit_psych(PSYCH);
  Result = PSYCH.Result;
elseif strmatch(flag,'next')
     it = findobj(gcf, 'Tag','TheList');
  CloseTag('Blocks');
     PSYCH.doblocks = [];
     n = get(it, 'value');
     if n < length(PSYCH.listids);
       n = n+1;
       set(it, 'value',n);
       PSYCH = doentry(PSYCH, PSYCH.listids(n));
     end
elseif strmatch(flag,'prev')
  CloseTag('Blocks');
     PSYCH.doblocks = [];
     it = findobj(gcf, 'Tag','TheList');
     n = get(it, 'value');
     if n > 1
       n = n-1;
       set(it, 'value',n);
       PSYCH = doentry(PSYCH, PSYCH.listids(n));
     end
elseif strmatch(flag,'closeblocks')
  CloseTag('Blocks');
  PSYCH.doblocks = [];
elseif strmatch(flag,'clearblocks')
  it = 1;
  j = 1;
  while ~isempty(it)
    it = findobj(gcf,'Tag',sprintf('Block%d',j));
    if(~isempty(it))
      set(it,'value',0);
    end
    j = j+1;
  end
elseif strmatch(flag,'checkblocks')
  it = 1;
  j = 1;
  while ~isempty(it)
    it = findobj(gcf,'Tag',sprintf('Block%d',j));
    if(~isempty(it))
      PSYCH.doblocks(j) = get(it,'value');
    end
    j = j+1;
  end
  it = findobj('Tag','TheList');
  n = get(it, 'value');
  PSYCH = doentry(PSYCH, PSYCH.listids(n),'doblocks',PSYCH.doblocks);
elseif strmatch(flag,'scroll')
     it = findobj('Tag','blockscroll');
     j = get(it,'value');
     scrollblocks(j);     
elseif strmatch(flag,'setentry')
  it = findobj('Tag','TheList');
  n = get(it, 'value');
  PSYCH = doentry(PSYCH, PSYCH.listids(n));
  CloseTag('Blocks');
  PSYCH.blocks = [];
  set(findobj('Tag','Addfile'),'string',PSYCH.files{PSYCH.listids(n)});
elseif(strncmpi(flag,'showoptions',5))
  showoptions(PSYCH);
elseif(strncmpi(flag,'showblocks',5))
  showblocks(PSYCH.blocklist, PSYCH);
elseif strmatch(flag,'touch')
    fprintf('%s: SD %.3g\n',PSYCH.labels{varargin{2}},PSYCH.fits(varargin{2}).sd);
elseif(strncmpi(flag,'Addfile',5))
  filename = get(findobj('Tag','Addfile'),'string');
  if isempty(strmatch(filename,PSYCH.files))
    nfiles = length(PSYCH.files)+1;
    PSYCH.files{nfiles} = filename;
    PSYCH.data{nfiles} = [];
    PSYCH.listids(nfiles) = nfiles;
    set(findobj('Tag','TheList'),'String',PSYCH.files);
  end
elseif(strncmpi(flag,'Update',4))
  PSYCH.opt.first = str2num(get(findobj('Tag','FirstTrial'),'string'));
  PSYCH.opt.last = str2num(get(findobj('Tag','LastTrial'),'string'));
  PSYCH.opt.blocksize = str2num(get(findobj('Tag','BlockSize'),'string'));
  PSYCH.opt.nmin = str2num(get(findobj('Tag','Mintrials'),'string'));
 
  nc = 0;
  x = get(findobj('Tag','Vara'),'string');
  PSYCH.opt.vars{1} = x;
  PSYCH.opt.conditions = {};
  if ~isempty(x)
      nc = nc+1;
      PSYCH.opt.conditions{nc} = x;
  end

  x = get(findobj('Tag','Varb'),'string');
  PSYCH.opt.vars{2} = x;
  if ~isempty(x)
      nc = nc+1;
      PSYCH.opt.conditions{nc} = x;
  end

  x = get(findobj('Tag','Varc'),'string');
  PSYCH.opt.vars{3} = x;
  if ~isempty(x)
      nc = nc+1;
      PSYCH.opt.conditions{nc} = x;
  end
  
  
  
  x = get(findobj('Tag','ShowInitial'),'value');
  if ~isempty(x)
    PSYCH.opt.showinit = x;
    PSYCH.opt.meanlabel = get(findobj('Tag','ShowMean'),'value');
    PSYCH.opt.sdlabel = get(findobj('Tag','ShowSD'),'value');
    PSYCH.opt.forceread = get(findobj('Tag','ForceRead'),'value');
    PSYCH.opt.blockexpts = get(findobj('Tag','BlockExpts'),'value');
    PSYCH.opt.collapsedata = get(findobj('Tag','CollapseData'),'value');
    PSYCH.opt.holdon = get(findobj('Tag','HoldOn'),'value');
    PSYCH.opt.prefix = get(findobj('Tag','Prefix'),'string');
    PSYCH.opt.counterrs = get(findobj('Tag','BadRate'),'value');
    PSYCH.plotsds = get(findobj('Tag','PlotSD'),'value');
    PSYCH.plotmeans = get(findobj('Tag','PlotMean'),'value');
    PSYCH.opt.pcorrect = get(findobj('Tag','PercentCorrect'),'value');
    PSYCH.explist = sscanf(get(findobj('Tag','Explist'),'string'),'%d,');
  end
  PSYCH.opt.xmin = str2num(get(findobj('Tag','Xmin'),'string'));
  PSYCH.opt.xmax = str2num(get(findobj('Tag','Xmax'),'string'));
  psfgui('setentry');
  figure(findobj('Tag',PSYCH.tag));
end



function PSYCH = doentry(PSYCH, id,varargin)

j = 1;
while j <= nargin-2;
  if(strncmpi(varargin{j},'doblocks',4))
    j = j+1;
    PSYCH.doblocks = varargin{j};
  end
  j = j+1;
end

Trials = [];
file = PSYCH.files{id};
if ~isempty(PSYCH.opt.prefix)
  file = [PSYCH.opt.prefix PSYCH.files{id}];
end



if isempty(PSYCH.data{id})  | PSYCH.opt.forceread
[Data, Trials, PSYCH.blocklist] = readpsychfile(PSYCH,file);
  PSYCH.data{id} = Data;
  if ~isempty(Trials)  %% This is a monkey serialout file
      PSYCH.trials{id} = Trials;
      PSYCH.explists{id} = PSYCH.blocklist;
  else
      PSYCH.trials{id}= [];
  end
else
  Data = PSYCH.data{id};
  if ~isempty(PSYCH.trials{id})
      Trials = PSYCH.trials{id};
  end
end

if ~isempty(Trials)
    if(PSYCH.opt.last > PSYCH.opt.first)
        lasttrial = PSYCH.opt.last;
    else
        lasttrial = length(Trials);
    end
  if PSYCH.opt.blockexpts
      Data = CountPsychTrials(Trials,'expts',PSYCH.blocklist,'sortby','sd','nmin', ...
          PSYCH.opt.nmin);
  elseif PSYCH.opt.blocksize > 10
      Data = CountPsychTrials(Trials,'expts',PSYCH.opt.first:PSYCH.opt.blocksize:lasttrial,'sortby',PSYCH.opt.vars{1},'nmin', ...
          PSYCH.opt.nmin);
  else
      Data = CountPsychTrials(Trials,'skip',PSYCH.opt.first,'last',PSYCH.opt.last,'sortby',PSYCH.opt.vars{1},'nmin', ...
          PSYCH.opt.nmin);
  end
end


if ~isempty(Data)
  idx = find([Data.n] >= PSYCH.opt.nmin);
  Data = Data(idx);
  it = findobj('Tag',PSYCH.figtag);
  figure(it);

  nfit = 1;
if ~PSYCH.opt.holdon
  hold off;
end

PSYCH.Result.Data = Data;
nresample = str2num(get(findobj('Tag','ResampleN'),'String'));

%explist = [17 18 19 20];
explist = [];
if isfield(PSYCH,'explist')
    explist = PSYCH.explist;
end
fits = [];
ex = unique([Data.expno]);
if PSYCH.opt.pcorrect %% convert choices in two directions to percent correct
    id = find([Data.expno] == ex(2));
    for j = id
        Data(j).x = Data(j).x;
        Data(j).expno = ex(1);
        Data(j).resp = Data(j).n - Data(j).resp;
        Data(j).p = Data(j).resp./Data(j).n;
    end
end
for ex = unique([Data.expno]);
  if isempty(explist) | ismember(ex,explist)
      if sum([Data.x] < 0) == 0
          nd = length(Data)+1;
          Data(nd).x = 0;
          Data(nd).n = 100000;
          Data(nd).p = 0.5;
          Data(nd).resp = Data(nd).n/2;
          Data(nd).expno = ex;
          fit = fitpsf(Data,'expno',ex,'xmin',PSYCH.opt.xmin,'xmax',PSYCH.opt.xmax);
      else
          fit = fitpsf(Data,'expno',ex,'xmin',PSYCH.opt.xmin,'xmax',PSYCH.opt.xmax);
      end
      if nresample
          fit = resamplepsf(fit,nresample,'xmin',PSYCH.opt.xmin,'xmax',PSYCH.opt.xmax);
      end
      PSYCH.Result.fits{nfit} = fit;
      
      if ~isnan(fit.fit(1))
          idx = find([Data.expno] == ex);
          if isfield(Data,'sname')
              labels{nfit} = Data(idx(1)).sname;
          else
              labels{nfit} = Data(idx(1)).name;
          end

          desc = [];
          for k = 1:length(PSYCH.opt.conditions)
              if isfield(Data,PSYCH.opt.conditions{k}) 
                  val = eval(['mean([Data(idx).' PSYCH.opt.conditions{k} ']);']);
              desc = [desc sprintf(' %s=%.2f',PSYCH.opt.conditions{k},val)];
              end
          end
          if ~isempty(desc)
              labels{nfit} = sprintf('%s %s',labels{nfit},desc);
          end
          if(PSYCH.opt.meanlabel)
              labels{nfit} = sprintf('%s PSE %.3g',labels{nfit},fit.fit(1));
          end
          if(PSYCH.opt.sdlabel)
              labels{nfit} = sprintf('%s %.3g',labels{nfit},abs(fit.fit(2)));
          end
          
          if nfit > length(PSYCH.colors)
              colorid = length(PSYCH.colors);
          else
              colorid = nfit;
          end
          
          h(nfit) = plotpsych(fit.data,fit.fit(1),fit.fit(2),'color', ...
              PSYCH.colors{colorid},'shown');
          if(PSYCH.opt.showinit)
              plotpsych(fit.data,fit.initial(1),fit.initial(2),'color', ...
                  PSYCH.colors{colorid});
          end
    fits(nfit).pse = fit.fit(1);
    fits(nfit).sd = fit.fit(2);
    fits(nfit).x = nfit;
    
    sb = PSYCH.opt.vars{1};
    if ~isempty(sb)
        fits(nfit).(sb) = mean([Data(idx).(sb)]);
    end
      end
    nfit = nfit+1;
      end
  end
end
legendpos = get(findobj('Tag','LegendPos'),'value')-1;
if (legendpos < 6)
  if nfit > 1
      PSYCH.fits = fits;
      PSYCH.labels = labels;
        legend(h,labels,legendpos);
    end 
    if nfit > 2 
        if ~isempty(sb)
            xv = [fits.(sb)];
            xlb = PSYCH.opt.vars{1};
        else
            xv = [fits.x];
            xlb = 'Expt ID';
        end
      if PSYCH.plotsds
        GetFigure(PSYCH.figbtag);
        hold off;
        plot(xv,[fits.sd]);
        hold on;
        ylabel('Threshold (SD)');
        xlabel(xlb);
        for j = 1:length(xv)
            plot(xv(j),fits(j).sd,'o','MarkerFaceColor','b','buttondownfcn',['psfgui(''touch'',' num2str(j) ');']);
        end
      elseif PSYCH.plotmeans
        GetFigure(PSYCH.figbtag);
        hold off;
        plot(xv,[fits.pse]);
      end
      
    end
end



function CloseTag(tag)
it = findobj('Tag',tag);
if ~isempty(it)
  close(it);
end


function exit_psych(PSYCH)
  CloseTag(PSYCH.tag);
  CloseTag(PSYCH.figtag);
  CloseTag(PSYCH.seltag);
  CloseTag(PSYCH.figbtag);
  CloseTag('Blocks');
  return;


function InitPSYCH(PSYCH)

scrsz = get(0,'Screensize');
  wsc = scrsz(3)/1000;  
  cntrl_box = figure('Position', [100 scrsz(4)-270*wsc 300*wsc 250*wsc],...
  'NumberTitle', 'off', 'Tag',PSYCH.tag,'Name','PSF Fitter');
  lst = uicontrol(gcf, 'Style','listbox','String',PSYCH.files,...
		'Callback', 'psfgui(''setentry'')','Tag','TheList',...
		'Position',[5 20 190 90]*wsc);


  SPACE = 3;
  VSPACE = 3;
  cw = 5;
  bp(1) = SPACE; bp(2) = 200; bp(3) = 40; bp(4) = 20;
  bp(1) = SPACE;
  bp(2) = bp(2) + bp(4) + SPACE;
  bp(3) = 30;
  uicontrol(gcf,'Style', 'text','String','Min N','Position', bp * wsc);
  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'edit','String','10','Position', bp * wsc,'Tag','Mintrials','Callback', ...
	    'psfgui(''Update'')','Backgroundcolor',[1 1 1]);

  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'text','String','Type','Position', bp * wsc);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 70;
  uicontrol(gcf,'style','pop','string','All|Split', ...
		    'Callback', 'psfgui(''setplot'')', 'Tag','sizetype',...
		    'position',bp);


  uicontrol(gcf,'Style', 'text','String','Resample N','Position', bp * wsc);
  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'edit','String','0','Position', bp * wsc,'Tag','ResampleN','Callback', ...
	    'psfgui(''Update'')','Backgroundcolor',[1 1 1]);

  bp(3) = 40;
  bp(2) = bp(2)-bp(4)-VSPACE;
  bp(1) = SPACE;
  uicontrol(gcf,'Style', 'text','String','First','Position', bp * wsc);
  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'edit','String','1','Position', bp * wsc,'Tag','FirstTrial','Callback', ...
	    'psfgui(''Update'')','Backgroundcolor',[1 1 1]);
  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'text','String','Last','Position', bp * wsc);
  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'edit','String','1','Position', bp * wsc,'Tag','LastTrial','Callback', ...
	    'psfgui(''Update'')','Backgroundcolor',[1 1 1]);

  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'text','String','Block','Position', bp * wsc);
  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'edit','String','0','Position', bp * wsc,'Tag','BlockSize','Callback', ...
	    'psfgui(''Update'')','Backgroundcolor',[1 1 1]);


    bp(2) = bp(2)-bp(4)-VSPACE;
  bp(1) = SPACE;
  uicontrol(gcf,'Style', 'text','String','Xmin','Position', bp * wsc);
  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'edit','String',num2str(PSYCH.opt.xmin),'Position', bp * wsc,'Tag','Xmin','Callback', ...
	    'psfgui(''Update'')','Backgroundcolor',[1 1 1]);
  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'text','String','Xmax','Position', bp * wsc);
  bp(1) = bp(1)+bp(3)+SPACE;
  uicontrol(gcf,'Style', 'edit','String','100','Position', bp * wsc,'Tag','Xmax','Callback', ...
	    'psfgui(''Update'')','Backgroundcolor',[1 1 1]);


  bp(1) = SPACE; bp(3) = 25; bp(2) = 110; bp(4) = 22;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'psfgui(''next'')',...
'String', '>>', 'Position', bp * wsc);

  bp(1) = bp(1) + bp(3)+SPACE;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'psfgui(''prev'')',...
'String', '<<', 'Position', bp * wsc);

  bp(1) = bp(1) + bp(3)+SPACE;
  bp(3) = 40;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'psfgui(''All'')',...
'String', 'All', 'Position', bp * wsc);
  bp(1) = bp(1) +bp(3) + SPACE;
  bp(3) = 70;
        
        uicontrol(gcf,'style','pop','string','Auto|UR|UL|LL|LR|Off fig|None', ...
		    'Callback', 'psfgui(''update'')', 'Tag','LegendPos',...
		    'position',bp*wsc);
  x = 10; y = 180;
  bp(1) = SPACE;
  bp(2) = bp(2) + bp(4)+VSPACE;
  bp(3) =30;
  uicontrol(gcf,'Style', 'text','String','File','Position', bp * wsc);
  bp(1) = bp(1) + bp(3)+SPACE;
  bp(3) = 250;
  uicontrol(gcf,'Style', 'edit','String','','Position', bp * wsc,'Tag','Addfile','Callback', ...
	    'psfgui(''Addfile'')','Backgroundcolor',[1 1 1]);
  bp(3) = 5*cw;
  bp(2) = bp(2) + bp(4)+VSPACE;
  bp(1) = SPACE;
  uicontrol(gcf,'Style', 'text','String','Expts','Position', bp * wsc);
  bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = 100;
  uicontrol(gcf,'Style', 'edit','String','','Position', bp * wsc,'Tag','Explist','Callback', ...
	    'psfgui(''update'')','Backgroundcolor',[1 1 1]);
    bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = cw * 8;
  uicontrol(gcf,'Style', 'text','String','Split By','Position', bp * wsc);
   
    bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = 4 * cw;
  uicontrol(gcf,'Style', 'edit','String',PSYCH.opt.vars{1},'Position', bp * wsc,'Tag','Vara','Callback', ...
	    'psfgui(''update'')','Backgroundcolor',[1 1 1]);
    
    bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = 4 * cw;
  uicontrol(gcf,'Style', 'edit','String',PSYCH.opt.vars{2},'Position', bp * wsc,'Tag','Varb','Callback', ...
	    'psfgui(''update'')','Backgroundcolor',[1 1 1]);

    bp(1) = bp(1)+bp(3)+SPACE;
  bp(3) = 4 * cw;
  uicontrol(gcf,'Style', 'edit','String',PSYCH.opt.vars{3},'Position', bp * wsc,'Tag','Varc','Callback', ...
	    'psfgui(''update'')','Backgroundcolor',[1 1 1]);

    hm = uimenu(gcf,'Label','File');
  uimenu(hm,'Label','Close','Callback','psfgui(''close'')');
  hm = uimenu(gcf,'Label','Options');
  uimenu(hm,'Label','Options','Callback','psfgui(''showoptions'')');
  uimenu(hm,'Label','Show Blocks','Callback','psfgui(''showblocks'')');
  uimenu(hm,'Label','Show Expts','Callback','psfgui(''showexpts'')');
  set(gcf,'Menubar','none');

fign = figure('Tag',PSYCH.figtag,'Position',[300 scrsz(4)-470 512 ...
		    400]);
PSYCH.id = 1;

function showexpts(blocklist, PSYCH)
  scrsz = get(0,'Screensize');
  wsc = scrsz(3)/1000;
  
  SPACE = 5;
  VSPACE = 2;
  BOXH = 15;
  
  len = length(blocklist.blocks);
  if (len < 20)
      h = (len+1) * (15 + VSPACE);
  else
      h = (20+1) * (15 + VSPACE);
  end

  CloseTag('Blocks');
  cntrl_box = figure('Position', [200 scrsz(4)-(h+50)*wsc 500 h*wsc], 'Menubar', 'none',...
       'NumberTitle', 'off', 'Tag','Blocks','Name','BlockList');
  
   bp = [30 0 500 15];
   for j = 1:length(blocklist.blocks);
      desc = [];
      for k = 1:length(PSYCH.opt.conditions)
          if isfield(blocklist,PSYCH.opt.conditions{k}) & eval(['~isempty(blocklist.' PSYCH.opt.conditions{k} '{j})']);
              desc = [desc ','];
              eval(['desc = [desc ' 'PSYCH.opt.conditions{k}' ' sprintf('' %d'',blocklist.' PSYCH.opt.conditions{k} '{j})];']);
          end
      end
      uicontrol(gcf,'Style', 'CheckBox','String',sprintf('%d %d %s %s %s',j,blocklist.n(j),blocklist.times{j},blocklist.names{j},desc),'Position', bp * wsc,...
      'Tag',sprintf('Block%d',j),'value',1);
      bp(2) = bp(2) + bp(4) + VSPACE;
  end
  bp(3) = 25;
  bp(1) = 1;
  bp(2) = h-20;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'psfgui(''checkblocks'')',...
'String', 'Go', 'Position', bp * wsc);
  bp(2) = bp(2)-bp(4);

  uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'psfgui(''clearblocks'')',...
'String', 'clr', 'Position', bp * wsc);
  bp(2) = bp(2)-bp(4);
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'psfgui(''closeblocks'')',...
'String', 'end', 'Position', bp * wsc);

  uicontrol(gcf, 'callback', 'psfgui(''scroll'')','style','slider','min',1','max',j,'value',1,...
      'position',[480 5 20 h-10],'Tag','blockscroll');
  

  
function showblocks(blocklist, PSYCH)
  scrsz = get(0,'Screensize');
  wsc = scrsz(3)/1000;
  
  SPACE = 5;
  VSPACE = 2;
  BOXH = 15;
  
  len = length(blocklist.blocks);
  if (len < 20)
      h = (len+1) * (15 + VSPACE);
  else
      h = (20+1) * (15 + VSPACE);
  end

  CloseTag('Blocks');
  cntrl_box = figure('Position', [200 scrsz(4)-(h+50)*wsc 500 h*wsc], 'Menubar', 'none',...
       'NumberTitle', 'off', 'Tag','Blocks','Name','BlockList');
  
   bp = [30 0 500 15];
   for j = 1:length(blocklist.blocks);
      desc = [];
      for k = 1:length(PSYCH.opt.conditions)
          if isfield(blocklist,PSYCH.opt.conditions{k}) & eval(['~isempty(blocklist.' PSYCH.opt.conditions{k} '{j})']);
              desc = [desc ','];
              eval(['desc = [desc ' 'PSYCH.opt.conditions{k}' ' sprintf('' %d'',blocklist.' PSYCH.opt.conditions{k} '{j})];']);
          end
      end
      uicontrol(gcf,'Style', 'CheckBox','String',sprintf('%d %d %s %s %s',j,blocklist.n(j),blocklist.times{j},blocklist.names{j},desc),'Position', bp * wsc,...
      'Tag',sprintf('Block%d',j),'value',1);
      bp(2) = bp(2) + bp(4) + VSPACE;
  end
  bp(3) = 25;
  bp(1) = 1;
  bp(2) = h-20;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'psfgui(''checkblocks'')',...
'String', 'Go', 'Position', bp * wsc);
  bp(2) = bp(2)-bp(4);

  uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'psfgui(''clearblocks'')',...
'String', 'clr', 'Position', bp * wsc);
  bp(2) = bp(2)-bp(4);
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'psfgui(''closeblocks'')',...
'String', 'end', 'Position', bp * wsc);

  uicontrol(gcf, 'callback', 'psfgui(''scroll'')','style','slider','min',1','max',j,'value',1,...
      'position',[480 5 20 h-10],'Tag','blockscroll');
  

function scrollblocks(listpos)

BOXH = 15;
  VSPACE = 2;
j = 1;
it = 1;
while ~isempty(it)
it = findobj('Tag',sprintf('Block%d',j));
if ~isempty(it)
    blocks(j) = it;
    j = j+1;
end
end
for j = 1:length(blocks)
    pos = get(blocks(j),'position');
    pos(2) = (j - (listpos)) * (BOXH + VSPACE);
    set(blocks(j),'position',pos);
end


function showoptions(PSYCH)

  scrsz = get(0,'Screensize');
  wsc = scrsz(3)/1000;

  legendvalues = [1 2 0];
if isempty(findobj('Tag',PSYCH.seltag))
  SPACE = 5;
  VSPACE = 2;
  scrsz = get(0,'Screensize');

  cntrl_box = figure('Position', [200 scrsz(4)-300 250 280], 'Menubar', 'none',...
       'NumberTitle', 'off', 'Tag',PSYCH.seltag,'Name','Section Criteria');
  bp = [0 0 100 15];
  uicontrol(gcf,'Style', 'CheckBox','String','Show Guess','Position', bp * wsc,...
      'Tag','ShowInitial','Callback','psfgui(''Update'')','value',PSYCH.opt.showinit);
  bp(2) = bp(2) + bp(4) + VSPACE;
  uicontrol(gcf,'Style', 'CheckBox','String','Label Mean','Position', bp * wsc,...
      'Tag','ShowMean','Callback','psfgui(''Update'')','value',PSYCH.opt.meanlabel);
  bp(2) = bp(2) + bp(4) + VSPACE;
  uicontrol(gcf,'Style', 'CheckBox','String','Label SD','Position', bp * wsc,...
      'Tag','ShowSD','Callback','psfgui(''Update'')','value',PSYCH.opt.sdlabel);
  bp(2) = bp(2) + bp(4) + VSPACE;
  uicontrol(gcf,'Style', 'CheckBox','String','Reread','Position', bp * wsc,...
      'Tag','ForceRead','Callback','psfgui(''Update'')','value',PSYCH.opt.forceread);
  bp(2) = bp(2) + bp(4) + VSPACE;
  uicontrol(gcf,'Style', 'CheckBox','String','Collapse','Position', bp * wsc,...
      'Tag','CollapseData','Callback','psfgui(''Update'')','value',PSYCH.opt.collapsedata);
  bp(2) = bp(2) + bp(4) + VSPACE;
  uicontrol(gcf,'Style', 'CheckBox','String','Hold','Position', bp * wsc,...
      'Tag','HoldOn','Callback','psfgui(''Update'')','value',PSYCH.opt.holdon);
 bp(2) = bp(2) + bp(4) + VSPACE;
  uicontrol(gcf,'Style', 'CheckBox','String','by Expt','Position', bp * wsc,...
      'Tag','BlockExpts','Callback','psfgui(''Update'')','value',PSYCH.opt.blockexpts);
   bp(2) = bp(2) + bp(4) + VSPACE;
  uicontrol(gcf,'Style', 'CheckBox','String','BadRate','Position', bp * wsc,...
      'Tag','BadRate','Callback','psfgui(''Update'')','value',PSYCH.opt.counterrs);
   bp(2) = bp(2) + bp(4) + VSPACE;
  uicontrol(gcf,'Style', 'CheckBox','String','Plot SD','Position', bp * wsc,...
      'Tag','PlotSD','Callback','psfgui(''Update'')','value',PSYCH.plotsds);
  bp(2) = bp(2) + bp(4) + VSPACE;
  uicontrol(gcf,'Style', 'CheckBox','String','Plot Mean','Position', bp * wsc,...
      'Tag','PlotMean','Callback','psfgui(''Update'')','value',PSYCH.plotmeans);
   bp(2) = bp(2) + bp(4) + VSPACE;
  uicontrol(gcf,'Style', 'CheckBox','String','Pecent Correct','Position', bp * wsc,...
      'Tag','PercentCorrect','Callback','psfgui(''Update'')','value',PSYCH.opt.pcorrect);
  bp(2) = bp(2) + bp(4) + VSPACE;
  uicontrol(gcf,'Style', 'text','String','Prefix','Position', bp * wsc);
  bp(1) = bp(1) + bp(3) + SPACE;
  uicontrol(gcf,'Style', 'edit','String','Prefix','Position', bp * wsc,...
      'Tag','Prefix','Callback','psfgui(''Update'')','String',PSYCH.opt.prefix);
end


function [Data, Trials, blocklist] = readpsychfile(PSYCH, name)

blocklist = [];
Trials = [];
argon = {};
matname = sprintf('%s.mat',name);
if exist(matname,'file') & ~PSYCH.opt.forceread
  load(matname);
  Data = CountPsychTrials(Trials,'sortby',PSCYH.opt.vars{1},'nmin', ...
			  PSYCH.opt.nmin);
  if exist('ExptTrialList','var')
      blocklist = ExptTrialList;
  end
  return;
end

if ~exist(name,'file')
  fprintf('Cant Read %s\n',name);
  Data = [];
  return;
end

fid = fopen(name,'r');
tline = fgetl(fid);
fclose(fid);
[a,b] = sscanf(tline,'%d %f');

if(strncmp(tline,'Reopened',8)) % this is a binoc output file
  if PSYCH.opt.counterrs
     argon = {argon{:} 'all'};
  end
  [Trials, blocklist] = readpsychtrials(name,'save',argon{:});
  if PSYCH.opt.counterrs
      for j = 1:length(Trials);
          Trials(j).score = 2 * (abs(Trials(j).score) -0.5);
      end
  end


  Data = CountPsychTrials(Trials,'sortby',PSYCH.opt.vars{1},'nmin', ...
			  PSYCH.opt.nmin,argon{:});
elseif length(a) & length(b)  %% a table of numbers, suggesting a online table
  tr = textread(name);
  id = find(tr(:,4) > 0 | tr(:,5) > 0)  %% made a response
  for j = 1:length(id)
      Trials(j).score = tr(id(j),4) - tr(id(j),5);
      Trials(j).sign = sign(tr(id(j),6));
      Trials(j).x = tr(id(j),6);
  end
  Data = CountPsychTrials(Trials,'nmin', ...
			  PSYCH.opt.nmin,'skip',PSYCH.opt.first,'last',PSYCH.opt.last);
else
  args = {};
  narg = 1;
  if ~isempty(PSYCH.doblocks)
    args{narg} = 'blocks';
    args{narg+1} = PSYCH.doblocks
    narg = narg+2;
  else
  if(PSYCH.opt.first > 1)
    args{narg} = 'skip';
    args{narg+1} = PSYCH.opt.first-1;
    narg = narg+2;
  end
  if(PSYCH.opt.last > 1)
    args{narg} = 'last';
    args{narg+1} = PSYCH.opt.last;
    narg = narg+2;
  end
  end
  if(PSYCH.opt.collapsedata)
    args{narg} = 'collapse';
    narg = narg+1;
  end
  
  nc = 0;
  if ~isempty(PSYCH.opt.vars{1})

    args{narg} = 'sort';
    narg = narg+1;
    args{narg} = PSYCH.opt.vars{1};
    nc = nc+1;
    PSYCH.opt.conditions{nc} = PSYCH.opt.vars{1};
    narg = narg+1;
  end
  if ~isempty(PSYCH.opt.vars{2})
    args{narg} = 'sort';
    narg = narg+1;
    args{narg} = PSYCH.opt.vars{2};
    nc = nc+1;
    PSYCH.opt.conditions{nc} = PSYCH.opt.vars{2};
    narg = narg+1;
  end
  if ~isempty(PSYCH.opt.vars{3})
    args{narg} = 'sort';
    narg = narg+1;
    args{narg} = PSYCH.opt.vars{3};
    nc = nc+1;
    PSYCH.opt.conditions{nc} = PSYCH.opt.vars{3};
    narg = narg+1;
  end
  
  [Data, blocklist] = readpsychsum(name,args{:});

end
