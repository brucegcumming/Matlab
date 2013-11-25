function DATA = PlotMap(varargin)

%
%  PlotMap('monkeyname')
%  reads in penetration data for a monkey, and plots RF locations, and
%  the use grid
%

name = 'RfMap';
% Set defaults before reading varargin
%
strings = [];
tag = 'RfMap'; %need to change these.
init = 0;
MapDefs;

if length(varargin) & isnumeric(varargin{1})
    toplevel = varargin{1};
    varargin = varargin(2:end);
elseif length(varargin) & strncmpi(varargin{1},'new',3)
    toplevel = [];
    varargin = varargin(2:end);
else
    toplevel = findobj('Tag',tag);
end

if ~isempty(toplevel)
  if strncmpi(varargin{1},'store',5)
      DATA = varargin{2};
    set(toplevel,'UserData',varargin{2});
    return;
  else
    DATA = get(toplevel,'UserData');
  end
else
    DATA.monkeynames = {'duf' 'ruf' 'ica' 'lem'};
    DATA.tags.toplevel = tag;
    DATA.plot.area = 1;
    DATA.plot.labelpts = 1;
    DATA.plot.popuprf = 0;
    DATA.plot.type = 1;
    DATA.plot.fontsiz = 10;
    DATA.plot.maxage = 0;
    DATA.plot.minage = 0;
    DATA.verbose = 0;
    DATA.plot.auto = 0;
    DATA.plot.hemisphere = 0;
    DATA.plot.selecttype = 0;
    DATA.reloadpens = 0;
    monkey = 1;
    init = 1;
end

j = 1;
if nargin 
    if ischar(varargin{1})
        if  strncmpi(varargin{1},'point',4)
            id = varargin{2};
            colors = mycolors;
            for j = 1:length(id)
                plotpen(DATA, DATA.px(id(j)), DATA.py(id(j)),'colors',colors);
            DATA.currentpen = [DATA.px(id(j)), DATA.py(id(j))];
            end
            PlotMap(DATA.top,'store',DATA);
        elseif  strncmpi(varargin{1},'cellpoint',6)
            id = varargin{2};
            pe = DATA.map.pen(id,1);
            pname = sprintf('/bgc/bgc/anal/icarus/pens/pen%d.log',pe);
            GetFigure('OnePen');
            subplot(1,1,1);
            hold off;
            type = get(findobj('Tag','PenPlot'),'value') -1;
            if type == 4
                PlotOnePen(pname,'plot',type,'cmtype',[0 1 4]);
            elseif type > 0
                PlotOnePen(pname,'plot',type,'cmtype',1);
            end
        elseif  strncmpi(varargin{1},'cell',4)
            id = varargin{2};
            colors = mycolors;
            for j = 1:length(id)
                fprintf('%s %.2f %.2f %.1f\n',DATA.map.cellname{id(j)},DATA.map.rf(id(j),1),...
                    DATA.map.rf(id(j),2),DATA.map.depth(id(j)));
            end
            PlotMap(DATA.top,'store',DATA);
        elseif  strncmpi(varargin{1},'excludepen',8)
            x = DATA.currentpen(1);
            y = DATA.currentpen(2);
            idx = find(DATA.map.pen(:,2) == x & DATA.map.pen(:,3) == y);
            DATA.excluded(idx) = 1;
            ReBuild(DATA);
        elseif  strncmpi(varargin{1},'missedpoint',8)
            id = find(DATA.map.pen(:,1) == varargin{2});
            fprintf('Missed on pen %d at %d,%d\n',varargin{2},DATA.map.pen(id,2), DATA.map.pen(id,3));
            PlotMap(DATA.top,'store',DATA);
        elseif  strncmpi(varargin{1},'gridpoint',6)
            id = varargin{2};
            plotpen(DATA, DATA.px(id), DATA.py(id),'grid');
            DATA.currentpen = [DATA.px(id), DATA.py(id)];
            PlotMap(DATA.top,'store',DATA);
        elseif  strncmpi(varargin{1},'penpoint',6)
            id = find(DATA.map.pen(:,1) == varargin{2});
            pe = DATA.map.pen(id(1),1);
            fprintf('Pen %d\n',pe);
            pname = sprintf('/bgc/bgc/anal/icarus/pens/pen%d.log',pe);
            GetFigure('OnePen');
            subplot(1,1,1);
            hold off;
            type = get(findobj('Tag','PenPlot'),'value') -1;
            if type == 4
                PlotOnePen(pname,'plot',type,'cmtype',[0 1 4]);
            elseif type > 0
                PlotOnePen(pname,'plot',type,'cmtype',1);
            end
        elseif  strncmpi(varargin{1},'depthpoint',4)
            ShowCell(DATA.map,varargin{2});
        elseif  strncmpi(varargin{1},'setmonkey',5)
            [DATA.map, DATA.monkey] = GetMap(DATA);
            DATA.selected = ones(size(DATA.map.area));
            DATA.excluded = zeros(size(DATA.map.area));
            DATA.marked = zeros(size(DATA.map.area));
            PlotMap(DATA.top,'store',DATA);
            RePlot(DATA);
        elseif  strncmpi(varargin{1},'savemap',7)
            SaveMap(DATA.map, DATA.map.monkeyname);
        elseif  strncmpi(varargin{1},'SetArea',5)
            DATA.plot.area = get(findobj('Tag','SetArea','Parent',toplevel),'value');
            set(toplevel,'UserData',DATA);
            PlotRFLoc(DATA);
        elseif  strncmpi(varargin{1},'getmap',5)
            return;
        elseif  strncmpi(varargin{1},'update',5)
            DATA = update(DATA);
            if DATA.plot.auto
                ReBuild(DATA);
            end
        elseif  strncmpi(varargin{1},'plottype',5)
            DATA = update(DATA);
            ReBuild(DATA);          
        elseif  strncmpi(varargin{1},'dufus',5)
            init = 1;
            monkey = 1;
        elseif  strncmpi(varargin{1},'rufus',5)
            init = 1;
            monkey = 2;
        elseif  strncmpi(varargin{1},'icarus',5)
            init = 1;
            monkey = 3;
        elseif  strncmpi(varargin{1},'lem',3)
            init = 1;
            monkey = 4;
        elseif  strncmpi(varargin{1},'dae',3)
            init = 1;
            monkey = 5;
        elseif  strncmpi(varargin{1},'testmap',5)
            init = 1;
            monkey = 5;
        elseif strncmpi(varargin{1},'scroll',5)
            it = findobj('Tag','blockscroll','Parent',DATA.top);
            j = get(it,'value');
            scrolllist(j);
        elseif strncmpi(varargin{1},'checkblocks',5)
            [a, iset, values] = GetBlockList('Block');
            id = find(iset == 0);
            DATA.excluded(values(id)) = 1;
            PlotMap(DATA.top,'store',DATA);
            ReBuild(DATA);    
        elseif strncmpi(varargin{1},'clearblocks',5)
            DATA.excluded(find(DATA.excluded ==1)) = 0;
            PlotMap(DATA.top,'store',DATA);
            ReBuild(DATA);
        elseif strncmpi(varargin{1},'closeblocks',8)
            CloseTag('Blocks');
        elseif strncmpi(varargin{1},'choosepens',5)
            ChoosePens(DATA);
        elseif strncmpi(varargin{1},'markpen',5)
            x = DATA.currentpen(1);
            y = DATA.currentpen(2);
            idx = find(DATA.map.pen(:,2) == x & DATA.map.pen(:,3) == y);
            DATA.marked(idx) = 1;
            idx = find(DATA.xc == x & DATA.yc == y);
            DATA.penmarked(idx) = 1;
            PlotMap(DATA.top,'store',DATA);
        elseif strncmpi(varargin{1},'unmarkpen',5)
            x = DATA.currentpen(1);
            y = DATA.currentpen(2);
            idx = find(DATA.map.pen(:,2) == x & DATA.map.pen(:,3) == y);
            DATA.marked(idx) = 0;
        elseif strncmpi(varargin{1},'close',5)
            CloseTag(DATA.tags.fig);
            close(DATA.top);
        elseif strncmpi(varargin{1},'plot',4)
            DATA = ReBuild(DATA);    
            RePlot(DATA);
        elseif strncmpi(varargin{j},'print',5)
            PrintPens(DATA,DATA.map);
        else
            j = 1;
            init = 1;
            while(j < nargin)
                if(strncmpi(varargin{j},'name',3))
                    j = j+1;
                    name = varargin{j};
                end
                j = j+1;
            end
        end
    elseif iscell(varargin{j})
        DATA.map = rflist2map(DATA, varargin{j});
        DATA.monkey = DATA.map.monkey;
        DATA.map = ReadAllPens(DATA.map,DATA.map.monkeyname,'reload');
        init = 1;
    end
end

showgrid = 0;
j = 2;
while j <= length(varargin)
    if strncmpi(varargin{j},'right',5) || strcmp(varargin{j},'R')
            DATA.plot.hemisphere = 1;
    elseif  strncmp(varargin{j},'grid',4)
        showgrid = 1;
    elseif  strncmpi(varargin{j},'maxage',4)
        j = j+1;
        DATA.plot.maxage = varargin{j};
    elseif  strncmp(varargin{j},'MT',2)
        DATA.plot.area = MT;
    elseif  strncmp(varargin{j},'reload',4)
        DATA.reloadpens =1;
    end
    j = j+1;
end
if init & isempty(toplevel)
    scrsz = get(0,'Screensize');
    wsc = scrsz(3) /1280;  %scale everything as is 1280 wide
    if scrsz(3) > 2400
        wsc = 1.2;
    end
    wsize(1) = 380 * wsc;
    wsize(2) = 200 * wsc;
    
    cntrl_box = figure('Position', [10 scrsz(4)-220*wsc 300*wsc 200*wsc],...
        'NumberTitle', 'off', 'Tag',tag,'Name',name);
    DATA.top = cntrl_box;
    top = num2str(DATA.top);
    DATA.tags.fig = ['RFMapPlot' top];
    DATA.tags.figb = ['RFGridPlot' top];
    if( ~isempty(strings))
        lst = uicontrol(gcf, 'Style','listbox','String', strings,...
            'Callback', ' runcp(''setentry'')','Tag',listtag,...
            'Position',[10 10 wsize(1) wsize(2)]);
        
    end
    HSPACE = 3;
    VSPACE = 2;
    cw = 9;
    ch = 18*wsc + VSPACE;
    
    bp(1) = HSPACE; bp(3) = 25; bp(2) = wsize(2)-ch; bp(4) = 22;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['PlotMap(' top ',''next'');'],...
        'String', '>>', 'Position', bp);
    bp(1) = bp(1) + bp(3) + HSPACE;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', [' PlotMap(' top ',''next'');'],...
        'String', '>>', 'Position', bp);
    
 %New row
    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    bp(3) = 6 * cw * wsc;
    uicontrol(gcf,'Style', 'checkbox','String', 'Label Pts', 'Tag', 'LabelPts', 'Callback', [' PlotMap(' top ',''update'');'],...
        'Position', bp,'value',DATA.plot.labelpts);
    bp(1) = bp(1) + bp(3)+ HSPACE;
    uicontrol(gcf,'Style', 'checkbox','String', 'Pop RF', 'Tag', 'PopupRF', 'Callback', [' PlotMap(' top ',''update'');'],...
        'Position', bp);
    
    bp(2) = bp(2) - ch;
    bp(3) = 6 * cw * wsc;
    bp(1) = HSPACE;
    uicontrol(gcf,'Style', 'pushbutton','String','Print','Callback',['PlotMap(' top ',''print'');'],'Position', bp);
    bp(1) = bp(1) + bp(3)+ HSPACE;
    uicontrol(gcf,'Style', 'pushbutton','String','Plot','Callback',['PlotMap(' top ',''plot'');'],'Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'style','pop','string','Map|Grid|GridVar|GridDepth|Contour|Xpcolor|Ypcolor|dZdX|Arrows|MeanArrow|GMdepth|GMsmooth', ...
           'Callback', [' PlotMap(' top ',''plottype'');'], 'Tag','plottype',...
        'position',bp,'value',1); bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'text','String','Monkey','Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'style','pop','string','Dufus|Rufus|Icarus|Lemieux|Daedalus|Test', 'Tag','RFMonkeyName','Callback',['PlotMap(' top ',''setmonkey'');'], 'position',bp,'value',monkey);
  
    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    uicontrol(gcf,'Style', 'text','String','Area','Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'style','pop','string',areanames, 'Tag','SetArea','Callback',['PlotMap(' top ',''SetArea'');'], ...
        'value',DATA.plot.area,'position',bp);
    bp(1) = bp(1) + bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'text','String','HemiSphere','Position', bp,'value',DATA.plot.hemisphere+1);
    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = cw*4;
    uicontrol(gcf,'style','pop','string','L|R','Tag','HemiSphere','Callback',['PlotMap(' top ',''update'');'], 'position',bp,...
        'value',DATA.plot.hemisphere+1);
    bp(1) = bp(1) + bp(3)+HSPACE;
    bp(3) = cw * 6;
    uicontrol(gcf,'Style', 'text','String','Fontsiz','Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'style','pop','string','4|6|8|10|12', 'value',4,'Tag','FontSize','Callback',['PlotMap(' top ',''update'');'], 'position',bp);
    

    
 
    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    uicontrol(gcf,'Style', 'text','String','Max Age','Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.2f',DATA.plot.maxage),'Position', bp,'Tag','MaxAge','Callback', ...
	    ['PlotMap(' top ',''Update'')'],'Backgroundcolor',[1 1 1]);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.2f',0.0),'Position', bp,'Tag','MinAge','Callback', ...
	    ['PlotMap(' top ',''Update'')'],'Backgroundcolor',[1 1 1]);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'CheckBox','String','Verbose','Tag','Verbose',...
        'Value',DATA.verbose,'Callback',' PlotMap(DATA.top,''update'')','Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'CheckBox','String','Auto','Tag','Auto',...
        'Value',DATA.plot.auto,'Callback',' PlotMap(DATA.top,''update'')','Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'style','pop','string','None|Depths|RFs|3DRF|Comments', 'value',4,'Tag','PenPlot', 'position',bp);

    
    
    
    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    uicontrol(gcf,'Style', 'text','String','Select','Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = cw*4;
    uicontrol(gcf,'style','pop','string','All|MultContact', 'value',1,'Tag','SelectType','Callback',['PlotMap(' top ',''update'');'], 'position',bp);
    
  hm = uimenu(gcf,'Label','File');
  top = DATA.top;
  uimenu(hm,'Label','Close','Callback',[' PlotMap(' num2str(DATA.top) ',''close'');']);
  uimenu(hm,'Label','Reload','Callback',[' PlotMap(' num2str(DATA.top) ',''setmonkey'');']);
  uimenu(hm,'Label','Choose','Callback',[' PlotMap(' top ',''choosepens'');']);
  uimenu(hm,'Label','Exclude','Callback',[' PlotMap(' top ',''excludepen'');']);
  uimenu(hm,'Label','Mark','Callback',[' PlotMap(' top ',''markpen'');']);
  uimenu(hm,'Label','UnMark','Callback',[' PlotMap(' top ',''unmarkpen'');']);
  uimenu(hm,'Label','Save Map','Callback',[' PlotMap(' num2str(top) ',''savemap'');']);
    set(gcf,'Menubar','none');
    if ~isfield(DATA,'map')
    [DATA.map, DATA.monkey] = GetMap(DATA);
    end
    DATA.selected = ones(size(DATA.map.area));
    DATA.excluded = zeros(size(DATA.map.area));
    DATA.marked = zeros(size(DATA.map.area));
    DATA = update(DATA);
    DATA = ReBuild(DATA);
    set(gcf,'UserData',DATA);
    if showgrid
        GetFigure(DATA.tags.figb);
        hold off;
        PlotGrid(DATA);
    end
end


function map = rflist2map(DATA, R)
MapDefs;
name2area = [1 1 2 1 10];

for j = 1:length(R)
    map.rf(j,:) = R{j}.rf;
    map.cellname{j} = R{j}.name;
    map.area(j) = 1;
    if isfield(R{j},'area')
        a = find(strcmp(R{j}.area,{'V1' 'Vd' 'V2' 'unknown' 'Vc'}));
        if ~isempty(a)
            map.area(j) = name2area(a);
        end
    end
    map.hemisphere(j) = 0;
    if isfield(R{j},'date')
        map.datestr{j} = datestr(R{j}.date);
    else
        map.datestr{j} = 'unknown';
    end
    map.depth(j) = R{j}.depth;
    monkey{j} = GetMonkeyName(R{j}.name);
    if strcmp(R{j}.electrode,'uProbe')
        map.types(j) = MULTICONTACT;
    else
        map.types(j) = NORMAL;
    end
end
[a,b] = Counts(monkey);
[c,d] = max(a);
map.monkey = b{d};
map.monkeyname = b{d};
txt = scanlines(['/bgc/anal/' map.monkey '/pens.err']);
for j = 1:length(txt)
    a = sscanf(txt{j},'%f');
    if length(a) > 2
        id = find(map.rf(:,6) == a(1));
        map.rf(id,7) = a(2);
        map.rf(id,8) = a(3);
    end
end
   
map.monkey = find(strcmp(map.monkeyname, DATA.monkeynames));
map.pen = map.rf(:,6:8);

function [map, monkey] = GetMap(DATA, varargin)
MapDefs;

global bgcfileprefix
monkey =  get(findobj('Tag','RFMonkeyName','Parent',DATA.top),'value');
nextfile = [];
if monkey == DUFUS
  mapfile = '/bgc/bgc/anal/dufus/dufv1.fixtab';
  missfile = '/bgc/bgc/anal/dufus/missed.pens';
  monkeyname = 'dufus';
elseif monkey == RUFUS
  mapfile = '/bgc//bgc/anal/rufus/rufv1.fixtab';
  missfile = '/bgc/bgc/anal/rufus/missed.pens';
  monkeyname = 'rufus';
elseif monkey == ICARUS
  mapfile = '/bgc//bgc/anal/icarus/icarus.fixtab';
  missfile = '/bgc/bgc/anal/icarus/missed.pens';
  nextfile = '/bgc/bgc/anal/icarus/penlist';
  monkeyname = 'icarus';
elseif monkey == LEMIEUX
  mapfile = '/bgc/bgc/anal/lem/lem.fixtab';
  missfile = '/bgc/bgc/anal/lem/missed.pens';
  nextfile = '/bgc/bgc/anal/lem/penlist';
  monkeyname = 'lem';
elseif monkey == DAE
  mapfile = '/bgc/bgc/anal/dae/dae.fixtab';
  missfile = '/bgc/bgc/anal/dae/missed.pens';
  nextfile = '/bgc/bgc/anal/dea/penlist';
  monkeyname = 'lem';

elseif monkey == TESTMAP
  mapfile = '/bgc//bgc/anal/maps/test.tab';
  missfile = '/bgc/bgc/anal/lem/missed.pens';
end
DATA.monkey = monkey;
if exist('bgcfileprefix','var')
%    mapfile = [bgcfileprefix mapfile];
end

map = loadmap(mapfile);
map.monkeyname = monkeyname'
if isfield(DATA,'map') & isfield(DATA.map, 'pens')
    map.pens = DATA.map.pens;
    map.missed = DATA.map.missed;
else
    args = {};
    if DATA.reloadpens
        args = {args{:} 'reload'};
    end
    map = ReadAllPens(map,map.monkeyname,args);
end

if monkey == RUFUS
    id = find(map.datenum < datenum('9/1/2006'));
    map.hemisphere(id) = 0;
    id = find(map.datenum > datenum('9/1/2006'));
    map.hemisphere(id) = 1;
end
if exist(missfile) & isnumeric(map.missed) & sum(map.missed) == 0
    [a,b] = textread(missfile,'%d %s');
    v = strmatch('V2',b);
    if ~isempty(v)
        missed{2} = a(v);
    end
    if isempty(a)
        missed{1} = [];
        missed{2} = [];
    end
        map.missed = missed;
end
if exist(nextfile)
    fid = fopen(nextfile,'r');
    a = textscan(fid,'%s','delimiter','\n');
    txt = a{1};
    for j = 1:length(txt)
        crd(j) = textscan(txt{j},'%n');
        nf(j) = length(crd{j});
        if nf(j) <= 2 | crd{j}(3) == 0
            unused(j) = 1;
        else
            unused(j) = 0;
        end
    end
    id = find(unused)
    if length(id)
        for j = 1:length(id)
        map.nextpens(1:length(crd{id(j)}),j) = crd{id(j)};
        end
    end
end

function DATA = update(DATA, varargin)
  DATA.plot.labelpts = get(findobj('Tag','LabelPts','Parent',DATA.top),'value');
  DATA.plot.popuprf = get(findobj('Tag','PopupRF','Parent',DATA.top),'value');
  DATA.plot.type = get(findobj('Tag','plottype','Parent',DATA.top),'value');
  DATA.plot.area = get(findobj('Tag','SetArea','Parent',DATA.top),'value');
  DATA.plot.selecttype = get(findobj('Tag','SelectType','Parent',DATA.top),'value');
  DATA.plot.hemisphere = get(findobj('Tag','HemiSphere','Parent',DATA.top),'value')-1;
  str = get(findobj('Tag','FontSize','Parent',DATA.top),'string');
  j = get(findobj('Tag','FontSize','Parent',DATA.top),'value');
  DATA.plot.fontsiz = str2num(str(j,:));
  DATA.verbose = get(findobj('Tag','Verbose','Parent',DATA.top),'value');
  DATA.plot.auto = get(findobj('Tag','Auto','Parent',DATA.top),'value');

  
  DATA.plot.maxage = str2num(get(findobj('Tag','MaxAge','Parent',DATA.top),'string'));
  DATA.plot.minage = str2num(get(findobj('Tag','MinAge','Parent',DATA.top),'string'));

  PlotMap(DATA.top,'store',DATA);

  
function CloseTag(tag)
it = findobj('Tag',tag);
if ~isempty(it)
    close(it);
end


function DATA = ReBuild(DATA);
    DATA = PlotRFLoc(DATA);
    
function DATA = RePlot(DATA)

MapDefs;
if DATA.plot.area == COMPAREV1V2
    gridlines = 0;
else
    gridlines = 1;
end
    
GetFigure(DATA.tags.fig);
hold off;
if DATA.plot.type == 2
    GetFigure(DATA.tags.figb);
    hold off;
    PlotGrid(DATA);
    return;
elseif DATA.plot.type == 3
    PlotGridVar(DATA);
    return;
elseif DATA.plot.type == 9
    PlotGridArrows(DATA);
    return;
elseif DATA.plot.type == 10
    PlotGridArrows(DATA,'mean');
    return;
elseif DATA.plot.type == 11
    PlotGMDepths(DATA);
    return;
elseif DATA.plot.type == 12
    PlotGMDepths(DATA, 'smooth', 0.3);
    return;
elseif DATA.plot.type == 4
    PlotGridDepths(DATA);
    return;
elseif ismember(DATA.plot.type, [5 6 7]) %contour/pcolor plot
    PlotRFContour(DATA, DATA.plot.type);
elseif DATA.plot.type == 8
    PlotZX(DATA);
else
    PlotRFMap(DATA, gridlines);
    offset = range(DATA.xc)/60;
    if isfield(DATA.map,'nextpens') & length(DATA.pid) > 10
        recent = (length(DATA.pid)-10):length(DATA.pid);
        id = find(~ismember(DATA.map.nextpens(1,:)+i*DATA.map.nextpens(2,:),DATA.ppx(recent)+i*DATA.ppy(recent)));
        for j = id;
        x = PredictRF(DATA, DATA.map.nextpens(:,j));
        plot(x(1),x(2),'go');
        text(x(1)+offset,x(2),sprintf('%d,%d',DATA.map.nextpens(1,j),DATA.map.nextpens(2,j)),'color','g');
        fprintf('Pen at %.1f,%.1f RFs should be %.1f,%.1f\n',DATA.map.nextpens(1,j),DATA.map.nextpens(2,j),x(1),x(2));
        end
    end
end
 
function px = PredictRF(DATA, x)
sd = 1;
C = (x(1) +i * x(2)) - (DATA.px + i * DATA.py);
w = exp(-(abs(C).^2)/2.*sd.^2);
px(1) = sum(w.*DATA.xc)./sum(w);
px(2) = sum(w.*DATA.yc)./sum(w);

function PlotZX(DATA)
plotdz = 6;
usefits = 1;
top = num2str(DATA.top);
pens = unique(DATA.map.pen(:,1));
dzs = [];
dxs = [];
dys = [];
deccs = [];
dangles = [];
ids = [];
for j = 1:length(pens)
    id = find(DATA.map.pen(:,1) == pens(j));
    dz = DATA.map.depth(id);
    [mz, mzi] = min(dz);
    dz = dz -mz;
    if usefits
        dx = DATA.map.rf(id,6)- DATA.map.rf(id(mzi),6);
        dy = DATA.map.rf(id,7)- DATA.map.rf(id(mzi),7);
        ecc = abs(DATA.map.rf(id,6)+ i * DATA.map.rf(id,7));
        angles = atan2(DATA.map.rf(id,6),DATA.map.rf(id,7));
    else
        dx = DATA.map.rf(id,1)- DATA.map.rf(id(mzi),1);
        dy = DATA.map.rf(id,2)- DATA.map.rf(id(mzi),2);
        ecc = abs(DATA.map.rf(id,1)+ i * DATA.map.rf(id,2));
        angles = atan2(DATA.map.rf(id,2),DATA.map.rf(id,1));
    end
    decc = ecc - ecc(mzi);
    dangle = angles - angles(mzi);
    for k = 1:length(id)
    if plotdz == 1
    h = plot(dz(k),dangle(k),'o');
    hold on;
    elseif plotdz == 6
        h = plot(dz(k),decc(k),'o');
        hold on;
    elseif plotdz == 5
        h = plot(dz(k),dx(k)-dy(k),'o');
        hold on;
    else
    h = plot(dz(k),dx(k),'o');
    hold on;
    h = plot(dz(k),dy(k),'ro');
    end
    set(h,'buttondownfcn',['PlotMap(' top ',''cell'',' num2str(id(k)) ',1);']);
    end
    if length(dz) > 1
        dangles = [dangles dangle(2:end)'];
        deccs = [deccs decc(2:end)'];
        ids = [ids id(2:end)'];
        dzs = [dzs dz(2:end)];
        dxs = [dxs dx(2:end)'];
        dys = [dys dy(2:end)'];
    end
end
[x,y] = meshgrid([min(dxs):0.01:max(dxs)],[min(dys):0.01:max(dys)]);
Z = Interpf(dxs, dys, dzs,x,y,1,0.01);
hold off;
if plotdz == 3
imagesc(x(1,:),y(:,1),Z);
elseif plotdz == 4
    plot3(dxs,dys,dzs,'o');
end
    

function DATA = PlotRFLoc(DATA, varargin)

MapDefs;
showall = 0;
holdon = 0;
map = DATA.map;
GetFigure(DATA.tags.fig);
area = DATA.plot.area;

j = 1;
while(j < nargin)
    if strncmpi(varargin{j},'hold',4)
        holdon = 1;
    elseif strncmpi(varargin{j},'all',3)
        showall = 1;
    elseif strncmpi(varargin{j},'area',3)
        area = varargin{j+1};
        j = j+1;
    end
    j = j+1;
end

gridlines = 1;
if area == COMPAREV1V2
    gridlines = 0;
    idx = find(map.area ~=1);
    idx = find(map.area == V2);
    idxall = [];
    for j = idx
        v1 = find(map.area' == V1 & map.pen(:,2) == map.pen(j,2) & map.pen(:,3) == map.pen(j,3));
        idxall = [idxall v1'];
    end
    DATA.selected = zeros(size(map.area));
    DATA.selected(unique(idxall)) = 1;
    DATA.selected(idx) = 1;
elseif area ~= ALL_AREAS
    idx = find(map.area == area & map.hemisphere == DATA.plot.hemisphere);
    DATA.selected = zeros(size(map.area));
    DATA.selected(idx) = 1;
    if area <= length(map.missed)
        if iscell(map.missed)
            DATA.missed = find(map.missed{area} == area);
        else
            DATA.missed = find(map.missed == area);
        end
    else
        DATA.missed = [];
    end
else
    idx = 1:length(map.area);
    DATA.selected = ones(size(map.area));
    DATA.missed = [map.missed{1} map.missed{2}];
end

if DATA.monkey == DUFUS
    idx = find(map.pen(:,1) < 2);
else
    idx = find(map.pen(:,1) < 1);
end

    
DATA.selected(idx) = 0;
idx = find(DATA.excluded ==1);
DATA.selected(idx) = 0;

if DATA.plot.maxage > 0
    id = find(map.age > DATA.plot.maxage);
    DATA.selected(id) = 0;
end
if DATA.plot.minage > 0
    id = find(map.age < DATA.plot.minage);
    DATA.selected(id) = 0;
end

if DATA.plot.selecttype ==MULTICONTACT
    id = find(map.types ~= MULTICONTACT)
    DATA.selected(id) = 0;
end
    
PlotMap(DATA.top,'store',DATA);

if DATA.plot.type == 4
    PlotGridDepth(DATA);
    return;
end

cellsonly = 1;

idx = find(DATA.selected > 0);
selected = idx;
if isempty(idx)
    fprintf('No data for Area %d\n',area);
    return;
end
tic;

for j = idx;
    xpos(j) = map.rf(j,1);
    ypos(j) = map.rf(j,2);
    if map.pen(j,2) > -99 & DATA.selected(j)
        all.px(j) = map.pen(j,2);
        all.py(j) = map.pen(j,3);
        all.area(j) = map.area(j);
        all.id(j) = map.pen(j,1);
        if regexp(map.cellname{j},'M[0-9][0-9][0-9]')
            all.ptype(j) = 1;
        elseif regexp(map.cellname{j},'S[0-9][0-9][0-9]')
            all.ptype(j) = 1;
        else
            all.ptype(j) = 0;
        end
    else
        all.px(j) = NaN;
        all.py(j) = NaN;
        all.area(j) = NaN;
        all.id(j) = NaN;
    end
    if strcmp(map.cellname{j},'Nocell')
        all.cell(j) = 0;
    else
        all.cell(j) = 1;
    end
end


toc;
idx = find(DATA.selected == 0);
all.px(idx) = NaN;
all.py(idx) = NaN;
all.cell(idx) = NaN;
all.id(idx) = NaN;

if ~holdon
    hold off;
end

if showall
    plot(xpos,ypos,'o');
    hold on;
end
np = 1;



%Collect data together by penetration location
DATA.xc = [];
DATA.yc = [];
DATA.px = [];
DATA.py = [];
DATA.rvar = [];
DATA.id = [];
DATA.ptype = []
%ids = DATA.id;
xvals = all.px(find(~isnan(all.px)));
yvals = all.py(find(~isnan(all.py)));
tic;
for x= unique(xvals)
    rowstart = np;
    for y= unique(yvals);
        idx = find(all.px == x & all.py == y & all.cell > 0);
        if isempty(idx) && ~cellsonly
            idx = find(all.px == x & all.py == y);
        end
        if(length(idx) > 0)
            
            if area == COMPAREV1V2
                v1 = find(all.area(idx) == 1);
                v2 = find(all.area(idx) ~= 1);
                DATA.xc(np) = mean(xpos(idx(v1)));
                DATA.yc(np) = mean(ypos(idx(v1)));
                DATA.axc(np) = mean(xpos(idx(v2)));
                DATA.ayc(np) = mean(ypos(idx(v2)));
            else
                DATA.xc(np) = mean(xpos(idx));
                DATA.yc(np) = mean(ypos(idx));
            end
                rad = sqrt((xpos(idx) - DATA.xc(np)) .^2 + (ypos(idx) - DATA.yc(np)).^2);
            if length(rad) > 1
                DATA.rvar(np) = std(rad);
            else
                DATA.rvar(np) = 0;
            end
%            DATA.id(np) = ids(idx);
            if sum(DATA.marked(idx)) > 0
                DATA.penmarked(np) = 1;
            else
                DATA.penmarked(np) = 0;
            end      
            DATA.id(np) = max(all.id(idx));
            DATA.allid{np} = unique(all.id(idx));
                
            DATA.px(np) = x;
            DATA.py(np) = y;
            if sum(map.types(idx) == MULTICONTACT) % some multicontact probes 
                DATA.ptype(np) = MULTICONTACT;
            else
                DATA.ptype(np) = 0;
            end
            np = np+1;
            if(DATA.verbose)
                fprintf('%.0f %.0f %d\n',x,y,np);
            end
        end
    end
%    plot(DATA.xc(rowstart:np-1),DATA.yc(rowstart:np-1),':');
end
np = 0;
for pid= unique(all.id)
    rowstart = np;
    idx = find(all.id == pid);
        if(length(idx) > 0) & pid> 0 & sum(all.cell(idx)) > 0
            np = np+1;
                DATA.pxc(np) = mean(xpos(idx));
                DATA.pyc(np) = mean(ypos(idx));
                DATA.pid(np) = pid;
            DATA.ppx(np) = mean(all.px(idx));
            DATA.ppy(np) = mean(all.py(idx));
            pxs = unique(all.px(idx));
            pys = unique(all.py(idx));
            if length(pxs) >1 | length(pys) > 1
                fprintf('Pen %d 2 Co-ordinates: px ',pid);
                fprintf('%.1f,',pxs);
                fprintf(' py: ');
                fprintf('%.1f,',pys);
                fprintf('\n');
            end
            DATA.ppid(np) = find(DATA.px == all.px(idx(1)) & DATA.py == all.py(idx(1)));;
        end
end
    
toc;


PlotMap(DATA.top,'store',DATA);

RePlot(DATA);


function PlotRFContour(DATA, type)
[X,Y] =meshgrid(floor(min(DATA.px)):ceil(max(DATA.px)),floor(min(DATA.py)):ceil(max(DATA.py)));
for j = 1:size(X,1)
    for k = 1:size(X,2)
        id = find(abs(DATA.px - X(j,k)) < 0.7 & abs(DATA.py - Y(j,k)) < 0.7);
        if ~isempty(id)
            Z(j,k) = mean(DATA.xc(id));
            ZY(j,k) = mean(DATA.yc(id));
        else
            Z(j,k) = NaN;
            ZY(j,k) = NaN;
        end
    end
end

if type == 5
yrange = [floor(min(min(Y))) ceil(max(max(Y)))];
xrange = [floor(min(min(X))) ceil(max(max(X)))];
[Xi, Yi] = meshgrid(xrange(1):0.1:xrange(2), yrange(1):0.1:yrange(2));
Zi = Interpf(X,Y, Z, Xi, Yi, 1, 0.5);
    [c, h] = contour(X,Y,Z,floor(min(min(Z))):ceil(max(max(Z))),'r');
clabel(c,h,'color','r');
hold on;
[c, h] = contour(X,Y,ZY,floor(min(min(ZY))):ceil(max(max(ZY))),'b');
clabel(c,h,'color','b');

elseif type == 6
    GetFigure('XPcolor','parent',DATA.top);
    [X,Y,Z] = fillpmesh(X,Y,Z)
    h = pcolor(X,Y,Z);
    set(h,'buttondownfcn',@HitImage);
    colorbar;
elseif type == 7
    GetFigure('YPcolor','parent',DATA.top);
    [X,Y,Z] = fillpmesh(X,Y,ZY)
    h = pcolor(X,Y,Z);
    set(h,'buttondownfcn',@HitImage);
    colorbar;
end

function HitImage(a,b)

DATA = GetDataFromFig(a);
c = get(gca,'CurrentPoint');
px = floor(c(1));
py = floor(c(1,2));
fprintf('At %d,%d Pos = %.1f\n',px,py,c(1,3));
plotpen(DATA,px,py);

function PlotRFMap(DATA,gridlines)


MapDefs;
GetFigure(DATA.tags.fig);
hold off;
area = DATA.plot.area;
top = num2str(DATA.top);

if gridlines
%Draw blue lines for rows, and Red lines connecting columns
for x= unique(DATA.px)
    idx = find(DATA.px == x);
    plot(DATA.xc(idx),DATA.yc(idx),':');
    hold on;
end
for y= unique(DATA.py)
    idx = find(DATA.py == y);
    plot(DATA.xc(idx),DATA.yc(idx),'r:');
end
end

offset = range(DATA.xc)/60;

for j = 1:length(DATA.xc)
    if area == COMPAREV1V2
        plot([DATA.xc(j) DATA.axc(j)],[DATA.yc(j) DATA.ayc(j)],'-');
        hold on;
        plot(DATA.axc(j),DATA.ayc(j),'ro','buttondownfcn',['PlotMap(' top ',''point'',' num2str(j) ',1);']);
    end
    h = plot(DATA.xc(j),DATA.yc(j),'o','buttondownfcn',['PlotMap(' top ',''point'',' num2str(j) ',1);']);
    if DATA.penmarked(j)
        set(h, 'MarkerFaceColor','b');
    elseif DATA.ptype(j) == MULTICONTACT
        set(h, 'MarkerFaceColor','r');
    end
    
    if DATA.plot.labelpts
        if mod(DATA.px(j),1) > 0.1 | mod(DATA.py(j),1) > 0.1
            text(DATA.xc(j)+offset,DATA.yc(j),sprintf('%.1f,%.1f',DATA.px(j),DATA.py(j)),'fontsiz',DATA.plot.fontsiz);
        else
            text(DATA.xc(j)+offset,DATA.yc(j),sprintf('%.0f,%.0f',DATA.px(j),DATA.py(j)),'fontsiz',DATA.plot.fontsiz);
        end
    end
end

axis('equal');



function PlotGrid(DATA)

MapDefs;

top = num2str(DATA.top);


for j = 1:length(DATA.xc)
    if DATA.plot.labelpts
        text(DATA.px(j),DATA.py(j),sprintf('%.1f,%.1f',DATA.xc(j),DATA.yc(j)),'fontsiz',DATA.plot.fontsiz,'HorizontalAlignment','center','VerticalAlign','bottom');
    end
    h = plot(DATA.px(j),DATA.py(j),'o','buttondownfcn',['PlotMap(' top ',''gridpoint'',' num2str(j) ',1);'],'MarkerFaceColor','b');
    if DATA.penmarked(j)
        set(h, 'MarkerFaceColor','m');
    elseif DATA.ptype(j) == MULTICONTACT
        set(h, 'MarkerFaceColor','r');
    end
    hold on;
end
[sorted, idx] = sort(DATA.pid);
colors = {[0 0 0], [0.3 0 0], [0.6 0.3 0] [0.6 0 0], [1 1 0 ], [1 0 0]} ;

k = 1;
if length(idx) > 5
    back = 5;
else
    back = length(idx)-1;
end
for j = idx(end-back:end)
    hs(k) = plot(DATA.ppx(j),DATA.ppy(j),'o','buttondownfcn',['PlotMap(' top ',''gridpoint'',' num2str(DATA.ppid(j)) ',1);'],'color',colors{k},'MarkerFaceColor',colors{k});
    labels{k} = sprintf('%d',DATA.pid(j));
    k = k+1;
end
legend(hs,labels);
for j =1:length(DATA.missed)
    k = find(DATA.map.pen(:,1) == DATA.missed(j));
    h = plot(DATA.map.pen(k,2),DATA.map.pen(k,3),'o','buttondownfcn',['PlotMap(' top ',''missedpoint'',' num2str(j) ',1);'],'color','k','MarkerFaceColor','none');
end


xr = get(gca,'Xlim');
yr = get(gca,'Ylim');
set(gca,'Xtick',[xr(1):xr(2)],'YTick',[yr(1):yr(2)],'Xgrid','on','Ygrid','on');



function PlotGMDepths(DATA, varargin)

np = 0;
smooth = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'smooth',4)
        j = j+1;
        smooth = varargin{j};
    end
    j = j+1;
end
for j = 1:length(DATA.map.pens)
    pen = DATA.map.pens{j};
    if isfield(pen,'cmtype')
    cid = find(pen.cmtype == 6);
    mtcells = [];
    n = 0;
    for k = 1:length(cid)
        if strfind(pen.comments{cid(k)},'MT')
            n = n+1;
            mtcells(n) = pen.cmtime(cid(k));
        end
    end
    if n
        id = find(pen.cmtype == 5);
        gmid = find(pen.cmtime(id) <= mtcells(1));
        if length(gmid) && now - pen.datenum < DATA.plot.maxage
            gmid = gmid(end);
            if strncmpi(pen.comments{id(gmid)},'Grey',4)
            np = np+1;
            penx(np) = pen.pos(1);
            peny(np) = pen.pos(2);
            penz(np) = pen.depths(pen.cmtime(id(gmid))) - pen.enterdepth;
            pid(np) = j;
            else
                np = np+1;
                penx(np) = pen.pos(1);
                peny(np) = pen.pos(2);
                penz(np) = pen.depths(pen.cmtime(cid(1))) - pen.enterdepth;
                pid(np) = j;
            end
        end
    end
    end
end
if smooth > 0
[X,Y] = meshgrid(linspace(min(penx),max(penx)),linspace(min(peny),max(peny)));
Z = Interpf(penx,peny,penz,X,Y,1,smooth);
pcolor(X,Y,Z);
shading('interp');
colorbar;
else
    for j = 1:length(penz)
        plot3(penx(j),peny(j),penz(j),'o','buttondownfcn',['PlotMap(' num2str(DATA.top) ',''penpoint'',' num2str(pid(j)) ',1);']);
        hold on;
    end
end




function PlotGridArrows(DATA, varargin)

MapDefs;


plotmean = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'mean',4)
        plotmean = 1;
    end
    j = j+1;
end

colors = {[0 0 0], [0.3 0 0], [0.6 0.3 0] [0.6 0 0], [1 1 0 ], [1 0 0]} ;
top = num2str(DATA.top);
scale = 10;
id = find(DATA.selected);
for k = 1:length(id)
    j = id(k);
    p = DATA.map.pen(j);
    px = DATA.map.pen(j,2);
    py = DATA.map.pen(j,3);
    xs(:,k) = [px px+DATA.map.rf(j,1)./scale];
    ys(:,k) = [py py+DATA.map.rf(j,2)./scale];
    arrow([px px+DATA.map.rf(j,1)./scale],[py py+DATA.map.rf(j,2)./scale],20,0.1);
    hold on;
    pe = find(DATA.id == p);
    if length(pe) == 1
    ahs(k) = plot(px,py,'o', 'markerfacecolor','b','buttondownfcn',['PlotMap(' top ',''gridpoint'',' num2str(pe) ',1);']);
    end
    hs(k) = plot(px+DATA.map.rf(j,1)./scale,py+DATA.map.rf(j,2)./scale,'o','buttondownfcn',['PlotMap(' top ',''cellpoint'',' num2str(j) ',1);']);
    text(0.1+px+DATA.map.rf(j,1)./scale,py+DATA.map.rf(j,2)./scale,DATA.map.cellname{j});
end

if plotmean
    hold off;
    pens = xs(1,:) + i* ys(1,:);
    up = unique(pens);
    for j = 1:length(up)
        id = find(pens == up(j));
        xm = mean(xs(2,id));
        ym = mean(ys(2,id));
        px = real(up(j));
        py = imag(up(j));
        ahs(k) = plot(px,py,'o', 'markerfacecolor','b','buttondownfcn',['PlotMap(' top ',''gridpoint'',' num2str(pe) ',1);']);
        arrow([px xm],[py ym],20,0.1);
        hold on;
    end
end

minx = floor(min(xs(:)))-1;
maxx = ceil(max(xs(:)))+1;
miny = floor(min(ys(:)))-1;
maxy = ceil(max(ys(:)))+1;
set(gca,'xlim',[minx maxx]);
set(gca,'ylim',[miny maxy]);
[sorted, idx] = sort(DATA.pid);
colors = {[0 0 0], [0.3 0 0], [0.6 0.3 0] [0.6 0 0], [1 1 0 ], [1 0 0]} ;

if length(idx) > 5
    back = 5;
else
    back = length(idx)-1;
end
k = 1;
hs = [];
for j = idx(end-back:end)
    hs(k) = plot(DATA.ppx(j),DATA.ppy(j),'o','buttondownfcn',['PlotMap(' top ',''gridpoint'',' num2str(DATA.ppid(j)) ',1);'],'color',colors{k},'MarkerFaceColor',colors{k});
    labels{k} = sprintf('%d',DATA.pid(j));
    k = k+1;
end
legend(hs,labels);
if ~isempty(DATA.missed)
id = find(DATA.missed == MT);
for j =1:length(id)
    k = find(DATA.map.pen(:,1) == id(j));
    h = plot(DATA.map.pen(k,2),DATA.map.pen(k,3),'co','buttondownfcn',['PlotMap(' top ',''missedpoint'',' num2str(id(j)) ',1);'],'color','c','MarkerFaceColor','c');
end
id = find(DATA.missed == STS);
for j =1:length(id)
    k = find(DATA.map.pen(:,1) == id(j));
    h = plot(DATA.map.pen(k,2),DATA.map.pen(k,3),'co','buttondownfcn',['PlotMap(' top ',''missedpoint'',' num2str(id(j)) ',1);'],'color','g','MarkerFaceColor','g');
end
end


xr = get(gca,'Xlim');
yr = get(gca,'Ylim');
set(gca,'Xtick',[xr(1):xr(2)],'YTick',[yr(1):yr(2)],'Xgrid','on','Ygrid','on');


function PlotGridVar(DATA)

for j = 1:length(DATA.xc)
    if DATA.plot.labelpts
        text(DATA.px(j),DATA.py(j),sprintf('%.1f,%.1f',DATA.xc(j),DATA.xc(j)),'fontsiz',DATA.plot.fontsiz,'HorizontalAlignment','center','VerticalAlign','bottom');
    end
    h = plot3(DATA.px(j),DATA.py(j),DATA.rvar(j),'o','buttondownfcn',['PlotMap(DATA.top,''gridpoint'',' num2str(j) ',1);'],'MarkerFaceColor','b');
    if DATA.penmarked(j)
        set(h,'MarkerFaceColor','r');
    end
    hold on;
end
xr = get(gca,'Xlim');
yr = get(gca,'Ylim');
set(gca,'Xtick',[xr(1):xr(2)],'YTick',[yr(1):yr(2)],'Xgrid','on','Ygrid','on');

function PlotGridDepth(DATA)

hold off;
idx = find(DATA.selected > 0);
selected = idx;
if isempty(idx)
    fprintf('No data for Area %d\n',area);
    return;
end

for j = idx;
    if(DATA.map.depth(j) > -2000)
    if DATA.plot.labelpts
        text(DATA.px(j),DATA.py(j),sprintf('%.1f,%.1f',DATA.xc(j),DATA.xc(j)),'fontsiz',DATA.plot.fontsiz,'HorizontalAlignment','center','VerticalAlign','bottom');
    end
    h = plot3(DATA.map.rf(j,1),DATA.map.rf(j,2),DATA.map.depth(j),'o','buttondownfcn',['PlotMap(DATA.top,''depthpoint'',' num2str(j) ',1);'],'MarkerFaceColor','b');
    if DATA.marked(j)
        set(h,'MarkerFaceColor','r');
    end
    hold on;
    end
end
xr = get(gca,'Xlim');
yr = get(gca,'Ylim');
set(gca,'Xtick',[xr(1):xr(2)],'YTick',[yr(1):yr(2)],'Xgrid','on','Ygrid','on');

function ChoosePens(DATA)
  scrsz = get(0,'Screensize');
  wsc = scrsz(3)/1000;
  SPACE = 5;   VSPACE = 2;   BOXH = 15;
  
  selid = find(DATA.selected > 0)
  len = length(selid);
  if (len < 20)
      h = (len+1) * (15 + VSPACE);
  else
      h = (20+1) * (15 + VSPACE);
  end

  CloseTag('Blocks');
  cntrl_box = figure('Position', [200 scrsz(4)-(h+50) 500 h], 'Menubar', 'none',...
       'NumberTitle', 'off', 'Tag','Blocks','Name','BlockList');
  
   bp = [30 0 500 15];
   for k = 1:length(selid);
       j = selid(k);
      desc = sprintf('%d: %d,%d %s %.1f,%.1f',DATA.map.pen(j,1),DATA.map.pen(j,2),DATA.map.pen(j,3),...
          DATA.map.cellname{j},DATA.map.rf(j,2),DATA.map.rf(j,3));
      val = ~(DATA.excluded(j));
      uicontrol(gcf,'Style', 'CheckBox','String',desc,'Position', bp * wsc,...
      'Tag',sprintf('Block%d',k),'value',val,'UserData',j);
      bp(2) = bp(2) + bp(4) + VSPACE;
  end
  bp(3) = 25;
  bp(1) = 1;
  bp(2) = h-100;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['PlotMap(' top ',''checkblocks'')'],...
'String', 'Go', 'Position', bp * wsc);
  bp(2) = bp(2)-bp(4);

  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['PlotMap(' top ',''clearblocks'')'],...
'String', 'clr', 'Position', bp * wsc);
  bp(2) = bp(2)-bp(4);
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['PlotMap(' top ',''close'')'],...
'String', 'end', 'Position', bp * wsc);

  uicontrol(gcf, 'callback', ['PlotMap(' top ',''scroll'')'],'style','slider','min',1','max',k,'value',1,...
      'position',[480 5 20 h-10],'Tag','blockscroll');
  

function scrolllist(listpos)

BOXH = 15;
  VSPACE = 2;

blocks = GetBlockList('Block');
for j = 1:length(blocks)
    pos = get(blocks(j),'position');
    pos(2) = (j - (listpos)) * (BOXH + VSPACE);
    set(blocks(j),'position',pos);
end


function [blocks, isset, data] = GetBlockList(prefix, varargin)

j = 1;
it = 1;
while ~isempty(it)
it = findobj('Tag',sprintf('%s%d',prefix,j));
if ~isempty(it)
    blocks(j) = it;
    isset(j) = get(it,'value');
    data(j) = get(it,'UserData');
    j = j+1;
end
end

function ShowCell(map, id)
fprintf('%s %d %s %.1f %.1f %.1f\n',map.cellname{id},map.pen(id,1),map.datestr{id},map.rf(id,1),map.rf(id,2),map.depth(id));
  

function PrintPens(DATA,map)

[idx, order] = sort(map.datenum)
for j = 1:length(order)
    if DATA.selected(order(j))
    fprintf('%s P%d %.1f %.1f %s\n',map.datestr{order(j)},map.pen(order(j),1),map.pen(order(j),2),map.pen(order(j),3),map.cellname{order(j)});
    end
end

function SaveMap(map, monkey, varargin)

mapfile = sprintf('/bgc/bgc/anal/%s/MapData.mat',monkey);
save(mapfile,'map');

function map = ReadAllPens(map, monkey, varargin)

pes = unique(map.pen(:,1));
reload = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'reload',4)
        reload = 1;
    end
    j = j+1;
end

penfile = sprintf('/bgc/bgc/anal/%s/pens/pendata.mat',monkey);
if exist(penfile,'file') && ~reload
    load(penfile);
    map.pens = pens;
    map.missed = missed;
    np = length(pens);
    for j = 1:length(pens)
        if isfield(pens{j},'num')
            oldpes(j) = pens{j}.num;
        end
    end
    id = find(oldpes < 0);
    oldpes(id) = 9000-oldpes(id); 
    id = find(pes < 0);
    pes(id) = 9000-pes(id); 
    newpes = setdiff(pes,oldpes);
else
    newpes = pes;
    np = 0;
end
for j = 1:length(newpes)
    pe = newpes(j);
    if pe >= 1
        name = sprintf('/bgc/anal/%s/pens/pen%d.log',monkey,pe);
    else
        name = sprintf('/bgc/anal/%s/pens/pen%d.log',monkey,9000-pe);
        pe = 9000-pe;
    end
    map.pens{pe} = ReadPen(name,'noplot');
    if ~isfield(map.pens{pe},'num')
        map.pens{pe}.num = pe;
    end
    if isfield(map.pens{pe},'missed')
        map.missed(pe) = map.pens{pe}.missed;
    else
        map.missed(pe) = 0;
    end
end
if length(newpes)
    pens = map.pens;
    missed = map.missed;
    save(penfile,'pens','missed');
end
