function binoc(varargin)
%binoc
%GUI for running binoclean via pipe.

if length(varargin) && ishandle(varargin{1})
    f = varargin{1};
    while ~isfigure(f)
        f = get(f,'parent');
    end
    DATA = get(f,'UserData');
    varargin = varargin(3:end);
else
TOPTAG = 'binoc';
it = findobj('Tag',TOPTAG,'type','figure');
if isempty(it)
    DATA.stimfile = 'test.in';
    DATA.name = 'Binoc';
    DATA.tag.top = 'Binoc';
    DATA = SetDefaults(DATA);
    DATA = InitInterface(DATA);
end
end
if length(varargin)
    if strncmpi(varargin{1},'close',5)
        f = fields(DATA.tag);
        for j = 1:length(f)
            CloseTag(DATA.tag.(f{j}));
        end
        return;
    end
end

function DATA = InterpretLine(DATA, s)

if s(1) == '#' %defines stim code/label
    [a,b] = sscanf(s,'#%d %s');
    a = a(1);
    id = find(s == ' ');
    DATA.comcodes(a).code = s(id(1)+1:id(2)-1);
    DATA.comcodes(a).label = s(id(2)+1:end);
    DATA.comcodes(a).const = a;
elseif strncmp(s,'Expts1',6)
    DATA.extypes{1} = sscanf(s(8:end),'%d');
elseif strncmp(s,'Expts2',6)
    DATA.extypes{2} = sscanf(s(8:end),'%d');
elseif strncmp(s,'Expts3',6)
    DATA.extypes{3} = sscanf(s(8:end),'%d');
elseif strncmp(s,'qe=',3)
    s = s(3:end);
    [a,b] = fileparts(s);
    n = length(DATA.quickexpts)+1;
    DATA.quickexpts(n).name = b;
    DATA.quickexpts(n).filename = s;
end

function DATA = SetDefaults(DATA)

DATA.outpipe = '/dev/binocin';
DATA.inpipe = '/dev/binocout';
DATA.outid = fopen(DATA.outpipe,'w');
DATA.inid = fopen(DATA.inpipe,'r');
DATA.incr = [0 0 0];
DATA.nstim = [0 0 0];
DATA.quickexpts = [];
DATA.stepsize = [20 10];
DATA.stepperpos = -2000;
DATA.tag.stepper = 'Stepper';
fid = fopen(DATA.stimfile,'r');
tline = fgetl(fid);
while ischar(tline)
    DATA = InterpretLine(DATA,tline);
    tline = fgetl(fid);
end
fclose(fid);


DATA.expstrs{1} = {};
DATA.expstrs{2} = {};
DATA.expstrs{3} = {};
for j = 1:length(DATA.comcodes)
    if ismember(DATA.comcodes(j).const,DATA.extypes{1})
        DATA.expstrs{1} = {DATA.expstrs{1}{:} DATA.comcodes(j).label};
    end
    if ismember(DATA.comcodes(j).const,DATA.extypes{2})
        DATA.expstrs{2} = {DATA.expstrs{2}{:} DATA.comcodes(j).label};
    end
    if ismember(DATA.comcodes(j).const,DATA.extypes{3})
        DATA.expstrs{3} = {DATA.expstrs{3}{:} DATA.comcodes(j).label};
    end
end


function DATA = InitInterface(DATA)

    scrsz = get(0,'Screensize');
    cntrl_box = figure('Position', [10 scrsz(4)-480 300 450],...
        'NumberTitle', 'off', 'Tag',DATA.tag.top,'Name',DATA.name,'menubar','none');
    DATA.toplevel = cntrl_box;
    lst = uicontrol(gcf, 'Style','edit','String', '',...
        'HorizontalAlignment','left',...
        'Callback', {@TextEntered}, 'Tag','NextButton',...
'units','norm', 'Position',[0.01 0 0.98 0.1]);
    DATA.txtui = lst;
    nr = 25;
     lst = uicontrol(gcf, 'Style','list','String', 'Commmand History',...
        'HorizontalAlignment','left',...
        'Max',10,'Min',0,...
        'Callback', {@TextEntered}, 'Tag','NextButton',...
'units','norm', 'Position',[0.01 1./nr 0.98 8/nr]);
   DATA.txtrec = lst;
    nc = 4;

    bp(2) = 9./nr;

    bp(1) = 0.12;
    bp(4) = 4/nr;
    bp(3) = 0.3;
    uicontrol(gcf,'style','list','string',num2str(DATA.nstim(1)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1List');
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','list','string',num2str(DATA.nstim(2)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt2List');
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','list','string',num2str(DATA.nstim(3)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt3List');
    
    
    bp(1) = 0.01;
    bp(2) = bp(2)+bp(4);
    bp(3) = 0.1;
    bp(4) = 1/nr;
    uicontrol(gcf,'style','text','string','N stim',  'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = 0.2;
    uicontrol(gcf,'style','edit','string',num2str(DATA.nstim(1)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1Nstim');
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.nstim(2)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt2Nstim');
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.nstim(3)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt3Nstim');

    bp(1) = 0.01;
    bp(2) = bp(2)+bp(4);
    bp(3) = 0.1;
    uicontrol(gcf,'style','text','string','Incr',  'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = 0.3;
    uicontrol(gcf,'style','edit','string',num2str(DATA.incr(1)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1Incr');
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.incr(2)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt2Incr');
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.incr(3)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt3Incr');

    bp(1) = 0.01;
    bp(2) = bp(2)+bp(4);
    bp(3) = 0.1;
    
    uicontrol(gcf,'style','pushbutton','string','Run', ...
        'Callback', {@NextButton, 'l'}, 'Tag','NextButton',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = 0.3;
    uicontrol(gcf,'style','pop','string',DATA.expstrs{1}, ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1List');
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pop','string',DATA.expstrs{3}, ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1List');
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pop','string',DATA.expstrs{3}, ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1List');
     
    bp(1) = 0.01;
    bp(2) = 0.99-1/nr;
    bp(3) = 1/nc;
    bp(4) = 1./nr;
    uicontrol(gcf,'style','checkbox','string','go', ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1List');
    bp(1) = bp(1)+bp(3);
    uicontrol(gcf,'style','checkbox','string','store', ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1List');
    
    
    hm = uimenu(cntrl_box,'Label','File','Tag','BinocFileMenu');
    uimenu(hm,'Label','Close','Callback',{@binoc, 'close'});
    hm = uimenu(cntrl_box,'Label','Quick','Tag','QuickMenu');
    for j = 1:length(DATA.quickexpts)
    uimenu(hm,'Label',DATA.quickexpts(j).name,'Callback',{@binoc, 'quick', DATA.quickexpts(j).filename});
    end
    hm = uimenu(cntrl_box,'Label','Pop','Tag','QuickMenu');
    uimenu(hm,'Label','Stepper','Callback',{@StepperPopup});

 set(DATA.toplevel,'UserData',DATA);
 
function StepperPopup(a,b);
  DATA = GetDataFromFig(a);
  cntrl_box = findobj('Tag',DATA.tag.stepper,'type','figure');
if ~isempty(cntrl_box)
    figure(cntrl_box);
    return;
end
scrsz = get(0,'Screensize');
cntrl_box = figure('Position', [10 scrsz(4)-480 300 450],...
        'NumberTitle', 'off', 'Tag',DATA.tag.stepper,'Name','Stepper','menubar','none');
    nr = 6;
    bp = [0.01 0.99-1/nr 0.1 1./nr]
    uicontrol(gcf,'style','pushbutton','string','+', ...
        'Callback', {@Stepper, 1, 1},...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pushbutton','string','-', ...
        'Callback', {@Stepper, -1, 1},...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = 0.3;
    uicontrol(gcf,'style','pop','string','10|20|50|100|200', ...
        'units', 'norm', 'position',bp,'value',1,'Tag','StepSize1');
    bp(1) = 0.01;
    bp(2) = bp(2) - bp(4)-0.01;
    uicontrol(gcf,'style','pushbutton','string','+', ...
        'Callback', {@Stepper, 1, 2},...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pushbutton','string','-', ...
        'Callback', {@Stepper, -1, 2},...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = 0.3;
    uicontrol(gcf,'style','edit','string','10', ...
        'units', 'norm', 'position',bp,'value',1,'Tag','StepSize2');
    bp(1) = 0.01;
    bp(2) = bp(2) - bp(4)-0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.stepperpos), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','StepperPosition');
    set(cntrl_box,'UserData',DATA.toplevel);
    
function Stepper(a,b, step, type)
    DATA = GetDataFromFig(a);
    sfig = get(a,'parent');
    it = findobj(sfig,'Tag','StepSize1');
    DATA.stepsize(1) = Menu2Val(it);
    it = findobj(sfig,'Tag','StepSize2');
    DATA.stepsize(2) = Text2Val(it);
    if step > 0
        s = sprintf('ed+%.3f\n',DATA.stepsize(type));
    else
        s = sprintf('ed-%.3f\n',DATA.stepsize(type));
    end
    fprintf('%s\n',s);
function val = Menu2Val(it)
val = NaN;
if isempty(it)
    return;
end
j = get(it(1),'value');
s = get(it(1),'string');
val = str2num(s(j,:));

function val = Text2Val(it)
val = NaN;
if isempty(it)
    return;
end
s = get(it,'string');
val = str2num(s);

        
        
        
 
function TextEntered(a,b)
    DATA = GetDataFromFig(a);
txt = get(a,'string');
fprintf('%s',txt);
set(a,'string','');

a =  get(DATA.txtrec,'string');
n = size(a,1);
a(n+1,1:length(txt)) = txt;
set(DATA.txtrec,'string',a);
set(DATA.txtrec,'listboxtop',n+1);
