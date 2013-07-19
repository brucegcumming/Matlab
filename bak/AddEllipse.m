function varargout = AddEllipse(F,varargin)
DATA = GetDataFromFig(F);
getresp = 0;
timeout = 10;
shape = 0;
plotargs = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'get',3)
        varargout{1} = DATA.elmousept;
        return;
    elseif strncmpi(varargin{j},'wait',4)
        getresp = 1;
    elseif strncmpi(varargin{j},'line',4)
        shape = 1;
    elseif strncmpi(varargin{j},'timeout',6)
        j = j+1;
        timeout = varargin{j};
    else
        plotargs = {plotargs{:} varargin{j}};
    end
    j = j+1;
end

DATA.elmousept.h= -1;
DATA.elmousept.shape= shape;
DATA.elmousept.down = 0;
DATA.elmousept.done = 0;
DATA.elmousept.plotargs = plotargs;
DATA.elmousept.dragfcn = get(F,'WindowButtonMotionFcn');
%should get old Fcns, then reset them after button release
set(F, 'WindowButtonDownFcn',@ButtonPressed);
set(F, 'WindowButtonMotionFcn',@ButtonDragged);
set(F, 'WindowButtonUpFcn',@ButtonReleased);
set(F, 'WindowScrollWheelFcn',@ScrollWheel);
set(F,'UserData',DATA);

if getresp
    waited = 0;
    while DATA.elmousept.done == 0
        DATA = GetDataFromFig(F);
        pause(0.05);
        waited = waited + 0.05;
        if waited > timeout
            varargout{1} = [];
            return;
        end
    end
    varargout{1} = DATA.elmousept;
end

function h= DrawEllipse(E,varargin)

if E.shape == 1
    h = DrawLine(E,varargin{:});
    return;
end
a = (E.pos(3)-E.pos(1))/2; %x radius
b = (E.pos(4)-E.pos(2))/2;
sn = 0;
cn = 1;
x = linspace(0,a);
y =  sqrt(b.^2 - (x.*b/a).^2);
x = [x fliplr(x) -x fliplr(-x)]+mean(E.pos([1 3]));
y = [y fliplr(-y) -y fliplr(y)]+mean(E.pos([2 4]));
if ishandle(E.h) 
    set(E.h,'Xdata',x,'Ydata',y);
    h = E.h;
else
    hold on;
    h = plot(real(x),real(y),varargin{:});
    hold off;
end

function h= DrawLine(E,varargin)

x = [E.pos(1) E.pos(3)];
y = [E.pos(2) E.pos(4)];
if ishandle(E.h) 
    set(E.h,'Xdata',x,'Ydata',y);
    h = E.h;
else
    hold on;
    h = plot(real(x),real(y),varargin{:});
    hold off;
end



function ButtonPressed(src, data)
DATA = GetDataFromFig(src);

start = get(gca,'CurrentPoint');
DATA.elmousept.down = 1;
DATA.elmousept.pos =[start(1,1) start(1,2) 0 0 ];
DATA.elmousept.axis = gca;
set(src,'UserData',DATA);

function ButtonReleased(src, data)
DATA = GetDataFromFig(src);
start = get(gca,'CurrentPoint');
DATA.elmousept.down = 0;
DATA.elmousept.done = 1;
p = DATA.elmousept.pos;
DATA.elmousept.xyr = [mean(p([1 3])) mean(p([2 4])) diff(p([1 3]))/2 diff(p([2 4]))/2]; 
if isfield(DATA.elmousept,'dragfcn')
set(src, 'WindowButtonMotionFcn',DATA.elmousept.dragfcn);
end
set(src,'UserData',DATA);

function ButtonDragged(src, data)
DATA = GetDataFromFig(src);

if isfield(DATA,'elmousept') && DATA.elmousept.down > 0
    start = get(gca,'CurrentPoint');
    DATA.elmousept.pos(3) = start(1,1);
    DATA.elmousept.pos(4) = start(1,2);
    DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
end
set(src,'UserData',DATA);

