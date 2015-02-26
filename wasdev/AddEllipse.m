function varargout = AddEllipse(F,varargin)
%adds an ellipse to a figure interactively
%AddEllipse(fig) returns immediately
%once drawn, adds appdata 'InteractiveEliipse' to the figure
%AddEllipse(fig, 'wait') waits for the user to draw the ellipse, then
%            returns a structure defining the ellipse. Time out after 10s
%AddEllipse(fig, 'wait', x) sets the timeout to x. If the timeout is NaN,
%it waits forever
%  


getresp = 0;
timeout = 10;
shape = 0;
plotargs = {};
j = 1;
appstr = 'InteractiveEllipse';
mousept = getappdata(F, appstr);
while j <= length(varargin)
    if strncmpi(varargin{j},'get',3)
        varargout{1} = mousept;
        return;
    elseif strncmpi(varargin{j},'wait',4)
        getresp = 1;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            timeout = varargin{j};
        end
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


mousept.starttime = now;
mousept.showtimes = 0;
if ~isfield(mousept,'h')
    mousept.h= -1;
end
mousept.shape= shape;
mousept.down = 0;
mousept.done = 0;
mousept.plotargs = plotargs;
mousept.dragfcn = get(F,'WindowButtonMotionFcn');
%should get old Fcns, then reset them after button release
set(F, 'WindowButtonDownFcn',@ButtonPressed);
set(F, 'WindowButtonMotionFcn',@ButtonDragged);
set(F, 'WindowButtonUpFcn',@ButtonReleased);
set(F, 'WindowScrollWheelFcn',@ScrollWheel);
setappdata(F,appstr, mousept);
if mousept.showtimes
            mytoc(mousept.starttime);
end
figure(F);
if getresp
    waited = 0;
    while mousept.done == 0
        mousept = getappdata(F,appstr);
        pause(0.05);
        waited = waited + 0.05;
        if waited > timeout
            varargout{1} = [];
            return;
        end
    end
    varargout{1} = mousept;
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
appstr = 'InteractiveEllipse';
mousept = getappdata(src,appstr);

if isfield(mousept,'h') && ishandle(mousept.h)
    delete(mousept.h);
    mousept.h = -1;
end


%if mousept.showtimes
%    mytoc(mousept.starttime);
%end
start = get(gca,'CurrentPoint');
mousept.down = 1;
mousept.pos =[start(1,1) start(1,2) 0 0 ];
mousept.axis = gca;
setappdata(src,appstr,mousept);
%if mousept.showtimes
%    mytoc(mousept.starttime);
%end

function ButtonReleased(src, data)
appstr = 'InteractiveEllipse';
mousept = getappdata(src,appstr);

start = get(gca,'CurrentPoint');
mousept.down = 0;
mousept.done = 1;
p = mousept.pos;
mousept.xyr = [mean(p([1 3])) mean(p([2 4])) diff(p([1 3]))/2 diff(p([2 4]))/2]; 
if isfield(mousept,'dragfcn')
set(src, 'WindowButtonMotionFcn',mousept.dragfcn);
end
setappdata(src,appstr, mousept);

function ButtonDragged(src, data)
appstr = 'InteractiveEllipse';
mousept = getappdata(src,appstr);

if isfield(mousept,'down') && mousept.down > 0
    start = get(gca,'CurrentPoint');
    mousept.pos(3) = start(1,1);
    mousept.pos(4) = start(1,2);
    mousept.h = DrawEllipse(mousept, mousept.plotargs{:});
end
setappdata(src,appstr, mousept);

