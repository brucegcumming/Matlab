function SubplotRedraw(np,npts)
% SubplotRedraw(np,npts)  create np x np subplots, and 
% drawy npts points in each, also draw a line in each.
% then just change the color of one line, and see how long redraw takes.
%adds a simple button drag functino that moves the line end when mouse is
%dragged. with np = 10 and npts = 100, there is a very noticeable delay
if nargin < 2
    npts = 1000;
end
optimize = 2;

ng = 0;
w = 0.9./np;
for j = 1:np;
    for k = 1:np;
        rnd = rand(1000,2);
        ng = ng+1;
        ax(ng) = subplot('Position',[(j-1)/np (k-1)/np w w]);
        hold off;
        plot(rnd(:,1),rnd(:,2),'.');
        hold on;
        set(gca,'xticklabel',[],'yticklabel',[],'xlim',[0 1],'ylim',[0 1]);
        set(gca,'xlimmode','manual','ylimmode','manual');
        h(ng) = plot([0 1],[0.5 0.5],'r-');
        set(gca,'userdata',h(ng));
        if optimize
            set(gca,'zlimmode','manual','climmode','manual');
            set(gca,'xtickmode','manual','ytickmode','manual','ztickmode','manual');
            set(gca,'xticklabelmode','manual','yticklabelmode','manual','zticklabelmode','manual');
        end
    end
end
if optimize > 1
set(gcf,'renderer','opengl');
end
set(gcf, 'WindowButtonMotionFcn', @ButtonDragged)
set(gcf, 'WindowButtonDownFcn',{@ButtonPressed, 1});
set(gcf, 'WindowButtonUpFcn',{@ButtonPressed, 0});
setappdata(gcf,'mouse_is_down', 0);
tic
drawnow;
toc
tic;
set(h(ng),'color','g');
drawnow;
toc

function ButtonPressed(src, data, down)
%just record if mouse is down/up, so Dragged works while button held down
setappdata(gcf,'mouse_is_down', down);

function ButtonDragged(src, data)
%only draw line if mouse is dragged with button down
down = getappdata(gcf,'mouse_is_down');
if down == 0
    return;
end
    h = get(gca,'userdata');
    x = get(h,'xdata');
    y = get(h,'ydata');
    start = get(gca,'CurrentPoint');
    x(2) = start(1,1);
    y(2) = start(1,2);
    set(h,'xdata',x,'ydata',y);
    tic
    drawnow;
    toc
