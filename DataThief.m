function DataThief(im, varargin);

%DataThief(im,....)
%use mouse to get data points from an image. 
%
%DataThief(im,[xl yl xu yu])
%specifies axis range. Then User clicks on extremes, and range for mouse
%clicks is set automatially
%get UserData from the figure to use data
DATA = [];
DATA.x = [];
DATA.y = [];
DATA.axptsx = [];
DATA.axptsy = [];
DATA.scaled = 0;
DATA.state.deleting = 0;
DATA.state.showduplicate = 1;
DATA.title = 'Data Image';

j = 1;
if isnumeric(varargin{1}) & length(varargin{1}) == 4
    ranges = varargin{1};
    j = 2;
end
while j <= length(varargin)
    if isstruct(varargin{j}) && isfield(varargin{j},'axptsx')
        DATA = varargin{j};
        DATA.scaled = 0;
    elseif strncmpi(varargin{j},'nodup',5)
        DATA.state.showduplicate = 0;
    end
    j = j+1;
end

GetFigure(DATA.title);
DATA.fig = gcf;
if isempty(findobj(gcf,'Tag','DataMenu'))
hm = uimenu(gcf,'Label','DataThief','Tag','DataMenu');
uimenu(hm,'Label','Redo axes','Callback',{@SetAxes});
uimenu(hm,'Label','Delete one point','Callback',{@DeletePoint});
uimenu(hm,'Label','Delete mulitple points','Callback',{@DeletePoint});
uimenu(hm,'Label','Make Histogram','Callback',{@MakeHist});
end

hold off;
DATA.image = flipdim(im,1);
if isfield(DATA,'ranges') && length(DATA.ranges) == 4
    ranges = DATA.ranges;
end
if length(ranges) == 4
imagesc([ranges(1) ranges(3)],[ranges(2) ranges(4)], DATA.image);
DATA.ranges = ranges;
else
imagesc(DATA.image);
end
set(gca,'ydir','normal');
if length(DATA.axptsx) == 4
    set(gcf,'UserData',DATA);
    set(gca,'xlim',DATA.imscale([1 2]),'ylim',DATA.imscale([ 3 4]));
    ButtonPressed(gcf,0);
else
title('Click on axes limits');
end
set(gcf, 'WindowButtonDownFcn',@ButtonPressed);
set(gcf,'UserData',DATA);
set(gcf, 'Pointer','crosshair');

function ButtonPressed(src, data)

DATA = get(gcf,'UserData');
pt = get(gca,'CurrentPoint');
if (length(DATA.axptsx) < 4 || DATA.scaled == 0) & isfield(DATA,'ranges')
    ranges = DATA.ranges;
DATA.axptsx = [DATA.axptsx pt(1,1)];
DATA.axptsy = [DATA.axptsy pt(1,2)];
hold on;
plot(pt(1,1),pt(1,2),'g+');
   if length(DATA.axptsx) >= 4
       if isfield(DATA,'imrange')
           imrange = DATA.imrange;
       else
           DATA.imscale = [get(gca,'xlim') get(gca,'ylim')];
           imrange = ranges([1 3 2 4]);
       end
           
       xsc = (ranges(3)-ranges(1))./range(DATA.axptsx);
       newxrange = (imrange(2)-imrange(1)) .* xsc;
       xoff = (ranges(1) - min(DATA.axptsx)) * xsc;
       ysc = (ranges(4)-ranges(2))./range(DATA.axptsy);
       newyrange = (imrange(4)-imrange(3))*ysc;
       yoff = (ranges(2) - min(DATA.axptsy)) * ysc;
       hold off;
       DATA.imrange = [imrange(1)+xoff imrange(1)+newxrange+xoff imrange(3)+yoff imrange(3)+newyrange+yoff];
       h= imagesc(DATA.imrange([1 2]),DATA.imrange([3 4]),DATA.image);
       dy = diff(ranges([2 4]))/20;
       dx = diff(ranges([1 3]))/20;
       set(gca,'ydir','normal','ylim',[ranges(2)-dy ranges(4)+dy],'xlim',[ranges(1)-dx ranges(3)+dx]);
       hold on;
       plot(ranges([1 1 3 3 1]),ranges([2 4 4 2 2]),'g');
%       set(gca,'xlim',[ranges(1) ranges(1)+newxrange],'ylim',[ranges(2) ranges(2)+newyrange]);
       title('Now set data pts');
       DATA.scaled = 1;
   end
elseif DATA.state.deleting
    d = abs((DATA.x+i*DATA.y)-(pt(1,1) + i * pt(1,2)));
    [a,b] = min(d);
    id = setdiff(1:length(DATA.x),b);
    DATA.x = DATA.x(id);
    DATA.y = DATA.y(id);
    delete(DATA.h(b));
    DATA.h = DATA.h(id);
 %state 2 = delete just deletes on point, so put state back. 
    if DATA.state.deleting == 2
    DATA.state.deleting = 0;
    end
else
DATA.x = [DATA.x pt(1,1)];
DATA.y = [DATA.y pt(1,2)];
hold on;
DATA.h(length(DATA.x)) = plot(pt(1,1),pt(1,2),'rx');
end
set(gcf,'UserData',DATA);
if length(DATA.x) & DATA.state.showduplicate
GetFigure('Duplicate');
plot(DATA.x,DATA.y,'o');
end

function MakeHist(a,b)
DATA = get(gcf,'UserData');
of = gcf;
GetFigure('Duplicate');
hold off;
y = round(DATA.y);
w = mean(diff(DATA.x));
x = min(DATA.x):mean(diff(DATA.x)):max(DATA.x);
bar(x,y,1);
set(gca,'xlim',DATA.ranges([1 3]),'ylim',DATA.ranges([2 4]));
pts = [];
py = [];
for j = 1:length(x)
    if y(j) > 0
    pts = [pts [x(j)-w/2+[1:y(j)] * w/y(j)]-w/(2*y(j))];
    py = [py ones(1,y(j)) .* y(j)]; 
    end
end
hist(pts,x);
hold on;
plot(pts,py,'ro');
set(gca,'xlim',DATA.ranges([1 3]),'ylim',DATA.ranges([2 4]));
DATA.histx = pts;
set(of,'UserData',DATA);

function SetAxes(a,b)
DATA = get(gcf,'UserData');
DATA.axptsx = [];
DATA.axptsy = [];
title('Click on axes limits to set scale');
set(gcf,'UserData',DATA);

function DeletePoint(a,b)
DATA = get(gcf,'UserData');
la = get(a,'label');

if strncmp('Delete one',la,10)
DATA.state.deleting = 2;
elseif strncmp('Delete multiple',la,10)
DATA.state.deleting = 1;
set(a,'Label','Stop Deleting');
elseif strncmp('Stop Deleting',la,10)
DATA.state.deleting = 0;
set(a,'Label','Delete multiple points');
end
set(gcf,'UserData',DATA);
