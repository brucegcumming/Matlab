function TestEllipse(ratio, varargin)

j = 1;
x = rand(100,2) - 0.5;
freeangle = 1;
x(:,2) = x(:,2).*ratio;
ne = 0;
guesspt = [NaN NaN NaN]
while j <= length(varargin)
    if strncmpi(varargin{j},'angle',4);
    j = j+1;
        ne = ne+1;
        angles = varargin{j};
    elseif strncmpi(varargin{j},'fixangle',4);
        freeangle = 0;
    elseif strncmpi(varargin{j},'Guess',4);
        j = j+1;
        guesspt = varargin{j};
    end
    j = j+1;
end

F = GetFigure('TestEllipse');
hold off;
plot(x(:,1),x(:,2),'.');
DATA.pts = x;
DATA.elmousept.h= -1;
shapes = [0 1 0];
        DATA.elmousept.shape= 0;
        DATA.elmousept.down = 0;
        DATA.elmousept.done = 0;
        DATA.elmousept.steps = 0;
        DATA.elmousept.angle = 0;
        DATA.elmousept.cluster = 1;
        DATA.elmousept.color = [1 0 0];
        DATA.elmousept.plotargs = {};
%        DATA.elmousept.dragfcn = get(F,'WindowButtonMotionFcn');
        %should get old Fcns, then reset them after button release
        set(F, 'WindowButtonDownFcn',@XYButtonPressed);
        set(F, 'WindowButtonMotionFcn',@XYButtonDragged);
        set(F, 'WindowButtonUpFcn',@XYButtonReleased);
        set(F, 'WindowScrollWheelFcn',@ScrollWheel);
        r(1) = 0.5;
        r(2) = ratio/4;
        ne = length(angles);
        for j = 1:ne
            DATA.elmousept.xyr = [0 0  r(1) r(2)];
            DATA.elmousept.angle = angles(j);
            [a,details] = DrawEllipse(DATA.elmousept,'color','g','linewidth',2);
            hold on;
            set(gca,'xlimmode','manual','ylimmode','manual');
            if 0
                DATA.elmousept.angle = atan(tan(angles(j))./ratio);
                rfix = tan(angles(j)).*sqrt(ratio);
                a = angles(j);
                a = atan(tan(a).*ratio.*1.);
                
                DATA.elmousept.xyr = [0 0 0.5 .* rfix  ratio/(4*rfix)];
                DATA.elmousept.xyr(3) = abs(cos(a)* r(1)  +sin(a) * r(1) * ratio *i);
                DATA.elmousept.xyr(4) = abs(sin(a)* r(2)/ratio  +cos(a) * r(2) *i) ;
                a = atan(tan(a).*ratio*1);
                
%         DrawEllipse(DATA.elmousept,'nonorm','color','r');
            end
%            DrawEllipse(DATA.elmousept,'nonorm','color','g','linewidth',2);
            [a,b] = min(details.x);
            yr(1) = details.y(b);
            [c,b] = max(details.x);
            yr(2) = details.y(b);
            yi = interp1([a c],yr,details.x);
            id = find(details.y <= yi);
            x = details.x(id);
            [y,b] = EllipseForX(DATA.elmousept,x);
            options = optimset('TolFun',1e-6);
            G = DATA.elmousept;
            xr = DATA.elmousept.xyr(3);
            yr = DATA.elmousept.xyr(4);
            G.angle =  atan(tan(DATA.elmousept.angle).*ratio);
            gangles(j) = G.angle;
            e = 1;
            ag = angles(j);
            ag =  atan(tan(DATA.elmousept.angle).*ratio);
            G.xyr(3) = (xr.^e *cos(ag).^e + xr.^e *ratio.^e *  sin(ag).^e).^(1/e);
            G.xyr(4) = (yr.^e *cos(ag).^e + xr.^e * sin(ag).^e/ratio.^e).^(1/e);
            G.xyr(3) = (xr *cos(ag) + xr *ratio *  sin(ag));
            G.xyr(4) = (yr *cos(ag) + yr* sin(ag)./ratio);
            G.xyr(3) =sqrt( ((xr *cos(ag)).^2 + (xr *ratio *  sin(ag)).^2));
            G.xyr(4) = sqrt(((yr *cos(ag)).^2 + (yr* sin(ag)./ratio)).^2);
            if ~isnan(guesspt(1))
                G.angle = guesspt(1);
            end
            if ~isnan(guesspt(2))
                G.xyr(3) = guesspt(2);
            end
            if ~isnan(guesspt(3))
                G.xyr(4) = guesspt(3);
            end
            fprintf('Guess %.2f, %.2f %.2f\n',G.angle,  G.xyr(3), G.xyr(4));
            gxy(j,:) = G.xyr([3 4]);
%            [y,b] = EllipseForX(G,x);
        end
        DATA.toplevel = F;
        set(F,'UserData',DATA);

    function [err, details] =  findellipse(xyr, x,yin)
        E.xyr(1) = 0;
        E.xyr(2) = 0;
        E.xyr(3) = xyr(1);
        E.xyr(4) = xyr(2);
        E.angle = xyr(3);
        [y,b] = EllipseForX(E,x);
        err = sum(abs((y(1,:) - yin).^2));
%        DrawEllipse(E,'nonorm','color','k');
        details.x = x;
        details.y =y;
        details.yin = yin;

     function [err, details] =  findellipseR(xyr, angle,x,yin)
        E.xyr(1) = 0;
        E.xyr(2) = 0;
        E.xyr(3) = xyr(1);
        E.xyr(4) = xyr(2);
        E.angle = angle;
        [y,b] = EllipseForX(E,x);
        err = sum(abs((y(1,:) - yin).^2));
%        DrawEllipse(E,'nonorm','color','k');
        details.x = x;
        details.y =y;
        details.yin = yin;

             function in = InGraph(pt, ax)
        xl = get(ax,'Xlim');
        yl = get(ax,'Ylim');
      in = pt(1,1) > xl(1) & pt(1,1) < xl(2) & pt(1,2) > yl(1) & pt(1,2) < yl(2);
        
            function XYButtonPressed(src, data)
DATA = GetDataFromFig(src);

DATA.ts = now;
start = get(gca,'CurrentPoint');
if InGraph(start,gca)
    mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
%    DATA.elmousept.mode = mode;
    distance = DistanceToEllipse(DATA.elmousept,start(1,1:2));
    DATA.elmousept.steps = 0;
    if ishandle(DATA.elmousept.h)
        delete(DATA.elmousept.h);
    end
    set(gca,'xlimmode','manual','ylimmode','manual');
    if  mode  == 2 %right button press
        DATA.elmousept.down = 3;
        DATA.elmousept.done = 0;
        if DATA.elmousept.shape == 0 %for line, don't move
            DATA.elmousept.start = start(1,1:2);
        end
    elseif distance < 1.05 %test fails for NaN
%pressed inside ellipse, so just tmove it. 
        DATA.elmousept.down = 2;
        DATA.elmousept.pos =[0  0 0 0 ];
        DATA.elmousept.start = start(1,1:2);
        DATA.elmousept.axis = gca;
%set up pos in case released with no drag        
     DATA.elmousept.pos(1) = DATA.elmousept.xyr(1)-DATA.elmousept.xyr(3);
    DATA.elmousept.pos(2) = DATA.elmousept.xyr(2)-DATA.elmousept.xyr(4);
    DATA.elmousept.pos(3) = DATA.elmousept.xyr(1)+DATA.elmousept.xyr(3);
    DATA.elmousept.pos(4) = DATA.elmousept.xyr(2)+DATA.elmousept.xyr(4);
    id = ClassifyPts(DATA.pts, DATA.elmousept);
    plot(DATA.pts(:,1),DATA.pts(:,2),'.');
    hold on;
    plot(DATA.pts(id,1),DATA.pts(id,2),'r.');
    DrawEllipse(DATA.elmousept);
    else
        DATA.elmousept.down = 1;
        DATA.elmousept.done = 0;
        DATA.elmousept.angle = 0;
        DATA.elmousept.pos =[start(1,1) start(1,2) 0 0 ];
        DATA.elmousept.axis = gca;
    end
set(DATA.toplevel,'UserData',DATA);
end


function id = ClassifyPts(x, E)
    yl = get(gca,'Ylim');
    xl = get(gca,'Xlim');
    ar = diff(yl)./diff(xl);
    xy = xyrotate(x(:,1) - E.xyr(1),(x(:,2) - E.xyr(2))./ar,E.angle);
    cn = cos(-E.angle);
    sn = sin(-E.angle);
    r =sqrt(( xy(:,1)./E.xyr(3)).^2 + (xy(:,2)./(E.xyr(4)/ar)).^2);
    id = find(r < 1);
    
function distance = DistanceToEllipse(E, pos);
   
if isempty(E) | ~isfield(E,'pos');
    distance = NaN;
    return;
end


a(1) = (E.pos(3)+E.pos(1))/2; %x center
a(2) = (E.pos(4)+E.pos(2))/2;


if E.shape == 1
angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
r(1) = abs(E.pos(3)-E.pos(1) + i * (E.pos(4)-E.pos(2)))/4;
r2 = r(1)/10;
elseif isfield(E,'xyr')
r(1) = E.xyr(3); %x radius
r(2) = E.xyr(4);
angle = E.angle;
a(1) = E.xyr(1);
a(2) = E.xyr(2);
else
r(1) = (E.pos(3)-E.pos(1))/2; %x radius
r(2) = (E.pos(4)-E.pos(2))/2;
angle = E.angle;
end
xy = pos - a;
xy = xy ./r;
cn = cos(-angle);
sn = sin(-angle);
p(1) = xy(1) * cn + xy(2) * sn;
p(2) = xy(2) * cn - xy(1) * sn;

distance = sum(p.^2);


function XYButtonReleased(src, data)
DATA = GetDataFromFig(src);
if DATA.elmousept.down == 0 
    return;
end
mode = DATA.elmousept.down;
DATA.elmousept.mode = mode;
start = get(gca,'CurrentPoint');
DATA.elmousept.done = 1;
p = DATA.elmousept.pos;
DATA.elmousept.down = 0;
DrawEllipse(DATA.elmousept,'nonorm','color','r');
set(DATA.toplevel,'UserData',DATA);
%touch inside ellispe to make cut. If drawing a line, make cut when
%released


function ScrollWheel(src, evnt)
DATA = GetDataFromFig(src);
DATA.elmousept.angle = DATA.elmousept.angle+0.02*evnt.VerticalScrollCount;
fprintf('Angle %.2f\n',DATA.elmousept.angle);
DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
set(DATA.toplevel,'UserData',DATA);

 function XYButtonDragged(src, data)
DATA = GetDataFromFig(src);
if isfield(DATA,'elmousept') &&  DATA.elmousept.down >0
%    fprintf('D%d,%d %.2f %.2f\n',DATA.elmousept.down,DATA.elmousept.steps,DATA.elmousept.pos(1),DATA.elmousept.pos(3));
    DATA.elmousept.steps = DATA.elmousept.steps +1;
if  DATA.elmousept.down == 1
    start = get(gca,'CurrentPoint');
    DATA.elmousept.pos(3) = start(1,1);
    DATA.elmousept.pos(4) = start(1,2);
    DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
%    drawnow;
%    mytoc(DATA.ts);
elseif  DATA.elmousept.down == 3 %set radius with R mouse button
    if DATA.elmousept.shape == 0 %set radius with R mouse button
    start = get(gca,'CurrentPoint');
    r = abs(start(1,1:2) - DATA.elmousept.xyr(1:2));
    sina = sin(DATA.elmousept.angle);
    cosa = cos(DATA.elmousept.angle);
    r = r * [cosa sina; -sina cosa];
    DATA.elmousept.xyr([3 4]) = r;
    DATA.elmousept.pos(1) = DATA.elmousept.xyr(1)-DATA.elmousept.xyr(3);
    DATA.elmousept.pos(2) = DATA.elmousept.xyr(2)-DATA.elmousept.xyr(4);
    DATA.elmousept.pos(3) = DATA.elmousept.xyr(1)+DATA.elmousept.xyr(3);
    DATA.elmousept.pos(4) = DATA.elmousept.xyr(2)+DATA.elmousept.xyr(4);
    DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
    else %for line R press doesnt move line. Just inverts cluster.sign
    end
elseif  DATA.elmousept.down == 2 %moving ellipse
    start = get(gca,'CurrentPoint');
    if DATA.elmousept.steps > 5
        start = get(gca,'CurrentPoint');
    end
    DATA.elmousept.pos(1) = start(1,1)-DATA.elmousept.start(1)+DATA.elmousept.xyr(1)-DATA.elmousept.xyr(3);
    DATA.elmousept.pos(2) = start(1,2)-DATA.elmousept.start(2)+DATA.elmousept.xyr(2)-DATA.elmousept.xyr(4);
    DATA.elmousept.pos(3) = start(1,1)-DATA.elmousept.start(1)+DATA.elmousept.xyr(1)+DATA.elmousept.xyr(3);
    DATA.elmousept.pos(4) = start(1,2)-DATA.elmousept.start(2)+DATA.elmousept.xyr(2)+DATA.elmousept.xyr(4);
    DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});

end
set(DATA.toplevel,'UserData',DATA);
end

function [h, details]= DrawEllipse(E,varargin)



norm = 1;
j = 1;
makenew = 0;
nused = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nonorm',4);
        norm = 0;
        makenew = 1;
        nused = nused+ 1;
    end
    j = j+1;
end
varargin = varargin(nused:end);
if norm
yl = get(gca,'ylim');
xl = get(gca,'xlim');
x = get(gca,'position');
fp = get(gcf,'position');
ar = abs(diff(x([1 3])).* diff(fp([1 3])) ./   (diff(x([2 4])) .*diff(fp([2 4]))));
ar = abs(diff(x([1 3])) ./   (diff(x([2 4]))));
ar=1;
E.nxyr([1]) = (E.xyr([1]) - mean(xl))./diff(xl);
E.nxyr([2]) = (E.xyr([2]) - mean(yl))./diff(yl);
E.nxyr([3]) = (E.xyr([3]) )./(diff(xl).*ar);
E.nxyr([4]) = (E.xyr([4]) )./diff(yl);
a = E.nxyr(3); %x radius
b = E.nxyr(4);
else
a = E.xyr(3); %x radius
b = E.xyr(4);
end

sn = 0;
cn = 1;
x = linspace(0,a);
y =  sqrt(b.^2 - (x.*b/a).^2);

sn = sin(E.angle);
cn = cos(E.angle);
x = [x fliplr(x) -x 0 fliplr(-x) 0];
y = [y fliplr(-y) -y 0 fliplr(y) 0];
xr = (x .* cn + y .*sn);
yr = (y .* cn - x .*sn);
if norm
    x = xr .*diff(xl).*ar+mean(E.xyr(1));
    y = yr.*diff(yl)+mean(E.xyr(2));
    fprintf('Center %.2f %.2f R%.2f,%.2f NR %.2f %.2f Aspect %.2f\n',E.xyr(1),E.xyr(2),E.xyr(3),E.xyr(4),E.nxyr(3),E.nxyr(4),ar);
else
    x = xr +mean(E.xyr(1));
    y = yr+mean(E.xyr(2));
end

if isfield(E,'h') && ishandle(E.h) && ~makenew
    set(E.h,'Xdata',real(x),'Ydata',real(y));
    h = E.h;
else
    hold on;
    h = plot(real(x),real(y),varargin{:});
    hold off;
end
details.x = real(x);
details.y = real(y);

    
 function [y, details]= EllipseForX(E,x,varargin)


norm = 0;
j = 1;
makenew = 0;
nused = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nonorm',4);
        norm = 0;
        makenew = 1;
        nused = nused+ 1;
    end
    j = j+1;
end
varargin = varargin(nused:end);
if norm
yl = get(gca,'ylim');
xl = get(gca,'xlim');
xp = get(gca,'position');
fp = get(gcf,'position');
ar = abs(diff(xp([1 3])).* diff(fp([1 3])) ./   (diff(xp([2 4])) .*diff(fp([2 4]))));
ar = abs(diff(xp([1 3])) ./   (diff(xp([2 4]))));
E.nxyr([1]) = (E.xyr([1]) - mean(xl))./diff(xl);
E.nxyr([2]) = (E.xyr([2]) - mean(yl))./diff(yl);
E.nxyr([3]) = (E.xyr([3]) )./(diff(xl).*ar);
E.nxyr([4]) = (E.xyr([4]) )./diff(yl);
a = E.nxyr(3); %x radius
b = E.nxyr(4);
else
a = E.xyr(3); %x radius
b = E.xyr(4);
end

sn = 0;
cn = 1;
x = x-E.xyr(1);
sn = sin(-E.angle);
cn = cos(-E.angle);
%ar > 1 means that Y range is bigger than X range
%so that in data space angle is bigger for cos(theta)
ar = 1;
cs = cos(atan(tan(E.angle * ar)));
ss = sin(atan(tan(E.angle/ar)));

aa = (cs^2 +(b^2/a^2)*ss^2);
d = aa;
bb = 2.*x.*cs;
cc = x.^2-b^2.*ss^2; 
xr = x.*cs + sqrt( d.*b^2*ss^2 - (x.^2.*ss^2.*b^2./a^2 ));
xr = xr./d;  %distance along one axes 0- would be x if 0 rotation
xr = (bb+sqrt(bb.^2-4.*aa.*cc))./(2.*aa);
xrb = (bb-sqrt(bb.^2-4.*aa.*cc))./(2.*aa);
xrc = xr;
xrd = xrb;
yr = b .*sqrt(1-(xr.^2/a^2));
yrb = b .*sqrt(1-(xrb.^2/a^2));
yrc = -b .*sqrt(1-(xrc.^2/a^2));
yrd = -b .*sqrt(1-(xrd.^2/a^2));
%xr = [xr xrb xrc xrd];
%yr = [yr yrb yrc yrd]; 
yall = cat(1,yr,yrb,yrc,yrd);
xall = cat(1,xr,xrb,xrc,xrd);

y = yall.*cn - xall .*sn;
y = min(real(y));
ya = cat(1,yr,yrc);
xa = cat(1,xr,xr);
ya = ya.*cn - xa .*sn;
xb = cat(1,xrb,xrb);
yb = cat(1,yrb,yrd);
yb = yb.*cn - xb .*sn;
xcheck(1,:) = (xr.*cn) + (yr.*sn);
xcheck(2,:) = (xr.*cn) - (yr.*sn);
xdiff = abs(xcheck -repmat(x,2,1));
[a,b] = min(xdiff);
id = sub2ind(size(xcheck),b,1:size(xcheck,2));

xcheck(1,:) = (xrb.*cn) + (yrb.*sn);
xcheck(2,:) = (xrb.*cn) - (yrb.*sn);
xdiff = abs(xcheck -repmat(x,2,1));
[a,b] = min(xdiff);
bid = sub2ind(size(xcheck),b,1:size(xcheck,2));
details.x = real([x]);
details.xcheck = xcheck(id);
y = [ya(id); yb(bid)];
fdetails.x = real(x);

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

