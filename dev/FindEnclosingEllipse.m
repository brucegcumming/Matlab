function fittedparams = FindEnclosingEllipse(xy, id, cluster, aspectratio, varargin)

cid = find(id == cluster);
nid = find(id ~= cluster);
c = mean(xy(cid,:));
guess(1:2) = c;
guess(3) = std(xy(cid,1)) .* 2;
guess(4) = std(xy(cid,2)) .* 2;
guess(5) = 0;
state.aspectratio = aspectratio;
maxiter = 1000;

hold off;
plot(xy(:,1),xy(:,2),'.');
hold on;
plot(xy(cid,1),xy(cid,2),'r.');
state.fig = gcf;
if isappdata(gcf,'params')
    rmappdata(gcf,'params');
end
options = optimset('MaxFunEvals',100000,'maxiter',maxiter,'display','off');
[fittedparams,fval,exitflag, output] = fminsearch(@Minimise,guess,options,xy,id,state);
E.xyr(1:4) = fittedparams(1:4);
E.angle = fittedparams(5);
E.color = 'r';


DrawEllipse(E);
return;

for j = 1:10
    a = j./pi;
    v = [cos(a) sin(a)];
    d = (xy-c) * v';
    r(j) = max(d(cid));
    or(j) = max(d(nid));
end

[ra,id] = max(r);
E.xyr(c1) = c(1);
E.xyr(2) = c(2);
E.xyr(3) = max(r);


function SSD = Minimise(params,xy,clst, state)

cx = params(1);
cy = params(2);
rx = params(3);
ry = params(4);
a = params(5);
xys = xyrotate(xy(:,1)-cx,(xy(:,2)-cy) ./state.aspectratio,a);
r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./state.aspectratio)).^2;

id = find(r < 1);
nid = find(r>1);
SSD = sum(clst(id) ~= 2) + sum(clst(nid) ==2);
if isfigure(state.fig)
    allparams = getappdata(state.fig,'params');
    allparams = [allparams; params];
    setappdata(state.fig,'params',allparams);
end

function h= DrawEllipse(E,varargin)

if isfield(E,'pos')
a = (E.pos(3)-E.pos(1))/2; %radius 
b = (E.pos(4)-E.pos(2))/2;
cntr(1) = mean(E.pos([1 3]));
cntr(2) = mean(E.pos([2 4]));
elseif isfield(E,'xyr')
    a = E.xyr(3);
    b = E.xyr(4);
    cntr(1) = E.xyr(1);
    cntr(2) = E.xyr(2);
end
    
sn = 0;
cn = 1;

x = linspace(0,a);
sn = sin(E.angle);
cn = cos(E.angle);
if isfield(E,'aspectratio') && E.aspectratio > 0
    b  = b ./E. aspectratio;
    y =  sqrt(b.^2 - (x.*b/a).^2);
else
    y =  sqrt(b.^2 - (x.*b/a).^2);
end

x = [x fliplr(x) -x fliplr(-x)];
y = [y fliplr(-y) -y fliplr(y)];
xr = (x .* cn + y .*sn);
yr = (y .* cn - x .*sn);
x = xr+cntr(1);
if isfield(E,'aspectratio') && E.aspectratio > 0
    y = yr.*E.aspectratio+cntr(2);
else
    y = yr+cntr(2);
end

if isfield(E,'h') & ishandle(E.h)  & get(E.h,'parent') == gca
%    fprintf('Using handle %f\n',get(E.h,'parent'));
    set(E.h,'Xdata',x,'Ydata',y,'color',E.color);
    h = E.h;
else
   hold on;
    h = plot(real(x),real(y),'color',E.color,varargin{:});
    hold off;
end
