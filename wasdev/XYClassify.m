function cluster = XYClassify(x,y,E, varargin)

theta = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
exy = xyrotate(E.pos([1 3]),E.pos([2 4]),theta);
cluster.y = exy(:,2);
crit = mean(exy(:,1));
cluster.angle = theta;
cluster.crit = crit;
xy = xyrotate(x,y,theta);
if mean(xy(:,1)) < 0
cluster.id = find(xy(:,1) > mean(exy(:,1)));
else
cluster.id = find(xy(:,1) < mean(exy(:,1)));
end
