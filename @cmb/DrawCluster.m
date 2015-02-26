function h =  DrawCluster(cluster, color)

h =0;
if isempty(cluster) | ~isfield(cluster,'x')
return;
end
if isfield(cluster,'h') & ishandle(cluster.h)
delete(cluster.h);
end
tmp.r = [cluster.x(2) cluster.y(2)];
tmp.c = [cluster.x(1) cluster.y(1)];
tmp.xrange = cluster.x(3);
tmp.yrange = cluster.y(3);
tmp.angle = -cluster.angle;
tmp.color = color;
tmp = cmb.myellipse(tmp);
h = tmp.lasth;


