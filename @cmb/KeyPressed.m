function KeyPressed(a, ks)
global mousept;

if ~isfield(mousept,'mode')
    return;
end
parent = get(a,'parent');
mousept.mode;
if strmatch(ks.Key,'delete') & mousept.mode == 5
cmb.DeleteCluster(mousept.cluster, a);
mousept.mode = 0;
mousept.angle = 0;
if mousept.lasth & ishandle(mousept.lasth)
delete(mousept.lasth);
end
elseif ks.Key == 'n'
cmb.NewCluster(a);
elseif strmatch(ks.Key,'add')
it = findobj(a,'Tag','Clusterid');
c = get(it,'value');
set(it,'value',c+1);
elseif strmatch(ks.Key,'subtract')
it = findobj(a,'Tag','Clusterid');
c = get(it,'value');
if c > 1 
set(it,'value',c-1);
end
elseif ks.Key == 'r'
mousept.angle = mousept.angle+0.02;
mousept.mode = 10;
if mousept.lasth & ishandle(mousept.lasth)
delete(mousept.lasth);
end
mousept = cmb.myellipse(mousept,[0 0; 0 0 ]);
end

