function DelClusterButton(caller,b)
it = findobj(get(caller,'parent'),'Tag','Clusterid');
c = get(it,'value');
cmb.DeleteCluster(c, caller);

