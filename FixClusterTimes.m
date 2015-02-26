function FixClusterTimes(Clusters, ClusterDetails)
%fix mismatch between Clusters{c}.times and ClusterDetails{c}
plottimes = 1;

for p = 1:length(Clusters)
    C = Clusters{p};
    CD = ClusterDetails{p};
    cid = find(CD.clst ==2);
    if length(C.times) ~= length(cid)
        fprintf('Checking probe %d\n',p);

        [a,mid,b,cmid] = MatchTimes(CD.t(cid),C.times,0.001);
        [a,id,b,c] = MatchTimes(C.times,CD.t(cid),0.0002);
        fprintf('%d not in clst, %d not in times\n',length(mid),length(id))
        if plottimes %plot results
            GetFigure('Missing events');
            hold off;
            plot(C.times(b),CD.t(cid(a)),'.');
            hold on;
            plot(C.times(c),CD.t(cid(id)),'r.');
            plot(C.times(mid),CD.t(cid(cmid)),'g.'); %in Clsuter.Times but not clst
            legend('matches','missing from clst','missing from Cluster.times');
            title(sprintf('Porbe%d',p));
        end
        Clusters{p}.times = sort([C.times CD.t(cid(id))]);
    end
end
