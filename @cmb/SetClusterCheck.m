function SetClusterCheck(DATA)
cmb.SetCheck('UseCluster1',ismember(1,DATA.spikelist),DATA.toplevel);
cmb.SetCheck('UseCluster2',ismember(2,DATA.spikelist),DATA.toplevel);
cmb.SetCheck('UseCluster3',ismember(3,DATA.spikelist),DATA.toplevel);
cmb.SetCheck('UseCluster4',ismember(4,DATA.spikelist),DATA.toplevel);



