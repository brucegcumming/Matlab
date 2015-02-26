function DATA = CalcOnlineShapes(DATA)
ts = now;
DATA = cmb.ClearExptClusters(DATA);
DATA = cmb.LoadOnlineClusters(DATA, cmb.ClusterFile(DATA,'getonline'));
DATA = cmb.ReClassifyAll(DATA,'mkmean');
cmb.SaveSpikeShape(DATA, strrep(DATA.meanspkfile,'spks','ospk'));
fprintf('Took %.2f sec\n',(now-ts)*(24*60*60));


