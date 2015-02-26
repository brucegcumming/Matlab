function OptimizeDprimeHit(a,b)
DATA = GetDataFromFig(a);
if isfield(DATA.plot,'useprobe') & sum(DATA.plot.useprobe)
DATA = cmb.GetAllProbeFig(DATA);
[nr, nc] = Nsubplots(length(DATA.probelist));
probes = find(DATA.plot.useprobe);
for j = probes;
subplot(nr,nc,j);
DATA.probe = DATA.probelist(j);
DATA.ispk = DATA.AllClusters(j).spklist;
C = cmb.OptimizeDprime(DATA);
DATA.cluster{DATA.currentcluster,DATA.probe} = C;
[DATA, newd, nc] = SetSpkCodes(DATA,DATA.ispk,DATA.probe,2);
cmb.DrawClusters(DATA,DATA.cluster, 0);
end
else

DATA.ispk = DATA.Expts{DATA.currentexpt(1)}.gui.spks;
C = cmb.OptimizeDprime(DATA);
DATA.cluster{DATA.currentcluster,DATA.probe} = C;
[DATA, newd, nc] = SetSpkCodes(DATA,DATA.ispk,DATA.probe,2);
cmb.DrawClusters(DATA,DATA.cluster, 0);
end
set(DATA.toplevel,'UserData',DATA);

