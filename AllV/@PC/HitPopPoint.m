function HitPopPoint(a,b, ex, p, cell)DATA = GetDataFromFig(a);Clusters = getappdata(DATA.toplevel,'Clusters');C = Clusters{ex}{p};fprintf('Expt %d Probe %d\n',C.exptno,p);DATA.selectprobe = zeros(length(Clusters),DATA.nprobes);DATA.selectprobe(ex,p) = 1;if nargin > 4 && cell > 0    DATA.currentcell = cell;endPC.ShowData(DATA, ex, p);