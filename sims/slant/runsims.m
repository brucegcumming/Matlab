%
% good mock up of gregs example cell
[res, details] = AbsSlant('rfp',[-1 -1],'dtscale',[0.4 60],'dxs',[0.25:0.05:0.45]);
AbsSlant(details,res,'plot')
