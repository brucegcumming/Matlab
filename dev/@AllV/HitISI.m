function HitISI(a,b, t)DATA = GetDATAFromFig(a);sid(2) = find(DATA.t == t);fprintf('t=%.2f event %d\n',t,sid(2));if size(DATA.t,2) == size(DATA.clst,2)id = find(DATA.clst == DATA.clst(sid(2)) & DATA.t < t);elseid = find(DATA.clst == DATA.clst(sid(2)) & DATA.t' < t);endsid(1) = id(end);if length(sid)AllV.PlotSpikes(DATA,sid);AllV.ReplotPCs(DATA,[],'setid',sid);AllV.SetFigure(DATA.tag.fullv, DATA);AllV.PlotFullV(DATA,[DATA.t(sid(1))-0.0005 DATA.t(sid(2))+0.0005]);end