function ClrSpkWin(caller,b)
%DATA = cmb.combine('getstate');
if isstruct(caller)
    DATA = caller;
else
    DATA = GetDataFromFig(caller);
end
cmb.SetFigure(DATA.tag.clusterxy,DATA);
ym = get(gca,'ylim');
xm = get(gca,'xlim');
hold off;
plot(0,0,'+');
set(gca,'ylim',ym);
set(gca,'xlim',xm);
hold on;
if ~isstruct(caller) % from a mouse button
    DATA = cmb.DrawClusters(DATA,DATA.cluster, 0);
    set(DATA.toplevel,'UserData',DATA);
end


