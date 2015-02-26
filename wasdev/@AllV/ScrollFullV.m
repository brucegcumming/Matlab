function ScrollFullV(src, evnt)
DATA = GetDataFromFig(src);

if src ~= gcf
    return;
end

xr = get(gca,'xlim');
if sign(evnt.VerticalScrollCount) > 0
    xr = xr + diff(xr)/10;
    set(gca,'xlim',xr);
elseif sign(evnt.VerticalScrollCount) < 0 
    xr = xr - diff(xr)/10;
    set(gca,'xlim',xr);
end
