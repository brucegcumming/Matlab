function FullVButtonPressed(a,b)

DATA = GetDataFromFig(a);


mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});

if mode == 2    
axdata.xrange = get(gca,'xlim');

pos = get(gca,'CurrentPoint');
x = pos(1,1);
set(gca,'xlim',[x-0.1 x+0.1]);
set(gca,'UserData',axdata);
else
    set(gca,'UserData',[]);
end