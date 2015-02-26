function FullVButtonReased(a,b)

DATA = GetDataFromFig(a);


axdata = get(gca,'UserData');
if isfield(axdata,'xrange')
    set(gca,'xlim',axdata.xrange);
end