function AddClusterIdButton(a,b)
DATA = GetDataFromFig(a);
it = findobj(DATA.xyfig,'Tag','OptimizeDprime');
fc = findobj(DATA.xyfig,'Tag','ForceClusterid');
if isempty(fc)
set(a,'label','Stop cluster ids forcing');
bp = get(it,'position');
bp(2) = bp(2)+bp(4);

uicontrol(DATA.xyfig,'style','pop','string','1|2|3|4|5|6|7','Position',bp,'Tag','ForceClusterid',...
'Callback',@cmb.Update);
else
set(a,'label','Force cluster ids');
delete(fc);
end




