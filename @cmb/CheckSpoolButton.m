function CheckSpoolButton(DATA)
if ~isfigure(DATA.xyfig)
return;
end
spoolbutton = findobj(DATA.xyfig,'Tag','SpoolSpikes');
if  isempty(spoolbutton)
return;
end
if DATA.spooling == 2 || DATA.firsttrial > 1
set(spoolbutton,'backgroundcolor','r');
else
set(spoolbutton,'backgroundcolor','w');
end

