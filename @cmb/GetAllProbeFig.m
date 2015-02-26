function DATA = GetAllProbeFig(DATA)
if DATA.plot.prettyfigs
GetFigure('AllProbeXYFig');
elseif isfigure(DATA.figs.allprobes)
GetFigure(DATA.figs.allprobes);
else
[DATA.figs.allprobes, isnew] = GetFigure(DATA.tag.allprobes);
dat.parentfigtag = DATA.tag.top;
set(DATA.figs.allprobes,'UserData',dat);
if isnew
bp = [10 10 DATA.plot.cw*7 20];
uicontrol(DATA.figs.allprobes,'Style', 'checkbox',...
'String', 'Density', 'Tag', 'AllDensity', 'Position', bp,'value',0,...
'Callback',@cmb.Update);

abp = bp;
bp(2) = bp(2)+DATA.plot.ch * 2.2;
uicontrol(DATA.figs.allprobes,'Style', 'pushbutton','String','Set','Position',bp,'Callback', @cmb.SetExptClusters);
bp(1) = bp(1)+bp(3);
bp(3) = DATA.plot.cw * 6;
uicontrol(DATA.figs.allprobes,'Style', 'pushbutton', 'Callback', @cmb.NextList, ...
'String', 'Next', 'Position',bp);
bp = abp;
for j = 1:length(DATA.probelist)
bp(1) = bp(1)+bp(3);
if j < 10
bp(3) = DATA.plot.cw * 3;
else
bp(3) = DATA.plot.cw * 4;
end
uicontrol(DATA.figs.allprobes,'Style', 'checkbox',...
'String', sprintf('%d',j), 'Tag', sprintf('MarkProbe%d',j), 'Position', bp,'value',0,...
'Callback',@cmb.Update);
end
bp(1) = bp(1)+bp(3);
bp(3) = DATA.plot.cw * 5;

uicontrol(DATA.figs.allprobes,'Style', 'pushbutton',...
'String', 'Mark', 'Position', bp,'value',0,...
'Callback',@cmb.MarkProbes,'UserData',dat);

bp(1) = bp(1)+bp(3);
bp(3) = DATA.plot.cw * 5;
uicontrol(DATA.figs.allprobes,'Style', 'pushbutton',...
'String', 'Hide SpikeV', 'Position', bp,'value',0,...
'Callback',@cmb.HideSpikes,'UserData',dat);
end
end


