function MakeCellList(a,b, varargin)
DATA = GetDataFromFig(a);

if ~isfield(DATA,'CellList')
DATA = cmb.LoadCellFile(DATA);
set(DATA.toplevel,'UserData',DATA);
end
wsc = DATA.wsc;
SPACE = 3 * wsc(1);
VSPACE = 5 * wsc(2);
ch = DATA.plot.ch;
cw = DATA.plot.cw;
tag = DATA.tag.celllist;

cntrl_box = findobj('Tag',tag,'Name','Cell List');
if ~isempty(cntrl_box)
figure(cntrl_box);
return;
end
dat.parentfigtag = DATA.tag.top;
% standard figure. Will Plot Map
cntrl_box = figure('Menubar', 'none',...
'NumberTitle', 'off', 'Tag',tag,'Name','Cell List','UserData',dat);

bp = [10 10 cw*4 ch*2];
uicontrol(gcf,'Style', 'pushbutton', 'Callback', {@cmb.AddCellToList, 'Tag' DATA.tag.top},...
'String', 'Add', 'Position', bp);

bbp = [10 bp(2)+ch*2 cw*4 ch*2];
uicontrol(gcf,'Style', 'pushbutton', 'Callback', @cmb.NextList,...
'String', 'Next', 'Position', bbp);
bbp(2) = bbp(2)+ch*2;

uicontrol(gcf,'Style', 'pushbutton', 'Callback', @cmb.MakeCellTemplate,...
'String', 'Tmplt', 'Position', bbp);

bbp(2) = bbp(2)+ch*2;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', @cmb.SaveCellList,...
'String', 'Save', 'Position', bbp);

bbp(2) = bbp(2)+ch*2;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', @cmb.AutoFillCellList,...
'String', 'Auto', 'Position', bbp);


np = length(DATA.probelist);
for j = 1:np
pstrs{j} = num2str(DATA.probelist(j));
pstrs{j+np} = sprintf('%d/2',DATA.probelist(j));
end
bp(1) = bp(1)+bp(3)+SPACE;
cellx = bp(1);

bp(1) = cellx;
bp(3) = cw*5;
uicontrol(gcf,'Style', 'Text','String', 'Quality','Position', bp);
for j = 1:DATA.state.listlen
bp(1) = bp(1)+bp(3)+SPACE;
uicontrol(gcf,'Style', 'pop',...
'String', 'Nothing|MU|MU+|Poor|OK|Good|VGood|Excellent|Automatic', 'Tag', sprintf('CellQuality%d',j), 'Position', bp,'value',1);
end

bp(1) = bp(1)+bp(3)+SPACE;
uicontrol(gcf,'Style', 'Text','String', 'Plot','Position', bp);
bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = cw*5;
uicontrol(gcf,'Style', 'pop',...
'String', 'Quality|Number|One Cell|Cells|CQuality|Both|dprime|Q-Dprime|CellDprime', 'Tag', 'PlotType', 'Position',bp, 'callback', @cmb.RePlotCellList);
bp(1) = bp(1)+bp(3)+SPACE;


bp(1) = cellx;
bp(3) = cw*5;
bp(2) = bp(2)+bp(4)+VSPACE;
uicontrol(gcf,'Style', 'Text','String', 'Probe','Position', bp);
for j = 1:DATA.state.listlen
bp(1) = bp(1)+bp(3)+SPACE;
uicontrol(gcf,'Style', 'pop',...
'String', pstrs, 'Tag', sprintf('CellProbeAdd%d',j), 'Position', bp,'value',DATA.probe);
end
bp(1) = cellx;
bp(2) = bp(2)+bp(4)+VSPACE;
bp(3) = cw*5;
uicontrol(gcf,'Style', 'Text','String', 'Cell','Position', bp);

for j = 1:DATA.state.listlen
bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = cw*5;
'String', num2str([0:40]), 'Tag', sprintf('CellNumber%d',j), 'Position', bp,'value',1,...
uicontrol(gcf,'Style', 'pop',...
'String', '0|1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25', 'Tag', sprintf('CellNumber%d',j), 'Position', bp,'value',1,...
'callback',@cmb.SetCellNum);
end
cmb.PlotCellList(DATA);
set(gca,'position',[0.1 0.2 0.85 0.78])

