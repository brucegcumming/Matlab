function AddProbeList(DATA)
% Add A list of probes for multiple probe recording.
%
ch = DATA.plot.ch;
cw = DATA.plot.cw;
wsiz = [DATA.plot.cw*40,DATA.plot.ch*40];

ival = 1;

it =  findobj(DATA.toplevel,'Tag','UseProbe1');
if ~isempty(it)
%    delete(it); %for debugging
return;
end
%
% if probes were selected on the command line, useprobe will exist.
% Otherwise not.
if isfield(DATA.plot,'useprobe')
up = DATA.plot.useprobe;
else
up = zeros(size(DATA.probelist));
end
for j = 1:length(DATA.probelist)
bp = [0.94 (j-1)./(length(DATA.probelist)+1) 0.07 0.05];
uicontrol(DATA.toplevel,'Style', 'checkbox',...
'String', num2str(j), 'Tag', sprintf('UseProbe%d',j),'units','norm', 'Position', bp,'value',up(j),...
'Callback',@cmb.Update);
end
bp = [0.94 (j)./(length(DATA.probelist)+1) 0.06 0.04];
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['cmb.combine(''combine'',''Tag'',''' DATA.tag.top ''');'],...
'String', 'Combine','Units','norm', 'Position', bp);
uicontrol(DATA.toplevel,'Style', 'pushbutton',...
'String', '+/-', 'Tag', 'ToggleProbe', 'units','norm','Position', bp,'value',0,...
'Callback',{@cmb.MainMenu,1});
ws = get(DATA.elst,'position');
ws(3) = 0.92;
set(DATA.elst,'position',ws);
ws = get(DATA.clst,'position');
ws(3) = 0.92;
set(DATA.clst,'position',ws);
ws = get(DATA.toplevel,'position');
ws(3) = ws(3) + cw * 4;
set(DATA.toplevel,'position',ws);


