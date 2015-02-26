
function runtf(flag)

list = '/d/bgc/data/dtsd/alllist';

global top fstrings fign

    if nargin > 1  | nargin < 1 | strcmp(flag,'start')
     fstrings = textread(list,'%s');
%%     top = dialog('WindowStyle','normal');
figure('Position', [0 800 200 200], 'Menubar', 'none',...
       'NumberTitle', 'off', 'Tag','TopLevel');
     lst = uicontrol(gcf, 'Style','listbox','String',fstrings,...
		'Callback', 'runtf(''set'')','Tag','TheList',...
		'Position',[20 20 170 90]);
uicontrol(gcf,'Style', 'checkbox', 'Callback', 'runtf(''smooth'')',...
'String', 'Interpolation', 'Tag', 'Interpolation', 'Position', [20 130 90 20]);
uicontrol(gcf,'Style', 'checkbox', 'Callback', 'runtf(''shade'')',...
'String', 'Shade', 'Tag', 'Shading', 'Position', [110 130 60 20]);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'runtf(''next'')',...
'String', 'Next', 'Position', [20 110 50 20]);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'runtf(''prev'')',...
'String', 'Prev', 'Position', [70 110 50 20]);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'runtf(''close'')',...
'String', 'Close', 'Position', [120 110 50 20]);
%fign = figure('Tag','Contours');
elseif strcmp(flag,'next')
     it = findobj(gcf, 'Tag','TheList');
     n = get(it, 'value');
     last = size(fstrings)
     if n < last(1)
     n = n+1;
     set(it, 'value',n);
     plotmember(n, fstrings, fign);
     end
elseif strcmp(flag,'prev')
     it = findobj(gcf, 'Tag','TheList');
     n = get(it, 'value');
     if(n > 1)
     n = n-1;
     set(it, 'value',n);
     plotmember(n, fstrings, fign);
     end
elseif strcmp(flag,'close')
     it = findobj(gcf, 'Tag','TopLevel');
     close(it);
     it = findobj(gcf, 'Tag','Contours');
     close(it);
else
     it = findobj(gcf, 'Tag','TheList');
     n = get(it, 'value');
     plotmember(n, fstrings);
end

function plotmember(n, fstrings)

file = fstrings(n);
thefile = file{1}
     plottf(thefile);

     
function plottf(file)

contents = textread(file,'%s')