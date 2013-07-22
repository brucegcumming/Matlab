function rundpsf(flag)

list = '/d/bgc/data/dpp/dpsflist';

global top fstrings fign

    if nargin > 1  | nargin < 1 | strcmp(flag,'start')
      allstrings = textread(list,'%s');
      fstrings= allstrings(strmatch('/usr',allstrings));
%%     top = dialog('WindowStyle','normal');
figure('Position', [0 800 200 200], 'Menubar', 'none',...
       'NumberTitle', 'off', 'Tag','TopLevel');
     lst = uicontrol(gcf, 'Style','listbox','String',fstrings,...
		'Callback', 'rundpsf(''set''),','Tag','TheList',...
		'Position',[20 20 170 90]);
uicontrol(gcf,'Style', 'checkbox', 'Callback', 'rundpsf(''smooth'')',...
'String', 'Force Read', 'Tag', 'ReRead', 'Position', [20 130 90 20]);
uicontrol(gcf,'style','pop','string','Pcolor|SF|dPhase|AC|Monoc', ...
		    'Callback', 'rundpsf(''smooth'')', 'Tag','plottype',...
		    'position',[20 150 50 20]);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'rundpsf(''all'')',...
'String', 'All', 'Position', [120 150 50 20]);
uicontrol(gcf,'Style', 'checkbox', 'Callback', 'rundpsf(''smooth'')',...
'String', 'Fit', 'Tag', 'ShowFit', 'Position', [70 150 50 20]);
uicontrol(gcf,'Style', 'checkbox', 'Callback', 'rundpsf(''shade'')',...
'String', 'Shade', 'Tag', 'Shading', 'Position', [110 130 60 20]);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'rundpsf(''next'')',...
'String', 'Next', 'Position', [20 110 50 20]);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'rundpsf(''prev'')',...
'String', 'Prev', 'Position', [70 110 50 20]);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'rundpsf(''close'')',...
'String', 'Close', 'Position', [120 110 50 20]);
fign = figure('Tag','Contours');
elseif strcmp(flag,'next')
     it = findobj('Tag','TheList');
     n = get(it, 'value');
     last = size(fstrings);
     if n < last(1)
       n = n+1;
     set(it, 'value',n);
     plotmember(n, fstrings);
     end
elseif strcmp(flag,'all')
     for j=1:length(fstrings)
       plotmember(j, fstrings);
     end
elseif strcmp(flag,'prev')
     it = findobj( 'Tag','TheList');
     n = get(it, 'value');
     if(n > 1)
     n = n-1;
     set(it, 'value',n);
     plotmember(n, fstrings);
     end
elseif strcmp(flag,'close')
     it = findobj(gcf, 'Tag','TopLevel');
     close(it);
     it = findobj(gcf, 'Tag','Contours');
     close(it);
    else %clicked on line
     it = findobj(gcf, 'Tag','TheList');
     n = get(it, 'value'); 
    plotmember(n, fstrings);
end

function plotmember(n, fstrings)

file = fstrings(n);
thefile = file{1};
it = findobj('Tag','Contours');
figure(it);
reread = get(findobj( 'Tag','ReRead'),'value');
acfile = strrep(thefile, 'grating','rds');
acfile = strrep(acfile, 'DPxSF','OXAC');

if ~exist(acfile)
  acfile = strrep(acfile, 'rds','rls');
  if ~exist(acfile)
    acfile = strrep(thefile, 'grating','rds');
    acfile = strrep(acfile, 'DPxSF','ODX');
    if ~exist(acfile)
      acfile = strrep(acfile, 'rds','rls');
    end
  end
end
    
DPP = Read_DPP_data(thefile,reread);
save('tmp');
subplot(1,1,1);
hold off;
plottype = get(findobj('Tag','plottype'),'value');
if plottype == 4 & exist(acfile)
    OXAC = Read_DPP_data(acfile,reread);
    AC = Plot_AC_data(OXAC);
    hold on;
else
  AC.disprange(1) = -1;
  AC.disprange(2) = 1;
end
Plot_DPSF_data(DPP,plottype,AC);

shade = get(findobj( 'Tag','Shading'),'value');
if(shade)
  shading('interp');
else
    shading('flat');
end

