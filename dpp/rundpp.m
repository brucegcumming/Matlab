function rundpp(flag)

list = '/d/bgc/data/dpp/alllist';

global top fstrings fign

    if nargin > 1  | nargin < 1 | strcmp(flag,'start')
     fstrings = textread(list,'%s');
%%     top = dialog('WindowStyle','normal');
figure('Position', [100 100 200 200], 'Menubar', 'none',...
       'NumberTitle', 'off', 'Tag','TopLevel','Name','DPP');
     lst = uicontrol(gcf, 'Style','listbox','String',fstrings,...
		'Callback', 'rundpp(''set''),','Tag','TheList',...
		'Position',[20 20 170 90]);
uicontrol(gcf,'Style', 'checkbox', 'Callback', 'rundpp(''smooth'')',...
'String', 'Force Read', 'Tag', 'ReRead', 'Position', [20 130 90 20]);
uicontrol(gcf,'style','pop','string','Pcolor|Phase1|Phase2|AC|DPSF-pcolor|DPSF:SF|DPSF:DP|DPSF_AC|DPSF', ...
		    'Callback', 'rundpp(''smooth'')', 'Tag','plottype',...
		    'position',[20 150 50 20]);
uicontrol(gcf,'Style', 'checkbox', 'Callback', 'rundpp(''smooth'')',...
'String', 'Fit', 'Tag', 'ShowFit', 'Position', [70 150 50 20]);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'rundpp(''all'')',...
'String', 'All', 'Tag', 'runall', 'Position', [140 150 50 20]);
uicontrol(gcf,'Style', 'checkbox', 'Callback', 'rundpp(''shade'')',...
'String', 'Shade', 'Tag', 'Shading', 'Position', [110 130 60 20]);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'rundpp(''next'')',...
'String', 'Next', 'Position', [20 110 50 20]);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'rundpp(''prev'')',...
'String', 'Prev', 'Position', [70 110 50 20]);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'rundpp(''close'')',...
'String', 'Close', 'Position', [120 110 50 20]);
fign = figure('Tag','Contours');
elseif strcmp(flag,'next')
     it = findobj('Tag','TheList');
     n = get(it, 'value');
     last = size(fstrings)
     if n < last(1)
       n = n+1;
     set(it, 'value',n);
     plotmember(n, fstrings);
     end
elseif strcmp(flag,'prev')
     it = findobj( 'Tag','TheList');
     n = get(it, 'value');
     if(n > 1)
     n = n-1;
     set(it, 'value',n);
     plotmember(n, fstrings);
     end
elseif strcmp(flag,'all')
     it = findobj( 'Tag','TheList');
     for n = 1:length(fstrings)
       fprintf('%s\n',fstrings{n});
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
acfile = strrep(thefile, '2grating','rds');
dpsffile = strrep(thefile, '2grating.DPP','grating.DPxSF');
acfile = strrep(acfile, 'DPP','OXAC');

if ~exist(acfile)
  acfile = strrep(acfile, 'rds','rls');
  if ~exist(acfile)
    acfile = strrep(thefile, '2grating','rds');
    acfile = strrep(acfile, 'DPP','ODX');
    if ~exist(acfile)
      acfile = strrep(acfile, 'rds','rls');
    end
  end
end
DPP = Read_DPP_data(thefile,reread);
plottype = get(findobj('Tag','plottype'),'value');

if plottype >= 4 & exist(acfile)
    hold on; 
    OXAC = Read_DPP_data(acfile,reread);
    if(plottype == 4)
      Plot_AC_data(OXAC);
      hold on;
    end
end

if plottype >= 5 & exist(dpsffile)
  DPSF = Read_DPP_data(dpsffile,reread); %Read_DPP works for DPSF too
  hold off;
  idx = findstr(OXAC.title,'dx');
  col = 1 + (idx -1)/3;
  dxs = OXAC.data(col,:);
  OXAC.disprange(1) = min(dxs);
  OXAC.disprange(2) = max(dxs);
  if(plottype == 8)
    Plot_AC_data(OXAC);
    hold on;
  end
  Plot_DPSF_data(DPSF,plottype-4,OXAC);
elseif plottype >= 5
  fprintf('No %s file\n',dpsffile);
else
  Plot_DPP_data(DPP,plottype);
end


shade = get(findobj( 'Tag','Shading'),'value');
if(shade)
  shading('interp');
else
    shading('flat');
end

