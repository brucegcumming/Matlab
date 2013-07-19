function AppDataPopup(a, b, fcn, datalabel)
% AppDataPopup(F, datalabel)
% Builds a popup to edit elements of a structre 
% stored as appdata with tag datalabel
%args 2 and 2
DATA = GetDataFromFig(a);
  fontsize = DATA.gui.fontsize;
  if ~isempty(b)
      datalabel = b;
      GetFigure(datalabel,DATA.toplevel);
      SetFigPos(DATA,datalabel);
      ClearPlot;
      S = getappdata(DATA.toplevel,b)
      f = fields(S);
      nr = length(f)+1;
      nc = 2;
      bp(3) = 0.95./nc;
      bp(4) =0.95./nr;
      for j = 1:length(f)
          bp(1) = 0.01;
          bp(2) = 0.01+(j).* 1./nr;
          uicontrol(gcf,'style','Text','String',f{j},...
              'Tag',f{j},...
              'fontsize',fontsize,...
              'units','normalized','position',bp);
          bp(1) = 0.5;
          uicontrol(gcf,'style','edit','String',num2str(DATA.gui.(f{j})),...
              'fontsize',fontsize,...
              'units','normalized','position',bp);
      end
      bp(2) = 0.01;
      uicontrol(gcf,'style','pushbutton','String','Apply',...
          'units','normalized',...
          'fontsize',DATA.gui.fontsize,...
          'Callback',{@AppDataPopup, 'apply' datalabel},'position',bp);
  elseif strcmp(fcn,'apply')
      S = getappdata(DATA.toplevel,'datalabel')
      f = fields(S);
      for j = 1:length(f)
          DATA.gui.(f{j}) = Text2Val(get(a,'parent'),f{j});
      end
      set(DATA.toplevel,'UserData',DATA);
  end
  
function value = Text2Val(F, tag)
    it = findobj(F,'tag',tag);
    if ~isempty(it)
        value = str2num(get(it,'string'));
    end
