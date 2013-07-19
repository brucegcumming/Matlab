function ClearPlot(a)
%ClearPlot() Clears a plot - deleting all axes 
%unlike clf does NOT remove menus etc that the user has added
%ClearPlot(h) Clears the plot in handle h


if nargin
  delete(allchild(a));
else
    subplot(1,1,1);
  delete(allchild(gca));
  lc = get(gcf,'children');
  for j=1:length(lc)
      if strcmp(get(lc(j),'type'),'axes')
          delete(lc(j));
      end
  end
end  

