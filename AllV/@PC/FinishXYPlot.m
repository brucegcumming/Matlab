function FinishXYPlot(ax, DATA, e,p)    axis('tight');    set(ax,'xtick',[],'ytick',[]);     tag = PC.GetFigureTag(ax);    if strcmp(tag,DATA.tag.allxy)        set(ax,'ButtonDownFcn',{@PC.HitXYPlot, e,p});    end    