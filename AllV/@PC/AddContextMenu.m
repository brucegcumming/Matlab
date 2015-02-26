function cmenu = AddContextMenu(DATA, type)    cmenu = uicontextmenu;        if strcmp(type,'subplot')        uimenu(cmenu,'label','spool','Callback',{@PC.SetCellFromLine, j-1,  'spool'});        uimenu(cmenu,'label','->FullV','Callback',{@PC.SetCellFromLine, j-1,  'getfullv'});        uimenu(cmenu,'label','cut cluster1 (&Line)','Callback',{@PC.XYCluster, 'lines', 'useaxdata'});        uimenu(cmenu,'label','cut ellipse&1','Callback',{@PC.XYCluster, 'ellipse1', 'useaxdata'});        uimenu(cmenu,'label','cut ellipse&2','Callback',{@PC.XYCluster, 'ellipse2', 'useaxdata'});        uimenu(cmenu,'label','cut ellipse&3','Callback',{@PC.XYCluster, 'ellipse3', 'useaxdata'});        uimenu(cmenu,'label','cut ellipse&4','Callback',{@PC.XYCluster, 'ellipse4', 'useaxdata'});        uimenu(cmenu,'label','cut ellipse&5','Callback',{@PC.XYCluster, 'ellipse5', 'useaxdata'});        uimenu(cmenu,'label','cut ellipse&6','Callback',{@PC.XYCluster, 'ellipse6', 'useaxdata'});        uimenu(cmenu,'label','&Flip Line','Callback',{@PC.FlipLineCrit});               uimenu(cmenu,'label','No cutting','Callback',{@PC.XYCluster, 'nocut'});        uimenu(cmenu,'label','save (&Quick)','Callback',{@PC.XYCluster, 'quicksave'});        uimenu(cmenu,'label','Optimize','Callback',{@PC.XYCluster, 'optimize'});        uimenu(cmenu,'label','test line','Callback',{@PC.XYCluster, 'testline'});    elseif strcmp(type,'cellplot')        tm = uimenu(cmenu,'label','Ex0 P0:','Tag','Title');        mum = uimenu(cmenu,'label','Ex0 P0:','Tag','Title');        uimenu(cmenu,'label','spool','Callback',{@PC.SetCellFromLine, 1,  'spool'});        uimenu(cmenu,'label','spool auto','Callback',{@PC.SetCellFromLine, 1,  'spoolauto'});        uimenu(cmenu,'label','->AllVPcs','Callback',{@PC.SetCellFromLine, 1,  'getfullv'});        pm = uimenu(cmenu,'label','AllVPcs');        uimenu(pm,'label','RefCut','Callback',{@PC.SetCellFromLine, 1,  'refcut'});        uimenu(pm,'label','Apply Last','Callback',{@PC.SetCellFromLine, 1,  'lastcut'});        uimenu(pm,'label','Use Fullv','Callback',{@PC.SetCellFromLine, 1,  'forcefullv'});        uimenu(cmenu,'label','Find Duplicates','Callback',{@PC.SetCellFromLine, 1,  'findduplicate'});        pm = uimenu(cmenu,'label','Plots');        uimenu(pm,'label','All XY This Probe','Callback',{@PC.SetCellFromLine, 1,  'allxyprobe'});        uimenu(pm,'label','All XY This Expt','Callback',{@PC.SetCellFromLine, 1,  'allxyexpt'});        uimenu(pm,'label','All Spks This Probe','Callback',{@PC.SetCellFromLine, 1,  'allspksprobe'});        uimenu(pm,'label','All Spks This Expt','Callback',{@PC.SetCellFromLine, 1,  'allspksexpt'});        uimenu(pm,'label','All XY This Cell','Callback',{@PC.SetCellFromLine, 0,  'allxycell'});        uimenu(pm,'label','All Spks This Cell','Callback',{@PC.SetCellFromLine, 0,  'allspkscell'});        uimenu(pm,'label','Rates For This Cell','Callback',{@PC.SetCellFromLine, 0,  'cellrates'});        pm = uimenu(cmenu,'label','Other');        uimenu(pm,'label','Comment','Callback',{@PC.TagMenu, 'comment'});        uimenu(pm,'label','Flip Line Crit','Callback',{@PC.TagMenu, 'flipline'});        cells = unique(DATA.CellList(:));        cells = cells(cells > 0);        if max(cells) >= 20            nc = max(cells)+4;        else            nc = 20;        end        for k = 1:nc            c = uimenu(tm,'label',sprintf('Cell %d',k),'Callback',{@PC.SetCellFromLine, 0,  k});            d = uimenu(mum,'label',sprintf('Cell %d is MU',k),'Callback',{@PC.SetCellFromLine, 0,  k, 'muonly'});            if ismember(k,cells)                if isfield(DATA,'cellcolors') && k <= length(DATA.cellcolors)                set(c,'foregroundcolor',DATA.cellcolors{k});                else                set(c,'foregroundcolor','r');                end                set(c,'checked','on');            end            if k == 1                set(c,'separator','on');            end        end    else        uimenu(cmenu,'label',sprintf('Spike%d',j-1),'foregroundcolor',DATA.colors{j});    for k = 1:50        c = uimenu(cmenu,'label',sprintf('Cell %d',k),'Callback',{@PC.SetCellFromLine, j-1,  k});        if k == 1            set(c,'separator','on');        end    end    end    