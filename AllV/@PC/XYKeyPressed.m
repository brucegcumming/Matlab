  function XYKeyPressed(src, ks, fcn)      if strcmp(ks.Key, 'control')          setappdata(0,'control_is_down',1);          return;      elseif strcmp(ks.Key,{'alt'})          setappdata(0,'alt_is_down',1);          return;      elseif sum(strcmp(ks.Key,{'shift' 'control' 'alt'}))          return;      end      DATA = GetDataFromFig(src);tag = get(src,'tag');bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});if sum(strcmp(ks.Key,{'uparrow' 'downarrow'})) & strcmp(ks.Modifier,'shift')    if DATA.NewCut.saved == 0        DATA = PC.SaveCluster(DATA, DATA.currentpoint,1);    end    DATA = PC.AddSelectedCells(DATA);endif strcmp(tag,DATA.tag.allxy) && DATA.allclustering    if strcmp(ks.Key,'l')        PC.XYCluster(src, [],'lines');    elseif strcmp(ks.Key,'1')        PC.XYCluster(src, [],'ellipse1');    elseif strcmp(ks.Key,'2')        PC.XYCluster(src, [],'ellipse2');    elseif strcmp(ks.Key,'3')        PC.XYCluster(src, [],'ellipse3');    elseif strcmp(ks.Key,'4')        PC.XYCluster(src, [],'ellipse4');    elseif strcmp(ks.Key,'5')        PC.XYCluster(src, [],'ellipse5');    elseif strcmp(ks.Key,'6')        PC.XYCluster(src, [],'ellipse6');    elseif strcmp(ks.Key,'f')        DATA = PC.FlipLineCrit(DATA);    elseif strcmp(ks.Key,'q')        PC.XYCluster(src, [],'quicksave');    elseif strcmp(ks.Key,'n')        PC.XYCluster(src, [],'nocut');    endelse    if strcmp(ks.Key,'f')        DATA = PC.FlipLineCrit(DATA);    elseif strcmp(ks.Key,'d')        DATA.plot.density = ~DATA.plot.density;        PC.ReplotXY(DATA,[],DATA.currentpoint(1),DATA.currentpoint(2),DATA.currentcluster);        set(DATA.toplevel,'UserData',DATA);    elseif strcmp(ks.Key,'q')               PC.XYCluster(src, [],'quicksave');    elseif strcmp(ks.Key,'add') | ks.Character == '+' %characte can be empty        DATA = PC.AddSelectedCells(DATA);        endendif strcmp(ks.Key,'rightarrow')    DATA = PC.NextButton(src, ks, 'r');elseif strcmp(ks.Key,'leftarrow')    DATA = PC.NextButton(src, ks, 'l');elseif strcmp(ks.Key,'downarrow')    DATA =  PC.NextButton(src, ks, 'd');elseif strcmp(ks.Key,'uparrow')    DATA = PC.NextButton(src, ks, 'u');elseif strcmp(ks.Key,'s')    PC.SpoolSpikes(DATA, DATA.currentpoint);endset(DATA.toplevel,'UserData',DATA);