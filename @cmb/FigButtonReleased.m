function FigButtonReleased(src, data)
global mousept;
mousept.mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
mousept.down = 0;
mousept.lasth = 0;
DATA = GetDataFromFig(src);
cmb.ExcludeTrials(DATA, mousept);

