function HitXcorrPt(src, b, id, j)
DATA = GetDataFromFig(src);
GetFigure('Counts','front');
X = DATA.xcorrs(id);
fprintf('Trial id %d\n',X.trialids(j));

