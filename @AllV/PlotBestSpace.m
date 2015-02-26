function DATA = PlotBestSpace(DATA, varargin)
%Call PlotOneXY with space that shows bestspace

C = DATA.cluster
if ~isfield(C,'bestisolation')
    if C.auto && strcmp(C.autocutmode,'ecker')
        cl = DATA.currentcluster+1;
        if size(C.eckercluster.bestspace,1) < cl
            return;
        end
        space = C.eckercluster.bestspace(cl,:);
    end
else
    space = C.bestisolation.space;
end
AllV.SetFigure(DATA.tag.tmplscore, DATA);
if space(1) == 3 %template space
    xname = DATA.TemplateLabels{space(2)};
    yname = DATA.TemplateLabels{space(3)};
elseif space(1) == 1
    xname = sprintf('PC%d',space(2));
    yname = sprintf('PC%d',space(3));
end
if exist('xname')
    DATA.xyplot.xy = {xname yname};
    xy = AllV.PlotOneXY(DATA,{ xname yname});
    SetData(DATA);
    [E, score] = FindEllipse(xy, DATA.clst,'cluster',DATA.currentcluster);
    DrawEllipse(E,'add');
end