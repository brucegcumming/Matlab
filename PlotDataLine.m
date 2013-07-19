function allh = PlotDataLine(x,y,sds, n, id, varargin)
color = 'r';
symb = 'o';
linestyle = '-';
fillsymbols = 1;
if isempty(id)
    id = find(~isnan(x));
end
j = 1;
while j <=length(varargin)
    if strncmpi(varargin{j},'color',5)
        j = j+1;
        color = varargin{j};
    end
    j = j+1;
end
er = errorbar(x(id),y(id),sds(id) ./ sqrt(max(n(id),1)),symb,'linestyle','-');
hold on;
set(er,'color',color,'linestyle',linestyle);
h = plot(x(id),y(id),symb,'color',color);
if fillsymbols
    set(h,'MarkerFaceColor',color);
end
allh = [h er];
