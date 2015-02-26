function [x, details] = CalcIsolation(xy, idlist, clnum, varargin)
%[x, details] = CalcIsolation(pts, idlist, clnum...   calculate islolation metrics from data points
%pts is an n (points) x m (dimension of space) matrix
%idlist is a vector of n classfifications
%without using GM fits
plottype = 'none';
parentfig = 0;
id = find(idlist == clnum);
nid = find(idlist ~= clnum);
details = [];
x(1:2) = NaN;  %default. Be sure length matches x in good returns
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'plothist',5)
        plottype = 'hist';
        if length(varargin) > j && isfigure(varargin{j+1})
            j = j+1;
            parentfig = varargin{j};
        end
    end
    j =j+1;
end

if length(id) < 3  || length(nid) < 3
    return;
end

%easy to get ill conditioned matrices if variables are highly correlated
wastate = warning('Off','MATLAB:illConditionedMatrix');  
wbstate = warning('Off','MATLAB:nearlySingularMatrix');  
%distance of cell points from MU dist
if length(nid) > size(xy,2) && length(id) > size(xy,2)
    dm = sqrt(mahal(xy(id,:),xy(nid,:)));
    %distance of mu points from MU dist
    dmm = sqrt(mahal(xy(nid,:),xy(nid,:)));
    %distance of MU points from SU dist
    ds = sqrt(mahal(xy(nid,:),xy(id,:)));
    %distance of cell points from SU dist
    dss = sqrt(mahal(xy(id,:),xy(id,:)));
    maxd = max([max(ds) max(dm)]);
    bins = linspace(0,maxd);
    x(1) = prctile(ds,1);
    x(1) = AllV.CalcDprime(ds,dss);
    x(2) = AllV.CalcDprime(dm,dmm);
else
    maxd = NaN;
end

warning(wastate);
warning(wbstate);


if strcmp(plottype,'hist') && ~isnan(maxd)
    if double(parentfig) > 0
        GetFigure('MahalHist','parent',parentfig);
    else
        GetFigure('MahalHist');
    end
      
    a = hist(dss,bins);
    b = hist(ds,bins);
    c = hist(dmm,bins);
    d = hist(dm,bins);
    n(1) = max(a+b);
    n(2) = max(c+d);
    hold off;
    plot(bins,a./n(1));
    hold on;
    plot(bins,b./n(1),'g-');
    plot(bins,(a+b)./n(1),'r');
    plot(bins,c./n(2),'--');
    plot(bins,(c+d)./n(2),'r--');
    plot(bins,d./n(2),'g--');
    title(sprintf('Dprimes for SU %.2f,MU %.2f',x));
end

