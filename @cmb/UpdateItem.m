function UpdateItem(a,b, varargin)
%handle change in checkbox for a single item

DATA = GetDataFromFig(a);
calltag = get(a,'Tag');
val = get(a,'value');
if strcmp(calltag,'UseLastLayout')
    DATA.optoins.uselastlayout = val;
elseif strcmp(calltag,'UseLastConfig')
    DATA.options.uselastconfig = val;
elseif strcmp(calltag,'CalcCCF')
    e = DATA.currentexpt(1);
    cp = get(a,'value');
    DATA = cmb.BuildAllTimes(DATA, 'expts', e);
    A = cmb.SpkCache(DATA,e,DATA.probe,'get');
    aid = find(A.codes ==1);
    B = cmb.SpkCache(DATA,e,cp,'get');
    bid = find(B.codes ==1);
    [xc, details] = xcorrtimes(A.times(aid)./10000,B.times(bid)./10000);
    GetFigure(DATA.tag.dataplot);
    hold off; 
    plot(details.xpts,xc);
    title(sprintf('Probe %d->%d. %d and %d Spikes',DATA.probe,cp,length(aid),length(bid)));
    xlabel('sec');
    ylabel('coincidences');
end
SetData(DATA);