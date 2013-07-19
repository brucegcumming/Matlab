function [disps, resp, aresp] = TwoSubunits(nreps, varargin)

p = [1 1];
plottype = 0;
dx = 0.5;
sd = 0.5;
baseline = 1;
coeff = [1 -1];
figlabel = {'TwoSubUnitsa', 'TwoSubUnitsb'};

j = 1;
while j < nargin
    if strncmpi(varargin{j},'base',2)
        j = j+1;
        baseline = varargin{j};
    elseif strncmpi(varargin{j},'dx',2)
        j = j+1;
        dx = varargin{j};
    elseif strncmpi(varargin{j},'coeffs',2)
        j = j+1;
        coeff = varargin{j};
    elseif strncmpi(varargin{j},'exp',3)
        j = j+1;
        p = varargin{j};
    elseif strncmpi(varargin{j},'sd',2)
        j = j+1;
        sd = varargin{j};
    elseif strncmpi(varargin{j},'show',3)
        plottype = 1;
    end
    j = j+1;
end

dispset = -1:0.1:1;
[disps, dresps, aresps, anti{1}, rnd] = rls(nreps,'sd',sd,'disps',dispset);
%[disps, dresps, bresps, rnd] = rls(nreps,'rnd',rnd,'dx',0.5);
[disps, dresps, bresps, anti{2}, rnd] = rls(nreps,'rnd',rnd,'dx',dx,'sd',sd,'disps',dispset);

aresps = aresps.^p(1);
bresps = bresps.^p(2);
anti{1} = anti{1}.^p(1);
anti{2} = anti{2}.^p(2);

if(baseline <= 0)
    tresp = aresps+coeff(2) * bresps;
    taresp = anti{1} +coeff(2) * anti{2};
else
    tresp = max(baseline+aresps+coeff(2) * bresps,zeros(size(bresps)));
    taresp = max(baseline+anti{1}+coeff(2) * anti{2},zeros(size(bresps)));
end
resp = mean(tresp,2);
aresp = mean(taresp,2);


if plottype == 1
    GetFigure(figlabel{1});
    hold off;
    sds = std(aresps');
    %errorbar(disps,mean(aresps,2),sds);
    plot(disps,mean(aresps,2));
    hold on;
    sds = std(anti{1}');
%    errorbar(disps,mean(anti{1},2),sds,':');
    plot(disps,mean(anti{1},2),':');
    plot(disps,mean(bresps,2),'r');
    plot(disps,mean(anti{2},2),'r:');
    GetFigure(figlabel{2});
    hold off;
    if(baseline <= 0)
        errorbar(disps,resp,std(tresp'));
        hold on;
        errorbar(disps,aresp,std(taresp'),'r');
    else
        plot(disps,resp);
        hold on;
        plot(disps,aresp,':');
    end
    
end