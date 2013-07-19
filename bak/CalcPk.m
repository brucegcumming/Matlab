function [kernel, details] = CalcPk(Expt,varargin)
% Calculate a psychophysical kernel from revcor Expt data
dcval = 0;
plottype = 0;
details = [];

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'dc',2)
        j = j+1;
        dcval= varargin{j};
    elseif strncmpi(varargin{j},'plot',4)
        plottype = 1;
        
    end
    j = j+1;
end

type = Expt.Stimvals.e2;
emean = mean([Expt.Trials.ori]);
details.signal = sort(unique([Expt.Trials.ori]));
details.signaln = hist([Expt.Trials.ori],details.signal);
for j = 1:length(dcval)
dcvals = [Expt.Trials.Dc] .* sign([Expt.Trials.ori] - emean); 
uid = find(dcvals == dcval(j) & [Expt.Trials.RespDir] > 0);
uframes = [Expt.Trials(uid).(type)];
did = find(dcvals == dcval(j) & [Expt.Trials.RespDir] < 0);
dframes = [Expt.Trials(did).(type)];

xv = unique([uframes(:); dframes(:)]);
un(:,j) = hist(uframes(:),xv);
dn(:,j) = hist(dframes(:),xv);
uk = un(:,j)'./length(uid);
dk = hist(dframes(:),xv)./length(did);
nstim(j) = length(uid)+length(did);
kernels(j,:) = (dk-uk) .* nstim(j);
end
kernel = sum(kernels,1)./sum(nstim);

if plottype
    plot(xv,kernel);
end
details.xv = xv;
details.n = [sum(un,2) sum(dn,2)];