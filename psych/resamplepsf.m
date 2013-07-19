function fit = resamplepdf(fit, nresample,varargin)

if nresample == 0
    return;
end

alpha = 0.05;
args = {};
j = 1;
while j < nargin -1
    if(strncmpi(varargin{j},'alpha',5))
        j = j+1;
        alpha = varargin{j};
    else
        args{j} = varargin{j};
    end
    j = j+1;
end


for j = 1:length(fit.data)
    fit.data(j).reresp = binornd(fit.data(j).n,fit.data(j).p,nresample,1);
end

redata = fit.data;

for res = 1:nresample;
    for j = 1:length(fit.data)
        redata(j).resp = fit.data(j).reresp(res);    
        redata(j).p = redata(j).resp / redata(j).n;    
    end
    refit = fitpsf(redata,args{:});
    refits(res,:) = refit.fit;
end
fit.refits = refits;


ll = round(nresample * alpha/2);
ul = nresample - ll;
ms = sort(refits(:,1));
fit.meanci = [ms(ll) ms(ul)];
ms = sort(refits(:,2));
fit.sdci = [ms(ll) ms(ul)];
fit.alpha = alpha;
