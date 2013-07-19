function corrs = Threshcorr(varargin)

nsmp = 500;
nloops = 500;
corr = 0.2;
th= 0.9;
j = 1;
while j <= nargin
    if strncmpi(varargin{j},'loops',3)
        j = j+1;
        nloops = varargin{j};
    elseif strncmpi(varargin{j},'th',2)
        j = j+1;
        th = varargin{j};
    elseif strncmpi(varargin{j},'corr',3)
        j = j+1;
        corr = varargin{j};
    elseif strncmpi(varargin{j},'nsmp',3)
        j = j+1;
        nsmp = varargin{j};
    end
    j = j+1;
end

crnd = randn(nloops,nsmp) .* corr;
rndb = ((1-corr) * randn(nloops,nsmp))+crnd;
rnda = ((1-corr) * randn(nloops,nsmp))+crnd;
prndb = (rndb > th) .* 1.0;
prnda = (rnda > th) .* 1.0;
for j = 1:nloops
    xc = corrcoef(prndb(j,:),prnda(j,:));
    corrs(j) = xc(1,2);
end
xc = corrcoef(sum(prndb),sum(prnda));
countcorr = xc(1,2);
xc = corrcoef(rndb,rnda);
membcorr = xc(1,2)
xc = corrcoef(prndb,prnda);
spkcorr = xc(1,2)
