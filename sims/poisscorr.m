function corrs = poisscorr(varargin)

nsmp = 500;
nloops = 500;
j = 1;
while j <= nargin
    if strncmpi(varargin{j},'loops',3)
        j = j+1;
        nloops = varargin{j};
    end
    j = j+1;
end
prndb = random('poiss',0.1,nloops,nsmp);
prnda = random('poiss',0.1,nloops,nsmp);
for j = 1:nloops
    xc = corrcoef(prndb(j,:),prnda(j,:));
    corrs(j) = xc(1,2);
end