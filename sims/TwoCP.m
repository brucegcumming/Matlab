function TwoCP(N)
%TwoCP(N) build a population of N nuerons that allows CP for two features
%If pooled accordingly
%
%
%    NL    NR     FL    FR
% NL 0.2   0.1    0.1   0.0
% NR 0.1   0.2   0.0    0.1
% FL 0.1   0.0   0.2   0.1
% FR 0.0   0.1   0.1   0.2

ra=0.2
rb=0.1;
ntrials = 100;
ai = 1:N/4;
bi = 1+N/4:N/2;
ci = 1+N/2:3*N/4;
di = 1+0.75*N:N;

R = ones(N,N) .* rb;
R(ai,ai) = ra;
R(bi,bi) = ra;
R(ci,ci) = ra;
R(di,di) = ra;
R(ai,di) = 0;
R(bi,ci) = 0;
R(ci,bi) = 0;
R(di,ai) = 0;

C = MatrixCorrCounts(R,100, size(R,1));
dvm = mean(C([ai ci],:)) - mean(C([bi di],:));
dvd = mean(C([ai bi],:)) - mean(C([ci di],:));
dpref = find(dvd > 0);
dnull = find(dvd <= 0);
mpref = find(dvm > 0);
mnull = find(dvm <= 0);

%dpref is a list of trials where a NEar choice is made
for j = 1:length(ai)
    cp(ai(j),1) = CalcCP(C(ai(j),dpref),C(ai(j),dnull));
    cp(ai(j),2) = CalcCP(C(ai(j),mpref),C(ai(j),mnull));
    cp(bi(j),1) = CalcCP(C(bi(j),dpref),C(bi(j),dnull));
    cp(bi(j),2) = CalcCP(C(bi(j),mnull),C(bi(j),mpref));
    cp(ci(j),1) = CalcCP(C(ci(j),dnull),C(ci(j),dpref));
    cp(ci(j),2) = CalcCP(C(ci(j),mpref),C(ci(j),mnull));
    cp(di(j),1) = CalcCP(C(di(j),dnull),C(di(j),dpref));
    cp(di(j),2) = CalcCP(C(di(j),mnull),C(di(j),mpref));
end
mean(cp(ai,:))
mean(cp(bi,:))
mean(cp(ci,:))
mean(cp(di,:))
