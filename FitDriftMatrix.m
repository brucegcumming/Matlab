function [P, details] = FitDriftMatrix(M, varargin)

maxiter = 20000;
quickfirst = 1;
%initial estimate is just measured steps between adjacent cells
for j = 2:size(M,1)
    dvals(j-1) = M(j,j-1);
end
%MakeGuess(M, dvals);

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'guess',5)
        j = j+1;
        guess = varargin{j};
        dvals(1:length(guess)) = guess;
        if sum(abs(dvals) > 0) > length(dvals)/2 %Not sparse
            quickfirst = 0;
        end
    elseif strncmpi(varargin{j},'maxiter',5)
        j = j+1;
        maxiter = varargin{j};
    end
    j  = j+1;
end

dvals(isnan(dvals)) = 0;
setappdata(gcf,'fittedparams',[]);
[ssq, guess] = TryMatrix(dvals,M,1:length(dvals));
options = optimset('MaxFunEvals',100000,'maxiter',maxiter,'display','off');
if quickfirst
    [a,b] = find(abs(dvals) > 0.8);
    if ~isempty(a)
        try
        [ssq,P, errs] = TryMatrix(a,M,b);
        [fittedparams,fval,exitflag, output] = fminsearch(@TryMatrix,a,options,M,b);
        dvals(b) = fittedparams;
        end
    end
end
b = 1:length(dvals);
[fittedparams,fval,exitflag, output] = fminsearch(@TryMatrix,dvals,options,M,b);
[ssq,P, errs] = TryMatrix(fittedparams,M,b);
details.jumps = fittedparams;
details.ssq = ssq;
details.errs = errs;
details.guess = dvals;
details.exitflag = exitflag;
details.output = output;
imagesc(errs);

function [ssq, P, errs] = TryMatrix(x, M, xi)

X = zeros(1,size(M,1)-1);
X(xi) = x;
n = length(X);
for j = 1:n
    P(j+1,j) = X(j);
end
P(n+1,n+1) = 0;
for j = 1:n
    for k = j:-1:1
        P(j+1,k) = sum(X(k:j));
        errs(j,k) = P(j+1,k)-M(j+1,k);
    end
end
params = getappdata(gcf,'fittedparams');
errs(isnan(errs)) = 0;
ssq = sum(errs(:).^2);
%params = cat(1,params, [ssq x]);
params = cat(1,params, [ssq sum(x)]);
if mod(size(params,1),1000) == 0
    fprintf('SSQ %.1f, %d evals %s\n',ssq, size(params,1),datestr(now));
end
setappdata(gcf,'fittedparams',params);


function MakeGuess(M, x)

for j = 1:length(x)
    P(j,j) = x(j);
end
for j = 1:size(P,1)
    for k = j-1:-1:1
        P(j,k) = sum(x(k:j));
        errs(j,k) = P(j,k)-M(j+1,k);
    end
end
errs(isnan(errs)) = 0;
for j = 2:size(P,1)
    rowerr(j-1) = mean(errs(j,1:j-1));
    colerr(j-1) = mean(errs(j-1:end,j-1));
end
[a,b] = max(rowerr);
[c,d] = max(colerr);

