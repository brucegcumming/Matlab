function results = fitpsf(probitin, varargin)
%
%function results = fitpsf(probitin, varargin)
%fits a culmulative Gaussian to binomial data. The data is in a structure
%array 'probit'. Each element must contain
%probit.x  the metameter
%probit.n  the number of stimulus presentations
%probit.resp the number of responses/detections
% the number of array elements = number of conditions (x values).
%
%An additional field, probit.expno can be used to define sepatate data
%sets.  fitpsf(probit,'expno',n) will fit only the data 
%for which probit.expno == n
%
%results.fit contains the result [mean sd] of fitted Gaussian
%fitpsf(psf, 'eval', x)   evaluates a returned fit (psf) at the points x
% fitpsf(psf, 'showfit')  plots data and fit
%


%initialize optional variables     
fitopt.twomeans = 0;
     fitopt.expnos(1) = 1;
     fitopt.expnos(2) = 2;
     fitopt.twosd = 0;
     fitopt.twofits = 0;
     fitopt.getsderr = 0;
     fitopt.findsdlim = 0;
     showinit = 0;
     coff = 0;

 %make a copy of the input data.    
    if sum([probitin.n]) == 0;
        results = [];
        return;
    end
     if length(probitin) ==1 & isfield(probitin,'resp')
         for j = 1:length(probitin.x)
             probit(j).x = probitin.x(j);
             probit(j).n = probitin.n(j);
             if ~isfield(probitin,'p')
                 probit(j).p = probitin.resp(j)./probitin.n(j);
             else
                 probit(j).p = probitin.p(j);
             end
             probit(j).resp = probitin.resp(j);
             if isfield(probitin,'expno')
             probit(j).expno = probitin.expno(j);
             else
                 probit(j).expno = 1;
             end
                 
         end
     else
         probit = probitin;
     end
         
     j = 1;
     while(j < nargin)
        if(strcmpi(varargin{j},'showinit'))
           showinit = 1;
       elseif(strcmpi(varargin{j},'expno'))
%if a single expno is given, extract this data from the full set by
%matching the expno field in the input data structure.
            j = j+1;
           expno = varargin{j};
           idx = find([probit.expno] == expno);
           probit = probit(idx);
        elseif(strncmpi(varargin{j},'coff',4))
            j = j+1;
            coff = varargin{j};
        elseif(strncmpi(varargin{j},'twomeans',6))
            fitopt.twomeans = 1;
        elseif(strncmpi(varargin{j},'eval',4))
            j = j+1;
            x = varargin{j};
            fit = probit;
            results = 0.5 + erf((x - fit.fit(1))/(fit.fit(2) * sqrt(2)))/2;
            return;
        elseif(strncmpi(varargin{j},'twosd',5))
            fitopt.twosd = 1;
        elseif(strncmpi(varargin{j},'sderr',5))
            fitopt.getsderr = 1;
        elseif(strncmpi(varargin{j},'twosd',5))
            fitopt.twosd = 1;
        elseif(strncmpi(varargin{j},'showfit',5))
        if length(varargin) > j && isstruct(varargin{j+1})
            j = j+1;
            infit = varargin{j};
            probit = infit.data;
        else
            infit = fitpsf(probit,varargin{1:j-1});
        end
        if fitopt.twomeans
            colors = mycolors;
            id = find([infit.data.expno] == 1);
            [results.fith(1), results.fitx, results.fity] = plotpsych(infit.data(id), infit.fit(1),infit.fit(2),'color',colors{1+coff},varargin{j+1:end});
            id = find([infit.data.expno] == 2);
            [results.fith(2) , results.fitx, results.fity]= plotpsych(infit.data(id), infit.fit(3),infit.fit(2),'color',colors{2+coff},varargin{j+1:end});
        else
            [results.fith, results.fitx, results.fity] = plotpsych(infit.data, infit.fit(1),infit.fit(2),varargin{j+1:end});
        end
        elseif(strncmpi(varargin{j},'xmax',4))
            j = j+1;
            fitopt.xmax = varargin{j};
            idx = find([probit.x] < fitopt.xmax);
            probit = probit(idx);
        elseif(strncmpi(varargin{j},'nmin',4))
            j = j+1;
            fitopt.nmin = varargin{j};
            idx = find([probit.n] >= fitopt.nmin);
            probit = probit(idx);
        elseif(strncmpi(varargin{j},'orbw',4))
        elseif(strncmpi(varargin{j},'xmin',4))
            j = j+1;
            fitopt.xmin = varargin{j};
            idx = find([probit.x] > fitopt.xmin);
            probit = probit(idx);
        elseif(strncmpi(varargin{j},'twofits',6))
            fitopt.twofits = 1;
        end
     j = j+1;
     end

% get rid of empty data 
id = find([probit.n] > 0);
probit = probit(id);

%calculate p
p = [probit.resp]./[probit.n];
for j = 1:length(probit)
    probit(j).p = p(j);
end
options = optimset([]);

if length(p) < 2
    results.fit(1) = NaN;
    results.fit(2) = NaN;
    results.exit = 0;
    return;
end
%make an initial guess. x(1) is the mean, x(2) is the sd.
w = exp(-(p -0.5).^2 ./0.1);
x(1) = sum(w .* [probit.x] .* [probit.n])/sum(w .*[probit.n]);
nsd = erfinv((p *2) -1);
idx = find(nsd > 2);
nsd(idx) = 2; 
idx = find(nsd < -2);
nsd(idx) = -2; 
if(max(nsd) > 0.1)
    idx = find(abs(nsd) > 0.1);
elseif(max(nsd) > 0.02)
    idx = find(abs(nsd) > 0.02);
else 
    idx = 1:length(nsd);
end
nsd(find(nsd == 0)) = mean(nsd)/10;
x(2) = sum((([probit(idx).x] - x(1))./nsd(idx)) .* [probit(idx).n])./sum([probit(idx).n]);
if x(2) > diff(minmax([probit.x]))
    x(2) = diff(minmax([probit.x]));
end

results.preg = polyfit([probit.x],nsd,1);
x(2) = 1/results.preg(1);
a = llike(x, probit, fitopt);
b = llike([x(1) -x(2)], probit, fitopt);
if(b < a)
    x(2) = -x(2);
end

%here is the initial guess. Calclate its log likelihood.
results.initial = x;
results.initlike = llike(x,probit,fitopt);
if(showinit)
  x
end
    
    
% if twomeans, twosd, or two fit options are set, First fit all the
% data together, then refit with two curves, and compare likelihoods
    if(fitopt.twomeans > 0)
        fitopt.twomeans = 0;
        [results.allfit, results.allllike] = fminsearch(@llike, x, options, probit, fitopt);
        x(3) = x(1);  
       fitopt.twomeans = 1;
    end
    
    if(fitopt.twosd > 0)
        fitopt.twosd = 0;
        [results.allfit, results.allllike] = fminsearch(@llike, x, options, probit, fitopt);
        x(3) = x(2);  
       fitopt.twosd = 1;
    end    

    
    if(fitopt.twofits > 0)
        fitopt.twofits = 0;
        [results.allfit, results.allllike] = fminsearch(@llike, x, options, probit, fitopt);
        x(3) = x(1);
        x(4) = x(2);
       fitopt.twofits = 1;
    end
    [results.fit, results.llike, results.exit] = fminsearch(@llike, x, options, probit, fitopt);

    results.data = probit;
    results.options = fitopt;

    if(fitopt.getsderr)
      fitopt.findsdlim = 1;
      fitopt.llike = results.llike;
      fitopt.mean = results.fit(1);
      x(2) = results.fit(2) * 1.2;
      [tmp, results.llikesd] = fminsearch(@llike, x, options, probit, ...
					  fitopt);
      results.sdlim(1) = tmp(2);
      x(2) = results.fit(2) * 0.8;
      [tmp, results.llikesd] = fminsearch(@llike, x, options, probit, fitopt);
      results.sdlim(2) = tmp(2);
    end
    
    
% comparing two log-likelihoods, where one has an extra paramenter, 2* the
% difference in log likelihood has a chisq distribution with 1 d.f.
    if(fitopt.twomeans > 0 | fitopt.twofits > 0 | fitopt.twosd > 0)
        results.chi = 2 * (results.allllike - results.llike);
        results.pval = 1 - chi2cdf(results.chi,1);
 %if there are two fits, there are two extra d.f. I am GUESSING that this
 %means that the chisq test should be for two degreees of freedom, but I
 %have not checked this...............
        if(fitopt.twofits > 0)
          results.pval = 1 - chi2cdf(results.chi,2);          
        end
    end

function result = llike(params, probit,fitopt);

mean = params(1);
sd = params(2);

% if params > 2, this means two sets are being fit, with the same SD
% but different means. The two datasets are
if length(params) > 2  & fitopt.twomeans > 0
     idx = find([probit.expno] == fitopt.expnos(1));
     means(idx) = deal(params(1));
     idx = find([probit.expno] == fitopt.expnos(2));
     means(idx) = deal(params(3));
     y = ([probit.x] - means) ./ sd;
 elseif  length(params) > 2  & fitopt.twosd > 0
     idx = find([probit.expno] == fitopt.expnos(1));
     sds(idx) = deal(params(2));
     idx = find([probit.expno] == fitopt.expnos(2));
     sds(idx) = deal(params(3));
     y = ([probit.x] - mean) ./ sds; 
 elseif  length(params) > 3  & fitopt.twofits > 0
     idx = find([probit.expno] == fitopt.expnos(1));
     sds(idx) = deal(params(2));
     means(idx) = deal(params(1));
     idx = find([probit.expno] == fitopt.expnos(2));
     sds(idx) = deal(params(4));
     means(idx) = deal(params(3));
     y = ([probit.x] - means) ./ sds;     
elseif fitopt.findsdlim
    y = ([probit.x] - fitopt.mean)./sd;
else
    y = ([probit.x] - mean)./sd;
end

p = (erf(y/sqrt(2)) +1)/2;
q = 1-p;

%can't have log(0).
p(find(p<1e-20)) = 1e-20;
q(find(q<1e-20)) = 1e-20;

result = sum((([probit.n] - [probit.resp]) .* log(q)) + (([probit.resp] .* log(p))));
result = -result;

if fitopt.findsdlim
   result = abs(result - (fitopt.llike+3.84));
end
   


function [fitline, x, y] = plotpsych(data, mean, sd, varargin)

color = 'r';
shown = 0;
j = 1;
while(j < nargin-2)
  if(strncmpi(varargin{j},'color',5))
    j = j+1;
    color = varargin{j};
  elseif strncmpi(varargin{j},'shown',4)
    shown = 1;
  end
  j = j+1;
end

if length(data) == 0
    fitline = 0;
    x = [];
    y = [];
    return;
end
   
h = errorbar([data.x],[data.p],sqrt([data.p] .* (1 - [data.p])./[data.n]),'o');
set(h,'color',color);
set(h,'MarkerFaceColor',color,'linewidth',2);
if length(data) < 2
    fitline = h;
    x = [];
    y = [];
    return;
end
hold on;
step = (max([data.x]) - min([data.x]))/100;
x = min([data.x]):step:max([data.x]);
y = 0.5 + erf((x - mean)/(sd * sqrt(2)))/2;
fitline = plot(x,y,'color',color,'linewidth',2);

if(shown)
for j = 1:length(data)
  text(data(j).x+2*step,data(j).p,sprintf('%d',data(j).n),'color',color);
end
end
