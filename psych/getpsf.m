function results = getpsf(probitin, varargin)
%
%function results = getpsf(probitin, varargin)
%fits a culmulative Gaussian to binomial data. The data is in a structure
%array 'probit'. Each element must contain
%probit.x  the metameter
%probit.n  the number of stimulus presentations
%probit.resp the number of responses/detections
%
%An addition field, probit.expno can be used to define sepatate data
%sets.  getpsf(probit,'expno',n) will fit only the data 
%for which probit.expno == n

%initialize optional variables     
fitopt.twomeans = 0;
     fitopt.expnos(1) = 1;
     fitopt.expnos(2) = 2;
     fitopt.twosd = 0;
     fitopt.twofits = 0;
     showinit = 0;

 %make a copy of the input data.    
     probit = probitin;
     j = 1;
     while(j < nargin)
        if(strcmpi(varargin{j},'showinit'))
           showinit = 1;
       elseif(strcmpi(varargin{j},'expno'))
%if a single expno is given, extract this data from the full set by
%matching the expno field in the input data structure.
            j = j+1;
           expno = varargin{j};
           idx = find([probitin.expno] == expno);
           probit = probitin(idx);
        elseif(strncmpi(varargin{j},'twomeans',6))
            fitopt.twomeans = 1;
        elseif(strncmpi(varargin{j},'twosd',5))
            fitopt.twosd = 1;
        elseif(strncmpi(varargin{j},'twofits',6))
            fitopt.twofits = 1;
        end
     j = j+1;
     end
     

p = [probit.resp]./[probit.n];
options = optimset([]);
%make an initial guess. x(1) is the mean, x(2) is the sd.
x(1) = mean(exp(-(p -0.5).^2 ./0.1) .* [probit.x]);
     nsd = erfinv((p *2) -1);
     idx = find(nsd > 2);
     nsd(idx) = 2; 
     idx = find(nsd < -2);
     nsd(idx) = -2; 
     x(2) = mean(nsd .* ([probit.x] - x(1)));
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
    [results.fit, results.llike] = fminsearch(@llike, x, options, probit, fitopt);
    results.data = probit;
    results.options = fitopt;
    
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
else
    y = ([probit.x] - mean)./sd;
end
p = (erf(y) +1)/2;
q = 1-p;
result = sum((([probit.n] - [probit.resp]) .* log(q)) + (([probit.resp] .* log(p))));
result = -result;
