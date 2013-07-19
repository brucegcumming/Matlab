function [dprime, details] = ExptCellQuality(Expt, varargin)

verbose = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'verbose',4)
        verbose = 1;
    end
    j = j+1;
end
errs = 0;
ps = Expt.probes;
blks = [Expt.Header.BlockStart Expt.Trials(end).Trial];
for j = 1:length(blks)-1;
    Trials = blks(j):blks(j+1);
    id = find(ismember([Expt.Trials.Trial],Trials));
    if length(id)
    p = median(Expt.probes(id));
    C = Expt.Header.Clusters{j}{1,p};
    if isfield(C,'dprime')
        dp(j) = C.dprime;
    else  %Shouldn't happen
        if verbose
            if isfield(Expt.Header,'Filename')
            fprintf('%s Block %d Missing cluster for probe %d\n',Expt.Header.Filename,j,p);
            else
            fprintf('%s.%s Block %d Missing cluster for probe %d\n',Expt.Header.Name,Expt.Header.expname,j,p);
            end
        end
        dp(j) = NaN;
        errs = errs+1;
    end
    n(j) = length(id);
    else
        dp(j) = 0;
        n(j) = 0;
    end
end
dprime = WeightedSum(dp,n);
details.dprimes = dp;
details.n = n;
details.errs = errs;