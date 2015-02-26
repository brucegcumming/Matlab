function SetSubplots(nr, nc, plots, flag, varargin)

%SetSubplots(nrow, ncol, plots, flag)
%iterates through a set of subplots to ensure that they all have
%the same value some property, detemined by flag
%SetSubplots(3,3,[1:3],'caxis')
%    makes sure that the firstthree plots all use the same colorscale
% if flag is a cell string, then the cells are passed as arguments to set(gca,
% for each subbplot


if iscell(flag)
     for j = plots;
        subplot(nr, nc, j);
        set(gca,flag{:});
     end
elseif strncmpi(flag,'pass',3)  %% pass on varargin
     for j = plots;
        subplot(nr, nc, j);
        set(gca,varargin{:});
     end
elseif strncmpi(flag,'caxis',3)  %% make caxis the same for all subplots
    crange = [NaN NaN];
    for j = plots;
        subplot(nr, nc, j);
        cx = caxis;
        crange(1) = min([cx(1) crange(1)]);
        crange(2) = max([cx(2) crange(2)]);
    end
    for j = plots;
        subplot(nr, nc, j);
        caxis(crange);
    end
end