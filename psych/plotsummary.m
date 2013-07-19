function plotsummary(pdata, varargin)

colors = mycolors;

j = 1;
narg = 1;
plotarg = {};
nmin = 0;
while j < nargin
    if strncmpi(varargin{j},'nmin',4)
        j = j + 1;
        nmin = varargin{j};
        plotarg{narg} = 'nmin';
        plotarg{narg+1} = nmin;
        narg = narg +2;
    elseif strncmpi(varargin{j},'xmin',4)
        j = j + 1;
        xmin = varargin{j};
        plotarg{narg} = 'xmin';
        plotarg{narg+1} = xmin;
        narg = narg +2;
    elseif strncmpi(varargin{j},'sds',3)
        j = j + 1;
        sdlist = varargin{j};
    elseif strncmpi(varargin{j},'mindiff',4)
        j = j + 1;
        mindiff = varargin{j};
    end
    j = j+1;
end


hold off;
for j = 1:length(pdata)
  fit = fitpsf(pdata(j).data,plotarg{:});
  h(j) = plotpsych(fit.data, fit.fit(1), fit.fit(2), 'color', ...
		   colors{j},'shown',plotarg{:});
  labels{j} = sprintf('dt %.0f xo %.1f fit %.2f %.4f',pdata(j).delay, pdata(j).xo,fit.fit(1),fit.fit(2));
  hold on;
end
legend(h,labels);