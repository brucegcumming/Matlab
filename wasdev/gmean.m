function u = gmean(x, varargin)

u = exp(mean(log(x,varargin{:})));
