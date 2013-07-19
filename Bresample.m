function resampled = Bresample(values,varargin)
% Resample(values,nrep,nout)
% Given a list of n values, this function replaces them with nout random resamples from the same list.
% Resample(values)
% If nout is not supplied, it is assumed to be n (so resampled and values have the same size)

n = length(values);
nrep = 1;
nv = length(varargin);
if nv==0
    nout = n;
elseif nv==1
    nrep = varargin{1};
    nout = n;
elseif nv==2
    nrep = varargin{1};
    nout = varargin{2};
else
    disp('Too many input arguments in Resample.m')
    stop
end

% Pick n random numbers between 1 and n: these will be the resampled values
newindices = unidrnd(n,nrep,nout);
resampled = values(newindices);
