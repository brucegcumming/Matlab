function [sx, sy] = xysmooth(x, y, w, varargin)
%[sx,xy] = xsmooth(x, y, w, ...) smooth data with boxcar, width w
%  if w is omitted, 5 is used. 
%smooth(x,y, w,'gauss') uses Gaussian smoothing with an SD of w

if nargin < 3
    w = 5
end
 if w == 0
     w = 1;
 end

j = 1;
kernel = ones(1,w)/w;
while j <= nargin-3
    if strncmpi(varargin{j},'gauss',4)
        kernel = Gauss(w,[-(w*2):(w*2)]);
        kernel = kernel./sum(kernel);
    elseif strncmpi(varargin{j},'kernel',4)
        j = j+1;
        kernel = varargin{j};
    end
    j = j + 1;
end

if length(x) == 0 || length(y) == 0
    sx = 0;
    sy = 0;
    return;
end
[X,a] = sort(x);
Y = y(a);
sx = smooth(X,w,varargin);
sy = smooth(Y,w,varargin);
