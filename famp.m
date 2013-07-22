function [a, c] = famp(x,y,f, varargin)

% [a, c] = famp(x,y,f)
% calculates spectral power of signal y, over time samples x, 
% at frequency f.
% returns  c  Amlitude (complex number)
% returns  a  abs(c) 

if size(x,1) == size(y,2)
    y = y';
end
pointwise = 0;

j=1;
while j < nargin - 2
    if strncmpi(varargin{j},'point',3)
        pointwise = 1;
    end
    j = j+1;
end
if pointwise
    sigw = mean(diff(x)) .* length(x);
    c = 2 * sum(y .* exp(2*pi*i*f.*x))/sigw;
else
    c = 2 * trapz(x, y .* exp(2*pi*i*f.*x))/(max(x) - min(x));
end
    a = abs(c);
