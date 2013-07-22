function S = parzenSurf(x, y, z, Xg, Yg, smoothParzen)
%
% parzenSurface -------------------------------------
% function S = parzenSurf(x, y, z, Xg, Yg, smoothParzen)
%
% make Parzen estimate of surface data in (x, y, z) triplets
% with gaussian kernel widths based on smoothParzen (0.1-3).
% output is S over Xg, Yg meshgrids
%
% 24may00 LM Optican

holes = 0;	% 1 = surface estimation with holes

S = zeros(size(Xg));	% data surface

if (holes)
	nS = S;			% normalization surface
end

uz = [];
ux = unique(x);
uy = unique(y);
sx = smoothParzen(1);
if length(smoothParzen) ==2
sy = smoothParzen(2);
else
    sy = sx;
end
%sy = sx * ( (uy(2) - uy(1)) / (ux(2) - ux(1)) ) / 2;

scale = (mean(diff(ux)) * mean(diff(uy))/(sx * sy))/(2 * pi);
%
% subtract out mean to eliminateedge effects. Add this back at the end;
%

zmean = mean(mean(z(find(~isnan(z)))));
z = z - zmean;

[ m , n ] = size(z);
for j = 1:n
for i = 1:m
%    x2 = ((Xg - x(i,j)) / sx).^2;
%	y2 = ((Yg - y(i,j)) / sy).^2;
% 
    x2 = ((Xg - x(i,j)) / sx).^2;
	y2 = ((Yg - y(i,j)) / sy).^2;
	R2 = x2 + y2;

%  normalization depends on both SD and area of data cell.
%  by sx, sy to normalize area under kernel. BGC Jan 2005.
    G = exp(-R2 / 2) * scale;
	if (holes)
		nS = nS + G;	% use for surface estimation
    end
    
    if ~isnan(z(i,j))
        S = S + G * z(i,j);
    end
end
end
scale = mean(mean(S));
if (holes)
	% normalize
	i = find(nS > 0);
	S(i) = S(i) ./ nS(i);
end
	[ m , n ] = size(z);
S = S + zmean;