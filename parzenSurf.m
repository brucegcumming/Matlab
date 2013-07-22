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
sx = smoothParzen;
sy = smoothParzen;
%sy = sx * ( (uy(2) - uy(1)) / (ux(2) - ux(1)) ) / 2;


n = length(z);
for i = 1:n
x2 = ((Xg - x(i)) / sx).^2;
	y2 = ((Yg - y(i)) / sy).^2;
	R2 = x2 + y2;

        G = exp(-R2 / 2);

	if (holes)
		nS = nS + G;	% use for surface estimation
	end

	S = S + G * z(i);
end

if (holes)
	% normalize
	i = find(nS > 0);
	S(i) = S(i) ./ nS(i);
end
