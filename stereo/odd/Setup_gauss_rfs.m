function [rf_left_a, rf_left_b, rf_right_a, rf_right_b] = Setup_gauss_rfs(rf_disparity, filtersize, sd)
%
%	[l_c, l_s, r_c, r_s] = Setup_gabor_rfs(rf_disparity, filtersize, sd, cycles)
%
%	cycles	:cycles per filtersize
%	returns 4 by filtersize matrix
%	[lefteye_cos;lefteye_sin;righteye_cos;righteye_sin]
%
%	left eye rf is fixed in the centre (x=filtersize/2)
%	right eye rf shifts, zero disparity = filtersize/2
%


        xf = 1:filtersize;
	dx = rf_disparity;
	lf = xf - filtersize/2;
	rf_left_a = exp(-(lf.^2)/sd.^2);
	rf_left_b = rf_left_a;
	lf = xf - filtersize/2+dx;
	rf_right_a = exp(-(lf.^2)/sd.^2);
	lf = xf - filtersize/2-dx;
	rf_right_b = -exp(-(lf.^2)/sd.^2);










