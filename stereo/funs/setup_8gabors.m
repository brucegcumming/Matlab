function [rf_left_cos, rf_left_sin, rf_right_cos, rf_right_sin, rf_left_cos2, rf_left_sin2, rf_right_cos2, rf_right_sin2] = Setup_gabor_rfs(rf_disparity, filtersize, sd, cycles)
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

%setup axis
	xf = 1:filtersize;
	centre = xf - filtersize/2;
	dblcycles = cycles *2;
	dblsd=sd/2;
% make and scale gabor filters (equal area with change of sd)
	rf_left_cos = Gabor(centre, sd, cycles, pi);
	rf_left_sin = Gabor(centre, sd, cycles, pi/2);
	rf_left_cos2 = Gabor(centre, dblsd, dblcycles, pi);
	rf_left_sin2 = Gabor(centre, dblsd, dblcycles, pi/2);
	gaussian = Gabor(xf - (filtersize/2),32,0,0);
	scale = sum(gaussian);
	rf_left_cos = rf_left_cos/scale;
	rf_left_sin = rf_left_sin/scale;
	rf_left_cos2 = rf_left_cos2/scale;
	rf_left_sin2 = rf_left_sin2/scale;

% now make right eye with appropriate disparity
	rf_right_cos = Gabor(centre + rf_disparity, sd, cycles, pi);
	rf_right_sin = Gabor(centre + rf_disparity, sd, cycles, pi/2);
	rf_right_cos2 = Gabor(centre + rf_disparity, dblsd, dblcycles, pi);
	rf_right_sin2 = Gabor(centre + rf_disparity, dblsd, dblcycles, pi/2);
	rf_right_cos = rf_right_cos/scale;
	rf_right_sin = rf_right_sin/scale;	
	rf_right_cos2 = rf_right_cos2/scale;
	rf_right_sin2 = rf_right_sin2/scale;	
