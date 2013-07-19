function response = Run_one_rds(itterations, stim_disparity, rf_disparity, filtersize, sd, freq)
%
%	response = Run_one_rds(itterations, stim_disparity, rf_disparity, filtersize, sd, freq)
%
%	returns 4 by itterations array of responses of 4 rf's to random contrast stimuli
%	[left_eye_cos, left_eye_sin, right_eye_cos, right_eye_sin]
%	... where the rf's are odd and even gabors shifted the appropriate disparity
%
%

% set appropriate rf's
% setup_4gabors returns 2 subunits for each eye (cos, sin)
% setup_8gabors returns 4 subunits for each eye, with different SF 
% Setup_odd_gabor RFs returns sin,cos for left, -sin,cos for right, for
% odd symetric tuning

	[rf_left_cos, rf_left_sin, rf_right_cos, rf_right_sin,rf_left_cos2, rf_left_sin2, rf_right_cos2, rf_right_sin2] = setup_8gabors(rf_disparity, filtersize, sd, freq);
% calculate response to different random 1d patterns
	for i = 1:itterations
		[leftstim, rightstim] = Stim_rds_1d(stim_disparity, filtersize);
		a(1) = sum(leftstim .* rf_left_cos);
		a(2) = sum(leftstim .* rf_left_sin);
		a(3) = sum(rightstim .* rf_right_cos);
		a(4) = sum(rightstim .* rf_right_sin);
		a(5) = sum(leftstim .* rf_left_cos2);
		a(6) = sum(leftstim .* rf_left_sin2);
		a(7) = sum(rightstim .* rf_right_cos2);
		a(8) = sum(rightstim .* rf_right_sin2);
		a(9) = sum(-leftstim .* rf_left_cos);
		a(10) = sum(-leftstim .* rf_left_sin);
		a(11) = sum(-leftstim .* rf_left_cos2);
		a(12) = sum(-leftstim .* rf_left_sin2);
		response = [response;a];
	end
