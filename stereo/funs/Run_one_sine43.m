function response = Run_one_sine43(stim_disparity, rf_disparity, filtersize, sd, freq, stimcycles)
%
%	response = Run_one_sine43(stim_disparity, rf_disparity, filtersize, sd, freq, stimcycles)
%
%	returns 4 by itterations array of responses of 4 rf's to 4f3f sine stimuli stimuli
%	[left_eye_cos, left_eye_sin, right_eye_cos, right_eye_sin]
%	... where the rf's are odd and even gabors shifted the appropriate disparity
%	... stimcycles is no of cycles of high freq sine per filtersize
%

% set appropriate rf's
	[rf_left_cos, rf_left_sin, rf_right_cos, rf_right_sin] = Setup_gabor_rfs(rf_disparity, filtersize, sd, freq);
% calculate response to different random 1d patterns
	for phase = 1:filtersize
		[leftstim, rightstim] = Stim_sine43(stim_disparity, filtersize, stimcycles, phase);
		a(1) = sum(leftstim .* rf_left_cos);
		a(2) = sum(leftstim .* rf_left_sin);
		a(3) = sum(rightstim .* rf_right_cos);
		a(4) = sum(rightstim .* rf_right_sin);
		response = [response;a];
	end





