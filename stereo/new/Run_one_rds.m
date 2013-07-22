function response = Run_one_rds(itterations, stim_disparity, rf, filtersize)
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
	response = [0 0 0 0];

% calculate response to different random 1d patterns
	for i = 1:itterations
		[leftstim, rightstim] = Stim_rds_1d(stim_disparity, filtersize);
		a(1) = sum(leftstim .* rf(1,:));
		a(2) = sum(leftstim .* rf(2,:));
		a(3) = sum(rightstim .* rf(3,:));
		a(4) = sum(rightstim .* rf(4,:));
		response = [response;a];
	end
