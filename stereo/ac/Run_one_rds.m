function response = Run_one_rds(itterations, stim_disparity, filtersize, lc, ls, rc, rs)
%
%	response = Run_one_rds(itterations, stim_disparity, rf_disparity, filtersize, sd, freq)
%
%	returns 4 by itterations array of responses of 4 rf's to random contrast stimuli
%	[left_eye_cos, left_eye_sin, right_eye_cos, right_eye_sin]
%	... where the rf's are odd and even gabors shifted the appropriate disparity
%
%
	offset = 0;
	response = [0 0 0 0 0 0 0 0];
% calculate response to different random 1d patterns
	for i = 1:itterations
		[leftstim, rightstim] = Stim_rds_1d(stim_disparity, filtersize);
		[uleftstim, urightstim] = Stim_rds_1d(stim_disparity, filtersize);
	aleftstim = -leftstim;
		a(1) = sum(leftstim .* lc) + offset;
		a(2) = sum(leftstim .* ls) + offset;
		a(3) = sum(rightstim .* rc) + offset;
		a(4) = sum(rightstim .* rs) + offset;
		a(5) = sum(aleftstim .* lc) + offset;
		a(6) = sum(aleftstim .* ls) + offset;
		a(7) = sum(uleftstim .* rc) + offset;
		a(8) = sum(uleftstim .* rs) + offset;
	for i = 1:8
		if a(i) < 0
			a(i) = 0;
		end
	end
		response = [response;a];
	end
