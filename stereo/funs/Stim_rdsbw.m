function [leftstim, rightstim]=Stim_rds_1d(stim_disparity, filtersize)
%
%
%	[leftstim, rightstim]=Stim_rds_1d(stim_disparity, filtersize)
%
%	returns 1d random contrast stimuli
%

%initialize an array
	x = stim_disparity + filtersize/2;
	disp = abs(stim_disparity);
	xf = 1:filtersize;
%only allow stim disparity if less than half filtersize
	if disp > (filtersize/2);
		error('Stim_Disparity too big');
	end
%take random sample for stimulus
	leftstim = Newsample(filtersize);	
%make right eye stim and copy over left eye stuff execpt for
%occluded regions - insert random stuff in there
	if stim_disparity == 0
		rightstim = leftstim;
	else 
		rightstim(1:disp) = Newsample(disp);
		rightstim((disp + 1):filtersize) = leftstim(1:(filtersize - disp));
%swap over images if negative disparity
		if stim_disparity < 0
			a = rightstim;
			rightstim = leftstim;
			leftstim = a;
		end
	end
% include next line if you want anticorrelated
%leftstim = -leftstim;
