function [leftstim, rightstim] = Stim_sine(stim_disparity, filtersize, highcycles, phase)
%
%	Stim_sine(stim_disparity, filtersize, cycles, phase)
%

	f = Cycles2const(highcycles, filtersize);
	
	buf = 1:filtersize;
	buf = buf+phase; %this gives the x value for the sin calculation at correct phase

%make leftstim 
	leftstim = sin(f .* buf);
%make rightstim
	rightstim = sin(f .* (buf - stim_disparity));
	
