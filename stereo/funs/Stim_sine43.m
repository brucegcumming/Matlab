function [leftstim, rightstim] = Stim_sine43(stim_disparity, filtersize, highcycles, phase)
%
%	Stim_sine43(stim_disparity, filtersize, highfreq, phase)
%
%	highfreq	:cycles/filtersize
%	phase		:integer i.e. 1:filtersize gives complete cycle 
%	[leftstim, rightstim]=stim_sine43sine_stim(stim_disparity, filtersize, phase)
%
%	returns 1d 3f+4f sinusoidal at required phase at required disparity
%	

	lowcycles = highcycles * 0.75;
	
	lf = Cycles2const(lowcycles, filtersize);
	hf = Cycles2const(highcycles, filtersize);
	
	buf = 1:filtersize;
	buf = buf+phase; %this gives the x value for the sin calculation at correct phase

%make leftstim 
	low = sin(lf .* buf);
	high = sin(hf .* buf);
	leftstim = (low + high);
%make rightstim
	low = sin(lf .* (buf - stim_disparity));
	high = sin(hf .* (buf - stim_disparity));
	rightstim = (low + high);
	
