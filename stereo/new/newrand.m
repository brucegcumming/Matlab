function newrand(niter)
% niter = number of iterations to run
%	1 dimentional random contrast stimulus
%
% This is a script which can be used to run the stereo rf modeling
%
% m files used run_one_set -> newstimuli -> newsample
%			   -> new_rf
%		analyze
%
% Run_one_rds.m calls setup_4gabors
%in funs		in funs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the parameters here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	const = .4;
	frequency = 4;
	filtersize = 32;
	sd = const * filtersize/frequency;
	answer = [0 0 0 0 0 0 0 0 0 0];
	subset = [0];
	response = [0 0 0 0];
	rf_disparity = 1; %range -filtersize/4 -> +filtersize/4
	itterations = 20;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if (niter == 0)
		load convolutions;
	else
		allresps = convolve(niter, rf_disparity, filtersize, sd, frequency);
		save convolutions allresps;
	end
	
% dimension 1 of allresps is the number of iterations in the file saved
% to disk
	msiz = size(allresps);
	for i = 1:filtersize
		a = bcanalyze(msiz(1)-1, allresps(:,:,i));
		answer = [answer; a];
	end
	save acmodel.dat answer -ascii
	even = figure;
	plot(answer(:,1:4));
	sum = figure;	
	plot(answer(:,9:10));
	odd= figure;
	plot(answer(:,5:8));
	xlabel('stimulus disparity');
	ylabel('response');
	legend('Corr','AC','Rect Corr','Rect AC');

%
%Convolve function runs a set of convolutsions and stores teh
%results to disk, then can load them in and do various calculations
%

function allresps = convolve(itterations, rf_disparity, filtersize, sd, frequency)
% set appropriate rf's
% setup_gabors returns 2 subunits for each eye (cos, sin)

	response = [0 0 0 0];

	subunits = setup_gabors(rf_disparity, filtersize, sd, frequency);
	allresps = [];

	n = 1;
	for i = -filtersize/2:filtersize/2
		response = Run_one_rds(itterations, i, subunits, filtersize);
		allresps(:,:,n) = response;
		n = n + 1;
	end


