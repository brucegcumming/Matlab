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

	clear;
	const = .8;
	frequency = 4;
	filtersize = 32;
	sd = const * filtersize/frequency;
	answer = [0 0 0 0 0];
	subset = [0];
	response = [0 0 0 0];
	itterations = 40; %thats per disparity increment
	rf_disparity = 0; %range -filtersize/4 -> +filtersize/4
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% set appropriate rf's
% setup_gabors returns 2 subunits for each eye (cos, sin)
	subunits = setup_gabors(rf_disparity, filtersize, sd, frequency);

	for i = -filtersize/2:filtersize/2
		response = Run_one_rds(itterations, i, subunits);
%		save test.mat response;
		a = bcanalyze(itterations, response);
		cx = a(5);
		subset = [subset; cx];
		answer = [answer; a];
	end		
	plot(answer);
	xlabel('stimulus disparity');
	ylabel('response');
	legend('on','off','+','-','complex');

