%	1 dimentional random contrast stimulus
%
% This is a script which can be used to run the stereo rf modeling
%
% m files used run_one_set -> newstimuli -> newsample
%			   -> new_rf
%		analyze
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the parameters here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	clear;
	const = .8;
	frequency = 4;
	filtersize = 32;
	sd = const * filtersize/frequency;

	
	itterations = 200; %thats per disparity increment
	
	rf_disparity = 0; %range -filtersize/4 -> +filtersize/4
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

	for i = -filtersize/2:filtersize/2
		response = Run_one_rds(itterations, i, rf_disparity, filtersize, sd, frequency);
		a = analyze(itterations, response);
		answer = [answer; a];
	end		
	plot(answer);
	xlabel('stimulus disparity');
	ylabel('response');
	legend('on','off','+','-','complex');
