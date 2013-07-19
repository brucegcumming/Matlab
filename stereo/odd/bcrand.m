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
	const = 2
	frequency = 4;
	filtersize = 32;
	sd = const * filtersize/frequency;
	response = [0 0 0 0];
	therf = [0 0 0 0];
	itterations = 4000; %thats per disparity increment

	
	
	rf_disparity = 0; %range -filtersize/4 -> +filtersize/4
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%	[left_cos, left_sin, right_cos, right_sin] = Setup_odd_gabor_rfs(rf_disparity, filtersize, sd, frequency);
	[left_cos, left_sin, right_cos, right_sin] = Setup_other_gabor_rfs(rf_disparity, filtersize, sd, frequency);

	for i = -(filtersize/2):filtersize/2
		response = Run_one_rds(itterations, i, filtersize, left_cos, left_sin, right_cos, right_sin);
		a = oddanalyze(itterations, response);
		if(i > -filtersize/2)
		  answer = [answer; a];
		else
		  answer = [a];
		end
	end		
	plot(answer);
	xlabel('stimulus disparity');
	ylabel('response');
	legend('on','off','+','-','complex');
