%	1 dimentional 4f+3f sine stimulus  stimulus
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
	clear;
	const = .6;
	frequency = 4;
	filtersize = 64;
	sd = const * filtersize/frequency;

	stim_high_cycles = 4; %number of cycles per filtersize - other will be * .75
	rf_disparity = 0; %range -filtersize/4 -> +filtersize/4
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

%find response of complex cell to sine components alone
%first the high freq component
	for disp = -filtersize/2:filtersize/2
	    response4 = Run_one_sine(disp, rf_disparity, filtersize, sd, frequency, stim_high_cycles);
	    a = analyze(1, response4); %analyzed once for each phase
	    ans_sin4 = [ans_sin4; a];
	end


	
%build a final answer matrix of only the complex cell responses
%	sin_average = (ans_sin4(:,5) + ans_sin3(:,5))/2;
%	answer = [ans_sin43(:,5), ans_sin4(:,5), ans_sin3(:,5), sin_average];
%	plot(answer);
	plot(ans_sin4);
	xlabel('stimulus disparity');
	ylabel('response');
	legend('lc','ls','rc','rs','cplx');


