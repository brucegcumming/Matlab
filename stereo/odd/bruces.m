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
	const = .3;
	frequency = 2;
	filtersize = 32;
	sd = 4;
	response = [0 0 0 0];
	therf = [0 0 0 0];
	itterations = 200; %thats per disparity increment

	
	
	rf_disparity = 4; %range -filtersize/4 -> +filtersize/4
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
        xf = 1:filtersize;
	dx = 0;
	lf = xf - filtersize/2;
	lefta = exp(-(lf.^2)/sd.^2);
	lf = xf - filtersize/2+dx;
	righta = exp(-(lf.^2)/sd.^2);
	lf = xf - filtersize/2-dx;

	sd = sd * 0.8;
	power = 2;
	rightb = (exp(-abs(lf ./ sd).^power));
	rightb = -1.69 * rightb;
	leftb = (exp(-abs(lf ./ sd).^power));
	leftb = 1.69 *leftb;

	uncorr = 0;
	for j = 1:10000
	for i = -(filtersize/2):filtersize/2
		response = Run_one_rds(itterations, i, filtersize, ...
				       lefta, leftb, righta, rightb);
		a = bruceanalyze(itterations, response);
		if(i > -filtersize/2)
		  answer = [answer; a];
		else
		  uncorr = uncorr + a(1)/2;
		  answer = [a];
		end
	end
	fprintf('uncorr %.1f\n',uncorr/j);
	if(j == 1)
	  suma = answer;
	else
	  suma = suma + answer;
	end
	subplot(2,1,1);
	avg = suma/j;
	plot(avg);
	xlabel(sprintf('%d iterations',j*itterations));
	ylabel('response');
	legend('sum','+','-');
	subplot(2,1,2);
%  This line plots R(x) vs R(-x)
if dx == 0
	plot(avg(:,2),avg(:,3));
else
	plot(avg(2:17,1),flipud(avg(17:32,1)));
end

	ratio = (max(avg(:,2)) - min(avg(:,2)))/(max(avg(:,3)) - ...
						 min(avg(:,3)));	
	legend(sprintf('Range ratio %3f',ratio));
	drawnow;
	end
