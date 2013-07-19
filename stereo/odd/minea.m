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
	  power = 1;
	  factor = 0.75;
	hold off;
%	for power = [1 1.5 2 2.5 3 3.5]
	for factor = [0.48 0.49 0.5 0.51]
	  datestr(now)
	  sd = 4;
	  lf = xf - filtersize/2;
	  lefta = exp(-(lf.^2)/sd.^2);
	  lf = xf - filtersize/2+dx;
	  righta = exp(-(lf.^2)/sd.^2);
	  lf = xf - filtersize/2-dx;

	  sd = sd * factor;
	  rightb = (exp(-abs(lf ./ sd).^power));
	  rightb = -1.69 * rightb;
	  leftb = (exp(-abs(lf ./ sd).^power));
	  leftb = 1.69 *leftb;

	  uncorr = 0;
	  for j = 1:500
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
	if(j == 1)
	  suma = answer;
	else
	  suma = suma + answer;
	end
	end %for j
	avg = suma/j;
	xlabel(sprintf('%d iterations power = %.1f SDR = %.2f',j*itterations,power,factor));
	ylabel('TI response');
	legend('sum','+','-');
%  This line plots R(x) vs R(-x)
if dx == 0
	ratio = (max(avg(:,2)) - min(avg(:,2)))/(max(avg(:,3)) - ...
						 min(avg(:,3)));	
	plot(avg(:,2),avg(:,3));
	text(avg(16,2),avg(16,3),sprintf('SDR %.2f RR %.3f',factor,ratio));
        hold on;
	plot([min(avg(:,2)) max(avg(:,2))],[max(avg(:,3)) min(avg(:,3))]);
else
	ratio = (max(avg(:,2)) - min(avg(:,2)))/(max(avg(:,3)) - ...
						 min(avg(:,3)));	
	plot(avg(2:17,1),flipud(avg(17:32,1)));
end
        hold on;
%	plot(xf,lefta);
	legend(sprintf('Range ratio %3f',ratio));
	drawnow;
	end
	
