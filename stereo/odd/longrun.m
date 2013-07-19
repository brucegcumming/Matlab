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

	clear all;
	const = .3;
	frequency = 2;
	filtersize = 32;
	sd = 4;
	response = [0 0 0 0];
	therf = [0 0 0 0];
	itterations = 100; %thats per disparity increment
	clear avg;

	
	
	rf_disparity = 4; %range -filtersize/4 -> +filtersize/4
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
        xf = 1:filtersize;
	dx = 3;
	lf = xf - filtersize/2;
	lefta = exp(-(lf.^2)/sd.^2);
	lf = xf - filtersize/2+dx;
	righta = exp(-(lf.^2)/sd.^2);
	lf = xf - filtersize/2-dx;

%1.15, 0.58 is good with gain at 2 give ratio of 1.03
%1.15, 0.56 is better with gain at 2 give ratio of 1.03 errmax 0.23
	factor = 0.58;
	factor = 0.56;
	power = 1.15;
	
	
	  sd = sd * factor;
	  rightb = (exp(-abs(lf ./ sd).^power));
	  rightb = -2 * rightb;
	  leftb = (exp(-abs(lf ./ sd).^power));
	  leftb = 2 *leftb;


	uncorr = 0;
	uncorra = 0;
	uncorrb = 0;
	for j = 1:100000
	for i = -(filtersize/2):filtersize/2
		response = Run_one_rds(itterations, i, filtersize, ...
lefta, leftb, righta, rightb);
		a = bruceanalyze(itterations, response);
		if(i > -filtersize/2)
		  answer = [answer; a];
		else
		  answer = [a];
		end
	end
	uncorr = uncorr + (answer(1,1) + answer(filtersize,1))/2;
	uncorra = uncorra + (answer(1,2) + answer(filtersize,2))/2;
	uncorrb = uncorrb + (answer(1,3) + answer(filtersize,3))/2;
	fprintf('uncorr %.1f %.1f %.1f\n',uncorr/j,uncorra/j,uncorrb/j);
	if(j == 1)
	  suma = answer;
	else
	  suma = suma + answer;
	end
	hold off;
	subplot(2,1,1);
	hold off;
	avg = suma/j;
        do = avg(:,1) - uncorr/j;				      ;
        da = avg(:,2) - uncorra/j;				      ;
        db = uncorrb/j - avg(:,3);
	plot(avg);
	hold on;
	if dx == 0
	  ratio = (max(avg(:,2)) - min(avg(:,2)))/(max(avg(:,3)) - ...
						 min(avg(:,3)));	
	  errs = da - db * ratio;
	  plot(abs(errs));
	  crit = (max(da) + max(db * ratio)) * 0.05;
	  plot([1 filtersize],[crit crit]);
	  errsum = sum(errs);
	  errmax = max(abs(errs));
	else
	  ratio = max(do)/(-min(do));
	  errs = do(1:15) - flipud(do(16:30));
	  errsum = sum(errs);
	  errmax = max(abs(errs));
	  plot(errs);
	end
	xlabel(sprintf('%d iterations err %.2f max %.2f',j*itterations,errsum,errmax));
	ylabel('response');
	legend('sum','+','-');
	subplot(2,1,2);
	hold off

	if dx == 0
	  plot(avg(:,2),avg(:,3));
	  hold on;
	  plot([min(avg(:,2)) max(avg(:,2))],[max(avg(:,3)) min(avg(:,3))],'r');
	else
	  %  This line plots R(x) vs R(-x)
	  plot(avg(1:15,1),flipud(avg(16:30,1)));
	  hold on;
%  This to plot the two components against each other (assuming symmetry)  
%	  plot(avg(:,2),flipud(avg(:,3)),'g');
	  plot([min(avg(1:15,1)) max(avg(1:15,1))],[max(avg(17:32,1)) min(avg(17:21,1))],'r');
	end

	legend(sprintf('power %.2f factor %.2f RRatio %3f',power,factor,ratio));
	drawnow;
	end
