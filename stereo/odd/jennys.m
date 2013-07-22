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
	const = 1;
	frequency = 2; %in cycles per filtersize
	filtersize = 32;
	sd = 4; %in pixels
	response = [0 0 0 0];
	therf = [0 0 0 0];
	itterations = 200; %thats per disparity increment

	
	
	rf_disparity = 4; %range -filtersize/4 -> +filtersize/4
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
        xf = 1:filtersize;
	lf = xf - filtersize/2;
	filtersize = size(xf,2);
	f = 2*pi*frequency/filtersize;
	phasel = pi/4;
	lefta = exp( - (lf.^2/(2 *sd.^2))) .* (cos(f .* lf + phasel));
	leftb = exp( - (lf.^2/(2 *sd.^2))) .* (cos(f .* lf + phasel));
	phaser = (3 *pi)/4;
	righta = exp( - (lf.^2/(2 *sd.^2))) .* (cos(f .* lf + phaser));
	rightb = exp( - (lf.^2/(2 * sd.^2))) .* (cos(f .* lf + phaser));
	dp = phasel-phaser;
	jcr = 1/(8 * pi * sd)* exp(-lf.^2 ./ (4 .* sd^2)) .* (cos(f ...
						  .* lf - dp) + ...
						  exp(-sd^2 * ...
						      f^2)*cos(phaser ...
						  + phasel));
	uncorr = 0;
	for j = 1:10000
	for i = lf;
		response = Run_one_rds(itterations, i, filtersize, ...
				       lefta, leftb, righta, rightb);
		a = odfanalyze(itterations, response);
		if(i > lf(1))
		  answer = [answer; a];
		else
		  uncorr = uncorr + a(1);
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
	plot(avg(:,1)-uncorr/j);
	hold on;
	plot((lefta .* uncorr/j) ,'r');
	plot((righta .* uncorr/j) ,'g');
	plot(((jcr) / (max(jcr) - min(jcr))) * (max(avg(:,1)) ...
						  - min(avg(:,1))),'k');
	hold off;
	xlabel(sprintf('%d iterations',j*itterations));
	ylabel('response');
	legend('sum','+','-');
	subplot(2,1,2);
	plot(avg(2:17,1),flipud(avg(17:32,1)));
	hold on;
	plot([min(avg(1:15,1)) max(avg(1:15,1))],[max(avg(17:32,1)) ...
		    min(avg(17:21,1))],'r');
	hold off;
	ratio = (max(avg(:,2)) - min(avg(:,2)))/(max(avg(:,3)) - ...
						 min(avg(:,3)));	
	legend(sprintf('Range ratio %3f',ratio));
	drawnow;
	end
