function answer=bcanalyze(trials, response)
%	
%	answer=analyze(trials)
%	
%	gets binocular responses from monocular inputs
%
%
	simple_ON_cos =  0;
	simple_OFF_cos = 0;
	simple_PLUS_sin = 0;
	simple_MINUS_sin = 0;
	complex = 0;
	
% binocular response is linear sum of monouclar responses for simple cells
% and is squared sum of all binocular simple cells for complex cell
% mututal inhibition - essential for spatial antagonism at flanks modelled by
% winner take all.
% response = [left cos, left sin, right cos, right sin]

	for i = 1:trials
	  lc = rc = ls = rs = 0;
	        if(response(i,1) > 0)
		  lc = response(i,1);
		end
	        if(response(i,3) > 0)
		  rc = response(i,3);
		end
	        if(response(i,2) > 0)
		  ls = response(i,2);
		end
	        if(response(i,4) > 0)
		  rs = response(i,4);
		end
		
		cossum = lc + rc;
		if cossum >= 0 
			simple_ON_cos = simple_ON_cos + cossum;
		else 
			simple_OFF_cos = simple_OFF_cos + cossum;
		end

		sinsum = ls + rs;
		if sinsum >= 0
			simple_PLUS_sin = simple_PLUS_sin + sinsum;
		else 
			simple_MINUS_sin = simple_MINUS_sin + sinsum;
		end

		complex = complex + (cossum^2 + sinsum^2);
	end
	answer = [simple_ON_cos, simple_OFF_cos, simple_PLUS_sin, simple_MINUS_sin, complex];
	answer = answer/trials;
	answer = complex/trials;

