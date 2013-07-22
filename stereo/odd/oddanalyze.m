function answer=oddanalyze(trials, response)
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
%setting thresh to a large negative number restores the
%original ODF model
	thresh = 0;
%	thresh = -100000;
	
% binocular response is linear sum of monouclar responses for simple cells
% and is squared sum of all binocular simple cells for complex cell
% mututal inhibition - essential for spatial antagonism at flanks modelled by
% winner take all.
% response = [left cos, left sin, right cos, right sin]

	for i = 1:trials
	        lc = max([thresh, response(i,1)]);
	        rs = max([thresh, response(i,4)]);
	        rc = max([thresh, response(i,3)]);
	        ls = max([thresh, response(i,2)]);
	        lcn = max([thresh, -1 * response(i,1)]);
	        rsn = max([thresh, -1 * response(i,4)]);
	        rcn = max([thresh, -1 * response(i,3)]);
	        lsn = max([thresh, -1 * response(i,2)]);
%
%       Standard ODF parins for phase of -pi/2
                cossum = lc + rs;
		cossumn = lcn + rsn;
		if cossum >= 0 
			simple_ON_cos = simple_ON_cos + cossum;
		else 
			simple_OFF_cos = simple_OFF_cos + cossum;
		end

		sinsum = ls +rc;
		sinsumn = lsn +rcn;
		if sinsum >= 0
			simple_PLUS_sin = simple_PLUS_sin + sinsum;
		else 
			simple_MINUS_sin = simple_MINUS_sin + sinsum;
		end

		complex = complex + (cossum^2 + sinsum^2 + cossumn^2 ...
				     + sinsumn^2);
	end
	answer = [simple_ON_cos, simple_OFF_cos, simple_PLUS_sin, simple_MINUS_sin, complex];
	answer = complex/trials;

