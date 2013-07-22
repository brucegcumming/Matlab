function answer=analyze(trials, response)
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
	simple_ON_cos2 =  0;
	simple_OFF_cos2 = 0;
	simple_PLUS_sin2 = 0;
	simple_MINUS_sin2 = 0;
	complex = 0;
	antisimple_ON_cos =  0;
	antisimple_OFF_cos = 0;
	antisimple_PLUS_sin = 0;
	antisimple_MINUS_sin = 0;
	antisimple_ON_cos2 =  0;
	antisimple_OFF_cos2 = 0;
	antisimple_PLUS_sin2 = 0;
	antisimple_MINUS_sin2 = 0;
	anticomplex = 0;
	
% binocular response is linear sum of monouclar responses for simple cells
% and is squared sum of all binocular simple cells for complex cell
% mututal inhibition - essential for spatial antagonism at flanks modelled by
% winner take all.
% response = [left cos, left sin, right cos, right sin, left cos2, left sin2, right cos2, right sin2]

	for i = 1:trials
		cossum = response(i,1) + response(i,3);
		if cossum >= 0 
			simple_ON_cos = simple_ON_cos + cossum;
		else 
			simple_OFF_cos = simple_OFF_cos + cossum;
		end

		cossum2 = response(i,5) + response(i,7);
		if cossum2 >= 0 
			simple_ON_cos2 = simple_ON_cos2 + cossum2;
		else 
			simple_OFF_cos2 = simple_OFF_cos2 + cossum2;
		end
		
		sinsum = response(i,2) + response(i,4);
		if sinsum >= 0
			simple_PLUS_sin = simple_PLUS_sin + sinsum;
		else 
			simple_MINUS_sin = simple_MINUS_sin + sinsum;
		end

		sinsum2 = response(i,6) + response(i,8);
		if sinsum2 >= 0
			simple_PLUS_sin2 = simple_PLUS_sin2 + sinsum2;
		else 
			simple_MINUS_sin2 = simple_MINUS_sin2 + sinsum2;
		end
		
		cosantisum = response(i,9) + response(i,3);
		if cosantisum >= 0 
			antisimple_ON_cos = antisimple_ON_cos + cosantisum;
		else 
			antisimple_OFF_cos = antisimple_OFF_cos + cosantisum;
		end

		cosantisum2 = response(i,11) + response(i,7);
		if cosantisum2 >= 0 
			antisimple_ON_cos2 = antisimple_ON_cos2 + cosantisum2;
		else 
			antisimple_OFF_cos2 = antisimple_OFF_cos2 + cosantisum2;
		end
		
		sinantisum = response(i,10) + response(i,4);
		if sinantisum >= 0
			antisimple_PLUS_sin = antisimple_PLUS_sin + sinantisum;
		else 
			antisimple_MINUS_sin = antisimple_MINUS_sin + sinantisum;
		end

		sinantisum2 = response(i,12) + response(i,8);
		if sinantisum2 >= 0
			antisimple_PLUS_sin2 = antisimple_PLUS_sin2 + sinantisum2;
		else 
			antisimple_MINUS_sin2 = antisimple_MINUS_sin2 + sinantisum2;
		end
		
		complex = complex + (cossum^2 + sinsum^2 + sinsum2^2 + cossum2^2);		
		anticomplex = anticomplex + (cosantisum^2 + sinantisum^2 + sinantisum2^2 + cosantisum2^2);
	end
	answer = [simple_ON_cos, simple_ON_cos2, simple_PLUS_sin, simple_PLUS_sin2, complex, anticomplex];
	answer = answer/trials;

