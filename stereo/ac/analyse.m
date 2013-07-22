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
	complex = 0;
	coscomplex = 0;
	antisimple_ON_cos =  0;
	antisimple_OFF_cos = 0;
	antisimple_PLUS_sin = 0;
	antisimple_MINUS_sin = 0;
	anticomplex = 0;
	cosanticomplex = 0;
	rl = 0; antirl = 0;
	sumrl = 0; sumantirl = 0; rlusum = 0;
	ucomplex = 0;
	
% binocular response is linear sum of monouclar responses for simple cells
% and is squared sum of all binocular simple cells for complex cell
% mututal inhibition - essential for spatial antagonism at flanks modelled by
% winner take all.
% response = [left cos, left sin, right cos, right sin, left cos2, left sin2, right cos2, right sin2]

	for i = 1:trials
		cossum = response(i,1) + response(i,3);
		ucossum = response(i,1) + response(i,7);
		rl = response(i,1) * response(i,3);
		antirl = response(i,5) * response(i,3);
		rlu = response(i,1) * response(i,7);
		if cossum >= 0 
			simple_ON_cos = simple_ON_cos + cossum;
		else 
			simple_OFF_cos = simple_OFF_cos + cossum;
		end

		sinsum = response(i,2) + response(i,4);
		usinsum = response(i,2) + response(i,8);
		if sinsum >= 0
			simple_PLUS_sin = simple_PLUS_sin + sinsum;
		else 
			simple_MINUS_sin = simple_MINUS_sin + sinsum;
		end

		cosantisum = response(i,5) + response(i,3);
		if cosantisum >= 0 
			antisimple_ON_cos = antisimple_ON_cos + cosantisum;
		else 
			antisimple_OFF_cos = antisimple_OFF_cos + cosantisum;
		end

		sinantisum = response(i,6) + response(i,4);
		if sinantisum >= 0
			antisimple_PLUS_sin = antisimple_PLUS_sin + sinantisum;
		else 
			antisimple_MINUS_sin = antisimple_MINUS_sin + sinantisum;
		end

		complex = complex + (cossum^2 + sinsum^2);		
		anticomplex = anticomplex + (cosantisum^2 + sinantisum^2);
		ucomplex = ucomplex + (ucossum^2);
		coscomplex = coscomplex + (cossum^2);		
		cosanticomplex = cosanticomplex + (cosantisum^2);
		sumrl = sumrl + rl;
		sumantirl = sumantirl + antirl;
		rlusum = rlusum + rlu;
	end
	answer = [coscomplex, cosanticomplex, ucomplex, sumrl, sumantirl, rlusum];
%	answer = [coscomplex, cosanticomplex,  complex, anticomplex];
	answer = answer/trials;

