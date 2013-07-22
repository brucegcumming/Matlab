function answer=analyze(trials, response)
%	
%	answer=analyze(trials)
%	
%	gets binocular responses from monocular inputs
%
%
%	      
	simple_ON_cos =  0;
	simple_OFF_cos = 0;
	simple_PLUS_sin = 0;
	simple_MINUS_sin = 0;
	rectresp = [0 0 0 0 0 0 0 0];
	arectresp = [0 0 0 0 0 0 0 0];
	rsinsum = [0 0];
	rcossum = [0 0];
	rodd = [0 0 0 0];
	rsodd = [0 0 0 0];
	acomplex = 0;
	rcomplex = 0;
	racomplex = 0;
	complex = 0;
	oddcx = 0;
	aoddcx = 0;
	raoddcx = 0;
	roddcx = 0;
	lsum = 0;
	rsum = 0;
	rasum = 0;
	sumx = 0;
	asumx = 0;
	
% binocular response is linear sum of monouclar responses for simple cells
% and is squared sum of all binocular simple cells for complex cell
% response = [left cos, left sin, right cos, right sin]

	for i = 1:trials
		for j = 1:4
			if(response(i,j) > 0)
				rectresp(i,j) = response(i,j);
				arectresp(i,j) = 0;
			else
				arectresp(i,j) = response(i,j) * -1;
				rectresp(i,j) = 0;
			end
		end
		for j = 5:8
			if(response(i,j-4) < 0)
				rectresp(i,j) = response(i,j-4) * -1; 
				arectresp(i,j) = 0;
			else
				arectresp(i,j) = response(i,j-4);
				rectresp(i,j) = 0;
			end
		end

		lsum = rectresp(i,1) + rectresp(i,2) + rectresp(i,5) + rectresp(i,6);
		rsum = rectresp(i,3) + rectresp(i,4) + rectresp(i,7) + rectresp(i,8);
		rasum = arectresp(i,3) + arectresp(i,4) + arectresp(i,7) + arectresp(i,8);
		acossum = response(i,1) - response(i,3);
		cossum = response(i,1) + response(i,3);
		rcossum(1) = rectresp(i,1) + rectresp(i,3);
		rcossum(2) = rectresp(i,5) + rectresp(i,7);
		racossum(1) = rectresp(i,1) + arectresp(i,3);
		racossum(2) = rectresp(i,5) + arectresp(i,7);
		if cossum >= 0 
			simple_ON_cos = simple_ON_cos + cossum;
		else 
			simple_OFF_cos = simple_OFF_cos + cossum;
		end

		odda = response(i,1) + response(i,4);
		oddb = response(i,2) - response(i,3);
		aodda = response(i,1) - response(i,4);
		aoddb = response(i,2) + response(i,3);
		rodd(1) = rectresp(i,1) + rectresp(i,8);
		rodd(2) = rectresp(i,5) + rectresp(i,4);
		rodd(3) = rectresp(i,2) + rectresp(i,3);
		rodd(4) = rectresp(i,6) + rectresp(i,7);
		raodd(1) = rectresp(i,1) + arectresp(i,8);
		raodd(2) = rectresp(i,5) + arectresp(i,4);
		raodd(3) = rectresp(i,2) + arectresp(i,3);
		raodd(4) = rectresp(i,6) + arectresp(i,7);

		sinsum = response(i,2) + response(i,4);
		asinsum = response(i,2) - response(i,4);
		rsinsum(1) = rectresp(i,2) + rectresp(i,4);
		rsinsum(2) = rectresp(i,6) + rectresp(i,8);
		rasinsum(1) = rectresp(i,2) + arectresp(i,4);
		rasinsum(2) = rectresp(i,6) + arectresp(i,8);

		if sinsum >= 0
			simple_PLUS_sin = simple_PLUS_sin + sinsum;
		else 
			simple_MINUS_sin = simple_MINUS_sin + sinsum;
		end

		aoddcx = aoddcx + (aodda^2 + aoddb^2);
		oddcx = oddcx + (odda^2 + oddb^2);
		complex = complex + (cossum^2 + sinsum^2);
		acomplex = acomplex + (acossum^2 + asinsum^2);
		rcomplex = rcomplex + (rcossum(1)^2 + rcossum(2)^2 + rsinsum(1)^2 + rsinsum(2)^2);
		racomplex = racomplex + (racossum(1)^2 + racossum(2)^2 + rasinsum(1)^2 + rasinsum(2)^2);
		raoddcx = raoddcx + (raodd(1)^2 + raodd(2)^2 + raodd(3)^2 + raodd(4)^2);
		roddcx = roddcx + (rodd(1)^2 + rodd(2)^2 + rodd(3)^2 + rodd(4)^2);
		sumx = sumx + (lsum + rsum)^2;
		asumx = asumx + (lsum + rasum)^2;
	end

	answer = [complex, acomplex, rcomplex, racomplex, oddcx, aoddcx, roddcx, raoddcx, sumx, asumx ];
	answer = answer/trials;

