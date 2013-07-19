function convole(itteration, rf_disparity, filtersize, sd, cycles)
% set appropriate rf's
% setup_gabors returns 2 subunits for each eye (cos, sin)

	response = [0 0 0 0];

	subunits = setup_gabors(rf_disparity, filtersize, sd, frequency);
	allresps = [];

	n = 1;
	for i = -filtersize/2:filtersize/2
		response = Run_one_rds(itterations, i, subunits, filtersize);
		allresps(:,:,n) = response;
		n = n + 1;
	end

	save convolutions allresps;
