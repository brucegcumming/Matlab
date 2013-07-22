function answer=bruceanalyze(trials, response)
%	
%	answer=analyze(trials)
%	
%	gets binocular responses from monocular inputs
%
%
	complex = 0;
	acomplex = 0;
	bcomplex = 0;
%setting thresh to a large negative number restores the
%original ODF model
	thresh = 0;
%	thresh = -100000;
	
% response = [left a, left b, right a, right b]

	for i = 1:trials
	  
%rectfied monocular inputs
          la = max([thresh, response(i,1)]);
	  lb = max([thresh, response(i,2)]);
	  ra = max([thresh, response(i,3)]);
	  rb = max([thresh, response(i,4)]);

%rectfied monocular inputs for "off" halves
	  lan = max([thresh, -1 * response(i,1)]);
	  lbn = max([thresh, -1 * response(i,2)]);
	  ran = max([thresh, -1 * response(i,3)]);
	  rbn = max([thresh, -1 * response(i,4)]);
	  asum = la + ra;
	  asumn = lan + ran;
	  bsum = lb+rb;
	  bsumn = lbn +rbn;
	  acomplex = acomplex + (asum^2 + asumn^2);
	  bcomplex = bcomplex + (bsum^2 + bsumn^2);
	  complex = complex + (asum^2 + asumn^2 + bsum^2 + bsumn^2);
	end
	answer = [complex/trials acomplex/trials bcomplex/trials];

