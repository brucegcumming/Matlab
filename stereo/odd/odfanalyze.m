function answer=odfanalyze(trials, response)
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
          la = response(i,1);
	  lb = response(i,2);
	  ra = response(i,3);
	  rb = response(i,4);

%rectfied monocular inputs for "off" halves
	  asum = la + ra;
	  bsum = lb+rb;
	  acomplex = acomplex + (asum^2);
	  bcomplex = bcomplex + (bsum^2);
	  complex = complex + (asum^2 + bsum^2);
	end
	answer = [complex/trials acomplex/trials bcomplex/trials];

