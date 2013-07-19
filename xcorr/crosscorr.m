function r = crosscorr()
%Makes 2d Gabor patch and performs 2d cross-correlation
%in the Fourier domain..

xsize = 100;
ysize = 100;


%Gabor Parameters.
%1  horz offset
%2  vert offset
%3  frequency in cpd.
%4  phase relative to mean position
%5  orientation clockwise from horizontal
%6  sd perpendicular to bars
%7  sd parallel to bars
%8  peak amplitude
%9  mean 

sf = 2;
sdh = 0.2;
sdv = 0.2;

asymfilter =  [0 0 sf 0    45 sdh sdv 0.5 0.0];%defines params for this Gabor
hsymfilter =  [0 0 sf 0    90 sdh sdv 0.5 0.0];	
bsymfilter =  [0.176 0 sf 0    45 sdh sdv 0.5 0.0];
csymfilter =  [0 0 sf -1.57    45 sdh sdv 0.5 0.0];
afiltershape = GFilter2 (xsize,ysize,asymfilter);%routine to create Gabor
afiltershape = afiltershape -mean(mean(afiltershape));
bfiltershape = GFilter2 (xsize,ysize,bsymfilter);
bfiltershape = bfiltershape -mean(mean(bfiltershape));
cfiltershape = GFilter2 (xsize,ysize,csymfilter);
cfiltershape = cfiltershape -mean(mean(cfiltershape));
hfiltershape = GFilter2 (xsize,ysize,hsymfilter);
hfiltershape = hfiltershape -mean(mean(hfiltershape));
fhandle = figure; %create new figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

afourtr = fft2 (afiltershape);						%2d fast fourier transform
hfourtr = fft2 (hfiltershape);						%2d fast fourier transform
bfourtr = fft2 (bfiltershape);						%2d fast fourier transform
cfourtr = fft2 (cfiltershape);						%2d fast fourier transform


fourtr_ccorr = hfourtr.*conj(hfourtr);
aucorr = ifft2 (fourtr_ccorr);		%inverse fft
aucorr = real (aucorr);			%sometimes v.small complex 
					%parts left from roundoff
aucorr = fftshift (aucorr);		%shift so zero displacement is
subplot (1,3,1);	
draw (aucorr,-99,-99);			%draw auto correlation
line([ 50 50], [0 100], 'Color', 'black');
line([ 0 100], [50 50], 'Color', 'black');
ylabel('Vertical Disparity');
title ('Zero Disparity');


fourtr_ccorr = afourtr.*conj(bfourtr);
aucorr = ifft2 (fourtr_ccorr);		%inverse fft
aucorr = real (aucorr);			%sometimes v.small complex 
					%parts left from roundoff
aucorr = fftshift (aucorr);		%shift so zero displacement is
subplot (1,3,2);	
draw (aucorr,-99,-99);			%draw auto correlation
line([ 50 50], [0 100], 'Color', 'black');
line([ 0 100], [50 50], 'Color', 'black');
xlabel('Horizontal Disparity');
title ('Postion Disparity');

fourtr_ccorr = afourtr.*conj(cfourtr);	%multiply by its own complex conjugate
aucorr = ifft2 (fourtr_ccorr);		%inverse fft
aucorr = real (aucorr);			%sometimes v.small complex 
					%parts left from roundoff
aucorr = fftshift (aucorr);		%shift so zero displacement is
subplot (1,3,3);	
draw (aucorr,-99,-99);							%draw auto correlation
line([ 50 50], [0 100], 'Color', 'black');
line([ 0 100], [50 50], 'Color', 'black');
title ('Phase Disparity');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%r = aucorr;								%return autocorrelation







