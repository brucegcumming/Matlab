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

asymfilter =  [0 0 sf  0 90 sdh sdv 0.5 0.0];%defines params for this Gabor
bsymfilter =  [0 0 sf 0 90 sdh 0.6 0.5 0.0];	
csymfilter =  [0 0 sf  0 90 sdh 0.1 0.5 0.0];	
dsymfilter =  [0 0 0.6  0 45 sdh 0.6 0.5 0.0];
esymfilter =  [0 0 sf  0 90 sdh sdv 0.5 0.0];
fsymfilter =  [0 0 0.2  0 90 sdh sdv 0.5 0.0];
barsymfilter =  [0 0 sf -1.57    45 sdh 20 0.5 0.0];
barasymfilter =  [0 0 sf 0    45 sdh 20 0.5 0.0];
afiltershape = GFilter2 (xsize,ysize,asymfilter);%routine to create Gabor
afiltershape = afiltershape -mean(mean(afiltershape));
bfiltershape = GFilter2 (xsize,ysize,bsymfilter);
bfiltershape = bfiltershape -mean(mean(bfiltershape));
cfiltershape = GFilter2 (xsize,ysize,csymfilter);
cfiltershape = cfiltershape -mean(mean(cfiltershape));
dfiltershape = GFilter2 (xsize,ysize,dsymfilter);
dfiltershape = dfiltershape -mean(mean(dfiltershape));
efiltershape = GFilter2 (xsize,ysize,esymfilter);
efiltershape = efiltershape -mean(mean(efiltershape));
fhandle = figure; %create new figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

afourtr = fft2 (afiltershape);	%2d fast fourier transform
bfourtr = fft2 (bfiltershape);	%2d fast fourier transform
cfourtr = fft2 (cfiltershape);	%2d fast fourier transform
dfourtr = fft2 (dfiltershape);	%2d fast fourier transform
efourtr = fft2 (efiltershape);	%2d fast fourier transform
afourtr = fftshift(real(afourtr));
bfourtr = fftshift(real(bfourtr));
cfourtr = fftshift(real(cfourtr));
dfourtr = fftshift(real(dfourtr));
efourtr = fftshift(real(efourtr));

subplot (2,4,1);	
draw (afiltershape,-99,-99);	
subplot (2,4,5);	
draw (afourtr,-99,-99);		

subplot (2,4,2);	
draw (bfiltershape,-99,-99);	
%contour(bfiltershape);
subplot (2,4,6);	
draw (bfourtr,-99,-99);		
subplot (2,4,3);	
draw (cfiltershape,-99,-99);	
subplot (2,4,7);	
draw (cfourtr,-99,-99);		

subplot (2,4,4);	
draw (dfiltershape,-99,-99);
subplot (2,4,8);	
draw (dfourtr,-99,-99);		








