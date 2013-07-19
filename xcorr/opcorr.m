function r = opcorr()
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

asymfilter =  [0 0 sf 0    -45 sdh sdv 0.5 0.0];%defines params for this Gabor
hsymfilter =  [0 0 sf 0    90 sdh sdv 0.5 0.0];	
vsymfilter =  [0 0 sf 0    0 sdh sdv 0.5 0.0];	
bsymfilter =  [0.176 0 sf 0    45 sdh sdv 0.5 0.0];
csymfilter =  [0 0 sf -1.57    45 sdh sdv 0.5 0.0];
barsymfilter =  [0 0 sf -1.57    45 sdh 20 0.5 0.0];
barasymfilter =  [0 0 sf 0    45 sdh 20 0.5 0.0];
afiltershape = GFilter2 (xsize,ysize,asymfilter);%routine to create Gabor
afiltershape = afiltershape -mean(mean(afiltershape));
bfiltershape = GFilter2 (xsize,ysize,bsymfilter);
bfiltershape = bfiltershape -mean(mean(bfiltershape));
cfiltershape = GFilter2 (xsize,ysize,csymfilter);
cfiltershape = cfiltershape -mean(mean(cfiltershape));
hfiltershape = GFilter2 (xsize,ysize,hsymfilter);
hfiltershape = hfiltershape -mean(mean(hfiltershape));
vfiltershape = GFilter2 (xsize,ysize,vsymfilter);
vfiltershape = vfiltershape -mean(mean(vfiltershape));
barfiltershape = GFilter2 (xsize*2,ysize*2,barsymfilter);
barfiltershape = barfiltershape -mean(mean(barfiltershape));
abarfiltershape = GFilter2 (xsize*2,ysize*2,barasymfilter);
abarfiltershape = abarfiltershape -mean(mean(abarfiltershape));
fhandle = figure; %create new figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

afourtr = fft2 (afiltershape);	%2d fast fourier transform
hfourtr = fft2 (hfiltershape);	%2d fast fourier transform
vfourtr = fft2 (vfiltershape);	%2d fast fourier transform
bfourtr = fft2 (bfiltershape);	%2d fast fourier transform
cfourtr = fft2 (cfiltershape);	%2d fast fourier transform
ofourtr = fft2 (barfiltershape);	%2d fast fourier transform
oafourtr = fft2 (abarfiltershape);	%2d fast fourier transform


fourtr_ccorr = vfourtr.*conj(vfourtr);
aucorr = ifft2 (fourtr_ccorr);		%inverse fft
aucorr = real (aucorr);			%sometimes v.small complex 
					%parts left from roundoff
aucorr = fftshift (aucorr);		%shift so zero displacement is

fourtr_ccorr = afourtr.*conj(afourtr);
aucorr = ifft2 (fourtr_ccorr);		%inverse fft
aucorr = real (aucorr);			%sometimes v.small complex 
					%parts left from roundoff
aucorr = fftshift (aucorr);		%shift so zero displacement is
subplot (1,2,1);	
draw (aucorr,-99,-99);			%draw auto correlation
set(gca,'Ytick',[]);
set(gca,'Xtick',[]);

fourtr_ccorr = vfourtr.*conj(vfourtr);
aucorr = ifft2 (fourtr_ccorr);		%inverse fft
aucorr = real (aucorr);			%sometimes v.small complex 
					%parts left from roundoff
aucorr = fftshift (aucorr);		%shift so zero displacement is
subplot (1,2,2);	
draw (aucorr,-99,-99);			%draw auto correlation
set(gca,'Ytick',[]);
set(gca,'Xtick',[]);









