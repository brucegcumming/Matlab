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
sdv = 0.1;

asymfilter =  [0 0 sf 0    0 sdh sdv 0.5 0.0];%defines params for this Gabor
bsymfilter =  [0 0 sf -1.0  0 sdh sdv 0.5 0.0];
csymfilter =  [0 0 sf -0.78  0 sdh sdv 0.5 0.0];
hsymfilter =  [0 0 sf 0    90 sdh sdv 0.5 0.0];	
vsymfilter =  [0 0 sf 0    0 sdh sdv 0.5 0.0];	
barsymfilter =  [0 0 sf -1.57    45 sdh 20 0.5 0.0];
barasymfilter =  [0 0 sf 0    45 sdh 20 0.5 0.0];
afiltershape = GFilter2 (xsize,ysize,asymfilter);%routine to create Gabor
afiltershape = afiltershape -mean(mean(afiltershape));
bfiltershape = GFilter2 (xsize,ysize,bsymfilter);
bfiltershape = bfiltershape -mean(mean(bfiltershape));
cfiltershape = GFilter2 (xsize,ysize,csymfilter);
cfiltershape = cfiltershape -mean(mean(cfiltershape));
fhandle = figure; %create new figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

afourtr = fft2 (afiltershape);	%2d fast fourier transform
bfourtr = fft2 (bfiltershape);	%2d fast fourier transform
cfourtr = fft2 (cfiltershape);	%2d fast fourier transform


fourtr_ccorr = afourtr.*conj(bfourtr);
aucorr = ifft2 (fourtr_ccorr);		%inverse fft
aucorr = real (aucorr);			%sometimes v.small complex 
					%parts left from roundoff
aucorr = fftshift (aucorr);		%shift so zero displacement is
subplot (2,1,1);	
draw (aucorr,-99,-99);			%draw auto correlation
set(gca,'Ytick',[]);
set(gca,'Xtick',[]);
line([ 50 50], [0 100], 'Color', 'black');
line([ 0 100], [50 50], 'Color', 'black');
ylabel('Vertical Disparity');
title ('Phase shift 1.0');
y=aucorr(50,:);

fourtr_ccorr = afourtr.*conj(cfourtr);
aucorr = ifft2 (fourtr_ccorr);		%inverse fft
aucorr = real (aucorr);			%sometimes v.small complex 
					%parts left from roundoff
aucorr = fftshift (aucorr);		%shift so zero displacement is
z=aucorr(50,:);
x=(1:length(y));
subplot (2,1,2);	
plot(x,y,x,z);		%plot cross section

hold on
fourtr_ccorr = bfourtr.*conj(cfourtr);
aucorr = ifft2 (fourtr_ccorr);		%inverse fft
aucorr = real (aucorr);			%sometimes v.small complex 
					%parts left from roundoff
aucorr = fftshift (aucorr);		%shift so zero displacement is
z=aucorr(50,:);
plot(x,z);
hold off

%r = aucorr;				%return autocorrelation







