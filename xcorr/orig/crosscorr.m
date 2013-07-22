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


dsymfilter =  [0 0 4 0    45 0.10 0.15 0.5 0.0];			%defines params for this Gabor
filtershape = GFilter2 (xsize,ysize,dsymfilter);			%routine to create Gabor
filtershape = filtershape -mean(mean(filtershape));			%remove any remaining DC.
fhandle = figure;							%create new figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1);								%use 2x2 subplot grid, posn1
draw(filtershape,-99,-99);						%calls MY routine draw
title ('Gabor Receptive Field');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2);
fourtr = fft2 (filtershape);						%2d fast fourier transform
draw(abs (fftshift(fourtr)),-99,-99);					%shifts it so dc is in middle
									%not corner and plots
title ('Fourier Transform');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fourtr_aucorr = fourtr.*conj(fourtr);					%multiply by its own complex
									%conjugate		
aucorr = ifft2 (fourtr_aucorr);						%inverse fft
aucorr = real (aucorr);							%sometimes v.small complex 
									%parts left from roundoff
subplot (2,2,3);	
draw (aucorr,-99,-99);							%draw auto correlation
title ('2d Autocorrelation function');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aucorr = fftshift (aucorr);						%shift so zero displacement is
									%in middle
subplot (2,2,4);							%draw shifted correlation
draw (aucorr,-99,-99);
title ('Shifted 2d Auto correlation function');



r = aucorr;								%return autocorrelation







