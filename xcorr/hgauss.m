%function r = hgauss()
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

dsymfilter =  [0 0 0.01 0    0 0.6 0.2 0.5 0.0];%defines params for this Gabor

dfiltershape = GFilter2 (xsize,ysize,dsymfilter);%routine to create Gabor
dfiltershape = dfiltershape' -mean(mean(dfiltershape));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



colormap('hot');
output = imagesc (dfiltershape, [-max(max(dfiltershape)), ...
		    max(max(dfiltershape))]);
axis('image');
colormap('hot');
set(gca,'Ytick',[]);
set(gca,'Xtick',[]);








