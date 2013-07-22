%function r = vhocorr()
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

sf = 1.3;
sdh = 0.15;
%sdh = 0.25;
sdv = 0.25;

asymfilter =  [0 0 sf 0    0 sdh sdv 0.5 0.0];%defines params for this Gabor
afiltershape = GFilter2 (xsize,ysize,asymfilter);%routine to create Gabor
afiltershape = afiltershape -mean(mean(afiltershape));

%fhandle = figure; %create new figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

afourtr = fft2 (afiltershape);	%2d fast fourier transform
ot = abs(fftshift(afourtr));
fourtr_ccorr = afourtr.*conj(afourtr);
aucorr = ifft2 (fourtr_ccorr);		%inverse fft
aucorr = real (aucorr);			%sometimes v.small complex 
					%parts left from roundoff
aucorr = fftshift (aucorr);		%shift so zero displacement

hold off;					
h = pcolor(aucorr);			%draw auto correlation
colormap('hot');
set(h,'EdgeColor', 'none');
rotate(h,[0 90],45);
hold on;
[C , h] = contour(aucorr,5);
set(h,'edgeColor', 'w','LineWidth',2);
rotate(h,[0 90],45);
axis('square');


set(gca,'fontsize',18);
set(gca,'Ytick',[]);
set(gca,'Xtick',[]);

arr = arrow([32 107], [45 120]);
text(-13,62,'Parallel Disparity','Rotation', 45, 'fontsize', 18);
arr = arrow([110 65], [125 50]);
text(55,120,'Orthogonal Disparity','Rotation', -45, 'fontsize', 18);
xlabel('Horizontal Disparity');
ylabel('Vertical Disparity');
axis('image');

x = 1:xsize;
y = 1:ysize;
[x2d,y2d]=meshgrid(x,y);

angles = [];
resps = [];
freq= sf * 100/pi;
freq =  sf * pi/100;
for angle=[0:0.1:2 * pi]
  xprime = x2d.*  cos(angle) - y2d.*sin(angle);
  grating = exp(2 * i *  pi .* freq .* xprime);
  resp = abs(sum(sum(grating.*afiltershape)));
  angles = [angles angle];
  resps = [resps resp];
  end
  angles = angles - pi/4;
resps = resps .^2 .* 0.0037;
[xx, yy] = pol2cart(angles, resps);
hold on;
plot(xx +50,yy +50,'linewidth',2);
hold off;









