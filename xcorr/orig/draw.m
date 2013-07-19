function r = draw(imagedata,fixmin,fixmax)
%function r = draw(imagedata,fixmin,fixmax)
%draws image in a window, with a scale
%if fix = -99 then it fixes the zero level to min and rescales
%the rest  by the data. Else it fixes at the given level






imagedata = imagedata';
max_value = max (max(imagedata));
min_value = min (min(imagedata));

if (fixmin==-99)
	fixmin = min_value;
end;

if (fixmax==-99)
	fixmax = max_value;
end;

colordata = [0:0.01:1 ; 0:0.01:1; 0:0.01:1]';
set (gcf,'colormap',colordata);



output = imagesc (imagedata,[fixmin fixmax]);
axis square;
axis off;

%colorbar('vert');




