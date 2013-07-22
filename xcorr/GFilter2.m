function r = GFilter2 (xdims,ydims,Gabor_Params)
%function r = GFilter2 (xdims,ydims,Gabor_Params)
%produces Gabor filter that can then be convolved defined in Gabor params
%60 pixels = 1 degree
%Gabor Params is a vector of components defining 2d Gabor parameters
%1  horz offset
%2  vert offset
%3  frequency in cpd.
%4  phase relative to mean position
%5  orientation clockwise from horizontal
%6  sd perpendicular to bars
%7  sd parallel to bars
%8  peak amplitude
%9  mean 

scaling_factor = 60;  %100


xposn = Gabor_Params(1);
yposn = Gabor_Params(2);
freq =  Gabor_Params(3);
phase = Gabor_Params(4);
or    = Gabor_Params(5);
sdx   = Gabor_Params(6);
sdy   = Gabor_Params(7);
amp   = Gabor_Params(8);
mean  = Gabor_Params(9);

%scaling parameters

xposn = xposn *scaling_factor;
yposn = yposn *scaling_factor;
freq = freq/scaling_factor;
sdx = sdx *scaling_factor;
sdy = sdy *scaling_factor;



%define arrays of x and y positions for convenience.

x= (1:xdims)';
x= x-xposn;
x2 = zeros (xdims,ydims);
for (count = 1:ydims);
	x2 (:,count) = x; 
end;

y= 1:ydims;
y = y+yposn;
y2 = zeros (xdims,ydims);
for (count1=1:xdims)
	y2(count1,:) = y;
end	


%centre x and y positions

x2 = x2- (xdims+1)/2;
y2 = y2- (ydims+1)/2;

%convert X2 and Y2 to dimensions of orientatation

or = or *pi/180;
x3 = cos (or).*x2 - sin (or).*y2;
y3 = sin (or).*x2 + cos (or).*y2;

% Calculate Filter

filter = cos (x3.*freq*2*pi + phase);
%filter = filter.*exp(-((x3-xposn).^2) / (2*sdx^2));
%filter = filter.*exp(-((y3-yposn).^2) / (2*sdy^2));
filter = filter.*exp(-((x3).^2) / (2*sdx^2));
filter = filter.*exp(-((y3).^2) / (2*sdy^2));

filter = filter * amp;
filter = filter +  mean;
filter;

%draw(filter,mean-amp,mean+amp);

r = filter;

%>> GFilter(temp,[0 0 5 0 0 0.5 0.5 1 1])
%GFilter(temp,[0 0 5 0 30 0.5 0.5 1 1]) 
