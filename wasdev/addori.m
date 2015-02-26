function fig = plotgrid(file)

OR =  dlmread(file,' ');
[rows, cols] = size(OR);
A = OR(:,2);
R = OR(:,1) * pi/180.0;

[XX, YY] = pol2cart(R,A);
%polar(R, A);
plot(XX, YY);
