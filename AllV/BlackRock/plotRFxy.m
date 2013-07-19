function PlotRFxy(rfX,rfY)



xs = [2:9 11:90 92:99];
xi = 2;
yi = 1;

[X,Y] = meshgrid(1:10,1:10);
Z(2:9) = rfX(1:8);
Z(11:90) = rfX(9:88);
Z(92:99) = rfX(89:96);
Z([1 10 91 100]) = NaN;
Z = reshape(Z,size(X));
pcolor(X,Y,Z);

X(1,1) = NaN;
X(1,10) = NaN;
X(1,10) = NaN;
X(10,1) = NaN;
Y(1,1) = NaN;
Y(1,10) = NaN;
Y(1,10) = NaN;
Y(10,1) = NaN;
