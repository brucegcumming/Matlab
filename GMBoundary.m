function XY  = GMBoundary(C, varargin)
showplot = 0;

j = 0;
while j <= length(varargin)
    if strncmpi(varargin{j},'plot',4)
        showplot = 1;
    end
    j = j+1;
end

xi = linspace(min(C.xy(:,1)),max(C.xy(:,1)));
yi = linspace(min(C.xy(:,2)),max(C.xy(:,2)));
[x,y] = meshgrid(xi,yi);

xy = cat(2,x(:),y(:));
X = posterior(C.gmfit2d,xy);
Xd = reshape(X(:,1),size(x));
Yd = reshape(X(:,2),size(x));
Z = abs(Xd-Yd);
[a,b] = min(Z');
if showplot
imagesc(minmax(C.xy(:,1)),minmax(C.xy(:,2)),Z);
hold on;
plot(xi(b),yi);
end
XY(:,1) = xi(b);
XY(:,2) = yi;
