function logrecip(N)
%
% what does y = 1/r look like on log and lin axes
 
 x = 0.1:0.1:30;
 
 
 x = repmat(x,1,20);
 y = 1./( x) + rand(size(x));
 plot(x(:),y(:),'o');

