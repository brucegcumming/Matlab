function [x, y] = circle(r,c, varargin)
%h = circle(r, c, args) plots a circle, with args{:} passed on to plot()
% retruns handle to the plotted function
%[x,y] = circle(r,c) (> 1 output argument) returns co-ordinates of the
%circumference, no plotting
% now superceded by ellipse.

j = [-pi:pi/100:pi];
if nargin == 1
    c = [0 0];
end
x = r.*cos(j) + c(1); 
y = r.*sin(j) + c(2); 
if nargout > 1
    return;
end
h = plot(x,y,varargin{:});
