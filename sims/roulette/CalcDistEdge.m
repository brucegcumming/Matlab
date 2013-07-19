function [e,r,p] = CalcDistEdge(x,y, varargin)


betw = 1;
angles = mod(x+pi,2*pi) - pi;
angles = atan2(sin(x),cos(x));
win = abs(angles) <  2 * (pi/37) * betw/2;
winsum = sum(y(find(win)));
loss = sum(y(find(~win)));
p = sum(y(find(win)))/sum(y);
e = (p * 36 - betw)/betw;

%e = ((winsum * (36-betw)) - (loss * betw))/((winsum+loss) * betw);

v = y .* sin(x) + i * y .* cos(x);
r = abs(sum(v))/sum(y);
