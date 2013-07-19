function y = TemporalImpulse(varargin)

 x = -100:100;
 y = Gauss([-50 10],x);
 z = Gauss([0 20],x);
y = y-z;
