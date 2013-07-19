function x = angle3d(a,b, varargin)
%
%x = angle3d(a,b)
%calculate the angle between two vectors
%



ma = sqrt(sum(a.^2));
mb = sqrt(sum(b.^2));
x = acos(dot(a,b)./(ma *mb));
