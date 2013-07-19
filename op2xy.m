function o = op2xy(x, angle,varargin)
%
% op2xy(op, angle,..)
% converts a location in orhtog distance, paralell distance into x,y
% (screen co-ordingate), given an RF orientation (in degrees).
%
% 0 degrees points down as in binoc (not right, as is common).

theta = -pi/2 + angle * pi/180;
j = 1;
while j < nargin-1;
    if strncmpi(varargin{j},'rad',3)
        theta = -pi/2 + angle;
    end
    j = j+1;
end



mx = [-cos(theta) -sin(theta); -sin(theta) cos(theta)];
o = x * mx;
