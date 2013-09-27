function o = op2xy(x, angle,varargin)
%
% op2xy(op, angle,..)
% converts a location in orhtog distance, paralell distance into x,y
% (screen co-ordingate), given an RF orientation (in degrees).
%
% 0 degrees points down as in binoc (not right, as is common).
% In binoc, if Ro = 90, Op+ -> -v X. Pp+ = +v y
%              Ro=-90 Op+ -> +v X. Pp+ = -v y
%              Ro=0 Op+ -> +v Y. Pp+ = +v x
%              Ro=180 Op+ -> -v Y. Pp+ = -v x
%
mirror = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'mirror',5)
        mirror = 1;
    end
    j = j+1;
end
theta = -pi/2 + angle * pi/180;
j = 1;
while j < nargin-1;
    if strncmpi(varargin{j},'rad',3)
        theta = -pi/2 + angle;
    end
    j = j+1;
end



mx = [-cos(theta) -sin(theta); -sin(theta) cos(theta)];
if mirror
    mx = [-cos(theta) -sin(theta); -sin(theta) cos(theta)];
end
if size(x,1) == 2
o = x' * mx';
else
o = x * mx;
end

