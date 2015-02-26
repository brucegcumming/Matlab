function o = xy2op(x, angle, varargin)
%
% xy2op(op, angle,..)
% converts a location in orhtog distance, paralell distance into x,y
% (screen co-ordingate), given an RF orientation (in degrees).
%
% 0 degrees points down (as in binoc), not left
% In binoc, if Ro = 90, Op+ -> -v X. Pp+ = +v y
%              Ro=-90 Op+ -> +v X. Pp+ = -v y
%              Ro=0 Op+ -> +v Y. Pp+ = +v x
%              Ro=180 Op+ -> -v Y. Pp+ = -v x
%
theta = -pi/2 + angle * pi/180;
j = 1;
while j < nargin-1;
    if strncmpi(varargin{j},'rad',3)
        theta = -pi/2 + angle;
    end
    j = j+1;
end
    ca = cos(theta);
    sa = sin(theta);
% need to recheck sign conventions here in replay...    

mx = [-cos(theta) -sin(theta); -sin(theta) cos(theta)];
if size(x,1) == 2
    o = x' * mx';
else
    o = x * mx;
end
