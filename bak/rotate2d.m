function o = rotate2d(x, theta,varargin)

j = 1;
while j < nargin-1;
    if strncmpi(varargin{j},'deg',3)
        theta = theta * pi/180;
    end
    j = j+1;
end


mx = [cos(theta) sin(theta); -sin(theta) cos(theta)];
o = x * mx;


