function XY = xyrotate(x,y,angle,varargin)
%XY = xyrotate(x,y,angle)

if min(size(x)) == 2 && length(y) == 1
    angle = y;
    xy = x;
else
    x = reshape(x,length(x),1);
    y = reshape(y,length(y),1);
    xy = cat(2,x,y);
end
cosa = cos(angle);
sina = sin(angle);
XY = xy * [cosa sina; -sina cosa];
