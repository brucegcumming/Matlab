function XY = xyrotate(x,y,angle,varargin)

cosa = cos(angle);
sina = sin(angle);
x = reshape(x,length(x),1);
y = reshape(y,length(y),1);
xy = cat(2,x,y);
XY = xy * [cosa sina; -sina cosa];
