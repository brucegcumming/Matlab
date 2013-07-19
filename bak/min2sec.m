function y = min2sec(x,varargin)


y = floor(x) * 60  + rem(x,1) * 100;