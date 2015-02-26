function Fn = FigureProp(F, prop, varargin)
%return figure number indepdent of matlab version
if nargin < 2
    prop = 'Number';
end
if isobject(F)
    Fn = get(F, prop)
else
    Fn = F;
end
