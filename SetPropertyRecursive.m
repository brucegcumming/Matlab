function SetPropertyRecursive(X, prop, value, varargin)
% SetPropertyRecursive(X, prop, ...
%itereates through all children of a graphics handle and sets property to
%value, if prop is a valid property
%SetPropertyRecursive(X, 'fontsize', 14)  sets all object that have a
%fontsize to fontsise

P = get(X);
if isfield(P,prop)
    set(X, prop, value);
end
if isfield(P,'Children')
    for j = 1:length(P.Children)
        SetPropertyRecursive(P.Children(j),prop,value);
    end
end