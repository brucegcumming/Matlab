function d = ArrayDistance(Array, x,y)



if isfield(Array,'Y')
    d = abs(diff(Array.X([x y]))+ ...
        i * diff(Array.Y([x y])));
else
    d = cid(j)-cid(k);
end
if isfield(Array,'separation')
    d = d * Array.separartion;
end