function distance = DistanceToCluster(C, pos);

if isempty(C) | ~isfield(C,'x');
distance = NaN;
return;
end
xy = pos - [C.x(1) C.y(1)];
xy = xy ./ [C.x(3) C.y(3)];
cn = cos(-C.angle);
sn = sin(-C.angle);
p(1) = xy(1) * cn + xy(2) * sn;
p(2) = xy(2) * cn - xy(1) * sn;

distance = p./[C.x(2)./C.x(3) C.y(2)./C.y(3)];
distance = sum(distance.^2);



