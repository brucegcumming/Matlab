function gasprice(varargin)
%when to fill up with gas coming back from chincoteage
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'mpg')
    j = j+1;
    mpg = varargin{j};
end
end

G = 0:5
D=200; %journey in miles
mpg=22;
d = 1:D;
gasprice = 400 + d.*0.25;
for j = 1:length(d)
cost(1,j) = 450 .*(D-d(j))./mpg + (d(j) * gasprice(d(j)))./mpg;
end
rcost(1,:) = cost(1,:)-min(cost(1,:));
for j = 2:length(G)
cost(j,:) = cost(1,:) + G(j) .* (gasprice) ;
rcost(j,:)  = cost(j,:) - min(cost(j,:));
end
plot(d,rcost);
