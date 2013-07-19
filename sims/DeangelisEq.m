function DeangleisEqn(b, varargin);

x = -1:0.1:1;
thcrit = 0.75;
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'xrange',5)
    j = j+1;
    x = varargin{j};
end
j = j+1;
end

thk = -log(1/thcrit -1);
if nargin == 0
b(1) = 0;
b(2) = 5;
end
y = 1./(1 + exp(-(b(1) +  b(2) .* x)));
plot(x,y);

t = thk./b(2);
hold on;
plot([x(1) t t],[0.75 0.75 0],'r-')

slope = (y(101)-y(99))./(x(101)-x(99))