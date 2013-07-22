function cp = GrandCP(ratio)

nl = 10
nt = 1000;
x = randn(2,nt);
y = randn(2,nt);
choices(1,:) = x(1,:)+y(1,:) > prctile(x(:),50);
choices(2,:) = x(2,:)+y(2,:) > prctile(x(:),50);

for j = 1:nl
    x = randn(2,nt);
    y = randn(2,nt);
    choices(1,:) = x(1,:)+y(1,:) > prctile(x(:),100*ratio);
    choices(2,:) = x(2,:)+y(2,:) > prctile(x(:),100*(1-ratio));
    cp(j,1) = CalcCP(x(1,choices(1,:) ==1), x(1,choices(1,:) ==0));
    cp(j,2) = CalcCP(x(2,choices(2,:) ==1), x(2,choices(2,:) ==0));
    cp(j,3) = CalcCP(x(choices ==1), x(choices ==0));
end
