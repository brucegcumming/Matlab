

clc

num=input('Enter a  positive number: ');
num2=num;

while (num<0)
    clc
    fprintf('Invalid entry!\n');
    a=input('Enter a positive number: ');
end

a=1;
c=num;
d=0;

while (c>=1)
    b=rem(c,2);
    c=c/2;
    c=floor(c);
    d=d+(b*a);
    a=a*10;
end

lsb=rem(num,2);
msb=b;
    
fprintf('\n\n%g is written %g in binary form.\n\n',num2,d)

