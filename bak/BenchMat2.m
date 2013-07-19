function [time1 time2] = BenchMat2(n, double)

% matlab benchmark meant to resemble LSR/Vision work
% returned number should be close to 1.2773

if nargin==1
    double=1;
end

if double
    tic;

    m1a=zeros(n);
    m2a=zeros(n);
    m4a=zeros(10,n,n/10);
    a1a=zeros(n,10);

    ta_alloc=toc;

    tic;

    m1a=randn(n);
    m2a=randn(n);
    m3a=m1a*m2a;

    for i=1:10
        m3a=sqrt(exp(sin(atan(m3a)./pi)));
    end

    ta_vector=toc;

    tic;

    for i=1:n/10
        for j=1:n
            for k=1:10
                m4a(k,j,i)=sqrt(exp(sin(mod(m3a(i*9+mod(j,10)-k+2,mod(j+i-k+1,n)+1),1))));
            end
        end
        a1a(i*10,:)=mean(reshape(m4a,n*n/10,10));
    end

    result_calculation=mean(a1a(find(a1a>0)));

    ta_loop = toc;
end

tic;

m1b=zeros(n);
m2b=zeros(n);
m4b=zeros(10,n,n/10);
a1b=zeros(n,10);

    for i=1:10000
        nn=fix(n/100+sin(i));
        mem1{i}=zeros(1,nn*nn-1);
        mem2{i}=zeros(nn);
        mem3{i}=zeros(fix(nn/10),fix(nn/10),i);
        mem1{i}=deal(sin(i));
        mem2{i}=deal(sin(i+1));
        mem3{i}=deal(sin(i+2));
    end
    
clear mem1; clear mem2; clear mem3;

tb_alloc=toc

tic;

m1b=randn(n);
m2b=randn(n);
m3b=m1b*m2b;

for i=1:10
    m3b=sqrt(exp(sin(atan(m3b)./pi)));
end

tb_vector=toc

tic;

for i=1:n/10
    for j=1:n
        for k=1:10
            m4b(k,j,i)=sqrt(exp(sin(mod(m3b(i*9+mod(j,10)-k+2,mod(j+i-k+1,n)+1),1))));
        end
    end
    a1b(i*10,:)=mean(reshape(m4b,n*n/10,10));
end

tmp=mean(a1b(find(a1b>0)));

tb_loop=toc

result_calculation=tmp

if double
    if abs((ta_alloc/tb_alloc)-1)>0.2
        warning('time diff between two runs of memory test greater than 20% in both runs');
    end
    if abs((ta_vector/tb_vector)-1)>0.2
        warning('time diff between two runs of vector computation test greater than 20% in both runs');
    end
    if abs((ta_loop/tb_loop)-1)>0.2
        warning('time diff between two runs of loop computation test greater than 20% in both runs');
    end
end
