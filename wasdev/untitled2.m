lambda=[1/6 2/6 1/6 1/6 1/6 0/6];
x=[1 2 3 4 5 6];
p_x=[1/6 1/6 1/6 1/6 1/6 1/6];
for i=1:10^5
    oldexpval=sum(x.*lambda);
    q_star=p_x*exp(oldexpval);
    normq=q_star/sum(q_star);
    newexpval=sum(x.*normq);
    epsilon=10^-5;
    lambda=lambda-epsilon*(newexpval-3.2);
end

    