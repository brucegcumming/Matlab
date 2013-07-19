
nreps = 8;
hold off;
for ysd = [1 4 16];
xsds = [];
means = [];
for xsd = 0:5:50;
    answers = [];
for j = 1:1000
    sd = 20;
x = [-100:10:100];
X = repmat(x,nreps,1);
gp = [0 20 40];
guess = [0 20 40];
noise = randn(nreps,length(x)) * xsd;
X = X + noise;
Y = gauss(gp, X);
noise = randn(nreps,length(x)) * ysd;
Y = Y + noise;
y = sum(Y) ./ nreps;
answer = nlinfit(x,y,'gauss',guess);
fity = gauss(answer, x);
guessy = gauss(guess, x);
answers = [answers, answer];
end
%plot(x,y);
%hold on;
%plot(x,fity,'r');
drawnow;
xsds = [xsds; xsd];
means = [means; sum(answers') ./ length(answers')];
end
means
plot(xsds,means(:,2));
hold on
plot(xsds,means(:,3),'r');
pred = sqrt((xsds .* xsds) + (20 * 20));
end
plot(xsds,pred,'g');
