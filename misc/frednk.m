function [contrasts, sumresps] = fred(varargin)


colors = mycolors;
j= 1;
maxc = 1;
expo = 2; 
pwr = 2;
cfifty = 0.5;
model = 1;
power = 4;


while(j < nargin+1)
    if(strncmpi(varargin{j},'maxc',4))
        j = j+1;
        maxc = varargin{j};
    elseif strncmpi(varargin{j},'test',4)
        atest(expo);
        return;
    elseif(strncmpi(varargin{j},'exp',3))
        j = j+1;
        expo = varargin{j};
    elseif(strncmpi(varargin{j},'model',3))
        j = j+1;
        model = varargin{j};
    elseif(strncmpi(varargin{j},'c50',3))
        j = j+1;
        cfifty = varargin{j};
    end
    j = j+1;
end


contrast = 0.01;
inc = 1.1;
contrasts = [];
while contrast < 1;
    contrasts = [contrasts contrast];
    contrast = contrast * inc;
end

lc = log(contrasts) - log(contrasts(1));
aresp = log(contrasts) - log(contrasts(1));
id = find(contrasts > maxc);
if ~isempty(id)
    aresp(id) = deal(min(aresp(id)));
end

aresp = contrasts .^ expo./(contrasts.^expo + cfifty.^expo);
bresp = -1 * aresp;

again = aresp(end);
bgain = bresp(end);
subplot(2,1,2);
hold off;
subplot(2,1,1);
hold off;
plot(contrasts,aresp);
hold on;
plot(contrasts,bresp,'r');
%plot(contrasts, contrasts.^expo .* bgain, 'r');
%plot(contrasts, contrasts.^expo .* again, 'b');
set(gca,'Xscale','log');
hold on;


%for j = 2:length(contrasts)
%The negative response is driven by testc, held constant.
%The contract of the component giving the positive response is varied.

j = 1;
for testc = [0.2 0.35 0.5 0.65 0.8]
ratios = contrasts/testc;
subplot(2,1,1);


%
% Model 1  linearly sum responses, but divisor is unsigned. I.e. both
% components noralize.
if model == 1
    sumresp = (contrasts .^ expo - testc.^expo) ./ ((contrasts.^expo + testc.^expo) + cfifty.^expo);
%Model 2 straightforward summation of responses
elseif model == 2
    sumresp = (contrasts .^ expo) ./ ((contrasts.^expo) + cfifty.^expo) - testc.^expo/(testc^expo + cfifty^expo);
%Model 3  sum as nrt(a^n + b^n)
% i.e. when n = 2, as if orthogonal vectors
elseif model == 3
    sumrespb = ((contrasts .^ expo) ./ ((contrasts.^expo) + cfifty.^expo)).^power - (testc.^expo/(testc^expo + cfifty^expo))^power;
    sumresp = abs(sumrespb).^(1/power) .* sign(sumrespb);
% Model 4 and 5, winner takes all style. In model 5, weight is 0 or 1,
% accoring to relative contrast (= real winner takes all). In model 4, the
% weigths are taken from raising the ratio to a power, giving a sigmoid
% weight function.
elseif model == 4 | model == 5
    scales = ratios.^10;
    if model == 5
        scales(find(ratios < 1)) = 0;
        scales(find(ratios > 1)) = 1;
    elseif model == 4
        scales = scales ./ (1 + scales);
    end
    sumresp = (contrasts .^ expo) ./ ((contrasts.^expo) + cfifty.^expo) .* scales - (testc.^expo/(testc^expo + cfifty^expo)) .* (1- scales);
end
    plot(contrasts,sumresp,'g');
    subplot(2,1,2);
    plot(ratios, sumresp, 'color',colors{j});
    hold on;
j = j+1;
end
set(gca,'Xscale','log');

function atest(expo)

lmin = 1;
x = 0.01:0.01:0.99;
lmax = lmin * (1+x)./(1-x);
a = 10;
b = -10;

r = a.*(x.^expo) ./ (x.^expo + 2);
s = (a.*x.^expo + b .* x(end) .^expo) ./ (x.^expo + x(end).^expo+2);

hold off;
plot(x,r);

hold on;
plot(x,s,'r');
set(gca,'Xscale','log');



