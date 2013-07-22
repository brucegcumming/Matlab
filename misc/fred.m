function fred(varargin)


colors = mycolors;
j= 1;
maxc = 1;
expo = 2; 
pwr = 2;
cfifty = 0.5;

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
    elseif(strncmpi(varargin{j},'c50',3))
        j = j+1;
        cfifty = varargin{j};
    end
    j = j+1;
end


contrast = 0.01;
inc = 1.5;
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

for j = 2:length(contrasts)
%sumresp = sign(aresp(j)) * aresp(j).^expo + sign(bresp) .* bresp.^expo;
%sumresp = sign(aresp(j)) * aresp(j) * contrasts(j).^expo + sign(bresp) .* bresp .* contrasts .^expo;
%sumresp =  again * contrasts(j).^expo + bgain * contrasts .^expo;
%sumresp = sumresp ./ (contrasts(j) .^ expo + contrasts .^ expo);
sumresp = (contrasts .^ expo - contrasts(j).^expo) ./ (contrasts.^expo + contrasts(j) .^expo + cfifty.^expo);
for k = 1:length(contrasts)
if sumresp(k) < 0
    sumresp(k) = -1 * (-1 * sumresp(k));
    sumresp(k) = (sumresp(k));
end
end
subplot(2,1,1);
plot(contrasts,sumresp,'g');
ratios = contrasts/contrasts(j);
%ratios = lc/lc(j);
subplot(2,1,2);
plot(ratios, sumresp, 'color',colors{j});
hold on;
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



