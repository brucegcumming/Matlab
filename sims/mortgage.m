function [capsum, growth, payment] = mortgage(growthrate,varargin)

taxrate = 0.25;
payment = 0;
capital = 400000;
borrow = 400000;
term = 30;
itaxrate = taxrate;
ctaxrate = 0;
income = 2302.61;
noise = 0;

j = 1;
while j < nargin
    if strncmpi(varargin{j},'Payment',5)
        j = j + 1;
        payment = varargin{j};
    elseif strncmpi(varargin{j},'taxrate',5)
        j = j + 1;
        taxrate = varargin{j};
    elseif strncmpi(varargin{j},'capital',5)
        j = j + 1;
        capital = varargin{j};
    elseif strncmpi(varargin{j},'borrow',5)
        j = j + 1;
        borrow = varargin{j};
    elseif strncmpi(varargin{j},'noise',5)
        j = j + 1;
        noise = varargin{j};
    elseif strncmpi(varargin{j},'income',5)
        j = j + 1;
        income = varargin{j};
    elseif strncmpi(varargin{j},'term',4)
        j = j + 1;
        term = varargin{j};
    elseif strncmpi(varargin{j},'taxmodel',5)
        j = j + 1;
        if strncmpi(varargin{j},'income',5)
            itaxrate = taxrate;
            ctaxrate = 0;
        elseif strncmpi(varargin{j},'capital',5)
            itaxrate = 0;
            ctaxrate = taxrate;
        elseif strncmpi(varargin{j},'none',4)
            itaxrate = 0;
            ctaxrate = 0;
        end
    end
    j = j+1;
end

% sucessfully replicates Marc Dorfmans fax, if payment is set to 2302.61,
% borrow 400,000

if(payment > 0)
    [ipay, left] = amortize(borrow, 5.625, term, 'Payment', payment);
else
    [ipay, left, payment] = amortize(borrow, 5.625, term);
end

%p is what you would pay if you borrowed 'borrow'. If you have borrowed
%less than thie, then each month you gain p - payment, which you will of
%course invest in the same way as your capital.
msave(1:term) = income - payment;
if term < 30
msave(term:30) = income;
ipay(term*12+1:1+30*12) = 0;
end

if(noise)
    growthrates = (mod([1:30],2) - 0.5)/10 + growthrate;
else
    growthrates = zeros(1,30) + growthrate;
end

for j = 1:30
    growth(j) = ((capital + msave(j) * 6) * growthrates(j)) * (1- itaxrate);
    k = j-1;
    capsum(j) = capital;
%    capital = capital * growth + payment * 12;
taxbreak(j) = sum(ipay(1+(k*12):1+(j*12))) * taxrate;
    capital = capital + growth(j)  + msave(j) * 12 + taxbreak(j);
end


capsum(31) = capsum(30) - ctaxrate * sum(growth);
