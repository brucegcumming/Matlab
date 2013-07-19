function [ipay, left, payment] = amortize(principal, interest, term, varargin)
%[ipay, left, payment] = amortize(principal, interest, term, varargin)

monthly = 1;
payment = 0;
printdebt = 0;
payed = [];
j = 1;
while j < nargin -2
  if strncmpi(varargin{j},'Payment',5)
        j = j + 1;
        payment = varargin{j};
  elseif strncmpi(varargin{j},'Paid',4)
      j = j+1;
      payed = varargin{j};
      payed = [0 payed];
  elseif strncmpi(varargin{j},'print',4)
      printdebt = 1;
    end
    j = j+1;
end

interest = interest/100;
if payment == 0
    intst = 1 + interest/12; %monthly interest
    mterm = term * 12;
    payment = (principal * intst^(mterm))/((1-intst^(mterm))/(1-intst));
end

left(1) = principal;
isum = 0;

mi = ((1+interest)^(1/12)) -1;
mi = interest/12;
mi = ((1+interest)^(1/12)) -1;
if(monthly)
months = term * 12;
if length(payed)
    months = length(payed);
else
    months = term * 12;
    payed = ones(1,months).*payment;
end
for j = 2:months
    ipay(j) = mi * left(j-1);
    left(j) = left(j-1) + ipay(j) - payed(j);
end
else
for j = 1:term
    ipay(j) = interest * left(j);
    ipay(j+1) = ipay(j);
    payed(j) = payment * 12;
    left(j+1) = left(j) + ipay(j) - payed(j+1);
end
end


if printdebt
    last = find(left > 0);
    last = last(end);
    fprintf('Month\tInterest\tPaid\tdebt\n');
    for j = 1:last
        fprintf('%d\t%.2f\t%.2f\t%.2f\n',j-1,ipay(j),payed(j),left(j));
    end
end
    end
