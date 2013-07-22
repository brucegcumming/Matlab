function mortgage(growthrate)

if ~exist('compare','var')
    compare = 1;
end

if ~exist('taxrate','var')
    taxrate = 0.25;
end

if ~exist('growthrate','var')
    growthrate = 0.05625;
end

% sucessfully replicates Marc Dorfmans fax, if payment is set to 
[ipay, left] = amortize(400000, 5.625, 30, 'Payment', 2302.63);
capital = 400000;
payment = 2302.63;

%Two scenarios when buying a house for 400K. 
%scenario 1 get mortgage, invest 400k with growth rate growth.
for j = 1:30
    capsum(1,j) = capital;
    growth = capital * growthrate;
    capital = capital + growth * (1-taxrate);
end

%scenario 2. Buy house with 400K, then have payment to invest each month.
capital = [0 0];
for j = 1:30
    k = j-1;
    capsum(2,j) = capital(1);
    capsum(3,j) = capital(2);
%    capital = capital * growth + payment * 12;
taxbreak(j) = sum(ipay(1+(k*12):1+(j*12)));
    growth(1) = (capital(1) + payment * 6) * growthrate + payment * 6 + taxbreak(j);
    growth(2) = (capital(2) + payment * 6) * growthrate + payment * 6;
    capital = capital + growth .* (1-taxrate);
end


if compare
    hold off;
    plot(capsum(1,:),'r');
    hold on;
    plot(capsum(2,:),'g');
    plot(capsum(3,:),'g');
    plot(taxbreak,'b');
else
    hold off;
    plot(1:length(left),left,'b');
    hold on;
    plot(1:length(ipay),cumsum(ipay),'r');
    plot((1:length(capsum)).*12,capsum,'g');
end