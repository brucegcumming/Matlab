in = [];
out = [];
[io, x, s] = textread('taxi.dat','%s%n%s%*[^\n]','delimiter',':');
id = strmatch('O',io);
fprintf 'Expenses\n';
totalout = 0;
for j = id'
    total = x(j);
    n = sscanf(io{j},'O%d');
    if isempty(n)
        totalout = totalout+total;
    else
        totalout = totalout+total*n;
    end
    if total
        fprintf('%.2f %.2f %s\n',totalout,total,s{j});
    end
end

fprintf 'Receipts\n';
totalin = 0;
id = strmatch('I',io);
for j = id'
    totalin = totalin+x(j);
    if x(j)
        fprintf('%.2f %.2f %s\n',totalout-totalin,x(j),s{j});
    end
end
fprintf('Total Receipts %.2f\n',totalin);
