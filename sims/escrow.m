balence = 4544;
payment = 651.49
taxes = 6585 + 1098;
months = 24;
for j = 1:months;
    balence = balence + payment;
    paid(j) = j * payment;
    if mod(j-1,6) == 0
        balence = balence - taxes/2;
    end
    total(j) = balence;
end

plot(1:months,total);