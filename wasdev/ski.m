price_per_visit = 106;
season_price = 855;
card_price = 149;

for j = 1:20
  if j == 1
    sum(j) = card_price + (price_per_visit * 0.6);
  elseif mod(j,6) ~= 0 | j == 1
    sum(j) = sum(j-1) + price_per_visit * 0.6;
  else
    sum(j) = sum(j-1);
  end
  plain(j) = price_per_visit * j;
end

plot(1:20,sum,'r');
hold on;
plot(1:20,plain,'g');
plot([1 20],[800 800],'b');