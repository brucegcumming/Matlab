clear Trials
for j = 1:10
  Trials(j).index = j;
end
for j = 12:100
  Trials(j).index = j;
end
for j = 1:99
  Trials(j).idx = j;
end


find([Trials.idx] == [Trials.index])
