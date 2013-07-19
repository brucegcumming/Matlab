function probit = consolidate(probit);

ndat = 1;
allprobit = [];

for j = unique([probit.expno])
  idx = find([probit.expno] == j);
  for val = unique([probit(idx).x])
    id = find([probit.x] == val & [probit.expno] == j);
    if ~isempty(id)
      allprobit(ndat).x = val;
      allprobit(ndat).expno = j;
      allprobit(ndat).n = sum([probit(id).n]);
      allprobit(ndat).resp = sum([probit(id).resp]);
      allprobit(ndat).p = allprobit(ndat).resp/allprobit(ndat).n;
      allprobit(ndat).name = probit(id(end)).name;
      ndat = ndat+1;
    end
  end
end
probit = allprobit;
