function res = DeTrendExpt(res, varargin)

tw = 20; %20 trials, all types
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'width',3)
        j = j+1;
        tw = varargin{j};
    end
    j = j+1;
end


Trials = res(1).Data.Trials;
e1 = res(1).type{1};
if length(res(1).type) > 1
e2 = res(1).type{2};
  
end

for j = 1:length(Trials)
if length(res(1).type) > 1 && ~strcmp(res(1).type{2},'e0')
   ids = find(res(1).x == Trials(j).(e1) & res(1).y == Trials(j).(e2));
else
   ids = find(res(1).x == Trials(j).(e1));
end
   idlist(j) = ids(1);
   nid(j) = length(ids);
   counts(j) = Trials(j).count;
end

means = cat(3,res.means);
means(isnan(means)) = 0;
ns = cat(3,res.n);
means = sum(means .* ns,3)./sum(ns,3);

predrates = means(idlist);
if diff(size(predrates)) < 0
    predrates = predrates';
end
adjrates = smooth(predrates -counts,tw);
for j = 1:length(Trials)
    Trials(j).count = Trials(j).count + adjrates(j);
end

res(1).Data.Trials = Trials;