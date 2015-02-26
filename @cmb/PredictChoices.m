function PredictChoices(DATA);
xc = DATA.xcorrs;
Expts = DATA.AllExpts;
ids = [];
for j = 1:length(Expts)
zid = find([Expts{j}.Trials.ob] > 120 & abs([Expts{j}.Trials.RespDir]) == 1);
zids{j} = [Expts{j}.Trials(zid).id];
ids = cat(2,ids,zids{j});
end
[a,b] = Counts(ids);
bar(b,a);
c = max(a);
smw = 10;
for j = 1:floor(c/2);
smw = 10;
t = smooth(a >= c,smw);
while max(t) > 0.999
smw = smw+1;
t = conv(a >= (c-j),ones(1,smw)./smw);
end
maxlen(j) = smw-1;
end
id = find(maxlen > max(maxlen)/2);
nc = id(1);
nt = maxlen(id(1));
t = conv(a >= (c-j),ones(1,smw));
w = ceil(smw/2);
[a,c] = max(t(w:end));
w = floor(nt/2)-1;
blist = b([c-w:c+w]);
allids = [];
for j = 1:length(zids)
if sum(ismember(zids{j},blist)) > nt*0.8 
includecell(j) = 1;
if isempty(allids)
allids = zids{j};
else
allids = intersect(zids{j},allids);
end
end
end
length(allids);    
includecell = find(includecell);            
for j = 1:length(includecell)
c = includecell(j);
id = find(ismember([Expts{c}.Trials.id],allids));
counts(:,j) = [Expts{c}.Trials(id).count];
choices(:,j) = [Expts{c}.Trials(id).RespDir];
end
C = classify(counts,counts,choices(:,1),'quadratic');
scores = sum(C == choices(:,j))./size(choices,1)





