function scores = CalcIsolation(DATA, cnum, varargin)
%try mulitple spaces to check isolation number

id = find(DATA.clst == cnum+1);
nid = find(DATA.clst ~= cnum+1 & DATA.clst > 0);
tmpspaces = { [2 10] [3 10] [2 5] [1 8] [2 8 10] [2 11] [2 12] [2 11 12]};
pcspaces = {[1 2] [1 2 3]  [1 2 3 4]};
for j = 1:length(tmpspaces)
    x = CalcIsolation(DATA.TemplateScores(:,tmpspaces{j}),DATA.clst,cnum+1);
    scores(j) = x(1);
end

for j = 1:length(pcspaces)
    x = CalcIsolation(DATA.pcs(:,pcspaces{j}),DATA.clst,cnum+1);
    scores(end+1) = x(1);
end
