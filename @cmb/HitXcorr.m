function HitXcorr(src, b, id)
DATA = GetDataFromFig(src);
F = GetFigure('Counts','front');
set(F,'UserData',DATA.toplevel);
hold off;
X = DATA.xcorrs(id);
Expt = DATA.AllExpts{X.cells(1)};
if isfield(Expt.Trials,'RespDir')
iid = find(ismember([Expt.Trials.id],X.trialids));
choices = [Expt.Trials(iid).RespDir];
else
choices = ones(size(X.counts,2));
end
for j = 1:size(X.counts,2)
h = plot(DATA.xcorrs(id).counts(1,j),DATA.xcorrs(id).counts(2,j),'o','buttondownfcn',{@cmb.HitXcorrPt, id, j});
if choices(j) < 0
set(h,'color','r');
end
hold on;
end

cid = find(abs(choices) == 1);
C = classify(X.counts(:,cid)',X.counts(:,cid)',choices(cid),'quadratic');
score = sum(C == choices(cid)')./length(cid)
cp(1) = cmb.ExptCP(DATA.AllExpts{X.cells(1)});
cp(2) = cmb.ExptCP(DATA.AllExpts{X.cells(2)});
xlabel(sprintf('Cell %d (%.1f) CP%.2f',DATA.AllExpts{X.cells(1)}.Header.cellnumber,X.probes(1),cp(1)));
ylabel(sprintf('Cell %d (%.1f) CP %.2f',DATA.AllExpts{X.cells(2)}.Header.cellnumber,X.probes(2),cp(2)));
if isfield(DATA.xcorrs(id).state,'choicexc');
GetFigure('CCF','front');
hold off;
plot(DATA.xcorrs(id).state.choicexc');
hold on;
plot(diff(DATA.xcorrs(id).state.choicexc),'r');
end

