function SetCellByTemplateScore(DATA, cell, probes, trials, tn)

scores = squeeze(DATA.TemplateScores(probes, tn, trials));
sid = find(sum(scores > 0) == length(probes)); %trials where all probes have scores
[ts, id] = max(scores(:,sid));
probe(sid) = probes(id);
id = find(probe == 0);
ip = interp1(sid,probe(sid),id);
probe(id) = round(ip);
GetFigure('TemplateScore');
hold off;
imagesc(trials,[1:size(DATA.TemplateScores,1)],squeeze(DATA.TemplateScores(:,tn,trials)));
hold on;
plot(trials,probe,'w-');
GetFigure(DATA.tag.celllist);
plot(trials,probe,'w-');
a = questdlg('Apply?','popup','Cancel','OK','OK');
if strcmp(a,'OK')
DATA.CellList(cell,trials) = probe;
% not clear how to set quatlity yet.  So leave as defined for this cell
%        DATA.CellQuality(cell,trials) = probe;
set(DATA.toplevel,'UserData',DATA);
end





