function DATA = SaveSpikeShape(DATA, outname)
for j = 1:length(DATA.Expts)
DATA.MeanSpike.eds(j) = DATA.Expts{j}.Stimvals.ed;
DATA.MeanSpike.Trials(j,:) = [DATA.Expts{j}.Trials([1 end]).Trial];
if isfield(DATA,'Header')
DATA.MeanSpike.Header = DATA.Header;
else
DATA.MeanSpike.probesep = 150; %default
end
end
MeanSpike = DATA.MeanSpike;
if isfield(DATA,'TrialVar')
TrialVar = DATA.TrialVar;
TemplateScores = DATA.TemplateScores;
TemplateInfo = DATA.TemplateInfo;
Templates = DATA.Templates;
fprintf('Saving MeanSpike to %s\n',outname);

save(outname,'MeanSpike','TrialVar','TemplateScores','Templates','TemplateInfo');
else
fprintf('Saving MeanSpike to %s\n',outname);
save(outname,'MeanSpike');
end


