function Expts = CountTxtSpikes(Expts, probe, cl)

if length(cl) > 2
cl = cl(1:2);
end
for j = 1:length(Expts)
if isfield(Expts{j},'Trials')
for k = 1:length(Expts{j}.Trials)

if size(Expts{j}.Trials(k).AllSpikes,2) < cl(1)+1
Expts{j}.Trials(k).Spikes = [];
elseif length(cl) > 1 && size(Expts{j}.Trials(k).AllSpikes,2) > cl(end)
Expts{j}.Trials(k).Spikes = union(Expts{j}.Trials(k).AllSpikes{probe,cl+1});
else
Expts{j}.Trials(k).Spikes = Expts{j}.Trials(k).AllSpikes{probe,cl(1)+1};
end
end
else
fprintf('Expt %d no trials\n');
end
Expts{j}.gui.counted = 1;
end

