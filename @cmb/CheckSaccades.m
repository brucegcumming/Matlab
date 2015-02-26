function Expts = CheckSaccades(Expts, emfile)

load(emfile);
saccth = 10;

for j = 1:length(Expts)
crate = Expt.Header.CRrates(1) * 10000;
for k = 1:length(Expts{j}.Trials)
if strfind(Expts{j}.Trials(k).OptionCode,'+2a')
bs = Expts{j}.Trials(k).TrialStart;
es = Expts{j}.Trials(k).End(end);
[tdiff, emt] = min(abs([Expt.Trials.Start] - bs));
if tdiff < 160
maxi = min([length(Expt.Trials(emt).lh) length(Expt.Trials(emt).rh) length(Expt.Trials(emt).lv) length(Expt.Trials(emt).rv)]);
ch = (Expt.Trials(emt).lh(1:maxi)+Expt.Trials(emt).rh(1:maxi))/2;
cv = (Expt.Trials(emt).lv(1:maxi)+Expt.Trials(emt).rv(1:maxi))/2;
dv = smooth(diff(cv),5);
dh = smooth(diff(ch),5);
v = 10000 * sqrt(dv.^2+dh.^2)./crate;
first = (bs -Expt.Trials(emt).ftime)/crate;
last =  (es -Expt.Trials(emt).ftime)/crate;
times = ftime + [1:length(cv)] .*crate;
saccs = find(v > 10);
id = find(saccs > last);
if ~isempty(id)
[rate, pk] = max(v(saccs(id)));
pk = id(pk);
sv(1) = mean(cv(first:saccs(pk)-10));
sv(2) = mean(cv(saccs(pk)+10:end));
sh(1) = mean(ch(first:saccs(pk)-10));
sh(2) = mean(ch(saccs(pk)+10:end));
saccsize = abs(diff(sh) + i* diff(sv));
sacdir = atan2(diff(sv),diff(sh));
end
end
end
end
end


