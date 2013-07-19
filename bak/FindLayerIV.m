function [lams, latency, blatency] = FindLayerIV(rc)

subbase = 0;
showplot = 0;
vars = squeeze(sum(var(rc.lfp,[],2),3));
preid = find(rc.lfptimes < 400);
b = std(vars(preid,:));
m = mean(vars(preid,:));
premax = max(vars(preid,:));

if subbase
for j = 1:size(vars,2)
    vars(:,j) = vars(:,j) - m(j);
end
end
ivars = interp1(vars,1:.1:size(vars,1));

[m, mid] = max(ivars);
for j = 1:size(ivars,2)
    for k = 1:10
    id = find(ivars(1:mid(j),j) > premax(j)*k);
    if isempty(id)
        latency(j,k) = NaN; % bad probes
    else
        latency(j,k) = id(1);
    end
    end
end
if showplot
plot(latency);
end
[a, lams(1)] = min(latency(:,2));

if isfield(rc,'lfpblank')
iblank = interp1(rc.lfpblank,1:.1:size(rc.lfpblank,1));
premax = max(abs(rc.lfpblank(preid,:)));
for j = 1:size(vars,2)
    id = find(abs(iblank(:,j)) > premax(j)*2);
    if isempty(id)
        blatency(j) = NaN; % bad probes
    else
        blatency(j) = id(1);
    end
end
if showplot
hold on;
plot(blatency,'k');
end
[a, lams(2)] = min(blatency);
end