function PlotQuickEMTrial(Trials, useid)

washold = ishold;

for j = 1:length(useid)
    tid = useid(j);
    T = Trials(tid);
    ns = T.neyesamples;
    h = [];
    v = [];
    h(1:ns) = NaN;
    h(1) = T.EyePos(1);
    v(1) = T.EyePos(2);
    h(ns) = T.EyePos(3);
    v(ns) = T.EyePos(4);
    for k = 1:size(T.Saccade,1)
        t = T.Saccade(k,1)+1;
        h(t) = T.Saccade(k,2);
        v(t) = T.Saccade(k,3);
        h(t-8) = T.Saccade(k,4);
        v(t-8) = T.Saccade(k,5);
    end
    t = find(~isnan(h));
    plot(t,h(t));
    hold on;
    plot(t,v(t),'r');
    hold on;
end

if washold
    hold on;
else
    hold off;
end