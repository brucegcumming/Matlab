function ManData(E, LFP)

th(1) = 20;
th(2) = 5;
el = length(E.ManRec(1).Eyevals.rh);
ll = size(LFP.ManRec(1).LFP,1);
    Saccades = [];
    times = [1:el] .* (E.Header.CRrates(1)); %ticks, not ms
    ch = (E.ManRec(1).Eyevals.rh + E.ManRec(1).Eyevals.lh)/2;
    cv = (E.ManRec(1).Eyevals.rv + E.ManRec(1).Eyevals.lv)/2;
    speed = sqrt(diff(cv).^2 + diff(ch).^2) .* 400;
    speed = smooth(speed,3);
    id = find(speed > th(1));
    if length(id)
    sid = find(speed > th(2)); 
    gaps = find(diff(id) > 10);
    ends = [id(gaps)' id(end)];
    starts = [id(1) id(gaps+1)'];
    ltrig = round(starts .* ll/el);
    lid = find(ltrig > 201 & ltrig < ll-201);
    for j = 1:length(lid)
        pts = ltrig(lid(j))-200:ltrig(lid(j))+200;
        lsum(j,:,:) = LFP.ManRec(1).LFP(pts,:);
    end
    lavg = squeeze(mean(lsum));
    imagesc(lavg);
    end