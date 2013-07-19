function result = CalcFixDisp(Expt, varargin)


result.fd = NaN;
j = 1;
while j < nargin
    j = j+1;
end

% Left monocular: me = -1

ev = [Expt.Trials.Eyevals];
lid = find([Expt.Trials.me] == -1);
rid = find([Expt.Trials.me] == 1);
bid = find([Expt.Trials.me] == 0);

rmo = mean([ev(rid).ro]);
lmo = mean([ev(lid).lo]);
result.fdo = mean([ev(bid).ro] - [ev(bid).lo]) - (rmo - lmo);

rmp = mean([ev(rid).rp]);
lmp = mean([ev(lid).lp]);
result.fdp = mean([ev(bid).rp] - [ev(bid).lp]) - (rmp - lmp);

rh = [ev.ch] + [ev.vh]/2;
lh = [ev.ch] - [ev.vh]/2;
rv = [ev.cv] + [ev.vv]/2;
lv = [ev.cv] - [ev.vv]/2;
rm = mean(rh(rid));
lm = mean(lh(lid));
result.fd = mean([rh(bid)] - [lh(bid)]) - (rm - lm);

rm = mean(rv(rid));
lm = mean(lv(lid));
result.fdv = mean([rv(bid)] - [lv(bid)]) - (rm - lm);
