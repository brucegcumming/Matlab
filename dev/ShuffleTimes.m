function [st, trt] = ShuffleTimes(t, Expt,varargin)
%st = ShuffleTimes(t, Expt,varargin)
%given a set of times, t, return a set of times from other equivalent
%trials for

checkhist = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'hist',4)
        checkhist = 1
    end
    j = j+1;
end
if ~isfield(Expt.Header,'preperiod')
    Expt.Header.preperiod = 1000;
end
Expt = FillTrials(Expt,'st');
checkfields = {'st' Expt.Stimvals.et Expt.Stimvals.e2 Expt.Stimvals.e3};
id = ~strcmp('e0',checkfields);
checkfields = checkfields(id);
for j = 1:length(Expt.Trials)
    ts(j) = (Expt.Trials(j).Start(1)-Expt.Header.preperiod)./10000;
    if Expt.Header.rc
        if j < length(Expt.Trials)
        shuffletrial(j) = j+1;
        else
            shuffletrial(j) = 1;
        end
    else
        id = find([Expt.Trials.(checkfields{1})] == Expt.Trials(j).(checkfields{1}));
        for k = 2:length(checkfields)
            xid = find([Expt.Trials.(checkfields{k})] == Expt.Trials(j).(checkfields{k}));
            id = intersect(id,xid);
        end
        if id(1) == j
            shuffletrial(j) = id(end);
        else
            x = find(id == j);
            shuffletrial(j) = id(x-1);
        end
    end
%    shuffletrial(j) = j;
end

xid = [];
sid = [];
for j = length(Expt.Trials):-1:1
    if j < length(Expt.Trials)
        tnext = ts(j+1);
    else
        tnext = 1+Expt.Trials(end).End(end)./10000;
    end
    tid = find(t > ts(j) & t < tnext);
    trt(tid) = t(tid) - Expt.Trials(j).Start(1)./10000;
    if shuffletrial(j) < length(Expt.Trials)
        snext = ts(shuffletrial(j)+1);
    else
        snext = 1+Expt.Trials(end).End(end)./10000;
    end
    st(tid) = ts(shuffletrial(j)) + trt(tid);
    id = find(st(tid)>snext); %if shufled spikes after trial end into a trial, 
    st(tid(id)) = st(tid(id))-0.05; %kludge
    st(tid(id)) = NaN; %kludge two
    xcounts(j) = length(id);
    xid = [xid tid(id)];
    sid = [sid find(st(tid)<=snext)];
end

if checkhist
    for j = 1:length(Expt.Trials)
        tid = find(t > Expt.Trials(j).Start(1)./10000 & t < Expt.Trials(j).End(end)./10000);
        a(tid) = t(tid) - Expt.Trials(j).Start(1)./10000;
        sid = find(st > Expt.Trials(j).Start(1)./10000 & st < Expt.Trials(j).End(end)./10000);
        b(sid) = st(sid) - Expt.Trials(j).Start(1)./10000;
        counts(j) = length(tid);
        scounts(j) = length(sid);
    end
    xv = 0:0.01:(0.1+mean([Expt.Trials.dur])./10000);
    y= histc(a,xv);
    sy = histc(b,xv);
    ty = histc(trt,xv);
    hold off; plot(xv,y); hold on; plot(xv,sy,'r');
end