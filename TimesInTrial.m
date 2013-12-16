function [uid, postid] = TimesInTrial(vt, Expt, preperiod, postperiod, varargin)
%TimesInTrial(t, Expt, preperiod, postperiod)
%retuns list of times t that are within Trials (-preperiod + postperiod) in Expt
% t is in sec, not ticks
iid = [];
piid = [];
for j = 1:length(Expt.Trials)
     oid = find(vt > Expt.Trials(j).Start(1)./10000 - preperiod & ...
                vt < Expt.Trials(j).End(end)/10000 + postperiod);
            pid = find(vt(oid) > Expt.Trials(j).End(end)/10000); %in postperiod
            iid = [iid oid];
            piid = [piid oid(pid)];
        end
ntrials = j;
uid = unique(iid);
ignoreid = setdiff(1:length(vt),uid);
[a, postid] = ismember(piid,uid);
