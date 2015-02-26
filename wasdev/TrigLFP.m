function [avg, details] = TrigLFP(Trials, times, samplerate, nch, varargin)
%
%TrigLFP(Trials, times, samplerate, nch)
%Builds a triggered LFP for a set of trials. From times(1) to times(2)
%relative to the trigger. Trigger times are relative to Trial Start, NOT
%abolute times. Trigger Times can be set in the Trials, or given as an
% optional cell array:%
% TrigLFP(Trials, times, samplerate, nch, 'triggers', trig)
%      where Trails(n).trigger is set to trig{n};
% sampletime is the duration of a single sampel in tics 
% I.e. the units of Header.LFPsamplrate
%
latency = 500;
nspk = 0;


%samplerate is in sec  per sample. 
% so ns * samplerate = duration in sec.
lfplen =  round(diff(times)./(samplerate * 10000));
details.lfptimes = ([1:lfplen] .* samplerate * 10000) + times(1);

lfpavg = zeros(lfplen,nch);
lfpavn = zeros(lfplen,1);
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'triggers',4)
        j = j+1;
        triggers = varargin{j};
        if iscell(triggers)
        for k = 1:min([length(Trials) length(triggers)])
            Trials(k).Trigger = triggers{k};
        end
        else
            for k = 1:min([length(Trials) length(triggers)])
                Trials(k).Trigger = triggers(k,:);
            end
        end
        
        
    end
    j = j+1;
end

if(~isfield(Trials,'Trigger'))
  [Trials.Trigger] = deal(0);
end
nt = length(cat(2,Trials.Trigger));


len = length(lfpavg);
lfpall = ones(len,nch,nt) .* NaN;
nt = 0;
for tr = 1:length(Trials);
    lfp = Trials(tr).LFP;
    if ~isempty(lfp)
        for j = 1:length(Trials(tr).Trigger)
            nt = nt+1;
            first = 1+round((Trials(tr).Trigger(j) +times(1)+Trials(tr).Start(1)-Trials(tr).lfptime)./(samplerate * 10000));
%            first = first+Trials(tr).lfpo;
        if first <= 0  %looking for samples that precede the data
            start = 1;
            last = round(lfplen+first);
            istart = 1-first;
            iend = lfplen;
            if last > size(lfp,1)
                last = size(lfp,1);
                iend = last+istart-1;
            end
        elseif first+lfplen > length(lfp) %need samples after data is over
            last = length(lfp);
            start = first;
            istart = 1;
            iend = 1+last-start;
        else
            start = first;
            last = first+lfplen-1;
            istart = 1;
            iend = lfplen;
        end
        if ~isnan(last) & ~isnan(start) & iend > istart
            if 0
            lfpall(istart:iend,:,nt) = lfp(start:last,:);
            else
                idx = find(~isnan(mean(lfp(start:last,:),2)));
            lfpavg(istart:iend,:) = lfpavg(istart:iend,:) + lfp(start:last,:);
            lfpavn(istart:iend) = lfpavn(istart:iend)+1;
            end
        end
%        lfpall(nlfp,:) = lfp(start:end);
        end
    end
end
if 0
    avg = nanmean(lfpall,3);
    lfpavn = sum(~isnan(sum(lfpall,2)),3);
else
    idx = find(lfpavn == 0);
    lfpavg(idx) = NaN;
    lfpavn(idx) = 0;
    avg = lfpavg ./ repmat(lfpavn,1,nch);
end
details.n = lfpavn;
details.ntrials = nspk;