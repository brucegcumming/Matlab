function res = CheckLFP(LFP, varargin)

plottype = 0;
fixtype = 0;
verbose = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'fix',3)
        fixtype = 1;
        verbose = 1;
    elseif strncmpi(varargin{j},'plot',4)
        plottype = 1;
    elseif strncmpi(varargin{j},'verbose',4)
        verbose = 1;
    end
    j = j+1;
end

for j = 1:length(LFP.Trials)
    if ~isempty(LFP.Trials(j).Start)
    if isfield(LFP.Trials,'lfptime')
    delays(j) = LFP.Trials(j).Start(1) - LFP.Trials(j).lfptime;
    elseif isfield(LFP.Trials,'lfptime')
    delays(j) = LFP.Trials(j).Start(1) - LFP.Trials(j).ftime;
    else
        delays(j) = NaN;
    end
    lens(j) = size(LFP.Trials(j).LFP,1);
    chs(j) = size(LFP.Trials(j).LFP,2);
    end
end
if plottype == 1
subplot(2,1,1);
hist(delays);
subplot(2,1,2);
hold off;
hist(lens);
title(sprintf('len %.1f',mean(lens)));
end
    nch = floor(prctile(chs,90));
if(verbose)
    bid = find(chs < nch);
    cid = find(lens < max(lens));
    if length(bid) | length(cid)
        fprintf('%d Short Trials. %d missing chans\n',length(cid),length(bid))
    end
end

if fixtype == 1
    minlen = prctile(lens,10);

    id = find(chs >= nch & lens >= minlen);
    if length(id) < length(LFP.Trials)
        fprintf('Omintting %d short trials\n',length(LFP.Trials)-length(id));
        LFP.Trials = LFP.Trials(id);
    end
    if minlen < max(lens)
        if minlen < max(lens)
            fprintf('Setting Length to %d (%d)\n',minlen,max(lens));
        end
        for j = 1:length(LFP.Trials)
            LFP.Trials(j).LFP = LFP.Trials(j).LFP(1:minlen,1:nch);
            LFP.Trials(j).FTlfp = fft(LFP.Trials(j).LFP);
        end
        LFP.Header.lfplen = minlen;
        if isfield(LFP.Header,'LFPtimes')
            LFP.Header.LFPtimes = LFP.Header.LFPtimes(1:minlen);
        end
    end
    if nch < max(chs)  %some trials with more LFP chans
        for j = 1:length(LFP.Trials)
            LFP.Trials(j).LFP = LFP.Trials(j).LFP(1:minlen,1:nch);
            LFP.Trials(j).FTlfp = fft(LFP.Trials(j).LFP);
        end
    end
    res = LFP;
end

res.delays = delays;
res.lens = lens;
res.nch = chs;


