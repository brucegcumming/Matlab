
function [LFP, avgs, pwrs] = FixLFPMains(LFP, tics, varargin)
% LFP = FixLFPMains(LFP, tics)
% remove mains artifact from LFP by calculating triggered averad
% Also scales channels 17:24 by 20 to match gain
% FixLFPMains(LFP,tics,'gain', xx) forces a different scale factor
%
%tics is the timestamps for mains phase, taken if necessary from the other
%.mat file

%before Sep19, had a gain difference between LFP chans of 20
%New boards on Sep 19 all have the same gain.

rebuild = 0;
maxtime = 0;
exptlist = [];
if ischar(tics)
    varargin = {tics varargin{:}};
end
for j = 1:length(varargin)
    if strcmp(varargin{j},'rebuild')
        rebuild = 1;
    elseif strcmp(varargin{j},'exptlist')
        j = j+1;
        exptlist = varargin{j};        
    elseif strcmp(varargin{j},'maxtime')
        j = j+1;
        maxtime = varargin{j};
    end
    j = j+1;
end

if ischar(LFP) && isdir(LFP) %find .lfp files in dir and process all
    d = mydir([LFP '/*.lfp.mat']);
    Headers = {};
    if strcmp(tics,'rebuild')
        rebuild = 1;
    end
    for j = 1:length(d);
        eid = GetExptNumber(d(j).name);
        if (isempty(exptlist) || ismember(eid,exptlist));
        fprintf('Loading %s\n',d(j).name);
        load(d(j).name);
        if ~isfield(LFP,'Trials')
            fprintf('No Trials in %s\n',d(j).name);
        elseif (~isfield(LFP.Header,'MainsNoise') ||  rebuild) 
            mfile = regexprep(d(j).name,'A(\.[0-9]*).lfp.mat','$1.mat');
            eid = GetExptNumber(d(j).name);
            if exist(mfile)  && (isempty(exptlist) || ismember(eid,exptlist));
                X = load(mfile);
                f = fields(X);
                if isfield(X,'Ch31')
                    maxtime = max(X.Ch31.times).*10000;
                end
                tics = [];
                for k = 1:length(f);
                    if isfield(X.(f{k}),'title') && strcmp(X.(f{k}).title,'Mains')
                        tics = X.(f{k}).times .*10000;
                        maxtime = max([max(tics) maxtime]);
                        LFP = FixLFPMains(LFP,tics,'maxtime',maxtime);
                        Headers{j} = LFP.Header;
                        if isfield(LFP.Header,'MainsNoise')
                            save(d(j).name,'LFP');
                        end
                    end
                end
                if isempty(tics)
                    cprintf('error','No Mains Record in %s\n',mfile);
                end
            end
        else
            Headers{j} = LFP.Header;
        end
        if isfield(LFP.Trials,'MainsGain')
            Headers{j}.MainsGain = [LFP.Trials.MainsGain];
        end
        end
    end
    LFP = Headers;
    return;
elseif ischar(LFP)
    cprintf('red','No Directory %s\n',LFP);
    return;
end

if datenum(LFP.Header.RecDate) > 733669
    forcegain = 1;
else
    forcegain = 20;
end
plotmean = 1;
trackpwr = 0;

rate = LFP.Header.CRsamplerate .* 10000; %tics/sample
len = ceil(333.3./rate);
blocks = [];
avgs = [];
pwrs = [];
fillchannel = 0;
nerr = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'blocks',3)
        j = j+1;
        blocks = varargin{j};
        if length(blocks) == 1
            blocks = [1:blocks:length(LFP.Trials)+1];
        end
    elseif strncmpi(varargin{j},'fillch',5)
        j = j+1;
        fillchannel = varargin{j};
    elseif strncmpi(varargin{j},'gain',3)
        j = j+1;
        forcegain = varargin{j};
    end
    j = j+1;
end

if isempty(blocks)
    blocks = [1 length(LFP.Trials)+1];
end

for j = 1:length(LFP.Trials)
    lens(j) = size(LFP.Trials(j).LFP,1);
    chs(j) = size(LFP.Trials(j).LFP,2);
    starts(j) = LFP.Trials(j).Start;    
end

%id Start AND id both step backwards, its the bug when spike2 failed to
%clear LFP
if isfield(LFP.Header,'Spike2ver')
    s2ver = LFP.Header.s2ver .* 100;
else
    s2ver = 0;
end
%breaks = find(diff([LFP.Trials.Start]) < 0 & diff([LFP.Trials.id]) <0);
breaks = find(diff([LFP.Trials.id]) <0);
if ~isempty(breaks)  && s2ver < 100;
    id = 1:(breaks(1));
    gid = [LFP.Trials(id).id];
    nid = setdiff([LFP.Trials.id],gid);
    if ~isempty(nid)
        nerr = nerr+1;
        err = sprintf('%s has ids %d - %d, but then "Hanging" Trials after  %d - %d',LFP.Header.Name,min(gid),max(gid),min(nid),max(nid));
        cprintf('red','%s\n',err);
        LFP.errs{nerr} = err;
    end
    LFP.Trials = LFP.Trials(id);
    chs = chs(id);
    lens = lens(id);    
end

nch = prctile(chs,90);
goodch = find(chs >= nch);
for j = goodch

    if fillchannel
        LFP.Trials(j).LFP(:,fillchannel) = mean(LFP.Trials(j).LFP(:,[fillchannel-1 fillchannel+1]),2);
    end
    start = LFP.Trials(j).ftime;
    last = start+ length(LFP.Trials(j).LFP) .* rate;
    if ~isempty(start)
        t = start+ [1:length(LFP.Trials(j).LFP)] .* rate;

        mid = find(tics > LFP.Trials(j).ftime & tics < last+170);
        if length(mid) < 2
            period = median(diff(tics([mid:mid+5])));
        else
            period = median(diff(tics(mid)));
        end

        kid = 1;
        id = [];
        if isempty(mid)
            id(1:size(LFP.Trials(j).LFP)) = NaN;
        else
            for k = 1:size(LFP.Trials(j).LFP,1)
                while t(k) > tics(mid(kid)) & kid < length(mid)
                    kid = kid+1;
                end
                if kid > length(mid)
                    fprintf('Trial longer than Mains pulses at %.3f\n', tics(mid(end))./10000);
                    id(k) =  NaN;
                else
                    p = len .* (t(k)-tics(mid(kid)))./period;
                    id(k) = len + ceil(p);
                    rem(k) = mod(p,1);
                end
            end
        end
        nid = find(id <=0);
        id(nid) = id(nid)+len;
        nid = find(id <=0);
        if ~isempty(nid)
            fprintf('Trial %.3f precedes Mains pulses at %.3f\n', t(k)./10000,tics(1)./10000);
            id = 1+mod(id-1,len);
        end
        nid = find(id > len);
        if ~isempty(nid)
            fprintf('Trial %d (%.3f - %.3f) %.3f long gap before %.3f\n', j,start./10000,last./10000,t(k)./10000,tics(mid(end))./10000);
            id = 1+mod(id-1,len);
        end
        ids{j} = id(find(~isnan(id)));
        if isempty(ids{j})
            nerr = nerr+1;
            err = sprintf('Trial %d No ids', j);
            cprintf('red','%s\n',err);
            LFP.errs{nerr} = err;
        end
        if j >= 717
            weighted = 0;
        end
        weighted = 0;

        if weighted
            rems{j} = rem(find((~isnan(id))));
            for k = 1:len
                kid = find(id==k);
                ns(k) = length(kid);
                if length(kid)
                    avg(k,:) = mean(LFP.Trials(j).LFP(kid,:));
                    aavg(k,:) = sum(LFP.Trials(j).LFP(kid,:) .* repmat(rems{j}(kid)',1,nch));
                    if k == 1
                        bavg(len,:) = sum(LFP.Trials(j).LFP(kid,:) .* (1-repmat(rems{j}(kid)',1,nch)));
                    else
                        bavg(k-1,:) = sum(LFP.Trials(j).LFP(kid,:) .* (1-repmat(rems{j}(kid)',1,nch)));
                    end
                else
                    aavg(k,:) = 0;
                    bavg(k,:) = 0;
                end
            end
            for k = 1:len
                avg(k,:) = (aavg(k,:) + bavg(k,:))./ ns(k);
            end
        elseif nch > 0
            for k = 1:len
                kid = find(id==k);
                if length(kid)
                    avg(k,:) = mean(LFP.Trials(j).LFP(kid,1:nch));
                else
                    avg(k,:) = zeros(1,nch);
                end
            end
        else
            avg(1:len,1) = 0;
        end
        for p = 1:size(avg,2)
            avgs(j,:,p) = avg(:,p) - mean(avg(:,p));
        end
        lens(j) = length(LFP.Trials(j).LFP);
    end
end
tlen = min(lens);

if isempty(avgs)
    return;
end
if plotmean && size(avgs,1) > 100
    k = 1;
    for j = 1:100:size(avgs,1)-100
        Z(k,:) = squeeze(mean(avgs(j:j+100,:,plotmean)));
        k = k+1;
    end
    imagesc(Z);
end

if nch > 1
    if isfield(LFP.Header,'MainsNoise') && trackpwr
        navg = LFP.Header.MainsNoise;
    else
        navg = squeeze(mean(avgs,1));
    end
    ugain = sum(navg(:).*navg(:));
    for b = 1:size(avgs,1)
        gains(b) = sum(sum(squeeze(avgs(b,:,:)) .* navg))./ugain;
        LFP.Trials(b).MainsGain = gains(b);
    end
    for b = 1:length(blocks)-1;
        id = find(goodch >= blocks(b) & goodch < blocks(b+1));
        id = goodch(id);
        if isfield(LFP.Header,'MainsNoise') && trackpwr
            navg = LFP.Header.MainsNoise;
        else
            navg = squeeze(mean(avgs(id,:,:),1));
        end
        ugain = sum(navg.*navg);
    if size(avgs,4) > 16
        aavg = mean(squeeze(mean(avgs(id,:,1:16),1)),2);
    bavg = mean(squeeze(mean(avgs(id,:,17:end),1)),2);
    dgain = std(bavg)./std(aavg);
    else
        aavg = mean(squeeze(mean(avgs(id,:,:),1)),2);
    end
    subplot(2,1,1);
    plot(squeeze(mean(avgs,1)));
    if forcegain
        dgain = forcegain;
    end
%could subtract average lumped over chans 1-16 and 17-24. But when using
%average from whole expt,(no blocking) there is no need, and this avoids
%problems with missing channels, mis-scaled channels, etc.
    if nch > 16
        ach = 16;
    else
        ach = nch;
    end
    for j = id
        if length(ids{j}) == size(LFP.Trials(j).LFP,1)
                pwr = LFP.Trials(j).LFP .* avg(ids{j},:);
                pwrs(j,:) = sum(pwr);
            if trackpwr
                for k = 1:ach
                    LFP.Trials(j).LFP(:,k) = LFP.Trials(j).LFP(:,k) - navg(ids{j},k);
                end
            else
        for k = 1:ach
            LFP.Trials(j).LFP(:,k) = LFP.Trials(j).LFP(:,k) - navg(ids{j},k);
        end
        for k = 17:nch
            LFP.Trials(j).LFP(:,k) = LFP.Trials(j).LFP(:,k) - navg(ids{j},k);
        end
            end
        end
        LFP.Trials(j).LFP(:,17:end) = LFP.Trials(j).LFP(:,17:end)./dgain;
        LFP.Trials(j).LFP = single(LFP.Trials(j).LFP); %save some disk space
        %   LFP.Trials(j).FTlfp = fft(LFP.Trials(j).LFP);
    end
    end
    subplot(2,1,2);
    plot(pwrs);
    xlabel('Trial');
    ylabel('Template Amplitude');
    LFP.Header.MainsNoise = navg;
    LFP.Header.adjgain = dgain;
else
    plot(mean(avgs));
    avg = mean(avgs(1:blocks(1),:));
    for j = 1:length(LFP.Trials)

        avg = avgs(j,:);
        LFP.Trials(j).LFP= LFP.Trials(j).LFP - avg(ids{j})';

        LFP.Trials(j).FTlfp = fft(LFP.Trials(j).LFP);
    end
end
