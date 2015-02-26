function [latency, details] = LFPLatency(lfp,times,varargin)

%[latency, details] = LFPLatency(LFP,times, varargin)
%calculate latency of LFP response in subspace map from
%difference between + an - resps
% returns an nprobe x w array
% latency(probe,1) is linear regression 10-50% of initial peak,intercept
% wiht base
% 
% latency(probe,2) is the time where variance exceeds preperiod max by channel
%  latency 3 is when varianec exceeds preperiod max in all channels

smoothresp = 1;
if ndims(lfp) < 3 % just one stim
    tid  = find(times < 400);
    if size(lfp,2) == length(times)
        lfp = lfp';
    end
    if smoothresp
        [gx,gy] = meshgrid(-3:0.5:3,-3:3);
        G = abs(gx + i*gy);
        G = exp(-(G.^2));
        lfp = conv2(lfp,G,'same');
    end
    maxr = max(lfp);
    premax = max(abs(lfp(tid,:)));
    critlevel = max(premax).*1.5;
    for j = 1:size(lfp,2)
        in = find(abs(lfp(:,j)) > max(premax)*2);
        if length(in)
            id = find(abs(lfp(1:in(1),j)) < critlevel);
            aid = id(end)-1:id(end)+1;
            %            aid = in(1)-1:in(1)+1;
            latency(j,3) = interp1(abs(lfp(aid,j)),times(aid),critlevel);
        else
            latency(j,3) = NaN;
        end
    end
    details.critlevel = critlevel;
    [details.maxr, details.maxt] = max(lfp);
    [a,details.maxprobe] = max(details.maxr);
else
    %find time of max variance
    v = squeeze(sum(var(lfp,[],2),3));
    [vmax, tm] = max(v);
    pret = find(times < 400);
    vbase = mean(v(pret,:));
    basesd = std(v(pret,:));
    halfmax = (vbase+vmax)/2;
    premax = max(v(pret,:));
    [details.maxr, details.maxt] = max(v);
    [a,details.maxprobe] = max(details.maxr);

    nt = sum(max(v) > max(premax)*2);
    if nt > 16
        critlevel = max(premax) * 2;
    else
        critlevel = max(premax);
    end

    for j = 1:size(lfp,4)
        a = [];
        nloop = 1;
        %if no resps reach vmax/4, lower threhold..
        while isempty(a)
            a = find(v(:,j) > vmax(j)/(4 * nloop));
            nloop = nloop+1;
            if nloop > 1000 %error
                fprintf('LFP Var is %.2f\n',mean((v(:))));
                return;
            end
        end
        t = a(1);
        a = find(diff(v(t:end,j)) < 0);
        t = t+a(1)-1;  %find first peak
        base = mean(mean(lfp(t,:,:,j)));
        [px, py] = find(squeeze(lfp(t,:,:,j) > base));
        [nx, ny] = find(squeeze(lfp(t,:,:,j) < base));
        if size(lfp,3) > 1
            nr = mean(mean(lfp(:,nx,ny,j),2),3);
            pr = mean(mean(lfp(:,px,py,j),2),3);
        else
            nr = mean(mean(lfp(:,ny,1,j),2),3);
            pr = mean(mean(lfp(:,py,1,j),2),3);
        end
        resp = pr-nr;
        id = find(resp(1:t) > max(resp(1:t))/2);
        sid = find(resp(1:t) < max(resp(1:t))/10);
        if length(sid)
            a = polyfit([sid(end):id(1)]',resp(sid(end):id(1)),1);
            latency(j,1) =  interp1(times,- a(2)/a(1));
        else
            latency(j,1) = NaN;
        end

        resps(:,j) = resp;
        %
        % latecny 2 is the time where variance exceeds preperiod max by channel
        %
        in = find(v(:,j) > premax(j) .* 4);
        if length(in)
            id = find(v(1:in(1),j) < premax(j));
            latency(j,2) = times(id(end));
        else
            latency(j,2) = NaN;
        end
        %
        %  latency 3 is when varianec exceeds preperiod max in all channels
        in = find(v(:,j) > critlevel);
        if length(in) > 1
            id = find(v(1:in(1),j) < max(premax));
            aid = id(end)-1:id(end)+1;
            aid = in(1)-1:in(1)+1;
            latency(j,3) = interp1(v(aid,j),times(aid),critlevel);
        else
            latency(j,3) = NaN;
        end

%       latecny 4 is linear regression for just after vbase to halfmax
        ts(j) = t;
        tid = find(v(1:tm(j),j) < halfmax(j));
        if isempty(tid)
            latency(j,4) = NaN;
        else
        halft = tid(end);
        tid = find(v(1:halft,j) < vbase(j)); %
        if isempty(tid)
            latency(j,4) = NaN;
        else
        hid = find(v(tid(end):halft,j) < vbase(j) + 2 * basesd(j));
        ts = tid(end)+hid(end)-1;
        if ts < halft
        lr = polyfit(times(ts:halft),v(ts:halft,j,1)'-vbase(j),1);
        latency(j,4) = -lr(2)/lr(1);
        else
            latency(j,4) = NaN;
        end
        end
        end
        if isnan(latency(j,1))
        id = find(times < 400);
        else
        id = find(times < latency(j,1));
        end
        presd(j) = std(v(id,j));
        premean(j) = mean(v(id,j));
    end
    plot(times,resps);
    details.times = times;
    details.resps = resps;
    % use std of samples prior to latency vs maxr to estimate reliability.
    details.presd = presd;
    details.premean = premean;
    details.maxsdr = (details.maxr-details.premean)./presd;
    details.var = v;
end