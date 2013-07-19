function [details] = GetFrameReversal(mlfp, lfptimes, probes)
%
%GetFrameReversal(mlfp, lfptimes, probes)
% finds laminae where frame-locked response reverses.
frameperiod = 166.6;
if ndims(mlfp) == 3
    mlfp = squeeze(sum(mlfp,2));
end

for j = probes
    [a,b] = famp(lfptimes, mlfp(:,j),1/frameperiod);
    phases(j) = angle(b);
    amps(j) = abs(b);
end

pdiff = diff(phases);
[a,b] = max(abs(pdiff));
if a >  pi
    if pdiff(b) < 0
        id = find(phases >= phases(b));
    else
        id = find(phases < phases(b+1));
    end
   phases(id) = phases(id) + 2 * pi;
end
details.amps = smooth(amps,2,'gauss');
[a,b] = max(details.amps);
details.maxprobe = probes(b);
[a,b] = min(details.amps(1:b));
details.minprobe = probes(b);


[y,x] =smhist(phases);
[a,b] = max(y);
details.phase = x(b); %dominant phase
sp = details.phase +pi/2; %90 deg ahead - halfway to 180;
details.phases = smooth(phases,2,'gauss');
id = find(details.phases > details.phase);
lastid = find(diff(id) > 1); %find break;
if length(lastid)
id = id(1:lastid(1));
end
[a,b] = max(abs(diff(details.phases(id))));
revp = probes(b);

if sp > max(details.phases)
    if length(id) <= 1
        revp = NaN;
    else
        revp = interp1(details.phases(id),probes(id),sp,'pchip','extrap');
        if revp < 0 || revp > max(probes)
            revp = NaN;
        end
    end
elseif length(id) > 1
    revp = interp1(details.phases(id),probes(id),sp);
    else
        revp = NaN;
end
details.revp = revp;