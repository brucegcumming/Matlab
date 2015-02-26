function [skipframes, derr, nskip] = CheckFrameDiffs(t, varargin)
%Gnve a set of frametime differences, find skips
Nf=401;
showplot = 0;
frameperiod = 10;
method = 2;
nskip = [];
    
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'frameperiod',6)
        j = j+1;
        frameperiod = varargin{j};
    elseif strncmpi(varargin{j},'plot',4)
        showplot = 1;
    end
    j = j+1;
end

skipframes = [];
derr = [];
if isempty(t)
    return;
end
cs = cumsum(t(:)) - [1:length(t)]' .* frameperiod;
crit = median(rem(cs+3*frameperiod/2,frameperiod)-frameperiod/2);
nskip = 0;
nomdur = (1+length(t)) .* frameperiod;

if method == 2
    smw = 1;
    ids = 1 + find(diff(cs) > frameperiod/2);
    while isempty(ids) && smw < 5
        smw = smw+1;
        xcs = downsample(smooth(cs,smw),smw);
        ids = 1 + find(diff(xcs) > frameperiod/2);
        ids = ids .* smw;
    end
    if isempty(ids)
        if cs(1) > 0
            skipframes(1) = 1;
        end
        derr = sum(t)-nomdur - frameperiod .* length(skipframes);
        return;
    end
    lastid = 1;
    for j = 1:length(ids)
        frametimes(j) = mean(cs(lastid:ids(j)-1));
        lastid = ids(j);
    end
    frametimes(end+1) = mean(cs(ids(end):end));
    x = find(diff(frametimes) > frameperiod/2);
    skipframes = ids(x);
    if frametimes(1) > 0
        skipframes = [1 skipframes(:)'];
    end
    if frametimes(end) > 1 && length(frametimes) > 1
        nskip = round(diff(frametimes));
        nskip = nskip(x);
    elseif  length(frametimes) > 2
        nskip = round(diff(frametimes));
        nskip = nskip(x);
    end
    derr = sum(t)-nomdur - frameperiod .* length(skipframes);    
else
skip =  -5;
sid = find(cs < skip);
lastskip = 1;
while isempty(sid)
    nskip = nskip+1;
    skip = skip+frameperiod;
    frametimes(nskip) = -frameperiod/2;
    sid = find(cs < skip);
    skipframes(nskip) = 1;
end
lid = find(cs > skip+5);
while ~isempty(sid) && ~isempty(lid)
    if lid(end) - sid(end) > 3 % there is a skip, but check it is here
        nskip = nskip+1;
        frametimes(nskip) = mean(cs(lastskip:sid(end)));
        if length(sid) > 1
            nid = find(cs(sid(end):end) - frametimes(nskip) > frameperiod/1.5);
        else
            nid = find(cs(sid(end):end) - cs(1) > frameperiod/1.5);
        end
        newtime = mean(cs(sid(end):sid(end)+nid(1)));
        skipframes(nskip) = sid(end) + nid(1)-1;
        frametimes(nskip) = mean(cs(lastskip:skipframes(nskip)));
        lastskip = skipframes(nskip);
    end
    skip =  skip+10;
    lid = find(cs > skip+5);
    sid = find(cs < skip);
end


if nskip
    lastskip = 1;
    for j = 1:nskip
        if skipframes(j) == 1
            frametimes(j) = -frameperiod/2;
        else
            frametimes(j) = mean(cs(lastskip:skipframes(j)));
        end
        lastskip = skipframes(j);
    end
    if length(cs) > skipframes(end)
        frametimes(end+1) = mean(cs(skipframes(end):end));
    end
    id = find(diff(frametimes) < frameperiod/2);
    skipframes(id) = [];
    derr = sum(t)-nomdur - frameperiod .* length(skipframes);
else
    derr = 0;
end
end    
if showplot
    plot(cs);
    hold on;
    yl = get(gca,'ylim');
    for j = 1:length(skipframes)
        line([skipframes(j) skipframes(j)], yl,'color','r');
    end
    ylabel('error in ms');
    xlabel('Frame');
end