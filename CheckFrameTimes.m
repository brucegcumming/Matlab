function E = CheckFrameTimes(E, varargin)
%Go through Trials from bnc file checking for dropped frames
%See also CheckStimDur, CheckFrameDiffs

usenf = [];
frameperiod = 10;
nohold = 1;
j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'framerate',6)
        j = j+1;
        frameperiod = 1000./varargin{j};
    elseif strncmp(varargin{j},'hold',4)
        nohold = 0;
    elseif strncmp(varargin{j},'nf',2)
        j = j+1;
        usenf = varargin{j};
    end
    j = j+1;
end


GetFigure('FrameTimes');
if iscell(E)
    for j = 1:length(E)
        CheckFrameTimes(E{j},'hold');
    end
    return;
end

if nohold
    hold off;
end

totaldrop = 0;
    dropframes = [];
    droptrials = [];
if isfield(E,'Trials') && isfield(E.Trials,'nframes') %REadSerialTrials
    if ~isempty(usenf)
        id = find([E.Trials.nf] == usenf);
        E.Trials = E.Trials(id);
    end
    if isfield(E.Trials,'fz')
        frameperiod = 1000./prctile([E.Trials.fz],50);
    end
    for t = 1:length(E.Trials)
    T = E.Trials(t);
    if isfield(T,'nf') && isempty(usenf)
        Nf = T.nf;
    else
        Nf = usenf;
    end
    E.Trials(t).ndrop = 0;
    if T.nframes >= Nf && T.duration > Nf * frameperiod/1000
        binocc(t) = length(T.mtFl);
        fd = diff(T.mtFi);
        plot(fd);
        cs = cumsum(T.mtFi(:)) - [1:length(T.mtFi)]' .* frameperiod;
        skipframes = CheckFrameDiffs(T.mtFi,'frameperiod',frameperiod);
        fdd = diff(fd);
        hold on;
        plot(fdd,'r');
        [a,b] =  sort(fd);
        fsum(t) = sum(T.mtFi);
        E.Trials(t).rptframes = skipframes;
        E.Trials(t).ndrop = length(E.Trials(t).rptframes);
        totaldrop = totaldrop+length(skipframes);
        dropframes = [dropframes(:)' skipframes(:)'];
        droptrials = [droptrials(:)' ones(size(skipframes(:)))' .* t];
    else
        nodata(t) = 1;
    end
    end
elseif isfield(E,'Trials') %ordinary Expt
    cs = [];
    for t = 1:length(E.Trials)
        T = E.Trials(t);
        if isfield(T,'mtFn')
            nf = length(T.mtFn);
            [skips, terr, nskip] = CheckFrameDiffs(diff(T.mtFn),'frameperiod',1);
            skipid = skips < nf-1; %'skip' on last frame is normal
            skipframes{t} = skips(skipid);
            nskips{t} = nskip(skipid);
            cs(t,1:nf) = T.mtFn(:)' - [1:nf];
            if isempty(skipframes{t}) && nf > 20
                baseline(t) = mean(cs(t,10:nf-01));
            end
        end
    end
    if ~isempty(cs);
        plot(cs');
        hold on;
    end
end

if totaldrop
    plot(droptrials,dropframes,'o');
    for j = 1:length(E.Trials)
        if isempty(E.Trials(j).ndrop)
            E.Trials(j).ndrop = 0;
        end
    end
end