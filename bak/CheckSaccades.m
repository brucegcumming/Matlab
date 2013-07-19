function Expt = CheckSaccades(Expt, varargin)

th(1) = 5;
th(2) = 2;
w = 10;  %number of presaccade samples to average for position
showplot = 0;
useeyes = 0;
fixedlen = 1;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'plot',4)
        showplot = 1;
    elseif strncmpi(varargin{j},'lmonoc',4)
        useeyes = 1;
        th(1) = 7;
    elseif strncmpi(varargin{j},'rmonoc',4)
        useeyes = 2;
        th(1) = 7;
    elseif strncmpi(varargin{j},'allsz',3)
        fixedlen = 0;
    elseif strncmpi(varargin{j},'threshold',4)
        j = j+1;
        th = varargin{j};
    elseif strncmpi(varargin{j},'trials',4)
        j = j+1;
        tid = varargin{j};
        Expt.Trials = Expt.Trials(tid);
    end
    j = j+1;
end


if ~isfield(Expt.Header,'emlen') && ~isfield(Expt.Trials,'EyeData')
    Expt = ConvertEM(Expt);
end


if fixedlen
len = Expt.Header.emlen;
if isfield(Expt.Header,'emtimes')
    times = Expt.Header.emtimes;
elseif Expt.Header.CRsamplerate < 0.002 %old rufus files need to fix
times = [1:Expt.Header.emlen] .* (Expt.Header.CRsamplerate *10000); 
else
times = [1:Expt.Header.emlen] .* 1/(Expt.Header.CRsamplerate); %ticks, not ms
end
erate = 10000./mean(diff(times));
end


if Expt.Header.CRrates(1) > 0
srate = 10 .* Expt.Header.CRrates(1);
else
    srate = 10000./597;
end

if ~isfield(Expt.Header,'CRsamplerate')
    cExpt.Header.CRsamplerate = 1/(10000 * Expt.Header.CRrates(1));
    Expt.Header.CRsamplerate = 1./(1000 * srate);
end
smoothw = 2;
for j = 1:length(Expt.Trials)
    if ~fixedlen
        len = min([size(Expt.Trials(j).Eyevals.rh,1) size(Expt.Trials(j).Eyevals.lh,1) size(Expt.Trials(j).Eyevals.rv,1) size(Expt.Trials(j).Eyevals.lv,1)]);
        times = [1:len] .* 1/Expt.Header.CRsamplerate;
        Expt.Trials(j).EyeData(:,1) = Expt.Trials(j).Eyevals.lh(1:len);
        Expt.Trials(j).EyeData(:,2) = Expt.Trials(j).Eyevals.rh(1:len);
        Expt.Trials(j).EyeData(:,3) = Expt.Trials(j).Eyevals.lv(1:len);
        Expt.Trials(j).EyeData(:,4) = Expt.Trials(j).Eyevals.rv(1:len);
        erate = 10000./mean(diff(times));
        if srate < 0.017 %KLUDGE to deal with raw em files from spike2. Why different with loaded Expts?
            erate = 10./srate;
        end
    end
    Saccades = [];
    if isempty(Expt.Trials(j).EyeData)
        ch = zeros(1,len);
        cv = zeros(1,len);
    elseif useeyes == 0
        ch = (Expt.Trials(j).EyeData(:,1) + Expt.Trials(j).EyeData(:,2))/2;
        cv = (Expt.Trials(j).EyeData(:,3) + Expt.Trials(j).EyeData(:,4))/2;
    elseif useeyes == 1
        ch = Expt.Trials(j).EyeData(:,1);
        cv = Expt.Trials(j).EyeData(:,3);
    elseif useeyes == 2
        ch = Expt.Trials(j).EyeData(:,2);
        cv = Expt.Trials(j).EyeData(:,4);
    end
        

    speed = sqrt(diff(cv).^2 + diff(ch).^2) .* erate;
    speed = smooth(speed,3);
    if median(speed) > th(1)/2
        cv = smooth(cv,smoothw,'gauss');
        ch = smooth(ch,smoothw,'gauss');
        speed = sqrt(diff(cv).^2 + diff(ch).^2) .* erate;
        speed = smooth(speed,3);
    end
    
    if fixedlen
    speeds(j,:) = speed;
    end
    id = find(speed > th(1));
    if length(id)
    sid = find(speed > th(2)); 
    gaps = find(diff(id) > 10);
    ends = [id(gaps)' id(end)];
    starts = [id(1) id(gaps+1)'];
    for k = 1:length(starts)
        id = find(sid <= starts(k));
        starts(k) = sid(id(end));
        id = find(sid >= ends(k));
        ends(k) = sid(id(1));
        Saccades(k).start = times(starts(k));
        Saccades(k).end = times(ends(k));
        [Saccades(k).peakv t] = max(speed([starts(k):ends(k)]));
        peakt = starts(k)+t-1;
        Saccades(k).peakt = times(peakt);
        ipt = max([starts(k)-w 1]);
        Saccades(k).pos(1,1) = mean(ch([ipt:starts(k)]));
        Saccades(k).pos(2,1) = mean(cv([ipt:starts(k)]));
        ept = min([len ends(k)+w]);
        Saccades(k).pos(1,2) = mean(ch([ends(k)+1:min([len ends(k)+w])]));
        Saccades(k).pos(2,2) = mean(cv([ends(k)+1:min([len ends(k)+w])]));
        sc = diff(Saccades(k).pos');
        Saccades(k).size = abs(sc(1) + i * sc(2));
        Saccades(k).dir = angle(sc(1) + i * sc(2));
        if showplot
            plot(times, ch, 'r');
            hold on;
            plot(times,cv,'m');
                plot([times(ipt) times(starts(k)) times(ends(k)) times(ept)],...
                    [Saccades(k).pos(1,[1 1]) Saccades(k).pos(1,[2 2])]);
                plot([times(ipt) times(starts(k)) times(ends(k)) times(ept)],...
                    [Saccades(k).pos(2,[1 1]) Saccades(k).pos(2,[2 2])]);
                plot(times(peakt), ch(peakt),'o');
                plot(times(peakt), cv(peakt),'ro');
            if showplot == 2
                plot(times(2:end), speed);
            end
        end
    end
    end
    Expt.Trials(j).Saccades = Saccades;
end

if ~isfield(Expt.Header,'emtimes')
    
    for j = 1:length(Expt.Trials)
    if isfield(Expt.Trials,'Eyevals')
        emlen(j) = length(Expt.Trials(j).Eyevals.lh);
    elseif isfield(Expt.Trials,'EyeData')
        emlen(j) = size(Expt.Trials(j).EyeData,1);
    else
        emlen(j) = length(Expt.Trials(j).lh);
    end
    end
    emlen = round(mean(emlen));

    
    if isfield(Expt.Trials,'ftime')
    pre = [Expt.Trials.Start] - [Expt.Trials.ftime];
    else
        pre = ones(size(Expt.Trials));
    end

    Expt.Header.emtimes = ([1:emlen] .* 1./Expt.Header.CRsamplerate) - mean(pre(find(~isnan(pre))));
    Expt.Header.emtimes = times;
end

function Expt = ConvertEM(Expt)

fields = {'rh' 'lh' 'rv' 'lv'};
for t = 1:length(Expt.Trials)
    for j = 1:length(fields)
        lens(j) = length(Expt.Trials(t).Eyevals.(fields{j}));
    end
    E = Expt.Trials(t).Eyevals;
    n = min(lens);
    Expt.Trials(t).EyeData = cat(2,E.rh(1:n), E.lh(1:n), E.rv(1:n), E.lv(1:n));
    alllen(t) = n;
end
Expt.Header.emlen = min(alllen);
Expt.Header.emtimes =  [1:Expt.Header.emlen] .* (Expt.Header.CRrates(1) *10000);
