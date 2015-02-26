function Expt = SetExptRC(Expt, varargin)
%Expt = SetExptRC(Expt, varargin) Make sure RC fields are set in Exot
%Useful if a manual Expt did not get set as RC by APlaySpkFile


keepframes = 1; %don't remove mtFn
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'short')
        keepframes = 0;
    end
    j = j+1;
end
Fr = GetEval(Expt,'Fr');
if ~(Fr > 1)
    Fr = 1;
end
fixseq = 0;
if isfield(Expt.Header,'rcvars')
    rcvars = unique({Expt.Header.rcvars{:} Expt.Stimvals.et Expt.Stimvals.e2});
else
    rcvars = {Expt.Stimvals.et Expt.Stimvals.e2};
end
rcvars = intersect(rcvars,fields(Expt.Trials));
rcvars = setdiff(rcvars,{'mtFl' 'mtFn' 'mtFi' 'e0'});
if isfield(Expt.Trials,'st') && length(Expt.Trials(1).st) >= length(Expt.Trials(1).(rcvars{1}))
    %rcvars = {rcvars{:} 'st'};
end
if isempty(rcvars)
    Expt = AddError(Expt,'No Fields to used for RC');
    return;
end
et = Expt.Stimvals.et;
if strfind(Expt.Header.Options,'+FN')
    allframes = 0;
else
    allframes = 1;
end
%check that all rcvars do in fact have the right number of samples
for j = 1:length(Expt.Trials)
    T = Expt.Trials(j);
    for k = 1:length(rcvars)
        nfs(j,k) = length(T.(rcvars{k}));
    end
end
nfc = sum(nfs);
gid = find(nfc > max(nfc * 0.9));
bid = find(nfc <= max(nfc * 0.9));
rcvars = rcvars(gid);
frameperiod = Expt.Header.frameperiod;
if min(nfs(:,1)) == max(nfs(:,1))-1
    nf = min(nfs(1,:));
    Expt = AddError(Expt,'Two Frame counts %d and %d',nf,nf+1);
end
for j = 1:length(Expt.Trials)
    T = Expt.Trials(j);
    nf = length(T.(rcvars{1}));
    if isfield(T,'Fr') && ~isempty(T.Fr)
        Fr = T.Fr;
    end

    dur = T.End(end)-T.Start(1);
    if nf * Expt.Header.frameperiod .* Fr > dur + Expt.Header.frameperiod %Need to take out repeats
        for k = 1:length(rcvars)
            et = rcvars{k};
        x{1} = T.(et)(1:Fr:end);
        x{2} = T.(et)(2:Fr:end);
        x{1} = x{1}(1:length(x{2}));
        nerr = sum(x{1}-x{2}~=0);
        if nerr
            cprintf('red','Trial %d(id%d) %s does not follow Fr properly\n',j,T.id,et);
        end
        Expt.Trials(j).(et) = x{1};
        end
        nf = nf./Fr;
    end
    if length(Expt.Trials(j).Start) < nf
        Expt.Trials(j).Start = Expt.Trials(j).Start(1) + [0:nf-1]' .* Expt.Header.frameperiod .* Fr;
        Expt.Trials(j).End = Expt.Trials(j).Start + Expt.Header.frameperiod.*Fr;
    end
    if length(Expt.Trials(j).Start) > nf
        Expt.Trials(j).Start = Expt.Trials(j).Start(1:nf);
        Expt.Trials(j).End = Expt.Trials(j).End(1:nf)
        for k = 1:length(rcvars)
            Expt.Trials(j).(rcvars{k}) = Expt.Trials(j).(rcvars{k})(1:nf);
        end
    end
    T = Expt.Trials(j);
    if isfield(T,'mtFn') && ~isempty(T.mtFn) 
        if isempty(strfind(T.OptionCode,'+FN'))
            if isfield(T,'rptframes') && length(T.rptframes) == length(T.mtFn)
                % think this only happens if op-FN. But should proabaly check for that and
                %hande op=+FN appropriately
                id = T.rptframes;
                sumdelay = 1:length(id).*Expt.Header.frameperiod;
                for k = 1:length(id)
                    if fixseq %danderous to do this with manual expts. Just let the use know about the delay
                        Expt.Trials(j).Start(id(k):end) = T.Start(id(k):end) + sumdelay(k) .* Expt.Header.frameperiod;
                        for f = 1:length(rcvars)
                            Expt.Trials(j).(rcvars{f})([id(k):end]+sumdelay(k)) = T.(rcvars{f})(id(k):end);
                        end
                    end
                end
                if length(id) > 1
                    id = T.rptframes;
                end
            else
                framets = round(T.mtFn(:)' - [1:length(T.mtFn)]);
                ds = diff(framets(1:end-1));
                id = find(ds > 0);
                sumdelay = cumsum(ds(ds > 0));
                Expt.Trials(j).framerpt = id;
                id = round(id./Fr);
                Expt.Trials(j).rptframes = id;
                %Add delays to T, which keeps record of originial times. So addition is not
                %recursive when there are two jumps
                for k = 1:length(id)
                    Expt.Trials(j).Start(id(k):end) = T.Start(id(k):end) + sumdelay(k) .* Expt.Header.frameperiod;
                end
            end
            Expt.Trials(j).ndrop = length(id);
        else %op=+FN  rptframes should already be detected
            if isfield(T,'rptframes') && ~isempty(T.rptframes)
                for k = T.rptframes(:)'
                    fx = ceil(k/Fr);
                    T.Start(fx+1:end) = T.Start(fx+1:end)+frameperiod;
                    T.End(fx+1:end) = T.End(fx+1:end)+frameperiod;
                end
                Expt.Trials(j) = T;
            end
        end
    end
end

if keepframes
    Expt.Trials = rmfields(Expt.Trials,{'mtFi' 'mtFl'}); %don't need these any more
else
    Expt.Trials = rmfields(Expt.Trials,{'mtFn' 'mtFi' 'mtFl'}); %don't need these any more
end