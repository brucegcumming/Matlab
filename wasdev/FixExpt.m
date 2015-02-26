function Expt = FixExpt(Expt,type)
%Expt = FixExpt(Expt,type)  Fixes misleading fields in Expt Structs
%type    'Trials'   makes sure that Expt.Trials.Trial is sequential, and
%                   that Header.BlockStart Matches

needpattern = 0;



if nargin == 1 || strncmp(type,'auto',4)
if isfield(Expt.Header,'rc') && Expt.Header.rc
    return;
end
    if isfield(Expt.Trials,'st') && isfield(Expt.Trials,'jv') && isfield(Expt.Trials,'nsf')
        if sum([Expt.Trials.st] == 18) && sum([Expt.Trials.st] == 15) %RLS and nsines
            Expt = FixExpt(Expt,'jv');
        end
    end
    f = setdiff(fields(Expt.Trials),{'Start' 'End' 'TrialsStart' 'dur' 'op' 'Trial' 'TrueEnd' 'delay' 'LFP' 'lfptime' 'FTlfp' 'st'});
    rdsf = intersect(f,{'c2' 'f2' 't2' 'a2' 'sf'}); %fields that do not affect rds 
    f2f = intersect(f,{'jx'}); %fields that do not affect grating/2gratin
    grf = intersect(f,{'nc'}); %fields that do not affect rds 
    if isfield(Expt.Trials,'st')
        bid = find([Expt.Trials.st] ==0);
        for j = 1:length(bid)
            for k = 1:length(f)
                Expt.Trials(bid(j)).(f{k}) = 0;
            end
        end
        bid = find([Expt.Trials.st] == 2);
        for j = 1:length(bid)
            for k = 1:length(rdsf)
                Expt.Trials(bid(j)).(rdsf{k}) = 0;
            end
        end
        bid = find([Expt.Trials.st] == 3);
        for j = 1:length(bid)
            for k = 1:length(grf)
                Expt.Trials(bid(j)).(grf{k}) = 0;
            end
        end
        bid = find([Expt.Trials.st] == 10);
        for j = 1:length(bid)
            for k = 1:length(f2f)
                Expt.Trials(bid(j)).(f2f{k}) = 0;
            end
            if isfield(Expt.Trials,'nph')
            if Expt.Trials(bid(j)).nph > 0 && ~isempty(strfind(Expt.Trials(bid(j)).OptionCode,'+rr'))
                Expt.Trials(bid(j)).t2 = 0;
            end
            end
        end
    end
    return;
    
end

if strcmp(type,'ed') & isfield(Expt.Trials,'ed')
    
    id = find(Expt.Trials.ed ~= 0);
    if length(id) > 1 & id(1) > 10 & diff(Expt.Trials.ed(id(1)-1:id(1))) > 0.1
        Expt.Trials.ed(1:id(1)) = Expt.Trials.ed(id(1));
    end
end

if strcmp(type,'dp')
    if isfield(Expt.Trials,'a2')
        needpattern = 1;
    end
    if needpattern
        for j = 1:length(Expt.Trials)
            or = Expt.Trials(j).or;
            a2 = or-Expt.Trials(j).a2;
            da = (or-a2) .* pi/180;
            angle = (or+a2) .*pi/180;
            dmag = Expt.Trials(j).dp .* cos(da);
            Expt.Trials(j).pdx = dmag .* cos(angle);
            Expt.Trials(j).pdy = dmag .* cos(angle);
        end
    end
end


if strcmp(type,'jv') && isfield(Expt.Trials,'jv') && isfield(Expt.Trials,'st')
    for j = 1:length(Expt.Trials)
        if isfield(Expt.Trials,'nsf')
            Expt.Trials(j).nsfs = length(Expt.Trials(j).nsf);
        end                        
        if Expt.Trials(j).st == 3
            Expt.Trials(j).jv = 0;
            if isfield(Expt.Trials,'nsf')
                Expt.Trials(j).nsf = Expt.Trials(j).sf;
                Expt.Trials(j).nsfs = 1;
            end
        elseif Expt.Trials(j).st == 18 && isfield(Expt.Trials,'nsf')
            Expt.Trials(j).nsfs = length(Expt.Trials(j).nsf);
            Expt.Trials(j).sf = Expt.Trials(j).nsf(1);
        elseif Expt.Trials(j).st == 15 && isfield(Expt.Trials,'nsf')
            Expt.Trials(j).nsfs = 0;
            if isfield(Expt.Trials,'sf')
                Expt.Trials(j).sf = 0;
            end
            if isfield(Expt.Trials,'tf')
                Expt.Trials(j).tf = 0;
            end
        elseif isfield(Expt.Trials,'nsf')
            Expt.Trials(j).nsfs = 0;
        end
    end
end

if strcmp(type,'jv') && isfield(Expt.Trials,'jv') && Expt.Stimvals.st ==18
    for j = 1:length(Expt.Trials)
        if isfield(Expt.Trials,'nsf')
            Expt.Trials(j).nsfs = length(Expt.Trials(j).nsf);
            Expt.Trials(j).sf = Expt.Trials(j).nsf(1);
        end                        
    end
end


if strncmp(type,'Trial',5)
    ts = [Expt.Trials.Trial];
    id = find(diff(ts) < 0);
    if ~isempty(id)
        toff = zeros(size(ts));
        toff(id+1) = abs(ts(id));
        toff = cumsum(toff);
        for j = 1:length(Expt.Trials)
            Expt.Trials(j).Trial = Expt.Trials(j).Trial + toff(j);
        end
        if isfield(Expt.Header,'BlockStart')
            if length(id) == length(Expt.Header.BlockStart) -1;
                Expt.Header.BlockStart(2:end) = Expt.Header.BlockStart(2:end) + toff(id+1);
            end
        end
    end
end