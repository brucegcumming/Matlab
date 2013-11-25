function Expt = FixExpt(Expt,type)

needpattern = 0;
if nargin == 1 || strncmp(type,'auto',4)
    if isfield(Expt.Trials,'st') && isfield(Expt.Trials,'jv') && isfield(Expt.Trials,'nsf')
        if sum([Expt.Trials.st] == 18) && sum([Expt.Trials.st] == 15) %RLS and nsines
            Expt = FixExpt(Expt,'jv');
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
