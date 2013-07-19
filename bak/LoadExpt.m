function Expt = LoadExpt(name);

if ~exist(name,'file')
    if(name(1) == '/')
        Expt = [];
        return;
    else
        load(name2path(name));
    end
else
    load(name);
end
if ~exist('Expt','var') && exist('cExpt','var')
    Expt = cExpt;
    return;
end
if ~exist('Expt','var') && exist('muExpt','var')
    Expt = muExpt;
    return;
end


idx = [];
if exist('AllExpt','var')
    AllExpt.Expt.Header.loadname = name;
    Expt = AllExpt;
    return;
end
if ~exist('Expt','var') | ~isfield(Expt,'Trials')
    if exist('cExpt','var') | ~isfield(cExpt,'Trials')
        Expt =  cExpt;
    else
    Expt = [];
    return;
    end
end
for j = 1:length(Expt.Trials)
    if ~isempty(Expt.Trials(j).Start)
        idx = [idx j];
    end
end

Expt.Header.loadname = name;

if isnumeric(Expt.Stimvals.et) 
    if strfind(name,'.AC.')
        Expt.Stimvals.et = 'dx';
        Expt.Stimvals.e2 = 'ce';
    elseif strfind(name,'OXAC.')
        Expt.Stimvals.et = 'dO';
        if ~isfield(Expt.Trials,'dO')
            Expt = FillTrials(Expt,'dO');
        end
        Expt.Stimvals.e2 = 'ce';
    else
    etvals = [6 117  9 132 119];
    fprintf('%s Expts %dX%d\n',name,Expt.Stimvals.et,Expt.Stimvals.e2);
    etnames = {'dx' 'ce' 'dx' 'ce' 'ce'};
    id = find(Expt.Stimvals.et == etvals);
    if length(id)
        Expt.Stimvals.et = etnames{id(1)};
    end
    id = find(Expt.Stimvals.e2 == etvals);
    if length(id)
        Expt.Stimvals.e2 = etnames{id(1)};
    end
    end
end

Expt.Trials = Expt.Trials(idx);
if(strfind(name,'rds.OxPD'))
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).dO = round(Expt.Trials(j).dO * 100)/100;
        Expt.Trials(j).dP = round(Expt.Trials(j).dP * 100)/100;
    end
end
for j = 1:length(Expt.Trials)
    if ~isfield(Expt.Trials,'dur') || isempty(Expt.Trials(j).dur)
        Expt.Trials(j).dur = Expt.Trials(j).End(end)-Expt.Trials(j).Start(1);
    end
    if ~isfield(Expt.Trials,'id') || isempty(Expt.Trials(j).id)
        Expt.Trials(j).id = j;
    end
end

if ~isfield(Expt.Stimvals,'sM')
    Expt.Stimvals.sM = 0;
end