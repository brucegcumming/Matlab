function Expt = LoadExpt(name, varargin);
%Expt = LoadExpt(name,...
%Loads an Expt File, checks for some inconsistencies, and returns;
%Expt = LoadExpt(name,'loadem') Adds eye position data to teh Trials
%Expt = LoadExpt(name,'loadlfp') Adds LFP data to teh Trials.  The length of Trials.LFP is
%   forced to be the same for all trials. Missing samples are filled with NaN
%Expt = LoadExpt(name,'loadlfp','zeropad') fills missing samples with 0

loadem = 0;
loadlfp = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'loadem',6)
        loadem = 1;
    elseif strncmpi(varargin{j},'loadlfp',6)
        loadlfp = 1;
    end
    j = j+1;
end

if ~exist(name,'file')
    Expt = [];
    if(name(1) == '/')
        return;
    else
        uname = regexprep(name,'\\','/');
        if exist(uname,'file')
            load(uname);
        else
            uname = name2path(name);
            if exist(uname,'file')
                load(uname);
            end
        end
    end
    if isempty(Expt)
        return;
    end
elseif isdir(name)
    Expt = [];
    return;
else
    load(name);
end
if ~exist('Expt','var') && exist('AllExpt','var')
    AllExpt.Expt.Header.loadname = name;
    if loadlfp
        AllExpt.Expt = LoadLFP(AllExpt.Expt);
    end
    Expt = AllExpt;
    return;
end
if ~exist('Expt','var') && exist('cExpt','var')
    Expt = cExpt;
end
if ~exist('Expt','var') && exist('muExpt','var')
    Expt = muExpt;
end
Expt.Header.loadname = name;


idx = [];
if ~exist('Expt','var') | ~isfield(Expt,'Trials')
    if exist('cExpt','var') & isfield(cExpt,'Trials')
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


if loadem
    Expt = LoadEmData(Expt);
end

if loadlfp
    Expt = LoadLFP(Expt,varargin{:});
end
