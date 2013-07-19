function Expt = FixSerialExpt(Expt, varargin)
%Expt = FixSerialExpt(Expt)
%If called with 'readexpt' will return a single expt file for  the whole
%serial file.  Make sure all trials belong to expt

details.exptype = [];
goodonly = 1;
j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})  && isfield(varargin{j},'exptype')
        details = varargin{j};
    elseif strcmp(vararing{j},'useall')
        goodonly = 0;
    end
    j=j+1;
end

%values to keep in each trial even if all the same
keepvals = {'rw'};
if strcmp(details.exptype,'ORBW')
    for j = 1:length(Expt.Trials)
        if isempty(Expt.Trials(j).Seedseq) || isempty(Expt.Trials(j).ob)
            good(j) = 0;
        else
            good(j) = 1;
        end
    end
    Expt.Trials = Expt.Trials(find(good));
    Expt.Trials = rmfields(Expt.Trials,'Stimseq');
elseif strcmp(Expt.Stimvals.et,'pR') || strcmp(details.exptype,'CR') || strcmp(details.exptype,'ICOD')
    xtype = 'pR';
    ytype = 'e0';
    if strcmp(details.exptype,'ICOD')
        xtype = 'ic';
        ytype = 'od'
    end
    Expt.Stimvals.et = xtype;
    Expt.Stimvals.e2 = ytype;
    for j = 1:length(Expt.Trials)
        if isempty(Expt.Trials(j).Stimseq) || isempty(Expt.Trials(j).(xtype))
            good(j) = 0;
        else
            good(j) = 1;
            if isfield(Expt.Trials,'or')
                Expt.Trials(j).ori = Expt.Trials(j).or;
            end
            Expt.Trials(j).or = Expt.Trials(j).Stimseq;
            %Expt.Trials(j).or(Expt.Trials(j).or == -1010);
        end
    end
    Expt.Trials = Expt.Trials(find(good));
    X = rmfields(Expt.Trials,'Stimseq','Seedseq');
    Expt.Trials = rmfields(Expt.Trials,'Stimseq','Seedseq');
elseif strcmp(Expt.Stimvals.et,'Dc') || strcmp(details.exptype,'DCOR') && isfield(Expt.Trials,'Dc') %need Stimseq
    for j = 1:length(Expt.Trials)
        if isempty(Expt.Trials(j).Stimseq) || isempty(Expt.Trials(j).Dc)
            good(j) = 0;
        else
            good(j) = 1;
            if strcmp(Expt.Stimvals.e2,'or')
            Expt.Trials(j).ori = Expt.Trials(j).or;
            Expt.Trials(j).or = Expt.Trials(j).Stimseq;
            %Expt.Trials(j).or(Expt.Trials(j).or == -1010);
            end
        end
    end
    Expt.Trials = Expt.Trials(find(good));
    Expt.Trials = rmfields(Expt.Trials,'Stimseq','Seedseq');
end
Expt.Trials = rmfields(Expt.Trials,'so');
f = fields(Expt.Trials);
for j = 1:length(f)
    if ~sum(strcmp(f{j},keepvals))
        uvals = unique(cat(1,Expt.Trials.(f{j})));
        if length(uvals) == 1 
            Expt.Stimvals.(f{j}) = uvals(1);
            Expt.Trials = rmfield(Expt.Trials,f{j});
        elseif isempty(uvals) % all empty
            Expt.Trials = rmfield(Expt.Trials,f{j});
        end
    end
end


if goodonly && isfield(Expt.Trials,'result')
    good = find([Expt.Trials.result] >= 0);
    Expt.Trials = Expt.Trials(good);
end
if ~isfield(Expt.Stimvals,'e2') || strcmp(Expt.Stimvals.e2,'sn')
    Expt.Stimvals.e2 = 'e0';
end