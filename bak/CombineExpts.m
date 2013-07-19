function Expt = CombineExpts(varargin)

j = 1;
nx = 0;
while j <= length(varargin)
    if iscell(varargin{j}) && isstruct(varargin{j}{1})
        Expts = varargin{j};
        nx = length(Expts);
    elseif isstruct(varargin{j}) & isfield(varargin{j},'Trials');
        nx = nx+1;
        Expts{nx} = varargin{j};
    end
    j = j+1;
end

goodx = [];
for j = 1:nx
    if ~isempty(Expts{j})
        goodx = [goodx j];
    end
end
nx = length(goodx);
Expts = Expts(goodx);

for k = 1:nx
    if ~isfield(Expts{k}.Trials,'Trial')
        fprintf('Filliing Trial\n');
        for j = 1:length(Expts{k}.Trials)
            Expts{k}.Trials(j).Trial = j;
        end
    end
end
fn = fields(Expts{1}.Stimvals);
diffs = zeros(size(fn));
for j = 1:length(fn)
    for k = 2:nx
        if isnumeric(Expts{1}.Stimvals.(fn{j})) & isfield(Expts{k}.Stimvals,fn{j})
        if Expts{k}.Stimvals.(fn{j}) ~= Expts{1}.Stimvals.(fn{j})
            diffs(j) = 1;
        end
        end
    end
end
id = find(diffs);
for j = 1:length(id)
    fprintf('Filling %s\n',fn{id(j)});
    for k = 1:nx
        Expts{k} = FillTrials(Expts{k},fn{id(j)});
    end
end


if ~isfield(Expts{1}.Header,'BlockStart')
    Expts{1}.Header.BlockStart = Expts{1}.Trials(1).Trial;
end
Expt = Expts{1};
for x = 2:nx;
    
    nt = Expt.Trials(end).Trial;
    ni = length(Expt.Trials);
    for j = 1:length(Expts{x}.Trials)
        Expts{x}.Trials(j).Trial = Expts{x}.Trials(j).Trial + nt;
        f = fields(Expts{x}.Trials);
        for k = 1:length(f)
        Expt.Trials(ni+j).(f{k}) = Expts{x}.Trials(j).(f{k});
        end
    end
%    Expt.Trials = [Expt.Trials Expts{x}.Trials];
     if isfield(Expts{x}.Header,'BlockStart')
         Expt.Header.BlockStart = [Expt.Header.BlockStart Expts{x}.Header.BlockStart+nt];
     else
         Expt.Header.BlockStart = [Expt.Header.BlockStart nt];
     end
end
Expt.Header.BlockStart(end) = length(Expt.Trials);


      