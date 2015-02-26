function Mx = CalcEfficacies(C, varargin)
%E= CalcEfficacies(C, varargin) Calc Efficacy for all pairs in a Cluster/Expt set
%E= CalcEfficacies(Expts,) works on expt set - these should he different
%cells,same expt
%E= CalcEfficacies(Expts,) does only adjacent pais

adjacentonly = 0;
verbose = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'adjacent',5)
        adjacentonly = 1;
    elseif strncmpi(varargin{j},'verbose',5)
        adjacentonly = 1;
    end
    j = j+1;
end

nc = 0;
if iscell(C) && isfield(C{1},'Trials') %expt list
    E = C;
    C = expt.spktimes(E);
elseif isfield(C,'Spikes') && iscell(C.Spikes) %AllExpt struct
    E = All2Expt(C,'all');
    C = expt.spktimes(E);
end
for j = 1:length(C)
    if isnumeric(C{j}) %cell array of timee
        nc = nc+1;
        spkt{nc} = C{j};
    else
        if isfield(C{j},'times')
            nc = nc+1;
            spkt{nc} = C{j}.times;
        end
    for c = 1:length(C{j}.next)
        if isfield(C{j}.next{c},'times') && ~isempty(C{j}.next{c}.times)
            nc = nc+1;
            spkt{nc} = C{j}.next{c}.times;
        end
    end
    end
end

if adjacentonly
for j = 1:length(spkt)-1
    e = CalcEfficacy(spkt{j},spkt{j+1});
    Mx(j,1) = e(1);
    Mx(j,2) = e(2);
end
else
for j = 1:length(spkt)
    if verbose
        fprintf('Cell %d\n',GetCellNumber(E{j}));
    end
    for k = 1:j-1
        e = CalcEfficacy(spkt{j},spkt{k});
        Mx(j,k) = e(1);
        Mx(k,j) = e(2);
    end
end
end