function res = CheckAllExpts(AllExpts,varargin)
dpmin = 2;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'dprime',2)
        j = j+1;
        dpmin= varargin{j};
    end
    j = j+1;
end

res = {};
for j = 1:length(AllExpts.exps)
    if isfield(AllExpts.exps{j},'spkres')
    for k = 1:length(AllExpts.exps{j}.spkres)
        res{j,k} = CheckSpkRes(AllExpts.exps{j}.spkres{k},AllExpts.names{j},dpmin);
    end
    end
end

function r = CheckSpkRes(spk, name, dpmin)

r = [];
if isempty(spk)
    return;
end

if isfield(spk,'Data')
    Expt = spk.Data;
    Header = Expt.Header;
else
    Header = spk.Header;
    Expt = [];
end
cell = spk.cellnumber;
Header.BlockStart = [Header.BlockStart spk.Trials(end)];
if isfield(Header, 'filename')
[a,expname] = fileparts(Header.filename);
else
    expname = name;
end

if length(spk.probes) == 1
    probes = ones(size(Header.Clusters)) .* spk.probes;
    tp(1:length(spk.Trials)) = spk.probes;
else
    k =1;
    for j = 1:length(spk.probestep)
        tp(k:spk.probestep(j)) = spk.probes(j);
        k = spk.probestep(j);
    end
    tp(k:length(spk.Trials)) = spk.probes(end);
    for j = 1:length(Header.Clusters)
        trials = Header.BlockStart(j):(Header.BlockStart(j+1)-1);
        probes(j)= mode(tp(ismember(spk.Trials,trials)));
    end
end
if isfield(Header,'BlockCount')
    id = find(~isnan(Header.BlockCount));
    if length(id)
    bc = Header.BlockCount(id);
    if max(bc) > min(bc) * 2
        [a,b] = max(bc);
        [c,d] = min(bc);
        fprintf('%s Counts %d:%.2f %d:%.2f, mean %.2f\n',expname,id(b),a,id(d),c,mean(bc));
    end
    end
end
for j = 1:length(Header.Clusters)
    C = Header.Clusters{j}{1,probes(j)};
    bt = Header.BlockStart(j):(Header.BlockStart(j+1)-1);
    if sum(ismember(spk.Trials,bt)) %only check blocks where cell is defined
    if isfield(C,'dprime')
    dprime(j) = C.dprime;
    if dprime(j) < dpmin || isnan(dprime(j))
        fprintf('%s Dprime %.2f Cell %d probe %d, Block %d(%d)\n',expname,dprime(j),cell,probes(j),j,Header.Combineids(j));
    end
    else
        dprime(j) = NaN;
        fprintf('No Dprime %s Cell %d probe %d, Block %d(%d)\n',expname,cell,probes(j),j,Header.Combineids(j));
    end
    end
end
r.dprime = dprime;
        
        