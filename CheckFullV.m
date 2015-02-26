function result = CheckFullV(FullV, varargin)
% CheckFullV(FullV, varargin)
% Check(FullV, Expt)
%Check Integrity of data in FullV Files. Looks for epochs were there are
%missing channels, based on low amplitudes in probes1-16, or 17-24% CheckFullV(FullV, varargin)
%Check Integrity of data in FullV Files. Looks for epochs were there are
%missing channels, based on low amplitudes in probes1-16, or 17-24
%
%If FullV is a directory name, checks all FullVs in the folder
%if FullV is a cell string array, checks each file/folder in arrar
%If FullV is a cell array of FullVs checks each
%
%Any Errors are recorded in Errors.mat 



result = fullv.Check(FullV, varargin{:});
return;

ispk = 22;
scales = 1;
spkrate = 100;
showth = 0;
filtershape = 0;
maxsize = NaN;
Expt = [];
Expts = {};
j = 1;

while j <= length(varargin)
    if isstruct(varargin{j})
        if isfield(varargin{j},'Header')
            Expt = varargin{j};
        end
    elseif iscell(varargin{j})
        if sum(CellToMat(varargin{j},'Header.frameperiod'))
            Expts = varargin{j};
        end
    elseif strncmpi(varargin{j},'maxsize',5)
        j = j+1;
        maxsize = varargin{j};
    elseif strncmpi(varargin{j},'odd',3)
        filtershape = 1;
    elseif strncmpi(varargin{j},'scale',4)
        j = j+1;
        scales = varargin{j};;
    elseif strncmpi(varargin{j},'tchan',4)
        j = j+1;
        ispk = varargin{j};
    end
    j = j+1;
end


if iscellstr(FullV)
    names = FullV;
    if ~isempty(Expts)
        expts = GetExptNumber(Expts);
    else
        expts = [];
    end
    for j = 1:length(names)
        if isdir(names{j})
            result{j} = CheckFullV(names{j},varargin{:});
        else
            go = 1;
            if maxsize > 0
                d = dir(names{j});
                if d.bytes > maxsize
                    fprintf('Ignoring %s size is %d\n',names{j},d.bytes);
                    go = 0;
                end
            end
            if go
                e = GetExptNumber(names{j});
                a = find(expts == e);
                if ~isempty(a)
                    Expt = Expts{a};
                else
                    Expt = [];
                end
                fprintf('Loading %s\n',names{j});
                FullV = LoadFullV(names{j});
                result{j} = CheckFullV(FullV,Expt,varargin{:});
                clear FullV;
            end
        end
    end
    return;
elseif iscell(FullV) %array of results
    if sum(CellToMat(FullV,'checktime'))
        ShowResults(FullV);
        return;
    end
elseif ischar(FullV)
    if isdir(FullV)
        d = mydir([FullV '/*FullV.mat']);
        Expts = ReadExptDir(FullV);
        if isempty(d)
            result.name = FullV;
            result.errs = 'No FullVs';
        else            
            result = CheckFullV({d.name},Expts);
        end
    else
        fprintf('Checking %s\n',FullV);
        X = LoadFullV(FullV);
        result{j} = CheckFullV(X,varargin{:});
    end
    return;
end

result.checktime = now;
if ~isfield(FullV,'blkend') && isfield(FullV,'blkstart')
for j = 1:length(FullV.blkstart)
    FullV.blkend(j) = FullV.blkstart(j)+ FullV.blklen(j).*FullV.samper;
end
end
if ~isfield(FullV,'V') || size(FullV.V,1) < 2 %may be FullVdata from cluster file or single channel from Utah
    if isfield(FullV,'blkstart')
        if ~isfield(FullV,'V')
            fprintf('No Voltages - checking blkstarts against Expt\n');
        elseif ~isfield(FullV,'t')
            fprintf('Building times list from blkstarts\n');
        end
        t = BuildFullVt(FullV);
        [a,result.missingids{1}] = FindMissingTrials(Expt, t);
        result.missingid = [];  % can't do individual probes
        if ~isempty(a)
            fprintf('missing ids%s\n',sprintf(' %d',result.missingids{1}));
        end
    else
        fprintf('No Voltages to check\n');
    end
    return;
end
result.name = FullV.loadname;
ends = cumsum(FullV.blklen);
result.exptno = FullV.exptno;
starts = [1 1+ends(1:end-1)];
errorfile = [];
if isfield(FullV,'loadname')
    errorfile = [fileparts(FullV.loadname) '/Errors.mat'];
    BackupFile(errorfile,'copy');
end
%When thre are large artifacts, std meanV can be bigger than individual
%channels. loa
if isfield(FullV,'meanV')
    crit = std(FullV.meanV)./10;
%add back in hte mean voltage for channels where variance is bigger
%in case the channel is actually blank
    if isfield(FullV,'meangain')
    tsd = std(FullV.V,[],2);
    for j = 1:size(FullV.V,1)
        if tsd(j) > crit
            FullV.V(j,:) = FullV.V(j,:) + FullV.meanV;
        end
    end
    end
else 
    crit = 0.01;
end
for j = 1:length(FullV.blklen)
    sds(j,:) = std(FullV.V(:,starts(j):ends(j))');
end
std(sds);
[a,b] = find(sum(sds(:,1:16) < crit,2) > 14);
err(1) = 0;
if ~isempty(Expt)
    [result.missingtrials, result.missingid] = FindMissingTrials(Expt,FullV.t);
    for j = 1:length(Expt.Trials)
        trials(j,1) = Expt.Trials(j).Start(1)./10000;
        trials(j,2) = Expt.Trials(j).End(end)./10000;
        trialids(j) = Expt.Trials(j).id;
    end
else
    trials = [];
end

if ~isempty(a)
    result = AddError(result,'-write',errorfile,'Missing probes 1:16 in Expt%d Blocks %s',FullV.exptno,sprintf(' %d',a));    
    err(1)  = 1;
    if ~isempty(Expt)
        result.missinggroup{1} = FindTrialsInBlocks(FullV,a, trials);
        result.missingids{1} = trialids(result.missinggroup{1});
    end
end
[a,b] = find(sum(sds(:,17:24) < crit,2) > 7);
if ~isempty(a)
    result = AddError(result,'-write',errorfile,'Missing probes 17:24 in Expt%d Blocks %s',FullV.exptno,sprintf(' %d',a));
    result.missinggroup{2} = FindTrialsInBlocks(FullV,a, trials);
    result.missingids{2} = trialids(result.missinggroup{2});
    err(2) = 1;
end
if sum(err) == 0 && sum(sds(:) < crit)
    [a,b] = find(sds < crit); %prob missing
    bad = unique(b); %blocks with errors
    for j = 1:length(bad)
        id = find(b == bad(j));
        result = AddError(result,'-write',errorfile,'FullV Flat in Blocks %s in Expt%d Probe %d',sprintf(' %d',a(id)),FullV.exptno,bad(j));
    end
    [a,b] = find(sds < crit); %low amp
    weak = unique(b); %blocks with errors
    weak = setdiff(weak,bad);
    for j = 1:length(weak)
        id = find(b == weak(j));
        result = AddError(result,'-write',errorfile,'FullV small in Blocks %s in Expt%d Probe %d',sprintf(' %d',a(id)),FullV.exptno,bad(j));
    end
end
fullv.SaveErrors(FullV); %also record any in the file

function tlist = FindTrialsInBlocks(FullV, bid, trials)

tlist = [];
if isempty(trials)
    return;
end
for j = 1:length(bid)
    tid = find(trials(:,1) > FullV.blkstart(bid(j)) & trials(:,2) < FullV.blkend(bid(j)));
    tlist = [tlist tid'];
end


function ShowResults(res)

for j = 1:length(res)
    if iscell(res{j})
        ShowResults(res{j});
    else
        if isfield(res{j},'errs') 
            for e = 1:length(res{j}.errs)
                fprintf('%s:%s\n',res{j}.name,res{j}.errs{e});
            end
        end
    end
end