function Data = CountPsychTrials(Trials, varargin)

Trials = Trials(find([Trials.score] ~= 0));
sortlist = ones(size(Trials));
listlen = [];
nmin = 1;
skip = 1;
last = length(Trials);
j = 1;
while j < nargin
    if(strncmpi(varargin{j},'sort',4))
        j = j+1;
        if ~isempty(varargin{j})
            sortby = varargin{j};
            sortlist = eval(['[Trials.' sortby ']']);
        end
    elseif(strncmpi(varargin{j},'nmin',4))
        j = j+1;
        nmin = varargin{j};
    elseif(strncmpi(varargin{j},'expts',4))
        j = j+1;
	explist = sort(unique(varargin{j}));
    id = find(explist > length(Trials));
    explist(id) = length(Trials);
	listlen(1,:) = explist;
	listlen(2,:) = [explist(2:end) length(Trials)];
    elseif(strncmpi(varargin{j},'skip',4))
        j = j+1;
        skip = varargin{j};
    elseif(strncmpi(varargin{j},'last',4))
        j = j+1;
        last = varargin{j};
	if(last <= 1 | last > length(Trials))
	  last = length(Trials);
	end
    end
    j = j+1;
end

if skip < 1
    skip = 1;
end
InTrials = Trials;

if isempty(listlen)
  listlen(1,1) = skip;
  listlen(2,1) = last;
end

n = 1;
expno = 1;
for nlist = 1:size(listlen,2)
Trials = InTrials(listlen(1,nlist):listlen(2,nlist));
respdir = [Trials.score] .* [Trials.sign];
for y = unique(sortlist);
    if exist('sortby','var')
        idx = eval(['find([Trials.' sortby '] == y)']);
        nTrials = Trials(idx);
        resps = respdir(idx);
    else
        nTrials = Trials;
        resps = respdir;
    end
    ntrials = length(nTrials);

 
for x = unique([nTrials.x])
    Data(n).x = x;
    idx = find([nTrials.x] == x);
    if ~isempty(idx)
        Data(n).n = length(idx);
        Data(n).resp = length(find(resps(idx) > 0));
        Data(n).p = Data(n).resp/Data(n).n;
        Data(n).expno = expno;
        if exist('sortby','var')
            eval(['Data(n).' sortby '= y;']);
            Data(n).name = sprintf('%s = %.3g %d Trials (%d - %d)',sortby,y,ntrials,listlen(1,nlist),listlen(2,nlist));
        else
            Data(n).name = 'x';
        end
        if length([nTrials(idx).xo]);
            Data(n).xo = mean([nTrials(idx).xo]);
            n = n+1;
        end
    end
end
    expno = expno+1;
end
end

Data = Data(find([Data.n] >= nmin));