function [nex, Exs] = FindExpt(idxfile, varargin)
%FindExpt(idxfile,'name1',name2,....)
%finds expts matching name in idx file. 
%if idxfile is a cell array, returns a cell array of answer, one for each
%file
%with no search patterns named. returns a list of Expts.

Exs = [];
nex = 0;
tofind = {};
listall = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'list',4)
        listall = 1;
    elseif ischar(varargin{j})
        tofind{length(tofind)+1} = varargin{j};
    end
    j = j+1;
end

if iscell(idxfile)
    Exs = {};
    for j = 1:length(idxfile)
        [nex(j) Exs{j}] = FindExpt(idxfile{j},varargin{:});
    end
    return;
end
if isempty(tofind)
    listall = 1;
end
load(idxfile);


if ~exist('ExptList','var')
    nex = NaN; %indicates no list, not no expts.
    if strmatch('image.orXob',tofind);
        oid = strmatch('or',(Expt.Trials.et));
        bid = strmatch('ob',(Expt.Trials.e2));
        nt = length(intersect(oid,bid));
        if nt > 40
            Exs(1).name = 'images.orXob';
            Exs(1).ids = 0;
            nex = 1;
        end
    end
    return;
end

names = unique({ExptList.expname});
if listall
    n = unique(names);
    for j = 1:length(n)
        fprintf('%s\n',n{j});
    end
end
for j = 1:length(names)
    Exs(j).name = names{j};
    Exs(j).ids = find(strmatch(names{j},{ExptList.expname}));
end
    
if length(tofind)
    ids = [];
    for j = 1:length(tofind)
        found = [];
        for k = 1:length(Exs)
            found(k) = length(strfind(Exs(k).name,tofind{j}));
        end
        ids = [ids find(found)];
    end
    Exs = Exs(ids);
end
nex = length(Exs);