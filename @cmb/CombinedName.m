function [outname, prefix] = CombinedName(DATA,eid,cluster, varargin)

probe = DATA.probe;
if DATA.listbycell == 1
    cell = DATA.probe;
end

if nargin < 3
    clsuter = 0;
end
modifier = '';

cell = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'cell',4)
        j = j+1;
        cell = varargin{j};
    elseif strncmpi(varargin{j},'modifier',5)
        j = j+1;
        modifier = varargin{j};
    elseif strncmpi(varargin{j},'probe',5)
        j = j+1;
        probe = varargin{j};
    end
    j = j+1;
end

suff = '';
if isempty(DATA.exptypelist{eid(1)})
    outname = [];
    prefix = [];
    return;
end
it = strmatch(DATA.exptypelist{eid(1)},DATA.expstrs,'exact');
%Changed Jan 2012. exptypelist doens't have 'RC' in, so this messes up
%how dit it ever work? ?Caused by Expt2Name
%it = strmatch(DATA.explist{eid(1)},DATA.expstrs,'exact');
if ~isempty(it)
    expname = DATA.expnames{it};
else
    it = strmatch(regexprep(DATA.exptypelist{eid(1)},'RC$',''),DATA.expstrs,'exact');
    if isempty(it)
        expname = DATA.exptypelist{eid(1)};
    else
        expname = [DATA.expnames{it} 'RC'];
    end
end
expname = [modifier expname];

a = strmatch(regexprep(DATA.explist{eid(1)},'\..*',''),DATA.stimnames);
if isempty(a)
    stimname = strrep(DATA.explist{eid(1)},['.' DATA.exptypelist{eid(1)}],'');
else
    stimname = DATA.stimnames{a};
end
if strfind(stimname,'CRC')
    stimname = strrep(stimname,'CRC','');
    suff = 'CRC';
elseif strfind(stimname,'RC')
    stimname = strrep(stimname,'RC','');
    suff = 'RC';
end
if isempty(strfind(expname,stimname))
    expname = [stimname '.' expname];
end

if DATA.listbycell == 1
    cell = DATA.probe;
end
if cell > 0
    cs = ['.cell' num2str(cell)];
elseif DATA.state.includeprobename
    cs = ['.p' num2str(probe) 'c' num2str(cluster)];
else
    cs = ['.c' num2str(cluster)];
end
prefix = DATA.datafilename;
if DATA.state.online
    outname = [DATA.datafilename '/' expname suff cs  '.mat'];
elseif DATA.bysuffix
    outname = [DATA.datafilename '/' DATA.prefix '.' expname suff cs  '.mat'];
    prefix = [DATA.datafilename '/' DATA.prefix];
else
    outname = [strrep(DATA.datafilename,'.mat',cs)  '.' expname suff '.mat'];
end


