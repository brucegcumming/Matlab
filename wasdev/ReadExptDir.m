function [AllExpts, Ex] = ReadExptDir(name, varargin)
% [Expts, Idx] = ReadExptDir(name, varargin)
% Read all SPike2 .mat files in a directory and combines the Expts lists
% into one list
%
% Called by AplaySpkFile  if first argument is a directory

state.relist = 0;
state.resort = 0; %redo SortExpts, but not whole listing
state.online = 0;
state.quick = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'online',4)
        state.online=1;
    elseif strncmpi(varargin{j},'resort',4)
        state.resort=1;
    elseif strncmpi(varargin{j},'quick',4)
        state.quick=1;
    elseif strncmpi(varargin{j},'relist',4)
        state.relist=1;
    end
    j=j+1;
end

[a,b,c] = GetMonkeyName(name);
outname  = [name '/' a c 'Expts.mat'];
if exist(outname) && state.resort == 0
    load(outname);
    for j = 1:length(Expts)
        if ~isfield(Expts{j}.Header,'suffix')
            Expts{j}.Header.suffix = GetExptNumber(Expts{j});
        end
        if isfield(Expts{j}.Header,'loadname')
            a= fileparts(Expts{j}.Header.loadname);
           loadname = strrep(Expts{j}.Header.loadname,a,name);
           Expts{j}.Header.loadname = loadname;
        end
    end
    AllExpts = Expts;
    Ex = Idx;
    return;
end

AllExpts = {};
%first sort numerically by suffix number
d = mydir([name '/*.mat']);
mnk =GetMonkeyName(name);
for j = 1:length(d)
    if state.online
        if ~isempty(regexp(d(j).filename,'Expt[0-9]*.mat'))
            suffixes(j) = str2double(regexprep(d(j).filename,'Expt([0-9]*).mat','$1'));
        end
    elseif ~isempty(regexp(d(j).name,[mnk '[M,\.,G,0-9]*.mat']))
        suffixes(j) = str2double(regexprep(d(j).filename,'.*[\.,A-z]([0-9]*).mat','$1'));        
    end
end
[a,b] =sort(suffixes);
sid = b(a> 0);


nex = 1;
for j = 1:length(sid)
    fprintf('Reading %s\n',d(sid(j)).name);
    [Ex{nex}, Expts] = APlaySpkFile(d(sid(j)).name,'nospikes','noerrs', varargin{:});
    if ~isempty(Expts)
        AllExpts = {AllExpts{:} Expts{:}};
    end
    if length(Expts) > 1
        fprintf('%d Expts in %s\n',length(Expts),d(sid(j)).name);
    end
    nex = length(AllExpts)+1; %so that it lines up with Expts,
end

for j = 1:length(Ex)
    if isfield(Ex{j},'errs') && ~isempty(Ex{j}.errs)
        cprintf('blue','Errors for %s\n',d(sid(j)).name);
        for k = 1:length(Ex{j}.errs)
            cprintf('red','%s\n',Ex{j}.errs{k});
        end
    end
end
[a,b,c] = GetMonkeyName(name);
outname  = [name '/' a c 'Expts.mat'];
Expts = AllExpts;
Idx = Ex;

if state.quick == 0 %don't write out Expts if didn't load everything
save(outname, 'Expts','Idx');
end

if nargout > 1 %? combine Ex{}
    
end