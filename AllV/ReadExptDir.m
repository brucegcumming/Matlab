function [AllExpts, Ex] = ReadExptDir(name, varargin)
% [Expts, Idx] = ReadExptDir(name, varargin)
% Read all SPike2 .mat files in a directory and combines the Expts lists
% into one list. Saves the result. nameExpts.mat
% ReadExptDir(name, 'relist') Calls APlaySpkFile and rebuilds all the
% idx.mat files
% ReadExptDir(name, 'resort') rebuilds Expt.Trials from the idx file
% 
% ReadExptDir(name, 'showerrs') prints out errors warnings from
% APlaySpkFile/CheckExpts
%
% ReadExptDir(name, 'checkstimtime') Checks duration in serial output file
% Called by AplaySpkFile  if first argument is a directory

state.relist = 0;
state.resort = 0; %redo SortExpts, but not whole listing
state.online = 0;
state.quick = 0;
state.checkstimtime =0;
        state.checkdur=0;
        state.checkframe=0;
        state.showerrs = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'online',4)
        state.online=1;
    elseif strncmpi(varargin{j},'resort',4)
        state.resort=1;
    elseif strncmpi(varargin{j},'checkall',4)
        state.checkdur=1;
        state.checkframe=1;
    elseif strncmpi(varargin{j},'checkstimstime',4)
        state.checkstimtime=1;
    elseif strncmpi(varargin{j},'quick',4)
        state.quick=1;
    elseif strncmpi(varargin{j},'showerrs',6)
        state.showerrs =1;
    elseif strncmpi(varargin{j},'relist',4)
        state.relist = 1;
        state.resort = 1;
    end
    j=j+1;
end

expdates = [];
%first sort numerically by suffix number
%And find modification times of individual expts files
d = mydir([name '/*.mat']);
if isempty(d)
    AllExpts = {};
    Ex = {};
    cprintf('red','No Matlab files in %s\n',name);
    return;
end
mnk =GetMonkeyName(name);
suffixes = [];
expfiles = {};

BinocFile = '';

for j = 1:length(d)
    if state.online
        if ~isempty(regexp(d(j).filename,'Expt[0-9]*.mat'))
            suffixes(j) = str2double(regexprep(d(j).filename,'Expt([0-9]*).mat','$1'));
        end
    elseif ~isempty(regexp(d(j).name,[mnk '[M,\.,G,0-9]*\.[0-9]+.mat']))
        suffixes(j) = str2double(regexprep(d(j).filename,'.*[\.,A-z,0-9]\.([0-9]+).mat','$1'));        
        BinocFile = regexprep(d(j).name,'\.[0-9]*.mat','.bnc.mat');
    elseif ~isempty(regexp(d(j).name,[mnk '[M,\.,G,0-9]+.mat'])) %no .N.mat files. 
%get suffix from first part of file name. Beware of extra files in a folter
%like jbeG068
        [a,b,c] = GetMonkeyName(d(j).name);
        if strfind(d(j).filename,c)
            suffixes(j) = str2double(regexprep(d(j).filename,'.*[\.,A-z]([0-9]+).mat','$1'));
        else
            suffixes(j) = 3000 + str2double(regexprep(d(j).filename,'.*[\.,A-z]([0-9]+).mat','$1'));
        end
    elseif ~isempty(regexp(d(j).name,[mnk '[M,\.,G,0-9]*Expts.mat']))
        e = GetExptNumber(d(j).name);
        if e > 0
            expdates(e) = d(j).datenum;
            expfiles{e} = d(j).name;
        end
    end
end

if isempty(suffixes) %if no suffixes, must be in single files. 
    ne = 0;
    for j = 1:length(d)
        if ~isempty(regexp(d(j).name,[mnk '[M,\.,G,0-9]*idx.mat']))
            e = GetExptNumber(d(j).name);
            ne = ne+1;
            suffixes(j) = ne;
        end
    end
end

[a,b] =sort(suffixes);
sid = b(a> 0 & a < 2000); %get rid of Utah ns5->.mat files 

[a,b,c] = GetMonkeyName(name);
outname  = [name '/' a c 'Expts.mat'];
errfile = [name '/ClusterErrors.mat'];

if state.relist && exist(BinocFile);
    fprintf('Loading %s\n',BinocFile);
    load(BinocFile);
    setappdata(0,'BinocExpt',BinocExpt);
elseif isappdata(0,'BinocExpt')
    rmappdata(0,'BinocExpt');
end


if state.checkstimtime
    CheckStimTimes({d(sid).name});
end
if exist(outname) && state.resort == 0
        cerrexlist = [];
        xd = dir(outname);
    load(outname);
    if exist(errfile)
        load(errfile);
    else
        errorlist = [];
    end
    if isfield(errorlist,'ex')
        cerrexlist = [errorlist.ex];
    end
    if ~exist('Idx','var') %old file
        state.resort = 1;
    else
    reloaded = 0;
    for j = 1:length(Expts)
        if ~isfield(Expts{j}.Header,'suffix')
            Expts{j}.Header.suffix = GetExptNumber(Expts{j});
        end
        e = Expts{j}.Header.suffix;
        if e > 0
            errid = find(ismember(cerrexlist,e));
            if ~isempty(errid)
                Expts{j}.clustererrs = errorlist(errid);
            end
            if e > 0 && length(expdates) >= e && expdates(e) > xd.datenum
                fprintf('Reloading %d\n',e);
                E = load(expfiles{e});
                if ~isempty(E.Expts)
                    Expts{j} = E.Expts{1};
                end
                Idx{j} = E.Tidx;
                reloaded = reloaded+1;
            end
            exptgood(e) = 1;
        else
            exptgood(j) = 1;
        end
        if isfield(Expts{j}.Header,'loadname')
            a= fileparts(Expts{j}.Header.loadname);
           loadname = strrep(Expts{j}.Header.loadname,a,name);
           Expts{j}.Header.loadname = loadname;
        end
        if ~isfield(Expts{j},'errs') && length(Idx) == length(Expts) && isfield(Idx{j},'errs')
            Expts{j}.errs = Idx{j}.errs;
        end
    end
    id = find(exptgood ==0);
    for j = 1:length(id) %Expts not saved in Expts for some reason
        e = id(j);
        if expdates(id(j)) > 0 
            cprintf('red','Adding Expt %d (missing) from %s\n',e,expfiles{e});
            X = load(expfiles{e});
            Expts(e+1:end+1) = Expts(e:end);
            Expts{e} = X.Expts{1};
            reloaded = reloaded+1;
        end
    end
    if reloaded
        [a,b, Expts] = CheckExpts(Expts,'quiet');
        fprintf('Saving Expts with reloaded suffixes\n');
        save(outname,'Expts','Idx');
    elseif state.checkframe
        [a,b, Expts] = CheckExpts(Expts,'quiet');
    else
        [a,b, Expts] = CheckExpts(Expts,'quick');
    end
    AllExpts = SetTrialOffsets(Expts);
    AllExpts = expt.AddComments(AllExpts);
    Ex = Idx;
    if state.showerrs
        ShowExptErrs(AllExpts);
    end
    return;
    end
end

AllExpts = {};
AllErrs = {};
combineexpts = 0;

if isempty(sid) %No files with Expt data
    return;
end
nex = 1;
for j = 1:length(sid)
    fprintf('Reading %s\n',d(sid(j)).name);
    [Ex{nex}, Expts] = APlaySpkFile(d(sid(j)).name,'nospikes','noerrs', varargin{:});
    if length(Expts) > 1 && length(sid) > 2 %Usually this is an error, and the first expt is bad.  
%Could add a test later to try and combine these if possible        
        fprintf('%d Expts in %s\n',length(Expts),d(sid(j)).name);
        if combineexpts
            AllExpts = {AllExpts{:} CombineExpts(Expts)};
        else
            AllExpts = {AllExpts{:} Expts{end}};
        end
    elseif ~isempty(Expts)
        AllExpts = {AllExpts{:} Expts{:}};
    end
    if isfield(Ex{nex},'errs') && ~isempty(Ex{nex}.errs)
        cprintf('blue','Errors for %s\n',d(sid(j)).name);
        for k = 1:length(Ex{nex}.errs)
            cprintf('red','%s\n',Ex{nex}.errs{k});
        end
        if isempty(Expts)
            AllErrs = {AllErrs{:} Ex{nex}.errs};
        elseif isfield(AllExpts{end},'errs')
            if isfield(AllExpts{end}.errs,'msg')
                for k = 1:length(Ex{nex}.errs)
                    AllExpts{end}.errs(end+1).msg = Ex{nex}.errs{k};
                    AllExpts{end}.errs(end+1).t = 0;
                end
            elseif isfield(AllExpts{end}.errs,'errmsg')
                for k = 1:length(Ex{nex}.errs)
                    AllExpts{end}.errs(end+1).errmsg = Ex{nex}.errs{k};
                    AllExpts{end}.errs(end+1).t = 0;
                end
            else
                AllExpts{end}.errs = {AllExpts{end}.errs{:} Ex{nex}.errs{:}};
            end
        else
            AllExpts{end}.errs = Ex{nex}.errs;
        end
    end
    nex = length(AllExpts)+1; %so that it lines up with Expts,
end
AllExpts = SetTrialOffsets(AllExpts);


for j = 1:length(AllExpts)
    e = GetExptNumber(AllExpts{j});
    if isfield(AllExpts{j},'Header') && e > 1
        AllExpts{j}.Header.exptno = e;
    end
    if e > 0 && e <= length(expfiles) && exist(expfiles{e})
        E = load(expfiles{e});
        if isfield(E,'Expts') && ~isempty(E.Expts)
            exptno = GetExptNumber(E.Expts{1});
            if exptno > 1
                E.Expts{1}.Header.exptno = exptno;
            end
            E.Expts{1}.Header.trialoffset = AllExpts{j}.Header.trialoffset;
            save(expfiles{e},'-struct','E');
        end       
    end
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

if ~exist(fileparts(outname))
    cprintf('error','Cannot save %s. Folder does not exist\n',outname);
elseif state.quick == 0 || ~exist(outname) %don't write out Expts if didn't load everything
    [a,b, Expts] = CheckExpts(Expts,'quiet');
    save(outname, 'Expts','Idx','AllErrs');
end

if nargout > 1 %? combine Ex{}
    
end


function CheckStimTimes(name)

if iscellstr(name)
    for j = 1:length(name)
        CheckStimTimes(name{j});
    end
    return;
end

aname = regexprep(name,'\.[0-9]+\.mat','A$0');
X = load(name);
f = fields(X);
for j = 1:length(f)
    V = X.(f{j});
    if isfield(V,'title') && strcmp(V.title,'StimOn')
        stimtimes = V.times;
    end
end
A = load(aname);
f = fields(A);
for j = 1:length(f)
    V = A.(f{j});
    if isfield(V,'title') && strcmp(V.title,'StimOn')
        astimtimes = V.times;
    end
end
lendiff = length(astimtimes) - length(stimtimes);
td = NaN;
if lendiff == 0
    xc = corrcoef(astimtimes, stimtimes);
    if xc(1,2) > 0.99
        td = mean(astimtimes-stimtimes);
    end
elseif lendiff == 1
    xc = corrcoef(astimtimes(2:end), stimtimes);
    if xc(1,2) > 0.99
        td = mean(astimtimes(2:end)-stimtimes);
    end
end
if abs(td) > 0.01
    fprintf('TIme mismatch %.3f in %s\n',td, name);
end