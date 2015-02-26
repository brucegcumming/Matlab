function [durs, details] = CheckStimDur(name,varargin)
% [durs, details] = CheckStimDur(name,varargin) Read online serial file and checks stimulus durations 
%and frame dropping
%
%CheckStimDur(details)  plots up dropped frames

splitfile = 0;
filelist = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'filelist',5)
        filelist = 1;
    elseif strncmpi(varargin{j},'split',4)
        splitfile = 1;
    end
    j = j+1;
end
details = [];


if iscellstr(name)
    if filelist
        for j = 1:length(name)
            [durs{j}, details{j}] = CheckStimDur(name{j});
        end
        return;
    else
        txt = name;
        filename = [];
    end
elseif iscell(name) %Results from multipple files
    CheckResults(name);
    return;
elseif isfield(name,'frametimes') %a single  result
    CheckResult(name);
    return;
elseif isfield(name,'uniquedur') %a single  result
    CheckResult(name);
    return;
elseif isfield(name,'datenum') %a dir result
    for j = 1:length(name)
        [durs{j}, details{j}] = CheckStimDur(name(j).name);
    end
    return;
    
elseif ischar(name)
    txt = scanlines(name);
    filename = name;
end

rundate = 0;
if splitfile
    x = {};
    monkey = 'ica';
    oid = find(strncmp('Reopened',txt,6));
    for j = 1:length(oid)
        dates(j) = datenum(txt{oid(j)}(14:end));
    end
    sameid = 1+find(diff(floor(dates)) < 0.5);
    oid(sameid) = [];  %now list of lines that are new days
    dates(sameid) = [];
    for j = 1:length(oid)
        fname = sprintf('%s%s',monkey,datestr(dates(j),'ddmmmyyyy'));
        x{j} = fname;
        fid = fopen(fname,'w');
        if j < length(oid)
            lastline = oid(j+1)-1;
        else
            lastline = length(txt);
        end
        
        if fid > 0
            for k = oid(j):lastline
                fprintf(fid,'%s\n',txt{k});
            end
            fclose(fid);
        end
         
    end
    durs = x;
    return;
end
    durs = [];
version = 0;
did = find(strncmp('#du',txt,3));
if isempty(did)
    fprintf('No durations in file %s\n',filename);
    return;
end
if regexp(txt{did(1)},'([0-9]+')  %binoc serial output
    fnid = find(strncmp('mtFn=',txt,5));
    ldiff = bsxfun(@minus,fnid,did');
    pid = find(strncmp('#Prep',txt,5));
    startid = find(strncmp('O 5',txt,3));
    for j = 1:length(did)
        a = sscanf(txt{did(j)},'#du%f(%f:%f)');
        if length(a) == 3
            x(j,1:3) = a;
        else %older versions have nominal dur, not frames
            x(j,[1 3]) = a;
        end
        id = find(ldiff(:,j)>0,1);
        if ~isempty(id)
            frametimes{j} = sscanf(txt{fnid(id)}(6:end),'%f');
        end
        a = find(pid < did(j));
        id = find(strncmp('id',txt(pid(a(end)):did(j)),2));
        if ~isempty(id)
            tid(j) = sscanf(txt{pid(a(end))+id(1)-1},'id%d');
        end
        id = find(startid < did(j));
        if ~isempty(id)
            starts(j) = sscanf(txt{startid(id(end))},'O 5 %f');
        end
    end
    [a,b] = Counts(x(:,2));
    [c,d] = max(a);
    id = find(x(:,2) == b(d));
    x = x(id,:);
    if size(x,2) > 2
    [a,b] = Counts(x(:,3)); %check nominal dur in case of change in framerate
    [c,d] = max(a);
    id = find(x(:,3) == b(d));
    x = x(id,:);
    end
    durs = x(:,1);
    frames = x(:,2);
    oid = find(strncmp('Reopened',txt,6));
    for j = 1:length(oid)
        dates(j) = datenum(txt{oid(j)}(14:end));
    end
    if isempty(oid)
        oid = find(strncmp('#SendAll at',txt,6));
        for j = 1:length(oid)
            dates(j) = datenum(txt{oid(j)}(17:end));
        end
    end
    
    veid = find(strncmp('by binoc',txt,6));
    for j = 1:length(veid)
        if strfind(txt{veid(j)},'binoclean')
            version(j) = sscanf(txt{veid(j)}(30:end),'%f');
        else
            version(j) = 0.001;
        end
        details.versionline = txt{veid(j)};
    end
    if ~isempty(oid)
        rundate = dates(1);
    end
else
    waits = {};
    a = textscan(char(txt(did))','#du%f');
    durs = a{1};
    oid = find(strncmp('Reopened',txt,6));
    for j = 1:length(oid)
        dates(j) = datenum(txt{oid(j)}(14:end));
    end
    veid = find(strncmp('by binoclean',txt,6));
    for j = 1:length(veid)
        version(j) = sscanf(txt{veid(j)}(30:end),'%f');
    details.versionline = txt{veid(j)};
    end
    if ~isempty(oid)
        rundate = dates(1);
    end
    fid = find(strncmp('Waits:',txt,6));
    for j = 1:length(fid)
        waits{j} = sscanf(txt{fid(j)}(9:end),'%f');
    end
    fid = find(strncmp('Testmode',txt,6));
    for j = 1:length(fid)
        s = regexprep(txt{fid(j)},'.*wait','');
        maxwaits(j) = sscanf(s,'%f');
    end

    fid = find(strncmp('Frames',txt,6));
    for j = 1:length(fid)
        frametimes{j} = sscanf(txt{fid(j)}(10:end),'%f');
        if frametimes{j}(1) > 1
            frametimess{j}(1) = 0;
        end
        frames(j) = length(frametimes{j});
        dframes{j} = diff(frametimes{j});
    end
    [a,b] =  find(cat(2,dframes{:}) > 2);
    for j = 1:length(a)
        line = fid(b(j));
        fprintf('Line %d %s\n',line+1,txt{line+1});
    end
    
    details.waits = waits;
end    
durs = durs(durs > 0);
[details.counts, details.uniquedur] = Counts(durs);
details.filename = filename;
details.date = dates(1);
details.nt = length(durs);
details.version = mean(version);
details.frames = frames;
details.durs = durs;
        details.frametimes = frametimes;
        details.trialids = tid;
details.starts = starts;
        hist(durs,100);


function CheckResult(X)

colors = mycolors('white');

if isfield(X,'frametimes')
    for j = 1:length(X.frametimes)
        c = 1 + mod(j-1,length(colors));
        plot(X.frametimes{j}-[1:length(X.frametimes{j})]','color',colors{c},'buttondownfcn',{@HitLine, j});
        hold on;
    end
    ylabel('lag (frames)');
elseif isfield(X,'frames') && iscell(X.frames)
end

function HitLine(a,b, trial)
fprintf('Trial %d\n',trial);


function CheckResults(X)

for j = 1:length(X)
    if isfield(X{j},'uniquedur')
        durrange(j,:) = [min(X{j}.uniquedur) max(X{j}.uniquedur)];
    end
    if isfield(X{j},'date')
        dates(j) = X{j}.date;
    end
    if isfield(X{j},'nt')
        ntrials(j) = X{j}.nt;
    end
    if isfield(X{j},'version')
        versions(j) = floor(X{j}.version.*10000)/10000;
    end
end
lid = find(durrange(:,2) > 2.002);
ds = diff(durrange,[],2);
bid = find(diff(durrange) > 0.02);
gid = find(diff(durrange) < 0.016);
id = find(dates > 0);
changedate = datenum('1 Feb 2014');
newid = find(dates > changedate);
oldid = find(dates < changedate & dates > 0);

newid = find(versions == 1.3991);
oldid = find(versions >= 1.3779 & versions < 1.3991);
oldid = find(versions >= 1.35 & versions < 1.3991);
newid = find(versions >= 1.3779);
oldid = find(versions > 0 & versions < 1.3779);

GetFigure('NTrials');
hold off;
plot(ntrials(newid),ds(newid),'o');
hold on;
plot(ntrials(oldid),ds(oldid),'ro');
%set(gca,'ylim',[0 0.2]);
GetFigure('Dates');
plot(dates(id),ds(id),'o');
datetick('x','mmyy');
