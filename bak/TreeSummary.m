function [res, details] = TreeSummary(path, varargin)
%[res, details] = TreeSummary(path, varargin)
res = [];
details = [];
plottype = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'noplot',5)
        plottype = 0;
    end
    j = j+1;
end

if ischar(path)
    res = TreeList(path);
    if plottype > 0
        [res, details] = FixData(res, varargin{:});
        PlotTree(res, details, varargin{:});
    end
elseif iscellstr(path)
    for j = 1:length(path)
        ress{j} = TreeList(path{j}, 'parent', [1 1+length(res) 1], res);
        res = cat(2,ress{:});
    end
    [res, details] = FixData(res, varargin{:});
    PlotTree(res, details, varargin{:});
else
    [res, details] = PlotTree(path, varargin{:});
    setappdata(gcf,'Files',res);
end


function [T, details] = PlotTree(T, varargin)

nlvl = max([T.level]);
details = [];
settree = 0;
j = 1;
while j <= length(varargin)
    if isstruct(varargin{j}) && isfield(varargin{j},'filetypes')
        details = varargin{j};
    elseif strncmpi(varargin{j},'settree',7)
        j = j+1;
        settree = varargin{j};
    end
    j = j+1;
end


if isempty(details)
    [T, details] = FixData(T, varargin);
end
    details.toplevel = gcf;
set(gcf,'UserData',details);
if settree
    details = SetTree(T, details, settree);
end
PlotFixedData(T, details);


function D = SetTree(T, D, tid)

allid = [];
newid = tid;
id = find(ismember(D.dirid,newid));
D.startlevel = D.levels(id(1));
D.topparent = D.parents(id);
allid = id;
while ~isempty(newid)
id = find(ismember(D.parents,newid));
allid = [allid id]; %allid is indedx in D structs
newid = D.dirid(id); %new id is index of T structs
end
D.sizes = D.sizes(allid,:);
f = {'dirids' 'levels' 'parents' 'dirid'};
for j = 1:length(f)
    D.(f{j}) = D.(f{j})(allid);
end


function [T, details] = FixData(T, varargin)
filetypes = { '.smr' '.ns5' 'FullV.mat'};

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'filetypes',7)
        j = j+1;
        filetypes = varargin{j};
    end
    j = j+1;
end
    
ntypes = 1+length(filetypes);
[T.dsize] = deal(zeros(1,ntypes));
[T.filetype] = deal(0);
nlvl = max([T.level]);
nd = 0;
tic;
for j = 1:length(T)
    for k = 1:length(filetypes)
        if regexp(T(j).names,filetypes{k})
            T(j).filetype = k;
            k = ntypes;  %stop at first match
        end
    end
end
%toc
dids = find([T.size] == -1);
for k = nlvl:-1:1
    lid = find([T.level] ==k);
    ids = unique([T(lid).dirid]);
    for j = 1:length(ids)
        zid = find([T(lid).dirid] == ids(j) & [T(lid).size] ~= 0);
        id = lid(zid);
        xid = find([T(dids).parent] == ids(j));
        xid = dids(xid);
        if ~isempty(id)
        nd = nd+1;
        details.sizes(nd,1) = sum([T(id).size]);
        for f = 1:length(filetypes)
            fid = find([T(id).filetype] == f);
            details.sizes(nd,f+1) = sum([T(id(fid)).size]);
        end
        if ~isempty(xid)
            details.sizes(nd,:) = details.sizes(nd,:) + sum(cat(1,T(xid).dsize),1);
        end
        details.dirids(nd) = ids(j);
        details.levels(nd) = k;
        details.parents(nd) = unique([T(id).parent]);
        did = find([T(dids).dirid] == ids(j));
        did = dids(did);
        if length(did) == 1
            T(did).dsize = details.sizes(nd,:);
            details.dirid(nd) = did;
        else
            fprintf('Missing dir %d\n',ids(j));
            details.dirid(nd) = 1;
        end
        end
    end
end
details.dids = dids;
details.filetypes = filetypes;
details.nlvl = nlvl;


function PlotFixedData(T, D, varargin)

if isfield(D,'startlevel')
    startlevel = D.startlevel;
else
startlevel = 1;
end

j = 1;

while j <= length(varargin)
    if strncmpi(varargin{j},'toplevel')
        j = j+1;
        startlevel = varargin{j};
    end
    j = j+1;
end

oldname = get(gcf,'name');
set(gcf,'name','Busy......');
drawnow;
nlvl = double(D.nlvl);
dids = D.dids;
plot(D.levels,D.sizes,'o');
for k = startlevel:nlvl
    id = find(D.levels == k);
    [a, ids] = sort(D.dirids(id));
    total = [0 cumsum(D.sizes(id(ids),1))'];
    totals(D.dirid((id(ids)))) = total(1:end-1);
    sumtotal(k) = total(end);
    for j = 1:length(ids)
        a = find([T(dids).nf] == T(D.dirid((id(ids(j))))).parent);
        if length(a) == 1
            p(j) = dids(a);
        else
            p(j) = 0;
        end
    end
    bars = [];
    for f = 1:size(D.sizes,2)
        lastp = 0;
        for j = 1:length(ids)
            if p(j) ~= lastp && p(j) > 0 %new parent
                h(1) = totals(p(j));
            elseif p(j) == 0 && j == 1
                h(1) = 0;
            else
                h(1) = h(2);
            end
            w=0.5;
            if f > 2
                h(1) = bars(j,1) + sum(D.sizes(id(ids(j)),2:f-1));
            elseif f == 2
                h(1) = bars(j,1);
            elseif f == 1
                starth(j) = h(1);
                w = 0.5;
            else
            end
            h(2) = h(1)+D.sizes(id(ids(j)),f);
            ph = patch([k k k+w k+w],h([1 2 2 1]),f./size(D.sizes,2));
            if f == 1
                bars(j,:) = h;
                set(ph,'edgecolor','r');
            end
            set(ph,'buttondownfcn',{@HitDir, D.dirid(id(ids(j)))});
            h(1) = h(2);
            lastp = p(j);
        end
    end
    lastp = 0;
    for j = 1:length(ids)
        h(1) = starth(j);
        h(2) = h(1)+D.sizes(id(ids(j)),1);
        th = text(k,mean(h),T(D.dirid(id(ids(j)))).names,...
            'color','k','verticalalignment','middle','horizontalalignment','right',...
            'buttondownfcn',{@HitDir, D.dirid(id(ids(j)))});
        if h(2) > h(1)
           ph = rectangle('Position',[k h(1) 0.5 h(2)-h(1)],'edgecolor','r','linewidth',2);
        end
        if k > startlevel
           ph = line([k k-1+w k],[h(1) mean(h) h(2)]); 
        end
    end
end

ph = getappdata(gcf,'BackupArrow');
if ~isempty(ph) && ph ~= 0 && ishandle(ph)
    delete(ph);
end

if startlevel == 1
    for k = 1:length(D.filetypes)
        h(1) = sumtotal(1) * 1.1;
        h(2) = sumtotal(1) * 1.2;
        patch([k k k+0.1 k+0.1],h([1 2 2 1]),(1+k)./size(D.sizes,2));
        text(k+0.2,mean(h),D.filetypes{k});
    end
else
    ph = annotation('textarrow', [0.2 0], [0.95 0.95]);
    set(ph,'String','Up','buttondownfcn',{@HitDir, D.topparent});
    setappdata(gcf,'BackupArrow',ph);
end
set(gcf,'name',oldname);
%hist(sizes);

function res = TreeList(path, varargin)
% [names, sizes, dates, pathnames] = TreeFind(path, varargin)
% returns a list of files under the tree starting in Path.
%
% TreeFind(path,'name','xxx')   finds files where strfind('xxx',pathname)
% is true
%
%
% TreeFind(path,'function','myfunc')   names a function to be executed. 
% Finds filew with return value >1 finds files where strfind('xxx',pathname)
% TreeFind(path,'function','myfunc',{args})   names a function to be executed. 
% calls myfunc with myfunc(file, args{:});
%
%   ....,'name','.*[0-9][0-9].mat') finds raw spk2 .mat files
name = [];
funcfcn = [];
funcargs = [];
printfiles = 0;
findzero = 1;
found = 0;
newer = 0;
level = int16(1);
res = {};
sizes = [];
dates = [];
extra = [];
parent = int32(0);
nf = int32(1);
thisdir = nf;
j = 1;
while j < nargin
    if strncmpi(varargin{j},'name',3)
        j = j+1;
        name = varargin{j};
    elseif strncmpi(varargin{j},'newer',3)
        j = j+1;
        newer = varargin{j};
    elseif strncmpi(varargin{j},'parent',3)
        j = j+1;
        parent = varargin{j}(1);
        nf = varargin{j}(2);
        level = varargin{j}(3);
        j = j+1;
        res = varargin{j};
        thisdir = length(res);
        nf = length(res)+1;
    elseif strncmpi(varargin{j},'function',3)
        j = j+1;
        funcfcn = varargin{j};
        if length(varargin) > j & iscell(varargin{j+1})
            j = j+1;
            funcargs = varargin{j};
        end
    elseif strncmpi(varargin{j},'print',3)
        printfiles = 1;
    end
    j = j+1;
end
ts = now;
d = dir(path);
if level == 1
    %res(nf).pathname = path;
    res(nf).names = path;
    res(nf).parent = int32(0);
    res(nf).date = 0;
    res(nf).nf = nf;
    res(nf).level = int16(1);
    res(nf).size = -1;
    res(nf).dirid = nf;
    nf = nf+1;
end

%read files then subdirectories, so that res(j) doesn't get set in a
%subdirectory
[a,b] = sort([d.isdir]);
d = d(b);
for j = 1:length(d)
%    res(nf).pathname = [path '/' d(j).name];
    res(nf).names = d(j).name;
    res(nf).parent = parent;
    res(nf).date = d(j).datenum;
    res(nf).nf = nf;
    res(nf).level = level; 
%    res(j).parents = [];
    if d(j).isdir 
        if level == 1 && strcmp('.',d(j).name)
            res(1).nf = thisdir;
            res(1).names = path;
            res(1).date = d(j).datenum;
        elseif ~strcmp('.',d(j).name) && ~strcmp('..',d(j).name)
            fprintf('Reading %s (%d/%d files, %.2f sec) parent is %d\n',[path '/' d(j).name],length(d),nf, mytoc(ts),thisdir);
            res(nf).size = -1;
            res(nf).dirid = nf;
            res(nf).level = level+1;
            res(nf).parent = thisdir;
            %        res(j).parents = [parent thisdir];
            res = TreeList([path '/' d(j).name],'parent',[thisdir nf level+1],res);
            if isfield(res,'nf')
                %res(nf).children = [res(end).nf];
            else
                fprintf('%s/%s empty\n',path,d(j).name);
            end
            nf = length(res)+1;
        else
            res(nf).size = 0;
            res(nf).dirid = nf;
%            res(nf).children = -1;
            nf = nf+1;
        end
    else
        res(nf).size = d(j).bytes;
        res(nf).dirid = thisdir;
   %     res(nf).children = 0;
        nf = nf+1;
    end
    if length(res) > nf
            res(nf).dirid = thisdir;
    end
end


function HitDir(a,b, did)


D = GetDataFromFig(a);
T = getappdata(D.toplevel,'Files');
delete(a);
details = SetTree(T, D, did);
PlotFixedData(T, details);
