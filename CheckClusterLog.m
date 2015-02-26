function errs = CheckClusterLog(name, varargin)
plottype = 'none';
%plottype = 'dates';

modes = {'findmanual'};

errs = [];
verbose = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'verbose')
        verbose = 1;
    end
    j = j+1;
end
if ischar(name) && exist(name,'dir')
    d = mydir([name '/ClusterLogExpt*.mat']);
    if isempty(d)
        d = mydir(name);
        id = find([d.isdir] & ~strncmp('.',{d.filename},1));
        for j = 1:length(id);
            fprintf('Checking %s\n',d(id(j)).name);
            errs{j} = CheckClusterLog(d(id(j)).name);
        end
    else
        for j = 1:length(d)
            errs{j} = CheckClusterLog(d(j).name);
            errs{j}.name = d(j).name;
        end
    end
    return;
elseif ischar(name) && exist(name,'file')
    load(name);
elseif iscell(name) && isfield(name{1},'space') % a clusterlog array
    ClusterLog = name;
elseif iscell(name)  % a previous result
    E = name;
    errs = {};
    for j = 1:length(E)
        if isfield(E{j},'nerr')
            if(E{j}.nerr > 0)
                errs = {errs{:} E{j}};
            end
        elseif isfield(E{j},'name')
        else
            x = CheckClusterLog(E{j});
            if ~isempty(x)
                errs = {errs{:} x{:}};
            end
        end
    end
    return;
elseif isempty(name)
    return;
else
    cprintf('errors','Can''t read %s\n',name);
    return;
end


exptno = 0;
p = CellToMat(ClusterLog,'probe');
[np, probes] = Counts(p);
id = find(np > 2);
errs.nerr = 0;
colors = mycolors;
eid = ClusterLog{end}.exptno;

for j = 1:length(id)
    recluster = [];
    shape = [];
    xyr = [];
    crit = [];
    dates = [];
    user = [];
    err = 0;
    pid = find(p == probes(id(j)));
    for k = length(pid):-1:1
        C = ClusterLog{pid(k)};
        recluster(k)=  C.recluster(1);
        shape(k) = C.shape;
        dates(k) = C.savetime(1);
        username{k} = C.user;
        user(k) = find(strcmp(C.user,{'bgc' 'bondya' 'pillardac' 'moeenya' 'expt'}));
        if C.shape == 0
            xyr(k,:) = C.xyr;
        else
            crit(k) = C.crit;
        end
    end
    uid = find(recluster == 0);
    if ~isempty(uid) && sum(strcmp('findmanual',modes))
        for u = 1:length(uid)
        fprintf('E%dP%d Manual cut by %s on %s\n',eid,probes(id(j)),username{uid(u)},datestr(dates(uid(u))));
        end
    end
    if recluster(end) == 2
        if shape(end) ~= shape(end-1)
            str = 'Shape Changed';
            err = 1;
        elseif shape(end) == 1
            ratio = crit(end)./crit(end-1);
            if abs(ratio-1) > 0.2
                err = 2;
            end
            str = [ 'Ratios:' sprintf(' %.2f',ratio)];
        elseif shape(end) == 0
            ratio = xyr(end,:)./xyr(end-1,:);
            if sum(abs(ratio-1) > 0.2) 
                err = 3;
            end
            str = [ 'Ratios:' sprintf(' %.2f',ratio)];
        else
            str = 'No Shape';
        end
        if err || verbose
        fprintf('E%dP%d Last Cut was reclassify %s\n',C.exptno,probes(id(j)),str);
        end
    end
    errs.errs(probes(id(j))) = err;
    errs.exptno = C.exptno;
    errs.nerr = sum(errs.errs(:) > 0);
    if strcmp(plottype,'dates')
        plot(dates,p(pid),'o','color',colors{j});
        hold on;
        uid = find(recluster == 0); %manual cuts
        plot(dates(uid),p(pid(uid)),'o','color',colors{j},'markerfacecolor',colors{j});        
        datetick;
    end
end