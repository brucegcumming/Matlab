function [errs, details] = ShowExptErrs(Expt, varargin)
%ShowExptErrs(Expt, varargin) Prints out errors in struct Expt 
%(or any struct with field errs, so works with anything that uses AdError
%Expt can also be a directoyr name
%see also AddError
markexpt = 'suffix';
j = 1;
Errors = {};
showall = 0;
rebuild = 0;
silent = 0;
details.type = {};

while j <= length(varargin)
    if strncmpi(varargin{j},'label',5)
        j = j+1;
        markexpt = varargin{j};
    elseif strncmpi(varargin{j},'showall',5)
        showall =1;
    elseif strncmpi(varargin{j},'silent',5)
        silent =1;
    elseif iscell(varargin{j})
        Errors = varargin{j};       
    elseif sum(strncmpi(varargin{j},{'check' 'relist' 'rebuild'},5))
        rebuild = 1;
    end
    j = j+1;
end

if iscellstr(Expt)
elseif iscell(Expt)
    allerrs = {};
    for j = 1:length(Expt)
        if isfield(Expt{j},'Header')
            if isempty(Errors) && isfield(Expt{j}.Header,'loadname');
                ErrorFile = [fileparts(Expt{j}.Header.loadname) '/Errors.mat'];
                if exist(ErrorFile)
                    load(ErrorFile);
                end
            end
            errs = ShowExptErrs(Expt{j},Errors,varargin{:});
            allerrs = {allerrs{:} errs{:}};
        elseif iscell(Expt{j}) %A previous set of results
            PrintErrors(Expt{j});
        elseif isfield(Expt{j},'s')
            fprintf('%s\n',Expt{j}.s);
        elseif isfield(Expt{j},'LFPerrs')
            ShowExptErrs(Expt{j});
        end
    end
    errs = allerrs;
    return;
elseif isstruct(Expt) %Drop down to main loop
elseif isdir(Expt)
    ErrorFile = [Expt '/Errors.mat'];
    CErrorFile = [Expt '/ClusterErrors.mat'];
    if exist(CErrorFile)
        load(CErrorFile);
        cid = [errorlist.ex];
    else
        cid = [];
        errorlist = [];
    end
    if exist(ErrorFile) && rebuild == 0
        E = load(ErrorFile);
        if ~silent
            PrintErrors(E.Errors);
        end
        if isempty(cid)
            errs = E.Errors;
        else
            ne = 1;
            for j = 1:length(E.Errors)
                id = find(cid < E.Errors{j}.eid);
                for k = 1:length(id)
                    for e = 1:length(errorlist(id(k)).errs)
                        Errors{ne}.err = errorlist(id(k)).errs{e};
                        Errors{ne}.eid = errorlist(id(k)).ex;
                        Errors{ne}.p = errorlist(id(k)).p;
                        Errors{ne}.type = 'cluster';
                        ne = ne+1;
                    end
                end
                Errors{ne} = E.Errors{j};
                Errors{ne}.type = 'matfile';
            end
            errs = Errors;
        end
    elseif rebuild ==1 % build if missing
        E = PlotErrors(Expt,'checkandsave');
        errs = E.Errors;
    else 
        errs = {};
    end
    details = MakeDetails(errs);
    return;
elseif exist(Expt,'file') % a filename
elseif ischar(Expt) 
    d = mydir(Expt);
    for j = 1:length(d)
        Header.name = d(j).name;
        errs{j} = ShowExptErrs(d(j).name,varargin{:});
        if ~isempty(errs{j}) && iscell(errs{j})
        errs{j} = {Header errs{j}{:}};
        else
            errs{j} = [];
        end
    end
    return;
end

errs = {};
if isfield(Expt,'errs')
    if isfield(Expt,'Header')
        H = Expt.Header;
    else
        H = [];
    end
    eid = GetExptNumber(Expt);
    if strcmp(markexpt,'suffix')
        if isfield(H,'suffix')
            elbl = sprintf('Ex%d:  ',H.suffix);
        else
            elbl = sprintf('Ex%d:  ',eid);
        end
    else
        elbl = [];
    end

    errmsg = {};
    if ~isempty(Errors) %previous list
        for j = 1:length(Errors)
            if isfield(Errors{j},'accepted')
                errmsg{j} = Errors{j}.s;
            end
        end
    end

    
    if iscell(Expt.errs)
        for j = 1:length(Expt.errs)
            errs{j}.s = deblank(Expt.errs{j});
            id = find(strncmp(errs{j}.s,errmsg,length(errs{j}.s)));
            if isempty(id)
                fprintf('%s%s\n',elbl,errs{j}.s);
            else
                errs{j}.accepted = Errors{id(1)}.accepted;
                errs{j}.user = Errors{id(1)}.user;
                if showall
                    cprintf('blue','%s%s:%s says %s\n',elbl,errs{j}.s,errs{j}.user,errs{j}.accepted);
                end
            end
            errs{j}.eid = GetExptNumber(Expt);
            errs{j}.t = 0;         
            errs{j}.type = 'matfile';
        end
    elseif isfield(Expt.errs,'msg')
        for j = 1:length(Expt.errs)
            if iscell(Expt.errs(j).msg)
                errs{j}.s = deblank(Expt.errs(j).msg{1});
            else
                errs{j}.s = deblank(Expt.errs(j).msg);
            end
            errs{j}.t = Expt.errs(j).t;
            errs{j}.eid = GetExptNumber(Expt);
            id = find(strncmp(errs{j}.s,errmsg,length(errs{j}.s)));
            if isempty(id)
                show = 1;
                exstr = [];
            else
                errs{j}.accepted = Errors{id(1)}.accepted;
                errs{j}.user = Errors{id(1)}.user;
                show = showall;
                exstr = sprintf('  %s says %s',errs{j}.user,errs{j}.accepted);
            end
            if show
            if Expt.errs(j).t > 0
                fprintf('%s%s at %.2f%s\n',elbl,deblank(Expt.errs(j).msg),errs{j}.t,exstr);
            else
                if iscell(Expt.errs(j).msg)
                    fprintf('%s%s%s\n',elbl,deblank(Expt.errs(j).msg{1}),exstr);
                elseif ~isempty(Expt.errs(j).msg)
                    fprintf('%s%s%s\n',elbl,deblank(Expt.errs(j).msg),exstr);
                end
            end
            end
            errs{j}.type = 'matfile';
        end
    end
elseif isfield(Expt,'LFPerrs')
    tid = [];
    eid = [];
    elbl =sprintf('E%d',GetExptNumber(Expt));
    for j = 1:length(Expt.LFPerrs)
        s = Expt.LFPerrs{j};
        id = regexp(s,'A.[0-9]+.lfp.mat');
        if ~isempty(id)
            eid(j) = sscanf(s(id(1)+2:end),'%d');
        end
        id = regexp(s,'Trial [0-9]+ ');
        if ~isempty(id)
            t = sscanf(s(id(1)+6:end),'%d');
            tid(j) = t(1); 
        end
    end
    if sum(tid)
        fprintf('%s Has missing LFP data for %d Trials over %d Expts\n',Expt.name, sum(tid > 0),sum(unique(eid) > 0)); 
    end
    for j = 1:length(Expt.LFPerrs)
        if tid(j) ==0
            fprintf('%s\n',Expt.LFPerrs{j});
        end
    end
end
if isfield(Expt,'clustererrs')
    for j = 1:length(Expt.clustererrs)
       
        for k = 1:length(Expt.clustererrs(j).errs)
            errs{end+1}.s = Expt.clustererrs(j).errs{k};
            errs{end}.eid = Expt.clustererrs(j).ex;
            errs{end}.p = Expt.clustererrs(j).p;
            errs{end}.type = 'cluster';
        end
    end
end
details = MakeDetails(errs);

function details = MakeDetails(E)
details = [];
for j = 1:length(E)
    details.type{j} = E{j}.type;
end

function PrintErrors(Errors)

name = '';
elbl = '';
for j = 1:length(Errors)
    if isfield(Errors{j},'eid')
        elbl = sprintf('E%d:',Errors{j}.eid);
    elseif ~isempty(name)
        elbl = ': ';
    else
        elbl = '';
    end
    if isfield(Errors{j},'accepted')
        cprintf('blue','%s%s%s. %s says %s\n',name,elbl,Errors{j}.s,Errors{j}.user,Errors{j}.accepted);
    elseif isfield(Errors{j},'s')
        cprintf('red','%s%s%s\n',name,elbl,Errors{j}.s);
    elseif isfield(Errors{j},'name')
        fprintf('%s\n',Errors{j}.name);
        name = path2name(Errors{j}.name);
    end
end
    