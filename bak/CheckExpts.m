function [good, res] = CheckExpts(Expts, varargin)

%CheckExpts(Expts, varargin)
%Checks an experiment, or cell struct list of expts, for timing accuracy,
%etc
%CheckExpts(Expt,'jxnz') returns 1 if jx > 0
%CheckExpts(Expt,'type') returns expt types
%

j = 1;
while j <= length(varargin)
    j = j+1;
end
if iscell(Expts)
    for j = 1:length(Expts)
        res{j} = CheckOneExpt(Expts{j},varargin);
        good(j) = res{j}.good;
    end
elseif isstruct(Expts) || ischar(Expts)
     res = CheckOneExpt(Expts,varargin);    
     good = res.good;
end




function res = CheckOneExpt(Expt, varargin)
plotdelays = 0;
checkmode = 1;
checkjx = 0;
res.good = 0;
verbose = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'plot',4)
        plotdelays = 1;
    elseif strncmpi(varargin{j},'clusters',4)
        checkmode = 4;
    elseif strncmpi(varargin{j},'jxnz',4)
        checkmode = 2;
    elseif strncmpi(varargin{j},'type',4)
        checkmode = 8;
    elseif strncmpi(varargin{j},'orbw',4)
        checkmode = 16;
    elseif strncmpi(varargin{j},'psign',4)
        checkmode = 32;
    end
    j = j+1;
end

if ischar(Expt)
    load(Expt)
end


if checkmode == 1
    id = find(~isnan([Expt.Trials.delay]));


    for j = id
        durs(j) = Expt.Trials(j).TrueEnd(1) - Expt.Trials(j).Start(1);
    end

    if plotdelays
        hist(durs(id));
    end
    fprintf('Delay %.0f - %.0f, dur %.0f - %.0f\n',min([Expt.Trials(id).delay]),max([Expt.Trials(id).delay]),...
        min(durs(id)),max(durs(id)));
end

if checkmode == 4
    res.nc = 0;
    res.np = 0;
    if isfield(Expt,'Cluster')
        res.nc = size(Expt.Cluster,1);
        res.np = size(Expt.Cluster,2);
        for j = 1:res.np
            nc = 0;
            for k = 1:res.nc
                if isfield(Expt.Cluster{k,j},'x')
                    nc = nc+1;
                end
            end
            res.pcls(j) = nc;
        end
    end
    if res.nc
        res.good = res.nc;
    end
    if verbose
        fprintf('%d,%d clusters/probes\n',res.nc,res.np);
    end
end
if bitand(checkmode,2)
    jx = GetEval(Expt,'jx');
    if verbose
    fprintf('%s jx is %.3f\n',Expt.Header.Name, jx);
    end
    res.jx = jx;
    if jx > 0
    res.good = 1;
    else
        res.good = 0;
    end
end

if bitand(checkmode,16)
    fprintf('%s SF%.2f,sz%.2f: ',Expt.Header.Name,GetEval(Expt,'sf'),GetEval(Expt,'sz'));
    if ~isfield(Expt.Trials,'or')
        fprintf(' No or in TrialsExpt %sX%s\n',Expt.Stimvals.et,Expt.Stimvals.e2);
        return;
    end
    [a,b] = Counts([Expt.Trials.or]);
    for j = 1:length(a)
        fprintf('%.0f(%d),',b(j),a(j));
    end
    if ~isfield(Expt.Trials,'ob')
        fprintf(' No ob in Trials: Expt %sX%s\n',Expt.Stimvals.et,Expt.Stimvals.e2);
        return;
    end
    n = sum([Expt.Trials.ob] > 120);
    fprintf('Inf %d\n',n)
    res.good = 1;
end

if bitand(checkmode, 32)
    et = Expt.Stimvals.et;
    nid = find([Expt.Trials.RespDir] < 0);
    pid = find([Expt.Trials.RespDir] > 0);
    nx = mean([Expt.Trials(nid).(et)]);
    px = mean([Expt.Trials(pid).(et)]);
    ori = GetEval(Expt,'or');
    d = -sin((ori+135) * pi/180);
    td = sign(d) .* sign(px);
    fprintf('Or %.2f (%.3f) RespDir 1 %s %.3f, Respdir -1 %s %.3f %.0f\n',ori,d,et,px, et,nx,td);
    
end