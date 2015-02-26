function [good, res, Expts] = CheckExpts(Expts, varargin)
%[good, res] CheckExpts(Expts, varargin)
%Checks an experiment, or cell struct list of expts, for timing accuracy,
%etc
%also fixes any errors in structure caused by old prog versions
%CheckExpts(Expt,'jxnz') returns 1 if jx > 0
%CheckExpts(Expt,'type') returns expt types
%
%returns a vector with 0/1 indicating inconsistent/consistent durations

plotdelays = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'plot',4)
        plotdelays = 1;
    end
    j = j+1;
end
errdata.errs = {};
errdata.errdata = [];

if plotdelays
    GetFigure('DurSeq');
    hold off;
end

if iscell(Expts)
    default.Header.nexpts = length(Expts);
    for j = 1:length(Expts)
        default.Header = CopyFields(default.Header,Expts{j}.Header,'DataType');
        [res{j}, Expts{j}] = CheckOneExpt(Expts{j},default,varargin{:});
        good(j) = res{j}.good;
        if isfield(res{j},'errs')
            errdata.errs = {errdata.errs{:} res{j}.errs{:}};
            errdata.errdata = [errdata.errdata res{j}.errdata];
        end
    end
elseif isstruct(Expts) || ischar(Expts)
     res = CheckOneExpt(Expts,varargin);    
     good = res.good;
end




function [res, Expt] = CheckOneExpt(Expt, default, varargin)
plotdelays = 0;
checkmode = 1;
checkjx = 0;
res.good = 0;
checkframes = 0;

verbose = 2;

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
    elseif strncmpi(varargin{j},'quiet',5) %errors only
        verbose = 1;
    elseif strncmpi(varargin{j},'silent',4) %no printing
        verbose = 0;
    elseif strncmpi(varargin{j},'quick',5) %just check structure integrity
        checkmode = 0;
    elseif strncmpi(varargin{j},'psign',4)
        checkmode = 32;
    elseif strncmpi(varargin{j},'readmethod',6)
        checkmode = 64;
    elseif strncmpi(varargin{j},'framedrop',6)
        checkframes = 1;
    end
    j = j+1;
end

if ischar(Expt)
    load(Expt)
end


if isfield(Expt,'Header') && isfield(Expt.Header,'exptno')
    eid = Expt.Header.exptno;
else
    eid = 0;
end
res.exptno = GetExptNumber(Expt);
if ~isfield(Expt.Header,'DataType')
    if isfield(default.Header,'DataType')
        Expt.Header.DataType = default.Header.DataType;
    end
end

if isfield(Expt,'Comments') && isfield(Expt.Comments,'times') && iscell(Expt.Comments.times)
    if length(Expt.Comments.times) == 1
        Expt.Comments.times = Expt.Comments.times{1};
    end
end


if checkmode == 1 && isfield(Expt.Trials,'delay')
    fz = Expt.Header.frameperiod;
    id = find(~isnan([Expt.Trials.delay]));
    delays = [Expt.Trials(id).delay];

    for j = id
        durs(j) = Expt.Trials(j).TrueEnd(1) - Expt.Trials(j).Start(1);
    end
    xstr = '';
    if diff(minmax(durs)) > fz %> 1 frame scatter
        d = median(durs);
        b = sum(abs(durs-d)>fz/2);
        xstr = sprintf(' %d outliers',b);
        Expt = AddError(Expt,'-silent','Expt%d Duration varies: %d/%d outliers',GetExptNumber(Expt),b,length(id));
    end
    s = sprintf('Expt%d (%s %d trials) Delay %.0f - %.0f, dur %.0f - %.0f%s\n',GetExptNumber(Expt),Expt2Name(Expt),length(id),min([Expt.Trials(id).delay]),max([Expt.Trials(id).delay]),...
            min(durs(id)),max(durs(id)),xstr);

        if plotdelays
            GetFigure('DurHist');
            hist(durs(id));
            GetFigure('DurSeq');
            plot(durs(id)-median(durs));
            hold on;
        end
        if verbose == 0
        elseif diff(minmax(delays)) > fz/2
        cprintf('blue',s);
    elseif diff(minmax(durs)) > fz
        cprintf('red',s);
    elseif verbose > 1
        fprintf('%s',s);
    end
end

if isfield(Expt.Trials,'Result') && sum([Expt.Trials.Result] ==0)
    bid = find([Expt.Trials.Result] == 0);
    gid = find(ismember([Expt.Trials.Result],[-1 1]));
    res.gooddur = sum([Expt.Trials(gid).dur]);
    res.baddur = sum([Expt.Trials(bid).dur]);
    res.badratio = res.baddur./(res.gooddur+res.baddur);
    
end

if checkframes
    colors = mycolors;
    ndrop = 0;
    for j = 1:length(Expt.Trials)
        adjust = 0;
        if isfield(Expt.Trials,'mtFi') && length(Expt.Trials(j).mtFi) > 10
            hold off;
            [a,b] = CheckFrameDiffs(Expt.Trials(j).mtFi,'plot');
            errs(j) = b;
            expectdur(j) = 10 +sum(Expt.Trials(j).mtFi);
            if isfield(Expt.Trials,'End')
                truedur(j) = (Expt.Trials(j).End(end) - Expt.Trials(j).Start(1))/10; %in ms
            else
                truedur(j) = Expt.Trials(j).duration .* 1000;
            end
            diffdur = truedur(j) - (1+length(Expt.Trials(j).mtFi)) .* 10;
            line(get(gca,'xlim'),[diffdur diffdur]);
            if isempty(a)
                firstframe(j) = 0;
            else
                firstframe(j) = a(1);
            end
        else
            expectdur(j) = NaN;
            truedur(j) = NaN;
        end
        if isfield(Expt.Trials,'rptframes') && ~isempty(Expt.Trials(j).rptframes)
            ndrop = ndrop+1;
            res.droptrials(ndrop) = j;
            res.firstdrop(ndrop) = Expt.Trials(j).rptframes(1);
            if Expt.Trials(j).rptframes(1) == 1
                adjust = -1;
            end
            res.dropids(ndrop) = Expt.Trials(j).id;
            res.dropcount(ndrop) = length(Expt.Trials(j).rptframes);
        elseif isfield(Expt.Trials,'mtFn') && ~isempty(Expt.Trials(j).mtFn)
            ndrop = ndrop+1;
            res.droptrials(ndrop) = j;
            res.dropids(ndrop) = Expt.Trials(j).id;
            ft = Expt.Trials(j).mtFn;
            frames = 1:length(ft);
            c = colors{1+mod(j,length(colors))};
            plot(frames(:)-ft(:),'color',c);
            hold on;
            if isfield(Expt.Trials,'framerpt')
                for k = 1:length(Expt.Trials(j).framerpt)
                    f = Expt.Trials(j).framerpt(k);
                    plot([f f],get(gca,'ylim'),'color',c);
                end
            end
        else
            res.goodtrials(j) = 1;
        end
        if isfield(Expt.Trials,'duration')
            res.dur(j) = Expt.Trials(j).duration;
        else
        res.dur(j) = Expt.Trials(j).End(end) - Expt.Trials(j).Start(1);        
        res.ndrops(j) = length(Expt.Trials(j).rptframes);
        res.expectdur(j) = (40000 + length(Expt.Trials(j).rptframes) .* 100);
        res.err(j) = res.dur(j) - res.expectdur(j);
        end
        res.adjust(j) = adjust;
    end
    if ndrop
        fprintf('%d/%d trials with dropped frames in %s\n',ndrop,length(Expt.Trials),Expt2Name(Expt));
    end
    title(sprintf('%d/%d Trials with Drops',ndrop, j));
    GetFigure('FrameDrop');
    id = find(~isnan(truedur));
    hold off; plot(truedur(id),expectdur(id),'o');
    refline(1);
    res.ndrop = ndrop;
    return;
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

if bitand(checkmode,64) %Check Read Mode
    if ~isfield(Expt.Header,'ReadMethod')
        mycprintf('red','No ReadMethod in %s E%d\n',Expt.Header.ReadMethod,eid)
    elseif Expt.Header.ReadMethod < 1
        mycprintf('red','Old ReadMethod in %s E%d\n',Expt.Header.ReadMethod,eid)
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