function fits = BuildRFFits(dirname, varargin)
%find OP/PP expts in a dir and fit x,Y, positions
%fits = BuildRFFits(dirname, varargin)
%currently crashed on M009
%if dirname is a cell array, fits is a cell array of fits
%See also PlotRFFits, CombineRFData
savefits = 0;
autosave = 1;
refit = 0;
parallel = 0;
relist = 0;
oldfits = [];
args ={};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'refit',4)
        refit = 1;
        args = {args{:} 'rebuild'};
    elseif strncmpi(varargin{j},'relist',4) %just run combinefits again
        relist = 1;
    elseif strncmpi(varargin{j},'check',4)
        j = j+1;
        oldfits = varargin{j};
    elseif strncmpi(varargin{j},'save',4)
        savefits = 1;
    elseif strncmpi(varargin{j},'parallel',6)
        parallel = 1;
    end
    j = j+1;
end

if iscellstr(dirname)
    if ~isempty(oldfits)
        for j = 1:length(dirname)
            if isempty(oldfits{j})
                fprintf('No fits for %s\n',dirname{j});
            end
        end
    elseif parallel
        parfor j = 1:length(dirname)
            fprintf('Worker %d Doing %s\n',mygetCurrentTask('number'),dirname{j})
            try
                fits{j} =BuildRFFits(dirname{j}, varargin{:});
            catch ME
                cprintf('red','Error in %s\n',dirname{j});
                fits{j}.dirname = dirname{j};
                [a,b,c, fits{j}.name] = GetMonkeyName(dirname{j});
                fits{j}.errstate = ME;                
            end
        end
        CheckExceptions(fits);
    else
        for j = 1:length(dirname)
            fits{j} =BuildRFFits(dirname{j}, varargin{:});
        end
    end
    return;
elseif isfield(dirname,'fits') && iscell(dirname.fits)
    fits = CombineFits(dirname);
    return;
elseif iscell(dirname)
    for j = 1:length(dirname)
        fits{j} = CombineFits(dirname{j});
    end
    return;
elseif strfind(dirname,'*')
    x = mydir(dirname);
    fits = BuildRFFits({x.name},varargin{:});
    if savefits && ~isempty(fits.fits)
    save(outname, 'fits');
    end
    return;
end

    outname = [dirname  '/rffits.mat'];
    if exist(outname) && ~refit
        load(outname);
        if ~isfield(fits,'spacing') || ~isfield(fits,'probelen')
            Array = GetArrayConfig(dirname);
            if isfield(Array,'spacing') && Array.spacing > 0
                fits.spacing = Array.spacing;
                fits.probelen = max(Array.Y) .* Array.spacing;
                save(outname,'fits');
            end
        end
        if ~isfield(fits,'rf') %Add RF description from UFL too
            therf = BuildRFData(dirname,args{:});
            fits = CopyFields(fits, therf,{'rf'});
        end
        if relist
            fits = CombineFits(fits);
        end
        return;
    end
    [a,b,c, name] = GetMonkeyName(outname);
    suffixes = {'/*PP.mat' '/*PP.cell*.mat' '/*PP.mu*.mat' '/*OP.mat' '/*OPRC*.mat' '/*OP.cell*.mat' '/*OP.mu*.mat' '/*PPRC.mat' ...
        '/*OPRC.cell*mat' '/*PPRC.cell*mat' '.*OPPP.mat'};
    suffixes = {'PP' 'OP' 'PPRC' 'OPRC' 'OPPP' 'XO' 'YO'};
    fits = {};
    for j = 1:length(suffixes)
        d = mydir({[dirname '/*' suffixes{j} '.mat'] [dirname '/*' suffixes{j} '.cell*.mat'] [dirname '/*' suffixes{j} '.mu*.mat']});
        a = FitExpts(d);
        fits = {fits{:} a{:}};
    end
    therf = BuildRFData(dirname,args{:});
    if isempty(fits)
        fits.dirname = dirname;
        fits.name = name;
        fits.nfits = 0;
        fits = CopyFields(fits, therf,{'StartDepth' 'electrode' 'area' 'date' 'depth' 'rf'});
        return;
    end
    X.fits = fits;
    clear fits;
    therf = BuildRFData(dirname,args{:});
    X = CopyFields(X, therf,{'StartDepth' 'electrode' 'area' 'date' 'depth'});
    Array = GetArrayConfig(dirname);
    X = CopyFields(X, Array,{'spacing'});
    if ~isfield(X,'spacing')
        X.spacing = 0;
    end
    if isfield(Array,'idstr')
        X.electrodeidstr = Array.idstr;
    end
    if length(therf.rf > 7)
        X.Pn = therf.rf(6);
        X.Xp = therf.rf(7);
        X.Yp = therf.rf(8);
    end
        
    X.dirname = dirname;
    X.nfits = length(X.fits);
    X.name = name;
    fits = CombineFits(X);
if autosave && ~exist(outname)
    savefits = 1;
end
if savefits && ~isempty(fits.fits)
    save(outname, 'fits');
end

function F = CombineFits(F)
probe = [];
if isempty(F)
    return;
end
fits = F.fits;
for j = 1:length(fits)
    if isfield(fits{j},'pvar')
        if isnan(fits{j}.depth) && isfield(F,'depth')
            F.fits{j}.depth = F.depth;
            fits{j}.depth = F.depth;
        end
        pvar(j) = fits{j}.pvar;
        
        probe(j) = fits{j}.probe;
        meanpos(j) = fits{j}.mean;
        fitamp(j) = fits{j}.amp;
        meansd(j) = fits{j}.sd;
        expts{j} = fits{j}.expt;
        depth(j) = fits{j}.depth;
        fxs(j) = fits{j}.fx;
        fys(j) = fits{j}.fy;
        if strfind(expts{j},'Pp')
            exptype(j) = 1;
        elseif strfind(expts{j},'Op')
            exptype(j) = 2;
        elseif strfind(expts{j},'sOXor') && ~isempty(fits{j}.RFpos)
            rfpos(j,:) = fits{j}.RFpos;
            exptype(j) = 3;
        else
            exptype(j) = 0;
        end
        if fits{j}.mean > max(fits{j}.xv) || fits{j}.mean < min(fits{j}.xv)
            pvar(j) = 0;
        end
        ros(j) = fits{j}.Ro;
        if isfield(fits{j},'rOp')
            rOps(j) = fits{j}.rOp;
        end
    else
        pvar(j) = NaN;
        exptype(j) = NaN;
        probe(j) = NaN;
    end
end

if isfield(F,'spacing')
    espc = F.spacing;
    if isfield(F,'probelen')
        probelen = F.probelen(1);
    else
        probelen = max(probe) .* espc;
    end
else
    espc = 0;
    probelen = 0;
end


pendepth = nanmean(depth);
if isnan(pendepth)
    pendepth = 0;
end
if ~isempty(probe)
    proberf = [];
    fitids = {};
    probes = unique(probe);
    nrf = 0;
    for j = 1:length(probes)
        pid = find(probe == probes(j) & ismember(exptype,[1 3]) & pvar > 0.5 & fitamp > 0);
        oid = find(probe == probes(j) & ismember(exptype,[2 3]) & pvar > 0.5 & fitamp > 0);
        a = depth([oid pid]);
        meandepth = mean(a(~isnan(a)));
        if isempty(meandepth)
            a = zeros(size(a));
        else
            a(isnan(a)) = meandepth;
        end
        depth([oid pid]) = a;
        bothexp = find(exptype ==3);
        if isnan(nanmean(depth(bothexp)))
            depth(bothexp) = pendepth;
        end
        if ~isempty(oid) && ~isempty(pid)
            d = sort(depth([oid pid]));
            bid = find(diff((d)) > 0.5); %500uM break in position
            bid = [bid length(d)];
            start = -2;
            for nd = 1:length(bid)
                ro = unique(ros([pid oid]));
                for o = 1:length(ro)
                    dpid = pid(find(depth(pid) > start & depth(pid)  <= d(bid(nd)) & ros(pid) == ro(o)));
                    doid = oid(find(depth(oid) > start & depth(oid)  <= d(bid(nd)) & ros(oid) == ro(o)));
                    if ~isempty(doid) & ~isempty(dpid)
                        ppos = mean(meanpos(dpid));
                        opos = mean(meanpos(doid));
                        id = intersect(dpid,bothexp);
                        if isempty(id)  %use separate O,P estimates
                            rf = op2xy([opos ppos],ro(o));
                        else %use data from OPPP expts
                            rf = mean(rfpos(id,:));
                        end
                        rf(3) = mean(meansd(doid));
                        rf(4) = mean(meansd(dpid));
                        rf(5) = ro(o);
                        rf(6) = F.Pn;
                        rf(7) = F.Xp;
                        rf(8) = F.Yp;
                        nrf = nrf+1;
                        rf(9) = mean(fxs([doid dpid]));
                        rf(10) = mean(fys([doid dpid]));
                        rf(11) = mean(depth([doid dpid])) - probelen + probes(j).*espc;
                        proberf(nrf,:) = [rf probes(j)];
                        fitids{nrf} = [dpid doid];
                    end
                end
                start = d(bid(nd));
            end
        end
    end
    F.proberf = proberf;
    F.fitids = fitids;
    if ~isempty(proberf)
        F.Ro = mean(proberf(:,5));
        F.Rx = mean(proberf(:,1));
        F.Ry = mean(proberf(:,2));
    else
        F.Ro = NaN;
    end
end

function fits = FitExpts(d)
fits = {};
nf = 0;
for j = 1:length(d)
    if isempty(strfind(d(j).name,'.lfp.'))
        Expt = LoadExpt(d(j).name);
        if isfield(Expt,'Stimvals')
        if Expt.Stimvals.ve < 4.8       
            Expt = FillTrials(Expt,'xo');
            Expt = FillTrials(Expt,'yo');
            if Expt.Header.rc == 0
                for t = 1:length(Expt.Trials)
                    op = xy2op([Expt.Trials(t).xo Expt.Trials(t).yo],Expt.Stimvals.Ro);
                    Expt.Trials(t).Op = op(1);
                    Expt.Trials(t).Pp = op(2);
                end
            else
            end
        end
        res = PlotExpt(Expt,'noplot','fbox');
        if ~isempty(res)
            E = res.Data;
            newfit = FitExpt(res);
            if E.Header.rc && ~isempty(newfit)
                expt = E.Header.expname;
                if strfind(expt,'Pp')
                    etype = 2;
                elseif strfind(expt,'Op')
                    etype =1 ;
                end 
                opm = 0;
                ppm = 0;
                if isfield(Expt.Trials,'Op')
                    ops = unique(cat(1,Expt.Trials.Op));
                    opm = median(ops(ops > -1000));                    
                end
                if isfield(Expt.Trials,'Pp')
                    ops = unique(cat(1,Expt.Trials.Pp));
                    ppm = median(ops(ops > -1000));
                end
                
                if (ppm == 0 && etype == 2) || (opm == 0 && etype == 1)
                    if isfield(Expt.Stimvals,'Rx') && abs(Expt.Stimvals.Rx)+abs(Expt.Stimvals.Ry) > 0
                        op = xy2op([Expt.Stimvals.Rx Expt.Stimvals.Ry],Expt.Stimvals.Ro);
                    else
                        op = xy2op([Expt.Stimvals.xo Expt.Stimvals.yo],Expt.Stimvals.Ro);
                    end
                    if strcmp(res.type{1},'Pp')
                        newfit.mean = newfit.mean+op(2);
                        newfit.xv = newfit.xv+op(2);
                        newfit.x = newfit.x+op(2);
                    elseif strcmp(res.type{1},'Op')
                        newfit.mean = newfit.mean+op(1);
                        newfit.xv = newfit.xv+op(1);
                        newfit.x = newfit.x+op(1);
                    end
                elseif Expt.Stimvals.ve < 4.8 
                    if strcmp(res.type{1},'Pp')
                        newfit.mean = newfit.mean;
                    elseif strcmp(res.type{1},'Op')
                        newfit.mean = -newfit.mean;
                        newfit.xv = -newfit.xv;
                        newfit.x = -newfit.x;
                    end
                end
            end
            if isfield(E.Header,'depth')
                depth = E.Header.depth;
            elseif isfield(E.Stimvals,'ed') && Expt.Stimvals.ed > 0
                depth = E.Stimvals.ed;
            elseif isfield(E.Trials,'ed')
                depth = mean([E.Trials.ed]);
            else
                depth = NaN;
            end
            
            if isfield(E.Header,'probe')
                p = E.Header.probe;
            else
                p = GetProbeFromName(d(j).name);
            end
            for k = 1:length(newfit)
                nf = nf+1;
                fits{nf} = newfit(k);
                fits{nf}.probe = p;
                fits{nf}.name = d(j).name;
                fits{nf}.depth = depth;
                fits{nf}.expt = E.Header.expname;
                fits{nf} = CopyFields(fits{nf},E.Stimvals,{'Ro' 'xo' 'yo' 'rPp' 'rOp' 'Rx' 'Ry' 'Op' 'Pp'});
                fits{nf}.fx = GetEval(E,'fx');
                fits{nf}.fy = GetEval(E,'fy');
                if strfind(fits{nf}.expt,'Pp')
                    fits{nf}.type = 2;
                elseif strfind(fits{nf}.expt,'Op')
                    fits{nf}.type = 1;
                end
            end
        end
        end
    end
end

